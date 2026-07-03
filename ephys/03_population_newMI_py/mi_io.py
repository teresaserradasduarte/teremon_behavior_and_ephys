"""mi_io.py - lean shared loader for exploratory MI work.
Parse the .mat ONCE into per-neuron dicts (spike times kept so counts can be
recomputed for any window), cache to pickle, iterate fast. Path-agnostic: callers
pass the .mat path and the struct-variable names (see config.py).

Supports both classic .mat (<v7.3, scipy) and v7.3 HDF5 .mat (h5py); detected
automatically from the file magic bytes.

reach_px lives in a SEPARATE cache file (reach_px.pkl), keyed by (mouse, sess).
Shape convention: (n_trials, 3, 201)  [trials × coord × frames]
  coord: 0=x, 1=y med-lat, 2=z
  Endpoint frames 147-167 MATLAB 1-idx → Python reach_px[:, :, 146:167]
Neuron dicts contain only spike times + labels (keeps main cache lean).
"""
import numpy as np, pickle
from pathlib import Path

LABELS = ["idx_reach_hit", "idx_reach_succ", "idx_reach_LCR", "idx_reach_PP", "idx_reach_cat"]

def _short(s):
    return s.replace("all_", "").replace("_sts", "")

def _is_hdf5(path):
    with open(path, "rb") as f:
        return f.read(10) == b"MATLAB 7.3"

# ---------------------------------------------------------------------------
# HDF5 (.mat v7.3) parser
# ---------------------------------------------------------------------------
def _str_from_ref(f, ref):
    return "".join(chr(c) for c in np.asarray(f[ref]).ravel())

def _spikes_from_ref(f, ref):
    out = []
    for tr in np.asarray(f[ref]).ravel():
        spk = np.asarray(f[tr], float).ravel()
        out.append(spk[np.isfinite(spk)])
    return out

def _label_from_ref(f, ref, n):
    v = np.asarray(f[ref], float).ravel()
    return v if v.size == n else np.full(n, np.nan)

def _parse_mat_h5(path, regions):
    import h5py
    neurons, reach_px = [], {}

    with h5py.File(str(path), "r") as f:
        for reg in regions:
            if reg not in f:
                print("missing:", reg); continue
            g = f[reg]
            n_neu = g["st_reach"].shape[1]

            for i in range(n_neu):
                if i % 100 == 0:
                    print(f"  {reg}: {i}/{n_neu}", flush=True)

                trials = _spikes_from_ref(f, g["st_reach"][0, i])
                n = len(trials)

                lab = {lf: (_label_from_ref(f, g[lf][0, i], n) if lf in g
                            else np.full(n, np.nan)) for lf in LABELS}

                mouse = _str_from_ref(f, g["mouse"][0, i])
                sess  = _str_from_ref(f, g["sess"][0, i])
                phyid = int(np.asarray(f[g["phyID"][0, i]], float).ravel()[0])

                # reach_px: one array per session, not per neuron
                key = (mouse, sess)
                if key not in reach_px:
                    rp = np.asarray(f[g["reach_px"][0, i]], float)
                    reach_px[key] = rp if (rp.ndim == 3 and rp.shape[0] == n) \
                                    else np.full((n, 3, 201), np.nan)

                neurons.append(dict(
                    uid=f"{_short(reg)}|{mouse}|{sess}|{phyid}",
                    region=_short(reg), mouse=mouse, sess=sess,
                    phyID=phyid, n_trials=n, spikes=trials, labels=lab))

    return neurons, reach_px

# ---------------------------------------------------------------------------
# Classic .mat (<v7.3) parser
# ---------------------------------------------------------------------------
def _spk(x):
    a = np.atleast_1d(np.asarray(x, float)).ravel(); return a[np.isfinite(a)]

def _trials_list(st):
    return list(st) if isinstance(st, np.ndarray) and st.dtype == object else [st]

def _parse_mat_scipy(path, regions):
    import scipy.io as sio
    mat = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    neurons, reach_px = [], {}
    for reg in regions:
        if reg not in mat:
            print("missing:", reg); continue
        for neu in np.atleast_1d(mat[reg]):
            trials = [_spk(t) for t in _trials_list(neu.st_reach)]
            n = len(trials)
            lab = {}
            for lf in LABELS:
                v = np.atleast_1d(np.asarray(getattr(neu, lf, np.nan), float)).ravel()
                lab[lf] = v if v.size == n else np.full(n, np.nan)
            key = (str(neu.mouse), str(neu.sess))
            if key not in reach_px:
                rp_raw = np.asarray(getattr(neu, "reach_px", np.full((201, 3, n), np.nan)), float)
                reach_px[key] = rp_raw.transpose(2, 1, 0) if rp_raw.ndim == 3 \
                                 else np.full((n, 3, 201), np.nan)
            neurons.append(dict(
                uid=f"{_short(reg)}|{neu.mouse}|{neu.sess}|{int(neu.phyID)}",
                region=_short(reg), mouse=str(neu.mouse), sess=str(neu.sess),
                phyID=int(neu.phyID), n_trials=n, spikes=trials, labels=lab))
    return neurons, reach_px

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def parse_mat(path, regions, **_kwargs):
    """Returns (neurons, reach_px).
    neurons: list of per-neuron dicts (spikes + labels, no reach_px).
    reach_px: dict keyed by (mouse, sess) → (n_trials, 3, 201) array.
    Save both with mi_io.save(); load both with mi_io.load()."""
    path = Path(path)
    if _is_hdf5(path):
        return _parse_mat_h5(path, regions)
    return _parse_mat_scipy(path, regions)

def endpoint_kin(neuron, reach_px, coord=1, frames=(147, 167)):
    """Mean kinematic value over endpoint window.
    reach_px: the dict returned by parse_mat / load.
    coord: 0=x, 1=y med-lat, 2=z.  frames: MATLAB 1-indexed."""
    rp = reach_px[(neuron["mouse"], neuron["sess"])]   # (n_trials, 3, 201)
    return rp[:, coord, frames[0]-1:frames[1]].mean(axis=1)

def counts(neuron, win, bin_size=0.0025):
    # bin_size kept for API compatibility; summing bins = counting spikes in window
    lo, hi = win
    return np.array([int(((s >= lo) & (s < hi)).sum()) for s in neuron["spikes"]])

def save(neurons, reach_px, path):
    """Save neurons list and reach_px dict as two separate files."""
    path = Path(path)
    with open(path, "wb") as f:
        pickle.dump(neurons, f, protocol=4)
    with open(path.with_suffix(".reach_px.pkl"), "wb") as f:
        pickle.dump(reach_px, f, protocol=4)

def load(path):
    """Returns (neurons, reach_px) — mirrors parse_mat."""
    path = Path(path)
    with open(path, "rb") as f:
        neurons = pickle.load(f)
    with open(path.with_suffix(".reach_px.pkl"), "rb") as f:
        reach_px = pickle.load(f)
    return neurons, reach_px
