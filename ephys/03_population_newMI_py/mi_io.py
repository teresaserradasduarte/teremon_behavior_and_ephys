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

# EMPTY ARRAYS IN A v7.3 .mat. HDF5 cannot hold MATLAB's 0x1 empty double as
# such, so MATLAB writes a small dataset containing the array's DIMENSIONS --
# the integers [0, 1] -- and marks it with the attribute MATLAB_empty = 1.
# Reading that dataset as numbers returns the SIZE of the empty array as if it
# were its content. Until 2026-10-12 _spikes_from_ref did exactly that, so every
# trial with no spikes came back as the spike times [0.0, 1.0]: 6-9% of CP
# trials, 0.3-0.4% of CB, and not one empty trial anywhere. The flag is now
# honoured, and load() repairs caches parsed before the fix.
# (19_spike_placeholder.py is the figure that shows it.)
def _is_matlab_empty(d):
    return bool(d.attrs.get("MATLAB_empty", 0))

def _spikes_from_ref(f, ref):
    cell = f[ref]
    if _is_matlab_empty(cell):                 # a neuron with no trials at all
        return []
    out = []
    for tr in np.asarray(cell).ravel():
        d = f[tr]
        if _is_matlab_empty(d):                # a trial with no spikes: keep the
            out.append(np.empty(0, float))     # TRIAL, with nothing in it
            continue
        spk = np.asarray(d, float).ravel()
        out.append(spk[np.isfinite(spk)])
    return out

def _label_from_ref(f, ref, n):
    d = f[ref]
    if _is_matlab_empty(d):                    # same trap: [0, 1] is a size
        return np.full(n, np.nan)
    v = np.asarray(d, float).ravel()
    return v if v.size == n else np.full(n, np.nan)

def repair_empty_trials(neurons):
    """Turn trials that are exactly [0.0, 1.0] back into EMPTY trials, in place.
    -> how many were repaired.

    For pickles parsed before _spikes_from_ref honoured MATLAB_empty: those hold
    the mis-read dimensions of an empty array where a silent trial should be.
    The trial is emptied, never removed, so every neuron keeps its trial count
    and order and trial k still pairs with reach_px[k] and labels[k].
    The value rule is safe: over the whole file it selects exactly the trials
    the MATLAB_empty flag marks, for st_init and st_reach, in CP and CB."""
    n = 0
    for neu in neurons:
        spk = neu["spikes"]
        for k, t in enumerate(spk):
            if len(t) == 2 and t[0] == 0.0 and t[1] == 1.0:
                spk[k] = np.empty(0, float); n += 1
    return n

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
    # A cache parsed before the MATLAB_empty fix is repaired here, in memory, so
    # every script gets empty trials as empty without re-parsing the .mat. A
    # cache parsed after the fix has nothing to repair and prints nothing.
    n_fix = repair_empty_trials(neurons)
    if n_fix:
        n_tr = sum(len(n["spikes"]) for n in neurons)
        print(f"mi_io.load: {n_fix}/{n_tr} trials ({100*n_fix/max(n_tr,1):.1f}%) "
              f"were empty trials mis-read as [0.0, 1.0] -> restored to EMPTY "
              f"(kept as trials). Re-run 01_load_data.py to fix the cache itself.")
    return neurons, reach_px


def describe(obj, name="obj", max_depth=3, max_items=12, _depth=0, _key=None):
    """Print the STRUCTURE of a nested dict/list/array -- what is in a pickle,
    without expanding the VS Code Variables pane one node at a time.

        import mi_io; mi_io.describe(out)
        mi_io.describe(WATER, max_items=3)

    Shows type, shape, dtype, and for numeric arrays the finite range and NaN
    fraction (the thing you usually actually want to check).
    """
    import numpy as np
    pad = "  " * _depth
    tag = f"{_key}: " if _key is not None else f"{name}: "

    if isinstance(obj, dict):
        print(f"{pad}{tag}dict({len(obj)} keys)")
        if _depth >= max_depth:
            return
        for i, (k, v) in enumerate(obj.items()):
            if i >= max_items:
                print(f"{pad}  ... {len(obj)-max_items} more keys")
                break
            describe(v, max_depth=max_depth, max_items=max_items,
                     _depth=_depth+1, _key=repr(k) if not isinstance(k, str) else k)
        return

    if isinstance(obj, (list, tuple)):
        print(f"{pad}{tag}{type(obj).__name__}(len {len(obj)})")
        if len(obj) and _depth < max_depth:
            describe(obj[0], max_depth=max_depth, max_items=max_items,
                     _depth=_depth+1, _key="[0]")
        return

    if isinstance(obj, np.ndarray):
        d = f"{pad}{tag}ndarray {obj.shape} {obj.dtype}"
        if obj.size and obj.dtype.kind in "fiu":
            a = obj[np.isfinite(obj)] if obj.dtype.kind == "f" else obj
            nan = float(np.mean(~np.isfinite(obj))) if obj.dtype.kind == "f" else 0.0
            if a.size:
                d += f"  range [{a.min():.4g}, {a.max():.4g}]"
            if nan:
                d += f"  NaN {nan*100:.1f}%"
        print(d)
        return

    v = repr(obj)
    print(f"{pad}{tag}{type(obj).__name__} = {v[:70]}{'...' if len(v) > 70 else ''}")
