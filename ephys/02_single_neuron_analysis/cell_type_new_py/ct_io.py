"""ct_io.py - reading the data (all read-only) and the small csv/pkl tables."""
import csv, os, pickle, re
import numpy as np
import h5py
import ct_config as C


# ── sessions and sorter folders ──────────────────────────────────────────────
def sessions(kinds=("CB", "BG"), mice=None, sess=None):
    """-> [(mouse, 'R3', kind, path of eg_neurons_<kind>.mat)] for the files that exist."""
    out = []
    for m in (mice or C.MICE):
        for s in range(1, C.N_SESS + 1):
            S = f"R{s}"
            if sess and S not in sess: continue
            for k in kinds:
                p = C.MAT_DIR / m / S / f"eg_neurons_{k}.mat"
                if p.is_file(): out.append((m, S, k, p))
    return out


_PATHS = None
def sorter_folder(mouse, sess, kind):
    """Curated sorter folder (kilosort4 / ibl_sorter_results...) of one session, from
    teremon_paths.mat -- the list s04_CB_cellType.m uses. None if not listed."""
    global _PATHS
    if _PATHS is None:
        from scipy.io import loadmat
        m = loadmat(str(C.PATHS_MAT), squeeze_me=True)
        _PATHS = {"CB": [str(x) for x in np.atleast_1d(m["files_curation_CB"])],
                  "BG": [str(x) for x in np.atleast_1d(m["files_curation_BG"])]}
    name = mouse.split("_", 1)[1].lower()
    n = sess[1:]
    hits = []
    for p in _PATHS[kind]:
        parts = [x for x in re.split(r"[\\/]+", p) if x]
        try: i = next(j for j, x in enumerate(parts) if x.lower() == mouse.lower())
        except StopIteration:
            if name not in p.lower(): continue
            i = -1
        rec = parts[i + 1] if i >= 0 else p                 # '04082023_CoteDor_Cer_S5_g0'
        if re.search(rf"S{n}(?=[A-Z_])", rec): hits.append(p)
    if len(hits) > 1:
        raise ValueError(f"{mouse} {sess} {kind}: {len(hits)} sorter folders match: {hits}")
    return hits[0] if hits else None


# ── eg_neurons_<kind>.mat (v7.3) ─────────────────────────────────────────────
def _s(f, ref):
    d = f[ref]
    if d.attrs.get("MATLAB_empty", 0): return ""
    a = np.asarray(d)
    return "".join(chr(int(c)) for c in a.ravel()) if a.dtype.kind in "ui" else ""


def _num(f, ref):
    d = f[ref]
    if d.attrs.get("MATLAB_empty", 0): return np.empty(0)
    return np.asarray(d, float)


def load_eg_neurons(path, regions=None, with_wf=True, with_st=True):
    """One session file -> list of dicts: phyID, brain_region, depth, st (s, sorted,
    whole session), W (n_channels x n_samples template). regions: brain_region filter."""
    out = []
    with h5py.File(str(path), "r") as f:
        g = f["eg_neurons"]
        n = g["phyID"].shape[1]
        for i in range(n):
            br = _s(f, g["brain_region"][0, i]).strip()
            if regions and br not in regions: continue
            d = dict(phyID=int(_num(f, g["phyID"][0, i]).ravel()[0]), brain_region=br,
                     depth=float(_num(f, g["depth"][0, i]).ravel()[0]) if "depth" in g else np.nan)
            if with_st:
                st = _num(f, g["st"][0, i]).ravel()
                d["st"] = np.sort(st[np.isfinite(st)])
            if with_wf:
                W = _num(f, g["templateWaveforms"][0, i])
                if W.ndim == 3:        # a unit merged from several templates: (n_templates,
                    wt = _num(f, g["templateWeight"][0, i]).ravel()    # channels, samples)
                    wt = wt / wt.sum() if wt.size == W.shape[0] and wt.sum() > 0 \
                        else np.full(W.shape[0], 1 / W.shape[0])
                    W = np.tensordot(wt, W, axes=1)                    # weighted mean template
                    d["n_templates"] = int(wt.size)
                if W.ndim == 2 and W.shape[0] < W.shape[1]: W = W.T    # -> channels x samples
                d["W"] = W
            out.append(d)
    return out


# ── s04's labels (unchanged, read only) ──────────────────────────────────────
def old_cb_labels():
    """-> {uid: label} from eg_neurons_CB_cellType.mat (all_CB_cellID). 'NaN'/numeric
    NaN/empty -> 'unlabelled'."""
    out = {}
    with h5py.File(str(C.OLD_CB_LABELS_MAT), "r") as f:
        g = f["all_CB_cellID"]
        for i in range(g["phyID"].shape[1]):
            ct = _s(f, g["cell_type"][0, i]).strip()
            if ct in ("", "NaN", "nan") or not ct.isprintable(): ct = C.UNLAB
            out[C.uid("CB", _s(f, g["mouse"][0, i]), _s(f, g["session"][0, i]),
                      _num(f, g["phyID"][0, i]).ravel()[0])] = ct
    return out


# ── C4 output of one sorter folder ───────────────────────────────────────────
def _tsv(path):
    if not os.path.isfile(C.lp(path)): return {}
    with open(C.lp(path), newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return {int(r[0]): r[1] for r in rows[1:] if len(r) > 1 and r[0].strip().isdigit()}


def read_c4(folder):
    """C4 files of a sorter folder -> dict:
      tsv      {cluster_id: dict(pred, prob, ratio, votes)}  the units C4 REPORTED
               (confidence ratio >= its threshold, 2 by default)
      classes  class names in the order of the ensemble's class axis
      ens      (n_units, n_classes, n_models) raw probabilities of EVERY unit C4 ran,
               including the ones below threshold; None if the file is missing
      ens_ids  cluster_id of each ensemble row, -1 where it could not be recovered.
    C4 does not save which unit each ensemble row is. Rows are in the order of the
    'good' units; the reported units are found by matching their saved probability /
    ratio, and an unreported row is identified only when the number of 'good',
    unreported cluster ids lying between its matched neighbours equals the number of
    unmatched rows there."""
    j = lambda *a: os.path.join(folder, *a)
    pred = _tsv(j("cluster_predicted_cell_type.tsv"))
    prob = _tsv(j("cluster_pred_probability.tsv"))
    ratio = _tsv(j("cluster_confidence_ratio.tsv"))
    votes = _tsv(j("cluster_model_votes.tsv"))
    def _v(s):
        try: a, b = s.split("/"); return float(a) / float(b)
        except Exception: return np.nan
    tsv = {k: dict(pred=pred[k], prob=float(prob.get(k, "nan")), ratio=float(ratio.get(k, "nan")),
                   votes=_v(votes.get(k, ""))) for k in pred}
    out = dict(tsv=tsv, classes=None, ens=None, ens_ids=None)
    pe = j("cell_type_classification", "ensemble_predictions_ncells_nclasses_nmodels.npy")
    pl = j("cell_type_classification", "label_correspondence.tsv")
    if not (os.path.isfile(C.lp(pe)) and os.path.isfile(C.lp(pl))): return out
    ens = np.load(C.lp(pe))
    with open(C.lp(pl)) as fh:
        rows = [r.split("\t") for r in fh.read().strip().splitlines()]
    lab = {int(a): b for a, b in zip(rows[0], rows[1])}
    classes = [lab.get(i, f"class{i}") for i in range(ens.shape[1])]
    mean = ens.mean(axis=2)
    srt = np.sort(mean, axis=1)
    e_prob, e_ratio = srt[:, -1], srt[:, -1] / np.maximum(srt[:, -2], 1e-12)
    ids = np.full(ens.shape[0], -1, int)
    r0 = 0
    for cid in sorted(tsv):                                 # same order as the rows
        t = tsv[cid]
        for r in range(r0, ens.shape[0]):
            if (classes[int(mean[r].argmax())] == t["pred"]
                    and abs(e_prob[r] - t["prob"]) < 1e-4 * max(1, t["prob"])
                    and abs(e_ratio[r] - t["ratio"]) < 1e-3 * max(1, t["ratio"])):
                ids[r] = cid; r0 = r + 1; break
    good = sorted(k for k, v in _tsv(j("cluster_group.tsv")).items() if v == "good" and k not in tsv)
    known = [(-1, -1)] + [(r, ids[r]) for r in range(ens.shape[0]) if ids[r] >= 0] \
            + [(ens.shape[0], np.inf)]
    for (ra, ca), (rb, cb) in zip(known[:-1], known[1:]):
        gap = list(range(ra + 1, rb))
        cand = [c for c in good if ca < c < cb]
        if gap and len(gap) == len(cand): ids[gap] = cand
    out.update(classes=classes, ens=ens, ens_ids=ids)
    return out


# ── small tables ─────────────────────────────────────────────────────────────
def write_csv(path, rows, cols=None):
    cols = cols or list(rows[0].keys())
    tmp = str(C.out_path(path)) + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    os.replace(tmp, str(path))


def read_csv(path):
    if not os.path.isfile(str(path)): return []
    with open(str(path), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def save_pkl(path, obj):
    tmp = str(C.out_path(path)) + ".tmp"
    with open(tmp, "wb") as fh: pickle.dump(obj, fh, protocol=4)
    os.replace(tmp, str(path))


def load_pkl(path):
    with open(str(path), "rb") as fh: return pickle.load(fh)


def load_features():
    """features.pkl (01) -> dict(rows=[dict per neuron], arrays={name: (n, ...)}, meta)."""
    if not C.FEATURES_PKL.is_file():
        raise FileNotFoundError(f"{C.FEATURES_PKL} not found -- run 01_extract_features.py")
    return load_pkl(C.FEATURES_PKL)


def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return np.nan
