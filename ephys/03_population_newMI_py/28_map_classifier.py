# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# CLASSIFY THE 2D MI MAP SHAPES FROM YOUR HAND LABELS (27_label_gui.py)
#
# Your labels are the ground truth. A simple, inspectable classifier learns them:
#   features per map  (the same for real and shuffle maps)
#       25  best template r per family (column, diagonal, blob)
#       23  hot cluster: has one, extent E, orientation O, dominance D, size
#       own strength  median z, top-2.5% smoothed z, share of cells z > 2
#           structure share of the map's variance explained by its column means
#                     (column), row means (row), lag-diagonal means (diagonal);
#                     number of hot clusters and the 2nd cluster's share
#           layout    6 x 6 block means of the smoothed z-map, divided by the map's
#                     top z (the shape without the strength)   [FEATURE_SET "full"]
#   model    multinomial logistic regression, standardised features, classes
#            weighted equally, regularisation chosen by an inner cross-validation
#   check    repeated stratified 5-fold cross-validation -> confusion matrix;
#            23's and 25's classes scored against the same labels
#   apply    every map: % per class, corrected for the classifier's own errors;
#            the shuffle maps give its false-positive rate
# "mixed/unsure" maps are left out of training.
import os, sys, time, pickle, csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, label as nd_label
from scipy.optimize import nnls
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import balanced_accuracy_score, cohen_kappa_score

_HERE = (os.path.dirname(os.path.abspath(__file__))
         if "__file__" in globals() else os.getcwd())
if not os.path.isfile(os.path.join(_HERE, "mi_io.py")):
    _sub = os.path.join("ephys", "03_population_newMI_py")
    _p, _found = os.getcwd(), None
    while _found is None:
        for _c in (_p, os.path.join(_p, _sub)):
            if os.path.isfile(os.path.join(_c, "mi_io.py")): _found = _c; break
        if _found is None:
            if _p == os.path.dirname(_p):
                raise ModuleNotFoundError("mi_io.py not found from " + os.getcwd())
            _p = os.path.dirname(_p)
    _HERE = _found
if _HERE not in sys.path: sys.path.insert(0, _HERE)
os.chdir(_HERE)
import config

# ── inputs (the same run as 26) ──────────────────────────────────────────────
TAG        = "x-y-z-speed-dist_water_nw70_cat1"
SRC_SET    = f"label_set_{TAG}.pkl"
SRC_LABELS = f"labels_{TAG}.csv"
SRC_12     = "td2D_v2_x-y-z-speed-dist_water_nw70_kw21_ks4_nall_s20_cat1_2Dcache.pkl"
SRC_23     = "shapes2D_x-y-z-speed-dist_water_nw70_ks4_cat1_nall_shufmap_cache.pkl"
SRC_23_RES = "shapes2D_x-y-z-speed-dist_water_nw70_ks4_cat1_nall_results.pkl"
SRC_25_RES = "tmpl2D_x-y-z-speed-dist_water_nw70_ks4_cat1_nall_results.pkl"

# ── classes and features ─────────────────────────────────────────────────────
CLASSES      = ("column", "diagonal", "tiling", "broad", "nothing")
# "2 tiles" is not a class: tiling maps get 1 or 2 hotspots instead (cell 3b).
MERGE_2TILES = False   # your "2 tiles" labels: True = trained as tiling, False = left out of training
FEATURE_SET  = "full"  # final model: "summary" (~20 features) or "full" (+ 6x6 layout)
N_BLOCK      = 6       # layout blocks per axis
SMOOTH_SIG   = (1.0, 0.5)   # (kin, neural) cells, as 23 / 26
Z_CLIP       = 20.0    # z clipped to +-Z_CLIP before the own features (as 25)
Z_HOT        = 2.0     # hot cells: smoothed z > Z_HOT (as 23)
MIN_CLUSTER  = 4       # cells, as 23
TOP_FRAC     = 0.025   # as 23
LAG_BIN      = 0.025   # s, lag bins for the diagonal structure
MIN_FINITE   = 0.6     # as 26: sparser maps are not classified

# ── cross-validation ─────────────────────────────────────────────────────────
N_FOLDS  = 5
N_REPS   = 10          # repeats of the 5-fold split (different random splits)
CS       = np.logspace(-3, 1, 9)   # regularisation grid, chosen inside each fold
SEED     = 0

# ── confidence (cell 3b) ─────────────────────────────────────────────────────
CONF_MIN      = 0.4    # calls below this probability -> "low confidence" (0: off); was 0.5
LOW = "low confidence"

# ── hotspots of the tiling maps (cell 3b) ────────────────────────────────────
HOT_REL       = 0.5    # hotspot level, as a fraction of the map's peak smoothed z
MIN_HOT_CELLS = 4
HOT_RANK      = "shape"  # "shape" (size x strength x roundness) or "peak" (highest peak z)
ROUND_PAD     = 0.3    # cells, added to both SDs so tiny patches are not called lines
SECOND_PEAK   = 0.7    # 2nd hotspot: share of the 1st one's peak z
SECOND_SIZE   = 0.5    # 2nd hotspot: share of the 1st one's cells
MIN_SEP       = 0.15   # s, minimum distance between the two centres
DIP_MAX       = 0.35   # the valley between them, as a share of the 2nd peak (0.5 let
                       # two lobes of one blob through)

# ── times (cells 8, 9) ───────────────────────────────────────────────────────
P_MIN         = 0.0    # extra cut for the times only: probability >= P_MIN (on top of CONF_MIN)
N_CHECK       = 8      # example maps per class
MIN_LAG_CELLS = 5      # lag diagonals with fewer defined cells are ignored (map corners)
STANDOUT_RANGE     = (0.4, 0.85)  # cell 9c: stand-out at or below / above these -> faintest / opaque
STANDOUT_STYLE     = "both"    # cell 9c: stand-out shown as "opacity", "size" or "both"
STANDOUT_SIZE      = (4, 40)   # cell 9c: dot size (points^2) of the faintest / most standing-out
STANDOUT_ALPHA_MIN = 0.25  # cell 9c: opacity of a column/diagonal hotspot that does not stand out

# ── outputs ──────────────────────────────────────────────────────────────────
N_EX       = 8         # example maps per class
SAVE_VARS  = True
SAVE_FIGS  = True
SHOW_FIGS  = True
OUT_SUBDIR = "MI_2D_classifier"
RUN_NAME   = "clf2D"

CLR = {"CP": np.array([106, 176, 134]) / 255, "CB": np.array([154, 123, 193]) / 255}
CLASS_CLR = {"column": "#3b82f6", "diagonal": "#f59e0b", "tiling": "#10b981",
             "2 tiles": "#ec4899", "broad": "#8b5cf6", "nothing": "#9ca3af",
             "mixed/unsure": "#a16207", "unsure": "#a16207"}
# 23's and 25's classes on the label scale (only to score them against your labels)
COARSE = {"column": "column", "diagonal": "diagonal", "tiling": "tiling", "broad": "broad",
          "not informative": "nothing", "noise-sized": "nothing", "no cluster": "nothing",
          "scattered": "nothing", "ambiguous": "unsure"}

T1, T2 = "tiling, 1 hotspot", "tiling, 2 hotspots"      # for the times (cells 8, 9)
CLASS_CLR[T1], CLASS_CLR[T2] = CLASS_CLR["tiling"], CLASS_CLR["2 tiles"]
K = len(CLASSES)


def _savefig_retry(fig, path, tries=6, **kw):
    """Dropbox briefly locks files it is uploading: retry instead of crashing."""
    for t in range(tries):
        try:
            return fig.savefig(path, **kw)
        except OSError:
            if t == tries - 1: raise
            time.sleep(0.5 * (t + 1))


def _long_path(p):
    """Windows refuses paths of 260+ characters unless they start with \\\\?\\"""
    p = os.path.abspath(str(p))
    if os.name == "nt" and len(p) >= 260 and not p.startswith("\\\\"):
        return "\\\\?\\" + p
    return p


RUN_TAG = f"{TAG}_{FEATURE_SET}{'_merged' if MERGE_2TILES else ''}"
fdir = config.fig_dir(OUT_SUBDIR, RUN_TAG) if SAVE_FIGS else None


def _save(fig, name):
    if SAVE_FIGS:
        _savefig_retry(fig, _long_path(fdir / f"{RUN_NAME}_{name}.png"), dpi=150,
                       bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 1. load: labels, maps, 23's and 25's results ───────────────────────────
# 1. cell --> your labels (first pass), 12's maps, 23's shuffle maps, 23 and 25.
with open(config.var_path(SRC_SET), "rb") as fh: LSET = pickle.load(fh)
with open(config.var_path(SRC_LABELS), newline="", encoding="utf-8") as fh:
    LROWS = list(csv.DictReader(fh))
with open(config.var_path(SRC_12), "rb") as fh: C12 = pickle.load(fh)
with open(config.var_path(SRC_23), "rb") as fh: C23 = pickle.load(fh)
if C23["key"].get("src") != SRC_12:
    raise ValueError(f"23's shuffle maps were not made from {SRC_12}")
with open(config.var_path(SRC_23_RES), "rb") as fh: R23 = pickle.load(fh)
with open(config.var_path(SRC_25_RES), "rb") as fh: R25 = pickle.load(fh)

win_centers, kin_centers = np.asarray(LSET["win_centers"]), np.asarray(LSET["kin_centers"])
if C12["mi"].shape[1:] != (len(kin_centers), len(win_centers)):
    raise ValueError("12's maps do not match the label set's grid")
Kc, Nc = np.meshgrid(kin_centers, win_centers, indexing="ij")
LAGc = Nc - Kc
extent = [win_centers[0], win_centers[-1], kin_centers[0], kin_centers[-1]]

uids = list(C12["uids"]); regions = np.asarray(C12["regions"]).astype(str)
ROW = {u: i for i, u in enumerate(uids)}
_row23 = {u: i for i, u in enumerate(C23["key"]["uids"])}
_i23 = {u: i for i, u in enumerate(R23["uids"])}
_i25 = {u: i for i, u in enumerate(R25["uids"])}
missing = [u for u in uids if u not in _row23 or u not in _i23 or u not in _i25]
print(f"{len(uids)} maps in 12; {len(missing)} not in 23's / 25's results (left out)")

FIRST = [r for r in LROWS if r["pass"] == "1"]
REPEAT_ROWS = [r for r in LROWS if r["pass"] == "repeat"]
print(f"labels: {len(FIRST)} first pass, {len(REPEAT_ROWS)} repeat pass")
import collections
print("  " + ", ".join(f"{k} {v}" for k, v in collections.Counter(r["label"] for r in FIRST).most_common()))


# %% ── 2. the features ───────────────────────────────────────────────────────
# 2. cell --> one feature vector per map, for every real map and its shuffle map.
def _smooth(Zm):
    ok = np.isfinite(Zm)
    if not ok.any(): return Zm.copy()
    num = gaussian_filter(np.where(ok, Zm, 0.0), SMOOTH_SIG)
    den = gaussian_filter(ok.astype(float), SMOOTH_SIG)
    with np.errstate(all="ignore"):
        return np.where(ok, num / den, np.nan)


def _r2_groups(Zm, groups):
    """share of the map's variance explained by the mean of each group of cells"""
    ok = np.isfinite(Zm)
    v, g = Zm[ok], groups[ok]
    if v.size < 10 or v.var() == 0: return 0.0
    _, inv = np.unique(g, return_inverse=True)
    means = np.bincount(inv, v) / np.bincount(inv)
    return float(means[inv].var() / v.var())


def _top(Zm):
    v = Zm[np.isfinite(Zm)]
    if v.size < 10: return np.nan
    return float(np.sort(v)[-max(1, int(round(TOP_FRAC * v.size))):].mean())


def _clusters(Zm):
    """-> number of hot clusters (>= MIN_CLUSTER cells), share of the 2nd largest in hot mass"""
    w = np.where(np.isfinite(Zm), Zm - Z_HOT, 0.0)
    lab, n = nd_label(w > 0, structure=np.ones((3, 3)))
    if n == 0: return 0, 0.0
    size = np.bincount(lab.ravel())[1:]
    mass = np.bincount(lab.ravel(), np.clip(w, 0, None).ravel())[1:]
    big = size >= MIN_CLUSTER
    if not big.any(): return 0, 0.0
    m = np.sort(mass[big])[::-1]
    return int(big.sum()), float(m[1] / m.sum()) if m.size > 1 else 0.0


_KB = np.array_split(np.arange(len(kin_centers)), N_BLOCK)
_NB = np.array_split(np.arange(len(win_centers)), N_BLOCK)
_COLG = np.broadcast_to(np.arange(len(win_centers))[None, :], Kc.shape)
_ROWG = np.broadcast_to(np.arange(len(kin_centers))[:, None], Kc.shape)
_LAGG = np.round(LAGc / LAG_BIN).astype(int)

SUMMARY_NAMES = ["r column (25)", "r diagonal (25)", "r blob (25)",
                 "has cluster (23)", "E (23)", "O (23)", "D (23)", "log size (23)",
                 "median z", "top z", "share z>2",
                 "R2 columns", "R2 rows", "R2 diagonals", "n clusters", "2nd cluster share"]
BLOCK_NAMES = [f"block k{a} n{b}" for a in range(N_BLOCK) for b in range(N_BLOCK)]


def _features(M, NU, SD, shape23, r25):
    """M: a map (real MI or 23's shuffle map), NU / SD: its null mean / SD;
    shape23: 23's cluster dict for that map; r25: 25's best r per family"""
    with np.errstate(all="ignore"):
        z = (M - NU) / np.where(SD > 0, SD, np.nan)
    z = np.clip(z, -Z_CLIP, Z_CLIP)
    zs = _smooth(z)
    has = np.isfinite(shape23["E"])
    ncl, m2 = _clusters(zs)
    top = _top(zs)
    with np.errstate(all="ignore"):
        summ = [r25["column"], r25["diagonal"], r25["blob"],
                float(has), shape23["E"] if has else 0.0, shape23["O"] if has else 0.0,
                shape23["D"] if has else 0.0, np.log1p(shape23["size"]),
                np.nanmedian(z), top, np.nanmean(z[np.isfinite(z)] > Z_HOT),
                _r2_groups(zs, _COLG), _r2_groups(zs, _ROWG), _r2_groups(zs, _LAGG), ncl, m2]
        blocks = [np.nanmean(zs[np.ix_(kb, nb)]) for kb in _KB for nb in _NB]
    blocks = np.nan_to_num(np.array(blocks) / max(top if np.isfinite(top) else 1.0, 1.0))
    return np.nan_to_num(np.array(summ, float)), blocks


OK = np.zeros(len(uids), bool)
FS = np.zeros((len(uids), len(SUMMARY_NAMES))); FB = np.zeros((len(uids), N_BLOCK ** 2))
FS_N, FB_N = np.zeros_like(FS), np.zeros_like(FB)
t0 = time.time()
for i, u in enumerate(uids):
    if u in missing: continue
    a, b = _i23[u], _i25[u]
    if R23["sparse"][a] or R25["sparse"][b] or str(R23["cls"][a]) == "too sparse": continue
    with np.errstate(all="ignore"):
        if np.isfinite(C12["mi"][i] - C12["null"][i]).mean() < MIN_FINITE: continue
    OK[i] = True
    FS[i], FB[i] = _features(C12["mi"][i], C12["null"][i], C12["null_std"][i], R23["shape"][a],
                             {f: R25["best_r"][f][b] for f in ("column", "diagonal", "blob")})
    FS_N[i], FB_N[i] = _features(C23["sh"][_row23[u]], C12["null"][i], C12["null_std"][i],
                                 R23["shape_null"][a],
                                 {f: R25["best_r_null"][f][b] for f in ("column", "diagonal", "blob")})
print(f"features for {int(OK.sum())} maps (+ their shuffle maps) in {time.time() - t0:.0f}s")


def _X(fs, fb, which):
    return fs if which == "summary" else np.hstack([fs, fb])


def _names(which):
    return SUMMARY_NAMES if which == "summary" else SUMMARY_NAMES + BLOCK_NAMES


# the labelled table
_lab = lambda l: "tiling" if (MERGE_2TILES and l == "2 tiles") else l
TR = []                       # dicts: row (in 12), kind, label, region, c23, c25
for r in FIRST:
    if r["label"] == "mixed/unsure": continue
    if r["label"] == "2 tiles" and not MERGE_2TILES: continue
    i = ROW.get(r["uid"])
    if i is None or not OK[i]:
        print(f"  labelled map not usable, skipped: {r['key']}"); continue
    TR.append(dict(i=i, kind=r["kind"], label=_lab(r["label"]), region=regions[i],
                   c23=COARSE.get(str(R23["cls"][_i23[r["uid"]]]), "?"),
                   c25=COARSE.get(str(R25["cls"][_i25[r["uid"]]]), "?"), uid=r["uid"]))
_shuf = np.array([t["kind"] == "shuffle" for t in TR])
_rows = np.array([t["i"] for t in TR])
Y = np.array([CLASSES.index(t["label"]) for t in TR])
XTR = {w: np.where(_shuf[:, None], _X(FS_N, FB_N, w)[_rows], _X(FS, FB, w)[_rows])
       for w in ("summary", "full")}
print(f"training set: {len(Y)} maps ({int(_shuf.sum())} shuffle), "
      f"{sum(r['label'] == 'mixed/unsure' for r in FIRST)} mixed/unsure left out"
      f"{'' if MERGE_2TILES else ', ' + str(sum(r['label'] == '2 tiles' for r in FIRST)) + ' 2 tiles left out'}")
for k, c in enumerate(CLASSES):
    print(f"  {c:>10}: {int((Y == k).sum()):3d}  (CP {int(((Y == k) & (np.array([t['region'] for t in TR]) == 'CP')).sum())},"
          f" CB {int(((Y == k) & (np.array([t['region'] for t in TR]) == 'CB')).sum())})")


# %% ── 3. cross-validation: how well are your labels reproduced? ─────────────
# 3. cell --> repeated stratified 5-fold, for both feature sets; 23 and 25 scored
# against the same labels (real maps; your "2 tiles" counted as their "tiling").
def _model():
    return Pipeline([("sc", StandardScaler()),
                     ("lr", LogisticRegressionCV(Cs=CS, cv=StratifiedKFold(4, shuffle=True, random_state=SEED),
                                                 class_weight="balanced", scoring="balanced_accuracy",
                                                 max_iter=5000))])


def _cv_proba(X, y):
    P = np.zeros((N_REPS, len(y), K))
    for r in range(N_REPS):
        cv = StratifiedKFold(N_FOLDS, shuffle=True, random_state=SEED + r)
        P[r] = cross_val_predict(_model(), X, y, cv=cv, method="predict_proba")
    return P


def _confusion(y, pred, n_pred):
    C = np.zeros((K, n_pred))
    for a, b in zip(y, pred): C[a, b] += 1
    return C


CV = {}
for w in ("summary", "full"):
    t0 = time.time()
    P = _cv_proba(XTR[w], Y)
    bacc = [balanced_accuracy_score(Y, p.argmax(1)) for p in P]       # per repeat
    CV[w] = dict(P=P.mean(0), bacc=np.array(bacc))
    print(f"{w:>8}: balanced accuracy {np.mean(bacc):.2f} (range over repeats "
          f"{np.min(bacc):.2f}-{np.max(bacc):.2f}; chance {1 / K:.2f})  [{time.time() - t0:.0f}s]")

# 23 and 25 on the real labelled maps
_real = ~_shuf
_y_old = np.array([CLASSES.index("tiling") if CLASSES[k] == "2 tiles" else k for k in Y])
OLD_COLS = list(CLASSES) + ["unsure"]
OLD = {}
for nm, key in (("23", "c23"), ("25", "c25")):
    pred = np.array([OLD_COLS.index(t[key]) if t[key] in OLD_COLS else len(OLD_COLS) - 1 for t in TR])
    OLD[nm] = dict(pred=pred[_real], acc=float(np.mean(pred[_real] == _y_old[_real])))
    print(f"      {nm}: agreement with your labels on the real maps {OLD[nm]['acc']:.2f}")

fig, axes = plt.subplots(1, 4, figsize=(22, 5.2), gridspec_kw=dict(width_ratios=[1, 1, 1.15, 1.15]))
panels = [(f"classifier, {w} features\nbalanced acc. {CV[w]['bacc'].mean():.2f}",
           _confusion(Y, CV[w]["P"].argmax(1), K), list(CLASSES)) for w in ("summary", "full")]
panels += [(f"script {nm} (real maps)\nagreement {OLD[nm]['acc']:.2f}",
            _confusion(_y_old[_real], OLD[nm]["pred"], len(OLD_COLS)), OLD_COLS) for nm in ("23", "25")]
for ax, (title, C, cols) in zip(axes, panels):
    R = C / np.maximum(C.sum(1, keepdims=True), 1)
    ax.imshow(R, cmap="Blues", vmin=0, vmax=1)
    for a in range(C.shape[0]):
        for b in range(C.shape[1]):
            if C[a, b]: ax.text(b, a, f"{int(C[a, b])}", ha="center", va="center", fontsize=8,
                                color="w" if R[a, b] > 0.6 else "k")
        if C[a].sum() == 0: continue
    for k in range(min(K, len(cols))): ax.add_patch(plt.Rectangle((k - .5, k - .5), 1, 1, fill=False,
                                                               ec="tab:red", lw=1.2))
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(K)); ax.set_yticklabels(CLASSES, fontsize=8)
    ax.set_xlabel("called"); ax.set_title(title, fontsize=9)
axes[0].set_ylabel("your label")
fig.suptitle("STEP 1 / cross-validated confusion (colour = share of the row, numbers = maps). "
             "The classifier never saw a map it is scored on.", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.92)); _save(fig, "step1_confusion")

PCV = CV[FEATURE_SET]["P"]; PRED_CV = PCV.argmax(1)
print(f"\nfinal feature set '{FEATURE_SET}', per class (cross-validated):")
for k, c in enumerate(CLASSES):
    rec = np.mean(PRED_CV[Y == k] == k) if (Y == k).any() else np.nan
    prec = np.mean(Y[PRED_CV == k] == k) if (PRED_CV == k).any() else np.nan
    print(f"  {c:>10}: recall {rec:.2f}  precision {prec:.2f}  (n = {int((Y == k).sum())})")


# %% ── 3b. hotspots of the tiling maps (1 or 2); confidence ────────────────────
# 3b. cell --> "2 tiles" is not a class: a tiling map has 1 or 2 hotspots, and 2 is
# given only when there are clearly two. A hotspot = connected patch of the smoothed
# z-map above max(Z_HOT, HOT_REL x the map's peak), >= MIN_HOT_CELLS cells.
# Ranking (HOT_RANK): "shape" = size x strength x roundness, where size x strength =
# the sum of (z - level) over the patch, and roundness = (minor SD + ROUND_PAD) /
# (major SD + ROUND_PAD) of the patch's cells (in cells, weighted by z - level):
# ~1 for a round blob, low for a 1-cell-wide line. Lines still count, they just rank
# lower. "peak" = highest peak z (the old ranking).
# The top-ranked one always counts; the next one counts as a 2nd only if ALL of
#   peak     its peak z >= SECOND_PEAK x the 1st one's      (both strong)
#   size     its cells  >= SECOND_SIZE x the 1st one's      (not a speck of noise)
#   distance centres >= MIN_SEP s apart (neural x kinematic plane)
#   valley   between the two peaks the map dips to <= DIP_MAX x the 2nd peak
#            (the saddle: the highest level at which the two are still connected)
# The tables check this against YOUR tiling and 2 tiles labels (your 2 tiles labels
# are not a class any more, but they are what "2 hotspots" should find).
# Confidence: a call whose probability is < CONF_MIN becomes "low confidence" and is
# left out of the class percentages, the examples and the times.
# (the settings HOT_REL ... DIP_MAX and CONF_MIN are in cell 0)


def _zs(i, shuffle=False):
    """smoothed z-map of map i (or of its shuffle map)"""
    M = C23["sh"][_row23[uids[i]]] if shuffle else C12["mi"][i]
    with np.errstate(all="ignore"):
        z = (M - C12["null"][i]) / np.where(C12["null_std"][i] > 0, C12["null_std"][i], np.nan)
    return _smooth(np.clip(z, -Z_CLIP, Z_CLIP))


def _hotspots_all(zs):
    """-> [(neural t, kinematic t, peak z, n cells, peak row, peak col, score)], best first
    (score: size x strength x roundness, or the peak z; see HOT_RANK)"""
    assert HOT_RANK in ("shape", "peak"), HOT_RANK
    if not np.isfinite(zs).any(): return []
    thr = max(Z_HOT, HOT_REL * np.nanmax(zs))
    w = np.where(np.isfinite(zs), zs - thr, 0.0)
    lab, n = nd_label(w > 0, structure=np.ones((3, 3)))
    out = []
    for c in range(1, n + 1):
        cl = lab == c
        if cl.sum() < MIN_HOT_CELLS: continue
        wc = w[cl]
        pr, pc = np.unravel_index(np.argmax(np.where(cl, zs, -np.inf)), zs.shape)
        if HOT_RANK == "shape":
            rr, cc = np.nonzero(cl)                          # same order as w[cl]
            ev = np.linalg.eigvalsh(np.atleast_2d(np.cov(np.vstack([rr, cc]), aweights=wc, bias=True)))
            sd = np.sqrt(np.clip(ev, 0, None))
            score = float(wc.sum() * (sd[0] + ROUND_PAD) / (sd[-1] + ROUND_PAD))
        else:
            score = float(zs[pr, pc])
        out.append((float((wc * Nc[cl]).sum() / wc.sum()), float((wc * Kc[cl]).sum() / wc.sum()),
                    float(zs[pr, pc]), int(cl.sum()), int(pr), int(pc), score))
    return sorted(out, key=lambda b: -b[6])


def _saddle(zs, a, b):
    """highest level at which cells a and b are still connected above it (bisection)"""
    Zf = np.where(np.isfinite(zs), zs, -np.inf)
    lo, hi = float(np.nanmin(zs)), float(min(Zf[a], Zf[b]))
    for _ in range(20):
        t = (lo + hi) / 2
        lab, _ = nd_label(Zf > t, structure=np.ones((3, 3)))
        if lab[a] and lab[a] == lab[b]: lo = t
        else: hi = t
    return lo


def _second_ok(zs, h1, h2, peak=None, size=None, sep=None, dip=None):
    peak = SECOND_PEAK if peak is None else peak; size = SECOND_SIZE if size is None else size
    sep = MIN_SEP if sep is None else sep; dip = DIP_MAX if dip is None else dip
    if h2[2] < peak * h1[2] or h2[3] < size * h1[3]: return False
    if np.hypot(h2[0] - h1[0], h2[1] - h1[1]) < sep: return False
    return dip >= 1 or _saddle(zs, h1[4:6], h2[4:6]) <= dip * h2[2]


def _hotspots(zs, **kw):
    """the hotspots that count: the top-ranked, plus the next one if it is clearly a 2nd"""
    hs = _hotspots_all(zs)
    if len(hs) < 2: return hs
    return hs[:2] if _second_ok(zs, hs[0], hs[1], **kw) else hs[:1]


def _conf(P, pred):
    """probability of each call"""
    return P[np.arange(len(pred)), pred]


t0 = time.time()
HOT, HOT_N = [None] * len(uids), [None] * len(uids)
for i in np.flatnonzero(OK):
    HOT[i], HOT_N[i] = _hotspots(_zs(i)), _hotspots(_zs(i, True))
NHOT = np.array([len(h) if h is not None else 0 for h in HOT])
NHOT_N = np.array([len(h) if h is not None else 0 for h in HOT_N])
print(f"hotspots of {int(OK.sum())} maps (+ shuffle maps) in {time.time() - t0:.0f}s")

# check against your labels: your tiling maps should get 1, your 2 tiles maps 2
_LAB_T = {c: [(ROW[r["uid"]], r["kind"] == "shuffle") for r in FIRST
              if r["label"] == c and r["uid"] in ROW and OK[ROW[r["uid"]]]] for c in ("tiling", "2 tiles")}


def _share2(c, **kw):
    return 100 * np.mean([len(_hotspots(_zs(i, sh), **kw)) >= 2 for i, sh in _LAB_T[c]]) \
        if _LAB_T[c] else np.nan


print(f"\nyour labels -> 2 hotspots (ranked by {HOT_RANK}; SECOND_PEAK {SECOND_PEAK},SECOND_SIZE {SECOND_SIZE}, "
      f"MIN_SEP {MIN_SEP} s, DIP_MAX {DIP_MAX}):")
print(f"  your tiling  ({len(_LAB_T['tiling'])}): {_share2('tiling'):3.0f}% given 2   (want low)")
print(f"  your 2 tiles ({len(_LAB_T['2 tiles'])}): {_share2('2 tiles'):3.0f}% given 2   (want high)")
print("one setting at a time (the others as above): tiling % | 2 tiles %")
for nm, vals in (("peak", (0.0, 0.6, 0.7, 0.8)), ("size", (0.0, 0.3, 0.5, 0.7)),
                 ("sep", (0.0, 0.1, 0.15, 0.2)), ("dip", (1.0, 0.6, 0.5, 0.4))):
    print(f"  {nm:>5}: " + "   ".join(f"{v:.2f}: {_share2('tiling', **{nm: v}):3.0f}|{_share2('2 tiles', **{nm: v}):3.0f}"
                                     for v in vals))

# confidence threshold, on the cross-validated calls: how many labelled maps are kept,
# and how often the kept ones agree with you
CONF_CV = _conf(PCV, PRED_CV)
print("\nconfidence threshold (cross-validated, labelled maps):")
for cm_ in (0.0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
    m = CONF_CV >= cm_
    if m.sum() < 2: continue
    per = "  ".join(f"{CLASSES[k]} {int((m & (Y == k)).sum())}/{int((Y == k).sum())}" for k in range(K))
    print(f"  CONF_MIN {cm_:.1f}: {100 * m.mean():3.0f}% kept, agree with you "
          f"{100 * np.mean(PRED_CV[m] == Y[m]):3.0f}%, balanced acc. "
          f"{balanced_accuracy_score(Y[m], PRED_CV[m]):.2f}{'   <- current' if np.isclose(cm_, CONF_MIN) else ''}"
          f"\n      kept per class (your label): {per}")
CONF_TR = CONF_CV >= CONF_MIN


# %% ── 4. the final model and what drives each class ─────────────────────────
# 4. cell --> fit on all labelled maps; standardised coefficients (summary
# features, and the 6 x 6 layout as an image per class).
MODEL = _model().fit(XTR[FEATURE_SET], Y)
COEF = MODEL.named_steps["lr"].coef_                       # (K, n_features)
print(f"chosen regularisation C = {MODEL.named_steps['lr'].C_[0]:.3g}")
ns = len(SUMMARY_NAMES)
fig = plt.figure(figsize=(16, 6 if FEATURE_SET == "summary" else 9))
gs = fig.add_gridspec(2 if FEATURE_SET == "full" else 1, K, height_ratios=[1.4, 1] if FEATURE_SET == "full" else [1])
ax = fig.add_subplot(gs[0, :])
lim = np.abs(COEF[:, :ns]).max()
ax.imshow(COEF[:, :ns], cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
ax.set_xticks(range(ns)); ax.set_xticklabels(SUMMARY_NAMES, rotation=40, ha="right", fontsize=8)
ax.set_yticks(range(K)); ax.set_yticklabels(CLASSES)
for a in range(K):
    for b in range(ns): ax.text(b, a, f"{COEF[a, b]:+.1f}", ha="center", va="center", fontsize=7)
ax.set_title("coefficients on standardised features (red: more of this feature -> more likely this class)",
             fontsize=9)
if FEATURE_SET == "full":
    lim = np.abs(COEF[:, ns:]).max()
    for k, c in enumerate(CLASSES):
        a2 = fig.add_subplot(gs[1, k])
        a2.imshow(COEF[k, ns:].reshape(N_BLOCK, N_BLOCK), origin="lower", cmap="RdBu_r",
                  vmin=-lim, vmax=lim, extent=extent, aspect="auto")
        d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
        a2.plot(d, d, "k--", lw=0.8)
        a2.set_title(f"{c}: layout weights", fontsize=9)
        a2.set_xlabel("neural (s)", fontsize=8)
        if k == 0: a2.set_ylabel("kinematic (s)", fontsize=8)
        a2.tick_params(labelsize=7)
plt.tight_layout(); _save(fig, "step2_coefficients")


# %% ── 5. apply to every map ─────────────────────────────────────────────────
# 5. cell --> class of every real map and every shuffle map; % per region, raw and
# corrected for the classifier's errors (cross-validated confusion matrix).
XA, XA_N = _X(FS, FB, FEATURE_SET)[OK], _X(FS_N, FB_N, FEATURE_SET)[OK]
PA, PA_N = MODEL.predict_proba(XA), MODEL.predict_proba(XA_N)
_pr, _pr_N = PA.argmax(1), PA_N.argmax(1)
CONF = np.full(len(uids), np.nan); CONF[OK] = _conf(PA, _pr)      # probability of the call
CONF_N = np.full(len(uids), np.nan); CONF_N[OK] = _conf(PA_N, _pr_N)
CLS = np.full(len(uids), "not classified", dtype=object)
CLS[OK] = np.where(CONF[OK] >= CONF_MIN, np.array(CLASSES)[_pr], LOW)
PROB = np.full((len(uids), K), np.nan); PROB[OK] = PA        # the model's own probabilities
CLS_MODEL = np.full(len(uids), "not classified", dtype=object)  # before the threshold
CLS_MODEL[OK] = np.array(CLASSES)[PA.argmax(1)]
CLS_N = np.full(len(uids), "not classified", dtype=object)
CLS_N[OK] = np.where(CONF_N[OK] >= CONF_MIN, np.array(CLASSES)[_pr_N], LOW)
print(f"low confidence (p < {CONF_MIN}): {int((CLS[OK] == LOW).sum())} of {int(OK.sum())} maps "
      f"(shuffle maps: {int((CLS_N[OK] == LOW).sum())})")

# P(called j | truly i) from the cross-validation, on the REAL labelled maps that pass
# the confidence threshold (the percentages below are of confident maps too)
CM = _confusion(Y[_real & CONF_TR], PRED_CV[_real & CONF_TR], K)
CM = CM / np.maximum(CM.sum(1, keepdims=True), 1)


def _corrected(p_obs):
    """the class mix p with CM^T p closest to the observed mix (p >= 0, sum 1)"""
    A = np.vstack([CM.T, 100 * np.ones(K)]); b = np.r_[p_obs, 100.0]
    p, _ = nnls(A, b)
    return p / p.sum()


PCT = {}
for reg in ("CP", "CB"):
    m_all = OK & (regions == reg)
    m = m_all & (CLS != LOW); mn = m_all & (CLS_N != LOW)     # confident maps only
    raw = np.array([np.mean(CLS[m] == c) for c in CLASSES])
    PCT[reg] = dict(raw=raw, corrected=_corrected(raw),
                    shuffle=np.array([np.mean(CLS_N[mn] == c) for c in CLASSES]), n=int(m.sum()),
                    n_low=int((m_all & ~m).sum()), n_all=int(m_all.sum()))
    print(f"{reg} ({int(m.sum())} confident maps of {int(m_all.sum())}; "
          f"{100 * np.mean(CLS[m_all] == LOW):.0f}% low confidence):")
    for k, c in enumerate(CLASSES):
        print(f"  {c:>10}: {100 * raw[k]:5.1f}%  corrected {100 * PCT[reg]['corrected'][k]:5.1f}%"
              f"   | shuffle maps {100 * PCT[reg]['shuffle'][k]:5.1f}%")
    _shape = (CLS_N[m_all] != LOW) & (CLS_N[m_all] != "nothing")
    print(f"  shuffle maps confidently called a shape (not 'nothing'): {100 * np.mean(_shape):.1f}% "
          f"of all {int(m_all.sum())}")
    mt = m & (CLS == "tiling")
    PCT[reg]["tiling_2hot"] = int((mt & (NHOT >= 2)).sum())
    print(f"  tiling maps: {int((mt & (NHOT < 2)).sum())} with 1 hotspot, "
          f"{PCT[reg]['tiling_2hot']} with 2 ({100 * np.mean(NHOT[mt] >= 2):.0f}%)")
_yn =np.array([t["label"] for t in TR])[_shuf]
print(f"your own labels on the {int(_shuf.sum())} labelled shuffle maps: "
      f"{100 * np.mean(_yn != 'nothing'):.0f}% called a shape")

fig, axes = plt.subplots(1, 2, figsize=(14, 4.5), sharey=True)
x = np.arange(K); wd = 0.27
for ax, reg in zip(axes, ("CP", "CB")):
    ax.bar(x - wd, 100 * PCT[reg]["raw"], wd, color=CLR[reg], label="classifier")
    ax.bar(x, 100 * PCT[reg]["corrected"], wd, color=CLR[reg], alpha=0.5, hatch="//",
           label="corrected for its errors")
    ax.bar(x + wd, 100 * PCT[reg]["shuffle"], wd, color="0.65", label="shuffle maps")
    ax.set_xticks(x); ax.set_xticklabels(CLASSES, rotation=30, ha="right")
    ax.set_title(f"{reg}: {PCT[reg]['n']} confident maps (+ {PCT[reg]['n_low']} with p < {CONF_MIN})",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("% of maps"); axes[0].legend(fontsize=8)
fig.suptitle("STEP 3 / class of every map", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.93)); _save(fig, "step3_percentages")


# %% ── 6. examples ───────────────────────────────────────────────────────────
# 6. cell --> per region, the most confident maps of each class (smoothed z), and
# the labelled maps where the cross-validated classifier disagrees with you.

def _show(ax, Zm, title, clr="k"):
    ok = np.isfinite(Zm).any()
    ax.imshow(Zm, origin="lower", aspect="auto", extent=extent, cmap="viridis",
              vmin=np.nanpercentile(Zm, 1) if ok else 0, vmax=np.nanpercentile(Zm, 99) if ok else 1)
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(d, d, "w--", lw=0.7)
    ax.set_title(title, fontsize=6.5, color=clr); ax.set_xticks([]); ax.set_yticks([])


EXAMPLES = []      # (group, uid, kind, comment): every map shown in an example figure,
                   # written as a list for 27's gallery (python 27_label_gui.py --relabel <file>)
for reg in ("CP", "CB"):
    fig, axes = plt.subplots(K, N_EX, figsize=(2.0 * N_EX, 1.8 * K), squeeze=False)
    for k, c in enumerate(CLASSES):
        idx = np.flatnonzero(OK & (regions == reg) & (CLS == c))
        idx = idx[np.argsort(-CONF[idx])][:N_EX]
        EXAMPLES += [(f"step 4 / {reg} called {c} (most confident)", uids[i], "real",
                      f"model: {c} ({CONF[i]:.2f})") for i in idx]
        for e in range(N_EX):
            ax = axes[k, e]
            if e < len(idx):
                i = idx[e]; _, mo, se, un = uids[i].split("|")
                _show(ax, _zs(i), f"{mo} {se} {un}  p={CONF[i]:.2f}")
            else: ax.axis("off")
        axes[k, 0].set_ylabel(c, fontsize=9, color=CLASS_CLR[c]); axes[k, 0].axis("on")
        axes[k, 0].set_xticks([]); axes[k, 0].set_yticks([])
    fig.suptitle(f"STEP 4 / {reg}: the most confident maps of each class (smoothed z)", fontsize=10)
    plt.tight_layout(rect=(0, 0, 1, 0.96)); _save(fig, f"step4_examples_{reg}")

wrong = np.flatnonzero(PRED_CV != Y)
wrong = wrong[np.argsort(-PCV[wrong, PRED_CV[wrong]])][:32]
nc = 8; nr = max(1, int(np.ceil(len(wrong) / nc)))
fig, axes = plt.subplots(nr, nc, figsize=(2.1 * nc, 2.0 * nr), squeeze=False)
for e, ax in enumerate(axes.ravel()):
    if e >= len(wrong): ax.axis("off"); continue
    t = TR[wrong[e]]; _, mo, se, un = t["uid"].split("|")
    EXAMPLES.append(("step 5 / labelled maps the model disagrees with", t["uid"], t["kind"],
                     f"you: {t['label']} | model: {CLASSES[PRED_CV[wrong[e]]]}"))
    _show(ax, _zs(t["i"], t["kind"] == "shuffle"),
          f"{t['region']} {mo} {se} {un}{' (shuffle)' if t['kind'] == 'shuffle' else ''}\n"
          f"you: {t['label']} | model: {CLASSES[PRED_CV[wrong[e]]]} "
          f"({PCV[wrong[e], PRED_CV[wrong[e]]]:.2f})")
fig.suptitle(f"STEP 5 / labelled maps the cross-validated classifier gets 'wrong' "
             f"({len(np.flatnonzero(PRED_CV != Y))} of {len(Y)}; most confident first). "
             "Some may be labels worth a second look.", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.95)); _save(fig, "step5_disagreements")


# %% ── 7. your agreement with yourself (repeat pass) ─────────────────────────
# 7. cell --> runs once the repeat pass of 27 is labelled: the ceiling for the
# classifier's accuracy.
if REPEAT_ROWS:
    first = {r["key"]: r["label"] for r in FIRST}
    pairs = [(first[r["key"]], r["label"]) for r in REPEAT_ROWS if r["key"] in first]
    a, b = map(np.array, zip(*pairs))
    print(f"repeat pass: {len(pairs)} maps labelled twice; same label {100 * np.mean(a == b):.0f}%, "
          f"Cohen's kappa {cohen_kappa_score(a, b):.2f}")
    allc = list(dict.fromkeys(list(CLASSES) + ["2 tiles", "mixed/unsure"]))
    C = np.zeros((len(allc), len(allc)))
    for p, q in pairs: C[allc.index(_lab(p) if p in allc else p), allc.index(q)] += 1
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.imshow(C / np.maximum(C.sum(1, keepdims=True), 1), cmap="Blues", vmin=0, vmax=1)
    for i in range(len(allc)):
        for j in range(len(allc)):
            if C[i, j]: ax.text(j, i, int(C[i, j]), ha="center", va="center", fontsize=8)
    ax.set_xticks(range(len(allc))); ax.set_xticklabels(allc, rotation=45, ha="right")
    ax.set_yticks(range(len(allc))); ax.set_yticklabels(allc)
    ax.set(xlabel="repeat pass", ylabel="first pass",
           title=f"your agreement with yourself: {100 * np.mean(a == b):.0f}%")
    plt.tight_layout(); _save(fig, "step6_self_agreement")
else:
    print("no repeat-pass labels yet (python 27_label_gui.py --repeat); this cell runs once there are")


# %% ── 8. where the centres are (check on example maps) ──────────────────────
# 8. cell --> one time per map, read from its smoothed z-map (positive part):
#   column    neural time of the peak of the column profile (mean over kinematic windows)
#   diagonal  lag (neural - kinematic) of the peak of the lag-diagonal profile
#   tiling    centre of its hotspot (cell 3b), in neural / kinematic time; tiling maps
#             with 2 hotspots (cell 3b) get both centres
# A peak's time = weighted centre of the above-half-maximum run around it.
# The figure marks them on random maps of each class, to check they sit right.
# (the settings P_MIN, N_CHECK, MIN_LAG_CELLS are in cell 0)
TIME_CLASSES = ("column", "diagonal", T1, T2)


def _peak_centre(prof, x):
    p = np.where(np.isfinite(prof), prof, -np.inf)
    i0 = int(np.argmax(p)); pk = p[i0]
    if not np.isfinite(pk) or pk <= 0: return np.nan
    lo, hi = i0, i0
    while lo > 0 and p[lo - 1] >= pk / 2: lo -= 1
    while hi < len(p) - 1 and p[hi + 1] >= pk / 2: hi += 1
    w = p[lo:hi + 1]
    return float((w * x[lo:hi + 1]).sum() / w.sum())


def _centres(zs):
    ok = np.isfinite(zs); zp = np.where(ok, np.clip(zs, 0, None), np.nan)
    with np.errstate(all="ignore"):
        col_t = _peak_centre(np.nanmean(zp, axis=0), win_centers)
    g = _LAGG[ok]; g0 = g.min()
    s, n = np.bincount(g - g0, zp[ok]), np.bincount(g - g0)
    with np.errstate(all="ignore"):
        prof = np.where(n >= MIN_LAG_CELLS, s / n, np.nan)
    diag_lag = _peak_centre(prof, (np.arange(len(prof)) + g0) * LAG_BIN)
    return col_t, diag_lag


POS = dict(col_t=np.full(len(uids), np.nan), diag_lag=np.full(len(uids), np.nan),
           blob1=np.full((len(uids), 2), np.nan), blob2=np.full((len(uids), 2), np.nan))
for i in np.flatnonzero(OK):
    POS["col_t"][i], POS["diag_lag"][i] = _centres(_zs(i))
    if HOT[i]: POS["blob1"][i] = HOT[i][0][:2]               # hotspots of cell 3b, top-ranked first
    if HOT[i] and len(HOT[i]) > 1: POS["blob2"][i] = HOT[i][1][:2]
PMAX = np.where(OK, np.nan_to_num(CONF, nan=0.0), 0.0)      # probability of the call (cell 5)
_base = OK & (PMAX >= P_MIN)
SEL = {"column": _base & (CLS == "column"), "diagonal": _base & (CLS == "diagonal"),
       T1: _base & (CLS == "tiling") & (NHOT < 2), T2: _base & (CLS == "tiling") & (NHOT >= 2)}
for c in TIME_CLASSES:
    print(f"{c:>18}: {int(SEL[c].sum())} maps (p >= {max(P_MIN, CONF_MIN)})")


def _mark(ax, i, c):
    kw = dict(color="r", lw=1.4)
    if c == "column" and np.isfinite(POS["col_t"][i]):
        ax.axvline(POS["col_t"][i], **kw)
    elif c == "diagonal" and np.isfinite(POS["diag_lag"][i]):
        L = POS["diag_lag"][i]; kk = np.array(extent[2:])
        ax.plot(kk + L, kk, **kw); ax.set_xlim(extent[:2]); ax.set_ylim(extent[2:])
    else:
        for b, mk in ((POS["blob1"][i], "+"), (POS["blob2"][i] if c == T2 else (np.nan,) * 2, "x")):
            if np.isfinite(b[0]): ax.plot(b[0], b[1], mk, color="r", ms=11, mew=2.2)


_rc = np.random.default_rng(3)
fig, axes = plt.subplots(len(TIME_CLASSES), N_CHECK, figsize=(2.1 * N_CHECK, 2.0 * len(TIME_CLASSES)),
                         squeeze=False)
for r, c in enumerate(TIME_CLASSES):
    idx = np.flatnonzero(SEL[c])
    idx = _rc.choice(idx, min(N_CHECK, len(idx)), replace=False) if len(idx) else idx
    EXAMPLES += [(f"step 7 / centres check: called {c}", uids[i], "real",
                  f"model: {c} ({PMAX[i]:.2f})") for i in idx]
    for e in range(N_CHECK):
        ax = axes[r, e]
        if e >= len(idx): ax.axis("off"); continue
        i = idx[e]; reg, mo, se, un = uids[i].split("|")
        val = (f"t={POS['col_t'][i]:+.3f}s" if c == "column" else
               f"lag={POS['diag_lag'][i]:+.3f}s" if c == "diagonal" else
               f"n={POS['blob1'][i][0]:+.2f} k={POS['blob1'][i][1]:+.2f}")
        _show(ax, _zs(i), f"{reg} {mo} {se} {un}\n{val}")
        _mark(ax, i, c)
    axes[r, 0].set_ylabel(c, fontsize=9, color=CLASS_CLR[c])
fig.suptitle("STEP 7 / where the centres are, on random maps of each class (red: column time, "
             "diagonal lag, blob centre +, 2nd blob x; smoothed z)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.95)); _save(fig, "step7_centres_check")


# %% ── 9. when: the times of the columns, diagonals and blobs ─────────────────
# 9. cell --> per region (rows): blob centres on the map plane (2 tiles: the two
# blobs joined), then the distributions of neural time (columns, blobs),
# kinematic time (blobs) and lag (diagonals, blobs). Dashed: medians.
from scipy.stats import mannwhitneyu


def _vals(reg, c, what):
    m = SEL[c] & (regions == reg)
    b = POS["blob1"][m] if c != T2 else np.vstack([POS["blob1"][m], POS["blob2"][m]])
    v = {"neural": POS["col_t"][m] if c == "column" else b[:, 0],
         "kinematic": b[:, 1],
         "lag": POS["diag_lag"][m] if c == "diagonal" else b[:, 0] - b[:, 1]}[what]
    return v[np.isfinite(v)]


PANELS = (("neural", "neural time (s)", ("column", T1, T2)),
          ("kinematic", "kinematic time (s)", (T1, T2)),
          ("lag", "lag = neural - kinematic (s)", ("diagonal", T1, T2)))
TIMES = {}
fig = plt.figure(figsize=(20, 8.5)); gs_ = fig.add_gridspec(2, 4)
for r, reg in enumerate(("CP", "CB")):
    ax = fig.add_subplot(gs_[r, 0])
    m1 = SEL[T1] & (regions == reg); m2 = SEL[T2] & (regions == reg)
    ax.plot(POS["blob1"][m1, 0], POS["blob1"][m1, 1], "o", ms=3, color=CLASS_CLR[T1], alpha=0.6,
            label=f"{T1} (n={int(m1.sum())})")
    ax.plot(POS["blob1"][m2, 0], POS["blob1"][m2, 1], "o", ms=3, color=CLASS_CLR[T2], alpha=0.7,
            label=f"{T2} (n={int(m2.sum())} maps, 2 centres each)")
    ax.plot(POS["blob2"][m2, 0], POS["blob2"][m2, 1], "o", ms=3, color=CLASS_CLR[T2], alpha=0.7)
    dd = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(dd, dd, color="k", lw=0.8, ls="--")
    ax.set(xlim=extent[:2], ylim=extent[2:], xlabel="neural time (s)", ylabel="kinematic time (s)",
           title=f"{reg}: blob centres")
    ax.legend(fontsize=7, loc="lower right")
    for c_, (what, xlab, classes) in enumerate(PANELS, start=1):
        ax = fig.add_subplot(gs_[r, c_])
        allv = np.concatenate([_vals(g, c, what) for g in ("CP", "CB") for c in classes] + [[0.0]])
        bins = np.linspace(np.percentile(allv, 0.5), np.percentile(allv, 99.5), 30)
        for c in classes:
            v = _vals(reg, c, what)
            TIMES[(reg, c, what)] = v
            if v.size < 2: continue
            ax.hist(np.clip(v, bins[0], bins[-1]), bins, histtype="step", lw=2, density=True,
                    color=CLASS_CLR[c], label=f"{c} (n={v.size}, median {np.median(v):+.3f})")
            ax.axvline(np.median(v), color=CLASS_CLR[c], ls="--", lw=1)
        ax.axvline(0, color="k", lw=0.6, ls=":")
        ax.set(xlabel=xlab, ylabel="density"); ax.legend(fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("STEP 8 / when: column times, diagonal lags and hotspot centres (classifier's classes, "
             f"p >= {max(P_MIN, CONF_MIN)}; tiling with 2 hotspots: both centres counted)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.95)); _save(fig, "step8_times")

print("median [25th, 75th pct] (s); CP vs CB: Mann-Whitney p")
for what, _, classes in PANELS:
    for c in classes:
        line = f"  {what:>9} {c:>18}:"
        for reg in ("CP", "CB"):
            v = TIMES.get((reg, c, what), np.array([]))
            line += (f"  {reg} {np.median(v):+.3f} [{np.percentile(v, 25):+.3f}, {np.percentile(v, 75):+.3f}]"
                     f" n={v.size}" if v.size else f"  {reg} --")
        a, b = TIMES.get(("CP", c, what), []), TIMES.get(("CB", c, what), [])
        if len(a) > 1 and len(b) > 1: line += f"  p={mannwhitneyu(a, b).pvalue:.2g}"
        print(line)


# %% ── 9b. the same times as counts, plus all classes together ─────────────────
# 9b. cell --> cell 9's histograms as COUNTS (not density), so regions and classes can
# be compared by size; the same y-axis for CP and CB in each column. Black: all the
# classes of that panel together. Tiling with 2 hotspots counts both centres, so the
# counts are centres, not maps.
fig, axes = plt.subplots(2, len(PANELS), figsize=(5.5 * len(PANELS), 8.5), sharey="col", squeeze=False)
for c_, (what, xlab, classes) in enumerate(PANELS):
    allv = np.concatenate([_vals(g, c, what) for g in ("CP", "CB") for c in classes] + [[0.0]])
    bins = np.linspace(np.percentile(allv, 0.5), np.percentile(allv, 99.5), 30)
    for r, reg in enumerate(("CP", "CB")):
        ax = axes[r, c_]
        vs = [_vals(reg, c, what) for c in classes]
        for c, v in zip(classes, vs):
            if v.size < 1: continue
            ax.hist(np.clip(v, bins[0], bins[-1]), bins, histtype="step", lw=2,
                    color=CLASS_CLR[c], label=f"{c} (n={v.size}, median {np.median(v):+.3f})")
            ax.axvline(np.median(v), color=CLASS_CLR[c], ls="--", lw=1)
        va = np.concatenate(vs)
        if va.size:
            ax.hist(np.clip(va, bins[0], bins[-1]), bins, histtype="step", lw=2.2, color="k",
                    label=f"all (n={va.size}, median {np.median(va):+.3f})")
            ax.axvline(np.median(va), color="k", ls="--", lw=1)
        ax.axvline(0, color="k", lw=0.6, ls=":")
        ax.set(xlabel=xlab, ylabel="count", title=f"{reg}: {what}")
        ax.legend(fontsize=7); ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("STEP 8b / when, as counts: column times, diagonal lags and hotspot centres (classifier's "
             f"classes, p >= {max(P_MIN, CONF_MIN)}; black: all classes of the panel; tiling with 2 "
             "hotspots: both centres counted)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.95)); _save(fig, "step8b_times_counts")


# %% ── 9c. hotspot centres of every class on the map plane ───────────────────
# 9c. cell --> the hotspots of cell 3b (centre of mass of each patch above the hotspot
# level; the top-ranked, plus a 2nd one only if it passes the same rules as for
# tiling), now for the column and diagonal maps too. Same maps as cell 9 (p >= the
# confidence cuts). Filled dot: 1st hotspot, open dot: 2nd. First panel: all classes
# together; then one panel per class. Dashed: zero lag.
# Column and diagonal dots: size and opacity = how much the hotspot stands out from the rest of
# its stripe. Stripe = the neural window through the hotspot's peak (column) or the lag
# diagonal through it (diagonal). Stand-out = 1 - (mean smoothed z of the stripe's
# cells outside the hotspot, floored at 0) / (hotspot peak z), in [0, 1]:
# 1 = a focal spot on a flat stripe (opaque), 0 = the stripe is as strong as the spot
# (small, faint). Size (STANDOUT_SIZE) and opacity (STANDOUT_ALPHA_MIN to 1) run over STANDOUT_RANGE (most values
# sit in ~0.5-0.75, so the full 0-1 range would make them all look alike).
# Tiling dots keep a fixed opacity.
HOT_PLANE_CLASSES = ("column", "diagonal", "tiling")
_selc = {"column": SEL["column"], "diagonal": SEL["diagonal"], "tiling": SEL[T1] | SEL[T2]}


def _standout(zs, h, c):
    """how much hotspot h stands out from the rest of its column / lag diagonal (0..1)"""
    pr, pc, pk = h[4], h[5], h[2]
    thr = max(Z_HOT, HOT_REL * np.nanmax(zs))
    lab, _ = nd_label(np.where(np.isfinite(zs), zs, -np.inf) > thr, structure=np.ones((3, 3)))
    hot = lab == lab[pr, pc]
    stripe = (_COLG == pc) if c == "column" else (_LAGG == _LAGG[pr, pc])
    rest = zs[stripe & ~hot & np.isfinite(zs)]
    if rest.size == 0 or pk <= 0: return 0.0              # the whole stripe is the hotspot
    return float(np.clip(1 - max(np.mean(rest), 0) / pk, 0, 1))


STANDOUT = np.full((len(uids), 2), np.nan)                 # 1st, 2nd hotspot
for c in ("column", "diagonal"):
    for i in np.flatnonzero(_selc[c]):
        z_ = _zs(i)
        for j, h in enumerate(HOT[i][:2]): STANDOUT[i, j] = _standout(z_, h, c)
for c in ("column", "diagonal"):
    for reg in ("CP", "CB"):
        v = STANDOUT[_selc[c] & (regions == reg), 0]; v = v[np.isfinite(v)]
        if v.size: print(f"stand-out of the 1st hotspot, {c:>8} {reg}: median {np.median(v):.2f} "
                         f"[{np.percentile(v, 25):.2f}, {np.percentile(v, 75):.2f}]  n={v.size}")


def _scaled(a):
    lo, hi = STANDOUT_RANGE
    return np.clip((np.nan_to_num(a, nan=0.0) - lo) / (hi - lo), 0, 1)


def _rgba(clr, a):
    from matplotlib.colors import to_rgb
    out = np.tile(np.r_[to_rgb(clr), 1.0], (len(a), 1))
    out[:, 3] = (STANDOUT_ALPHA_MIN + (1 - STANDOUT_ALPHA_MIN) * _scaled(a)
                 if STANDOUT_STYLE in ("opacity", "both") else 0.7)
    return out


def _size(a):
    if STANDOUT_STYLE not in ("size", "both"): return np.full(len(a), 14.0)
    return STANDOUT_SIZE[0] + (STANDOUT_SIZE[1] - STANDOUT_SIZE[0]) * _scaled(a)


assert STANDOUT_STYLE in ("opacity", "size", "both"), STANDOUT_STYLE


dd = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
fig, axes = plt.subplots(2, 1 + len(HOT_PLANE_CLASSES), figsize=(5 * (1 + len(HOT_PLANE_CLASSES)), 9.5),
                         squeeze=False)
for r, reg in enumerate(("CP", "CB")):
    for p_, show in enumerate([HOT_PLANE_CLASSES] + [(c,) for c in HOT_PLANE_CLASSES]):
        ax = axes[r, p_]
        for c in show:
            m = _selc[c] & (regions == reg)
            m2 = m & np.isfinite(POS["blob2"][:, 0])
            clr = CLASS_CLR[c]
            lab_ = f"{c} (n={int(m.sum())} maps, {int(m2.sum())} with a 2nd)"
            if c == "tiling":
                ax.plot(POS["blob1"][m, 0], POS["blob1"][m, 1], "o", ms=3.5, color=clr, alpha=0.6, label=lab_)
                ax.plot(POS["blob2"][m2, 0], POS["blob2"][m2, 1], "o", ms=4, mfc="none", mec=clr, alpha=0.8)
            else:                                            # opacity = stand-out
                ax.scatter(POS["blob1"][m, 0], POS["blob1"][m, 1], s=_size(STANDOUT[m, 0]),
                           c=_rgba(clr, STANDOUT[m, 0]), edgecolors="none")
                ax.scatter(POS["blob2"][m2, 0], POS["blob2"][m2, 1], s=_size(STANDOUT[m2, 1]),
                           facecolors="none", edgecolors=_rgba(clr, STANDOUT[m2, 1]))
                ax.plot([], [], "o", ms=3.5, color=clr, label=lab_ + f", {STANDOUT_STYLE} = stand-out")
        ax.plot(dd, dd, color="k", lw=0.8, ls="--")
        ax.axvline(0, color="k", lw=0.5, ls=":"); ax.axhline(0, color="k", lw=0.5, ls=":")
        ax.set(xlim=extent[:2], ylim=extent[2:], xlabel="neural time (s)", ylabel="kinematic time (s)",
               title=f"{reg}: {'all classes' if len(show) > 1 else show[0]}")
        ax.legend(fontsize=7, loc="lower right")
fig.suptitle("STEP 8c / hotspot centres of column, diagonal and tiling maps (classifier's classes, "
             f"p >= {max(P_MIN, CONF_MIN)}; filled: 1st hotspot, open: 2nd; column/diagonal {STANDOUT_STYLE}: "
             "how much the hotspot stands out from the rest of its stripe)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.95)); _save(fig, "step8c_hotspot_centres")


# %% ── 10. save ─────────────────────────────────────────────────────────────
if SAVE_VARS:
    p = config.var_path(f"{RUN_NAME}_{RUN_TAG}_results.pkl")
    with open(p, "wb") as fh:
        pickle.dump(dict(uids=uids, regions=regions, classified=OK, cls=CLS, proba=PROB, conf=CONF,
                         cls_shuffle=CLS_N, classes=CLASSES, pct=PCT, confusion_cv=CM,
                         positions=POS, p_min_times=P_MIN, n_hotspots=NHOT, n_hotspots_shuffle=NHOT_N,
                         cls_model_only=CLS_MODEL,
                         cv={w: dict(bacc=CV[w]["bacc"]) for w in CV}, old_acc={k: OLD[k]["acc"] for k in OLD},
                         feature_names=_names(FEATURE_SET), coef=COEF, model=MODEL,
                         train=dict(uids=[t["uid"] for t in TR], kind=[t["kind"] for t in TR], y=Y,
                                    proba_cv=PCV),
                         params=dict(FEATURE_SET=FEATURE_SET, MERGE_2TILES=MERGE_2TILES, N_BLOCK=N_BLOCK,
                                     HOT_REL=HOT_REL, MIN_HOT_CELLS=MIN_HOT_CELLS, HOT_RANK=HOT_RANK,
                                     ROUND_PAD=ROUND_PAD, SECOND_PEAK=SECOND_PEAK,
                                     SECOND_SIZE=SECOND_SIZE, MIN_SEP=MIN_SEP, DIP_MAX=DIP_MAX,
                                     CONF_MIN=CONF_MIN,
                                     SMOOTH_SIG=SMOOTH_SIG, Z_CLIP=Z_CLIP, Z_HOT=Z_HOT,
                                     MIN_CLUSTER=MIN_CLUSTER, TOP_FRAC=TOP_FRAC, LAG_BIN=LAG_BIN,
                                     MIN_FINITE=MIN_FINITE, N_FOLDS=N_FOLDS, N_REPS=N_REPS, CS=CS,
                                     SEED=SEED, SRC_SET=SRC_SET, SRC_LABELS=SRC_LABELS, SRC_12=SRC_12,
                                     SRC_23=SRC_23, SRC_23_RES=SRC_23_RES, SRC_25_RES=SRC_25_RES)),
                    fh, protocol=4)
    print(f"results -> {p}")
    # the labelled maps the cross-validated model disagrees with, as a relabel list
    # for 27 (python 27_label_gui.py --relabel <this file>); delete lines you trust
    p = config.var_path(f"relabel_{TAG}.txt")
    order = np.flatnonzero(PRED_CV != Y)
    order = order[np.argsort(-PCV[order, PRED_CV[order]])]
    with open(p, "w", encoding="utf-8") as fh:
        fh.write("# maps where the cross-validated classifier disagrees with your label\n"
                 "# one map per line; delete the ones you are sure of, add any others\n")
        for k in order:
            t = TR[k]
            fh.write(f"{t['uid'].replace('|', ' ')}{' (shuffle)' if t['kind'] == 'shuffle' else ''}"
                     f"   # you: {t['label']} | model: {CLASSES[PRED_CV[k]]} ({PCV[k, PRED_CV[k]]:.2f})\n")
    print(f"relabel list ({len(order)} maps) -> {p}")
    # every map shown in the example figures (steps 4, 5, 7), grouped as in the figures:
    # open them all in 27's gallery and click the ones to relabel
    p = config.var_path(f"examples_{TAG}.txt")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write("# the maps of 28's example figures; '## ...' lines are the gallery sections\n")
        grp = None
        for g, u, kind, com in EXAMPLES:
            if g != grp: fh.write(f"## {g}\n"); grp = g
            fh.write(f"{u.replace('|', ' ')}{' (shuffle)' if kind == 'shuffle' else ''}   # {com}\n")
    print(f"example maps ({len(EXAMPLES)}) -> {p}")

# %%
