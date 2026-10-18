# %% ── 0. imports & constants ─────────────────────────────────────────────────
# Multi-dimensional kinematics: instead of the single endpoint-y scalar used in
# 04/08, the kinematic variable here is the VECTOR [PC1, PC2, speed, accel] and
# the landscape is the JOINT MI( kin_vector ; spike count in window ).
#
# WHY A CUSTOM ESTIMATOR: sklearn's mutual_info_regression(X, y) with a (N, d)
# feature matrix returns d SEPARATE per-feature MIs -- not MI(X; y). Feeding it a
# 4-column X and summing (or taking the max) is a different, wrong quantity. This
# script uses a true multivariate KSG-1 (_ksg_md) that treats X as one joint
# variable under the Chebyshev norm. Validated in cell 3: at d=1 it reproduces
# sklearn exactly, and at d=2..4 it matches analytic Gaussian MI.
#
# KSG only (Ross needs a discrete variable); spike COUNT by default.
import os
import numpy as np
import matplotlib.pyplot as plt
import pickle, time, importlib
import mi_io, config, mi_rates
importlib.reload(config); importlib.reload(mi_rates)
from sklearn.neighbors import NearestNeighbors, KDTree
from sklearn.feature_selection import mutual_info_regression
from scipy.special import digamma
from scipy.signal import find_peaks, savgol_filter
from joblib import Parallel, delayed

if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

# EVERYTHING downstream works in MILLIMETRES, not pixels. This is not cosmetic:
# the px->mm ratio DIFFERS per axis (0.1, 0.07, 0.1), so in pixels the y axis is
# compressed relative to x and z. Any quantity that mixes axes -- speed, accel,
# the PCA, and the distance to water -- is therefore geometrically WRONG in px.
# (x, y, z on their own are unaffected: a per-axis rescale is monotone and _prep
# z-scores each column, so their MI is identical either way.)
PX2MM    = np.asarray(config.PX2MM, float)
REACH_MM = {k: v * PX2MM[None, :, None] for k, v in REACH_PX.items()}

CLR     = {"CP": np.array([106,176,134])/255, "CB": np.array([154,123,193])/255}
NAT2BIT = 1 / np.log(2)
rng     = np.random.default_rng(0)

# ── analysis parameters ───────────────────────────────────────────────────────
WIN_WIDTH   = 0.100       # sliding window width (s)
WIN_STEP    = 0.025       # step size (s)
WIN_RANGE   = (-1.0, 1.0) # landscape range (s), re reach
K_KSG       = 3           # k for KSG
CAT_KEEP    = {1, 2}      # reach categories to include
MIN_TRIALS  = 20

NEURAL_VAR  = "count"     # "count" (default) | "isi" | "conv" -- see mi_rates.py

# ── kinematic feature construction ────────────────────────────────────────────
FRAMES      = (147, 167)  # endpoint window, MATLAB 1-indexed (same as 04/08)
FPS         = 400.0       # camera frame rate
SAVGOL_WIN  = 11          # frames (27.5 ms) for the Savitzky-Golay derivative.
SAVGOL_ORD  = 3           # Differentiating the raw trace amplifies tracking jitter
                          # badly (naive diff gave speeds up to 1.2e5 px/s and an
                          # accel sd of 6.4e5, measured back when this ran in px);
                          # smoothing first cuts that ~5x. Everything is mm now.
# PROJECT DEFAULT  ->  FEATURES = ("x", "y", "z", "speed")
# Currently set to the PC variant so you can compare; restore the line above when done.
FEATURES    = ("x","y","z","speed","dist_water")
                          # Any subset/order of the features below works; everything downstream
                          # (PCA on/off, figures, output paths) follows FEATURES.
                          # Available: "x", "y", "z"  = endpoint-window mean of the raw
                          #            paw coordinate (y = med-lat, the 04/08 scalar)
                          #            "speed"        = ||d(pos)/dt|| over that window
                          #            "accel"        = d(speed)/dt over that window
                          #            "PC1", "PC2"   = trajectory PCs (see PCA_ON);
                          #                             fitted per PCA_SCOPE pool, so with
                          #                             the default "mouse" each animal
                          #                             gets its own axes
                          #            "dist_water"   = |paw - water| in MM over the window
                          #            "speed_water"  = d(dist_water)/dt in mm/s; NEGATIVE
                          #                             means closing on the water
                          #        The last two need 10_behavior_water.py to have been run
                          #        (they come from its pickle, WATER_PKL below). Every
                          #        feature here is in mm; what makes these two special is
                          #        that they know WHERE THE REWARD WAS -- x/y/z live in a
                          #        fixed camera frame, so they confound "where the paw went"
                          #        with "where the target was" across L/C/R trials.
                          #
                          # Measured peak MI, CP+CB, 180-neuron sample:
                          #   y only (= 04)      0.149   |  x+y              0.257
                          #   x+y+speed          0.280   |  x+y+z            0.297
                          #   x+y+z+speed        0.304   |  PC1+PC2          0.236
                          #   PC1+PC2+speed      0.247   |  +accel           0.218
                          #
                          # WHY THE DEFAULT IS RAW COORDINATES, NOT PCs: PCA maximises
                          # variance, not information about spiking. With PCA_ON=
                          # "trajectory" this went badly wrong -- PC1 was 45-79% reach-
                          # CATEGORY identity and scored pre-reach posture, separating at
                          # -300 ms and collapsing before the endpoint window. PCA_ON=
                          # "endpoint" (now the setting) fixes that by construction: the
                          # PCs are a rotation of [mean x, mean y, mean z] WITHIN the
                          # window, so nothing outside it can leak in.
                          # But note what that implies -- PC1+PC2 is then a 2-D rotation
                          # of a 3-D space, i.e. x+y+z with one dimension DROPPED, and MI
                          # is rotation-invariant. So PC1+PC2 can at best match x+y+z and
                          # will usually lose a little; the PCs buy interpretability of
                          # the axes, not information. Raw x/y/z need no PCA at all.
                          #
                          # WHY "accel" IS OUT OF THE DEFAULT: it carries only 0.006-0.009
                          # bits above its own shuffle null (vs ~0.10 for x/y), no filtering
                          # rescues it (pre-smoothing changes nothing, single-pass deriv=2
                          # halves the noise without adding information), and including it
                          # REDUCES the joint MI because a near-noise dimension inflates the
                          # joint Chebyshev distances. It stays selectable for exactly this
                          # kind of check.
                          #
                          # "z" earns its place: it gains more than speed does
                          # (x+y -> +z is +0.040, x+y -> +speed is +0.024), and the full
                          # x+y+z+speed set is the best measured (0.304 vs 0.149 for the
                          # 04 baseline of endpoint-y alone -- roughly double).

# ── trajectory-PCA settings (only used when PC1/PC2 are in FEATURES) ─────────
PCA_ON      = "endpoint"  # "endpoint"   -> each trial is [mean(x), mean(y), mean(z)] over
                          #     the endpoint window (3 dims). PCs are then a rotation
                          #     WITHIN that window, so they cannot be contaminated by
                          #     posture outside it. THIS IS THE WORKING CHOICE.
                          # "trajectory" -> each trial is the FLATTENED path (3 coords
                          #     x n_frames = 603 dims), so PCs describe the reach SHAPE.
                          #     REJECTED on this data: the per-frame projections showed
                          #     PC1/PC2 separating maximally at -300 to -250 ms and
                          #     collapsing to ~0 before the endpoint window opens, i.e.
                          #     they scored PRE-REACH POSTURE, not endpoint kinematics.
PCA_SCOPE   = "mouse"     # pool the PCs are fitted on: "mouse" | "session" | "global".
                          # "mouse" chosen from the data: only 2-5% of trajectory variance
                          # is BETWEEN sessions, and fitting 603 dims per session
                          # (330-1000 trials) is badly conditioned. Verified on the worst
                          # mouse (4_Milka) that the JOINT landscape is unchanged either
                          # way (<=0.001 bits), because MI is invariant to rotation of the
                          # kinematic basis; only the standalone PC2 curve is scope-
                          # dependent there.
TRAJ_FRAMES = None        # None = whole 201-frame trace; or (start, stop) to restrict
                          # the PCA to part of the reach.

WATER_PKL   = "behavior_water.pkl"   # written by 10_behavior_water.py; only read when
                          # "dist_water"/"speed_water" appear in FEATURES.

PER_FEATURE = True        # also compute the 1-D landscape for each feature on its
                          # own, so the joint can be judged against its parts

# ── per-window shuffle null (same convention as 04/07) ───────────────────────
DO_NULL     = True        # draw the null on the mean-landscape figures
N_SHUF_NULL = 10          # permutations per window per neuron. Costs one extra
                          # landscape x N_SHUF_NULL, so it roughly 10x's the joint
                          # run -- drop to 5, or set DO_NULL=False, while iterating.
NULL_PER_FEATURE = True   # also null the 1-D per-feature landscapes. Matters MORE
                          # here than in 04: KSG bias grows with d, so the joint and
                          # the 1-D curves do NOT share a baseline and cannot be
                          # compared on raw MI alone (see cell 9).

# ── peak / significance (same meaning as 04) ─────────────────────────────────
N_SHUF_LEADLAG    = 50
SMOOTH_WIN        = 3
PEAK_PROM_RATIO   = 0.5
PEAK_RUNNERUP_MAX = 0.7
PEAK_RANGE        = None

N_TEST      = None          # neurons PER REGION; None for the full run. Full run at
                          # d=4 is ~15 min joint + ~5.5 min per feature, parallelised.
N_JOBS      = -1
SHOW_FIGS   = True
SAVE_FIGS   = True
SAVE_VARS   = True

# Cells 5-7 are the expensive ones (joint landscape, per-feature landscapes,
# nulls, significance). CACHE_MODE decides whether to recompute them:
#   "load"  -> if <RUN_NAME_F>_ksg_multiD.pkl exists AND was written with the same
#              parameters, reuse it and SKIP cells 5-7 entirely. Cell 2 still runs
#              (it defines run_idx / regions_run / uids, which cells 8-10 need).
#   "rerun" -> always recompute, then overwrite the pickle if SAVE_VARS.
# Any parameter mismatch falls back to recomputing and says which key differed,
# so you can never silently plot a landscape from the wrong settings.
CACHE_MODE  = "load"

OUT_SUBDIR  = "multiD_kin"
RUN_NAME    = "multiD100"

# The feature set is appended to BOTH names, so switching FEATURES writes to its
# own folder and its own .pkl instead of overwriting the previous run:
#   multiD_kin_x-y-z-speed/        multiD100_x-y-z-speed_ksg_multiD.pkl
#   multiD_kin_PC1-PC2-speed-accel/  multiD100_PC1-PC2-speed-accel_ksg_multiD.pkl
FEAT_TAG    = "-".join(FEATURES)
OUT_DIR_F   = f"{OUT_SUBDIR}_{FEAT_TAG}"
RUN_NAME_F  = f"{RUN_NAME}_{FEAT_TAG}"

# defined here too, so cells 7-9 still run if you skip the null cell (6b)
null_joint, NULL_FEAT = None, {}

win_starts  = np.arange(WIN_RANGE[0], WIN_RANGE[1] - WIN_WIDTH + WIN_STEP/2, WIN_STEP)
win_centers = win_starts + WIN_WIDTH / 2
n_wins      = len(win_starts)
print(f"windows: {n_wins}  width={WIN_WIDTH*1000:.0f} ms  step={WIN_STEP*1000:.0f} ms")
print(f"kinematic features: {list(FEATURES)}  |  neural variable: {NEURAL_VAR}")
print(f"outputs -> {OUT_DIR_F}/  |  vars -> {RUN_NAME_F}_ksg_multiD.pkl")


# %% ── 1. helpers ────────────────────────────────────────────────────────────
def _traj_matrix(rp, frames=TRAJ_FRAMES):
    """(n_trials, 3*n_frames): each trial's x, y and z traces laid end to end.
    This is the row the trajectory PCA sees."""
    sub = rp if frames is None else rp[:, :, frames[0]:frames[1]]
    return sub.reshape(sub.shape[0], -1)


def _pca_input(rp):
    if PCA_ON == "trajectory": return _traj_matrix(rp)
    if PCA_ON == "endpoint":   return rp[:, :, EP].mean(axis=2)
    raise ValueError(f"PCA_ON must be 'trajectory' or 'endpoint', got {PCA_ON!r}")


def fit_pca(mats, n_comp=2):
    """Centre and SVD. Returns (mean, components, explained-variance-ratio).
    Centres but does NOT rescale: x, y, z are all in MM (a real, isotropic
    physical unit), so rescaling the axes would distort the geometry of the path.
    This only holds because the input is mm -- in pixels the axes have different
    scales and un-rescaled PCA would silently under-weight y."""
    X  = np.vstack(mats)
    mu = np.nanmean(X, axis=0)
    Z  = X - mu
    ok = np.isfinite(Z).all(axis=1)
    _, S, Vt = np.linalg.svd(Z[ok], full_matrices=False)
    return mu, Vt[:n_comp], S**2 / np.sum(S**2)


def _prep(A, r):
    """sklearn-style preprocessing: scale each column to unit sd (so no axis
    dominates the Chebyshev norm) and add the same 1e-10 jitter it uses to break
    ties. Without the scaling a 4-D max-norm is decided by whichever feature has
    the largest units (speed in mm/s vs positions in mm)."""
    A = np.atleast_2d(A).astype(float)
    A = A / np.maximum(A.std(axis=0, keepdims=True), 1e-12)
    return A + 1e-10 * np.maximum(1, np.abs(A).mean(axis=0, keepdims=True)) \
                     * r.standard_normal(A.shape)

def _ksg_md(X, y, k=3, seed=0):
    """TRUE multivariate KSG-1: MI(X ; y) for X (N,d) continuous and y (N,)
    continuous. Mirrors sklearn's _compute_mi_cc but keeps X d-dimensional, so
    at d=1 it is numerically identical to mutual_info_regression."""
    r = np.random.default_rng(seed)
    X = _prep(np.asarray(X, float).reshape(len(y), -1), r)
    y = _prep(np.asarray(y, float).reshape(-1, 1), r)
    n = X.shape[0]
    nn = NearestNeighbors(metric="chebyshev", n_neighbors=k).fit(np.hstack([X, y]))
    radius = np.nextafter(nn.kneighbors()[0][:, -1], 0)
    nx = np.array(KDTree(X, metric="chebyshev").query_radius(X, radius, count_only=True)) - 1.0
    ny = np.array(KDTree(y, metric="chebyshev").query_radius(y, radius, count_only=True)) - 1.0
    mi = digamma(n) + digamma(k) - np.mean(digamma(nx + 1)) - np.mean(digamma(ny + 1))
    return max(0.0, mi) * NAT2BIT

def _landscape_one(kin_m, C_m, k, seed):
    """MI landscape for ONE neuron -> (n_wins,). Top-level so joblib can ship it."""
    out = np.full(C_m.shape[1], np.nan)
    for wi in range(C_m.shape[1]):
        c_w = C_m[:, wi]
        if np.unique(c_w).size < 2: continue
        out[wi] = _ksg_md(kin_m, c_w, k, seed + wi)
    return out

def _null_one_md(kin_m, C_m, n_shuf, seed, k):
    """Per-window shuffle null for ONE neuron -> (n_wins,), the MEAN joint MI over
    `n_shuf` trial-permutations of the kinematic matrix at every window. Same
    convention as 04/07: whole kinematic ROWS are permuted, so the trial pairing
    is destroyed while the internal structure of the kinematic vector survives --
    the null therefore keeps the d-dimensional KSG bias that makes the raw joint
    MI look large, which is exactly the offset we want to subtract.
    The same permutations are reused across windows on purpose: each window's null
    stays unbiased and the curve is smoother window-to-window than redrawing."""
    r    = np.random.default_rng(seed)
    shuf = [r.permutation(kin_m) for _ in range(n_shuf)]
    out  = np.full(C_m.shape[1], np.nan)
    for wi in range(C_m.shape[1]):
        c_w = C_m[:, wi]
        if np.unique(c_w).size < 2: continue
        out[wi] = np.mean([_ksg_md(ks, c_w, k, seed + wi) for ks in shuf])
    return out

def _sig_one(kin_m, C_m, mi_arr, k, n_shuf, seed, smooth_win,
             prom_ratio, runnerup_max):
    """Peak-window shuffle test for ONE neuron -> bool. The whole kinematic ROW is
    permuted, so the trial pairing is destroyed while the internal structure of
    the kinematic vector is preserved."""
    if np.all(np.isnan(mi_arr)): return False
    pk = _find_peak(_smooth(mi_arr, smooth_win), prom_ratio, runnerup_max, smooth_win)
    if pk is None: return False
    wi = pk[0]; real = mi_arr[wi]
    if real <= 0: return False
    c_pk = C_m[:, wi]
    if np.unique(c_pk).size < 2: return False
    r = np.random.default_rng(seed)
    null = np.array([_ksg_md(kin_m[r.permutation(len(kin_m))], c_pk, k, seed + j)
                     for j in range(n_shuf)])
    return bool(real > np.percentile(null, 95))

def _smooth(arr, win):
    if win <= 1: return arr.copy()
    n = len(arr); out = np.full(n, np.nan); half = win // 2
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        seg = arr[lo:hi]
        if np.any(~np.isnan(seg)): out[i] = np.nanmean(seg)
    return out

def _find_peak(arr, prom_ratio=PEAK_PROM_RATIO, runnerup_max=PEAK_RUNNERUP_MAX,
               smooth_win=SMOOTH_WIN):
    valid = ~np.isnan(arr)
    if valid.sum() < 3: return None
    a = np.where(valid, arr, -np.inf)
    idx_peaks, props = find_peaks(a, prominence=0, distance=max(1, smooth_win))
    if len(idx_peaks) == 0: return None
    proms = props["prominences"]; heights = a[idx_peaks]
    order = np.argsort(heights)[::-1]; top = order[0]
    top_h, top_prom = heights[top], proms[top]
    runner_h = heights[order[1]] if len(order) > 1 else 0.0
    if top_h <= 0 or top_prom / top_h < prom_ratio or runner_h > runnerup_max * top_h:
        return None
    return int(idx_peaks[top]), float(top_prom / top_h), float(top_h)


# %% ── 2. build the kinematic feature matrix ─────────────────────────────────
EP = slice(FRAMES[0] - 1, FRAMES[1])          # MATLAB 1-idx -> python

if PCA_SCOPE == "session":
    POOLS = {k: [k] for k in REACH_PX}
elif PCA_SCOPE == "mouse":
    POOLS = {}
    for k in REACH_PX: POOLS.setdefault(k[0], []).append(k)
elif PCA_SCOPE == "global":
    POOLS = {"ALL": list(REACH_PX)}
else:
    raise ValueError(f"PCA_SCOPE must be mouse|session|global, got {PCA_SCOPE!r}")

NEEDS_PCA = any(f in ("PC1", "PC2") for f in FEATURES)
PCA_FIT  = {pool: fit_pca([_pca_input(REACH_MM[k]) for k in keys])
            for pool, keys in POOLS.items()} if NEEDS_PCA else {}
KEY_POOL = {k: pool for pool, keys in POOLS.items() for k in keys}

if NEEDS_PCA:
    evrs = np.vstack([PCA_FIT[p][2][:3] for p in POOLS])
    print(f"PCA on {PCA_ON!r} fitted per {PCA_SCOPE!r}: {len(POOLS)} pool(s), "
          f"input dim = {_pca_input(REACH_MM[list(REACH_MM)[0]]).shape[1]}")
    print(f"  explained variance PC1..3 (mean over pools): {np.round(evrs.mean(0), 3)}"
          f"   PC1+PC2 = {evrs.mean(0)[:2].sum()*100:.1f}%")
    if PCA_ON == "trajectory":
        print("  NOTE trajectory PC1 is 45-79% reach-category identity on this data "
              "and scores PRE-REACH posture; prefer PCA_ON='endpoint'.")
    else:
        print("  PCs are a rotation of [mean x, mean y, mean z] within the endpoint "
              "window -> no leakage from outside it, but PC1+PC2 drops one of the "
              "three dimensions, so expect <= x+y+z (MI is rotation-invariant).")
else:
    print("no PC in FEATURES -> PCA skipped entirely (raw coordinates need none)")

# --- water-referenced features (built by 10_behavior_water.py) ---------------
NEEDS_WATER = any(f in ("dist_water", "speed_water") for f in FEATURES)
WATER = {}
if NEEDS_WATER:
    _wp = config.var_path(WATER_PKL)
    if not _wp.is_file():
        raise FileNotFoundError(
            f"{_wp} not found -- run 10_behavior_water.py first (it builds the "
            f"water position, distance and speed-to-water per reach).")
    with open(_wp, "rb") as _f: WATER = pickle.load(_f)["sessions"]
    _miss = [k for k in REACH_PX if k not in WATER]
    if _miss: raise KeyError(f"sessions missing from {WATER_PKL}: {_miss}")
    # the reach arrays must be the SAME reaches, or the water is on the wrong rows
    for _k in REACH_PX:
        if WATER[_k]["dist_water"].shape[0] != REACH_PX[_k].shape[0]:
            raise ValueError(f"{_k}: {WATER[_k]['dist_water'].shape[0]} reaches in "
                             f"{WATER_PKL} vs {REACH_PX[_k].shape[0]} in REACH_PX")
    _n = sum(int(np.isfinite(WATER[k]["water_mm"]).all(1).sum()) for k in REACH_PX)
    _t = sum(REACH_PX[k].shape[0] for k in REACH_PX)
    print(f"water features loaded: {_n}/{_t} reaches have a water position "
          f"({_n/_t*100:.1f}%); the rest become NaN and are masked out")

def session_features(key):
    """-> (n_trials, len(FEATURES)) plus the smoothed speed/accel traces.
    PC scores project this session's trials onto the components fitted on its POOL,
    so every session of a mouse shares one set of axes."""
    rp = REACH_MM[key]
    if NEEDS_PCA:
        mu, comp, _ = PCA_FIT[KEY_POOL[key]]
        pcs = (_pca_input(rp) - mu) @ comp.T
    else:
        pcs = np.full((rp.shape[0], 2), np.nan)
    vel = savgol_filter(rp, SAVGOL_WIN, SAVGOL_ORD, deriv=1, delta=1/FPS, axis=2)
    spd = np.linalg.norm(vel, axis=1)
    acc = savgol_filter(spd, SAVGOL_WIN, SAVGOL_ORD, deriv=1, delta=1/FPS, axis=1)
    cols = {"x": rp[:, 0, EP].mean(axis=1), "y": rp[:, 1, EP].mean(axis=1),
            "z": rp[:, 2, EP].mean(axis=1),
            "PC1": pcs[:, 0], "PC2": pcs[:, 1],
            "speed": spd[:, EP].mean(axis=1), "accel": acc[:, EP].mean(axis=1)}
    if NEEDS_WATER:                       # in mm / mm-per-s, NOT pixels
        W = WATER[key]
        cols["dist_water"]  = W["dist_water"][:, EP].mean(axis=1)
        cols["speed_water"] = W["speed_water"][:, EP].mean(axis=1)
    return np.column_stack([cols[f] for f in FEATURES]), spd, acc

SESS_FEAT, SESS_SPD, SESS_ACC = {}, {}, {}
for key in REACH_PX:
    SESS_FEAT[key], SESS_SPD[key], SESS_ACC[key] = session_features(key)

uids    = [n["uid"] for n in NEURONS]
regions = [n["region"] for n in NEURONS]
kin_list, mask_list = [], []
for n in NEURONS:
    F   = SESS_FEAT[(n["mouse"], n["sess"])]
    cat = n["labels"]["idx_reach_cat"]
    kin_list.append(F)
    mask_list.append(np.isfinite(F).all(axis=1) & np.isin(cat, list(CAT_KEEP)))

cp_idx = [i for i, r in enumerate(regions) if r == "CP"]
cb_idx = [i for i, r in enumerate(regions) if r == "CB"]
run_idx = (cp_idx[:N_TEST] + cb_idx[:N_TEST]) if N_TEST is not None           else list(range(len(NEURONS)))
regions_run = [regions[i] for i in run_idx]
print(f"{'TEST' if N_TEST else 'FULL'} RUN: {len(run_idx)} neurons | "
      f"feature matrix per neuron: (n_trials, {len(FEATURES)})")


# %% ── 3. INSPECT the kinematic features ────────────────────────────────────
# One figure: what each feature is read from -- the raw coordinate / derived
# signal it averages, or (for PCs) the trajectory pattern it encodes.
# ROWS = mouse, COLUMNS = feature.
n_frames = REACH_MM[list(REACH_MM)[0]].shape[2]
tms    = (np.arange(n_frames) - EP.start) / FPS * 1000
epw    = (EP.stop - EP.start) / FPS * 1000
COORD  = ("paw x", "paw y (med-lat)", "paw z")   # all in mm
UNITS  = {"x": "mm", "y": "mm", "z": "mm", "PC1": "mm", "PC2": "mm",
          "speed": "mm/s", "accel": "mm/s^2",
          "dist_water": "mm", "speed_water": "mm/s"}

def _quartiles(v):
    ok = np.isfinite(v)
    q1, q3 = np.nanpercentile(v[ok], [25, 75])
    return ok & (v <= q1), ok & (v >= q3)

# --- what each feature is computed from, ONE ROW PER MOUSE ---
# Columns = the FEATURES; rows = mice, so idiosyncratic trajectories are visible
# side by side. Within a panel, trials are split into the bottom/top quartile OF
# THAT COLUMN'S FEATURE, computed within that mouse.
MICE = sorted({k[0] for k in REACH_PX})       # rows, independent of PCA_SCOPE

fig, axes = plt.subplots(len(MICE), len(FEATURES),
                         figsize=(5*len(FEATURES), 3.4*len(MICE)), sharex=True)
# reshape, NOT atleast_2d: with a single FEATURE, subplots returns shape
# (n_mice,) and atleast_2d would make it (1, n_mice) -- transposed, so axes[r, c]
# then indexes the wrong panel and crashes for the second mouse.
axes = np.asarray(axes).reshape(len(MICE), len(FEATURES))
for r, mouse in enumerate(MICE):
    mkeys = [k for k in REACH_PX if k[0] == mouse]
    rp_m  = np.vstack([REACH_MM[k]  for k in mkeys])
    F_m   = np.vstack([SESS_FEAT[k] for k in mkeys])
    spd_m = np.vstack([SESS_SPD[k]  for k in mkeys])
    acc_m = np.vstack([SESS_ACC[k]  for k in mkeys])
    for c, fn in enumerate(FEATURES):
        ax  = axes[r, c]
        src = {"x": rp_m[:, 0, :], "y": rp_m[:, 1, :], "z": rp_m[:, 2, :],
               "speed": spd_m, "accel": acc_m}.get(fn)
        if src is None and NEEDS_WATER and fn in ("dist_water", "speed_water"):
            src = np.vstack([WATER[k][fn] for k in mkeys])
        if src is None and fn in ("PC1", "PC2"):
            # Show the PC as a SIGNAL OVER TIME, like every other column, rather
            # than as a paw x/y/z reconstruction: project each frame onto the
            # loading, so the curve is the score unrolled in time and the blue/red
            # gap is what the feature measures.
            #   PCA_ON="endpoint"   -> one 3-vector loading, applied at every frame;
            #                          the curve's MEAN over the endpoint window is
            #                          the trial's score.
            #   PCA_ON="trajectory" -> a loading per frame; the SUM over frames is
            #                          the score.
            j = 0 if fn == "PC1" else 1
            # evr is per POOL, i.e. per mouse with PCA_SCOPE="mouse", so each row
            # of this figure gets its own number -- hence it goes on the y label,
            # not the shared column title.
            mu, comp, _evr = PCA_FIT[KEY_POOL[mkeys[0]]]
            pc_evr = float(_evr[j]) * 100
            if PCA_ON == "endpoint":          # mu, comp[j] are length-3
                base = mu.reshape(3, 1); load = comp[j].reshape(3, 1)
                nf   = rp_m.shape[2]
                src  = np.einsum("ncf,c->nf", rp_m - base[None], load[:, 0])
            else:                             # length 3*n_frames
                base = mu.reshape(3, -1); load = comp[j].reshape(3, -1)
                nf   = base.shape[1]
                src  = np.einsum("ncf,cf->nf", rp_m[:, :, :nf] - base[None], load)
        if src is not None:
            is_pc = fn in ("PC1", "PC2")
            tt = tms[:src.shape[1]]               # TRAJ_FRAMES may shorten a PC trace
            lo_m, hi_m = _quartiles(F_m[:, c])
            for m, col, tag in ((lo_m, "#3b6ea5", "bottom 25%"),
                                (hi_m, "#c1442e", "top 25%")):
                mn = np.nanmean(src[m], axis=0)
                se = np.nanstd(src[m], axis=0)/np.sqrt(max(m.sum(), 1))
                ax.plot(tt, mn, color=col, lw=1.8, label=f"{tag} (n={m.sum()})")
                ax.fill_between(tt, mn-se, mn+se, color=col, alpha=0.25)
            if is_pc: ax.axhline(0, color="k", lw=0.5)
            if is_pc:
                ax.set_ylabel(f"{fn} per-frame projection ({UNITS.get(fn,'')})"
                              f"\n{pc_evr:.1f}% of variance")
            else:
                ax.set_ylabel(f"{fn} ({UNITS.get(fn,'')})")
            if r == 0:
                ax.set_title(f"{fn}  --  " + ("the score, unrolled in time"
                             if is_pc else "the signal it averages"),
                             fontsize=12, pad=10)
        else:                                     # PC1 / PC2
            j = 0 if fn == "PC1" else 1
            mu, comp, _ = PCA_FIT[KEY_POOL[mkeys[0]]]
            base = mu.reshape(3, -1); load = comp[j].reshape(3, -1)
            sc   = 2*np.nanstd(F_m[:, c]); tt = tms[:base.shape[1]]
            pc_evr = float(PCA_FIT[KEY_POOL[mkeys[0]]][2][j]) * 100
            for k2, (nm, col) in enumerate(zip(COORD, ("#222", "#777", "#bbb"))):
                ax.plot(tt, base[k2] + sc*load[k2], color=col, lw=1.8,
                        label=f"{nm}: mean +2sd")
                ax.plot(tt, base[k2] - sc*load[k2], color=col, lw=1.4, ls="--",
                        label=f"{nm}: mean -2sd")
            ax.set_ylabel(f"paw position (mm)\n{fn}: {pc_evr:.1f}% of variance")
            if r == 0:
                ax.set_title(f"{fn}  --  the trajectory pattern it encodes",
                             fontsize=12, pad=10)
        ax.legend(fontsize=6, ncol=1 if src is not None else 2)
        ax.axvspan(0, epw, color="k", alpha=0.12)
        if r == len(MICE)-1: ax.set_xlabel("time re endpoint-window start (ms)")
        ax.spines[["top", "right"]].set_visible(False)
    axes[r, 0].annotate(mouse, xy=(-0.26, 0.5), xycoords="axes fraction",
                        fontsize=13, fontweight="bold", ha="center", va="center",
                        rotation=90, color="#333")
fig.suptitle("What each feature is computed from  --  ROWS = mouse, "
             "COLUMNS = kinematic feature\n"
             "blue = that feature's bottom quartile within the mouse, red = top "
             "quartile;  grey band = endpoint window", fontsize=13, y=1.002)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir(OUT_DIR_F) / "kin_feature_sources.png", dpi=140,
                bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

allF = np.vstack([SESS_FEAT[k] for k in SESS_FEAT])
print(f"{'feature':8s} {'median':>12s} {'IQR':>26s} {'range':>26s}")
for j, nm in enumerate(FEATURES):
    v = allF[:, j]; v = v[np.isfinite(v)]
    print(f"{nm:8s} {np.median(v):12.2f}  "
          f"[{np.percentile(v,25):11.2f},{np.percentile(v,75):11.2f}]  "
          f"[{v.min():11.1f},{v.max():11.1f}]")


# %% ── 4. validate the multivariate estimator (run once, keep it honest) ─────
# d=1 must reproduce sklearn exactly; d>1 must track analytic Gaussian MI.
_r = np.random.default_rng(0)
print("d=1 vs sklearn (must be identical):")
for rho in (0.3, 0.6, 0.85):
    x = _r.normal(size=400); yv = rho*x + np.sqrt(1-rho**2)*_r.normal(size=400)
    sk = mutual_info_regression(x.reshape(-1,1), yv, n_neighbors=3, random_state=0)[0]*NAT2BIT
    print(f"   rho={rho}: sklearn={sk:.4f}  _ksg_md={_ksg_md(x.reshape(-1,1), yv, 3):.4f}")
print("d>1 vs analytic Gaussian MI:")
for d, rho in ((2, 0.6), (4, 0.7)):
    A = _r.normal(size=(d,)); X = _r.normal(size=(4000, d))
    yv = X @ A * rho/np.linalg.norm(A) + np.sqrt(1-rho**2)*_r.normal(size=4000)
    S = np.cov(np.hstack([X, yv.reshape(-1,1)]).T)
    ana = 0.5*np.log2(np.linalg.det(S[:d,:d])*S[d,d]/np.linalg.det(S))
    print(f"   d={d} rho={rho}: analytic={ana:.4f}  _ksg_md={_ksg_md(X, yv, 3):.4f}")


# %% ── 5. joint MI landscape: MI([PC1,PC2,speed,accel] ; count) ──────────────
# Try the cache before doing any work. _LOADED gates cells 6, 6b and 7 as well.
_CACHE_P = config.var_path(f"{RUN_NAME_F}_ksg_multiD.pkl")
_LOADED  = False
if CACHE_MODE == "load" and os.path.exists(_CACHE_P):
    with open(_CACHE_P, "rb") as f: _c = pickle.load(f)
    # Compare only what actually changes the numbers. Figure/style parameters are
    # deliberately NOT compared, so a plot can be restyled without a recompute.
    _want = dict(features=list(FEATURES), WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP,
                 K=K_KSG, NEURAL_VAR=NEURAL_VAR, FRAMES=tuple(FRAMES),
                 CAT_KEEP=list(CAT_KEEP), SAVGOL_WIN=SAVGOL_WIN,
                 SAVGOL_ORD=SAVGOL_ORD, N_TEST=N_TEST)
    _bad = [k for k, v in _want.items() if list(np.atleast_1d(_c.get(k))) != list(np.atleast_1d(v))]
    if _bad or len(_c["win_centers"]) != n_wins:
        print(f"cache MISMATCH on {_bad or ['n_wins']} -> recomputing")
        for k in _bad:
            print(f"    {k}: cached={_c.get(k)!r}  now={_want[k]!r}")
    else:
        mi_joint   = _c["mi"]
        sig_joint  = _c["sig"]
        MI_FEAT    = dict(_c.get("mi_per_feature") or {})
        null_joint = _c.get("null")
        NULL_FEAT  = dict(_c.get("null_per_feature") or {})
        _LOADED    = True
        print(f"loaded cache <- {_CACHE_P}")
        print(f"  {mi_joint.shape[0]} neurons x {mi_joint.shape[1]} windows | "
              f"null={'yes' if null_joint is not None else 'no'} "
              f"(S={_c.get('N_SHUF_NULL')}) | per-feature={list(MI_FEAT)}")
        print("  -> SKIP cells 6, 6b and 7; go straight to cell 8")
elif CACHE_MODE == "load":
    print(f"no cache at {_CACHE_P} -> computing")

payload, keep = [], []
if not _LOADED:
    payload = []
    for pi, i in enumerate(run_idx):
        mask = mask_list[i]
        if mask.sum() < MIN_TRIALS:
            payload.append(None); continue
        C_m = np.asarray(mi_rates.windowed(NEURONS[i]["spikes"], win_starts, WIN_WIDTH,
                                           NEURAL_VAR), float)[mask]
        payload.append((kin_list[i][mask], C_m))
    keep = [pi for pi, it in enumerate(payload) if it is not None]
    print(f"\njoint {len(FEATURES)}-D landscape for {len(keep)} neurons "
          f"(n_jobs={N_JOBS})...")

    t0 = time.time()
    res_j = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
        delayed(_landscape_one)(payload[pi][0], payload[pi][1], K_KSG, 100 + pi)
        for pi in keep)
    mi_joint = np.full((len(run_idx), n_wins), np.nan)
    for pi, v in zip(keep, res_j): mi_joint[pi] = v
    print(f"joint landscape done in {time.time()-t0:.0f}s")


# %% ── 6. per-feature 1-D landscapes (so the joint can be judged vs its parts) ─
# Guarded: re-running this after a cache load would overwrite loaded MI_FEAT.
if not _LOADED:
    MI_FEAT = {}
    if PER_FEATURE:
        for fi, fname in enumerate(FEATURES):
            t0 = time.time()
            res_f = Parallel(n_jobs=N_JOBS, backend="loky", verbose=0)(
                delayed(_landscape_one)(payload[pi][0][:, [fi]], payload[pi][1],
                                        K_KSG, 200 + pi)
                for pi in keep)
            arr = np.full((len(run_idx), n_wins), np.nan)
            for pi, v in zip(keep, res_f): arr[pi] = v
            MI_FEAT[fname] = arr
            print(f"  1-D {fname:6s} done in {time.time()-t0:.0f}s")


# %% ── 6b. per-window shuffle nulls ──────────────────────────────────────────
# Guarded: re-running this after a cache load would overwrite the loaded null.
if not _LOADED:
    null_joint = None
    NULL_FEAT  = {}
    if DO_NULL:
        t0 = time.time()
        res_n = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
            delayed(_null_one_md)(payload[pi][0], payload[pi][1], N_SHUF_NULL,
                                  1000 + pi, K_KSG)
            for pi in keep)
        null_joint = np.full((len(run_idx), n_wins), np.nan)
        for pi, v in zip(keep, res_n): null_joint[pi] = v
        print(f"joint null done in {time.time()-t0:.0f}s (S={N_SHUF_NULL})")

        if PER_FEATURE and NULL_PER_FEATURE:
            for fi, fname in enumerate(FEATURES):
                t0 = time.time()
                res_nf = Parallel(n_jobs=N_JOBS, backend="loky", verbose=0)(
                    delayed(_null_one_md)(payload[pi][0][:, [fi]], payload[pi][1],
                                          N_SHUF_NULL, 2000 + pi, K_KSG)
                    for pi in keep)
                arr = np.full((len(run_idx), n_wins), np.nan)
                for pi, v in zip(keep, res_nf): arr[pi] = v
                NULL_FEAT[fname] = arr
                print(f"  1-D null {fname:6s} done in {time.time()-t0:.0f}s")

        # How much of the raw joint MI is just d-dimensional KSG bias?
        for reg in ("CP", "CB"):
            idx = [i for i, r in enumerate(regions_run) if r == reg]
            if not idx: continue
            o = np.nanmean(mi_joint[idx], axis=0); nu = np.nanmean(null_joint[idx], axis=0)
            print(f"  {reg}: peak observed={np.nanmax(o):.4f}  null there="
                  f"{nu[np.nanargmax(o)]:.4f}  excess={np.nanmax(o)-nu[np.nanargmax(o)]:.4f} bits")


# %% ── 7. peak significance for the joint landscape ──────────────────────────
# Guarded: sig_joint comes from the cache when one was loaded, and there is
# nothing new to write, so the pickle is left untouched.
if not _LOADED:
    t0 = time.time()
    res_s = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
        delayed(_sig_one)(payload[pi][0], payload[pi][1], mi_joint[pi], K_KSG,
                          N_SHUF_LEADLAG, 300 + pi, SMOOTH_WIN,
                          PEAK_PROM_RATIO, PEAK_RUNNERUP_MAX)
        for pi in keep)
    sig_joint = np.zeros(len(run_idx), dtype=bool)
    for pi, v in zip(keep, res_s): sig_joint[pi] = v
    print(f"significance: {sig_joint.sum()}/{len(run_idx)} ({time.time()-t0:.0f}s)")

if SAVE_VARS and not _LOADED:
    out = dict(method="ksg_multiD", features=list(FEATURES), mi=mi_joint,
               sig=sig_joint, mi_per_feature={k: v for k, v in MI_FEAT.items()},
               null=null_joint, null_per_feature={k: v for k, v in NULL_FEAT.items()},
               N_SHUF_NULL=(N_SHUF_NULL if DO_NULL else None),
               win_centers=win_centers, win_starts=win_starts,
               WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP, K=K_KSG,
               uids=[uids[i] for i in run_idx], regions=regions_run,
               FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), NEURAL_VAR=NEURAL_VAR,
               SAVGOL_WIN=SAVGOL_WIN, SAVGOL_ORD=SAVGOL_ORD,
               N_SHUF_LEADLAG=N_SHUF_LEADLAG, N_TEST=N_TEST)
    p = config.var_path(f"{RUN_NAME_F}_ksg_multiD.pkl")
    with open(p, "wb") as f: pickle.dump(out, f, protocol=4)
    print(f"saved -> {p}")


# %% ── 8. mean landscape + lead/lag (same style as 04) ───────────────────────
fdir = config.fig_dir(OUT_DIR_F)
lab  = f"KSG k={K_KSG}, {NEURAL_VAR}, joint [{', '.join(FEATURES)}]"

# Say so LOUDLY rather than silently drawing a figure with no null in it.
if null_joint is None:
    print("!! no shuffle null on the mean landscape -- run cell 6b "
          f"(DO_NULL is currently {DO_NULL})")
else:
    print(f"mean landscape includes the shuffle null (S={N_SHUF_NULL})")

fig, ax = plt.subplots(figsize=(8, 4))
for reg in ("CP", "CB"):
    idx = [i for i, r in enumerate(regions_run) if r == reg]
    if not idx: continue
    m  = mi_joint[idx]; ok = ~np.all(np.isnan(m), axis=1); m = m[ok]
    if m.shape[0] == 0: continue
    mn  = np.nanmean(m, axis=0)
    sem = np.nanstd(m, axis=0) / np.sqrt(np.maximum(np.sum(~np.isnan(m), axis=0), 1))
    ax.plot(win_centers, mn, color=CLR[reg], lw=2, label=f"{reg} (n={ok.sum()})")
    ax.fill_between(win_centers, mn-sem, mn+sem, color=CLR[reg], alpha=0.2)
    if null_joint is not None:
        nm   = null_joint[idx][ok]
        nmn  = np.nanmean(nm, axis=0)
        nsem = np.nanstd(nm, axis=0) / np.sqrt(
            np.maximum(np.sum(~np.isnan(nm), axis=0), 1))
        ax.plot(win_centers, nmn, color=CLR[reg], lw=1.4, ls="--", alpha=0.55,
                label=f"{reg} shuffle null")
        ax.fill_between(win_centers, nmn-nsem, nmn+nsem, color=CLR[reg], alpha=0.08)
    pk = win_centers[np.nanargmax(mn)]
    ax.axvline(pk, color=CLR[reg], lw=1.5, ls="--", alpha=0.85,
               label=f"{reg} peak = {pk*1000:.0f} ms")
ax.axvline(0, color="k", lw=0.8, ls=":", label="reach endpoint")
ax.axhline(0, color="k", lw=0.4)
ax.set(xlabel="window center re reach (s)", ylabel="mean MI ± SEM (bits)",
       title=f"{lab}  —  mean landscape")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.legend(fontsize=7)
plt.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "mean_landscape_joint.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

lo, hi = PEAK_RANGE if PEAK_RANGE is not None else (win_centers[0], win_centers[-1])
win_mask = (win_centers >= lo) & (win_centers <= hi)
bin_edges = np.linspace(lo, hi, 21)
fig, ax = plt.subplots(figsize=(6, 4))
for reg in ("CP", "CB"):
    peaks = []
    for pi, r in enumerate(regions_run):
        if r != reg or not sig_joint[pi]: continue
        row = mi_joint[pi]
        if not np.any(~np.isnan(row)): continue
        pk = _find_peak(np.where(win_mask, _smooth(row, SMOOTH_WIN), np.nan))
        if pk is not None: peaks.append(win_centers[pk[0]])
    if not peaks: continue
    ax.hist(peaks, bins=bin_edges, color=CLR[reg], alpha=0.6, edgecolor="white",
            label=f"{reg} (n={len(peaks)})")
    med = np.median(peaks)
    ax.axvline(med, color=CLR[reg], lw=1.5, ls="--", alpha=0.9,
               label=f"{reg} median = {med*1000:.0f} ms")
ax.axvline(0, color="k", lw=0.8, ls=":", label="reach endpoint")
ax.set(xlabel="time of peak MI re reach (s)", ylabel="# neurons",
       title=f"{lab}  —  lead / lag")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.legend(fontsize=8)
plt.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "lead_lag_joint.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 9. does the joint buy anything over the single features? ──────────────
if PER_FEATURE:
    # cycle, so FEATURES longer than the palette still gets a colour for every
    # entry (zip() silently truncated and raised KeyError on the 5th feature)
    _PAL = ["#d95f02", "#1f78b4", "#7570b3", "#66a61e", "#e7298a", "#a6761d", "#666666"]
    FCLR = {f: _PAL[i % len(_PAL)] for i, f in enumerate(FEATURES)}
    # THE COMPARISON THAT NEEDS THE NULL. The joint curve is d-dimensional and the
    # feature curves are 1-D, and KSG bias grows with d -- so on RAW MI the joint
    # is flattered by construction. With nulls available, subtract each curve's own
    # null and the comparison becomes like-for-like.
    have_nulls = bool(NULL_FEAT) and null_joint is not None
    J   = (mi_joint - null_joint) if have_nulls else mi_joint
    FMI = ({f: MI_FEAT[f] - NULL_FEAT[f] for f in FEATURES} if have_nulls
           else MI_FEAT)
    ylab = "mean Excess MI (bits)" if have_nulls else "mean MI (bits)"

    fig, axes = plt.subplots(1, 2, figsize=(13, 4), sharey=True)
    for ax, reg in zip(axes, ("CP", "CB")):
        idx = [i for i, r in enumerate(regions_run) if r == reg]
        if not idx: continue
        m = J[idx]; ok = ~np.all(np.isnan(m), axis=1)
        ax.plot(win_centers, np.nanmean(m[ok], axis=0), color="k", lw=2.4,
                label=f"JOINT {len(FEATURES)}-D")
        for fname in FEATURES:
            mf = FMI[fname][idx]
            ax.plot(win_centers, np.nanmean(mf[ok], axis=0), lw=1.5, alpha=0.85,
                    color=FCLR[fname], label=fname)
        ax.axvline(0, color="k", lw=0.8, ls=":")
        if have_nulls: ax.axhline(0, color="k", lw=0.8)
        ax.set(xlabel="window center re reach (s)", title=f"{reg} (n={ok.sum()})")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.legend(fontsize=8)
    axes[0].set_ylabel(ylab)
    plt.suptitle(f"Joint vs single-feature "
                 f"{'EXCESS ' if have_nulls else ''}MI  —  {NEURAL_VAR}, "
                 f"KSG k={K_KSG}", fontsize=11)
    plt.tight_layout()
    if SAVE_FIGS:
        fig.savefig(fdir / ("joint_vs_features_excess.png" if have_nulls
                            else "joint_vs_features.png"),
                    dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

    hdr = "peak excess" if have_nulls else "peak MI"
    print(f"\n{'region':7s} {'variable':10s} {hdr:>12s} {'peak t (ms)':>12s}"
          + ("  " + f"{'raw MI':>9s} {'null':>9s}" if have_nulls else ""))
    for reg in ("CP", "CB"):
        idx = [i for i, r in enumerate(regions_run) if r == reg]
        if not idx: continue
        ok = ~np.all(np.isnan(mi_joint[idx]), axis=1)
        rows = [("JOINT", J, mi_joint, null_joint)] + \
               [(f, FMI[f], MI_FEAT[f], NULL_FEAT.get(f)) for f in FEATURES]
        for nm, arr, raw, nul in rows:
            mn = np.nanmean(arr[idx][ok], axis=0)
            j  = int(np.nanargmax(mn))
            extra = ""
            if have_nulls and nul is not None:
                extra = (f"  {np.nanmean(raw[idx][ok], axis=0)[j]:9.4f} "
                         f"{np.nanmean(nul[idx][ok], axis=0)[j]:9.4f}")
            print(f"{reg:7s} {nm:10s} {np.nanmax(mn):12.4f} "
                  f"{win_centers[j]*1000:12.0f}{extra}")
print(f"\nfigures -> {fdir}")



# %% ── 10. per-mouse panels (same layout as 08) ──────────────────────────────
# Does the pooled result hold in every animal, or is it carried by one mouse?
# Grid of one panel per mouse plus a pooled "ALL MICE" panel last.
MIN_NEURONS_MOUSE = 20        # skip a region within a mouse below this

# recomputed here (not inherited from cell 8) so this cell runs on its own
lo, hi    = PEAK_RANGE if PEAK_RANGE is not None else (win_centers[0], win_centers[-1])
win_mask  = (win_centers >= lo) & (win_centers <= hi)
bin_edges = np.linspace(lo, hi, 21)

mice_run = [NEURONS[i]["mouse"] for i in run_idx]
MOUSE_LIST = sorted(set(mice_run))

def _peaks_for(sel):
    """Peak times of the significant, genuinely-peaked neurons in sel."""
    out = []
    for pi in sel:
        if not sig_joint[pi]: continue
        row = mi_joint[pi]
        if not np.any(~np.isnan(row)): continue
        pk = _find_peak(np.where(win_mask, _smooth(row, SMOOTH_WIN), np.nan))
        if pk is not None: out.append(win_centers[pk[0]])
    return out

def _draw_mean_landscape(ax, sel, title, fontsize=6, xlabel=True, ylabel=True):
    for reg in ("CP", "CB"):
        idx = [pi for pi in sel if regions_run[pi] == reg]
        if len(idx) < MIN_NEURONS_MOUSE: continue
        m = mi_joint[idx]; ok = ~np.all(np.isnan(m), axis=1); m = m[ok]
        if m.shape[0] == 0: continue
        mn  = np.nanmean(m, axis=0)
        sem = np.nanstd(m, axis=0)/np.sqrt(np.maximum(np.sum(~np.isnan(m), axis=0), 1))
        ax.plot(win_centers, mn, color=CLR[reg], lw=2, label=f"{reg} (n={ok.sum()})")
        ax.fill_between(win_centers, mn-sem, mn+sem, color=CLR[reg], alpha=0.2)
        if null_joint is not None:
            nm   = null_joint[idx][ok]
            nmn  = np.nanmean(nm, axis=0)
            nsem = np.nanstd(nm, axis=0)/np.sqrt(
                np.maximum(np.sum(~np.isnan(nm), axis=0), 1))
            ax.plot(win_centers, nmn, color=CLR[reg], lw=1.4, ls="--", alpha=0.55,
                    label=f"{reg} shuffle null")
            ax.fill_between(win_centers, nmn-nsem, nmn+nsem, color=CLR[reg], alpha=0.08)
        pk = win_centers[np.nanargmax(mn)]
        ax.axvline(pk, color=CLR[reg], lw=1.5, ls="--", alpha=0.85,
                   label=f"{reg} peak = {pk*1000:.0f} ms")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.axhline(0, color="k", lw=0.4)
    ax.set_title(title)
    if xlabel: ax.set_xlabel("window center re reach (s)")
    if ylabel: ax.set_ylabel("mean MI ± SEM (bits)")
    ax.spines[["top", "right"]].set_visible(False)
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=fontsize)

def _draw_lead_lag(ax, sel, title, fontsize=7, xlabel=True, ylabel=True):
    for reg in ("CP", "CB"):
        peaks = _peaks_for([pi for pi in sel if regions_run[pi] == reg])
        if not peaks: continue
        ax.hist(peaks, bins=bin_edges, color=CLR[reg], alpha=0.6,
                edgecolor="white", label=f"{reg} (n={len(peaks)})")
        med = np.median(peaks)
        ax.axvline(med, color=CLR[reg], lw=1.5, ls="--", alpha=0.9,
                   label=f"{reg} median = {med*1000:.0f} ms")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.set_title(title)
    if xlabel: ax.set_xlabel("time of peak MI re reach (s)")
    if ylabel: ax.set_ylabel("# neurons")
    if PEAK_RANGE is not None: ax.set_xlim(PEAK_RANGE)
    ax.spines[["top", "right"]].set_visible(False)
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=fontsize)

PANELS = [(ms, [pi for pi in range(len(run_idx)) if mice_run[pi] == ms])
          for ms in MOUSE_LIST]
PANELS.append(("ALL MICE", list(range(len(run_idx)))))
NCOL = 3
NROW = int(np.ceil(len(PANELS) / NCOL))

# mean landscape grid -- y SHARED, the panels are meant to be compared
fig, axes = plt.subplots(NROW, NCOL, figsize=(6*NCOL, 3.6*NROW),
                         sharex=True, sharey=True)
axes = np.atleast_1d(axes).ravel()
for k, (name, sel) in enumerate(PANELS):
    _draw_mean_landscape(axes[k], sel, name,
                         xlabel=(k // NCOL == NROW-1), ylabel=(k % NCOL == 0))
for ax in axes[len(PANELS):]: ax.set_visible(False)
plt.suptitle(f"Mean MI landscape per mouse  —  {lab}", fontsize=12)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(fdir / "panels_mean_landscape.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

# lead/lag grid -- y NOT shared: neuron counts differ a lot between mice, and the
# pooled panel would otherwise flatten every individual one
fig, axes = plt.subplots(NROW, NCOL, figsize=(6*NCOL, 3.6*NROW), sharex=True)
axes = np.atleast_1d(axes).ravel()
for k, (name, sel) in enumerate(PANELS):
    _draw_lead_lag(axes[k], sel, name,
                   xlabel=(k // NCOL == NROW-1), ylabel=(k % NCOL == 0))
for ax in axes[len(PANELS):]: ax.set_visible(False)
plt.suptitle(f"Peak-MI (lead/lag) distribution per mouse  —  {lab}", fontsize=12)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(fdir / "panels_lead_lag.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

print(f"\n{'mouse':16s} {'region':7s} {'n':>5s} {'peak (ms)':>10s} "
      f"{'n sig':>6s} {'median peak (ms)':>17s}")
for name, sel in PANELS:
    for reg in ("CP", "CB"):
        idx = [pi for pi in sel if regions_run[pi] == reg]
        if len(idx) < MIN_NEURONS_MOUSE: continue
        mn = np.nanmean(mi_joint[idx], axis=0)
        pks = _peaks_for(idx)
        med = f"{np.median(pks)*1000:17.0f}" if pks else f"{'-':>17s}"
        print(f"{name:16s} {reg:7s} {len(idx):5d} "
              f"{win_centers[np.nanargmax(mn)]*1000:10.0f} {len(pks):6d}{med}")

# %%
