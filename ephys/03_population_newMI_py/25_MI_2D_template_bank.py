# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# WHAT MAKES UP THE BANDS OF THE MEAN 2D MAP? -- by TEMPLATE MATCHING
#
# Same question as 23: the mean map is a vertical band in CP and a diagonal band
# in CB. Is that carried by single neurons that are themselves a column / a
# diagonal, or by neurons that each carry ONE blob and together tile the band?
#
# 23 read the shape of each map's hot cluster, which needs a hot-cell threshold
# and a connected cluster (speckled real maps break it up). Here each map is
# instead compared with a BANK OF TEMPLATES that covers the whole map:
#     column   (t_n, s)            exp(-(N - t_n)^2 / 2s^2), full kinematic height
#     diagonal (L, s)              exp(-((N - K) - L)^2 / 2s^2)
#     blob     (t_n, t_k, s_n, s_k) a 2D Gaussian, shorter than the kinematic range
# Score = Pearson r between the map and the template over the map's cells (= the
# best fit of that one template, amplitude + offset; r^2 = variance explained).
# A template sums over its whole shape (a matched filter): no hot-cell threshold,
# no smoothing, no cluster. Neighbouring templates are SPACING x their width apart,
# so any shape of that width is matched with r >= ~0.94.
#
# NOISE FLOOR: the best r of each family on the region's SHUFFLE maps (23's, one
# per neuron, all trials). Each family has its own floor, since a family with
# more templates (blobs) finds a higher best r in pure noise.
#     no family above its floor                -> not informative
#     elevated everywhere (median z > BROAD_Z) -> broad
#     the family with the highest r            -> column / diagonal / tiling
#     top two families within MARGIN           -> ambiguous
# The maps are 12's z-maps, z = (MI - mean of 20 shuffles) / their SD, as in 23.
import os, sys, time, pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

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

# ── which maps: 12's run (must match 12's cache key exactly) ────────────────
REACH_CAT     = (1, )
FEATURES_2D   = ("x", "y", "z", "speed", "dist_water")
KIN_WIN       = 21     # frames (52.5 ms)
KIN_STEP      = 4      # frames (10 ms)
KIN_RANGE     = (-0.34, 0.10)
NEU_WIN_WIDTH = 0.07   # (s)
NEU_WIN_STEP  = 0.025  # (s)
NEU_WIN_RANGE = (-0.5, 0.5)
N_SHUF_12     = 20
MIN_TRIALS    = 20
SRC_12 = "td2D_v2_x-y-z-speed-dist_water_nw70_kw21_ks4_nall_s20_cat1_2Dcache.pkl"
SRC_23 = "shapes2D_x-y-z-speed-dist_water_nw70_ks4_cat1_nall_shufmap_cache.pkl"  # shuffle maps
SRC_23_RES = "shapes2D_x-y-z-speed-dist_water_nw70_ks4_cat1_nall_results.pkl"   # 23's classes
SRC_22 = "split2D_x-y-z-speed-dist_water_nw70_ks8_cat1_nall_results.pkl"        # reliable set
SHUF_SEED_23  = 50000  # 23's seed: checked against its cache key

N_TEST        = None   # neurons PER REGION for a test run (fixed seed); None = all

# ── the template bank ────────────────────────────────────────────────────────
COL_SIG    = (0.025, 0.05)    # (s) widths (SD) of the column templates
DIAG_SIG   = (0.025, 0.05)    # (s) widths (SD) of the diagonal templates
BLOB_SIG_N = (0.025, 0.05)    # (s) blob width along neural time
BLOB_SIG_K = (0.03, 0.06)     # (s) blob width along kinematic time. Keep well below
                              # the kinematic range (0.43 s, SD of a full stripe
                              # 0.124 s): wider and a blob becomes a column
SPACING    = 1.0              # neighbouring templates are SPACING x their SD apart
                              # (1.0 -> worst-case r 0.95; 2.0 -> 0.82)
MIN_MASS   = 0.5              # a template needs this share of its full mass on the map

# ── classification ───────────────────────────────────────────────────────────
MIN_FINITE = 0.6       # maps defined in fewer cells are "too sparse" (set aside)
Z_CLIP     = 20.0      # z clipped to +-Z_CLIP before scoring, real AND shuffle maps:
                       # a few cells have a near-zero null SD and huge z (shuffle
                       # z-maps: skew ~490), and one such cell can dominate a Pearson r.
                       # Real peaks are ~z 14 (median max). None = no clipping
NULL_PCT   = 95        # floors: (100 - NULL_PCT)% of shuffle maps pass ANY family
MARGIN     = 0.05      # top two families closer than this in r -> "ambiguous"
BROAD_Z    = 2.0       # median z of the map above this -> "broad" (elevated everywhere;
                       # ~10% of real maps). None = no broad class

# ── synthetic maps of known type (no data, no MI), as 23 ─────────────────────
SYNTHETIC     = False
SYN_NOISE_FREE = False         # True: the real maps carry no noise (the pure populations)
SYN_POP       = {"CP": {"time": 40, "tiling": 40, "broad": 20, "noise": 40},
                 "CB": {"kinematics": 40, "tiling_diag": 40, "broad": 20, "noise": 40}}
SYN_AMP_RANGE = (3.0, 12.0)
SYN_WIDTH     = 0.04
SYN_NOISE_SMOOTH = (0.9, 0.5)
SYN_SEED      = 0
EXPECTED = {"time": "column", "kinematics": "diagonal", "tiling": "tiling",
            "tiling_diag": "tiling", "broad": "broad", "noise": "not informative"}

# ── outputs ──────────────────────────────────────────────────────────────────
SAVE_FIGS     = True
SAVE_VARS     = True
SHOW_FIGS     = True
N_EXAMPLES    = 3      # example neurons per class and region
OUT_SUBDIR    = "MI_2D_templates"
RUN_NAME      = "tmpl2D"
CMAP          = "viridis"
CLR = {"CP": np.array([106, 176, 134]) / 255, "CB": np.array([154, 123, 193]) / 255}
FAMS = ("column", "diagonal", "blob")
FAM_CLR = {"column": "tab:orange", "diagonal": "tab:red", "blob": "tab:blue"}
CLASSES = ("column", "diagonal", "tiling", "ambiguous", "broad", "not informative", "too sparse")
CLASS_CLR = {"column": "tab:orange", "diagonal": "tab:red", "tiling": "tab:blue",
             "ambiguous": "tab:purple", "broad": "0.45", "not informative": "0.85",
             "too sparse": "0.7"}
FAM2CLASS = {"column": "column", "diagonal": "diagonal", "blob": "tiling"}


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


def _save(fig, name):
    if SAVE_FIGS:
        _savefig_retry(fig, _long_path(fdir / f"{RUN_NAME}_{name}.png"), dpi=150,
                       bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 1. the grid, 12's maps, 23's shuffle maps, 22's reliable set ──────────
# 1. cell --> the same grid as 12; 12's z-maps and 23's shuffle z-maps (matched
# by uid). SYNTHETIC: maps of known type instead (23's generator).
FPS = 400.0
win_starts  = np.arange(NEU_WIN_RANGE[0], NEU_WIN_RANGE[1] - NEU_WIN_WIDTH + NEU_WIN_STEP / 2,
                        NEU_WIN_STEP)
win_centers = win_starts + NEU_WIN_WIDTH / 2
kin_starts  = np.arange(0, config.REACH_N_FRAMES - KIN_WIN + 1, KIN_STEP)
kin_centers = (kin_starts + KIN_WIN / 2 - config.REACH_ENDPOINT_FRAME) / FPS
_in = (kin_centers >= KIN_RANGE[0]) & (kin_centers <= KIN_RANGE[1])
kin_starts, kin_centers = kin_starts[_in], kin_centers[_in]
n_kin, n_wins = len(kin_starts), len(win_starts)
Kc, Nc = np.meshgrid(kin_centers, win_centers, indexing="ij")     # (n_kin, n_wins), s
LAGc   = Nc - Kc
extent = [win_centers[0], win_centers[-1], kin_centers[0], kin_centers[-1]]

CAT_TAG  = "".join(str(c) for c in sorted(REACH_CAT))
NW_TAG   = f"nw{int(round(NEU_WIN_WIDTH * 1000))}"
FEAT_TAG = "-".join(FEATURES_2D)
OUT_DIR_C = (f"{OUT_SUBDIR}_SYNTHETIC{'_noisefree' if SYN_NOISE_FREE else ''}" if SYNTHETIC
             else f"{OUT_SUBDIR}_{FEAT_TAG}_{NW_TAG}_cat{CAT_TAG}"
             + ("" if N_TEST is None else f"_test{N_TEST}"))
fdir = config.fig_dir(OUT_DIR_C)


def _synthetic():
    """23's generator -> Z (real), ZN (null) z-maps, uids, regions, true types"""
    rs = np.random.default_rng(SYN_SEED)
    g = lambda x, w=SYN_WIDTH: np.exp(-0.5 * (x / w) ** 2)
    def noise():
        v = gaussian_filter(rs.normal(size=Kc.shape), SYN_NOISE_SMOOTH)
        return v / v.std()
    Z, ZN, U, R, T = [], [], [], [], []
    for reg, pop in SYN_POP.items():
        for typ, n in pop.items():
            T0 = rs.permutation(np.linspace(-0.15, 0.25, n))
            for j in range(n):
                kw = kin_centers[0] + (j + rs.uniform()) / n * (kin_centers[-1] - kin_centers[0])
                if typ == "time":          M = g(Nc - T0[j])
                elif typ == "kinematics":  M = g(LAGc - rs.uniform(-0.1, 0.1))
                elif typ == "tiling":      M = g(Nc - T0[j]) * g(Kc - kw)
                elif typ == "tiling_diag": M = g(Nc - (kw + rs.uniform(-0.1, 0.1))) * g(Kc - kw)
                elif typ == "broad":       M = 0.6 + 0.4 * g(Nc - T0[j], 0.3)
                elif typ == "noise":       M = np.zeros_like(Kc)
                else: raise ValueError(f"unknown synthetic type {typ!r}")
                nz = noise()
                Z.append(rs.uniform(*SYN_AMP_RANGE) * M + (0 * nz if SYN_NOISE_FREE else nz))
                ZN.append(noise())
                U.append(f"{reg}|SYN_{typ}|#{j}"); R.append(reg); T.append(typ)
    return np.array(Z), np.array(ZN), U, np.array(R), np.array(T)


if SYNTHETIC:
    Z, ZN, uids, regions_run, TRUE_TYPE = _synthetic()
    EXC = Z.copy()                                  # no MI here: the class means use z
    REL22 = np.zeros(len(uids), bool)
    print(f"SYNTHETIC{' (noise-free)' if SYN_NOISE_FREE else ''}: {len(uids)} maps of known type")
else:
    TRUE_TYPE = None
    t0 = time.time()
    with open(config.var_path(SRC_12), "rb") as fh: C12 = pickle.load(fh)
    _want = dict(features=list(FEATURES_2D), KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP,
                 KIN_RANGE=list(KIN_RANGE), n_kin=n_kin, n_wins=n_wins,
                 N_SHUF_2D=N_SHUF_12, WIN_WIDTH=NEU_WIN_WIDTH, WIN_STEP=NEU_WIN_STEP,
                 WIN_RANGE=list(NEU_WIN_RANGE), CAT_KEEP=sorted(REACH_CAT),
                 MIN_TRIALS=MIN_TRIALS, NEURAL_VAR="count")
    _bad = {k: (C12["key"].get(k), v) for k, v in _want.items() if C12["key"].get(k) != v}
    if _bad:
        raise ValueError(f"12's cache does not match these parameters: {_bad}")
    with open(config.var_path(SRC_23), "rb") as fh: C23 = pickle.load(fh)
    if C23["key"].get("src") != SRC_12 or C23["key"].get("seed") != SHUF_SEED_23:
        raise ValueError(f"23's shuffle maps were not made from {SRC_12} with seed {SHUF_SEED_23}")
    _row23 = {u: i for i, u in enumerate(C23["key"]["uids"])}
    uids_all = list(C12["uids"]); regions_all = np.asarray(C12["regions"])
    _ok = np.isfinite(C12["mi"]).any(axis=(1, 2)) & np.array([u in _row23 for u in uids_all])
    _rn = np.random.default_rng(1)
    sel = []
    for reg in ("CP", "CB"):
        c = np.flatnonzero((regions_all == reg) & _ok)
        sel += sorted(c if N_TEST is None or N_TEST >= len(c)
                      else _rn.choice(c, N_TEST, replace=False).tolist())
    sel = np.array(sel)
    uids, regions_run = [uids_all[i] for i in sel], regions_all[sel]
    MI, NU, NSD = C12["mi"][sel], C12["null"][sel], C12["null_std"][sel]
    SH = C23["sh"][[_row23[u] for u in uids]]
    del C12, C23
    with np.errstate(all="ignore"):
        _sd = np.where(NSD > 0, NSD, np.nan)
        Z  = (MI - NU) / _sd
        ZN = (SH - NU) / _sd
    EXC = MI - NU                                   # excess MI (bits), for the class means
    REL22 = np.zeros(len(uids), bool)
    _p22 = config.var_path(SRC_22)
    if _p22.is_file():
        with open(_p22, "rb") as fh: _r22 = pickle.load(fh)
        _rel = {u for u, r in zip(_r22["uids"], _r22["reliable"]) if r}
        REL22 = np.array([u in _rel for u in uids])
        print(f"22's reliable set: {len(_rel)} neurons; {int(REL22.sum())} of this run's")
    for reg in ("CP", "CB"):
        print(f"{reg}: {int((regions_run == reg).sum())} neurons "
              f"(of {int(((regions_all == reg) & _ok).sum())} with a map in 12 and 23)")
    print(f"maps loaded ({time.time() - t0:.0f}s)")
print(f"grid: {n_kin} kinematic x {n_wins} neural windows  |  outputs -> {OUT_DIR_C}/")


# %% ── 2. the template bank ──────────────────────────────────────────────────
# 2. cell --> all templates on the grid; how many; how much the families overlap.
g_ = lambda x, s: np.exp(-0.5 * (x / s) ** 2)


def _centres(lo, hi, step):
    """evenly spaced centres from lo to hi, at most `step` apart"""
    n = max(1, int(np.ceil((hi - lo) / step)) + 1)
    return np.linspace(lo, hi, n)


def _bank():
    """-> T (n_templates, n_kin * n_wins), list of parameter dicts"""
    T, P = [], []
    # reference masses: the same template centred in the middle of the map, on a
    # grid 3x wider (so nothing is cut off)
    def _keep(M, full):
        return M.sum() >= MIN_MASS * full
    dn, dk = NEU_WIN_STEP, kin_centers[1] - kin_centers[0]
    for s in COL_SIG:
        full = n_kin * np.sqrt(2 * np.pi) * s / dn
        for t in _centres(win_centers[0], win_centers[-1], SPACING * s):
            M = g_(Nc - t, s)
            if _keep(M, full): T.append(M.ravel()); P.append(dict(fam="column", s=s, t_n=t))
    for s in DIAG_SIG:
        full = n_kin * np.sqrt(2 * np.pi) * s / dn
        for L in _centres(LAGc.min(), LAGc.max(), SPACING * s):
            M = g_(LAGc - L, s)
            if _keep(M, full): T.append(M.ravel()); P.append(dict(fam="diagonal", s=s, lag=L))
    for sn in BLOB_SIG_N:
        for sk in BLOB_SIG_K:
            full = 2 * np.pi * sn * sk / (dn * dk)
            for tk in _centres(kin_centers[0], kin_centers[-1], SPACING * sk):
                for tn in _centres(win_centers[0], win_centers[-1], SPACING * sn):
                    M = g_(Nc - tn, sn) * g_(Kc - tk, sk)
                    if _keep(M, full):
                        T.append(M.ravel())
                        P.append(dict(fam="blob", s=sn, s_k=sk, t_n=tn, t_k=tk, lag=tn - tk))
    return np.array(T), P


TB, TP = _bank()
TFAM = np.array([p["fam"] for p in TP])
print("template bank: " + ", ".join(f"{f} {int((TFAM == f).sum())}" for f in FAMS)
      + f"  (total {len(TP)})")


def _corr_rows(A, B):
    """Pearson r between every row of A and every row of B (no NaN)"""
    A = A - A.mean(1, keepdims=True); B = B - B.mean(1, keepdims=True)
    A /= np.linalg.norm(A, axis=1, keepdims=True); B /= np.linalg.norm(B, axis=1, keepdims=True)
    return A @ B.T


# overlap: for each template of family A, its best r with any template of family B;
# the median over A's templates. Off the diagonal of this table = how easily a
# shape of one family is mistaken for another.
G = _corr_rows(TB, TB)
OVL = np.array([[np.median(G[np.ix_(TFAM == a, TFAM == b)].max(axis=1)) for b in FAMS]
                for a in FAMS])
fig = plt.figure(figsize=(17, 6.5))
gs_ = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 1, 1.1])
_ex = {f: np.flatnonzero(TFAM == f) for f in FAMS}
for r_, f in enumerate(("column", "diagonal")):
    for c_, j in enumerate(np.linspace(0, _ex[f].size - 1, 4).astype(int)):
        ax = fig.add_subplot(gs_[r_, c_]); i = _ex[f][j]
        ax.imshow(TB[i].reshape(n_kin, n_wins), origin="lower", aspect="auto", extent=extent, cmap=CMAP)
        ax.set_title(f"{f} " + ", ".join(f"{k} {v:+.3f}" for k, v in TP[i].items() if k != "fam"),
                     fontsize=7)
        ax.tick_params(labelsize=6)
ax = fig.add_subplot(gs_[0, 4])
_bl = _ex["blob"][np.linspace(0, _ex["blob"].size - 1, 9).astype(int)]
ax.imshow(TB[_bl].max(axis=0).reshape(n_kin, n_wins), origin="lower", aspect="auto",
          extent=extent, cmap=CMAP)
ax.set_title("9 blob templates (overlaid)", fontsize=8); ax.tick_params(labelsize=6)
ax = fig.add_subplot(gs_[1, 4])
im = ax.imshow(OVL, cmap="Greys", vmin=0, vmax=1)
for a in range(3):
    for b in range(3):
        ax.text(b, a, f"{OVL[a, b]:.2f}", ha="center", va="center",
                color="w" if OVL[a, b] > 0.6 else "k")
ax.set_xticks(range(3), FAMS); ax.set_yticks(range(3), FAMS)
ax.set_xlabel("best match in family"); ax.set_ylabel("template of family")
ax.set_title("overlap: median best r", fontsize=8)
fig.suptitle(f"STEP 0 / the template bank ({', '.join(f'{f} {_ex[f].size}' for f in FAMS)}); "
             f"x: neural window, y: kinematic window (s)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.94)); _save(fig, "step0_bank")
print("overlap (median best r, rows = template family, cols = best match in family):")
for a, f in enumerate(FAMS):
    print(f"  {f:9s}" + "".join(f"{OVL[a, b]:7.2f}" for b in range(3)))


# %% ── 3. scores ─────────────────────────────────────────────────────────────
# 3. cell --> r of every map with every template (NaN cells left out, exactly),
# best r per family, and each family's floor from the shuffle maps.
def _scores(Zs, T):
    """Zs (n_maps, n_kin, n_wins) with NaN; T (n_templates, cells)
    -> r, amplitude (n_maps, n_templates): Pearson r over each map's finite cells,
    and the fitted amplitude (z per unit of template)."""
    X = Zs.reshape(len(Zs), -1)
    M = np.isfinite(X).astype(float); X = np.where(M > 0, X, 0.0)
    n = M.sum(1, keepdims=True)
    Sx, Sxx = X.sum(1, keepdims=True), (X * X).sum(1, keepdims=True)
    St, Stt, Sxt = M @ T.T, M @ (T * T).T, X @ T.T
    with np.errstate(all="ignore"):
        cov = Sxt - Sx * St / n
        vt  = Stt - St ** 2 / n
        vx  = Sxx - Sx ** 2 / n
        r   = cov / np.sqrt(vx * vt)
        a   = cov / vt
    r[~np.isfinite(r) | (vt < 1e-9)] = np.nan
    return r, a


def _best(r):
    """-> {family: (best r, index of its template)} per map"""
    out = {}
    for f in FAMS:
        idx = np.flatnonzero(TFAM == f)
        rr = np.where(np.isfinite(r[:, idx]), r[:, idx], -np.inf)
        j = rr.argmax(axis=1)
        b = rr[np.arange(len(r)), j]
        out[f] = (np.where(np.isfinite(b), b, np.nan), idx[j])
    return out


t0 = time.time()
_clip = (lambda A: A) if Z_CLIP is None else (lambda A: np.clip(A, -Z_CLIP, Z_CLIP))
R_REAL, A_REAL = _scores(_clip(Z), TB)
R_NULL, _      = _scores(_clip(ZN), TB)
BEST, BESTN = _best(R_REAL), _best(R_NULL)
print(f"scores: {R_REAL.shape[0]} maps x {R_REAL.shape[1]} templates, real + shuffle "
      f"({time.time() - t0:.1f}s)")

FINITE = np.array([min(np.isfinite(a).mean(), np.isfinite(b).mean()) for a, b in zip(Z, ZN)])
SPARSE = FINITE < MIN_FINITE
# The floors are set JOINTLY: one percentile q for all three families, chosen so
# that (100 - NULL_PCT)% of the shuffle maps pass ANY family. Separate 95th
# percentiles would let ~3 x 5% of noise maps through (12% on the synthetic noise).
def _floors(m):
    B = np.column_stack([BESTN[f][0][m] for f in FAMS])
    B = B[np.isfinite(B).all(axis=1)]
    for q in np.arange(NULL_PCT, 100.0, 0.1):
        fl = np.percentile(B, q, axis=0)
        if np.mean((B > fl).any(axis=1)) <= 1 - NULL_PCT / 100: break
    return {f: float(v) for f, v in zip(FAMS, fl)}, q

FLOOR, FLOOR_Q = {}, {}
for reg in ("CP", "CB"):
    m = (regions_run == reg) & ~SPARSE
    if m.any(): FLOOR[reg], FLOOR_Q[reg] = _floors(m)
for reg, fl in FLOOR.items():
    print(f"{reg} floors (best r on shuffle maps, {FLOOR_Q[reg]:.1f}th pct per family -> "
          f"{100 - NULL_PCT:.0f}% of shuffles pass any): "
          + "  ".join(f"{f} {v:.3f}" for f, v in fl.items()))

fig, axes = plt.subplots(2, 3, figsize=(15, 7), squeeze=False)
for r_, reg in enumerate(FLOOR):
    m = (regions_run == reg) & ~SPARSE
    for c_, f in enumerate(FAMS):
        ax = axes[r_, c_]
        b = np.linspace(-0.05, 1, 60)
        ax.hist(BESTN[f][0][m], b, color="0.6", alpha=0.7, label="shuffle maps")
        ax.hist(BEST[f][0][m], b, color=CLR[reg], alpha=0.6, label="real maps")
        ax.axvline(FLOOR[reg][f], color="k", ls="--", lw=1, label=f"floor ({FLOOR_Q[reg]:.1f}th pct)")
        ax.set(xlabel=f"best r, {f} templates", ylabel="neurons",
               title=f"{reg} {f}: {int((BEST[f][0][m] > FLOOR[reg][f]).sum())}/{int(m.sum())} above")
        ax.legend(fontsize=7); ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("STEP 1 / best match per template family, real vs shuffle maps", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.94)); _save(fig, "step1_best_r")


# %% ── 4. classification ─────────────────────────────────────────────────────
# 4. cell --> one class per neuron; example maps with their best templates; the
# decision plane with every neuron as a dot; the percentages.
RB = np.column_stack([BEST[f][0] for f in FAMS])           # (n, 3) best r per family
RBN = np.column_stack([BESTN[f][0] for f in FAMS])
MEDZ = np.array([np.nanmedian(z) if np.isfinite(z).any() else np.nan for z in Z])


def _classify(i):
    if SPARSE[i]: return "too sparse"
    fl = FLOOR[regions_run[i]]
    above = [f for f in FAMS if RB[i, FAMS.index(f)] > fl[f]]
    if not above: return "not informative"
    if BROAD_Z is not None and MEDZ[i] > BROAD_Z: return "broad"
    rs_ = sorted(((RB[i, FAMS.index(f)], f) for f in above), reverse=True)
    if len(rs_) > 1 and rs_[0][0] - rs_[1][0] < MARGIN: return "ambiguous"
    return FAM2CLASS[rs_[0][1]]


CLASS = np.array([_classify(i) for i in range(len(uids))])
WIN_FAM = np.array([FAMS[int(np.nanargmax(RB[i]))] if np.isfinite(RB[i]).any() else ""
                    for i in range(len(uids))])        # the top family, whatever the class


def _draw(ax, Zm, title, tmpl=()):
    ok = np.isfinite(Zm).any()
    ax.imshow(Zm, origin="lower", aspect="auto", extent=extent, cmap=CMAP,
              vmin=min(0, np.nanpercentile(Zm, 1)) if ok else 0,
              vmax=np.nanpercentile(Zm, 99) if ok else 1)
    for j, col, lw in tmpl:                            # template outlines at half height
        ax.contour(win_centers, kin_centers, TB[j].reshape(n_kin, n_wins), [0.5],
                   colors=[col], linewidths=lw)
    dd = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(dd, dd, color="w", lw=0.8, ls="--", alpha=0.6)
    ax.set_title(title, fontsize=7); ax.tick_params(labelsize=6)


_rx = np.random.default_rng(2)
_show = ("column", "diagonal", "tiling", "ambiguous", "broad")
fig, axes = plt.subplots(2 * N_EXAMPLES, len(_show), figsize=(4 * len(_show), 2.4 * 2 * N_EXAMPLES),
                         squeeze=False)
for c, cls in enumerate(_show):
    for r0, reg in enumerate(("CP", "CB")):
        idx = np.flatnonzero((CLASS == cls) & (regions_run == reg))
        pick = _rx.choice(idx, min(N_EXAMPLES, idx.size), replace=False) if idx.size else []
        for j in range(N_EXAMPLES):
            ax = axes[r0 * N_EXAMPLES + j, c]
            if j >= len(pick): ax.axis("off"); continue
            i = pick[j]
            tt = f" [true {TRUE_TYPE[i]}]" if TRUE_TYPE is not None else ""
            _draw(ax, Z[i], f"{cls.upper()}  {uids[i]}{tt}\n"
                  + "  ".join(f"{f[:4]} {RB[i, k]:.2f}" for k, f in enumerate(FAMS)),
                  [(BEST[f][1][i], FAM_CLR[f], 2.0 if f == WIN_FAM[i] else 0.8) for f in FAMS])
fig.suptitle("STEP 2 / example maps per class (z, unsmoothed). Outlines: the best template of each "
             "family at half height (thick = the best family); orange column, red diagonal, blue blob",
             fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.97)); _save(fig, "step2_examples")

# the decision plane: x = r_column - r_diagonal, y = r_blob - max(r_column, r_diagonal)
fig, axes = plt.subplots(1, 2, figsize=(14, 6), squeeze=False)
for ax, reg in zip(axes[0], ("CP", "CB")):
    m = (regions_run == reg) & ~SPARSE
    if not m.any(): continue
    xn = RBN[m, 0] - RBN[m, 1]; yn = RBN[m, 2] - np.maximum(RBN[m, 0], RBN[m, 1])
    ax.scatter(xn, yn, s=6, color="0.75", alpha=0.5, lw=0, label="shuffle maps")
    for cls in CLASSES:
        mc = m & (CLASS == cls)
        if not mc.any(): continue
        ax.scatter(RB[mc, 0] - RB[mc, 1], RB[mc, 2] - np.maximum(RB[mc, 0], RB[mc, 1]), s=10,
                   color=CLASS_CLR[cls], alpha=0.8, lw=0, label=f"{cls} ({int(mc.sum())})",
                   marker="x" if cls == "not informative" else "o")
    ax.axhline(0, color="k", lw=0.8); ax.axvline(0, color="k", lw=0.8)
    for v in (-MARGIN, MARGIN):
        ax.axhline(v, color="k", lw=0.5, ls=":"); ax.axvline(v, color="k", lw=0.5, ls=":")
    ax.set(xlabel="r_column - r_diagonal   (> 0 column-like)",
           ylabel="r_blob - max(r_column, r_diagonal)   (> 0 blob-like)", title=reg)
    ax.text(0.98, 0.02, "column", transform=ax.transAxes, ha="right", color=CLASS_CLR["column"])
    ax.text(0.02, 0.02, "diagonal", transform=ax.transAxes, color=CLASS_CLR["diagonal"])
    ax.text(0.5, 0.97, "blob (tiling)", transform=ax.transAxes, ha="center", va="top",
            color=CLASS_CLR["tiling"])
    ax.legend(fontsize=7, loc="upper left"); ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"STEP 3 / every neuron by its best r per family (dotted: MARGIN {MARGIN}). The class also "
             f"needs the family above its floor; 'broad' = median z > {BROAD_Z}", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.94)); _save(fig, "step3_plane")

PCT = {}
for subset, msk in (("all", np.ones(len(uids), bool)), ("reliable (22)", REL22)):
    if not msk.any(): continue
    print(f"\n{subset}: % of neurons per class")
    print(f"{'':7s}{'n':>6s}" + "".join(f"{c:>17s}" for c in CLASSES))
    for reg in ("CP", "CB"):
        m = msk & (regions_run == reg)
        if not m.any(): continue
        PCT[(subset, reg)] = {c: 100 * float(np.mean(CLASS[m] == c)) for c in CLASSES}
        print(f"{reg:7s}{int(m.sum()):6d}" + "".join(f"{PCT[(subset, reg)][c]:16.1f}%" for c in CLASSES))

fig, axes = plt.subplots(1, max(1, len({k[0] for k in PCT})), figsize=(13, 4.5), squeeze=False)
for ax, subset in zip(axes[0], sorted({k[0] for k in PCT})):
    regs = [r for r in ("CP", "CB") if (subset, r) in PCT]
    left = np.zeros(len(regs))
    for c in CLASSES:
        v = np.array([PCT[(subset, r)][c] for r in regs])
        ax.barh(regs, v, left=left, color=CLASS_CLR[c], label=c)
        for y, (l, w) in enumerate(zip(left, v)):
            if w >= 4: ax.text(l + w / 2, y, f"{w:.0f}", ha="center", va="center", fontsize=8)
        left += v
    ax.set(xlabel="% of neurons", title=subset, xlim=(0, 100))
axes[0, -1].legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
fig.suptitle("STEP 4 / class percentages (template matching)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.92)); _save(fig, "step4_percentages")

if TRUE_TYPE is not None:
    types = list(dict.fromkeys(TRUE_TYPE))
    print("\nSYNTHETIC confusion (rows true type, columns called class, % of the type)")
    print(f"{'':12s}" + "".join(f"{c[:11]:>12s}" for c in CLASSES) + "   expected")
    for t in types:
        m = TRUE_TYPE == t
        print(f"{t:12s}" + "".join(f"{100 * np.mean(CLASS[m] == c):11.0f}%" for c in CLASSES)
              + f"   {EXPECTED[t]} {100 * np.mean(CLASS[m] == EXPECTED[t]):.0f}%")


# %% ── 5. does the class mix explain the mean map? ───────────────────────────
# 5. cell --> per region: the mean excess-MI map of each class, and each class's
# SHARE of the region's mean map inside its band (top BAND_FRAC cells), as 23.
BAND_FRAC = 0.25
SHARE = {}
_cls5 = ("column", "diagonal", "tiling", "ambiguous", "broad", "not informative")
fig, axes = plt.subplots(2, 1 + len(_cls5), figsize=(4 * (1 + len(_cls5)), 6.5), squeeze=False)
for r, reg in enumerate(("CP", "CB")):
    m = (regions_run == reg) & ~SPARSE
    if not m.any(): continue
    MEAN = np.nanmean(EXC[m], axis=0)
    band = MEAN >= np.nanpercentile(MEAN, 100 * (1 - BAND_FRAC))
    _draw(axes[r, 0], MEAN, f"{reg} mean map, all {int(m.sum())} neurons\n(white: its top "
                            f"{BAND_FRAC:.0%} = the band)")
    axes[r, 0].contour(win_centers, kin_centers, band.astype(float), [0.5], colors="w", linewidths=1)
    tot = np.nansum(np.where(band, np.nansum(EXC[m], axis=0), 0))
    SHARE[reg] = {}
    for c, cls in enumerate(_cls5, start=1):
        mc = m & (CLASS == cls)
        SHARE[reg][cls] = (float(np.nansum(np.where(band, np.nansum(EXC[mc], axis=0), 0)) / tot)
                           if mc.any() and tot != 0 else 0.0)
        if mc.any():
            _draw(axes[r, c], np.nanmean(EXC[mc], axis=0),
                  f"{cls}: mean of {int(mc.sum())} neurons\nshare of the band: {100 * SHARE[reg][cls]:.0f}%")
        else:
            axes[r, c].axis("off")
    print(f"{reg}: share of the mean map's band -- " + "  ".join(
        f"{k} {100 * v:.0f}%" for k, v in SHARE[reg].items()))
fig.suptitle("STEP 5 / the mean map split by class: who carries the band? (each panel its own colour "
             "scale; x: neural window, y: kinematic window, s)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.94)); _save(fig, "step5_class_means")


# %% ── 6. where the templates sit ────────────────────────────────────────────
# 6. cell --> the position of each neuron's best template: blob centres of the
# tilers (map, neural time, kinematic time, lag); the column time of the column
# neurons; the lag of the diagonal neurons.
def _par(fam, key):
    return np.array([TP[BEST[fam][1][i]].get(key, np.nan) for i in range(len(uids))])

POS = {"col_t": _par("column", "t_n"), "diag_lag": _par("diagonal", "lag"),
       "blob_tn": _par("blob", "t_n"), "blob_tk": _par("blob", "t_k"), "blob_lag": _par("blob", "lag")}
fig = plt.figure(figsize=(18, 8))
gs_ = fig.add_gridspec(2, 4)
for r, reg in enumerate(("CP", "CB")):
    mt = (regions_run == reg) & (CLASS == "tiling")
    ax = fig.add_subplot(gs_[r, 0])
    ax.plot(POS["blob_tn"][mt], POS["blob_tk"][mt], "o", ms=3, color=CLR[reg], alpha=0.6)
    dd = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(dd, dd, color="k", lw=0.8, ls="--")
    ax.set(xlim=extent[:2], ylim=extent[2:], xlabel="neural time (s)", ylabel="kinematic time (s)",
           title=f"{reg} tiling: blob centres (n={int(mt.sum())})")
    for c, (lab, series) in enumerate((
            ("neural time (s)", (("tiling", "blob_tn"), ("column", "col_t"))),
            ("kinematic time (s)", (("tiling", "blob_tk"),)),
            ("lag = neural - kinematic (s)", (("tiling", "blob_lag"), ("diagonal", "diag_lag")))),
            start=1):
        ax = fig.add_subplot(gs_[r, c])
        for cls, key in series:
            mm = (regions_run == reg) & (CLASS == cls) & np.isfinite(POS[key])
            if mm.sum() > 1:
                ax.hist(POS[key][mm], 25, histtype="step", lw=2, density=True,
                        color=CLR[reg] if cls == "tiling" else CLASS_CLR[cls],
                        label=f"{cls} (n={int(mm.sum())})")
        ax.set(xlabel=lab, ylabel="density")
        if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("STEP 6 / where the best templates sit: tilers' blob centres; column neurons' "
             "neural time; diagonal neurons' lag", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.95)); _save(fig, "step6_positions")


# %% ── 7. agreement with 23 ──────────────────────────────────────────────────
# 7. cell --> this script's class against 23's class for the same neurons (uid):
# rows 25, columns 23, % of the 25 class.
if not SYNTHETIC and config.var_path(SRC_23_RES).is_file():
    with open(config.var_path(SRC_23_RES), "rb") as fh: _r23 = pickle.load(fh)
    _c23 = dict(zip(_r23["uids"], _r23["cls"]))
    C23 = np.array([_c23.get(u, "") for u in uids])
    cls23 = [c for c in dict.fromkeys(_r23["cls"])]
    for reg in ("CP", "CB"):
        m = (regions_run == reg) & (C23 != "")
        if not m.any(): continue
        print(f"\n{reg}: rows = 25's class, columns = 23's class (% of the row)")
        print(f"{'':16s}{'n':>6s}" + "".join(f"{c[:10]:>11s}" for c in cls23))
        for c in CLASSES:
            mc = m & (CLASS == c)
            if not mc.any(): continue
            print(f"{c:16s}{int(mc.sum()):6d}" + "".join(f"{100 * np.mean(C23[mc] == d):10.0f}%"
                                                        for d in cls23))
elif not SYNTHETIC:
    print(f"23's results not found ({SRC_23_RES}) -- no comparison")


# %% ── 8. save ──────────────────────────────────────────────────────────────
if SAVE_VARS and not SYNTHETIC:
    p = config.var_path(f"{RUN_NAME}_{FEAT_TAG}_{NW_TAG}_ks{KIN_STEP}_cat{CAT_TAG}"
                        f"_n{'all' if N_TEST is None else N_TEST}_results.pkl")
    with open(p, "wb") as fh:
        pickle.dump(dict(uids=uids, regions=regions_run, kin_centers=kin_centers,
                         win_centers=win_centers, sparse=SPARSE, finite=FINITE,
                         templates=TP, best_r={f: BEST[f][0] for f in FAMS},
                         best_idx={f: BEST[f][1] for f in FAMS},
                         best_r_null={f: BESTN[f][0] for f in FAMS}, floor=FLOOR,
                         median_z=MEDZ, cls=CLASS, top_family=WIN_FAM, positions=POS,
                         pct=PCT, share=SHARE, reliable22=REL22,
                         params=dict(COL_SIG=COL_SIG, DIAG_SIG=DIAG_SIG, BLOB_SIG_N=BLOB_SIG_N,
                                     BLOB_SIG_K=BLOB_SIG_K, SPACING=SPACING, MIN_MASS=MIN_MASS,
                                     MIN_FINITE=MIN_FINITE, NULL_PCT=NULL_PCT, MARGIN=MARGIN,
                                     BROAD_Z=BROAD_Z, Z_CLIP=Z_CLIP, SRC_12=SRC_12, SRC_23=SRC_23, SRC_22=SRC_22)),
                    fh, protocol=4)
    print(f"results -> {p}")

# %%
