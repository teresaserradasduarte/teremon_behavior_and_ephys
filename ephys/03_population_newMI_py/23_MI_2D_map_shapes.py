# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# WHAT MAKES UP THE BANDS OF THE MEAN 2D MAP? -- from ALL trials
#
# The mean map (12) is a vertical band in CP (locked to neural time) and a
# diagonal band in CB (locked to kinematics). Each can come from single neurons
# that are themselves a stripe / a ridge, or from neurons that each carry ONE
# blob and together tile the band -- or from a mix. This script puts a number
# on it: % column-locked, % diagonal-locked, % tiling, and where the tilers' blobs
# sit. (22 asked the same with a trial split; half the trials was too noisy.)
#
# THE MAPS are 12's, computed on all trials (no new MI): per cell,
#     z = (MI - mean of 20 shuffles) / SD of the 20 shuffles
# lightly smoothed. ONE extra shuffle map per neuron is computed here (all
# trials, new permutations) and goes through exactly the same steps: it is the
# reference for every statistic ("what would this look like with no information").
#
# THE SHAPE of a neuron's map is read from its HOT CLUSTER: the connected cells
# above Z_HOT that contain its strongest cell. No peak-picking, no model fit:
# only the weighted spread of that cluster, in seconds,
#     SD_n   along neural time      SD_k   along kinematic time
#     SD_lag along lag = neural time - kinematic time
# gives two numbers per neuron,
#     E  extent      SD_k / SD_k of a full-height stripe   (0 = blob, 1 = spans
#                    every kinematic window)
#     O  orientation (SD_lag - SD_n) / (SD_lag + SD_n)      (+1 = vertical =
#                    column-locked, -1 = diagonal = locked to kinematics)
# and a centre of mass (neural time, kinematic time, lag).
#     E >= E_THR, O >=  O_THR   column-locked
#     E >= E_THR, O <= -O_THR   diagonal-locked
#     E >= E_THR, |O| < O_THR   broad
#     E <  E_THR, D >= D_MIN    tiling (one blob that dominates the map)
#     E <  E_THR, D <  D_MIN    scattered (several similar specks)
#     E <  E_MIN                noise-sized (no wider than the shuffle maps' clusters)
# D = share of the map's hot mass that is in the main cluster.
# Every neuron with a cluster is plotted in (O, E), so the boundaries can be judged.
import os, sys, re, time, pickle
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from scipy.ndimage import gaussian_filter, label as nd_label

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
import mi_io, config, mi_rates

# ── which maps: 12's run (must match 12's cache key exactly) ────────────────
REACH_CAT     = (1, )
TRIAL_SELECT  = None
FEATURES_2D   = ("x", "y", "z", "speed", "dist_water")
KIN_WIN       = 21     # frames (52.5 ms)
KIN_STEP      = 4      # frames (10 ms)
KIN_RANGE     = (-0.34, 0.10)
NEU_WIN_WIDTH = 0.07   # (s)
NEU_WIN_STEP  = 0.025  # (s)
NEU_WIN_RANGE = (-0.5, 0.5)
N_SHUF_12     = 20     # shuffles 12 used for its null mean / SD
MIN_TRIALS    = 20     # 12's minimum number of trials
SRC_12 = "td2D_v2_x-y-z-speed-dist_water_nw70_kw21_ks4_nall_s20_cat1_2Dcache.pkl"
SRC_22 = "split2D_x-y-z-speed-dist_water_nw70_ks8_cat1_nall_results.pkl"   # reliable set

# ── the neurons ──────────────────────────────────────────────────────────────
N_TEST        = None   # neurons PER REGION for a test run (fixed seed); None = all
N_JOBS        = -1
SHUF_SEED     = 50000  # the extra shuffle map uses new permutations (12 used 100+i)

# ── map shape ────────────────────────────────────────────────────────────────
SMOOTH_SIG    = (1.0, 0.5)  # (kin, neural) cells: Gaussian smoothing of the z-maps
                       # (10 ms x 12.5 ms: well inside the windows' own overlap)
MIN_FINITE    = 0.6    # maps defined in fewer cells are "too sparse" (set aside)
TOP_FRAC      = 0.025  # informative test: mean of the map's top 2.5% cells ...
INF_PCT       = 95     # ... above this percentile of the same on the shuffle maps
HOT_REF       = "zero" # the hot cluster is the cells with z above
                       # "zero":   Z_HOT. A map high everywhere (~10% of real maps
                       #           sit above z=2 almost everywhere) is one big
                       #           cluster -> "broad"
                       # "median": the map's own median + Z_HOT (shape ABOVE its
                       #           baseline). Caution: in SYNTHETIC mode it calls
                       #           70% of broad maps "tiling" -- the bump left on
                       #           top of a high baseline looks like a blob
Z_HOT         = 2.0
MIN_CLUSTER   = 4      # cells; a smaller cluster has no measurable shape
E_THR         = 0.55   # extent threshold (blob below)
O_THR         = 0.25   # orientation threshold (|O| below = broad)
E_MIN         = None   # extent floor between noise and tiling: a cluster no wider
                       # (in kinematic time) than the shuffle maps' clusters is
                       # "noise-sized". None = the E_MIN_PCT percentile of E over the
                       # region's shuffle maps; or a number to fix it
E_MIN_PCT     = 95
D_MIN         = 0.5    # a small cluster is "tiling" only if it DOMINATES the map:
                       # D = its hot mass / the hot mass of ALL cells above the
                       # threshold. Lower D = several similar specks -> "scattered"

# ── synthetic maps of known type (no data, no MI) ────────────────────────────
SYNTHETIC     = False
SYN_POP       = {"CP": {"time": 40, "tiling": 40, "broad": 20, "noise": 40},
                 "CB": {"kinematics": 40, "tiling_diag": 40, "broad": 20, "noise": 40}}
                       #   time        vertical stripe     -> column-locked
                       #   kinematics  diagonal ridge      -> diagonal-locked
                       #   tiling      one blob, neural time independent of its
                       #               kinematic window (mean map: vertical band)
                       #   tiling_diag one blob on the diagonal (mean map: diagonal band)
                       #   broad       elevated almost everywhere
                       #   noise       nothing
SYN_AMP_RANGE = (3.0, 12.0)   # peak z of a synthetic neuron, uniform in this range
                       # (real maps: top-2.5% mean z 4.6 / 6.9 / 15 at pct 10/50/90)
SYN_WIDTH     = 0.04   # (s) SD of the stripe, ridge and blob
SYN_NOISE_SMOOTH = (0.9, 0.5)  # (kin, neural) cells: smoothness of the z noise
                       # (from the real maps' neighbour correlations, 0.72 / 0.38)
SYN_SEED      = 0
EXPECTED = {"time": "column", "kinematics": "diagonal", "tiling": "tiling",
            "tiling_diag": "tiling", "broad": "broad", "noise": "not informative"}

# ── outputs ──────────────────────────────────────────────────────────────────
CACHE_MODE    = "load"
SAVE_CACHE    = True
SAVE_FIGS     = True
SAVE_VARS     = True
SHOW_FIGS     = True
N_EXAMPLES    = 3      # example neurons per class and region
OUT_SUBDIR    = "MI_2D_map_shapes"
RUN_NAME      = "shapes2D"
CMAP          = "viridis"
CLR = {"CP": np.array([106, 176, 134]) / 255, "CB": np.array([154, 123, 193]) / 255}
CLASSES = ("column", "diagonal", "tiling", "scattered", "broad", "noise-sized", "no cluster",
           "not informative")
CLASS_CLR = {"column": "tab:orange", "diagonal": "tab:red", "tiling": "tab:blue",
             "scattered": "tab:purple", "noise-sized": "tab:olive",
             "broad": "0.45", "no cluster": "0.7", "not informative": "0.85"}

TIMES = {}
class _step:
    def __init__(self, name): self.name = name
    def __enter__(self): self.t = time.time(); return self
    def __exit__(self, *a):
        dt = time.time() - self.t
        TIMES[self.name] = TIMES.get(self.name, 0.0) + dt
        print(f"    [{self.name}] {dt:.1f}s" if dt < 90 else
              f"    [{self.name}] {dt/60:.1f} min", flush=True)
        return False


def _savefig_retry(fig, path, tries=6, **kw):
    """Dropbox briefly locks files it is uploading: retry instead of crashing."""
    for t in range(tries):
        try:
            return fig.savefig(path, **kw)
        except OSError:
            if t == tries - 1: raise
            time.sleep(0.5 * (t + 1))


def _long_path(p):
    """Windows refuses paths of 260+ characters unless they start with \\\\?\\ --
    the Dropbox output folder alone is ~210 characters"""
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


# %% ── 1. the grid, 12's maps, 22's reliable set ─────────────────────────────
# 1. cell --> the same grid as 12; load 12's maps (all trials) and 22's reliable
# neurons. SYNTHETIC: the grid only.
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
SD_K_FULL = float(np.std(kin_centers))     # SD_k of a stripe covering every kinematic window
extent = [win_centers[0], win_centers[-1], kin_centers[0], kin_centers[-1]]

CAT_TAG  = "".join(str(c) for c in sorted(REACH_CAT))
NW_TAG   = f"nw{int(round(NEU_WIN_WIDTH * 1000))}"
FEAT_TAG = "-".join(FEATURES_2D)
OUT_DIR_C = (f"{OUT_SUBDIR}_SYNTHETIC" if SYNTHETIC else
             f"{OUT_SUBDIR}_{FEAT_TAG}_{NW_TAG}_cat{CAT_TAG}"
             + ("" if N_TEST is None else f"_test{N_TEST}"))
fdir = config.fig_dir(OUT_DIR_C)

if not SYNTHETIC:
    with _step("read 12's maps"):
        with open(config.var_path(SRC_12), "rb") as fh: C12 = pickle.load(fh)
    _want = dict(features=list(FEATURES_2D), KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP,
                 KIN_RANGE=list(KIN_RANGE), n_kin=n_kin, n_wins=n_wins,
                 N_SHUF_2D=N_SHUF_12, WIN_WIDTH=NEU_WIN_WIDTH, WIN_STEP=NEU_WIN_STEP,
                 WIN_RANGE=list(NEU_WIN_RANGE), CAT_KEEP=sorted(REACH_CAT),
                 MIN_TRIALS=MIN_TRIALS, NEURAL_VAR="count")
    _bad = {k: (C12["key"].get(k), v) for k, v in _want.items() if C12["key"].get(k) != v}
    if _bad:
        raise ValueError(f"12's cache does not match these parameters: {_bad}")
    uids_all = list(C12["uids"]); regions_all = np.asarray(C12["regions"])
    _ok = np.isfinite(C12["mi"]).any(axis=(1, 2))           # neurons 12 could compute
    _rn = np.random.default_rng(1)
    sel = []
    for reg in ("CP", "CB"):
        c = np.flatnonzero((regions_all == reg) & _ok)
        sel += sorted(c if N_TEST is None or N_TEST >= len(c)
                      else _rn.choice(c, N_TEST, replace=False).tolist())
    sel = np.array(sel)
    uids, regions_run = [uids_all[i] for i in sel], regions_all[sel]
    MI, NU, NSD = C12["mi"][sel], C12["null"][sel], C12["null_std"][sel]
    del C12
    # 22's reliable neurons (split-half test), matched by uid
    REL22 = np.zeros(len(uids), bool)
    _p22 = config.var_path(SRC_22)
    if _p22.is_file():
        with open(_p22, "rb") as fh: _r22 = pickle.load(fh)
        _rel = {u for u, r in zip(_r22["uids"], _r22["reliable"]) if r}
        REL22 = np.array([u in _rel for u in uids])
        print(f"22's reliable set: {len(_rel)} neurons; {int(REL22.sum())} of this run's")
    else:
        print(f"22's results not found ({_p22.name}) -- no 'reliable from 22' subset")
    TRUE_TYPE = None
    for reg in ("CP", "CB"):
        print(f"{reg}: {int((regions_run == reg).sum())} neurons "
              f"(of {int(((regions_all == reg) & _ok).sum())} with a map in 12)")
print(f"grid: {n_kin} kinematic x {n_wins} neural windows  |  outputs -> {OUT_DIR_C}/")


# %% ── 2. one shuffle map per neuron (all trials) ────────────────────────────
# 2. cell --> 12's trials and windows, kinematics permuted across trials with
# NEW permutations; same estimator. Cached. SYNTHETIC: maps of known type instead.
def _shuffle_map(kin_stack, C_m, k, seed):
    """one KSG map with the trials' kinematics permuted (as 12's _landscape2d_one)"""
    out = np.full((kin_stack.shape[1], C_m.shape[1]), np.nan)
    p = np.random.default_rng(seed).permutation(kin_stack.shape[0])
    ok_neu = [wi for wi in range(C_m.shape[1]) if np.unique(C_m[:, wi]).size >= 2]
    for ki in range(kin_stack.shape[1]):
        X = kin_stack[:, ki, :]
        if not np.isfinite(X).all(): continue
        for wi in ok_neu:
            out[ki, wi] = _ksg_md(X[p], C_m[:, wi], k, seed + ki * 97 + wi)
    return out


def _synthetic():
    """-> Z (real), ZN (null) z-maps, uids, regions, true types"""
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
                Z.append(rs.uniform(*SYN_AMP_RANGE) * M + noise()); ZN.append(noise())
                U.append(f"{reg}|SYN_{typ}|#{j}"); R.append(reg); T.append(typ)
    return np.array(Z), np.array(ZN), U, np.array(R), np.array(T)


if SYNTHETIC:
    Z_RAW, ZN_RAW, uids, regions_run, TRUE_TYPE = _synthetic()
    EXC = Z_RAW.copy()                       # no MI here: the class means use z
    REL22 = np.zeros(len(uids), bool)
    print(f"SYNTHETIC: {len(uids)} maps of known type  |  "
          + "  ".join(f"{r}: " + ", ".join(f"{t} {n}" for t, n in p.items()) for r, p in SYN_POP.items()))
else:
    KEY = dict(src=SRC_12, uids=uids, seed=SHUF_SEED, v="v1")
    CACHE_PATH = config.var_path(f"{RUN_NAME}_{FEAT_TAG}_{NW_TAG}_ks{KIN_STEP}_cat{CAT_TAG}"
                                 f"_n{'all' if N_TEST is None else N_TEST}_shufmap_cache.pkl")
    SH = None
    if CACHE_MODE == "load" and CACHE_PATH.is_file():
        with open(CACHE_PATH, "rb") as fh: _c = pickle.load(fh)
        if _c.get("key") == KEY:
            SH = _c["sh"]; print(f"LOADED {CACHE_PATH.name} -- no recompute")
        else:
            print("cache exists but parameters differ -> recomputing")
    if SH is None:
        # 09's cells 0-2: data, kinematics in mm, water features, the estimator
        _src9 = open(os.path.join(_HERE, "09_MI_multiD_kin.py"), encoding="utf-8").read()
        _head = _src9[:re.search(r"(?m)^# %% ── 3\.", _src9).start()]
        for _nm, _val in (("FEATURES", tuple(FEATURES_2D)), ("N_TEST", None),
                          ("CAT_KEEP", set(REACH_CAT))):
            _head, _n = re.subn(rf"(?m)^{_nm}\s*=.*$", f"{_nm} = {_val!r}", _head, count=1)
            if not _n: raise RuntimeError(f"could not find {_nm} in 09 -- renamed?")
        # run in its OWN namespace: 09 defines its own grid, uids, ... which must
        # not overwrite this script's; only the data and the estimator are taken
        import io as _io, contextlib as _ctx
        _buf, _g9 = _io.StringIO(), {"__name__": "09_head"}
        with _step("load data + 09 setup"):
            try:
                with _ctx.redirect_stdout(_buf):
                    exec(compile(_head, "09_head", "exec"), _g9)
            except BaseException:
                print(_buf.getvalue()); raise
        NEURONS, REACH_MM, WATER = _g9["NEURONS"], _g9["REACH_MM"], _g9["WATER"]
        _ksg_md, K_KSG, NEURAL_VAR = _g9["_ksg_md"], _g9["K_KSG"], _g9["NEURAL_VAR"]
        savgol_filter, SAVGOL_WIN, SAVGOL_ORD = (_g9["savgol_filter"], _g9["SAVGOL_WIN"],
                                                 _g9["SAVGOL_ORD"])
        if NEURAL_VAR != "count":
            raise ValueError(f"09's NEURAL_VAR is {NEURAL_VAR!r}; 12's maps used 'count'")
        if TRIAL_SELECT:
            raise NotImplementedError("TRIAL_SELECT: 12's cache used here has none")
        # the kinematic feature of every window, per session (12's cell 2)
        KIN = {}
        for key in REACH_MM:
            rp  = REACH_MM[key]
            vel = savgol_filter(rp, SAVGOL_WIN, SAVGOL_ORD, deriv=1, delta=1/FPS, axis=2)
            spd = np.linalg.norm(vel, axis=1)
            acc = savgol_filter(spd, SAVGOL_WIN, SAVGOL_ORD, deriv=1, delta=1/FPS, axis=1)
            d = {"x": rp[:, 0], "y": rp[:, 1], "z": rp[:, 2], "speed": spd, "accel": acc,
                 "dist_water": WATER[key]["dist_water"], "speed_water": WATER[key]["speed_water"]}
            KIN[key] = np.stack([np.stack([d[f][:, s:s + KIN_WIN].mean(axis=1) for s in kin_starts],
                                          axis=1) for f in FEATURES_2D], axis=2)
        _byuid = {n["uid"]: i for i, n in enumerate(NEURONS)}
        jobs = []
        for pi, u in enumerate(uids):
            n = NEURONS[_byuid[u]]; K = KIN[(n["mouse"], n["sess"])]
            cat = np.asarray(n["labels"]["idx_reach_cat"]).ravel()
            msk = np.isfinite(K).all(axis=(1, 2)) & np.isin(cat, list(REACH_CAT))
            C = np.asarray(mi_rates.windowed(n["spikes"], win_starts, NEU_WIN_WIDTH, NEURAL_VAR),
                           float)[msk]
            jobs.append((K[msk], C, SHUF_SEED + pi))
        _calls = len(jobs) * n_kin * n_wins
        print(f"ETA: ~{_calls:,} KSG calls; at ~2500 calls/s that is ~{_calls/2500/60:.0f} min")
        with _step("shuffle maps"):
            SH = np.stack(Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
                delayed(_shuffle_map)(K_, C_, K_KSG, s_) for K_, C_, s_ in jobs))
        if SAVE_CACHE:
            with open(CACHE_PATH, "wb") as fh:
                pickle.dump(dict(key=KEY, sh=SH), fh, protocol=4)
            print(f"cache -> {CACHE_PATH}")
    if SH.shape != MI.shape:
        raise ValueError(f"shuffle maps {SH.shape} do not match 12's maps {MI.shape}")
    with np.errstate(all="ignore"):
        _sd = np.where(NSD > 0, NSD, np.nan)
        Z_RAW  = (MI - NU) / _sd
        ZN_RAW = (SH - NU) / _sd
    EXC = MI - NU                                # excess MI (bits), for the class means


# %% ── 3. z-maps, sparse maps, informative test ──────────────────────────────
# 3. cell --> smooth (NaN-aware), set sparse maps aside, and ask: is the map
# stronger than a shuffle map? Statistic = mean of its top TOP_FRAC cells;
# threshold = INF_PCT percentile of the same on the region's shuffle maps.
def _smooth(Zm):
    ok = np.isfinite(Zm)
    if not ok.any(): return Zm.copy()
    num = gaussian_filter(np.where(ok, Zm, 0.0), SMOOTH_SIG)
    den = gaussian_filter(ok.astype(float), SMOOTH_SIG)
    with np.errstate(all="ignore"):
        return np.where(ok, num / den, np.nan)


def _top(Zm):
    v = Zm[np.isfinite(Zm)]
    if v.size < 10: return np.nan
    return float(np.sort(v)[-max(1, int(round(TOP_FRAC * v.size))):].mean())


ZS  = np.stack([_smooth(z) for z in Z_RAW])
ZSN = np.stack([_smooth(z) for z in ZN_RAW])
FINITE = np.array([min(np.isfinite(a).mean(), np.isfinite(b).mean()) for a, b in zip(Z_RAW, ZN_RAW)])
SPARSE = FINITE < MIN_FINITE
TOP_R  = np.array([_top(z) for z in ZS]); TOP_N = np.array([_top(z) for z in ZSN])
INF_THR = {reg: float(np.nanpercentile(TOP_N[(regions_run == reg) & ~SPARSE], INF_PCT))
           for reg in ("CP", "CB") if ((regions_run == reg) & ~SPARSE).any()}
INFORM = np.array([(not SPARSE[i]) and TOP_R[i] > INF_THR.get(regions_run[i], np.inf)
                   for i in range(len(uids))])
for reg in INF_THR:
    m = regions_run == reg
    print(f"{reg}: {int(SPARSE[m].sum())}/{int(m.sum())} too sparse  |  informative "
          f"{int(INFORM[m].sum())}/{int((m & ~SPARSE).sum())} "
          f"({100*INFORM[m & ~SPARSE].mean():.0f}%; threshold top-{TOP_FRAC:.1%} z > {INF_THR[reg]:.2f})")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
for ax, reg in zip(axes, ("CP", "CB")):
    m = (regions_run == reg) & ~SPARSE
    if not m.any(): continue
    _v = np.r_[TOP_N[m], TOP_R[m]]; _v = _v[np.isfinite(_v)]
    b = np.linspace(np.percentile(_v, 0.5), np.percentile(_v, 99.5), 50)
    ax.hist(np.clip(TOP_N[m], b[0], b[-1]), b, color="0.6", alpha=0.7, label="shuffle maps")
    ax.hist(np.clip(TOP_R[m], b[0], b[-1]), b, color=CLR[reg], alpha=0.6, label="real maps")
    ax.axvline(INF_THR[reg], color="k", ls="--", lw=1, label=f"{INF_PCT}th pct of shuffles")
    ax.set(xlabel=f"mean smoothed z of the map's top {TOP_FRAC:.1%} cells", ylabel="neurons",
           title=f"{reg}: {int(INFORM[m].sum())}/{int(m.sum())} informative")
    ax.legend(fontsize=7); ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("STEP 1 / is the map stronger than a shuffle map? (all trials)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.93)); _save(fig, "step1_informative")


# %% ── 4. the hot cluster and its shape ──────────────────────────────────────
# 4. cell --> the connected cells above the hot threshold that contain the
# strongest cell; its weighted spread -> extent E, orientation O, centre of mass.
def _cluster(Zm):
    """-> boolean mask of the hot cluster (or None), and the weights"""
    if not np.isfinite(Zm).any(): return None, None
    ref = np.nanmedian(Zm) if HOT_REF == "median" else 0.0
    w = Zm - (ref + Z_HOT)
    hot = np.isfinite(w) & (w > 0)
    if not hot.any(): return None, None
    lab, _ = nd_label(hot, structure=np.ones((3, 3)))
    i0 = np.unravel_index(np.nanargmax(np.where(hot, Zm, -np.inf)), Zm.shape)
    cl = lab == lab[i0]
    if cl.sum() < MIN_CLUSTER: return None, None
    return cl, np.where(cl, w, 0.0)


def _shape(Zm):
    """-> dict(E, O, com_n, com_k, com_lag, size) of the hot cluster (NaNs if none)"""
    cl, w = _cluster(Zm)
    nan = dict(E=np.nan, O=np.nan, D=np.nan, com_n=np.nan, com_k=np.nan, com_lag=np.nan, size=0)
    if cl is None: return nan
    ws = w.sum()
    ref = np.nanmedian(Zm) if HOT_REF == "median" else 0.0
    w_all = np.clip(np.nan_to_num(Zm - (ref + Z_HOT), nan=0.0), 0, None)
    m_n, m_k, m_l = (float((w * A).sum() / ws) for A in (Nc, Kc, LAGc))
    sd = lambda A, mu: float(np.sqrt((w * (A - mu) ** 2).sum() / ws))
    sd_n, sd_k, sd_l = sd(Nc, m_n), sd(Kc, m_k), sd(LAGc, m_l)
    O = (sd_l - sd_n) / (sd_l + sd_n) if (sd_l + sd_n) > 0 else np.nan
    return dict(E=sd_k / SD_K_FULL, O=O, D=float(ws / w_all.sum()), com_n=m_n, com_k=m_k,
                com_lag=m_l, size=int(cl.sum()))


SHAPE  = [_shape(z) for z in ZS]
SHAPEN = [_shape(z) for z in ZSN]
_a = lambda S, k: np.array([d[k] for d in S], float)
E, O, E_N, O_N = _a(SHAPE, "E"), _a(SHAPE, "O"), _a(SHAPEN, "E"), _a(SHAPEN, "O")
D, D_N = _a(SHAPE, "D"), _a(SHAPEN, "D")
HAS_CL = np.isfinite(E)


# the floor between noise and tiling, per region (from the shuffle maps unless fixed)
E_MIN_R = {reg: (float(E_MIN) if E_MIN is not None else
                 float(np.nanpercentile(E_N[(regions_run == reg) & ~SPARSE], E_MIN_PCT))
                 if np.isfinite(E_N[(regions_run == reg) & ~SPARSE]).any() else 0.0)
           for reg in ("CP", "CB")}


def _classify(e, o, d, informative, sparse, e_min):
    if sparse:          return "too sparse"
    if not informative: return "not informative"
    if not np.isfinite(e): return "no cluster"
    if e < e_min:       return "noise-sized"
    if e < E_THR:       return "tiling" if d >= D_MIN else "scattered"
    return "column" if o >= O_THR else ("diagonal" if o <= -O_THR else "broad")


CLASS = np.array([_classify(E[i], O[i], D[i], INFORM[i], SPARSE[i], E_MIN_R[regions_run[i]])
                  for i in range(len(uids))])
print(f"\nhot cluster: z above {'the map median + ' if HOT_REF == 'median' else ''}{Z_HOT}, "
      f">= {MIN_CLUSTER} cells  |  E_THR {E_THR}, O_THR {O_THR}, D_MIN {D_MIN}  |  E_MIN "
      + ", ".join(f"{r} {v:.3f}" for r, v in E_MIN_R.items())
      + (f"  ({E_MIN_PCT}th pct of the shuffle maps' E)" if E_MIN is None else ""))
for reg in ("CP", "CB"):
    m = regions_run == reg
    if m.any():
        print(f"  {reg}: cluster found in {int((HAS_CL & m).sum())}/{int(m.sum())} maps "
              f"(shuffle maps: {int((np.isfinite(E_N) & m).sum())})")


def _draw(ax, Zm, title, cl=None, d=None):
    lim = np.nanpercentile(Zm, 99) if np.isfinite(Zm).any() else 1
    ax.imshow(Zm, origin="lower", aspect="auto", extent=extent, cmap=CMAP,
              vmin=min(0, np.nanpercentile(Zm, 1)) if np.isfinite(Zm).any() else 0, vmax=lim)
    if cl is not None:
        ax.contour(win_centers, kin_centers, cl.astype(float), [0.5], colors="w", linewidths=1)
    if d is not None and np.isfinite(d["E"]):
        ax.plot(d["com_n"], d["com_k"], "r+", ms=10, mew=2)
    dd = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(dd, dd, color="w", lw=0.8, ls="--", alpha=0.6)
    ax.set_title(title, fontsize=7); ax.tick_params(labelsize=6)


# FIGURE: example neurons of each class, with the cluster (white) and its centre (+)
_rx = np.random.default_rng(2)
_show_cls = ("column", "diagonal", "tiling", "scattered", "broad")
fig, axes = plt.subplots(2 * N_EXAMPLES, len(_show_cls), figsize=(4 * len(_show_cls), 2.3 * 2 * N_EXAMPLES),
                         squeeze=False)
for c, cls in enumerate(_show_cls):
    for r0, reg in enumerate(("CP", "CB")):
        idx = np.flatnonzero((CLASS == cls) & (regions_run == reg))
        pick = _rx.choice(idx, min(N_EXAMPLES, idx.size), replace=False) if idx.size else []
        for j in range(N_EXAMPLES):
            ax = axes[r0 * N_EXAMPLES + j, c]
            if j < len(pick):
                i = pick[j]; cl, _ = _cluster(ZS[i]); d = SHAPE[i]
                tt = f" [true {TRUE_TYPE[i]}]" if TRUE_TYPE is not None else ""
                _draw(ax, ZS[i], f"{cls.upper()}  {uids[i]}{tt}\nE {d['E']:.2f}  O {d['O']:+.2f}  D {d['D']:.2f}", cl, d)
            else:
                ax.axis("off")
fig.suptitle("STEP 2 / example maps per class (smoothed z; white = hot cluster, red + = its centre "
             "of mass; dashed = zero lag). x: neural window (s), y: kinematic window (s)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.97)); _save(fig, "step2_examples")


# %% ── 5. the (O, E) plane and the class percentages ─────────────────────────
# 5. cell --> EVERY neuron with a cluster as a dot (informative or not), the
# shuffle maps' dots for reference, and the threshold lines. Then % per class,
# for all neurons and for 22's reliable ones.
TYPE_STYLE = {"time": ("tab:orange", "s"), "kinematics": ("tab:red", "^"),
              "tiling": ("tab:blue", "o"), "tiling_diag": ("tab:cyan", "D"),
              "broad": ("0.4", "P"), "noise": ("0.6", "x")}
def _dots(ax, X, Y, XN, YN, m):
    """every neuron with a cluster: shuffle maps grey; real maps by region (or by
    true type in SYNTHETIC mode); open = not informative"""
    ax.plot(XN[m], YN[m], ".", ms=3, color="0.8", alpha=0.6, label="shuffle maps")
    if TRUE_TYPE is not None:
        for typ, (col_, mk) in TYPE_STYLE.items():
            for inf, fill in ((True, None), (False, "none")):
                mm = m & (TRUE_TYPE == typ) & (INFORM == inf)
                if mm.any():
                    ax.plot(X[mm], Y[mm], mk, ms=5, color=col_, mfc=fill, alpha=0.8,
                            label=f"true {typ}" + ("" if inf else " (not informative)"))
    else:
        ax.plot(X[m & ~INFORM], Y[m & ~INFORM], "o", ms=3, color=CLR[reg], mfc="none", alpha=0.4,
                label="not informative")
        ax.plot(X[m & INFORM], Y[m & INFORM], "o", ms=3.5, color=CLR[reg], alpha=0.6,
                label="informative")
    ax.spines[["top", "right"]].set_visible(False)


fig, axes = plt.subplots(2, 2, figsize=(14, 11), squeeze=False)
for c, reg in enumerate(("CP", "CB")):
    m = (regions_run == reg) & ~SPARSE
    # top: orientation vs extent
    ax = axes[0, c]
    _dots(ax, O, E, O_N, E_N, m)
    ax.axhline(E_THR, color="k", lw=0.8, ls="--")
    ax.plot([O_THR, O_THR], [E_THR, 1.3], color="k", lw=0.8, ls="--")
    ax.plot([-O_THR, -O_THR], [E_THR, 1.3], color="k", lw=0.8, ls="--")
    ax.axhline(E_MIN_R[reg], color="tab:olive", lw=1, ls="--")
    for x_, y_, t_ in ((0.75, 1.15, "column"), (-0.75, 1.15, "diagonal"), (0, 1.15, "broad"),
                       (0, (E_MIN_R[reg] + E_THR) / 2, "tiling / scattered (see below)")):
        ax.text(x_, y_, t_, ha="center", fontsize=9, color="0.3")
    ax.text(-0.98, E_MIN_R[reg] / 2, "noise-sized", fontsize=8, color="tab:olive", va="center")
    ax.set(xlim=(-1, 1), ylim=(0, 1.3), xlabel="orientation O  (-1 diagonal ... +1 vertical)",
           ylabel="extent E  (SD along kinematic time / that of a full stripe)",
           title=f"{reg}: {int((m & HAS_CL).sum())} maps with a cluster")
    ax.legend(fontsize=7, loc="center right")
    # bottom: dominance vs extent (only the small clusters are split by D)
    ax = axes[1, c]
    _dots(ax, D, E, D_N, E_N, m)
    ax.axhline(E_THR, color="k", lw=0.8, ls="--")
    ax.plot([D_MIN, D_MIN], [E_MIN_R[reg], E_THR], color="k", lw=0.8, ls="--")
    ax.axhline(E_MIN_R[reg], color="tab:olive", lw=1, ls="--")
    ax.text(0.25, (E_MIN_R[reg] + E_THR) / 2, "scattered", ha="center", fontsize=9, color="0.3")
    ax.text(0.8, (E_MIN_R[reg] + E_THR) / 2, "tiling", ha="center", fontsize=9, color="0.3")
    ax.text(0.02, E_MIN_R[reg] / 2, "noise-sized", fontsize=8, color="tab:olive", va="center")
    ax.set(xlim=(0, 1.02), ylim=(0, 1.3),
           xlabel="dominance D  (share of the map's hot mass in the main cluster)",
           ylabel="extent E")
fig.suptitle("STEP 3 / shape of every map's hot cluster -- every neuron with a cluster is a dot "
             "(dashed: the class boundaries)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.94)); _save(fig, "step3_OE_plane")

# the percentages
_subsets = [("all", np.ones(len(uids), bool))] + ([("reliable (22)", REL22)] if REL22.any() else [])
PCT = {}
print("\nclass percentages (of the neurons that are not too sparse):")
for sname, sset in _subsets:
    for reg in ("CP", "CB"):
        m = (regions_run == reg) & ~SPARSE & sset
        if not m.any(): continue
        PCT[(sname, reg)] = {c: 100 * float(np.mean(CLASS[m] == c)) for c in CLASSES}
        n_inf = int((m & INFORM).sum())
        print(f"  {sname:14s} {reg} n={int(m.sum()):4d} | " + "  ".join(
            f"{c} {PCT[(sname, reg)][c]:4.1f}%" for c in CLASSES))
        print(f"  {'':14s}       of the {n_inf} informative: " + "  ".join(
            f"{c} {100*np.mean(CLASS[m & INFORM] == c):4.1f}%" for c in CLASSES[:6]))

fig, axes = plt.subplots(1, len(_subsets), figsize=(6 * len(_subsets), 4.5), squeeze=False)
for ax, (sname, sset) in zip(axes[0], _subsets):
    regs = [r for r in ("CP", "CB") if (sname, r) in PCT]
    bottom = np.zeros(len(regs))
    for c in CLASSES:
        v = np.array([PCT[(sname, r)][c] for r in regs])
        ax.bar(range(len(regs)), v, bottom=bottom, color=CLASS_CLR[c], label=c, width=0.6)
        for x_, (b_, v_) in enumerate(zip(bottom, v)):
            if v_ >= 4: ax.text(x_, b_ + v_ / 2, f"{v_:.0f}%", ha="center", va="center", fontsize=8)
        bottom += v
    ax.set_xticks(range(len(regs))); ax.set_xticklabels(regs)
    ax.set(ylabel="% of neurons (not too sparse)", title=f"{sname} neurons")
    ax.spines[["top", "right"]].set_visible(False)
axes[0, -1].legend(fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
fig.suptitle("STEP 4 / map class per neuron", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.93)); _save(fig, "step4_class_pct")

if TRUE_TYPE is not None:
    print("\nSYNTHETIC check -- % of each TRUE type called as each class:")
    print(f"  {'true type':12s} {'expected':16s} | " + " ".join(f"{c[:10]:>10s}" for c in CLASSES))
    for typ in TYPE_STYLE:
        mm = (TRUE_TYPE == typ) & ~SPARSE
        if not mm.any(): continue
        print(f"  {typ:12s} {EXPECTED[typ]:16s} | " + " ".join(
            f"{100*np.mean(CLASS[mm] == c):9.0f}%" for c in CLASSES)
              + f"   | median D {np.nanmedian(D[mm]):.2f}")


# %% ── 6. does the class mix explain the mean map? ───────────────────────────
# 6. cell --> per region: the mean excess-MI map of each class, and each class's
# SHARE of the region's mean map (sum over the class's neurons / all neurons),
# inside the band where the mean map is strongest (its top BAND_FRAC cells).
BAND_FRAC = 0.25
SHARE = {}
fig, axes = plt.subplots(2, 8, figsize=(28, 6.5), squeeze=False)
for r, reg in enumerate(("CP", "CB")):
    m = (regions_run == reg) & ~SPARSE
    if not m.any(): continue
    MEAN = np.nanmean(EXC[m], axis=0)
    band = MEAN >= np.nanpercentile(MEAN, 100 * (1 - BAND_FRAC))
    _draw(axes[r, 0], MEAN, f"{reg} mean map, all {int(m.sum())} neurons\n(white: its top "
                            f"{BAND_FRAC:.0%} = the band)", band)
    tot = np.nansum(np.where(band, np.nansum(EXC[m], axis=0), 0))
    SHARE[reg] = {}
    for c, cls in enumerate(("column", "diagonal", "tiling", "scattered", "broad", "noise-sized",
                             "not informative"), start=1):
        mc = m & (CLASS == cls)
        SHARE[reg][cls] = (float(np.nansum(np.where(band, np.nansum(EXC[mc], axis=0), 0)) / tot)
                           if mc.any() and tot != 0 else 0.0)
        if mc.any():
            _draw(axes[r, c], np.nanmean(EXC[mc], axis=0),
                  f"{cls}: mean of {int(mc.sum())} neurons\nshare of the band: {100*SHARE[reg][cls]:.0f}%")
        else:
            axes[r, c].axis("off")
    print(f"{reg}: share of the mean map's band -- " + "  ".join(
        f"{k} {100*v:.0f}%" for k, v in SHARE[reg].items()))
fig.suptitle("STEP 5 / the mean map split by class: who carries the band? (each panel its own colour "
             "scale; x: neural window, y: kinematic window, s)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.94)); _save(fig, "step5_class_means")


# %% ── 7. where the blobs and stripes sit ───────────────────────────────────
# 7. cell --> tilers: centres of mass on the map, and their distributions in
# neural time, kinematic time and lag. Column-locked: their neural time;
# diagonal-locked: their lag.
com_n, com_k, com_l = _a(SHAPE, "com_n"), _a(SHAPE, "com_k"), _a(SHAPE, "com_lag")
fig = plt.figure(figsize=(18, 8))
gs_ = fig.add_gridspec(2, 4)
for r, reg in enumerate(("CP", "CB")):
    m = (regions_run == reg) & (CLASS == "tiling")
    ax = fig.add_subplot(gs_[r, 0])
    ax.plot(com_n[m], com_k[m], "o", ms=3, color=CLR[reg], alpha=0.6)
    dd = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(dd, dd, color="k", lw=0.8, ls="--")
    ax.set(xlim=extent[:2], ylim=extent[2:], xlabel="neural time (s)", ylabel="kinematic time (s)",
           title=f"{reg} tiling: blob centres (n={int(m.sum())})")
    for c, (vals, lab, cls) in enumerate(((com_n, "neural time (s)", "tiling"),
                                          (com_k, "kinematic time (s)", "tiling"),
                                          (com_l, "lag = neural - kinematic (s)", "tiling")), start=1):
        ax = fig.add_subplot(gs_[r, c])
        for cls_, col_, ls in (("tiling", CLR[reg], "-"), ("column", "tab:orange", "-"),
                               ("diagonal", "tab:red", "-")):
            if c == 1 and cls_ == "diagonal": continue      # neural time: tilers and columns
            if c == 2 and cls_ != "tiling": continue        # kinematic time: tilers only
            if c == 3 and cls_ == "column": continue        # lag: tilers and diagonals
            mm = (regions_run == reg) & (CLASS == cls_) & np.isfinite(vals)
            if mm.sum() > 1:
                ax.hist(vals[mm], 25, histtype="step", lw=2, color=col_, density=True,
                        label=f"{cls_} (n={int(mm.sum())})")
        ax.set(xlabel=lab, ylabel="density"); ax.legend(fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("STEP 6 / where the tilers' blobs sit (centre of mass of the hot cluster); "
             "column-locked: their neural time; diagonal-locked: their lag", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.95)); _save(fig, "step6_centres")


# %% ── 8. save ──────────────────────────────────────────────────────────────
if SAVE_VARS and not SYNTHETIC:
    p = config.var_path(f"{RUN_NAME}_{FEAT_TAG}_{NW_TAG}_ks{KIN_STEP}_cat{CAT_TAG}"
                        f"_n{'all' if N_TEST is None else N_TEST}_results.pkl")
    with open(p, "wb") as fh:
        pickle.dump(dict(uids=uids, regions=regions_run, kin_centers=kin_centers,
                         win_centers=win_centers, sparse=SPARSE, finite=FINITE,
                         top=TOP_R, top_null=TOP_N, inf_thr=INF_THR, informative=INFORM,
                         shape=SHAPE, shape_null=SHAPEN, cls=CLASS, pct=PCT, share=SHARE,
                         reliable22=REL22,
                         params=dict(SMOOTH_SIG=SMOOTH_SIG, TOP_FRAC=TOP_FRAC, INF_PCT=INF_PCT,
                                     HOT_REF=HOT_REF, Z_HOT=Z_HOT, MIN_CLUSTER=MIN_CLUSTER,
                                     E_THR=E_THR, O_THR=O_THR, D_MIN=D_MIN, E_MIN=E_MIN_R, MIN_FINITE=MIN_FINITE,
                                     SHUF_SEED=SHUF_SEED, SRC_12=SRC_12, SRC_22=SRC_22)),
                    fh, protocol=4)
    print(f"results -> {p}")

# %%
