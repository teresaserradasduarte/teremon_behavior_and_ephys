# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# WHAT MAKES UP THE STRIPES OF THE MEAN 2D MAP?
#
# The mean CP map (12) is a vertical stripe: information at one NEURAL time,
# whatever kinematic window is read. That can come from three kinds of neuron:
#   locked to neural time   each neuron is itself a vertical stripe: informative
#                           about EVERY kinematic window at its neural time
#   locked to kinematics    each neuron is a diagonal ridge (as the CB mean map)
#   tiling                  each neuron is a single blob -- one kinematic window,
#                           one neural time -- and only the POPULATION covers the
#                           whole kinematic axis
# This script tells them apart without fitting any model, from two profiles per
# neuron, read on HELD-OUT trials:
#   column profile     MI down the neuron's best COLUMN (fixed neural time), as a
#                      function of kinematic window
#   diagonal profile   MI along the neuron's best DIAGONAL (fixed lag), as a
#                      function of kinematic window
#                           column profile   diagonal profile
#   locked to neural time   flat (full)      peaked
#   locked to kinematics    peaked           flat (full)
#   tiling                  peaked           peaked
#
# THE SPLIT. Each session's trials are split at random into two halves, A and B,
# and every neuron's map is computed on each half separately (excess MI = MI -
# the mean of a few trial shuffles, in each half). Everything that involves a
# CHOICE -- which column, which diagonal, the preferred kinematic window used to
# sort the neurons -- is taken from half A; what is PLOTTED comes from half B.
# Sorting noisy profiles by their own maximum makes a diagonal band even out of
# pure noise; sorted by the other half, noise makes none.
#
# RELIABILITY. Some maps are noise. A neuron's map is RELIABLE when its half-A
# and half-B maps are correlated beyond what noise gives. The noise reference
# comes for free from the shuffles each half needs anyway: one shuffle map minus
# the mean of the others is a map with no information in it, processed exactly
# like the real one. The threshold is the 95th percentile of those null
# correlations, per region. No peak is involved.
#
# The grid is coarser than 12's (KIN_STEP 20 ms) and there are few shuffles per
# half: this is about the SHAPE of single-neuron maps, and adjacent 52.5 ms
# kinematic windows 10 ms apart overlap by 80% anyway.
import os, sys, re, time, pickle
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

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

# ── which reaches ────────────────────────────────────────────────────────────
REACH_CAT    = (1, )   # reach categories (idx_reach_cat)
TRIAL_SELECT = None    # None = every reach of the category, or e.g. {"succ": 1},
                       # {"hit": 0} -- the same convention as 12

# ── kinematic features (joint KSG over all of them) ──────────────────────────
FEATURES_2D = ("x", "y", "z", "speed", "dist_water")

# ── the grid (coarser than 12) ───────────────────────────────────────────────
KIN_WIN       = 21     # frames per kinematic window (52.5 ms), as in 12
KIN_STEP      = 8      # frames between kinematic windows (20 ms; 12 uses 10 ms)
KIN_RANGE     = (-0.34, 0.10)   # (s) re the reach endpoint, window centres
NEU_WIN_WIDTH = 0.07   # (s) spike-count window
NEU_WIN_STEP  = 0.025  # (s)
NEU_WIN_RANGE = (-0.5, 0.5)     # (s) re the reach endpoint, window centres

# ── the split, the shuffles, the neurons ─────────────────────────────────────
N_SHUF        = 5      # shuffles per half: one is kept as the reliability null,
                       # the mean of the other N_SHUF-1 is subtracted (from the
                       # real map AND from that null map, so the two are treated
                       # identically)
MIN_TRIALS_HALF = 12   # a neuron needs at least this many usable trials in EACH half
SPLIT_SEED    = 0      # the A/B split is drawn once per SESSION (all neurons of a
                       # session share it) from this seed
N_TEST        = 30     # neurons PER REGION for a test run, drawn at random from
                       # all mice/sessions (fixed seed). None = every neuron.
N_JOBS        = -1

# ── reliability and profiles ─────────────────────────────────────────────────
REL_STAT      = "hotspot"   # how reliability is measured (both are always computed):
                       # "hotspot"  the cells where half A is strongest (the top
                       #            REL_TOP_FRAC of the map) -- is half B high there
                       #            too? Sensitive to ANY structure, however small:
                       #            a stripe, a ridge or a single blob all put their
                       #            top cells on themselves.
                       # "corr"     correlation of the whole half-A and half-B maps.
                       #            Blind to small structures: a blob covering 2% of
                       #            the map is swamped by the noise everywhere else
                       #            (in SYNTHETIC mode only 20% of tiling neurons pass).
REL_TOP_FRAC  = 0.025  # "hotspot": fraction of the map's cells taken from half A
                       # (2.5% of 22 x 38 cells = 21 cells: about one stripe, one
                       # ridge, or a blob of ~2 x 2 windows' worth)
REL_PCT       = 95     # reliable = the statistic above this percentile of the same
                       # statistic on the null maps of the same region
PREF_Z        = 1.0    # the preferred kinematic window of a profile (used only to
                       # SORT the neurons) = centre of mass of the profile above
                       # z = PREF_Z, the same rule as 12's robust peak

# ── sparse maps ──────────────────────────────────────────────────────────────
MIN_FINITE  = 0.6      # a neuron's map must be defined in at least this fraction of
                       # its cells in BOTH halves. A cell is undefined when the
                       # spike count does not vary across the half's trials (a very
                       # sparse neuron); such maps are mostly empty and are set aside
                       # ("too sparse") and counted, not classified.

# ── synthetic data: see what each kind of neuron gives here ─────────────────
SYNTHETIC   = False    # True = no data and no MI at all: half-A/B maps of made-up
                       # neurons of KNOWN type are built directly, and steps 3-6 run
                       # on them unchanged (step 1, the split, is skipped)
SYN_POP     = {"CP": {"time": 25, "tiling": 25, "noise": 25},
               "CB": {"kinematics": 35, "tiling_diag": 25, "noise": 15}}
                       # neurons per type in each "region". Types:
                       #   time        vertical stripe (locked to neural time)
                       #   kinematics  diagonal ridge  (locked to kinematics)
                       #   tiling      one blob; the blobs' kinematic windows are
                       #               spread over the whole kinematic range, their
                       #               neural times independent of it -> the MEAN
                       #               map is a vertical band (like "time")
                       #   tiling_diag one blob, placed along the diagonal (its
                       #               neural time follows its kinematic window)
                       #               -> the MEAN map is a diagonal band (like
                       #               "kinematics")
                       #   noise       no information at all
SYN_AMP     = 0.06     # (bits) peak excess MI of a synthetic neuron
SYN_WIDTH   = 0.04     # (s) width (SD) of the stripe, ridge and blob
SYN_NOISE   = 0.025    # (bits) SD of the noise of each half-map
SYN_SMOOTH  = 1.0      # (cells) the noise is smoothed over neighbouring cells:
                       # real maps are smooth because neighbouring windows overlap
SYN_SEED    = 0

# ── outputs ──────────────────────────────────────────────────────────────────
CACHE_MODE  = "load"   # "load" = reuse the split maps if every parameter matches
SAVE_CACHE  = True
SAVE_FIGS   = True
SAVE_VARS   = True
SHOW_FIGS   = True
N_EXAMPLES  = 4        # example neurons per region in the step figures
OUT_SUBDIR  = "MI_2D_split_profiles"
RUN_NAME    = "split2D"
CMAP        = "viridis"
CLR = {"CP": np.array([106, 176, 134]) / 255, "CB": np.array([154, 123, 193]) / 255}

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


# %% ── 1. data, features, estimator (from 09) ───────────────────────────────
# 1. cell --> 09's cells 0-2 give the reaches in mm, the water features and the
# KSG estimator, with this script's FEATURES / categories (as 12 does). Skipped
# in SYNTHETIC mode, which needs only the grid.
if not SYNTHETIC:
    _src9 = open(os.path.join(_HERE, "09_MI_multiD_kin.py"), encoding="utf-8").read()
    _head = _src9[:re.search(r"(?m)^# %% ── 3\.", _src9).start()]
    for _nm, _val in (("FEATURES", tuple(FEATURES_2D)), ("N_TEST", None),
                      ("CAT_KEEP", set(REACH_CAT))):
        _head, _n = re.subn(rf"(?m)^{_nm}\s*=.*$", f"{_nm} = {_val!r}", _head, count=1)
        if not _n: raise RuntimeError(f"could not find {_nm} in 09 -- renamed?")
    _MINE = {k: globals()[k] for k in ("OUT_SUBDIR", "RUN_NAME", "SAVE_FIGS", "SAVE_VARS",
                                       "SHOW_FIGS", "N_JOBS", "N_TEST", "CLR", "CMAP")}
    # 09 prints a description of its OWN run; none of it describes this script
    import io as _io, contextlib as _ctx
    _buf = _io.StringIO()
    with _step("load data + 09 setup"):
        try:
            with _ctx.redirect_stdout(_buf):
                exec(compile(_head, "09_head", "exec"))
        except BaseException:
            print(_buf.getvalue()); raise
    globals().update(_MINE)
    if any(f in ("PC1", "PC2") for f in FEATURES_2D):
        raise ValueError("PC features are not supported here; use 12 for those")
else:
    FPS = 400.0

_LABEL_KEY = {"hit": "idx_reach_hit", "succ": "idx_reach_succ",
              "LCR": "idx_reach_LCR", "PP": "idx_reach_PP", "cat": "idx_reach_cat"}


def _select_mask(labels):
    n = len(np.asarray(labels["idx_reach_cat"]).ravel())
    m = np.ones(n, bool)
    for k, v in (TRIAL_SELECT or {}).items():
        a = np.asarray(labels[_LABEL_KEY[k]], float).ravel()
        m &= np.isin(a, np.asarray(np.atleast_1d(v), float))
    return m


SEL_TAG  = "" if not TRIAL_SELECT else "-".join(
    f"{k}{'+'.join(str(x) for x in np.atleast_1d(v))}" for k, v in sorted(TRIAL_SELECT.items()))
CAT_TAG  = "".join(str(c) for c in sorted(REACH_CAT))
NW_TAG   = f"nw{int(round(NEU_WIN_WIDTH * 1000))}"
FEAT_TAG = "-".join(FEATURES_2D)
OUT_DIR_C = (f"{OUT_SUBDIR}_SYNTHETIC" if SYNTHETIC else
             f"{OUT_SUBDIR}_{FEAT_TAG}_{NW_TAG}_cat{CAT_TAG}"
             + (f"_{SEL_TAG}" if SEL_TAG else "")
             + ("" if N_TEST is None else f"_test{N_TEST}"))
fdir = config.fig_dir(OUT_DIR_C)

# ── the grids (the reaches are cut to config.REACH_N_FRAMES frames) ──────────
win_starts  = np.arange(NEU_WIN_RANGE[0], NEU_WIN_RANGE[1] - NEU_WIN_WIDTH + NEU_WIN_STEP / 2,
                        NEU_WIN_STEP)
win_centers = win_starts + NEU_WIN_WIDTH / 2
kin_starts  = np.arange(0, config.REACH_N_FRAMES - KIN_WIN + 1, KIN_STEP)
kin_centers = (kin_starts + KIN_WIN / 2 - config.REACH_ENDPOINT_FRAME) / FPS
_in = (kin_centers >= KIN_RANGE[0]) & (kin_centers <= KIN_RANGE[1])
kin_starts, kin_centers = kin_starts[_in], kin_centers[_in]
n_kin, n_wins = len(kin_starts), len(win_starts)

if not SYNTHETIC:
    # the kinematic feature of every window, per session (as 12's cell 2)
    NEEDS_WATER_2D = any(f in ("dist_water", "speed_water") for f in FEATURES_2D)
    KIN = {}
    for key in REACH_MM:
        rp  = REACH_MM[key]
        vel = savgol_filter(rp, SAVGOL_WIN, SAVGOL_ORD, deriv=1, delta=1/FPS, axis=2)
        spd = np.linalg.norm(vel, axis=1)
        acc = savgol_filter(spd, SAVGOL_WIN, SAVGOL_ORD, deriv=1, delta=1/FPS, axis=1)
        d = {"x": rp[:, 0], "y": rp[:, 1], "z": rp[:, 2], "speed": spd, "accel": acc}
        if NEEDS_WATER_2D:
            d["dist_water"], d["speed_water"] = WATER[key]["dist_water"], WATER[key]["speed_water"]
        KIN[key] = np.stack([np.stack([d[f][:, s:s + KIN_WIN].mean(axis=1) for s in kin_starts],
                                      axis=1) for f in FEATURES_2D], axis=2)   # (n, n_kin, d)
    print(f"features {list(FEATURES_2D)} (joint KSG, k={K_KSG})  |  reach cat {sorted(REACH_CAT)}"
          + (f"  |  trials: {SEL_TAG}" if SEL_TAG else ""))
else:
    print("SYNTHETIC data: made-up neurons of known type, no MI computed  |  "
          + "  ".join(f"{r}: " + ", ".join(f"{t} {n}" for t, n in pop.items())
                      for r, pop in SYN_POP.items()))
print(f"grid: {n_kin} kinematic windows ({kin_centers[0]*1000:.0f}..{kin_centers[-1]*1000:.0f} ms, "
      f"step {KIN_STEP/FPS*1000:.0f} ms) x {n_wins} neural windows "
      f"({win_centers[0]*1000:.0f}..{win_centers[-1]*1000:.0f} ms, {NEU_WIN_WIDTH*1000:.0f} ms wide)")
print(f"outputs -> {OUT_DIR_C}/")


# %% ── 2. the A/B split, and which neurons ──────────────────────────────────
# 2. cell --> one random A/B split per SESSION, over the usable trials (finite
# kinematics in every window, the right category and trial selection). Then the
# neurons: those with at least MIN_TRIALS_HALF usable trials in each half; in a
# test run, N_TEST per region drawn at random across all sessions.
if not SYNTHETIC:
    _rs = np.random.default_rng(SPLIT_SEED)
    USABLE, HALF = {}, {}
    for key in sorted(KIN):
        n0 = next(n for n in NEURONS if (n["mouse"], n["sess"]) == key)
        ok = (np.isfinite(KIN[key]).all(axis=(1, 2))
              & np.isin(n0["labels"]["idx_reach_cat"], list(REACH_CAT))
              & _select_mask(n0["labels"]))
        idx = np.flatnonzero(ok)
        perm = _rs.permutation(idx)
        h = np.full(ok.size, -1)                     # -1 not used, 0 = A, 1 = B
        h[perm[: idx.size // 2]] = 0
        h[perm[idx.size // 2:]] = 1
        USABLE[key], HALF[key] = ok, h

    cand = {"CP": [], "CB": []}
    for i, n in enumerate(NEURONS):
        if n["region"] not in cand: continue
        h = HALF[(n["mouse"], n["sess"])]
        if min((h == 0).sum(), (h == 1).sum()) >= MIN_TRIALS_HALF:
            cand[n["region"]].append(i)
    _rn = np.random.default_rng(1)
    run_idx = []
    for reg in ("CP", "CB"):
        c = cand[reg]
        run_idx += sorted(c if N_TEST is None or N_TEST >= len(c)
                          else _rn.choice(c, N_TEST, replace=False).tolist())
    uids        = [NEURONS[i]["uid"] for i in run_idx]
    regions_run = np.array([NEURONS[i]["region"] for i in run_idx])
    TRUE_TYPE   = None
    for reg in ("CP", "CB"):
        print(f"{reg}: {int((regions_run == reg).sum())} neurons "
              f"(of {len(cand[reg])} with >= {MIN_TRIALS_HALF} usable trials in each half)")

    # FIGURE: how many trials each half has, per session
    fig, ax = plt.subplots(figsize=(max(7, 0.32 * len(HALF)), 3.6))
    _k = sorted(HALF)
    _a = [int((HALF[k] == 0).sum()) for k in _k]; _b = [int((HALF[k] == 1).sum()) for k in _k]
    ax.bar(np.arange(len(_k)) - 0.2, _a, width=0.4, color="0.35", label="half A")
    ax.bar(np.arange(len(_k)) + 0.2, _b, width=0.4, color="0.7", label="half B")
    ax.axhline(MIN_TRIALS_HALF, color="tab:red", lw=1, ls="--",
               label=f"MIN_TRIALS_HALF = {MIN_TRIALS_HALF}")
    ax.set_xticks(np.arange(len(_k)))
    ax.set_xticklabels([f"{k[0].split('_')[-1]}|{k[1]}" for k in _k], rotation=90, fontsize=6)
    ax.set(ylabel="usable trials", title="STEP 1 / the A/B split: usable trials in each half, per session")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)
    plt.tight_layout(); _save(fig, "step1_split")


# %% ── 3. the maps of each half ─────────────────────────────────────────────
# 3. cell --> per neuron and half: the MI map and N_SHUF trial-shuffled maps.
#   excess = MI       - mean(shuffles 2..N)
#   null   = shuffle1 - mean(shuffles 2..N)    (a map with no information,
#                                               processed exactly like excess)
# The shuffles permute whole kinematic rows WITHIN the half.
# In SYNTHETIC mode the four maps are built directly: the neuron's true map plus
# independent smoothed noise in each half, and pure noise for the null maps.
def _half_maps(K, C, k, seed, n_shuf):
    """K (n, n_kin, d), C (n, n_wins) -> excess, null, each (n_kin, n_wins)"""
    rs = np.random.default_rng(seed)
    perms = [rs.permutation(K.shape[0]) for _ in range(n_shuf)]
    mi = np.full((K.shape[1], C.shape[1]), np.nan)
    sh = np.full((n_shuf,) + mi.shape, np.nan)
    ok_w = [w for w in range(C.shape[1]) if np.unique(C[:, w]).size >= 2]
    for ki in range(K.shape[1]):
        X = K[:, ki, :]
        for w in ok_w:
            c = C[:, w]
            mi[ki, w] = _ksg_md(X, c, k, seed + ki * 97 + w)
            for s, p in enumerate(perms):
                sh[s, ki, w] = _ksg_md(X[p], c, k, seed + ki * 97 + w)
    base = np.nanmean(sh[1:], axis=0)
    return (mi - base).astype(np.float32), (sh[0] - base).astype(np.float32)


def _neuron_maps(KA, CA, KB, CB, k, seed, n_shuf):
    eA, nA = _half_maps(KA, CA, k, seed, n_shuf)
    eB, nB = _half_maps(KB, CB, k, seed + 50_000, n_shuf)
    return eA, eB, nA, nB


def _synthetic():
    """-> EA, EB, NA, NB, uids, regions, true types, for the SYN_POP population"""
    from scipy.ndimage import gaussian_filter
    rs = np.random.default_rng(SYN_SEED)
    Kc, Nc = np.meshgrid(kin_centers, win_centers, indexing="ij")
    g = lambda x: np.exp(-0.5 * (x / SYN_WIDTH) ** 2)
    noise = lambda: SYN_NOISE * (gaussian_filter(rs.normal(size=Kc.shape), SYN_SMOOTH)
                                 / max(gaussian_filter(rs.normal(size=(200, 200)),
                                                       SYN_SMOOTH).std(), 1e-9))
    maps, U, R, T = [], [], [], []
    for reg, pop in SYN_POP.items():
        for typ, n in pop.items():
            for j in range(n):
                t0 = rs.uniform(-0.15, 0.15)                    # where in neural time
                # tiling: kinematic windows spread evenly, one per n-th of the range
                kw = kin_centers[0] + (j + rs.uniform()) / n * (kin_centers[-1] - kin_centers[0])
                if typ == "time":         M = g(Nc - t0)
                elif typ == "kinematics": M = g((Nc - Kc) - rs.uniform(-0.1, 0.1))
                elif typ == "tiling":     M = g(Nc - t0) * g(Kc - kw)
                elif typ == "tiling_diag":
                    M = g(Nc - (kw + rs.uniform(-0.1, 0.1))) * g(Kc - kw)
                elif typ == "noise":      M = np.zeros_like(Kc)
                else: raise ValueError(f"unknown synthetic type {typ!r}")
                M = SYN_AMP * M
                maps.append((M + noise(), M + noise(), noise(), noise()))
                U.append(f"{reg}|SYN_{typ}|#{j}"); R.append(reg); T.append(typ)
    EA, EB, NA, NB = (np.stack([m[i] for m in maps]).astype(np.float32) for i in range(4))
    return EA, EB, NA, NB, U, np.array(R), np.array(T)


if SYNTHETIC:
    EA, EB, NA, NB, uids, regions_run, TRUE_TYPE = _synthetic()
    run_idx = list(range(len(uids)))
    print(f"{len(uids)} synthetic neurons built")
else:
    KEY = dict(features=list(FEATURES_2D), REACH_CAT=sorted(REACH_CAT), select=SEL_TAG,
               KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP, KIN_RANGE=list(KIN_RANGE),
               NEU_WIN_WIDTH=NEU_WIN_WIDTH, NEU_WIN_STEP=NEU_WIN_STEP,
               NEU_WIN_RANGE=list(NEU_WIN_RANGE), N_SHUF=N_SHUF,
               MIN_TRIALS_HALF=MIN_TRIALS_HALF, SPLIT_SEED=SPLIT_SEED, K=K_KSG,
               NEURAL_VAR=NEURAL_VAR, uids=uids, v="v1")
    CACHE_PATH = config.var_path(
        f"{RUN_NAME}_{FEAT_TAG}_{NW_TAG}_ks{KIN_STEP}_cat{CAT_TAG}"
        + (f"_{SEL_TAG}" if SEL_TAG else "")
        + f"_s{N_SHUF}_n{'all' if N_TEST is None else N_TEST}_cache.pkl")
    _c = None
    if CACHE_MODE == "load" and CACHE_PATH.is_file():
        with open(CACHE_PATH, "rb") as fh: _c = pickle.load(fh)
        if _c.get("key") != KEY:
            print("cache exists but parameters differ -> recomputing"); _c = None
    if _c is not None:
        EA, EB, NA, NB = _c["EA"], _c["EB"], _c["NA"], _c["NB"]
        print(f"LOADED {CACHE_PATH.name} -- no recompute")
    else:
        _calls = len(run_idx) * 2 * n_kin * n_wins * (1 + N_SHUF)
        print(f"ETA: ~{_calls:,} KSG calls; at ~2500 calls/s that is ~{_calls/2500/60:.0f} min")
        jobs = []
        for pi, i in enumerate(run_idx):
            n = NEURONS[i]; key = (n["mouse"], n["sess"]); h = HALF[key]
            C = np.asarray(mi_rates.windowed(n["spikes"], win_starts, NEU_WIN_WIDTH,
                                             NEURAL_VAR), float)
            jobs.append((KIN[key][h == 0], C[h == 0], KIN[key][h == 1], C[h == 1], 100 + pi))
        with _step("split maps (A, B, nulls)"):
            res = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
                delayed(_neuron_maps)(KA, CA, KB, CB, K_KSG, sd, N_SHUF)
                for KA, CA, KB, CB, sd in jobs)
        EA, EB, NA, NB = (np.stack([r[j] for r in res]) for j in range(4))
        if SAVE_CACHE:
            with open(CACHE_PATH, "wb") as fh:
                pickle.dump(dict(key=KEY, EA=EA, EB=EB, NA=NA, NB=NB), fh, protocol=4)
            print(f"cache -> {CACHE_PATH}")

# a map undefined in too many cells (a very sparse neuron) is set aside
FINITE = np.array([min(np.isfinite(EA[pi]).mean(), np.isfinite(EB[pi]).mean())
                   for pi in range(len(run_idx))])
SPARSE = FINITE < MIN_FINITE
for reg in ("CP", "CB"):
    m = regions_run == reg
    if m.any():
        print(f"{reg}: {int(SPARSE[m].sum())}/{int(m.sum())} maps too sparse "
              f"(defined in < {MIN_FINITE:.0%} of cells) -> set aside")

extent = [win_centers[0], win_centers[-1], kin_centers[0], kin_centers[-1]]


def _draw_map(ax, M, title, lim=None):
    lim = lim if lim is not None else (np.nanpercentile(np.abs(M), 99) if np.isfinite(M).any() else 1)
    im = ax.imshow(M, origin="lower", aspect="auto", extent=extent, cmap=CMAP, vmin=-lim / 3, vmax=lim)
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(d, d, color="w", lw=0.8, ls="--", alpha=0.7)
    ax.set_title(title, fontsize=8); ax.tick_params(labelsize=6)
    return im


# example neurons: spread over each region (fixed seed)
_re = np.random.default_rng(2)
EX = {reg: sorted(_re.choice(np.flatnonzero(regions_run == reg),
                             min(N_EXAMPLES, int((regions_run == reg).sum())), replace=False))
      for reg in ("CP", "CB") if (regions_run == reg).any()}

# FIGURE: what each half gives, for example neurons
_ex = [pi for reg in EX for pi in EX[reg]]
fig, axes = plt.subplots(len(_ex), 4, figsize=(13, 2.4 * len(_ex)), squeeze=False)
for r, pi in enumerate(_ex):
    lim = np.nanpercentile(np.abs(np.r_[EA[pi].ravel(), EB[pi].ravel()]), 99)
    for c, (M, t) in enumerate(((EA[pi], "half A: excess MI"), (EB[pi], "half B: excess MI"),
                                (NA[pi], "half A: null (no information)"),
                                (NB[pi], "half B: null"))):
        _draw_map(axes[r, c], M, t, lim)
    axes[r, 0].set_ylabel(f"{uids[pi]}\nkinematic (s)", fontsize=7)
fig.suptitle("STEP 2 / each neuron's map on the two halves of its trials, and the matching "
             "null maps (one shuffle minus the mean of the others)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.97)); _save(fig, "step2_split_maps")


# %% ── 4. reliability ───────────────────────────────────────────────────────
# 4. cell --> is the map the same on two independent halves of the trials? Two
# statistics, both computed identically on the real maps and on the null maps:
#   hotspot  mean of half B over the top REL_TOP_FRAC cells of half A (in units
#            of half B's own SD, so maps of different overall size compare)
#   corr     correlation of the whole half-A and half-B maps
# Reliable = the REL_STAT statistic above the REL_PCT percentile of the null
# values of the same region. Sparse maps are never reliable.
def _corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5 or np.std(a[m]) == 0 or np.std(b[m]) == 0: return np.nan
    return float(np.corrcoef(a[m], b[m])[0, 1])


def _hotspot(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 10 or np.std(b[m]) == 0: return np.nan
    av, bv = a[m], b[m]
    k = max(1, int(round(REL_TOP_FRAC * av.size)))
    top = np.argsort(av)[-k:]
    return float((bv[top].mean() - bv.mean()) / bv.std())


STATS = {}
for nm, fn in (("hotspot", _hotspot), ("corr", _corr)):
    STATS[nm] = (np.array([fn(EA[pi], EB[pi]) for pi in range(len(run_idx))]),
                 np.array([fn(NA[pi], NB[pi]) for pi in range(len(run_idx))]))
R_REAL, R_NULL = STATS[REL_STAT]
THR = {reg: float(np.nanpercentile(R_NULL[(regions_run == reg) & ~SPARSE], REL_PCT))
       for reg in ("CP", "CB") if ((regions_run == reg) & ~SPARSE).any()}
RELIABLE = np.array([(not SPARSE[pi]) and np.isfinite(R_REAL[pi])
                     and R_REAL[pi] > THR.get(regions_run[pi], np.inf)
                     for pi in range(len(run_idx))])
# PEAK-TO-NOISE of every map: its 99.5th-percentile cell (half A) over the SD of
# its null map. Detectability depends on it AND on the shape: in SYNTHETIC mode a
# stripe or a ridge is found reliably from a ratio of ~3, a single blob (tiling)
# only from ~4-6, because a blob spreads its information over far fewer cells.
# Compare the real ratios here with that, before reading the fractions of types.
SNR = np.array([np.nanpercentile(EA[pi], 99.5) / np.nanstd(NA[pi])
                if np.isfinite(NA[pi]).sum() > 10 and np.nanstd(NA[pi]) > 0 else np.nan
                for pi in range(len(run_idx))])
print("\npeak-to-noise of the maps (half A 99.5th-pct cell / SD of the null map):")
for reg in ("CP", "CB"):
    m = (regions_run == reg) & ~SPARSE
    if m.any():
        q = np.nanpercentile(SNR[m], [25, 50, 75, 90])
        print(f"  {reg}: median {q[1]:.1f}  (IQR {q[0]:.1f}-{q[2]:.1f}, 90th pct {q[3]:.1f})")
print(f"\nreliability ({REL_STAT}): real maps vs the {REL_PCT}th percentile of the null maps")
for reg in THR:
    m = regions_run == reg
    _alt = "corr" if REL_STAT == "hotspot" else "hotspot"
    _rr, _rn = STATS[_alt]
    _thr_alt = np.nanpercentile(_rn[m & ~SPARSE], REL_PCT)
    print(f"  {reg}: threshold {THR[reg]:.3f}  |  reliable {int(RELIABLE[m].sum())}/{int(m.sum())} "
          f"({100*RELIABLE[m].mean():.0f}%; {int(SPARSE[m].sum())} too sparse)  |  for comparison, "
          f"by '{_alt}': {int(((_rr > _thr_alt) & ~SPARSE & m).sum())}")

# FIGURE: the distributions, and maps across the range of r
fig = plt.figure(figsize=(15, 8.5))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.1], hspace=0.45)
for c, reg in enumerate(("CP", "CB")):
    if reg not in THR: continue
    ax = fig.add_subplot(gs[0, 2 * c:2 * c + 2])
    m = regions_run == reg
    _v = np.r_[R_NULL[m], R_REAL[m]]; _v = _v[np.isfinite(_v)]
    b = np.linspace(np.percentile(_v, 0.5), np.percentile(_v, 99.5), 40) if _v.size else np.linspace(-1, 1, 40)
    ax.hist(np.clip(R_NULL[m], b[0], b[-1]), bins=b, color="0.6", alpha=0.7, label="null maps")
    ax.hist(np.clip(R_REAL[m], b[0], b[-1]), bins=b, color=CLR[reg], alpha=0.7, label="real maps")
    ax.axvline(THR[reg], color="k", ls="--", lw=1.2,
               label=f"threshold ({REL_PCT}th pct of null): {THR[reg]:.2f}")
    ax.set(xlabel=("half B at half A's top cells (SD of half B above its mean)"
                   if REL_STAT == "hotspot" else "split-half correlation (half A vs half B map)"),
           ylabel="neurons",
           title=f"{reg}: {int(RELIABLE[m].sum())}/{int(m.sum())} reliable")
    ax.title.set_fontsize(10); ax.legend(fontsize=7); ax.spines[["top", "right"]].set_visible(False)
# examples: lowest, median, highest r (all neurons)
_o = np.argsort(np.nan_to_num(R_REAL, nan=-9))
_pick = [_o[0], _o[len(_o) // 2], _o[-1]] if len(_o) >= 3 else list(_o)
sub = gs[1, :].subgridspec(2, len(_pick), hspace=0.35)
for c, pi in enumerate(_pick):
    lim = np.nanpercentile(np.abs(np.r_[EA[pi].ravel(), EB[pi].ravel()]), 99)
    _draw_map(fig.add_subplot(sub[0, c]), EA[pi],
              f"{uids[pi]}  {REL_STAT} = {R_REAL[pi]:.2f}\nhalf A", lim)
    _draw_map(fig.add_subplot(sub[1, c]), EB[pi], "half B", lim)
fig.suptitle(f"STEP 3 / reliability ({REL_STAT}): is the map the same on two independent halves "
             f"of the trials? (bottom: the least, median and most reliable neuron)", fontsize=10)
_save(fig, "step3_reliability")


# %% ── 5. the two profiles of every neuron ──────────────────────────────────
# 5. cell --> from half A: the best COLUMN (neural window with the highest
# column-mean MI) and the best DIAGONAL (lag with the highest diagonal-mean MI).
# Along each, the profile over kinematic windows -- read from half A to choose
# the preferred window (for sorting), and from half B to be PLOTTED.
LAG     = win_centers[None, :] - kin_centers[:, None]           # (n_kin, n_wins)
_lb     = NEU_WIN_STEP
LAG_IDX = np.round(LAG / _lb).astype(int); _l0 = LAG_IDX.min(); LAG_IDX -= _l0
lag_axis = (np.arange(LAG_IDX.max() + 1) + _l0) * _lb
# only lags that every kinematic window reaches: a diagonal crossing just a
# corner of the map would give a profile with most of it missing
_full_lag = ((lag_axis >= LAG.min(axis=1).max() - _lb / 2) &
             (lag_axis <= LAG.max(axis=1).min() + _lb / 2))


def _diag_mean(M):
    ok = np.isfinite(M)
    s = np.bincount(LAG_IDX[ok], M[ok], minlength=lag_axis.size)
    n = np.bincount(LAG_IDX[ok], minlength=lag_axis.size)
    return np.where(n > 0, s / np.maximum(n, 1), np.nan)


def _diag_profile(M, lag):
    """M along the diagonal of lag `lag`: one value per kinematic window (the
    neural window whose centre is nearest kin_centre + lag)"""
    out = np.full(n_kin, np.nan)
    for k in range(n_kin):
        j = int(np.argmin(np.abs(win_centers - (kin_centers[k] + lag))))
        if abs(win_centers[j] - (kin_centers[k] + lag)) <= NEU_WIN_STEP / 2 + 1e-9:
            out[k] = M[k, j]
    return out


def _pref(prof):
    """preferred kinematic window (s) of a profile: centre of mass above z=PREF_Z"""
    v = np.asarray(prof, float)
    if np.isfinite(v).sum() < 3 or np.nanstd(v) == 0: return np.nan
    z = (v - np.nanmean(v)) / np.nanstd(v)
    w = np.clip(np.nan_to_num(z - PREF_Z), 0, None)
    if w.sum() == 0: w = (np.nan_to_num(z) == np.nanmax(z)).astype(float)
    return float(np.sum(w * kin_centers) / w.sum())


def _profiles(MA, MB):
    """choose column and diagonal on MA, read both profiles on MA and MB"""
    with np.errstate(all="ignore"):
        cm = np.nanmean(MA, axis=0) if np.isfinite(MA).any() else np.full(n_wins, np.nan)
    dm = _diag_mean(MA); dm[~_full_lag] = np.nan
    if not (np.isfinite(cm).any() and np.isfinite(dm).any()):
        # map undefined everywhere (a very sparse neuron): no profiles. It is
        # SPARSE, so it is left out of every later step anyway
        nan = np.full(n_kin, np.nan)
        return dict(col=None, lag=None, colA=nan, colB=nan, diaA=nan, diaB=nan)
    col = int(np.nanargmax(cm))
    lag = float(lag_axis[int(np.nanargmax(dm))])
    return dict(col=col, lag=lag,
                colA=MA[:, col], colB=MB[:, col],
                diaA=_diag_profile(MA, lag), diaB=_diag_profile(MB, lag))


PROF = [_profiles(EA[pi], EB[pi]) for pi in range(len(run_idx))]       # real
PNUL = [_profiles(NA[pi], NB[pi]) for pi in range(len(run_idx))]       # null, same procedure
for P in (PROF, PNUL):
    for d in P:
        d["prefCol"], d["prefDia"] = _pref(d["colA"]), _pref(d["diaA"])

# FIGURE: how the profiles are taken, for example neurons (those with a map)
_ex4 = [pi for pi in _ex if PROF[pi]["col"] is not None]
fig, axes = plt.subplots(max(len(_ex4), 1), 3, figsize=(13, 2.6 * max(len(_ex4), 1)), squeeze=False,
                         gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
for r, pi in enumerate(_ex4):
    d = PROF[pi]; reg = regions_run[pi]
    _draw_map(axes[r, 0], EA[pi], f"{uids[pi]} -- half A, with the chosen column and diagonal")
    axes[r, 0].axvline(win_centers[d["col"]], color="tab:orange", lw=1.5)
    _kk = np.array(kin_centers[[0, -1]])
    axes[r, 0].plot(_kk + d["lag"], _kk, color="tab:red", lw=1.5)
    axes[r, 0].set_xlim(extent[:2])
    for c, (a_, b_, t, pf) in enumerate(((d["colA"], d["colB"], f"column at {win_centers[d['col']]*1000:+.0f} ms",
                                          d["prefCol"]),
                                         (d["diaA"], d["diaB"], f"diagonal at lag {d['lag']*1000:+.0f} ms",
                                          d["prefDia"])), start=1):
        ax = axes[r, c]
        ax.plot(kin_centers * 1000, a_, color="0.6", lw=1.2, label="half A (chooses)")
        ax.plot(kin_centers * 1000, b_, color=CLR[reg], lw=2, label="half B (plotted)")
        if np.isfinite(pf): ax.axvline(pf * 1000, color="k", ls=":", lw=1, label="preferred (A)")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title(t, fontsize=8, color="tab:orange" if c == 1 else "tab:red")
        ax.tick_params(labelsize=6); ax.spines[["top", "right"]].set_visible(False)
        if r == 0: ax.legend(fontsize=6)
        if r == len(_ex4) - 1: ax.set_xlabel("kinematic window centre (ms re reach endpoint)", fontsize=7)
fig.suptitle("STEP 4 / the two profiles: MI down the best COLUMN (orange) and along the best "
             "DIAGONAL (red), chosen on half A, read on half B", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.97)); _save(fig, "step4_profiles")


# %% ── 6. sorted profiles: the figure that tells the options apart ──────────
# 6. cell --> per region, every neuron's half-B profile as a row (z-scored, so
# only its SHAPE shows), sorted by its preferred kinematic window from half A.
#   column profiles:   full rows = locked to neural time; a band = peaked
#   diagonal profiles: full rows = locked to kinematics;  a band = peaked
#   a neuron peaked in BOTH views is a single blob = tiling
# The null row shows what this procedure makes out of maps with no information.
def _z(v):
    v = np.asarray(v, float); s = np.nanstd(v)
    return (v - np.nanmean(v)) / s if s > 0 else v * np.nan


def _sorted_block(P, idx, which):
    pref = np.array([P[pi]["prefCol" if which == "col" else "prefDia"] for pi in idx])
    o = np.argsort(np.nan_to_num(pref, nan=np.inf))
    rows = np.stack([_z(P[pi]["colB" if which == "col" else "diaB"]) for pi in np.asarray(idx)[o]]) \
        if len(idx) else np.zeros((0, n_kin))
    return rows, pref[o]


# rows of each figure: reliable neurons, every defined map ("all"), and the null maps
def _set_rows(reg):
    m = (regions_run == reg) & ~SPARSE
    return [("reliable", PROF, np.flatnonzero(m & RELIABLE)),
            ("all", PROF, np.flatnonzero(m)),
            ("NULL maps", PNUL, np.flatnonzero(m))]


fig, axes = plt.subplots(4, 4, figsize=(16, 14), squeeze=False)
for c0, reg in enumerate(("CP", "CB")):
    for c1, which in enumerate(("col", "dia")):
        c = 2 * c0 + c1
        for r, (sname, P, idx) in enumerate(_set_rows(reg)):
            ax = axes[r, c]
            rows, pref = _sorted_block(P, idx, which)
            if rows.shape[0]:
                ax.imshow(rows, aspect="auto", cmap=CMAP, vmin=-1.5, vmax=2.5, origin="upper",
                          extent=[kin_centers[0] * 1000, kin_centers[-1] * 1000, rows.shape[0], 0],
                          interpolation="nearest")
                ax.plot(pref * 1000, np.arange(rows.shape[0]) + 0.5, "w.", ms=2, alpha=0.6)
            ax.set_title(f"{reg} {sname} (n={len(idx)})\n{'COLUMN' if which == 'col' else 'DIAGONAL'} "
                         f"profiles, sorted by preferred window (half A)", fontsize=8)
            ax.set_ylabel("neurons", fontsize=7); ax.tick_params(labelsize=6)
        # bottom: the average shape of the RELIABLE neurons, against the null maps
        ax = axes[3, c]
        (_, _, idx_rel), _, (_, _, idx_nul) = _set_rows(reg)
        for P, idx, col_, ls, lab in ((PROF, idx_rel, CLR[reg], "-", "reliable"),
                                      (PNUL, idx_nul, "0.5", "--", "null")):
            rows, _ = _sorted_block(P, idx, which)
            if rows.shape[0]:
                ax.plot(kin_centers * 1000, np.nanmean(rows, axis=0), color=col_, ls=ls, lw=2, label=lab)
        ax.axhline(0, color="k", lw=0.5)
        ax.set(xlabel="kinematic window centre (ms re reach endpoint)",
               ylabel="mean z-scored profile", title="average shape, reliable neurons (unsorted)")
        ax.title.set_fontsize(8); ax.legend(fontsize=7); ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("STEP 5 / sorted profiles -- half-B profiles, each row z-scored, sorted by the half-A "
             "preferred window (white dots).  Rows: reliable neurons, all defined maps, null maps\n"
             "COLUMN view: full rows = locked to neural time.  DIAGONAL view: full rows = locked to "
             "kinematics.  A band in BOTH views = tiling.", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.94)); _save(fig, "step5_sorted_profiles")


# %% ── 6b. sorted profiles in NOISE units (height kept) ──────────────────────
# 6b. cell --> the same rows as step 5, but each divided by the SD of that
# neuron's own half-B NULL map instead of z-scored, so a row keeps its HEIGHT:
#   locked to neural time   full bright rows in the COLUMN view
#   locked to kinematics    full bright rows in the DIAGONAL view
#   tiling                  a short band in both views
#   noise                   dark (about 0) everywhere
# (z-scoring a flat-but-high row turns it into noise -- that is why step 5
# shows a band for every type.) SYNTHETIC: one figure with the true types in
# separate blocks, then one with all neurons sorted together.
NOISE_VLIM = (-1, 3)   # colour range, in noise SDs


def _nsd(pi):
    v = NB[pi]; v = v[np.isfinite(v)]
    s = np.std(v) if v.size else np.nan
    return s if s > 0 else np.nan


def _noise_rows(P, idx, which):
    """half-B profiles / noise SD, sorted by the half-A preferred window"""
    idx = np.asarray(idx, int)
    pref = np.array([P[pi]["prefCol" if which == "col" else "prefDia"] for pi in idx])
    o = np.argsort(np.nan_to_num(pref, nan=np.inf))
    rows = (np.stack([np.asarray(P[pi]["colB" if which == "col" else "diaB"], float) / _nsd(pi)
                      for pi in idx[o]]) if idx.size else np.zeros((0, n_kin)))
    return rows, pref[o]


def _noise_panel(ax, blocks, which, title):
    """blocks = [(label, rows, pref)], stacked top to bottom with a line between"""
    y0 = 0
    for lab, rows, pref in blocks:
        if not rows.shape[0]: continue
        ax.imshow(rows, aspect="auto", cmap=CMAP, vmin=NOISE_VLIM[0], vmax=NOISE_VLIM[1],
                  origin="upper", interpolation="nearest",
                  extent=[kin_centers[0] * 1000, kin_centers[-1] * 1000, y0 + rows.shape[0], y0])
        ax.plot(pref * 1000, y0 + np.arange(rows.shape[0]) + 0.5, "w.", ms=2, alpha=0.6)
        if lab:
            ax.axhline(y0, color="w", lw=1.5)
            ax.text(kin_centers[0] * 1000, y0 + 0.5, f" {lab} ({rows.shape[0]})", color="w",
                    va="top", fontsize=7, fontweight="bold")
        y0 += rows.shape[0]
    ax.set_ylim(max(y0, 1), 0)
    ax.set_title(f"{title}\n{'COLUMN' if which == 'col' else 'DIAGONAL'} profiles", fontsize=8)
    ax.set_ylabel("neurons", fontsize=7); ax.tick_params(labelsize=6)


TYPE_ORDER = ("time", "kinematics", "tiling", "tiling_diag", "noise")
_layouts = [("together", None)]
if TRUE_TYPE is not None:
    _layouts = [("by_type", TYPE_ORDER), ("together", None)]
for lname, types in _layouts:
    fig, axes = plt.subplots(3, 4, figsize=(16, 14), squeeze=False)
    for c0, reg in enumerate(("CP", "CB")):
        for c1, which in enumerate(("col", "dia")):
            c = 2 * c0 + c1
            for r, (sname, P, idx) in enumerate(_set_rows(reg)):
                if types is None or sname == "NULL maps":
                    blocks = [("", *_noise_rows(P, idx, which))]
                else:
                    blocks = [(t, *_noise_rows(P, idx[TRUE_TYPE[idx] == t], which)) for t in types]
                _noise_panel(axes[r, c], blocks, which, f"{reg} {sname} (n={len(idx)})")
            axes[2, c].set_xlabel("kinematic window centre (ms re reach endpoint)", fontsize=7)
    if axes[1, -1].images:
        _cb = fig.colorbar(axes[1, -1].images[0], ax=axes, shrink=0.5, pad=0.01)
        _cb.set_label("half-B MI / SD of the neuron's null map", fontsize=8)
    fig.suptitle(f"STEP 5b / sorted profiles in NOISE units{' (true types in blocks)' if types else ''}"
                 f" -- the height is kept.  Rows: reliable neurons, all defined maps, null maps\n"
                 f"COLUMN view: full bright rows = locked to neural time.  DIAGONAL view: full "
                 f"bright rows = locked to kinematics.  A short band in both = tiling.", fontsize=10)
    _save(fig, "step5b_noise" + ("_types" if types else ""))


# %% ── 7. breadth of the two profiles ───────────────────────────────────────
# 7. cell --> one number per profile, peak-free: the PARTICIPATION RATIO of its
# positive part, (sum p)^2 / (N sum p^2) -- 1 when the excess MI is spread
# evenly over all kinematic windows, 1/N when it sits in a single one. From
# half B. Each neuron is a point (column breadth, diagonal breadth):
#   top-left      broad on the column, narrow on the diagonal -> locked to neural time
#   bottom-right  narrow on the column, broad on the diagonal -> locked to kinematics
#   bottom-left   narrow on both                               -> tiling
def _breadth(v):
    p = np.clip(np.nan_to_num(np.asarray(v, float)), 0, None)
    n = np.isfinite(np.asarray(v, float)).sum()
    return float(p.sum() ** 2 / (n * (p ** 2).sum())) if (p ** 2).sum() > 0 and n else np.nan


BR = np.array([[_breadth(d["colB"]), _breadth(d["diaB"])] for d in PROF])
BN = np.array([[_breadth(d["colB"]), _breadth(d["diaB"])] for d in PNUL])
fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), sharex=True, sharey=True)
TYPE_STYLE = {"time": ("tab:orange", "s"), "kinematics": ("tab:red", "^"),
              "tiling": ("tab:blue", "o"), "tiling_diag": ("tab:cyan", "D"),
              "noise": ("0.5", "x")}
for ax, reg in zip(axes, ("CP", "CB")):
    m = (regions_run == reg) & ~SPARSE
    ax.plot(BN[m, 1], BN[m, 0], "o", ms=3, color="0.8", alpha=0.6, label="null maps")
    if TRUE_TYPE is not None:
        # SYNTHETIC: every neuron in the colour of the type it was built as;
        # filled = called reliable, open = not
        for typ, (col_, mk) in TYPE_STYLE.items():
            for rel, fill in ((True, None), (False, "none")):
                mm = m & (TRUE_TYPE == typ) & (RELIABLE == rel)
                if mm.any():
                    ax.plot(BR[mm, 1], BR[mm, 0], mk, ms=6, color=col_, mfc=fill, alpha=0.85,
                            label=f"true {typ}" + ("" if rel else " (not reliable)"))
    else:
        ax.plot(BR[m & ~RELIABLE, 1], BR[m & ~RELIABLE, 0], "o", ms=4, color=CLR[reg],
                alpha=0.3, mfc="none", label="not reliable")
        ax.plot(BR[m & RELIABLE, 1], BR[m & RELIABLE, 0], "o", ms=5, color=CLR[reg], alpha=0.8,
                label="reliable")
    ax.plot([0, 1], [0, 1], color="k", lw=0.6, ls=":")
    ax.set(xlim=(0, 1.02), ylim=(0, 1.02), aspect="equal",
           xlabel="DIAGONAL profile breadth  (high = locked to kinematics)",
           ylabel="COLUMN profile breadth  (high = locked to neural time)",
           title=f"{reg}  (n={int(m.sum())} defined maps, reliable {int((m & RELIABLE).sum())})")
    ax.text(0.03, 0.97, "locked to\nneural time", transform=ax.transAxes, va="top", fontsize=8, color="0.4")
    ax.text(0.97, 0.03, "locked to\nkinematics", transform=ax.transAxes, ha="right", fontsize=8, color="0.4")
    ax.text(0.03, 0.03, "tiling", transform=ax.transAxes, fontsize=8, color="0.4")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=7, loc="upper right")
fig.suptitle("STEP 6 / breadth of each neuron's two profiles (half B; participation ratio of the "
             "positive excess MI: 1 = spread over every kinematic window)", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.94)); _save(fig, "step6_breadth")
if TRUE_TYPE is not None:
    print(f"\nSYNTHETIC check -- by TRUE type: reliable fraction, and median breadth "
          f"(column / diagonal) of the reliable ones")
    for typ in TYPE_STYLE:
        mm = TRUE_TYPE == typ
        if not mm.any(): continue
        r_ = RELIABLE[mm]
        cb = np.nanmedian(BR[mm & RELIABLE, 0]) if (mm & RELIABLE).any() else np.nan
        db = np.nanmedian(BR[mm & RELIABLE, 1]) if (mm & RELIABLE).any() else np.nan
        print(f"  {typ:11s} n={int(mm.sum()):3d}  reliable {100*r_.mean():3.0f}%   "
              f"column {cb:.2f} / diagonal {db:.2f}")


# %% ── 7b. PCA check: do types show up without a breadth formula? ──────────
# 7b. cell --> a cross-check of step 6 with no formula chosen in advance. Each
# neuron is described by its two half-B profiles in noise units (as step 5b),
# each SORTED from high to low: that keeps "flat and high" vs "one peak" and
# drops WHERE the peak is (shifted stripes would otherwise need many PCs).
# PCA over the defined maps of both regions together (same axes for CP and CB);
# the null maps are projected onto the same PCs for reference.
def _desc(d, pi):
    out = []
    for w in ("colB", "diaB"):
        v = np.asarray(d[w], float) / _nsd(pi)
        out.append(np.nan_to_num(-np.sort(-v)))        # high to low, NaNs (missing) -> 0
    return np.concatenate(out)


_ok    = np.flatnonzero(~SPARSE)
D_REAL = np.stack([_desc(PROF[pi], pi) for pi in _ok])
D_NULL = np.stack([_desc(PNUL[pi], pi) for pi in _ok])
_mu    = D_REAL.mean(axis=0)
_U, _S, PCS = np.linalg.svd(D_REAL - _mu, full_matrices=False)
EXPL   = _S ** 2 / np.sum(_S ** 2)
SC_REAL, SC_NULL = (D_REAL - _mu) @ PCS[:2].T, (D_NULL - _mu) @ PCS[:2].T
print(f"\nPCA of the sorted profiles: PC1 {EXPL[0]:.0%}, PC2 {EXPL[1]:.0%} of the variance")

fig, axes = plt.subplots(1, 3, figsize=(16, 5), gridspec_kw=dict(width_ratios=[1.2, 1, 1]))
ax = axes[0]
ax.plot(SC_NULL[:, 0], SC_NULL[:, 1], "o", ms=3, color="0.8", alpha=0.6, label="null maps")
_rg, _rel = regions_run[_ok], RELIABLE[_ok]
if TRUE_TYPE is not None:
    _tt = TRUE_TYPE[_ok]
    for typ, (col_, mk) in TYPE_STYLE.items():
        for rel, fill in ((True, None), (False, "none")):
            mm = (_tt == typ) & (_rel == rel)
            if mm.any():
                ax.plot(SC_REAL[mm, 0], SC_REAL[mm, 1], mk, ms=6, color=col_, mfc=fill, alpha=0.85,
                        label=f"true {typ}" + ("" if rel else " (not reliable)"))
else:
    for reg in ("CP", "CB"):
        for rel, fill in ((True, None), (False, "none")):
            mm = (_rg == reg) & (_rel == rel)
            if mm.any():
                ax.plot(SC_REAL[mm, 0], SC_REAL[mm, 1], "o", ms=5, color=CLR[reg], mfc=fill,
                        alpha=0.8 if rel else 0.4, label=f"{reg}" + ("" if rel else " (not reliable)"))
ax.set(xlabel=f"PC1 ({EXPL[0]:.0%})", ylabel=f"PC2 ({EXPL[1]:.0%})",
       title="each neuron's sorted profiles, projected on PC1/PC2")
ax.title.set_fontsize(9); ax.legend(fontsize=7); ax.spines[["top", "right"]].set_visible(False)
# what PC1 and PC2 mean: their loadings on the two sorted profiles
_rank = np.arange(1, n_kin + 1)
for c, k in enumerate((0, 1), start=1):
    ax = axes[c]
    ax.plot(_rank, PCS[k, :n_kin], color="k", lw=2, label="COLUMN profile part")
    ax.plot(_rank, PCS[k, n_kin:], color="k", lw=2, ls="--", label="DIAGONAL profile part")
    ax.axhline(0, color="k", lw=0.5)
    ax.set(xlabel="rank within the profile (1 = its highest window)", ylabel="loading",
           title=f"PC{k + 1} loadings  (+ score = more of this shape)")
    ax.title.set_fontsize(9); ax.legend(fontsize=7); ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("STEP 6b / PCA of the sorted profiles (half B, noise units): do the map types separate "
             "without a breadth formula?", fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.93)); _save(fig, "step6b_pca_sorted_profiles")
if TRUE_TYPE is not None:
    print("SYNTHETIC check -- median PC1 / PC2 score by TRUE type (defined maps):")
    for typ in TYPE_STYLE:
        mm = TRUE_TYPE[_ok] == typ
        if mm.any():
            print(f"  {typ:11s} PC1 {np.median(SC_REAL[mm, 0]):+6.2f}   PC2 {np.median(SC_REAL[mm, 1]):+6.2f}")


# %% ── 8. save ──────────────────────────────────────────────────────────────
if SAVE_VARS and not SYNTHETIC:
    p = config.var_path(f"{RUN_NAME}_{FEAT_TAG}_{NW_TAG}_ks{KIN_STEP}_cat{CAT_TAG}"
                        + (f"_{SEL_TAG}" if SEL_TAG else "")
                        + f"_n{'all' if N_TEST is None else N_TEST}_results.pkl")
    with open(p, "wb") as fh:
        pickle.dump(dict(key=KEY, uids=uids, regions=regions_run,
                         kin_centers=kin_centers, win_centers=win_centers,
                         rel_stat=REL_STAT, stats=STATS, r_real=R_REAL, r_null=R_NULL,
                         threshold=THR, reliable=RELIABLE,
                         sparse=SPARSE, finite=FINITE, synthetic=SYNTHETIC, true_type=TRUE_TYPE,
                         snr=SNR,
                         profiles=PROF, null_profiles=PNUL,
                         breadth=BR, breadth_null=BN), fh, protocol=4)
    print(f"saved -> {p}")

# %%
