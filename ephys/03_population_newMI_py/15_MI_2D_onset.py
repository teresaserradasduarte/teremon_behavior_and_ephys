# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# TIME x TIME multi-dimensional MI, ALIGNED TO REACH ONSET.
#
# This is 12_MI_2D_kin_x_neuron.py with one change: both axes are relative to
# each trial's own REACH ONSET instead of the max-reach frame. Everything after
# the kinematic stack (cells 3-7) is the same code, so the two scripts can be
# compared panel for panel.
#
#   rows    = kinematic window centre re ONSET   (when the MOVEMENT is read)
#   columns = neural window centre re ONSET      (when the SPIKES are counted)
#
# WHAT t=0 MEANS HERE, and why this script exists:
#   st_reach t=0 is the MAX-REACH frame, index 150 of the 201 tracked frames
#   (step1_find_reaches..m: reach_goBack_frames=150, max_reach=151 MATLAB-1-idx;
#   confirmed by tm_w in the behaviour mat, which runs -375..+125 ms with zero at
#   index 150). 09/12 instead use frame 147 as the kinematic zero, which is 10 ms
#   EARLIER than the true max-reach frame. Here both axes use the onset, taken
#   per trial, so that offset cannot leak in.
#
# ONSET comes from the MATLAB behaviour pipeline, not re-detected here:
#   reaches.start_forw.start_forw   start of the forward movement (frame index
#                                   into the same 201-frame window, 1-based)
#   mapped onto the ephys reaches with bhv.reach_inVec_idx (774 -> 565 mask;
#   verified: reach_mat[mask] equals REACH_PX exactly).
#   Measured on 1_CoteDor R2: median -107 ms re max reach, 5-95% -202..-75 ms,
#   0% missing. The alternatives (start_realStop 32% missing, start_stopX) are
#   selectable through ONSET_VAR.
#
# CONSEQUENCE OF PER-TRIAL ALIGNMENT: a kinematic window near the edge of
# KIN_RANGE exists for some trials and not for others, because the tracked trace
# is a fixed 201 frames around max reach while the onset moves. Each kinematic
# window is therefore estimated on the trials that HAVE it (see cell 3), and cell
# 1c prints and plots the coverage per window. A global "drop any trial with a
# missing window" rule would throw away every late-onset trial, which is a biased
# sample, not a conservative one.
import os
import re
import time
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import h5py
import scipy.io as sio
from scipy.signal import savgol_filter
from scipy.stats import mannwhitneyu, wilcoxon
from joblib import Parallel, delayed
import mi_io, config, mi_rates

# ── timing: every expensive step records itself, summarised at the end ───────
TIMES = {}
class _step:
    """with _step("name"): ...   -> records and prints the elapsed time."""
    def __init__(self, name): self.name = name
    def __enter__(self): self.t = time.time(); return self
    def __exit__(self, *a):
        dt = time.time() - self.t
        TIMES[self.name] = TIMES.get(self.name, 0.0) + dt
        print(f"    [{self.name}] {dt:.1f}s" if dt < 90 else
              f"    [{self.name}] {dt/60:.1f} min", flush=True)
        return False

# The output tree is Dropbox-synced: Dropbox briefly locks a file while it
# uploads it, and overwriting one then raises OSError 22. Retry, do not crash.
SAVE_RETRIES = {"n": 0}
def _savefig_retry(fig, path, tries=6, **kw):
    for t in range(tries):
        try:
            return fig.savefig(path, **kw)
        except OSError:
            if t == tries - 1: raise
            SAVE_RETRIES["n"] += 1
            time.sleep(0.5 * (t + 1))

# ── the kinematic feature vector (joint MI is taken over ALL of these) ───────
FEATURES_2D = ("x", "y", "z", "speed", "dist_water")
# PC1/PC2 are NOT supported here: 12 refits them per kinematic window across
# trials, which needs every trial to share that window. With per-trial onsets
# they do not, so the basis would be fitted on a different subset per window.

# ── which reaches ────────────────────────────────────────────────────────────
REACH_CAT = (1, )      # reach categories (idx_reach_cat) to include.
                       # The onset is only well computed for CATEGORY 1, so this
                       # script refuses anything else unless you override below.
ONSET_CAT_ONLY = True  # True -> raise if REACH_CAT is not exactly (1,)
# ── which reaches WITHIN that category ───────────────────────────────────────
# Outcome filter on top of REACH_CAT. None = every reach of the category, which
# is what every run so far did.
#   {"hit": 1}            only reaches that hit
#   {"hit": 0}            only misses
#   {"succ": 1}           only successful reaches
#   {"succ": 0}           only failures
#   {"succ": 1, "hit": 1} both at once (AND)
#   {"LCR": (1, 3)}       a set of values (left or right spout)
# Keys map onto the per-reach labels: hit -> idx_reach_hit, succ ->
# idx_reach_succ, LCR -> idx_reach_LCR, PP -> idx_reach_PP, cat -> idx_reach_cat.
# A reach whose label is NaN never matches, so it is dropped by any filter that
# names that label -- cell 3 prints how many that is.
# The choice goes into the cache key, the cache filename and the figure folder,
# so a filtered run cannot overwrite or silently load an unfiltered one.
TRIAL_SELECT = None


# ── ONSET ────────────────────────────────────────────────────────────────────
ONSET_VAR   = "start_forw"      # "start_forw" | "start_realStop" | "start_stopX"
ONSET_CACHE_MODE = "rerun"       # "load" = reuse the onset pickle if it matches
                                # "rerun" = re-read the behaviour mats

# ── the KINEMATIC / BEHAVIOUR window (rows), now re ONSET ───────────────────
KIN_WIN   = 21         # frames per kinematic window (21 @ 400 fps = 52.5 ms)
KIN_STEP  = 4          # frames between kinematic windows (4 = 10 ms)
KIN_RANGE = (-0.15, 0.20)       # (start, stop) in SECONDS re ONSET, for the
                       # window CENTRES. The tracked trace is only 201 frames
                       # around max reach, so the far end of this range exists
                       # for a shrinking fraction of trials -- cell 1c shows it.

# ── the NEURAL window (columns), re ONSET ───────────────────────────────────
NEU_WIN_WIDTH = 0.07  # (s) width of the spike-count window
NEU_WIN_STEP  = 0.025  # (s) step between neural windows
NEU_WIN_RANGE = (-0.5, 0.5)   # (s) re ONSET, for the window CENTRES

# ── how much to run ──────────────────────────────────────────────────────────

# Minimum reaches a neuron needs to be estimated at all. Set here, and it wins
# over 09's value (this script's settings are restored after 09's head is
# exec'd). It matters most for a filtered run (TRIAL_SELECT), which has far
# fewer reaches per session, and it decides how many sessions survive.
# Measured for cat 1 + succ==1:
#   20 -> 25/29 sessions, 2542 neurons | 25 -> 23/29, 2429 | 30 -> 22/29, 2245
#   40 -> 18/29, 1831    | 50 -> 13/29, 1273
# The excess MI stays unbiased even at 15 trials (the per-cell null removes the
# bias); what shrinks is sensitivity, so this trades neurons against how much of
# a real effect each neuron recovers. It is part of the cache key.
MIN_TRIALS  = 20

N_TEST      = None     # neurons PER REGION. None = whole population.
N_JOBS      = -1
DO_NULL     = True     # per-CELL shuffle null -> Excess MI (obs - null)
N_SHUF_2D   = 20       # permutations per cell
MIN_WIN_TRIALS = None  # minimum trials for a kinematic window to be estimated;
                       # None = MIN_TRIALS inherited from 09 (20)

# ── DIAGONAL vs NEURAL TIME ──────────────────────────────────────────────────
# Identical to 12: f(neural) = the map averaged down each COLUMN, h(lag) = along
# each DIAGONAL; each is built on half of the kinematic windows and scored on the
# other half, row-centred. dR2 = R2_lag - R2_neural. Read it as a population
# tendency: a noisy map and a ridge whose slope is not 1 both land near 0.
LAG_BIN    = None      # (s) lag resolution; None = the neural step
CV_GAP     = None      # kinematic windows dropped from the building half next to
                       # the split. None = the overlap, KIN_WIN // KIN_STEP
FIG_SUBSET = "all"     # "all" or "structured"

# ── caching ──────────────────────────────────────────────────────────────────
CACHE_MODE  = "load"   # "load" = reuse a cached run whose parameters MATCH
SAVE_CACHE  = True

# ── what the figures report ──────────────────────────────────────────────────
N_EXAMPLES   = 40      # individual landscapes PER REGION; None = every neuron
EXAMPLE_MODE = "spread"     # "spread" | "top" | "random"
PER_MOUSE_FILES = True      # also write one mean-landscape figure per animal
PER_MOUSE_MIN_N = 5         # skip a region/mouse panel with fewer neurons
EQUAL_AXES  = False
ALIGN_FIG_RANGE = (-0.25, 0.25)   # (s) re onset, the span of the check figure

SHOW_FIGS   = True
SAVE_FIGS   = True
SAVE_VARS   = True
CMAP        = "viridis"
OUT_SUBDIR  = "MI_2D_kin_x_neuron_ONSET_v2"
RUN_NAME    = "td2Donset_v2"


# %% ── 1. inherit 09's parameters, data and estimator ────────────────────────
# 1. cell --> same trick as 11/12: exec 09's cells 0-2 with FEATURES / N_TEST /
# CAT_KEEP overridden, then put this script's own settings back.
src  = open("09_MI_multiD_kin.py", encoding="utf-8").read()
_m3  = re.search(r"(?m)^# %% ── 3\.", src)
head = src[:_m3.start()]
head, n1 = re.subn(r"(?m)^FEATURES\s*=.*$", f"FEATURES    = {tuple(FEATURES_2D)!r}",
                   head, count=1)
head, n2 = re.subn(r"(?m)^N_TEST\s*=.*$",   f"N_TEST      = {N_TEST!r}", head, count=1)
head, n3 = re.subn(r"(?m)^CAT_KEEP\s*=.*$", f"CAT_KEEP    = {set(REACH_CAT)!r}",
                   head, count=1)
if not (n1 and n2 and n3):
    raise RuntimeError("could not find FEATURES / N_TEST / CAT_KEEP in 09 -- renamed?")
if ONSET_CAT_ONLY and set(REACH_CAT) != {1}:
    raise ValueError(
        f"REACH_CAT={REACH_CAT}: {ONSET_VAR} is only well computed for reach "
        f"category 1. Set ONSET_CAT_ONLY=False only if you have checked the "
        f"onset for the other categories yourself.")
if any(f in ("PC1", "PC2") for f in FEATURES_2D):
    raise NotImplementedError(
        "PC1/PC2 need one basis per kinematic window fitted across ALL trials; "
        "with per-trial onsets the trials differ per window. Use raw coordinates.")

_MINE = {k: globals()[k] for k in
         ("OUT_SUBDIR", "RUN_NAME", "SAVE_FIGS", "SAVE_VARS", "SHOW_FIGS", "N_JOBS",
          "MIN_TRIALS")}
with _step("load data + 09 setup"):
    # 09's cells print a description of 09's OWN run (its windows, its output
    # folder, its PCA, its neuron count). None of that describes this script, so
    # it is not shown: everything printed below states what THIS run uses. If 09
    # fails, its output is printed, so the error is not hidden.
    import io as _io, contextlib as _ctx
    _buf = _io.StringIO()
    try:
        with _ctx.redirect_stdout(_buf):
            exec(compile(head, "09_head", "exec"))
    except BaseException:
        print(_buf.getvalue()); raise
globals().update(_MINE)

# NOT a parameter: the reaches were cut around this frame in MATLAB, and
# st_reach t=0 is that same frame. It comes from config with the data, and the
# onset loader below re-checks it against tm_w in every session's behaviour mat.
ALIGN_FRAME = config.REACH_ENDPOINT_FRAME


# ── trial selection (TRIAL_SELECT) ───────────────────────────────────────────
_LABEL_KEY = {"hit": "idx_reach_hit", "succ": "idx_reach_succ",
              "LCR": "idx_reach_LCR", "PP": "idx_reach_PP", "cat": "idx_reach_cat"}

def _select_mask(labels, select=None):
    """-> boolean mask over a session's reaches for the TRIAL_SELECT dict.
    NaN labels never match: np.isin is False for them, so a reach with an
    unknown outcome is excluded whenever the filter names that label."""
    select = TRIAL_SELECT if select is None else select
    n = len(np.asarray(labels["idx_reach_cat"]).ravel())
    m = np.ones(n, bool)
    if not select: return m
    for k, v in select.items():
        if k not in _LABEL_KEY:
            raise KeyError(f"TRIAL_SELECT key {k!r} unknown; use {sorted(_LABEL_KEY)}")
        a = np.asarray(labels[_LABEL_KEY[k]], float).ravel()
        m &= np.isin(a, np.asarray(np.atleast_1d(v), float))
    return m

SEL_TAG = "" if not TRIAL_SELECT else "-".join(
    f"{k}{'+'.join(str(x) for x in np.atleast_1d(v))}"
    for k, v in sorted(TRIAL_SELECT.items()))
if SEL_TAG:
    print(f"TRIAL_SELECT {TRIAL_SELECT} -> tag {SEL_TAG!r}")

CAT_TAG   = "".join(str(c) for c in sorted(REACH_CAT))
OUT_DIR_C = f"{OUT_SUBDIR}_cat{CAT_TAG}" + (f"_{SEL_TAG}" if SEL_TAG else "")
MIN_WIN_N = MIN_TRIALS if MIN_WIN_TRIALS is None else MIN_WIN_TRIALS

# --- NEURAL axis: this script's grid, now relative to ONSET -------------------
WIN_WIDTH   = NEU_WIN_WIDTH
# The neural window WIDTH goes into every output name (figure folder, cache,
# final pickle), built from the value in use, so runs with different widths are
# different folders and files and cannot overwrite each other (same as 12).
NW_TAG    = f"nw{int(round(WIN_WIDTH * 1000))}"
OUT_DIR_C = (f"{OUT_SUBDIR}_{NW_TAG}_cat{CAT_TAG}" + (f"_{SEL_TAG}" if SEL_TAG else ""))
print(f"outputs for THIS run -> {OUT_DIR_C}/   (neural window {WIN_WIDTH*1000:.0f} ms)")
WIN_STEP    = NEU_WIN_STEP
WIN_RANGE   = tuple(NEU_WIN_RANGE)
win_starts  = np.arange(WIN_RANGE[0], WIN_RANGE[1] - WIN_WIDTH + WIN_STEP/2, WIN_STEP)
win_centers = win_starts + WIN_WIDTH / 2
n_wins      = len(win_starts)

# --- KINEMATIC axis: window centres in FRAMES relative to each trial's onset --
kin_off  = np.arange(int(round(KIN_RANGE[0]*FPS)), int(round(KIN_RANGE[1]*FPS)) + 1,
                     KIN_STEP)                      # centre offsets, frames re onset
kin_centers = kin_off / FPS                          # s re onset
kin_start_off = kin_off - KIN_WIN // 2               # window START offsets, frames
n_kin    = len(kin_off)

print(f"\n2D landscape ALIGNED TO ONSET ({ONSET_VAR}): {n_kin} kinematic windows "
      f"x {n_wins} neural windows = {n_kin*n_wins} cells per neuron")
print(f"  KIN : win {KIN_WIN} fr ({KIN_WIN/FPS*1000:.1f} ms), step {KIN_STEP} fr "
      f"({KIN_STEP/FPS*1000:.1f} ms), centres "
      f"{kin_centers[0]*1000:.0f}..{kin_centers[-1]*1000:.0f} ms re onset")
print(f"  NEUR: win {WIN_WIDTH*1000:.0f} ms, step {WIN_STEP*1000:.0f} ms, centres "
      f"{win_centers[0]*1000:.0f}..{win_centers[-1]*1000:.0f} ms re onset")
print(f"features: {list(FEATURES_2D)} (joint KSG, k={K_KSG})  |  neural variable: "
      f"{NEURAL_VAR}  |  neurons: {len(run_idx)}  |  reach cat {sorted(CAT_KEEP)}"
      + (f"  |  trials: {SEL_TAG}" if SEL_TAG else ""))


# %% ── 1b. reach onset per trial, from the MATLAB behaviour pipeline ─────────
# 1b. cell --> ONSET_FRAME[key] = 0-based frame index of the onset within the 201
# tracked frames; ONSET_T[key] = the same in seconds relative to MAX REACH, which
# is what the spike times are aligned to, so shifting spikes by it puts them on
# the onset clock.
_ON_PKL = config.var_path(f"reach_onsets_{ONSET_VAR}.pkl")
ONSET_FRAME, ONSET_T = {}, {}
_loaded = False
if ONSET_CACHE_MODE == "load" and _ON_PKL.is_file():
    with open(_ON_PKL, "rb") as f: _c = pickle.load(f)
    if (_c.get("onset_var") == ONSET_VAR and _c.get("align_frame") == ALIGN_FRAME
            and set(_c["frames"]) == set(REACH_MM)):
        ONSET_FRAME = _c["frames"]; ONSET_T = _c["t"]; _loaded = True
        print(f"loaded onsets <- {_ON_PKL.name}")
if not _loaded:
    with _step("read onsets from the behaviour mats"):
        for key in REACH_MM:
            mouse, sess = key
            b = sio.loadmat(config.bhv_fundamentals_mat(mouse, sess),
                            squeeze_me=True, struct_as_record=False)["bhv"]
            sel = np.asarray(b.reach_inVec_idx).ravel().astype(bool)
            with h5py.File(config.behav_session_mat(mouse, sess), "r") as f:
                st = np.asarray(f["reaches"][ONSET_VAR][ONSET_VAR][:], float).ravel()
                tmw = np.asarray(f["reaches"]["tm_w"][:], float).ravel()
            if sel.size != st.size:
                raise ValueError(f"{key}: reach_inVec_idx {sel.size} vs "
                                 f"{ONSET_VAR} {st.size}")
            if sel.sum() != REACH_MM[key].shape[0]:
                raise ValueError(f"{key}: mask selects {sel.sum()} reaches but "
                                 f"REACH_MM has {REACH_MM[key].shape[0]}")
            if abs(tmw[ALIGN_FRAME]) > 1e-9:
                raise ValueError(f"{key}: tm_w[{ALIGN_FRAME}] = {tmw[ALIGN_FRAME]:.4f} s, "
                                 f"expected 0 -- this session was cut around a "
                                 f"different frame than config.REACH_ENDPOINT_FRAME")
            fr = st[sel] - 1.0                      # MATLAB 1-based -> python 0-based
            ONSET_FRAME[key] = fr
            ONSET_T[key]     = (fr - ALIGN_FRAME) / FPS      # s re max reach
    if SAVE_VARS:
        with open(_ON_PKL, "wb") as f:
            pickle.dump(dict(onset_var=ONSET_VAR, align_frame=ALIGN_FRAME,
                             frames=ONSET_FRAME, t=ONSET_T), f, protocol=4)
        print(f"onsets -> {_ON_PKL}")

_all_t = np.concatenate([v[np.isfinite(v)] for v in ONSET_T.values()])
_miss  = np.mean(~np.isfinite(np.concatenate(list(ONSET_T.values())))) * 100
print(f"{ONSET_VAR}: {_all_t.size} reaches with an onset ({_miss:.1f}% missing)")
print(f"  onset re max reach: median {np.median(_all_t)*1000:+.0f} ms, "
      f"5-95% {np.percentile(_all_t,5)*1000:+.0f}..{np.percentile(_all_t,95)*1000:+.0f} ms")


# %% ── 1c. IS THE ALIGNMENT ANY GOOD? ────────────────────────────────────────
# 1c. cell --> the check figure. If the onsets are right, the speed traces that
# look smeared when aligned to max reach should line up when aligned to onset,
# and the rise should start at 0. Also shows how many trials each kinematic
# window actually has, which is the cost of per-trial alignment.
_LBL = {}                                  # per session: the per-reach labels
for n in NEURONS:
    _LBL.setdefault((n["mouse"], n["sess"]), n["labels"])

SPD = {}                                   # per session: speed traces (trials, frames)
for key in REACH_MM:
    vel = savgol_filter(REACH_MM[key], SAVGOL_WIN, SAVGOL_ORD, deriv=1,
                        delta=1/FPS, axis=2)
    SPD[key] = np.linalg.norm(vel, axis=1)

_pre  = int(round(abs(ALIGN_FIG_RANGE[0]) * FPS))
_post = int(round(ALIGN_FIG_RANGE[1] * FPS))
_t_on = np.arange(-_pre, _post + 1) / FPS
_t_mx = (np.arange(REACH_MM[list(REACH_MM)[0]].shape[2]) - ALIGN_FRAME) / FPS

_MX, _ON, _ONS = [], [], []                # max-reach aligned / onset aligned
for key in REACH_MM:
    cat = _LBL[key]["idx_reach_cat"]
    keep_tr = (np.isin(cat, list(CAT_KEEP)) & np.isfinite(ONSET_FRAME[key])
               & _select_mask(_LBL[key]))
    S, fr   = SPD[key], ONSET_FRAME[key]
    for ti in np.flatnonzero(keep_tr):
        _MX.append(S[ti]); _ONS.append(ONSET_T[key][ti])
        a, b = int(fr[ti]) - _pre, int(fr[ti]) + _post + 1
        row = np.full(_t_on.size, np.nan)
        lo, hi = max(a, 0), min(b, S.shape[1])
        if hi > lo: row[lo - a:hi - a] = S[ti, lo:hi]
        _ON.append(row)
_MX, _ON, _ONS = np.array(_MX), np.array(_ON), np.array(_ONS)
print(f"alignment check on {_MX.shape[0]} reaches (cat {sorted(CAT_KEEP)})")

# coverage: how many trials have each kinematic window at all
_cov = np.zeros(n_kin)
_tot = 0
for key in REACH_MM:
    fr = ONSET_FRAME[key]; n_fr = REACH_MM[key].shape[2]
    cat = _LBL[key]["idx_reach_cat"]
    ok0 = np.isin(cat, list(CAT_KEEP)) & np.isfinite(fr) & _select_mask(_LBL[key])
    _tot += ok0.sum()
    for wi, off in enumerate(kin_start_off):
        s0 = fr[ok0] + off
        _cov[wi] += np.sum((s0 >= 0) & (s0 + KIN_WIN <= n_fr))
_cov_frac = _cov / max(_tot, 1)

_af = plt.figure(figsize=(17, 9))
gs  = _af.add_gridspec(2, 3, hspace=0.38, wspace=0.28)

ax = _af.add_subplot(gs[0, 0])
ax.hist(_ONS * 1000, bins=40, color="0.45", edgecolor="white")
ax.axvline(np.median(_ONS) * 1000, color="crimson", lw=2, ls="--",
           label=f"median {np.median(_ONS)*1000:+.0f} ms")
ax.axvline(0, color="k", lw=1.0, ls=":", label="max reach")
ax.set(xlabel="onset time re reach endpoint (ms)", ylabel="# reaches",
       title=f"{ONSET_VAR}: when the onset falls")
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

ax = _af.add_subplot(gs[0, 1])
ax.plot(_t_mx * 1000, np.nanmean(_MX, axis=0), color="0.25", lw=2)
ax.fill_between(_t_mx * 1000,
                np.nanmean(_MX, 0) - np.nanstd(_MX, 0) / np.sqrt(_MX.shape[0]),
                np.nanmean(_MX, 0) + np.nanstd(_MX, 0) / np.sqrt(_MX.shape[0]),
                color="0.25", alpha=0.25)
ax.axvline(0, color="k", lw=1.0, ls=":")
ax.set(xlabel="time re reach endpoint (ms)", ylabel="speed ± SEM (mm/s)",
       title="aligned to MAX REACH\n(onset is smeared across trials)")
ax.spines[["top", "right"]].set_visible(False)

ax = _af.add_subplot(gs[0, 2])
ax.plot(_t_on * 1000, np.nanmean(_ON, axis=0), color="crimson", lw=2)
_n_ok = np.sum(np.isfinite(_ON), axis=0)
ax.fill_between(_t_on * 1000,
                np.nanmean(_ON, 0) - np.nanstd(_ON, 0) / np.sqrt(np.maximum(_n_ok, 1)),
                np.nanmean(_ON, 0) + np.nanstd(_ON, 0) / np.sqrt(np.maximum(_n_ok, 1)),
                color="crimson", alpha=0.25)
ax.axvline(0, color="k", lw=1.0, ls=":", label="onset")
ax.set(xlabel="time re reach onset (ms)", ylabel="speed ± SEM (mm/s)",
       title="aligned to ONSET\n(the rise should start at 0)")
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

_ord = np.argsort(_ONS)
_vmax = np.nanpercentile(_MX, 99)
ax = _af.add_subplot(gs[1, 0])
ax.imshow(_MX[_ord], aspect="auto", origin="lower", cmap="magma", vmin=0, vmax=_vmax,
          extent=[_t_mx[0]*1000, _t_mx[-1]*1000, 0, _MX.shape[0]])
ax.plot(_ONS[_ord] * 1000, np.arange(_ONS.size), color="w", lw=1.0)
ax.axvline(0, color="w", lw=1.0, ls=":")
ax.set(xlabel="time re reach endpoint (ms)", ylabel="reach (sorted by onset)",
       title="speed, aligned to MAX REACH\n(white line = onset: it slopes)")

ax = _af.add_subplot(gs[1, 1])
ax.imshow(_ON[_ord], aspect="auto", origin="lower", cmap="magma", vmin=0, vmax=_vmax,
          extent=[_t_on[0]*1000, _t_on[-1]*1000, 0, _ON.shape[0]])
ax.axvline(0, color="w", lw=1.2, ls=":")
ax.set(xlabel="time re reach onset (ms)", ylabel="reach (sorted by onset)",
       title="speed, aligned to ONSET\n(the onset edge should be VERTICAL)")

ax = _af.add_subplot(gs[1, 2])
ax.plot(kin_centers * 1000, _cov_frac * 100, "o-", color="teal", lw=2, ms=4)
ax.axhline(100, color="k", lw=0.6, ls=":")
ax.set(xlabel="kinematic window centre re onset (ms)",
       ylabel="% of reaches with this window", ylim=(0, 105),
       title=f"coverage of KIN_RANGE {KIN_RANGE}\n"
             f"(the trace is a fixed 201 frames around max reach)")
ax.spines[["top", "right"]].set_visible(False)

_af.suptitle(f"Reach-onset alignment check  —  {ONSET_VAR}, cat {sorted(CAT_KEEP)}, "
             f"{_MX.shape[0]} reaches", fontsize=12)
fdir_chk = config.fig_dir(OUT_DIR_C)
if SAVE_FIGS:
    _savefig_retry(_af, fdir_chk / f"{RUN_NAME}_alignment_check.png", dpi=150,
                   bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(_af)

print(f"coverage of the kinematic windows: "
      f"{_cov_frac.min()*100:.0f}%..{_cov_frac.max()*100:.0f}% of reaches "
      f"(window centres {kin_centers[0]*1000:.0f}..{kin_centers[-1]*1000:.0f} ms)")
_thin = kin_centers[_cov_frac < 0.5]
if _thin.size:
    print(f"  !! under 50% coverage at {_thin[0]*1000:.0f}..{_thin[-1]*1000:.0f} ms "
          f"-- those windows rest on a minority of reaches, narrow KIN_RANGE")


# %% ── 2. kinematic features per ONSET-RELATIVE window, per session ──────────
# 2. cell --> KIN[key] = (n_trials, n_kin, n_feats), NaN where the window falls
# outside the tracked trace. Derivatives are computed on the FULL trace and then
# averaged inside each window, exactly as in 12.
NEEDS_WATER_2D = any(f in ("dist_water", "speed_water") for f in FEATURES_2D)

_t_kin = time.time()
KIN = {}
for key in REACH_MM:
    rp  = REACH_MM[key]
    spd = SPD[key]
    acc = savgol_filter(spd, SAVGOL_WIN, SAVGOL_ORD, deriv=1, delta=1/FPS, axis=1)
    cols = {"x": rp[:, 0], "y": rp[:, 1], "z": rp[:, 2], "speed": spd, "accel": acc}
    if NEEDS_WATER_2D:
        W = WATER[key]
        cols["dist_water"], cols["speed_water"] = W["dist_water"], W["speed_water"]
    fr   = ONSET_FRAME[key]
    n_tr, n_fr = rp.shape[0], rp.shape[2]
    out  = np.full((n_tr, n_kin, len(FEATURES_2D)), np.nan)
    for wi, off in enumerate(kin_start_off):
        s0  = fr + off
        okt = np.isfinite(s0) & (s0 >= 0) & (s0 + KIN_WIN <= n_fr)
        if not okt.any(): continue
        rows = np.flatnonzero(okt)
        idx  = s0[rows, None].astype(int) + np.arange(KIN_WIN)[None, :]
        for fi, f in enumerate(FEATURES_2D):
            out[rows, wi, fi] = np.take_along_axis(cols[f][rows], idx, axis=1).mean(1)
    KIN[key] = out
TIMES["kinematic stack"] = time.time() - _t_kin
print(f"kinematic stack built in {TIMES['kinematic stack']:.0f}s: "
      f"{KIN[list(KIN)[0]].shape} per session (trials, kin_wins, features)")

_cells = len(run_idx) * n_kin * n_wins * (1 + (N_SHUF_2D if DO_NULL else 0))
print(f"\nETA: ~{_cells:,} KSG calls; at ~2500 calls/s parallel that is "
      f"~{_cells/2500/60:.0f} min")


# %% ── 3. the 2D landscape, per neuron ───────────────────────────────────────
# 3. cell --> compute the MI landscape for every neuron, and the nulls if asked.
def _landscape2d_one(kin_stack, C_m, k, seed, n_shuf=0, min_n=MIN_WIN_N):
    """kin_stack (n_trials, n_kin, d), C_m (n_trials, n_neu)
    -> (mi, null_mean, null_std)

    DIFFERENT FROM 12: the kinematic windows are placed relative to each trial's
    own onset, so a window near the edge of KIN_RANGE exists for some trials and
    not others. Each kinematic window is therefore estimated on the trials that
    HAVE it, provided there are at least min_n of them. The permutations are
    drawn per window (the subset size changes) but reused across neural windows,
    so the null surface stays smooth along a row rather than speckled."""
    n_kin_, n_neu_ = kin_stack.shape[1], C_m.shape[1]
    mi = np.full((n_kin_, n_neu_), np.nan)
    mu = np.full_like(mi, np.nan)
    sd = np.full_like(mi, np.nan)
    for ki in range(n_kin_):
        X = kin_stack[:, ki, :]
        rows = np.isfinite(X).all(axis=1)
        if rows.sum() < min_n: continue
        Xr = X[rows]
        rs = np.random.default_rng(seed + 7919 * ki)
        perms = [rs.permutation(Xr.shape[0]) for _ in range(n_shuf)]
        SH = np.full((n_shuf, n_neu_), np.nan) if n_shuf else None
        for wi in range(n_neu_):
            c = C_m[rows, wi]
            if np.unique(c).size < 2: continue
            mi[ki, wi] = _ksg_md(Xr, c, k, seed + ki*97 + wi)
            for s, p in enumerate(perms):
                SH[s, wi] = _ksg_md(Xr[p], c, k, seed + ki*97 + wi)
        if n_shuf:
            mu[ki] = np.nanmean(SH, axis=0)
            sd[ki] = np.nanstd(SH, axis=0, ddof=1)
    return mi, mu, sd


payload, keep = [], []
_t_pay = time.time()
for pi, i in enumerate(run_idx):
    n   = NEURONS[i]
    key = (n["mouse"], n["sess"])
    K, ons = KIN[key], ONSET_T[key]
    cat = n["labels"]["idx_reach_cat"]
    # a trial is usable when it has an onset, the right category, and at least
    # one kinematic window inside the tracked trace
    msk = (np.isfinite(ons) & np.isin(cat, list(CAT_KEEP))
           & np.isfinite(K).all(axis=2).any(axis=1)
           & _select_mask(n["labels"]))
    if msk.sum() < MIN_TRIALS: payload.append(None); continue
    # spikes onto the onset clock: st_reach is re MAX REACH, so subtract the
    # onset time of that same trial
    rows = np.flatnonzero(msk)
    spk  = [np.asarray(n["spikes"][t], float) - ons[t] for t in rows]
    C    = np.asarray(mi_rates.windowed(spk, win_starts, WIN_WIDTH, NEURAL_VAR), float)
    payload.append((K[msk], C)); keep.append(pi)
TIMES["spike counting + masks"] = time.time() - _t_pay
print(f"{len(keep)}/{len(run_idx)} neurons usable "
      f"({TIMES['spike counting + masks']:.0f}s)")
if SEL_TAG:
    _kept = sum(int(payload[pi][0].shape[0]) for pi in keep)
    print(f"TRIAL_SELECT {TRIAL_SELECT}: {_kept} reaches kept, summed over neurons")

# ── cache: this is the expensive step, so never repeat it for free ───────────
CACHE_KEY = dict(features=list(FEATURES_2D), KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP,
                 KIN_RANGE=list(KIN_RANGE), n_kin=n_kin, n_wins=n_wins,
                 N_TEST=N_TEST, N_SHUF_2D=(N_SHUF_2D if DO_NULL else 0),
                 WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP, WIN_RANGE=list(WIN_RANGE),
                 K_KSG=K_KSG, NEURAL_VAR=NEURAL_VAR, CAT_KEEP=sorted(CAT_KEEP),
                 MIN_TRIALS=MIN_TRIALS, MIN_WIN_N=MIN_WIN_N, n_neurons=len(run_idx),
                 align="onset", onset_var=ONSET_VAR, align_frame=ALIGN_FRAME,
                 v="onset", **({"select": SEL_TAG} if SEL_TAG else {}))
CACHE_PATH = config.var_path(
    f"{RUN_NAME}_{'-'.join(FEATURES_2D)}_{NW_TAG}_kw{KIN_WIN}_ks{KIN_STEP}"
    f"_n{'all' if N_TEST is None else N_TEST}_s{N_SHUF_2D}_cat{CAT_TAG}"
    f"{'_' + SEL_TAG if SEL_TAG else ''}_{ONSET_VAR}_2Dcache.pkl")

_cached = None
if CACHE_MODE == "load" and CACHE_PATH.is_file():
    with _step("cache read"):
        with open(CACHE_PATH, "rb") as f: _c = pickle.load(f)
    if _c.get("key") == CACHE_KEY:
        _cached = _c
        print(f"LOADED cache {CACHE_PATH.name} "
              f"({CACHE_PATH.stat().st_size/1e6:.0f} MB) -- no recompute")
    else:
        diff = [k for k in CACHE_KEY if _c.get("key", {}).get(k) != CACHE_KEY[k]]
        print(f"cache exists but parameters differ ({', '.join(diff)}) -> recomputing")

if _cached is not None:
    MI2, NU2, NSD2 = _cached["mi"], _cached["null"], _cached["null_std"]
    keep = list(_cached["keep"])
else:
    _t_ls = time.time()
    res = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
        delayed(_landscape2d_one)(payload[pi][0], payload[pi][1], K_KSG, 100+pi,
                                  N_SHUF_2D if DO_NULL else 0)
        for pi in keep)
    MI2  = np.full((len(run_idx), n_kin, n_wins), np.nan)
    NU2  = np.full_like(MI2, np.nan)
    NSD2 = np.full_like(MI2, np.nan)
    for pi, (m, mu, sd) in zip(keep, res):
        MI2[pi] = m; NU2[pi] = mu; NSD2[pi] = sd
    TIMES["2D landscapes (MI + nulls)"] = time.time() - _t_ls
    print(f"2D landscapes done in {TIMES['2D landscapes (MI + nulls)']/60:.1f} min")
    if SAVE_CACHE:
      with _step("cache write"):
        with open(CACHE_PATH, "wb") as f:
            pickle.dump(dict(key=CACHE_KEY, mi=MI2, null=NU2, null_std=NSD2,
                             keep=keep, uids=[uids[i] for i in run_idx],
                             regions=regions_run), f, protocol=4)
        print(f"cache -> {CACHE_PATH} ({CACHE_PATH.stat().st_size/1e6:.0f} MB)")

EXC2 = MI2 - NU2 if DO_NULL else MI2


# %% ── 4. mean 2D landscape per region ───────────────────────────────────────
# 4. cell --> average the 2D landscapes over all neurons in each region, and plot.
ASPECT = "equal" if EQUAL_AXES else "auto"
fdir   = config.fig_dir(OUT_DIR_C)
extent = [win_centers[0], win_centers[-1], kin_centers[0], kin_centers[-1]]
lab    = (f"KSG k={K_KSG}, {NEURAL_VAR}, joint [{', '.join(FEATURES_2D)}], "
          f"cat {sorted(CAT_KEEP)}{', ' + SEL_TAG if SEL_TAG else ''}, "
          f"aligned to {ONSET_VAR}{' — EXCESS' if DO_NULL else ''}")

def _draw2d(ax, M, title, vmin=None, vmax=None):
    im = ax.imshow(M, origin="lower", aspect=ASPECT, extent=extent, cmap=CMAP,
                   vmin=vmin, vmax=vmax)
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(d, d, color="w", lw=1.0, ls="--", alpha=0.7)     # simultaneity
    ax.axvline(0, color="w", lw=0.8, ls=":", alpha=0.7)
    ax.axhline(0, color="w", lw=0.8, ls=":", alpha=0.7)
    ax.set(xlabel="NEURAL window centre re reach onset (s)",
           ylabel="KINEMATIC window centre re reach onset (s)", title=title)
    return im

MEAN2 = {}
for reg in ("CP", "CB"):
    idx = [pi for pi in keep if regions_run[pi] == reg]
    MEAN2[reg] = np.nanmean(EXC2[idx], axis=0) if idx else None

vals = np.concatenate([m.ravel() for m in MEAN2.values() if m is not None])
vmax = np.nanpercentile(vals, 99); vmin = min(0.0, np.nanpercentile(vals, 1))
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, reg in zip(axes, ("CP", "CB")):
    if MEAN2[reg] is None: ax.set_visible(False); continue
    n_ = len([pi for pi in keep if regions_run[pi] == reg])
    im = _draw2d(ax, MEAN2[reg], f"{reg} (n={n_})", vmin, vmax)
    plt.colorbar(im, ax=ax, label="mean MI (bits)")
plt.suptitle(f"{lab}  —  mean 2D landscape\n"
             f"dashed = simultaneity; ABOVE it the kinematics happen AFTER the "
             f"spikes (predictive), BELOW it the spikes lag the movement",
             fontsize=11)
plt.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / f"{RUN_NAME}_mean2D.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

# CB - CP difference map: where do the two regions actually differ?
if all(MEAN2[r] is not None for r in ("CP", "CB")):
    D = MEAN2["CB"] - MEAN2["CP"]
    lim = np.nanpercentile(np.abs(D), 99)
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(D, origin="lower", aspect=ASPECT, extent=extent, cmap="RdBu_r",
                   norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim))
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(d, d, color="k", lw=1.0, ls="--", alpha=0.6)
    ax.axvline(0, color="k", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set(xlabel="NEURAL window centre re reach onset (s)",
           ylabel="KINEMATIC window centre re reach onset (s)",
           title=f"CB - CP  ({lab})")
    plt.colorbar(im, ax=ax, label="Δ mean MI (bits)")
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(fdir / f"{RUN_NAME}_diff2D.png", dpi=150,
                              bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 5. DIAGONAL vs NEURAL TIME ──────────────────────────────────────────
# 5. cell --> per neuron: does the map vary along the DIAGONAL (lag) or along the
# NEURAL axis? Both descriptions are built on one half of the kinematic windows
# and scored on the other half, row-centred. See the block in cell 0.
_lag_bin = WIN_STEP if LAG_BIN is None else LAG_BIN
_gap     = KIN_WIN // KIN_STEP if CV_GAP is None else CV_GAP
if n_kin // 2 - _gap < 2:
    raise ValueError(f"only {n_kin} kinematic windows: with CV_GAP={_gap} each fold "
                     f"trains on {n_kin//2 - _gap}. Widen KIN_RANGE or lower CV_GAP.")

LAG      = win_centers[None, :] - kin_centers[:, None]    # neural - kinematic (s)
LAG_IDX  = np.round(LAG / _lag_bin).astype(int)
lag_axis = np.arange(LAG_IDX.min(), LAG_IDX.max() + 1) * _lag_bin
LAG_IDX -= LAG_IDX.min()
n_lag    = lag_axis.size
# Lags that EVERY kinematic window reaches. Outside this range a diagonal average
# is built from only some kinematic windows (the corners of the map), so the ends
# of h(lag) are thinner -- drawn shaded in the figures.
lag_full = ((lag_axis >= LAG.min(axis=1).max() - _lag_bin / 2) &
            (lag_axis <= LAG.max(axis=1).min() + _lag_bin / 2))


def _rowc(M):
    """subtract each kinematic row's own mean, ignoring NaN. This is what makes
    the comparison independent of how strong the map is at that kinematic time."""
    ok = np.isfinite(M); ct = ok.sum(axis=1, keepdims=True)
    mu = np.where(ct > 0,
                  np.where(ok, M, 0).sum(axis=1, keepdims=True) / np.maximum(ct, 1),
                  np.nan)
    return M - mu


def _lag_profile(E, rows=None):
    """h(lag): mean of the map along each diagonal (optionally only some rows)."""
    Er = E if rows is None else E[rows]
    Lr = LAG_IDX if rows is None else LAG_IDX[rows]
    ok = np.isfinite(Er)
    sm = np.bincount(Lr[ok], Er[ok], n_lag)
    ct = np.bincount(Lr[ok], minlength=n_lag)
    return np.where(ct > 0, sm / np.maximum(ct, 1), np.nan)


def _neu_profile(E, rows=None):
    """f(neural): mean of the map down each column (optionally only some rows)."""
    Er = E if rows is None else E[rows]
    ok = np.isfinite(Er)
    ct = ok.sum(axis=0)
    return np.where(ct > 0, np.where(ok, Er, 0).sum(axis=0) / np.maximum(ct, 1), np.nan)


def _cv_axes(E):
    """-> (R2_neural, R2_lag) on HELD-OUT kinematic windows.
    Fold 1 builds both profiles from the first half of the kinematic windows and
    scores them on the second; fold 2 the reverse. _gap windows next to the split
    are left out of the building half. Both are scored on exactly the same cells,
    and data and predictions alike are ROW-CENTRED first (see _rowc), so neither
    is asked to predict how strong the map is at that kinematic time."""
    h = n_kin // 2
    folds = [(np.arange(0, h - _gap), np.arange(h, n_kin)),
             (np.arange(h + _gap, n_kin), np.arange(0, h))]
    sse_n = sse_l = sst = 0.0
    for tr, te in folds:
        pn = _rowc(np.broadcast_to(_neu_profile(E, tr), (te.size, n_wins)).copy())
        pl = _rowc(_lag_profile(E, tr)[LAG_IDX[te]])
        y  = _rowc(E[te])
        ok = np.isfinite(y) & np.isfinite(pn) & np.isfinite(pl)
        if ok.sum() < 5: continue
        y = y[ok]
        sst   += np.sum((y - y.mean()) ** 2)
        sse_n += np.sum((y - pn[ok]) ** 2)
        sse_l += np.sum((y - pl[ok]) ** 2)
    if sst <= 0: return np.nan, np.nan
    return 1 - sse_n / sst, 1 - sse_l / sst


R2N, R2L = {}, {}
for pi in keep:
    R2N[pi], R2L[pi] = _cv_axes(EXC2[pi])
DR2      = {pi: R2L[pi] - R2N[pi] for pi in keep}
STRUCT   = {pi: bool(np.nanmax([R2N[pi], R2L[pi]]) > 0) for pi in keep}
PROF_NEU = {pi: _neu_profile(EXC2[pi]) for pi in keep}      # bits, per neural window
PROF_LAG = {pi: _lag_profile(EXC2[pi]) for pi in keep}      # bits, per lag bin

print(f"\nlag bin {_lag_bin*1000:.0f} ms, {n_lag} lags "
      f"({lag_axis[0]*1000:.0f}..{lag_axis[-1]*1000:.0f} ms; every kinematic window "
      f"reaches {lag_axis[lag_full][0]*1000:.0f}..{lag_axis[lag_full][-1]*1000:.0f} ms)")
print(f"CV: 2 blocked folds over {n_kin} kinematic windows, gap {_gap} "
      f"-> each fold builds on {n_kin//2 - _gap}; scored row-centred")


def _med_iqr(v):
    q = np.percentile(v, [25, 75])
    return f"{np.median(v):+.3f} [{q[0]:+.2f},{q[1]:+.2f}]"


for sub_name, sub_ok in (("ALL neurons  <- the headline", lambda pi: True),
                         ("only those whose shape repeats (secondary, power-limited)",
                          lambda pi: STRUCT[pi])):
    print(f"\n{sub_name}")
    print(f"{'region':6s} {'n':>5s} {'R2_neural':>10s} {'R2_lag':>8s} "
          f"{'dR2 median [IQR]':>24s} {'% diagonal':>11s} {'dR2 vs 0':>10s}")
    print("-" * 80)
    _d = {}
    for reg in ("CP", "CB"):
        idx = [pi for pi in keep if regions_run[pi] == reg and sub_ok(pi)]
        d = np.array([DR2[pi] for pi in idx]); m = np.isfinite(d); d = d[m]
        _d[reg] = d
        if d.size == 0:
            print(f"{reg:6s} {0:5d}"); continue
        rn = np.array([R2N[pi] for pi in idx])[m]; rl = np.array([R2L[pi] for pi in idx])[m]
        p0 = wilcoxon(d).pvalue if d.size >= 10 else np.nan
        print(f"{reg:6s} {d.size:5d} {np.median(rn):+10.3f} {np.median(rl):+8.3f} "
              f"{_med_iqr(d):>24s} {100*np.mean(d > 0):10.0f}% {p0:10.1e}")
    if _d.get("CP", np.array([])).size and _d.get("CB", np.array([])).size:
        _u, _p = mannwhitneyu(_d["CP"], _d["CB"], alternative="two-sided")
        print(f"  CB - CP median dR2: {np.median(_d['CB']) - np.median(_d['CP']):+.3f}"
              f"   (Mann-Whitney p = {_p:.1e})")
print("\n  dR2 < 0: MI is set by WHEN the spikes are (vertical stripes)")
print("  dR2 > 0: MI is set by the LAG between spikes and movement (diagonal)")
print("  GRAIN OF SALT: a ridge whose slope is not 1 fits neither description, and")
print("  a noisy map fits neither either, so both land near dR2 = 0 for different")
print("  reasons. Read this as a population tendency, not a per-neuron verdict.")


# %% ── 5b. individual 2D landscapes ────────────────────────────────────────
# One figure per neuron: its map, f(neural) underneath, h(lag) below that, the
# kinematic marginal on the right. These are for LOOKING at the neuron; the dR2
# numbers are printed small at the bottom, deliberately not in the title, because
# on a single map they are weak (a ridge whose slope is not 1, or simply a noisy
# map, fits neither description).
# Ordered by dR2 so a run of them sweeps from diagonal-like to neural-time-like.
EX_DIR = config.fig_dir(OUT_DIR_C, "individual_landscapes")


def _shade_partial(ax, scale=1000.0):
    """grey where a diagonal average uses only SOME kinematic windows"""
    lo_, hi_ = lag_axis[lag_full][[0, -1]] * scale
    ax.axvspan(lag_axis[0] * scale, lo_, color="0.92", zorder=0, lw=0)
    ax.axvspan(hi_, lag_axis[-1] * scale, color="0.92", zorder=0, lw=0)


_rng_ex = np.random.default_rng(0)
examples = []
for reg in ("CP", "CB"):
    idx = [pi for pi in keep if regions_run[pi] == reg]
    if not idx: continue
    order = sorted(idx, key=lambda pi: -DR2[pi] if np.isfinite(DR2[pi]) else np.inf)
    if N_EXAMPLES is None:
        sel = order                                     # every neuron
    else:
        k_ = min(N_EXAMPLES, len(order))
        if EXAMPLE_MODE == "top":      sel = order[:k_]
        elif EXAMPLE_MODE == "random": sel = list(_rng_ex.choice(order, k_, replace=False))
        else:                          sel = [order[int(round(t))]
                                              for t in np.linspace(0, len(order)-1, k_)]
    for pi in sel:
        examples.append((reg, pi, order.index(pi) + 1, len(order)))
print(f"drawing {len(examples)} individual landscapes"
      f"{' (ALL neurons)' if N_EXAMPLES is None else ''}")

_t_ex = time.time()
for reg, pi, rank, ntot in examples:
    M   = EXC2[pi]
    uid = uids[run_idx[pi]]
    col = CLR[reg]
    fig = plt.figure(figsize=(11, 11))
    gs  = fig.add_gridspec(3, 2, width_ratios=[3, 1], height_ratios=[3, 1.1, 1.1],
                           wspace=0.5, hspace=0.35)
    axm = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[0, 1], sharey=axm)     # kinematic marginal (viewing)
    axb = fig.add_subplot(gs[1, 0], sharex=axm)     # f(neural), aligned with the map
    axl = fig.add_subplot(gs[2, 0])                 # h(lag), its own time axis

    lim = np.nanpercentile(np.abs(M), 99) if np.isfinite(M).any() else 1.0
    im = axm.imshow(M, origin="lower", aspect=ASPECT, extent=extent, cmap=CMAP,
                    vmin=min(0, np.nanpercentile(M, 1)), vmax=lim)
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    axm.plot(d, d, color="w", lw=1.0, ls="--", alpha=0.7)
    axm.axvline(0, color="w", lw=0.7, ls=":", alpha=0.6)
    axm.axhline(0, color="w", lw=0.7, ls=":", alpha=0.6)
    axm.set_ylabel("KINEMATIC window centre re reach onset (s)")
    axm.set_xlabel("NEURAL window centre re reach onset (s)")

    # kinematic marginal: the map averaged ACROSS each row (over neural windows).
    _ok = np.isfinite(M)
    g_kin = np.where(_ok.any(axis=1),
                     np.where(_ok, M, 0).sum(axis=1) / np.maximum(_ok.sum(axis=1), 1),
                     np.nan)
    axr.plot(g_kin, kin_centers, color=col, lw=1.8)
    axr.axvline(0, color="k", lw=0.6); axr.axhline(0, color="k", lw=0.6, ls=":")
    axr.set_xlabel("mean over neural\nwindows (bits)")
    axr.set_title("kinematic marginal\n(viewing only)", fontsize=8)
    plt.setp(axr.get_yticklabels(), visible=False)
    fig.colorbar(im, cax=axr.inset_axes([1.12, 0, 0.1, 1]), label="MI (bits)")

    axb.plot(win_centers, PROF_NEU[pi], color=col, lw=1.8)
    axb.axhline(0, color="k", lw=0.6); axb.axvline(0, color="k", lw=0.6, ls=":")
    axb.set(xlabel="NEURAL window centre re reach onset (s)", ylabel="f(neural)\n(bits)",
            title="f(neural): averaged down each COLUMN")
    axb.title.set_fontsize(9)

    _shade_partial(axl, scale=1.0)
    axl.plot(lag_axis, PROF_LAG[pi], color=col, lw=1.8)
    axl.axhline(0, color="k", lw=0.6); axl.axvline(0, color="k", lw=0.6, ls=":")
    axl.set(xlabel="lag = neural - kinematic (s)     < 0: spikes BEFORE the movement",
            ylabel="h(lag)\n(bits)",
            title="h(lag): averaged along each DIAGONAL  "
                  "(grey: only some kinematic windows reach this lag)")
    axl.title.set_fontsize(9)
    for a_ in (axb, axl, axr): a_.spines[["top", "right"]].set_visible(False)

    fig.suptitle(f"{uid}  |  {reg}", fontsize=11)
    # secondary, small and grey: a single map does not settle this
    fig.text(0.5, 0.055,
             f"held-out R2  neural-time {R2N[pi]:+.2f} / lag {R2L[pi]:+.2f}"
             f"   ->  dR2 {DR2[pi]:+.2f}   (order {rank}/{ntot})"
             f"{'' if STRUCT[pi] else '   neither repeats across the kinematic halves'}",
             ha="center", va="top", fontsize=7.5, color="0.45")
    if SAVE_FIGS:
        _savefig_retry(fig, EX_DIR / f"{uid.replace('|','_')}.png", dpi=130,
                       bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

TIMES["individual figures"] = time.time() - _t_ex
print(f"{len(examples)} individual landscapes -> {EX_DIR} "
      f"({TIMES['individual figures']:.0f}s)")


# %% ── 6. DIAGONAL vs NEURAL TIME: the figure ────────────────────────────────
# Top    : the per-neuron comparison (dR2), its two ingredients, its distribution
# Bottom : WHAT the two profiles look like, averaged per region
_sub_ok  = (lambda pi: True) if FIG_SUBSET == "all" else (lambda pi: STRUCT[pi])
_sub_lab = ("all neurons" if FIG_SUBSET == "all"
            else "only neurons whose shape repeats")
_lf = plt.figure(figsize=(18, 10))
gs  = _lf.add_gridspec(2, 3, hspace=0.45, wspace=0.28)

# A: dR2 per neuron, as a proportion of each region
ax = _lf.add_subplot(gs[0, 0])
_bins = np.linspace(-1.5, 1.5, 41)
for reg in ("CP", "CB"):
    d = np.array([DR2[pi] for pi in keep if regions_run[pi] == reg and _sub_ok(pi)])
    d = d[np.isfinite(d)]
    if not d.size: continue
    ax.hist(np.clip(d, _bins[0], _bins[-1]), bins=_bins,
            weights=np.full(d.size, 1 / d.size), color=CLR[reg], alpha=0.55,
            edgecolor="white", label=f"{reg} (n={d.size})")
    ax.axvline(np.median(d), color=CLR[reg], lw=1.8, ls="--",
               label=f"{reg} median {np.median(d):+.2f}")
ax.axvline(0, color="k", lw=1.0, ls=":")
ax.set(xlabel="dR2 = R2_lag - R2_neural   (beyond ±1.5 piled in the end bins)",
       ylabel="proportion of neurons",
       title=f"diagonal vs neural time -- {_sub_lab}\n"
             f"< 0 neural-time-like   |   > 0 diagonal-like")
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=7)

# B: the two ingredients, every neuron
ax = _lf.add_subplot(gs[0, 1])
_lo, _hi = -1.0, 1.0
for reg in ("CP", "CB"):
    idx = [pi for pi in keep if regions_run[pi] == reg]
    x = np.clip([R2N[pi] for pi in idx], _lo, _hi)
    y = np.clip([R2L[pi] for pi in idx], _lo, _hi)
    ax.plot(x, y, "o", ms=3, alpha=0.3, color=CLR[reg], label=f"{reg} ({len(idx)})")
ax.plot([_lo, _hi], [_lo, _hi], color="k", lw=1.0, ls="--")
ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
ax.set(xlim=(_lo, _hi), ylim=(_lo, _hi), aspect="equal",
       xlabel="held-out R2, neural-time  f(neural)",
       ylabel="held-out R2, lag  h(lag)",
       title="the two ingredients, every neuron\nabove the dashed line = diagonal wins")
ax.legend(fontsize=7, loc="lower right")

# C: cumulative dR2
ax = _lf.add_subplot(gs[0, 2])
for reg in ("CP", "CB"):
    d = np.sort([DR2[pi] for pi in keep if regions_run[pi] == reg and _sub_ok(pi)
                 and np.isfinite(DR2[pi])])
    if d.size:
        ax.plot(d, np.arange(1, d.size + 1) / d.size, color=CLR[reg], lw=1.8,
                label=f"{reg} (n={d.size})")
ax.axvline(0, color="k", lw=1.0, ls=":")
ax.set(xlim=(-2, 1.5), xlabel="dR2", ylabel="cumulative proportion of neurons",
       title="same, cumulative\ncurve further RIGHT = more diagonal-like")
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=7)

# D: mean h(lag) per region
ax = _lf.add_subplot(gs[1, 0])
_shade_partial(ax)
for reg in ("CP", "CB"):
    idx = [pi for pi in keep if regions_run[pi] == reg and _sub_ok(pi)]
    if not idx: continue
    Pm = np.stack([PROF_LAG[pi] for pi in idx])
    mn = np.nanmean(Pm, axis=0)
    se = np.nanstd(Pm, axis=0) / np.sqrt(np.maximum(np.isfinite(Pm).sum(0), 1))
    ax.plot(lag_axis * 1000, mn, color=CLR[reg], lw=2, label=f"{reg} (n={len(idx)})")
    ax.fill_between(lag_axis * 1000, mn - se, mn + se, color=CLR[reg], alpha=0.2)
ax.axvline(0, color="k", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.5)
ax.set(xlabel="lag = neural - kinematic (ms)     < 0: spikes BEFORE the movement",
       ylabel="mean excess MI ± SEM (bits)",
       title="h(lag): the map averaged along each DIAGONAL\n"
             "(grey: only some kinematic windows reach this lag)")
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

# E: mean f(neu) per region
ax = _lf.add_subplot(gs[1, 1])
for reg in ("CP", "CB"):
    idx = [pi for pi in keep if regions_run[pi] == reg and _sub_ok(pi)]
    if not idx: continue
    Pm = np.stack([PROF_NEU[pi] for pi in idx])
    mn = np.nanmean(Pm, axis=0)
    se = np.nanstd(Pm, axis=0) / np.sqrt(np.maximum(np.isfinite(Pm).sum(0), 1))
    ax.plot(win_centers * 1000, mn, color=CLR[reg], lw=2, label=f"{reg} (n={len(idx)})")
    ax.fill_between(win_centers * 1000, mn - se, mn + se, color=CLR[reg], alpha=0.2)
ax.axvline(0, color="k", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.5)
ax.set(xlabel="NEURAL window centre re reach onset (ms)",
       ylabel="mean excess MI ± SEM (bits)",
       title="f(neural): the map averaged down each COLUMN\n"
             "(whatever kinematic window is read)")
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

# F: how to read it, caveats included
ax = _lf.add_subplot(gs[1, 2]); ax.axis("off")
ax.text(0, 1, "\n".join([
    "Two 1-D summaries of each neuron's 2-D map:",
    "",
    "  f(neural)  MI depends on WHEN the spikes are,",
    "             whatever kinematics is read",
    "             -> vertical stripes",
    "  h(lag)     MI depends on neural - kinematic time",
    "             -> diagonal stripes",
    "",
    "Each is built from half of the kinematic windows and",
    "scored (R2) on the OTHER half, row by row, so how",
    "STRONG the map is at a kinematic time never enters,",
    "and h(lag)'s extra flexibility buys nothing on noise.",
    "",
    "  dR2 = R2_lag - R2_neural",
    "",
    "WITH A GRAIN OF SALT: a ridge whose slope is not 1",
    "fits neither summary, and so does a noisy map, for",
    "completely different reasons. This is a population",
    "tendency, not a per-neuron verdict.",
    "",
    f"grid {n_kin} kin x {n_wins} neural windows, lag bin "
    f"{_lag_bin*1000:.0f} ms, CV gap {_gap}",
]), va="top", ha="left", fontsize=9, family="monospace", transform=ax.transAxes)

_lf.suptitle(f"Diagonal vs neural time  —  {lab}", fontsize=12)
if SAVE_FIGS:
    _savefig_retry(_lf, fdir / f"{RUN_NAME}_diagonal_vs_neural.png", dpi=150,
                   bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(_lf)


print(f"\nfigures -> {fdir}")

# %% ── 7. per-mouse mean landscapes ──────────────────────────────────────────
# 7. cell --> the same three maps as cell 4 (CP, CB, CB - CP), one row per mouse,
# plus an ALL MICE row at the bottom for reference. Saves one grid figure and,
# if PER_MOUSE_FILES, one figure per mouse as well.

MOUSE_OF = {pi: uids[run_idx[pi]].split("|")[1] for pi in keep}
MICE_RUN = sorted({MOUSE_OF[pi] for pi in keep})
ROWS = [(m, [pi for pi in keep if MOUSE_OF[pi] == m]) for m in MICE_RUN]
ROWS.append(("ALL MICE", list(keep)))

def _mean_maps(idx):
    """-> {'CP': map or None, 'CB': ..., 'n': {'CP': n, 'CB': n}}"""
    out, ns = {}, {}
    for reg in ("CP", "CB"):
        sel = [pi for pi in idx if regions_run[pi] == reg]
        ns[reg] = len(sel)
        out[reg] = (np.nanmean(EXC2[sel], axis=0)
                    if len(sel) >= PER_MOUSE_MIN_N else None)
    out["n"] = ns
    return out

MAPS = {name: _mean_maps(idx) for name, idx in ROWS}
print(f"\nper-mouse mean landscapes ({len(MICE_RUN)} mice)")
for name, _ in ROWS:
    print(f"  {name:12s} CP n={MAPS[name]['n']['CP']:4d}   CB n={MAPS[name]['n']['CB']:4d}")

# One colour scale for CP/CB across ALL panels, and one for every difference
# panel, so the rows can be compared with each other rather than each being
# stretched to its own range.
_v = [m[r].ravel() for m in MAPS.values()
      for r in ("CP", "CB") if m[r] is not None]
if not _v:
    raise ValueError(f"no mouse has >= PER_MOUSE_MIN_N ({PER_MOUSE_MIN_N}) neurons "
                     f"in either region -- lower it, or run more neurons")
_v = np.concatenate(_v)
VMAX = np.nanpercentile(_v, 99); VMIN = min(0.0, np.nanpercentile(_v, 1))
_d = [m["CB"] - m["CP"] for m in MAPS.values()
      if m["CP"] is not None and m["CB"] is not None]
DLIM = np.nanpercentile(np.abs(np.concatenate([d.ravel() for d in _d])), 99) if _d else 1.0


def _panel(ax, M, title, diff=False, vmin=None, vmax=None, dlim=None):
    """vmin/vmax/dlim default to the figure-wide scales (VMIN/VMAX/DLIM) that the
    GRID uses, where every mouse must share one scale to be comparable. The
    single-mouse files pass their own, so each is scaled to itself."""
    if M is None:
        ax.text(0.5, 0.5, "too few neurons", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="0.5")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
        return None
    if diff:
        _d = DLIM if dlim is None else dlim
        im = ax.imshow(M, origin="lower", aspect=ASPECT, extent=extent, cmap="RdBu_r",
                       norm=TwoSlopeNorm(vcenter=0, vmin=-_d, vmax=_d))
        c = "k"
    else:
        im = ax.imshow(M, origin="lower", aspect=ASPECT, extent=extent, cmap=CMAP,
                       vmin=VMIN if vmin is None else vmin,
                       vmax=VMAX if vmax is None else vmax)
        c = "w"
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(d, d, color=c, lw=1.0, ls="--", alpha=0.7)
    ax.axvline(0, color=c, lw=0.7, ls=":", alpha=0.7)
    ax.axhline(0, color=c, lw=0.7, ls=":", alpha=0.7)
    ax.set_title(title, fontsize=9)
    return im


_t_pm = time.time()
_H = 3.6 * len(ROWS)
fig, axes = plt.subplots(len(ROWS), 3, figsize=(15, _H), squeeze=False)
for r, (name, _) in enumerate(ROWS):
    M = MAPS[name]
    D = (M["CB"] - M["CP"]) if (M["CP"] is not None and M["CB"] is not None) else None
    for c, (key, Mi, ttl, isd) in enumerate((
            ("CP", M["CP"], f"CP  (n={M['n']['CP']})", False),
            ("CB", M["CB"], f"CB  (n={M['n']['CB']})", False),
            ("D",  D,       "CB - CP", True))):
        ax = axes[r, c]
        im = _panel(ax, Mi, ttl, diff=isd)
        if im is not None:
            plt.colorbar(im, ax=ax, label=("Δ mean MI (bits)" if isd else "mean MI (bits)"),
                         fraction=0.046, pad=0.02)
        if c == 0:
            ax.set_ylabel(f"{name}\nKINEMATIC centre (s)", fontsize=9)
        if r == len(ROWS) - 1:
            ax.set_xlabel("NEURAL window centre re reach onset (s)")
fig.suptitle(f"Mean 2D landscape per mouse  —  {lab}\n"
             f"dashed = simultaneity; ONE colour scale for every mouse, so the "
             f"rows compare (the single-mouse files are scaled to themselves)",
             fontsize=12)
# leave room for the two-line suptitle: the figure grows with the number of
# mice, so the margin has to be a fraction of the actual height, not a constant
plt.tight_layout(rect=(0, 0, 1, 1 - 0.75 / _H))
if SAVE_FIGS:
    _savefig_retry(fig, fdir / f"{RUN_NAME}_mean2D_per_mouse.png", dpi=150,
                   bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

if PER_MOUSE_FILES:
    for name, _ in ROWS:
        M = MAPS[name]
        D = (M["CB"] - M["CP"]) if (M["CP"] is not None and M["CB"] is not None) else None
        # scale each single-mouse figure to ITS OWN maps: CP and CB share one
        # scale so the regions stay comparable within the animal, the difference
        # gets its own. (The grid above keeps one scale across all mice.)
        _mv = [m.ravel() for m in (M["CP"], M["CB"]) if m is not None]
        if _mv:
            _mv = np.concatenate(_mv)
            _vmax = np.nanpercentile(_mv, 99); _vmin = min(0.0, np.nanpercentile(_mv, 1))
        else:
            _vmin, _vmax = VMIN, VMAX
        _dlim = np.nanpercentile(np.abs(D), 99) if D is not None else DLIM
        f1, ax1 = plt.subplots(1, 3, figsize=(15, 4.2))
        for c, (Mi, ttl, isd) in enumerate((
                (M["CP"], f"CP  (n={M['n']['CP']})", False),
                (M["CB"], f"CB  (n={M['n']['CB']})", False),
                (D,       "CB - CP", True))):
            im = _panel(ax1[c], Mi, ttl, diff=isd,
                        vmin=_vmin, vmax=_vmax, dlim=_dlim)
            if im is not None:
                plt.colorbar(im, ax=ax1[c],
                             label=("Δ mean MI (bits)" if isd else "mean MI (bits)"),
                             fraction=0.046, pad=0.02)
            ax1[c].set_xlabel("NEURAL window centre re reach onset (s)")
        ax1[0].set_ylabel("KINEMATIC window centre re reach onset (s)")
        f1.suptitle(f"{name}  —  {lab}\ncolour scale fitted to this "
                    f"mouse (the per-mouse GRID shares one scale across mice)",
                    fontsize=11)
        plt.tight_layout()
        if SAVE_FIGS:
            _savefig_retry(f1, fdir / f"{RUN_NAME}_mean2D_{name.replace(' ', '_')}.png",
                           dpi=150, bbox_inches="tight")
        if SHOW_FIGS: plt.show()
        plt.close(f1)

TIMES["per-mouse figures"] = time.time() - _t_pm
print(f"per-mouse landscapes -> {fdir} ({TIMES['per-mouse figures']:.0f}s)")

# ── where the time went ─────────────────────────────────────────────────────
def _fmt(v):
    return f"{v:.1f}s" if v < 90 else f"{v/60:.1f} min"

_tot = sum(TIMES.values())
print(f"\n{'step':32s} {'time':>10s} {'% of total':>11s}")
print("-" * 56)
for _k, _v in sorted(TIMES.items(), key=lambda kv: -kv[1]):
    print(f"{_k:32s} {_fmt(_v):>10s} {100*_v/max(_tot, 1e-9):10.1f}%")
print("-" * 56)
print(f"{'TOTAL (timed steps)':32s} {_fmt(_tot):>10s}")
if "2D landscapes (MI + nulls)" in TIMES:
    _per = TIMES["2D landscapes (MI + nulls)"] / max(len(keep), 1)
    print(f"\n{_per:.2f}s per neuron over {len(keep)} neurons x {n_kin*n_wins} cells "
          f"x {1+N_SHUF_2D} (obs+nulls)")
    print(f"scaling to the full population (~{len(NEURONS)} neurons) would take "
          f"~{_per*len(NEURONS)/60:.0f} min")

if SAVE_VARS:
    p = config.var_path(f"{RUN_NAME}_{'-'.join(FEATURES_2D)}_{NW_TAG}_cat{CAT_TAG}"
                        f"{'_' + SEL_TAG if SEL_TAG else ''}_2D.pkl")
    with open(p, "wb") as f:
        pickle.dump(dict(mi=MI2, null=NU2 if DO_NULL else None,
                         null_std=NSD2 if DO_NULL else None,
                         # diagonal vs neural time, per neuron (cross-validated)
                         R2_neural=R2N, R2_lag=R2L, dR2=DR2, structured=STRUCT,
                         prof_neural=PROF_NEU, prof_lag=PROF_LAG,
                         lag_axis=lag_axis, lag_full=lag_full,
                         LAG_BIN=_lag_bin, CV_GAP=_gap,
                         REACH_CAT=sorted(CAT_KEEP),
                         features=list(FEATURES_2D),
                         kin_centers=kin_centers, win_centers=win_centers,
                         KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP,
                         KIN_RANGE=KIN_RANGE,
                         WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP,
                         WIN_RANGE=WIN_RANGE, K=K_KSG,
                         uids=[uids[i] for i in run_idx], regions=regions_run,
                         N_TEST=N_TEST, N_SHUF_2D=N_SHUF_2D,
                         NEURAL_VAR=NEURAL_VAR), f, protocol=4)
    print(f"saved -> {p}")

# %%

# %%
