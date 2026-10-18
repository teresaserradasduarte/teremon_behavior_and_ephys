# %% ── 0. imports & parameters ────────────────────────────────────────────────
# Kinematics around the LEVER PUSH / PULL, extracted and CHECKED.
#
# Everything so far has been aligned to reaches, for which MATLAB already cut a
# 201-frame record per reach. There is no such cut for the push/pull, only the
# TIME of each one (bhv.pp_times). This script builds the equivalent: a window
# of BOTH paws' kinematics around every push/pull, plus the plots needed to
# believe the alignment. It writes a pickle a later MI script can consume the
# way 09/12 consume REACH_PX.
#
# WHAT t=0 IS (measured, cell 1 re-checks it every run): pp_times marks the END
# of the push/pull -- the moment the lever reaches criterion. It is followed by
# water delivery after a FIXED +222 ms (5-95%: +220..+224 ms, i.e. hardware, not
# behaviour). The paw movement that earns it runs from about -0.6 s to 0: speed
# peaks just before 0 and the dominant paw's x and z reach their minimum AT 0.
# So the push/pull itself lives BEFORE t=0, and the window after 0 holds the
# return plus the subsequent reach (median +1.05 s). Plot and model accordingly.
#
# THE EVENT CODES (bhv.valPushPull_invalPushPull_idx), as defined upstream:
#       0 = valid push     1 = valid pull     2 = invalid push     3 = invalid pull
# A VALID push/pull is the one that starts a new trial. The invalid ones still
# have an alignment problem upstream, so this script keeps only 0 and 1 by
# default. (bhv.init_invPush_invPull_idx merges 0 and 1 into a single "init"
# code, which is why it cannot separate push from pull -- it is not used here.)
#
# HOW THE PIECES FIT (verified in cell 2, not assumed):
#   bhv.pp_times                   absolute harp timestamps of each push/pull
#   trials.sync.frames_timestamps  (frames, trials) the same clock, per frame
#   trials.size_trial              real frame count per trial (arrays are padded
#                                  to the widest trial; see _trials_of below --
#                                  which trials make up the session trace is
#                                  MEASURED, not assumed)
#   session.<paw>_sess             the continuous trace = the trials concatenated
# So: timestamp -> (trial, frame) -> global sample -> slice the continuous trace.
#
# WHY THE CONTINUOUS TRACE: a valid push/pull STARTS a trial, so it sits exactly
# on a trial boundary -- measured, the median event is the last frame of the
# trial containing it. Cutting inside one trial gave windows with nothing after
# t=0 (77% coverage before, 23% after). The recording is continuous, so cutting
# across the boundary is correct and gives 100% coverage.
#
# THE LOAD CELL is extracted too (behavior.inputs.paws_all in the RAW tree:
# absolute harp timestamp + force in mV/V at ~1 kHz, same clock as pp_times, so
# no conversion is involved). The force IS the signal that defines a push/pull --
# the paw can reach the same criterion by different paths, which is why single
# paw traces look ragged while the force traces are stereotyped. Cell 5 plots it
# and it is saved next to the kinematics, resampled onto the camera grid.
#
# WHAT IS NOT HERE: spikes. st_init (spike times re initiation) exists in the
# ephys .mat but mi_io does not parse it into the cache, so aligning spikes to
# the push/pull is the next step, not this one.
import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.io as sio
from scipy.signal import savgol_filter

# The VS Code interactive window does not always start in this folder, and then
# "import mi_io" fails even though mi_io.py sits right next to this file. Find
# the folder that holds the project modules and work from there -- data paths
# all come from config, so moving the working directory touches nothing else.
_HERE = (os.path.dirname(os.path.abspath(__file__))
         if "__file__" in globals() else os.getcwd())
if not os.path.isfile(os.path.join(_HERE, "mi_io.py")):
    # the cwd may be the repo root, in which case the modules are in a SUBfolder
    _sub = os.path.join("ephys", "03_population_newMI_py")
    _p, _found = os.getcwd(), None
    while _found is None:
        for _c in (_p, os.path.join(_p, _sub)):
            if os.path.isfile(os.path.join(_c, "mi_io.py")): _found = _c; break
        if _found is None:
            if _p == os.path.dirname(_p):
                raise ModuleNotFoundError(
                    f"mi_io.py not found from the current directory ({os.getcwd()}). "
                    "Open this script's folder in VS Code, or set the interactive "
                    "window's working directory to it.")
            _p = os.path.dirname(_p)
    _HERE = _found
if _HERE not in sys.path: sys.path.insert(0, _HERE)
os.chdir(_HERE)
import mi_io, config

# ── the window ───────────────────────────────────────────────────────────────
PP_WIN    = (-3.0, 3.0)     # (s) EXTRACTED and saved, so other windows stay available
PLOT_WIN  = (-1.0, 0.5)     # (s) what the figures show; must sit inside PP_WIN.
                            # The movement is BEFORE 0 (see the header), so a
                            # window starting at -0.1 hides it.
FPS       = 400.0
UNITS     = "mm"            # "mm" (applies config.PX2MM) or "px"
PAWS      = ("dom", "nondom")   # which paws to extract; both is the honest check,
                                # since the pushing paw need not be the reaching one

# ── which events ─────────────────────────────────────────────────────────────
PP_CODE  = {0: "valid push", 1: "valid pull", 2: "invalid push", 3: "invalid pull"}
PP_KEEP  = (0, 1)           # valid only: the invalid ones are still mis-aligned
SPLIT    = (0, 1)           # codes drawn as separate conditions in every figure

# ── trial boundaries ─────────────────────────────────────────────────────────
STAY_IN_TRIAL = False       # True = clamp the window to the event's own trial
                            # (leaves the post-event half empty: see the header)

# ── load cell ────────────────────────────────────────────────────────────────
LOADCELL   = True           # extract + plot the force (needs the RAW behaviour tree)
# The alignment check on the force needs NO clock conversion: the force must
# DEFLECT at t=0, upward for a push and downward for a pull. (An earlier version
# checked pp_times against the raw init_time instead, which meant solving a
# session-seconds offset; that matched in only 22/29 sessions for reasons that
# have nothing to do with the extraction -- every session's events do fall
# inside the raw load-cell span -- so it was dropped in favour of this.)
LC_BASE = (-1.00, -0.60)    # (s) baseline window, before any movement
LC_PEAK = (-0.05,  0.20)    # (s) where the deflection is measured
LC_MIN_N = 10               # a condition with fewer events than this is not
                            # judged: the sign of a 6-event median means little

# ── scope ────────────────────────────────────────────────────────────────────
SESSIONS   = None           # None = every session in the cache, or [(mouse, sess), ...]
EXAMPLE    = None           # (mouse, sess) for the SELF-TEST session; None = first
N_EXAMPLES = 6              # single events drawn per condition, taken from AS MANY
                            # DIFFERENT MICE/SESSIONS as possible (round-robin over
                            # mice, then sessions) so they are representative
VALIDATE_N = 40             # cell 2 re-cuts this many REACHES with the same code
                            # and PLOTS them on top of MATLAB's own cut, so the
                            # alignment can be seen. Nothing from it enters PP or
                            # the pickle -- reaches are only the known answer.

# ── output ───────────────────────────────────────────────────────────────────
OUT_SUBDIR = "pushpull_kinematics"
RUN_NAME   = "pp"
PP_PKL     = "pushpull_kin.pkl"
CACHE_MODE = "load"         # "load" = reuse the pickle when it matches; "rerun"
SAVE_VARS  = True
SAVE_FIGS  = True
SHOW_FIGS  = True
# the project's push/pull colours (as in the MATLAB figures; /256 as written there)
CLR = {0: np.array([132, 160, 124]) / 256,     # valid push   (sage green)
       1: np.array([26, 94, 99]) / 256,        # valid pull   (dark teal)
       2: np.array([132, 160, 124]) / 256 * 0.55 + 0.45,   # invalid push, lighter
       3: np.array([26, 94, 99]) / 256 * 0.55 + 0.45}      # invalid pull, lighter
PAW_STYLE = {"dom": dict(ls="-", lw=2.0), "nondom": dict(ls="--", lw=1.4)}

N_PRE   = int(round(abs(PP_WIN[0]) * FPS))
N_POST  = int(round(PP_WIN[1] * FPS))
T_WIN   = np.arange(-N_PRE, N_POST + 1) / FPS
PLOT_M  = (T_WIN >= PLOT_WIN[0]) & (T_WIN <= PLOT_WIN[1])
PX2MM   = np.asarray(config.PX2MM, float)
SCALE   = PX2MM[:, None] if UNITS == "mm" else np.ones((3, 1))
PAW_VAR = {"dom": "ppDom_paw", "nondom": "ppNonDom_paw"}
print(f"extract {PP_WIN} s ({T_WIN.size} frames at {FPS:.0f} fps), plot {PLOT_WIN} s "
      f"| paws {PAWS} in {UNITS}")
print(f"keeping codes {PP_KEEP} = {[PP_CODE[c] for c in PP_KEEP]}")


# %% ── 1. extract ────────────────────────────────────────────────────────────
# 1. cell --> per session: the push/pull times, their codes, and a window of each
# paw around every event. Reads only the slices it needs.
if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())
KEYS = sorted({(n["mouse"], n["sess"]) for n in NEURONS}) if SESSIONS is None \
       else [tuple(k) for k in SESSIONS]
print(f"{len(KEYS)} session(s)")


def _trials_of(siz, n_samples):
    """Which trials make up the continuous session trace, decided by MATCHING
    the sample count rather than by a heuristic.
    Measured over the 29 sessions: most concatenate EVERY trial, but 9 carry a
    dummy trial 0 that is left out. The dummy is not identifiable by its size
    (in 4_Milka R6/R7 it is not the widest trial) nor by its timestamps (in
    1_CoteDor R1 the LAST trial has no timestamps at all, yet its samples ARE in
    the trace). Trying both candidate sets and keeping the one whose sizes add
    up is the only rule that held on all of them.
    -> (label, list of trial numbers)"""
    cands = {"all trials":      [t for t in range(siz.size) if siz[t] > 0],
             "trial 0 dropped": [t for t in range(1, siz.size) if siz[t] > 0]}
    hit = [(k, v) for k, v in cands.items() if sum(siz[t] for t in v) == n_samples]
    if len(hit) != 1:
        raise ValueError(
            f"cannot tell which trials build the session trace: "
            f"{ {k: sum(siz[t] for t in v) for k, v in cands.items()} } vs "
            f"{n_samples} samples")
    return hit[0]


def _trial_clock(f, real):
    """-> per-trial timestamp arrays. A trial can be entirely NaN here (one
    session's last trial has no sync at all); such a trial simply never matches
    an event, while its samples still occupy their place in the trace."""
    siz = f["trials"]["size_trial"][:].ravel().astype(int)
    return [np.asarray(f["trials"]["sync"]["frames_timestamps"][:siz[t], t], float)
            for t in real]


def _spans(TS):
    """first/last finite timestamp per trial; NaN for a trial without any."""
    lo = np.array([np.nanmin(t) if np.isfinite(t).any() else np.nan for t in TS])
    hi = np.array([np.nanmax(t) if np.isfinite(t).any() else np.nan for t in TS])
    return lo, hi


LOCATE_TOL = 0.5 / FPS      # (s) half a frame


def _locate(s0, s1, TS, t_abs):
    """absolute timestamp -> (position in the real-trial list, frame index),
    or (None, None) when it falls between trials.

    The rig starts a new camera trial AT the trial initiation, so the timestamp
    of a VALID push/pull is the last frame of one camera trial and the first
    frame of the next. Compared exactly, 14% of valid events (531 over the 29
    sessions) fell in the sub-millisecond crack between those two frames, were
    placed in neither, and got no paws and no force at all -- which is what 21
    reported as "untracked". Half a frame of tolerance, used only as a
    fallback, places them; events the exact test places are unchanged."""
    with np.errstate(invalid="ignore"):
        c = np.flatnonzero((s0 <= t_abs) & (t_abs <= s1))
        if not c.size:
            # only when the exact test finds nothing: everything it placed
            # before stays exactly as it was (the two frames at a boundary carry
            # the same timestamp, so re-choosing between them would move those
            # events by one frame)
            c = np.flatnonzero((s0 - LOCATE_TOL <= t_abs) & (t_abs <= s1 + LOCATE_TOL))
    if not c.size: return None, None
    j = int(c[0])
    return j, int(np.nanargmin(np.abs(TS[j] - t_abs)))


def _cut(TRACE, g, lo_lim=None, hi_lim=None, n_pre=N_PRE, n_post=N_POST):
    """window around GLOBAL sample g of a continuous (3, N) trace.
    -> (3, n_pre+n_post+1) in UNITS, plus the fraction of finite samples."""
    out = np.full((3, n_pre + n_post + 1), np.nan)
    lo, hi = g - n_pre, g + n_post
    a = max(lo, 0 if lo_lim is None else lo_lim)
    b = min(hi, (TRACE.shape[1] - 1) if hi_lim is None else hi_lim)
    if b >= a:
        out[:, a - lo:b - lo + 1] = TRACE[:, a:b + 1]
    return out * SCALE, float(np.mean(np.isfinite(out[0])))


_P = config.var_path(PP_PKL)
PP = None
if CACHE_MODE == "load" and _P.is_file():
    with open(_P, "rb") as fh: _c = pickle.load(fh)
    if (_c.get("win") == list(PP_WIN) and _c.get("units") == UNITS
            and tuple(_c.get("paws", ())) == tuple(PAWS)
            and _c.get("stay_in_trial") == STAY_IN_TRIAL
            and _c.get("loadcell") == LOADCELL
            and _c.get("v") == "v2"        # v2: events at a trial boundary placed
            and set(_c["sessions"]) >= set(KEYS)):
        PP = _c["sessions"]; print(f"loaded {_P.name} -- no re-extraction")
    else:
        print("pickle exists but was built with other settings -> re-extracting")

if PP is None:
    PP = {}
    for key in KEYS:
        mouse, sess = key
        b = sio.loadmat(config.bhv_fundamentals_mat(mouse, sess), squeeze_me=True,
                        struct_as_record=False)["bhv"]
        pp_t = np.asarray(b.pp_times, float).ravel()
        code = np.asarray(b.valPushPull_invalPushPull_idx, float).ravel()
        rt   = np.asarray(b.reach_times_inVec, float).ravel()
        with h5py.File(config.behav_session_mat(mouse, sess), "r") as f:
            siz  = f["trials"]["size_trial"][:].ravel().astype(int)
            TR   = {p: np.asarray(f["session"][PAW_VAR[p] + "_sess"][:], float)
                    for p in PAWS}
            _n   = TR[PAWS[0]].shape[1]
            rule, real = _trials_of(siz, _n)
            TS   = _trial_clock(f, real)
            s0, s1 = _spans(TS)
            OFF  = np.concatenate([[0], np.cumsum([siz[t] for t in real])])
            for p, A in TR.items():
                if A.shape[1] != _n:
                    raise ValueError(f"{key}: {PAW_VAR[p]}_sess has {A.shape[1]} "
                                     f"samples but {PAW_VAR[PAWS[0]]}_sess has {_n}")
            K    = {p: np.full((pp_t.size, 3, T_WIN.size), np.nan) for p in PAWS}
            frac = np.zeros(pp_t.size)
            trl  = np.full(pp_t.size, -1)
            edge = np.full((pp_t.size, 2), np.nan)      # (s) to trial start / end
            for e, t_abs in enumerate(pp_t):
                j, i = _locate(s0, s1, TS, t_abs)
                if j is None: continue
                trl[e] = real[j]
                g = OFF[j] + i
                edge[e] = [-i / FPS, (siz[real[j]] - 1 - i) / FPS]
                for p in PAWS:
                    K[p][e], fr = _cut(TR[p], g,
                                       OFF[j] if STAY_IN_TRIAL else None,
                                       OFF[j + 1] - 1 if STAY_IN_TRIAL else None)
                    if p == PAWS[0]: frac[e] = fr
        # ---- load cell: the signal that defines the event ------------------
        FRC = np.full((pp_t.size, T_WIN.size), np.nan)
        lc_defl = {}                 # per code: median force deflection at t=0
        if LOADCELL:
            _raw = sio.loadmat(config.behav_raw_session_mat(mouse, sess),
                               squeeze_me=True, struct_as_record=False)["behavior"]
            _pa  = np.asarray(_raw.inputs.paws_all, float)
            _o   = np.argsort(_pa[:, 0])
            _tlc, _f = _pa[_o, 0], _pa[_o, 1]        # absolute harp, mV/V
            # same clock as pp_times, so the window is just an interpolation.
            # Resampled onto the camera grid T_WIN so force and kinematics share
            # one time base (native rate is ~1 kHz).
            for e, t_abs in enumerate(pp_t):
                if trl[e] < 0: continue
                tt = t_abs + T_WIN
                v  = np.interp(tt, _tlc, _f, left=np.nan, right=np.nan)
                v[(tt < _tlc[0]) | (tt > _tlc[-1])] = np.nan
                FRC[e] = v
            # CHECK (no clock conversion): the force must deflect at t=0, up for
            # a push and down for a pull. This uses only the cut window, so it
            # tests the alignment that the rest of the script depends on.
            _bs = (T_WIN >= LC_BASE[0]) & (T_WIN <= LC_BASE[1])
            _pk = (T_WIN >= LC_PEAK[0]) & (T_WIN <= LC_PEAK[1])
            for _c in PP_KEEP:
                _m = (code == _c) & (trl >= 0)
                if not _m.any(): continue
                _A = FRC[_m]
                _b0 = np.nanmean(_A[:, _bs], axis=1)
                _dd = _A[:, _pk] - _b0[:, None]
                # the signed extreme: whichever of min/max is further from baseline
                _ex = np.where(np.nanmax(_dd, axis=1) >= -np.nanmin(_dd, axis=1),
                               np.nanmax(_dd, axis=1), np.nanmin(_dd, axis=1))
                lc_defl[_c] = (float(np.nanmedian(_ex)), int(_m.sum()))

        # Two independent timing checks, both re-measured every run:
        #   next_reach  the reach must FOLLOW the push/pull (median ~1.05 s)
        #   next_water  water delivery follows at a FIXED latency, which is what
        #               identifies pp_times as the lever-criterion moment
        nxt = np.array([(rt[rt > t][0] - t) if np.any(rt > t) else np.nan
                        for t in pp_t])
        with h5py.File(config.behav_session_mat(mouse, sess), "r") as f:
            _dl = np.asarray(f["trials"]["sync"]["timestamp_delivered"][:], float).ravel()
        _dl = _dl[np.isfinite(_dl)]
        wat = np.array([(_dl[_dl > t][0] - t) if np.any(_dl > t) else np.nan
                        for t in pp_t])
        PP[key] = dict(t=pp_t, code=code, trial=trl, frac=frac, next_reach=nxt,
                       next_water=wat, force=FRC.astype(np.float32),
                       lc_defl=lc_defl,
                       trial_edge=edge,      # next_reach is a CHECK field (the
                       # only reach-derived quantity saved); the kinematics are
                       # push/pull only
                       **{f"kin_{p}": K[p].astype(np.float32) for p in PAWS})
        _in = trl >= 0
        print(f"  {mouse}|{sess}: [{rule}, {len(real)} trials] {pp_t.size} events, "
              f"{int(_in.sum())} located, coverage {100*frac[_in].mean():.0f}%  |  "
              + ", ".join(f"{PP_CODE[c]}: {int(np.sum((code == c) & _in))}"
                          for c in sorted(PP_CODE)))
        _w = wat[_in & np.isin(code, list(PP_KEEP))]
        _w = _w[np.isfinite(_w)]
        if _w.size:
            print(f"      pp -> water: median {np.nanmedian(_w)*1000:+.0f} ms "
                  f"(5-95% {np.nanpercentile(_w,5)*1000:+.0f}.."
                  f"{np.nanpercentile(_w,95)*1000:+.0f}) -- a tight range means "
                  f"pp_times is the lever event, not a behavioural estimate")
        if LOADCELL and lc_defl:
            # judged only where there are enough events; an absent condition
            # is "not applicable", not a failure
            _bads = [c for c, (v, n) in lc_defl.items()
                     if n >= LC_MIN_N and ((c == 0 and v <= 0) or (c == 1 and v >= 0))]
            print("      load cell deflection at t=0: "
                  + ", ".join(f"{PP_CODE[c]} {v:+.0f} (n={n})"
                              for c, (v, n) in sorted(lc_defl.items()))
                  + " mV/V"
                  + ("" if not _bads
                     else "   <-- UNEXPECTED SIGN for " + ", ".join(PP_CODE[c] for c in _bads)))


def _keep(S, codes=PP_KEEP):
    return (S["trial"] >= 0) & np.isin(S["code"], list(codes))

print()
for c in sorted(PP_CODE):
    n = sum(int(np.sum((S["code"] == c) & (S["trial"] >= 0))) for S in PP.values())
    print(f"  code {c} = {PP_CODE[c]:14s} {n:5d} events"
          f"{'   <- kept' if c in PP_KEEP else ''}")


# %% ── 2. SELF-TEST of the cutting code, with the proof PLOTTED ────────────
# 2. cell --> NOT an analysis: nothing here enters PP or the pickle.
#
# There is no ground truth for a push/pull window. If the chain
# timestamp -> trial -> frame -> slice were off by a few frames, every figure
# below would still look plausible and nothing would give it away. But MATLAB
# already cut the REACHES off the same clock, so those are a known answer:
# re-cutting them with THIS script's code must reproduce REACH_PX exactly
# (150 frames before the max-reach frame, 50 after). The figure overlays the two
# so the match is visible, not merely asserted, and the same verified code then
# cuts the push/pull windows.
key0 = tuple(EXAMPLE) if EXAMPLE else KEYS[0]
b0   = sio.loadmat(config.bhv_fundamentals_mat(*key0), squeeze_me=True,
                   struct_as_record=False)["bhv"]
rt0  = np.asarray(b0.reach_times_inVec, float).ravel()
with h5py.File(config.behav_session_mat(*key0), "r") as f:
    siz = f["trials"]["size_trial"][:].ravel().astype(int)
    DOM = np.asarray(f["session"]["ppDom_paw_sess"][:], float)
    rule, real = _trials_of(siz, DOM.shape[1])
    TS  = _trial_clock(f, real)
    s0, s1 = _spans(TS)
    OFF = np.concatenate([[0], np.cumsum([siz[t] for t in real])])
    cuts, refs, which, skip = [], [], [], 0
    for k in range(0, rt0.size, max(1, rt0.size // VALIDATE_N)):
        j, i = _locate(s0, s1, TS, rt0[k])
        if j is None: skip += 1; continue
        lo = OFF[j] + i - config.REACH_ENDPOINT_FRAME
        hi = lo + config.REACH_N_FRAMES
        if lo < 0 or hi > DOM.shape[1]: skip += 1; continue
        cuts.append(DOM[:, lo:hi]); refs.append(REACH_PX[key0][k]); which.append(k)
CUT = np.asarray(cuts); REF = np.asarray(refs)        # (n, 3, 201) in PIXELS
ERR = np.nanmax(np.abs(CUT - REF), axis=(1, 2))
print(f"\nSELF-TEST on {key0[0]}|{key0[1]} [{rule}]: re-cut {CUT.shape[0]} reaches "
      f"({skip} skipped)")
print(f"  max |difference| vs REACH_PX: {np.nanmax(ERR):.3g}")
print("  identical: clock, frame indexing and px units confirmed"
      if np.nanmax(ERR) < 1e-9 else
      "  !! NOT identical -- the alignment is wrong, stop here")

T_R = (np.arange(config.REACH_N_FRAMES) - config.REACH_ENDPOINT_FRAME) / FPS * 1000
fig, axes = plt.subplots(2, 3, figsize=(18, 8))
for n_, ax in enumerate(axes[0]):
    if n_ >= CUT.shape[0]: ax.set_visible(False); continue
    for c, nm in enumerate("xyz"):
        ax.plot(T_R, CUT[n_, c], lw=3.0, alpha=0.45,
                label=f"{nm}: re-cut here" if n_ == 0 else None)
        ax.plot(T_R, REF[n_, c], lw=1.0, ls="--", color="k",
                label="MATLAB REACH_PX" if (n_ == 0 and c == 0) else None)
    ax.axvline(0, color="crimson", lw=1.2, ls=":")
    ax.set(xlabel="time re reach endpoint (ms)", ylabel="paw (px)",
           title=f"reach {which[n_]}   max diff {ERR[n_]:.3g}")
    ax.spines[["top", "right"]].set_visible(False)
    if n_ == 0: ax.legend(fontsize=7)

ax = axes[1, 0]
ax.plot(T_R, (CUT - REF).reshape(-1, CUT.shape[-1]).T, lw=0.8, color="0.3")
ax.axhline(0, color="crimson", lw=1.0)
ax.set(xlabel="time re reach endpoint (ms)", ylabel="re-cut − REACH_PX (px)",
       title=f"difference, all {CUT.shape[0]} reaches x 3 coords\n"
             f"(a flat line at 0 is the pass)")
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1, 1]
ax.plot(which, ERR, "o", ms=5, color="teal")
ax.axhline(0, color="crimson", lw=1.0)
ax.set(xlabel="reach index in the session", ylabel="max |difference| (px)",
       title="worst sample per reach\n(every point on zero = exact)")
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1, 2]
_v = savgol_filter(CUT, 11, 3, deriv=1, delta=1/FPS, axis=2)
_sp = np.linalg.norm(_v * PX2MM[None, :, None], axis=1)
m = np.nanmean(_sp, axis=0)
sd = np.nanstd(_sp, axis=0) / np.sqrt(_sp.shape[0])
ax.plot(T_R, m, color="k", lw=2); ax.fill_between(T_R, m - sd, m + sd, color="k", alpha=0.2)
ax.axvline(0, color="crimson", lw=1.2, ls=":", label="reach endpoint (frame "
           f"{config.REACH_ENDPOINT_FRAME})")
ax.set(xlabel="time re reach endpoint (ms)", ylabel="speed ± SEM (mm/s)",
       title="mean speed of the re-cut reaches\n(the movement must end at 0)")
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

fig.suptitle(f"SELF-TEST of the cutting code  —  {key0[0]}|{key0[1]}  "
             f"(reaches used as the known answer; none of this is saved)",
             fontsize=12)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir(OUT_SUBDIR) / f"{RUN_NAME}_selftest_recut.png",
                dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)
del CUT, REF, cuts, refs            # the test fixture is not kept


# %% ── 3. figures: single events, split push vs pull ─────────────────────────
# 3. cell --> one panel per event, both paws' speed, push in one column and pull
# in the other. If the alignment is right the speed should lift at t=0 in every
# panel, and the pushing paw should be the one that moves.
fdir = config.fig_dir(OUT_SUBDIR)
TP   = T_WIN[PLOT_M]

def _speed(K):
    """(..., 3, T) -> (..., T) speed in UNITS/s, NaN-safe."""
    flat = np.asarray(K, float).reshape(-1, 3, K.shape[-1])
    out  = np.full((flat.shape[0], K.shape[-1]), np.nan)
    for i in range(flat.shape[0]):
        ok = np.isfinite(flat[i]).all(axis=0)
        if ok.sum() < 11: continue
        v = savgol_filter(flat[i][:, ok], 11, 3, deriv=1, delta=1/FPS, axis=1)
        out[i, ok] = np.linalg.norm(v, axis=0)
    return out.reshape(K.shape[:-2] + (K.shape[-1],))

def _pick_examples(c, n):
    """n events of condition c, spread over as many MICE and SESSIONS as there
    are: round-robin over mice, then over that mouse's sessions, taking a middle
    event each time. One session's idiosyncrasies cannot dominate the panel."""
    by_mouse = {}
    for k in KEYS:
        e = np.flatnonzero(_keep(PP[k], (c,)))
        if e.size: by_mouse.setdefault(k[0], []).append((k, list(e)))
    out, mice = [], sorted(by_mouse)
    if not mice: return out
    r = 0
    while len(out) < n and any(any(v for _, v in sess) for sess in by_mouse.values()):
        m = mice[r % len(mice)]
        for k, evs in by_mouse[m]:
            if evs:
                out.append((k, evs.pop(len(evs) // 2))); break
        r += 1
        if r > 500: break
    return out[:n]

EXAMPLES = {c: _pick_examples(c, N_EXAMPLES) for c in SPLIT}
for c in SPLIT:
    print(f"{PP_CODE[c]} examples: "
          + ", ".join(f"{k[0]}|{k[1]}#{e}" for k, e in EXAMPLES[c]))

fig, axes = plt.subplots(N_EXAMPLES, len(SPLIT), figsize=(7 * len(SPLIT), 2.0 * N_EXAMPLES),
                         sharex=True, squeeze=False)
for col, c in enumerate(SPLIT):
    sel = EXAMPLES[c]
    for row in range(N_EXAMPLES):
        ax = axes[row, col]
        if row >= len(sel):
            ax.set_visible(False); continue
        k, e = sel[row]; S = PP[k]
        for p in PAWS:
            v = _speed(S[f"kin_{p}"][e][None])[0]
            ax.plot(TP, v[PLOT_M], color=CLR[c], label=p, **PAW_STYLE[p])
        ax.axvline(0, color="k", lw=1.2, ls="--")
        nr = S["next_reach"][e]
        if np.isfinite(nr) and PLOT_WIN[0] <= nr <= PLOT_WIN[1]:
            ax.axvline(nr, color="0.3", lw=1.0, ls=":")
        ax.set_ylabel(f"{k[0].split('_')[0]}|{k[1]} #{e}\n({UNITS}/s)", fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        if row == 0:
            ax.set_title(f"{PP_CODE[c]}", fontsize=11, color=CLR[c])
            ax.legend(fontsize=7, ncol=len(PAWS))
    axes[-1, col].set_xlabel("time re push/pull (s)")
fig.suptitle("single events, paw speed  —  one per mouse/session in rotation  "
             "(dotted = next reach)", fontsize=12)
plt.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / f"{RUN_NAME}_single_events.png", dpi=150,
                          bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 4. figures: averages, split push vs pull ──────────────────────────────
# 4. cell --> ROWS = push / pull (one condition per axis), COLUMNS = x / y / z /
# speed, and BOTH PAWS drawn together in every panel (dominant solid, non-
# dominant dashed). Each column shares its y axis across the rows, so push and
# pull can be compared directly.
ALL = {c: {p: np.concatenate([np.asarray(PP[k][f"kin_{p}"][_keep(PP[k], (c,))], float)
                              for k in KEYS]) for p in PAWS} for c in SPLIT}
_wall = np.concatenate([PP[k]["next_water"][_keep(PP[k])] for k in KEYS])
_W_MED = float(np.nanmedian(_wall)) if np.isfinite(_wall).any() else np.nan
print(f"water delivery at {_W_MED*1000:+.0f} ms after the push/pull "
      f"(dotted line in the averages)")
for c in SPLIT:
    print(f"{PP_CODE[c]}: {ALL[c][PAWS[0]].shape[0]} events pooled")

fig, axes = plt.subplots(len(SPLIT), 4, figsize=(19, 4.0 * len(SPLIT)),
                         sharex=True, sharey="col", squeeze=False)
for r, c in enumerate(SPLIT):
    for cc, (q, nm) in enumerate(((0, "x"), (1, "y"), (2, "z"), (None, "speed"))):
        ax = axes[r, cc]
        for p in PAWS:
            A = ALL[c][p]
            if not A.size: continue
            X = _speed(A) if q is None else A[:, q]
            m = np.nanmean(X, axis=0)[PLOT_M]
            n = np.sum(np.isfinite(X), axis=0)[PLOT_M]
            se = (np.nanstd(X, axis=0)[PLOT_M]) / np.sqrt(np.maximum(n, 1))
            ax.plot(TP, m, color=CLR[c], label=f"{p} paw", **PAW_STYLE[p])
            ax.fill_between(TP, m - se, m + se, color=CLR[c], alpha=0.18)
        ax.axvline(0, color="k", lw=1.0, ls="--")
        if np.isfinite(_W_MED) and PLOT_WIN[0] <= _W_MED <= PLOT_WIN[1]:
            ax.axvline(_W_MED, color="0.4", lw=1.0, ls=":")
        ax.set_title(f"{PP_CODE[c]} (n={ALL[c][PAWS[0]].shape[0]})  —  {nm}",
                     fontsize=10, color=CLR[c])
        ax.set_ylabel(f"{nm} ({UNITS}/s)" if q is None else f"{nm} ({UNITS})")
        ax.spines[["top", "right"]].set_visible(False)
        if r == len(SPLIT) - 1: ax.set_xlabel("time re push/pull (s)")
        if r == 0 and cc == 0: ax.legend(fontsize=8)
plt.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / f"{RUN_NAME}_averages.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 5. figures: every event as a heatmap, plus the timing checks ──────────
fig = plt.figure(figsize=(18, 4.6 * len(PAWS)))
gs  = fig.add_gridspec(len(PAWS), len(SPLIT) + 2, hspace=0.35, wspace=0.3)
for r, p in enumerate(PAWS):
    V = {c: _speed(ALL[c][p]) for c in SPLIT}
    vm = np.nanpercentile(np.concatenate([v[np.isfinite(v)] for v in V.values()]), 99)
    for cc, c in enumerate(SPLIT):
        ax = fig.add_subplot(gs[r, cc])
        nr = np.concatenate([PP[k]["next_reach"][_keep(PP[k], (c,))] for k in KEYS])
        o  = np.argsort(np.where(np.isfinite(nr), nr, np.inf))
        ax.imshow(V[c][o][:, PLOT_M], aspect="auto", origin="lower", cmap="magma",
                  vmin=0, vmax=vm, extent=[TP[0], TP[-1], 0, V[c].shape[0]])
        ax.axvline(0, color="w", lw=1.2, ls="--")
        ax.set(xlabel="time re push/pull (s)", ylabel=f"{p} paw\nevent (sorted)",
               title=f"{PP_CODE[c]}  (n={V[c].shape[0]})")
    if r == 0:
        ax = fig.add_subplot(gs[:, len(SPLIT)])
        for c in SPLIT:
            nr = np.concatenate([PP[k]["next_reach"][_keep(PP[k], (c,))] for k in KEYS])
            nr = nr[np.isfinite(nr)]
            ax.hist(nr[nr < 6], bins=60, color=CLR[c], alpha=0.55, edgecolor="white",
                    label=f"{PP_CODE[c]}: median {np.median(nr):.2f} s")
        ax.set(xlabel="push/pull -> next reach (s)", ylabel="# events",
               title="latency to the next reach\n(the reach must FOLLOW)")
        ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

        ax = fig.add_subplot(gs[:, len(SPLIT) + 1])
        cov = np.mean(np.isfinite(np.concatenate(
            [ALL[c][PAWS[0]][:, 0] for c in SPLIT])), axis=0) * 100
        ax.plot(T_WIN, cov, color="teal", lw=2)
        ax.axvspan(PLOT_WIN[0], PLOT_WIN[1], color="0.85", zorder=0, lw=0)
        ax.axvline(0, color="k", lw=1.0, ls="--")
        ax.set(xlabel="time re push/pull (s)", ylabel="% of events with data",
               ylim=(0, 105), title=f"coverage of the EXTRACTED window\n"
               f"(grey = what the figures show; STAY_IN_TRIAL={STAY_IN_TRIAL})")
        ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"Every event, and the timing checks  —  {len(KEYS)} session(s)",
             fontsize=12)
if SAVE_FIGS: fig.savefig(fdir / f"{RUN_NAME}_events_and_checks.png", dpi=150,
                          bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)
print(f"figures -> {fdir}")


# %% ── 5b. LOAD CELL: the signal that defines the push/pull ────────────────
# 5b. cell --> the force around every event, cut at pp_times on the SAME absolute
# clock (no conversion). This is the check that matters: the paw can reach the
# lever criterion by different paths, so single paw traces look ragged, but the
# force must be stereotyped -- rising to a peak at 0 for a push, falling through
# 0 for a pull. If that is not what this shows, the alignment is wrong.
# Right-hand panel: the same means per MOUSE, so one animal cannot carry it.
if LOADCELL:
    FR = {c: np.concatenate([np.asarray(PP[k]["force"][_keep(PP[k], (c,))], float)
                             for k in KEYS]) for c in SPLIT}
    MOUSE_OF = {k: k[0] for k in KEYS}
    fig, axes = plt.subplots(1, len(SPLIT) + 1, figsize=(6.2 * (len(SPLIT) + 1), 4.4))
    for ax, c in zip(axes, SPLIT):
        A = FR[c][:, PLOT_M]
        if A.size:
            ax.plot(TP, A[:80].T, color=CLR[c], lw=0.4, alpha=0.25)
            ax.plot(TP, np.nanmean(A, axis=0), color=CLR[c], lw=2.8)
        ax.axvline(0, color="orange", lw=1.6, ls="--")
        ax.axhline(0, color="k", lw=0.6)
        if np.isfinite(_W_MED) and PLOT_WIN[0] <= _W_MED <= PLOT_WIN[1]:
            ax.axvline(_W_MED, color="0.4", lw=1.0, ls=":")
        ax.set(xlabel="time from trial initiation (s)", ylabel="force (mV/V)",
               title=f"{PP_CODE[c]}  (n={A.shape[0]})")
        ax.spines[["top", "right"]].set_visible(False)
    ax = axes[-1]
    for c in SPLIT:
        for m in sorted({MOUSE_OF[k] for k in KEYS}):
            A = np.concatenate([np.asarray(PP[k]["force"][_keep(PP[k], (c,))], float)
                                for k in KEYS if MOUSE_OF[k] == m])
            if A.shape[0] < 5: continue
            ax.plot(TP, np.nanmean(A[:, PLOT_M], axis=0), color=CLR[c], lw=1.4,
                    alpha=0.8, label=f"{PP_CODE[c]} {m.split('_')[0]}")
    ax.axvline(0, color="orange", lw=1.6, ls="--"); ax.axhline(0, color="k", lw=0.6)
    ax.set(xlabel="time from trial initiation (s)", ylabel="force (mV/V)",
           title="mean per mouse\n(every animal should show the same sign)")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=6, ncol=2)
    _d = {k: PP[k].get("lc_defl", {}) for k in KEYS}
    _bad = [k for k, v in _d.items()
            if any(n >= LC_MIN_N and ((c == 0 and x <= 0) or (c == 1 and x >= 0))
                   for c, (x, n) in v.items())]
    fig.suptitle("LOAD CELL around the push/pull  —  cut at pp_times, same clock, "
                 "no conversion\n"
                 f"force deflection at t=0 has the expected sign in "
                 f"{len(KEYS)-len(_bad)}/{len(KEYS)} sessions "
                 f"(push up, pull down)", fontsize=12)
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(fdir / f"{RUN_NAME}_loadcell.png", dpi=150,
                              bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)
    if _bad: print("  !! sessions whose force does not deflect as expected:", _bad)


# %% ── 6. save ───────────────────────────────────────────────────────────────
if SAVE_VARS:
    out = dict(sessions=PP, win=list(PP_WIN), t=T_WIN, fps=FPS, units=UNITS, v="v2",
               paws=list(PAWS), stay_in_trial=STAY_IN_TRIAL, codes=PP_CODE,
               loadcell=LOADCELL, force_units="mV/V")
    with open(_P, "wb") as fh: pickle.dump(out, fh, protocol=4)
    print(f"saved -> {_P}  ({_P.stat().st_size/1e6:.0f} MB)")
    print("  sessions[(mouse, sess)] = dict(t, code, trial, frac, next_reach, "
          "next_water, trial_edge, lc_defl, force (n,T) mV/V, "
          + ", ".join(f"kin_{p} (n,3,T) float32" for p in PAWS) + ")")

# %%
