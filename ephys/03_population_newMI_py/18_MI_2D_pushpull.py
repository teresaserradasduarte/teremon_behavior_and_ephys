# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# TIME x TIME multi-dimensional MI, aligned to the PUSH/PULL (trial initiation).
#
# 12 does this around the REACH: kinematics from reach_px, spikes from st_reach,
# t=0 = the reach endpoint. Here the same question is asked around the LEVER
# event: kinematics from 16's pushpull_kin.pkl, spikes from st_init, t=0 = the
# moment the push/pull criterion is crossed (pp_times).
#
#        MI( kin_features(kin_window) ; spike_count(neural_window) )
#
#   rows    = kinematic window centre   (when the MOVEMENT is read)
#   columns = neural window centre      (when the SPIKES are counted)
#
# Both axes share t=0 = the push/pull time, which is the zero of st_init and the
# zero of 16's extraction window, so no offset is involved.
#
# NOTE on where the movement is: 16 established that pp_times marks the END of
# the push/pull (water follows at a fixed +222 ms), so the paw movement happens
# BEFORE 0. KIN_RANGE is set accordingly -- a range starting at -0.1 s would
# miss the movement entirely.
#
# WHY A SEPARATE SCRIPT (and not a switch in 12): nothing is shared with the
# reach pipeline except the estimator and the analysis. Different spikes
# (st_init), different kinematics (push/pull windows, 2401 frames), different
# per-trial labels (idx_init_vPP_invPP, not idx_reach_*), different t=0. The
# ANALYSIS, though, must be identical for the two alignments to be comparable,
# so cells 5-7 are 12's own cells, exec'd (same trick 12 uses on 09) -- they are
# never copied here, and cannot drift.
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import re, time, pickle
from scipy.signal import savgol_filter
from scipy.stats import mannwhitneyu, wilcoxon
from joblib import Parallel, delayed

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
import mi_io, config, mi_rates

# ── the kinematic feature vector (joint MI is taken over ALL of these) ───────
# "force" is the LOAD CELL, resampled by 16 onto the camera grid, so it is one
# more column of the same (n_events, n_kin) shape as x/y/z/speed. Units do not
# matter: _prep scales every column to unit SD before the Chebyshev norm.
FEATURES_PP = ("x", "y", "z", "speed")
                       # any of: x, y, z, speed, accel, force
PAW         = "dom"   # which paw gives x/y/z/speed ("dom" or "nondom").
                       # The pushing paw need not be the reaching one, so this is
                       # worth running both ways; it is part of the cache key.

# ── which events ─────────────────────────────────────────────────────────────
# Codes from idx_init_vPP_invPP: 0 = valid push, 1 = valid pull,
#                                2 = invalid push, 3 = invalid pull.
# The invalid ones (2, 3) can be included: 21 removes those detected while the
# mouse was reaching (USE_CLEAN) and gives every event a force onset (ALIGN).
#
# Both valid codes together is the default: it is the full set of trial
# initiations and it keeps the event count per session high (median 66 push + 77
# pull, so ~140 rather than ~70).
# ONE THING TO KEEP IN MIND when reading a pooled run WITH "force" among the
# features: the force sign is what distinguishes the two conditions (push
# deflects +1800..+4200 mV/V, pull -4700..-11400 across the 29 sessions), so part
# of that MI is information about WHICH rule the trial was, not about graded
# force. A single-code run (PP_KEEP = (0,) or (1,)) removes that by construction
# and is the control to run if the pooled maps look force-dominated. The
# alignment figures in cells 2b-2d are always drawn PER CODE, so the size of the
# push/pull difference in every feature is visible before any MI is computed.
PP_KEEP = (0, 1, 2, 3)   # any subset of (0, 1, 2, 3)

# ── where t = 0 is ───────────────────────────────────────────────────────────
# "cue"    the timestamp of the event as the rig recorded it = the moment it was
#          DETECTED, which is when the cue fires (valid: sound + bar vibration;
#          invalid: blue light). The right zero for anything driven by the cue.
#          Valid and invalid are detected at DIFFERENT points of the force rise
#          (valid waits for a 50 ms hold; invalid does not), so pooling them on
#          this zero mixes two phases of the movement.
# "onset"  the FORCE ONSET from 21: when the force starts rising, defined the
#          same way for every event. Valid and invalid line up on it (21's
#          alignment figure). Events without an onset are left out.
# Both the kinematics/force AND the spikes are moved to the chosen zero.
ALIGN     = "onset"      # "cue" | "onset"
USE_CLEAN = True       # use 21's cleaning: invalid events detected while the
                       # mouse was reaching are left out. Valid events are never
                       # removed by it (they started a trial). Only matters when
                       # PP_KEEP contains 2 or 3.
CLEAN_PKL = "pushpull_clean.pkl"     # written by 21

# ── which mice ───────────────────────────────────────────────────────────────
EXCLUDE_TOBLERONE = False   # True = leave out EVERY session of Toblerone, which
                       # cannot push (its sessions hold 1-63 valid pushes against
                       # 38-111 pulls, and the pushes it does make are not the
                       # same movement). It is removed from the neurons, from the
                       # alignment figures and from the per-mouse maps alike, and
                       # the choice is in the folder, the cache name and the
                       # cache key, so a run with and a run without it never mix.

# ── the KINEMATIC / BEHAVIOUR window (rows) ──────────────────────────────────
KIN_WIN   = 21         # frames per kinematic window (21 @ 400 fps = 52.5 ms,
                       # the same width 09 and 12 use)
KIN_STEP  = 8          # frames between kinematic windows (8 = 20 ms).
                       # Cost is LINEAR in the number of kin windows; 12 uses 4
                       # (10 ms) over a much shorter range.
KIN_RANGE = (-0.20, 0.50)       # (start, stop) in SECONDS re t = 0 (see ALIGN),
                       # for the window CENTRES. 16 extracted [-3, 3] s, so this
                       # can be widened a long way -- but every extra row
                       # multiplies the run. Re the CUE the movement sits before
                       # 0; re the ONSET it starts AT 0 and the force peaks
                       # ~100-200 ms later, so with ALIGN = "onset" something like
                       # (-0.30, 0.40) is the comparable range. Cell 2b shows it.

# ── the NEURAL window (columns) ──────────────────────────────────────────────
WIN_WIDTH   = 0.07    # (s) width of the spike-count window
WIN_STEP    = 0.025    # (s) step between neural windows
WIN_RANGE   = (-0.50, 0.50)     # (s) re the PUSH/PULL, for the window CENTRES.
                       # Cost: the defaults give 40 kinematic x 37 neural windows
                       # = 1480 cells x 21 (obs + 20 nulls) = ~31k KSG calls per
                       # neuron, so ~75M over the ~2400 usable neurons, ~8 h at
                       # the ~2500 calls/s this machine does -- the same order as
                       # 12's full run (which is 31k per neuron too). Cell 3
                       # prints the ETA before the loop starts: read it.

# ── estimator (same as 09/12; not a knob to change casually) ─────────────────
K_KSG       = 3
NEURAL_VAR  = "count"  # "count" | "isi" | "conv" -- see mi_rates.py
FPS         = 400.0
SAVGOL_WIN  = 11       # frames for the Savitzky-Golay derivative (speed/accel)
SAVGOL_ORD  = 3

# ── how much to run ──────────────────────────────────────────────────────────
MIN_TRIALS  = 20       # minimum events a neuron needs to be estimated at all.
                       # One code is only part of a session, so this bites harder
                       # here than in 12. Measured over the 29 sessions (valid
                       # events per session, verified 2026-10-09):
                       #   push  0..245, median 66  -> 25/29 sessions at 20
                       #   pull  0..125, median 77  -> 26/29 sessions at 20
                       # The sessions that fail are the ones in a block of the
                       # other rule: 1_CoteDor|R7 has 0 pushes, 3_Toblerone|R6
                       # has 1, 4_Milka|R7 has 0 pulls. Cell 3 prints which
                       # survive. As in 12 the excess MI stays unbiased at low N;
                       # what shrinks is sensitivity.
N_TEST      = None     # neurons PER REGION. None = whole population.
SESSIONS    = None     # None = every session, or [(mouse, sess), ...] to restrict
                       # BOTH the st_init parse and the kinematics. With a subset
                       # the st_init cache is neither read nor written, so a quick
                       # one-session pass cannot leave a partial cache behind.
N_JOBS      = -1

DO_NULL     = True     # per-CELL shuffle null -> Excess MI (obs - null)
N_SHUF_2D   = 20       # permutations per cell. MULTIPLIES the whole run.

# ── diagonal vs neural time (cell 5, inherited from 12) ──────────────────────
LAG_BIN    = None      # (s) lag resolution; None = WIN_STEP
CV_GAP     = None      # kinematic windows dropped next to the CV split;
                       # None = the overlap, KIN_WIN // KIN_STEP
FIG_SUBSET = "all"     # "all" or "structured"

# ── caching ──────────────────────────────────────────────────────────────────
CACHE_MODE  = "load"   # "load" = reuse a cached run whose parameters MATCH
                       # "rerun" = always recompute
SAVE_CACHE  = True

# ── empty trials ─────────────────────────────────────────────────────────────
# A trial with no spikes is an EMPTY trial here, and it is KEPT as a trial (it
# counts towards MIN_TRIALS and gives a count of 0 in every window -- silence is
# informative, so silent trials are never dropped).
# Until 2026-10-12 the reader returned them as the spike times [0.0, 1.0] (the
# dimensions of MATLAB's 0x1 empty array, see mi_io._spikes_from_ref and
# 19_spike_placeholder.py). That is fixed in mi_io, which this script reads
# st_init through, so there is no switch for it here.

# ── the alignment checks (cells 2b-2d, all cheap, all before any MI) ─────────
CHECK_FIGS  = True     # draw them at all
EX_EVENTS   = 6        # single events per code in the single-event figure, taken
                       # from AS MANY DIFFERENT MICE/SESSIONS as possible (round
                       # robin over mice, then sessions), so one animal's quirk
                       # cannot pass for the alignment
PETH_BIN    = 0.010    # (s) bin for the spike PETH around the push/pull
PETH_WIN    = (-1.5, 1.5)       # (s) span of the PETH
PETH_BASE   = (-1.5, -1.0)      # (s) baseline for the per-neuron z-score

# ── what the figures report ──────────────────────────────────────────────────
N_EXAMPLES   = 40      # individual landscapes PER REGION. None = every neuron.
EXAMPLE_MODE = "spread"        # "spread" | "top" | "random"
PER_MOUSE_FILES = True
PER_MOUSE_MIN_N = 5
SAVE_PDF_ALL_MICE = True        # also save the ALL MICE mean-landscape figure as
                                # PDF (text kept editable), as 12 does
EQUAL_AXES   = False

SHOW_FIGS   = True
SAVE_FIGS   = True
SAVE_VARS   = True
CMAP        = "viridis"
OUT_SUBDIR  = "MI_2D_pushpull"
RUN_NAME    = "pp2D"

PP_PKL      = "pushpull_kin.pkl"      # written by 16
INIT_PKL    = "parsed_init.pkl"       # st_init cache, written by cell 1b below

CLR     = {"CP": np.array([106, 176, 134]) / 255,
           "CB": np.array([154, 123, 193]) / 255}
PP_CODE = {0: "valid push", 1: "valid pull", 2: "invalid push", 3: "invalid pull"}
PP_CLR  = {0: np.array([132, 160, 124]) / 256,     # the project's push colour
           1: np.array([26, 94, 99]) / 256}        # and pull colour
PP_CLR[2] = PP_CLR[0] * 0.55 + 0.45                 # invalid: the same, lighter
PP_CLR[3] = PP_CLR[1] * 0.55 + 0.45
PP_LS   = {0: "-", 1: "--", 2: ":", 3: "-."}
if ALIGN not in ("cue", "onset"):
    raise ValueError(f"ALIGN must be 'cue' or 'onset', got {ALIGN!r}")
ZERO_TXT = "the cue (detection)" if ALIGN == "cue" else "the force onset"
_INVALID = any(c in (2, 3) for c in PP_KEEP)
CLEAN_ON = USE_CLEAN and _INVALID
NAT2BIT = 1 / np.log(2)

PP_TAG  = "pp" + "".join(str(c) for c in sorted(PP_KEEP)) + f"_{PAW}"
NW_TAG  = f"nw{int(round(WIN_WIDTH * 1000))}"      # the neural window width is in
MICE_TAG = "_noTob" if EXCLUDE_TOBLERONE else ""
ALIGN_TAG = "_onset" if ALIGN == "onset" else ""    # "cue" keeps the old names
CLEAN_TAG = "_clean" if CLEAN_ON else ""
# the feature set is in the folder name too (as it already is in the cache and
# results names), so a run with other FEATURES_PP never writes over these figures
FEAT_TAG  = "-".join(FEATURES_PP)
OUT_DIR_C = (f"{OUT_SUBDIR}_{FEAT_TAG}_{NW_TAG}_{PP_TAG}{MICE_TAG}{ALIGN_TAG}"
             f"{CLEAN_TAG}")                       # every output name, as in 12

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

_bad = [c for c in PP_KEEP if c not in PP_CODE]
if _bad:
    raise ValueError(f"PP_KEEP has code(s) {_bad} that do not exist. The codes are "
                     + ", ".join(f"{c} = {v}" for c, v in PP_CODE.items())
                     + " -- e.g. PP_KEEP = (0, 1, 2, 3) for all four.")
print(f"push/pull 2D MI  |  codes {[PP_CODE[c] for c in PP_KEEP]}  |  paw {PAW}")
print(f"features {list(FEATURES_PP)}  |  neural var {NEURAL_VAR}, k={K_KSG}"
      + ("  |  Toblerone EXCLUDED" if EXCLUDE_TOBLERONE else ""))
print(f"t = 0 at {ZERO_TXT}"
      + (f"  |  invalid events cleaned by 21" if CLEAN_ON else
         ("  |  invalid events NOT cleaned (USE_CLEAN = False)" if _INVALID else "")))
if len(PP_KEEP) > 1 and "force" in FEATURES_PP:
    print("  note: push AND pull pooled with force among the features -- the "
          "force sign\n  separates the two conditions, so part of this MI is "
          "information about the\n  trial TYPE. A single-code run is the "
          "control for that.")


# %% ── 1. the estimator (same code as 09, imported not copied) ──────────────
# 1. cell --> pull _prep/_ksg_md out of 09 so the estimator is literally the same
# function 09 and 12 use. Only 09's cell 1 is exec'd: it defines helpers and
# touches no data, so nothing about the reach pipeline is loaded.
from sklearn.neighbors import NearestNeighbors, KDTree
from scipy.special import digamma

_src9 = open(os.path.join(_HERE, "09_MI_multiD_kin.py"), encoding="utf-8").read()
_m = re.search(r"(?m)^def _prep\(", _src9)
_e = re.search(r"(?m)^def _landscape_one\(", _src9)
if not (_m and _e):
    raise RuntimeError("could not find _prep.._ksg_md in 09 -- renamed?")
exec(compile(_src9[_m.start():_e.start()], "09_estimator", "exec"))
print("estimator: 09's _prep + _ksg_md (KSG-1, Chebyshev, true multivariate)")


# %% ── 1b. st_init spike times ──────────────────────────────────────────────
# 1b. cell --> spikes aligned to the PUSH/PULL.
# mi_io's parser keeps only st_reach (the reach cache would have to be rebuilt to
# add st_init, invalidating every existing analysis), so st_init is read here
# into its OWN pickle. mi_io's private readers are reused, not reimplemented, so
# the two caches are parsed by identical code.
#   INIT[uid] = dict(region, mouse, sess, phyID, spikes=[array per event],
#                    code=idx_init_vPP_invPP)
_IP = config.CACHE_DIR / INIT_PKL


def _parse_init(path, regions, keys=None):
    """-> {uid: dict(...)}. Mirrors mi_io._parse_mat_h5/_parse_mat_scipy, reading
    st_init + idx_init_vPP_invPP instead of st_reach + idx_reach_*.
    keys: restrict to these (mouse, sess); None = all. The mouse/sess strings are
    read FIRST so a restricted parse never touches the other sessions' spikes --
    that is what makes a one-session pass cheap on a 9 GB file."""
    out = {}
    keys = None if keys is None else {tuple(k) for k in keys}
    n_emp = [0, 0]                     # [empty trials, trials seen]

    def _clean(trials):
        """count the empty trials; the list itself is returned unchanged"""
        n_emp[0] += sum(1 for t in trials if t.size == 0)
        n_emp[1] += len(trials)
        return trials

    if mi_io._is_hdf5(path):
        import h5py
        with h5py.File(str(path), "r") as f:
            for reg in regions:
                if reg not in f: print("missing:", reg); continue
                g = f[reg]
                if "st_init" not in g:
                    raise KeyError(f"{reg} has no st_init -- fields: {list(g)}")
                n_neu = g["st_init"].shape[1]
                for i in range(n_neu):
                    if i % 200 == 0: print(f"  {reg}: {i}/{n_neu}", flush=True)
                    if keys is not None and (
                            mi_io._str_from_ref(f, g["mouse"][0, i]),
                            mi_io._str_from_ref(f, g["sess"][0, i])) not in keys:
                        continue
                    spk = _clean(mi_io._spikes_from_ref(f, g["st_init"][0, i]))
                    cd  = (mi_io._label_from_ref(f, g["idx_init_vPP_invPP"][0, i],
                                                 len(spk))
                           if "idx_init_vPP_invPP" in g
                           else np.full(len(spk), np.nan))
                    mouse = mi_io._str_from_ref(f, g["mouse"][0, i])
                    sess  = mi_io._str_from_ref(f, g["sess"][0, i])
                    pid   = int(np.asarray(f[g["phyID"][0, i]], float).ravel()[0])
                    r = mi_io._short(reg)
                    out[f"{r}|{mouse}|{sess}|{pid}"] = dict(
                        region=r, mouse=mouse, sess=sess, phyID=pid,
                        spikes=spk, code=cd)
    else:
        import scipy.io as sio
        mat = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
        for reg in regions:
            if reg not in mat: print("missing:", reg); continue
            for neu in np.atleast_1d(mat[reg]):
                if keys is not None and (str(neu.mouse), str(neu.sess)) not in keys:
                    continue
                if not hasattr(neu, "st_init"):
                    raise KeyError(f"{reg} neurons have no st_init field")
                spk = _clean([mi_io._spk(t)
                              for t in mi_io._trials_list(neu.st_init)])
                cd  = np.atleast_1d(np.asarray(
                    getattr(neu, "idx_init_vPP_invPP", np.nan), float)).ravel()
                if cd.size != len(spk): cd = np.full(len(spk), np.nan)
                r = mi_io._short(reg)
                out[f"{r}|{neu.mouse}|{neu.sess}|{int(neu.phyID)}"] = dict(
                    region=r, mouse=str(neu.mouse), sess=str(neu.sess),
                    phyID=int(neu.phyID), spikes=spk, code=cd)
    print(f"  empty trials: {n_emp[0]}/{n_emp[1]} "
          f"({100*n_emp[0]/max(n_emp[1],1):.1f}%) -- kept, as trials with "
          f"no spikes")
    return out


_FULL = SESSIONS is None          # only a full parse may use/write the cache
if CACHE_MODE == "load" and _FULL and _IP.is_file():
    with _step("st_init cache read"):
        with open(_IP, "rb") as fh: INIT = pickle.load(fh)
    print(f"loaded {_IP.name}: {len(INIT)} neurons -- no re-parse")
else:
    with _step("parse st_init from the .mat"):
        INIT = _parse_init(config.DATA_MAT, config.struct_names(), SESSIONS)
    if _FULL:
        config.CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with open(_IP, "wb") as fh: pickle.dump(INIT, fh, protocol=4)
        print(f"parsed {len(INIT)} neurons -> {_IP}")
    else:
        print(f"parsed {len(INIT)} neurons for SESSIONS={SESSIONS} "
              f"(not cached: a subset must not masquerade as the full parse)")

_nev = {}
for v in INIT.values(): _nev.setdefault((v["mouse"], v["sess"]), len(v["spikes"]))
print(f"{len(INIT)} neurons, {len(_nev)} sessions, "
      f"{min(_nev.values())}..{max(_nev.values())} init events per session")


# %% ── 1c. push/pull kinematics from 16, and the JOIN CHECK ─────────────────
# 1c. cell --> load 16's pickle and verify that its events and st_init's events
# are THE SAME events in THE SAME order. This is the one assumption the whole
# script rests on: 16 indexes events by bhv.pp_times (behaviour mats) while the
# spikes are indexed by st_init (ephys mat). Nothing guarantees the two agree, so
# it is checked per session and printed, not assumed.
# VERIFIED 2026-10-09: all 29 sessions agree index-for-index -- same number of
# events AND the same idx_init_vPP_invPP / valPushPull_invalPushPull_idx vector,
# NaNs included. The check stays in because that is a property of the current
# mats, not a guarantee.
_PK = config.var_path(PP_PKL)
if not _PK.is_file():
    raise FileNotFoundError(
        f"{_PK} not found. Run 16_pushpull_kinematics.py with SESSIONS = None "
        f"and SAVE_VARS = True first -- that is what writes the push/pull "
        f"kinematics this script reads.")
with _step("pushpull_kin.pkl read"):
    with open(_PK, "rb") as fh: _PPK = pickle.load(fh)
PP      = _PPK["sessions"]
T_PP    = np.asarray(_PPK["t"], float)          # (2401,) s re the push/pull
PP_UNIT = _PPK["units"]
if f"kin_{PAW}" not in next(iter(PP.values())):
    raise KeyError(f"16's pickle has no kin_{PAW}; it holds "
                   f"{[k for k in next(iter(PP.values())) if k.startswith('kin_')]}")
print(f"16's extraction: {len(PP)} sessions, {T_PP.size} frames "
      f"{T_PP[0]:+.2f}..{T_PP[-1]:+.2f} s, paws in {PP_UNIT}")

_Z = int(np.argmin(np.abs(T_PP)))               # the frame where t = 0
if abs(T_PP[_Z]) > 0.5 / FPS:
    raise ValueError(f"16's time axis has no zero frame (closest {T_PP[_Z]:+.4f} s)")

# ── the join check ───────────────────────────────────────────────────────────
print(f"\n{'session':26s} {'st_init':>8s} {'16 events':>10s} {'codes':>8s} "
      f"{'kept':>6s}")
print("-" * 68)
JOIN, _bad = {}, []
for key in sorted(set(_nev) & set(PP)):   # _nev already honours SESSIONS
    n_i, n_p = _nev[key], PP[key]["t"].size
    c_p = np.asarray(PP[key]["code"], float).ravel()
    c_i = next(np.asarray(v["code"], float).ravel() for v in INIT.values()
               if (v["mouse"], v["sess"]) == key)
    same = (n_i == n_p and c_i.size == c_p.size
            and np.array_equal(np.nan_to_num(c_i, nan=-9),
                               np.nan_to_num(c_p, nan=-9)))
    kept = int(np.isin(c_p, list(PP_KEEP)).sum())
    JOIN[key] = same
    if not same: _bad.append(key)
    print(f"{key[0]+'|'+key[1]:26s} {n_i:8d} {n_p:10d} "
          f"{'same' if same else 'DIFFER':>8s} {kept:6d}")
_miss = sorted(set(_nev) - set(PP))
if _miss:
    print(f"\n{len(_miss)} session(s) in the ephys mat but not in 16's pickle: {_miss}")
if _bad:
    print(f"\n!! {len(_bad)}/{len(JOIN)} session(s) where st_init and 16 do NOT "
          f"describe the same events. Those are EXCLUDED below -- the spikes could "
          f"not be paired with the kinematics. Fix before trusting anything here.")
else:
    print(f"\nall {len(JOIN)} shared sessions: st_init and 16 agree on the number "
          f"of events AND on every code -> the pairing is index-for-index.")
USE_KEYS = [k for k, ok in JOIN.items() if ok]
if EXCLUDE_TOBLERONE:
    _tob = [k for k in USE_KEYS if "toblerone" in k[0].lower()]
    USE_KEYS = [k for k in USE_KEYS if k not in _tob]
    print(f"EXCLUDE_TOBLERONE: {len(_tob)} session(s) left out "
          f"({', '.join(k[0] + '|' + k[1] for k in _tob) or 'none found'}); "
          f"{len(USE_KEYS)} sessions from "
          f"{len({k[0] for k in USE_KEYS})} mice remain")

# ── which events are usable, and where t = 0 is (both from 21) ───────────────
# EV_OK[key] marks, per event, whether it may be used: 21's cleaning (if on) and,
# for the onset alignment, an onset that could be measured. SHIFT[key] is where
# the new zero sits re the cue (0 for the cue alignment).
EV_OK = {k: np.ones(np.asarray(PP[k]["code"]).size, bool) for k in USE_KEYS}
SHIFT = {k: np.zeros(np.asarray(PP[k]["code"]).size) for k in USE_KEYS}
if ALIGN == "onset" or CLEAN_ON:
    _CP = config.var_path(CLEAN_PKL)
    if not _CP.is_file():
        raise FileNotFoundError(f"{_CP} not found -- run 21_pushpull_clean_align.py "
                                f"(it is needed for ALIGN='onset' and for USE_CLEAN)")
    with open(_CP, "rb") as fh: _CL = pickle.load(fh)["sessions"]
    _miss = [k for k in USE_KEYS if k not in _CL]
    if _miss:
        raise KeyError(f"{CLEAN_PKL} has no entry for {_miss} -- was 21 run with "
                       f"EXCLUDE_TOBLERONE or SESSIONS? Run it on every session.")
    for k in USE_KEYS:
        c21, c16 = _CL[k], PP[k]
        # 21 must describe THESE events: same times, same codes, same order
        if not (np.array_equal(np.asarray(c21["code"]), np.asarray(c16["code"]))
                and np.allclose(np.asarray(c21["t"]), np.asarray(c16["t"]))):
            raise ValueError(f"{k}: {CLEAN_PKL} and {PP_PKL} describe different events "
                             f"-- rerun 21 after 16")
        if CLEAN_ON:          EV_OK[k] &= np.asarray(c21["keep"], bool)
        if ALIGN == "onset":
            sh = np.asarray(c21["onset_shift"], float)
            EV_OK[k] &= np.isfinite(sh)
            SHIFT[k] = np.where(np.isfinite(sh), sh, 0.0)
    print(f"\nevents per code (kept sessions): total / usable"
          + (" (21's cleaning)" if CLEAN_ON else "")
          + (" (with an onset)" if ALIGN == "onset" else ""))
    for c in sorted(PP_KEEP):
        _n  = sum(int((np.asarray(PP[k]["code"]) == c).sum()) for k in USE_KEYS)
        _ok = sum(int(((np.asarray(PP[k]["code"]) == c) & EV_OK[k]).sum()) for k in USE_KEYS)
        print(f"  {PP_CODE[c]:13s} {_ok:5d} / {_n:5d}")

if ALIGN == "onset":
    # Move every trace so that t = 0 is the event's onset: new(t) = old(t + shift).
    # 21 measures the onset ON 16's frame grid, so the shift is a whole number of
    # frames and nothing is interpolated. Done on a COPY: _PPK keeps the original.
    def _shift_rows(A, fr):
        """A (n, ..., T); fr (n,) int frames -> A moved so that index i holds
        the old index i + fr (NaN where that falls outside the window)."""
        out = np.full(A.shape, np.nan, dtype=float)
        T_ = A.shape[-1]
        for e, s_ in enumerate(fr):
            if s_ >= 0: out[e, ..., :T_ - s_] = A[e, ..., s_:]
            else:       out[e, ..., -s_:]     = A[e, ..., :T_ + s_]
        return out
    _fps = float(_PPK["fps"])
    PP = {k: dict(v) for k, v in PP.items()}
    for k in USE_KEYS:
        fr = np.round(SHIFT[k] * _fps).astype(int)
        if not np.allclose(fr / _fps, SHIFT[k], atol=1e-6):
            raise ValueError(f"{k}: onset shifts are not on 16's frame grid")
        for nm in [f for f in PP[k] if f.startswith("kin_")] + ["force"]:
            PP[k][nm] = _shift_rows(np.asarray(PP[k][nm], float), fr)
    print(f"kinematics and force moved to the force onset "
          f"(median shift {np.median(np.concatenate([SHIFT[k][EV_OK[k]] for k in USE_KEYS]))*1000:+.0f} ms re the cue)")


def _spikes(u):
    """the neuron's spike times per event, with t = 0 at the chosen zero. The
    list keeps one entry per event (unusable events included), so event k still
    pairs with kinematics k."""
    n = INIT[u]; sh = SHIFT[(n["mouse"], n["sess"])]
    if ALIGN == "cue": return n["spikes"]
    return [np.asarray(t, float) - sh[e] for e, t in enumerate(n["spikes"])]


# %% ── 2. kinematic features for EVERY window, per session ───────────────────
# 2. cell --> KIN[key] = (n_events, n_kin, n_feats), built once; the MI loop only
# indexes it. Derivatives are computed on the FULL 6-s trace and then averaged
# inside each window (never differentiated per-window), so a short window does
# not distort the Savitzky-Golay filter -- same rule as 09/12.
kin_starts = np.arange(0, T_PP.size - KIN_WIN + 1, KIN_STEP)
# t=0 is frame _Z of 16's axis, exactly (16 cut N_PRE frames before the event),
# so a window's centre is its mid-frame measured from there. No fitted shift and
# no inherited endpoint frame is involved: this alignment has one zero only.
kin_centers = (kin_starts + KIN_WIN / 2 - _Z) / FPS
if KIN_RANGE is not None:
    _in = (kin_centers >= KIN_RANGE[0]) & (kin_centers <= KIN_RANGE[1])
    if not _in.any():
        raise ValueError(f"KIN_RANGE={KIN_RANGE} s selects no window; 16 extracted "
                         f"{T_PP[0]:+.2f}..{T_PP[-1]:+.2f} s")
    kin_starts, kin_centers = kin_starts[_in], kin_centers[_in]
n_kin = len(kin_starts)

win_starts  = np.arange(WIN_RANGE[0], WIN_RANGE[1] - WIN_WIDTH + WIN_STEP/2, WIN_STEP)
win_centers = win_starts + WIN_WIDTH / 2
n_wins      = len(win_starts)

print(f"\n2D landscape: {n_kin} kinematic windows x {n_wins} neural windows "
      f"= {n_kin*n_wins} cells per neuron")
print(f"  KIN : win {KIN_WIN} fr ({KIN_WIN/FPS*1000:.1f} ms), step {KIN_STEP} fr "
      f"({KIN_STEP/FPS*1000:.1f} ms), centres "
      f"{kin_centers[0]*1000:.0f}..{kin_centers[-1]*1000:.0f} ms")
print(f"  NEUR: win {WIN_WIDTH*1000:.0f} ms, step {WIN_STEP*1000:.0f} ms, centres "
      f"{win_centers[0]*1000:.0f}..{win_centers[-1]*1000:.0f} ms")


def _savgol_nan(A, deriv):
    """Savitzky-Golay along the LAST axis, NaN-SAFE.

    Unlike reach_px, these traces DO contain NaN: 16's +-3 s window runs off the
    end of the session trace for events near its edges, and the paw can be
    untracked. scipy fits a polynomial to the edges of the array and raises
    (LinAlgError: SVD did not converge) the moment a NaN is in it, so the filter
    is applied separately to each contiguous run of finite samples that is at
    least SAVGOL_WIN long. Everything else stays NaN, which is honest: a
    derivative across a gap is not measured, it is invented. The MI mask then
    drops those windows."""
    A    = np.asarray(A, float)
    flat = A.reshape(-1, A.shape[-1])
    out  = np.full_like(flat, np.nan)
    for i in range(flat.shape[0]):
        ok = np.isfinite(flat[i])
        if not ok.any(): continue
        edge = np.flatnonzero(np.diff(np.r_[0, ok.astype(np.int8), 0]) != 0)
        for a, b in zip(edge[::2], edge[1::2]):
            if b - a >= SAVGOL_WIN:
                out[i, a:b] = savgol_filter(flat[i, a:b], SAVGOL_WIN, SAVGOL_ORD,
                                            deriv=deriv, delta=1/FPS)
    return out.reshape(A.shape)


# Only the frames the selected windows need are filtered, plus SAVGOL_WIN of
# margin on each side so every USED frame still gets a full filter window. 16
# extracted 2401 frames; the windows typically span a few hundred, and the
# NaN-safe filter is a per-trace python loop, so this is what keeps cell 2 fast.
_f0 = int(max(0, kin_starts.min() - SAVGOL_WIN))
_f1 = int(min(T_PP.size, kin_starts.max() + KIN_WIN + SAVGOL_WIN))
_ks = kin_starts - _f0                       # window starts within the slice
print(f"filtering frames {_f0}..{_f1} of {T_PP.size} "
      f"({T_PP[_f0]:+.3f}..{T_PP[_f1-1]:+.3f} s) -- the windows plus "
      f"{SAVGOL_WIN} frames of margin")


def _win_mean(A, s):          # A (n, T) -> mean over frames [s, s+KIN_WIN)
    return A[:, s:s+KIN_WIN].mean(axis=1)


def _stack(d, feats):
    """-> (n, n_kin, len(feats)) of window means"""
    return np.stack([np.stack([_win_mean(d[f], s) for s in _ks], axis=1)
                     for f in feats], axis=2)


_t_kin = time.time()
KIN, NAN_FRAC = {}, {}
for key in USE_KEYS:
    P   = PP[key]
    rp  = np.asarray(P[f"kin_{PAW}"], float)[:, :, _f0:_f1]    # (n, 3, T') in mm
    vel = _savgol_nan(rp, 1)
    spd = np.linalg.norm(vel, axis=1)
    acc = _savgol_nan(spd, 1)
    frc = np.asarray(P["force"], float)[:, _f0:_f1]            # (n, T') in mV/V
    d = {"x": rp[:, 0], "y": rp[:, 1], "z": rp[:, 2], "speed": spd, "accel": acc,
         "force": frc}
    miss = [f for f in FEATURES_PP if f not in d]
    if miss: raise KeyError(f"unknown feature(s) {miss}; available {sorted(d)}")
    KIN[key] = _stack(d, FEATURES_PP)                          # (n, n_kin, d)
    NAN_FRAC[key] = {f: float(np.mean(~np.isfinite(_stack(d, (f,)))))
                     for f in FEATURES_PP}
TIMES["kinematic stack"] = time.time() - _t_kin
print(f"kinematic stack built in {TIMES['kinematic stack']:.0f}s: "
      f"{KIN[USE_KEYS[0]].shape} per session (events, kin_wins, features)")

# Per-feature NaN, so a session whose LOAD CELL failed is visible BEFORE it
# silently removes every one of its neurons (the MI mask needs all features
# finite, so one all-NaN feature drops the whole session).
print(f"\n{'session':26s} " + " ".join(f"{f:>10s}" for f in FEATURES_PP)
      + "    (% of windows NOT finite)")
print("-" * (27 + 11*len(FEATURES_PP)))
for key in USE_KEYS:
    row = " ".join(f"{100*NAN_FRAC[key][f]:9.0f}%" for f in FEATURES_PP)
    flag = "   <-- a feature is entirely missing" \
           if max(NAN_FRAC[key].values()) > 0.99 else ""
    print(f"{key[0]+'|'+key[1]:26s} {row}{flag}")


# %% ── 2b. ALIGNMENT CHECK 1: what the kinematic rows read ──────────────────
# 2b. cell --> before paying for the MI, look at the kinematics and the force the
# kinematic windows are reading, split by code. Left block: pooled over every
# session (mean +- SEM). Right block: the same means PER MOUSE, so one animal
# cannot carry it. Bottom right: how much of each window is actually tracked.
# The shaded band is KIN_RANGE and the ticks are the window centres -- if the
# movement and the force deflection do not sit inside that band, move KIN_RANGE
# rather than running anything.
_codes = sorted(set(PP_KEEP) | {0, 1})
_show  = (T_PP >= min(KIN_RANGE[0] - 0.4, -1.0)) & (T_PP <= max(KIN_RANGE[1] + 0.4, 0.5))
_TS    = T_PP[_show]
_QUANT = (("pos",   lambda rp, fr: rp[:, 0],                        f"{PAW} paw x ({PP_UNIT})"),
          ("speed", lambda rp, fr: np.linalg.norm(_savgol_nan(rp, 1), axis=1),
                                                                    f"{PAW} paw speed ({PP_UNIT}/s)"),
          ("force", lambda rp, fr: fr,                              "load cell (mV/V)"))


def _traces(key, c, fn):
    """one session's traces for code c, already sliced to the plotted range.
    Sliced BEFORE filtering on purpose: the NaN-safe filter is a python loop over
    traces, so filtering all 2401 frames of every event of every session just to
    draw 2 s of it is the one thing that would make this check slow."""
    m = (np.asarray(PP[key]["code"], float).ravel() == c) & EV_OK[key]
    if not m.any(): return None
    return fn(np.asarray(PP[key][f"kin_{PAW}"], float)[m][:, :, _show],
              np.asarray(PP[key]["force"], float)[m][:, _show])


def _band(ax, ticks=True):
    ax.axvspan(KIN_RANGE[0], KIN_RANGE[1], color="0.9", zorder=0, lw=0)
    if ticks:
        ax.plot(kin_centers, np.full(n_kin, ax.get_ylim()[1]), "|", color="0.35",
                ms=6, clip_on=False)
    ax.axvline(0, color="orange", lw=1.6, ls="--")
    # the zero line is drawn only where the data reach zero. Drawing it always
    # forces the axis down to 0, which squashes a trace that lives far from it
    # (a paw position around 20 mm) into the top of the panel.
    # (judged on the DATA range, not on the axis limits, which carry a margin)
    _lo, _hi = ax.dataLim.intervaly
    if _lo <= 0 <= _hi:
        _yl = ax.get_ylim()
        ax.axhline(0, color="k", lw=0.6)
        ax.set_ylim(_yl)
    ax.spines[["top", "right"]].set_visible(False)


MICE_PP = sorted({k[0] for k in USE_KEYS})
_MSTY   = ["-", "--", ":", "-.", (0, (3, 1, 1, 1))]     # one per mouse

if CHECK_FIGS:
    fig, axes = plt.subplots(2, 3, figsize=(18, 8.4))
    for col, (nm, fn, lb) in enumerate(_QUANT):
        # ── row 0: pooled over sessions ──────────────────────────────────────
        ax = axes[0, col]
        for c in _codes:
            A = [a for key in USE_KEYS if (a := _traces(key, c, fn)) is not None]
            if not A: continue
            A  = np.concatenate(A)
            mn = np.nanmean(A, axis=0)
            se = np.nanstd(A, axis=0) / np.sqrt(np.maximum(np.isfinite(A).sum(0), 1))
            ax.plot(_TS, mn, color=PP_CLR[c], lw=2,
                    label=f"{PP_CODE[c]} (n={A.shape[0]})")
            ax.fill_between(_TS, mn - se, mn + se, color=PP_CLR[c], alpha=0.2)
        ax.set(xlabel=f"time from {ZERO_TXT} (s)", ylabel=lb,
               title="all sessions pooled (mean ± SEM)")
        ax.title.set_fontsize(9)
        _band(ax); ax.legend(fontsize=7)

        # ── row 1: one line per MOUSE ────────────────────────────────────────
        ax = axes[1, col]
        for c in _codes:
            for mi, m in enumerate(MICE_PP):
                A = [a for key in USE_KEYS if key[0] == m
                     and (a := _traces(key, c, fn)) is not None]
                if not A: continue
                A = np.concatenate(A)
                if A.shape[0] < 5: continue
                ax.plot(_TS, np.nanmean(A, axis=0), color=PP_CLR[c], lw=1.4,
                        ls=_MSTY[mi % len(_MSTY)], alpha=0.85,
                        label=f"{PP_CODE[c]} {m.split('_')[-1]}")
        ax.set(xlabel=f"time from {ZERO_TXT} (s)", ylabel=lb,
               title="mean per mouse (every animal should agree in SIGN)")
        ax.title.set_fontsize(9)
        _band(ax); ax.legend(fontsize=6, ncol=2)
    fig.suptitle(f"ALIGNMENT 1 / what the KINEMATIC rows read  —  grey = KIN_RANGE "
                 f"{KIN_RANGE} s, ticks = the {n_kin} window centres\n"
                 + (f"t=0 is the cue (detection), i.e. the END of the push/pull, so "
                    f"the movement is BEFORE it and the force must deflect AT it "
                    f"(up for a push, down for a pull)" if ALIGN == "cue" else
                    f"t=0 is the force ONSET: the force must start rising AT it, in "
                    f"the same way for valid and invalid events (up for a push, "
                    f"down for a pull)"), fontsize=11)
    plt.tight_layout(rect=(0, 0, 1, 0.94))
    if SAVE_FIGS:
        _savefig_retry(fig, config.fig_dir(OUT_DIR_C) / f"{RUN_NAME}_align1_means.png",
                       dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 2c. ALIGNMENT CHECK 2: single events ──────────────────────────────────
# 2c. cell --> means can be stereotyped for the wrong reason, so look at single
# events. One panel per event: the paw trace in the code's colour on the left
# axis, the LOAD CELL in grey on the right axis. The force is the signal that
# DEFINES the event, so in every panel it must deflect at 0 -- that, and not the
# paw, is what says the alignment is right. Events are drawn from as many
# different mice and sessions as possible.
def _pick_events(c, n):
    """-> [(key, event index)] for code c, round-robin over MICE then sessions,
    so the examples are spread over animals instead of being the first n events
    of the first session. Deterministic (seed 0)."""
    rng, pool = np.random.default_rng(0), {}
    for key in USE_KEYS:
        ok = np.flatnonzero((np.asarray(PP[key]["code"], float).ravel() == c)
                            & EV_OK[key] & np.isfinite(KIN[key]).all(axis=(1, 2)))
        for i in ok: pool.setdefault(key[0], []).append((key, int(i)))
    for m in pool:
        pool[m] = [pool[m][j] for j in rng.permutation(len(pool[m]))]
    out, r = [], 0
    while len(out) < n and any(len(v) > r for v in pool.values()):
        for m in sorted(pool):
            if len(pool[m]) > r and len(out) < n: out.append(pool[m][r])
        r += 1
    return out


if CHECK_FIGS:
    _EV = {c: _pick_events(c, EX_EVENTS) for c in _codes}
    _nr = max(len(v) for v in _EV.values()) or 1
    fig, axes = plt.subplots(_nr, len(_codes),
                             figsize=(7.0 * len(_codes), 2.1 * _nr),
                             squeeze=False, sharex=True)
    for col, c in enumerate(_codes):
        for r in range(_nr):
            ax = axes[r, col]
            if r >= len(_EV[c]):
                ax.axis("off"); continue
            key, i = _EV[c][r]
            rp = np.asarray(PP[key][f"kin_{PAW}"], float)[i][:, _show]
            fr = np.asarray(PP[key]["force"], float)[i][_show]
            ax.plot(_TS, rp[0], color=PP_CLR[c], lw=1.6)
            ax.set_ylabel(f"{key[0].split('_')[-1]}|{key[1]} #{i}\n"
                          f"x ({PP_UNIT})", fontsize=7)
            ax.tick_params(labelsize=7)
            a2 = ax.twinx()
            a2.plot(_TS, fr, color="0.45", lw=1.1)
            a2.axhline(0, color="0.75", lw=0.5)
            a2.set_ylabel("force", fontsize=6, color="0.45")
            a2.tick_params(labelsize=6, colors="0.45")
            ax.axvspan(KIN_RANGE[0], KIN_RANGE[1], color="0.92", zorder=0, lw=0)
            ax.axvline(0, color="orange", lw=1.4, ls="--")
            ax.spines[["top"]].set_visible(False)
            if r == 0:
                ax.set_title(f"{PP_CODE[c]}  —  paw x (colour) + load cell (grey)",
                             fontsize=10)
            if r == _nr - 1: ax.set_xlabel(f"time from {ZERO_TXT} (s)")
    fig.suptitle(f"ALIGNMENT 2 / single events, {EX_EVENTS} per code from "
                 f"different mice and sessions\n"
                 f"the FORCE (grey) is what defines the event: it must "
                 + ("deflect at 0" if ALIGN == "cue" else "start rising at 0")
                 + " in every panel, whatever the paw is doing", fontsize=11)
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    if SAVE_FIGS:
        _savefig_retry(fig, config.fig_dir(OUT_DIR_C) / f"{RUN_NAME}_align2_single.png",
                       dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

    # the per-session force deflection 16 already measured, as a table figure:
    # 29 numbers that each have to have the right SIGN
    _D = {k: PP[k].get("lc_defl", {}) for k in USE_KEYS}
    fig, ax = plt.subplots(figsize=(max(7, 0.42 * len(USE_KEYS)), 4.2))
    _x = np.arange(len(USE_KEYS))
    for c in (0, 1):
        v = [(_D[k].get(c, (np.nan, 0))[0]
              if isinstance(_D[k].get(c), tuple) else _D[k].get(c, np.nan))
             for k in USE_KEYS]
        ax.bar(_x + (c - 0.5) * 0.4, v, width=0.4, color=PP_CLR[c],
               label=PP_CODE[c])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(_x)
    ax.set_xticklabels([f"{k[0].split('_')[-1]}|{k[1]}" for k in USE_KEYS],
                       rotation=90, fontsize=6)
    ax.set(ylabel="median force deflection at t=0 (mV/V)",
           title="ALIGNMENT 2b / the force deflection 16 measured, per session\n"
                 "push must be POSITIVE and pull NEGATIVE in every session")
    ax.title.set_fontsize(10)
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)
    plt.tight_layout()
    if SAVE_FIGS:
        _savefig_retry(fig, config.fig_dir(OUT_DIR_C) / f"{RUN_NAME}_align2b_force.png",
                       dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 2d. ALIGNMENT CHECK 3: the SPIKES share that zero ────────────────────
# 2d. cell --> the checks above are all about the kinematics, which come from 16.
# The spikes come from a different file (st_init in the ephys mat), and the whole
# 2D map assumes the two share t=0. This is the figure that tests it: the
# population PETH around the push/pull, per region and per code. Cell 1c already
# showed the two files describe the same EVENTS; this shows they describe them on
# the same CLOCK -- a modulation locked to 0 is the evidence, and a PETH whose
# features sit tens of ms off would mean the zeros differ.
# Left: mean rate in Hz. Right: per-neuron z-score against PETH_BASE, so a few
# high-rate CB cells cannot set the shape.
if CHECK_FIGS:
    _pb   = np.arange(PETH_WIN[0], PETH_WIN[1] + PETH_BIN / 2, PETH_BIN)
    _pc   = _pb[:-1] + PETH_BIN / 2
    _base = (_pc >= PETH_BASE[0]) & (_pc <= PETH_BASE[1])
    _t_pe = time.time()
    PETH = {}                      # (region, code) -> (n_neurons, n_bins) in Hz
    for reg in ("CP", "CB"):
        for c in _codes:
            rows = []
            for u, v in INIT.items():
                if v["region"] != reg: continue
                if (v["mouse"], v["sess"]) not in set(USE_KEYS): continue
                m = np.flatnonzero((np.asarray(v["code"], float).ravel() == c)
                                   & EV_OK[(v["mouse"], v["sess"])])
                if m.size < MIN_TRIALS: continue
                _sp = _spikes(u)
                spk = np.concatenate([_sp[i] for i in m]) \
                      if m.size else np.array([])
                h, _ = np.histogram(spk, bins=_pb)
                rows.append(h / (m.size * PETH_BIN))          # Hz
            PETH[(reg, c)] = np.asarray(rows) if rows else np.zeros((0, _pc.size))
    print(f"PETH built in {time.time()-_t_pe:.0f}s: "
          + ", ".join(f"{r}/{PP_CODE[c]} n={PETH[(r,c)].shape[0]}"
                      for r in ("CP", "CB") for c in _codes))
    # The REACH follows the push/pull by about a second (16 saves next_reach per
    # event), and that is a second movement with its own neural correlate. It is
    # marked so a bump out there is recognised rather than read as part of the
    # initiation -- and so it can be checked that it falls OUTSIDE WIN_RANGE.
    _nr = np.concatenate([np.asarray(PP[k]["next_reach"], float).ravel()
                          for k in USE_KEYS])
    _nr = float(np.nanmedian(_nr))

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.6))
    _drop = {}
    for ax, mode in zip(axes, ("hz", "z")):
        for reg in ("CP", "CB"):
            for c in _codes:
                A = PETH[(reg, c)]
                if not A.size: continue
                if mode == "z":
                    mu = A[:, _base].mean(axis=1, keepdims=True)
                    sd = A[:, _base].std(axis=1, keepdims=True)
                    # A neuron that is SILENT throughout the baseline has sd = 0
                    # and no z-score at all: dividing by a floor turns its first
                    # spike into a value of 1e9 and the panel shows nothing else.
                    # Those neurons are dropped from this panel only -- they are
                    # in the Hz panel, and in the MI, like any other.
                    _fin = (sd[:, 0] > 0)
                    _drop[(reg, c)] = int((~_fin).sum())
                    if not _fin.any(): continue
                    A = (A[_fin] - mu[_fin]) / sd[_fin]
                mn = np.nanmean(A, axis=0)
                se = np.nanstd(A, axis=0) / np.sqrt(max(A.shape[0], 1))
                ax.plot(_pc, mn, color=CLR[reg], lw=2,
                        ls=PP_LS[c],
                        label=f"{reg} {PP_CODE[c]} (n={A.shape[0]})")
                ax.fill_between(_pc, mn - se, mn + se, color=CLR[reg], alpha=0.18)
        ax.axvspan(KIN_RANGE[0], KIN_RANGE[1], color="0.9", zorder=0, lw=0,
                   label="KIN_RANGE" if mode == "hz" else None)
        ax.axvspan(WIN_RANGE[0], WIN_RANGE[1], color="none", zorder=0,
                   ec="0.55", ls=":", lw=1.2)
        ax.axvline(0, color="orange", lw=1.6, ls="--")
        if PETH_WIN[0] <= _nr <= PETH_WIN[1]:
            ax.axvline(_nr, color="0.4", lw=1.2, ls="-.",
                       label=f"median reach, +{_nr*1000:.0f} ms"
                             if mode == "hz" else None)
        ax.set(xlabel=f"time from {ZERO_TXT} (s)",
               ylabel="population mean rate (Hz)" if mode == "hz"
                      else f"z re {PETH_BASE} s",
               title="mean firing rate" if mode == "hz"
                     else "per-neuron z-score (shape, not amplitude)")
        ax.title.set_fontsize(10)
        ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=7)
    if any(_drop.values()):
        print("  z panel: dropped "
              + ", ".join(f"{v} {r} {PP_CODE[c].split()[1]}"
                          for (r, c), v in sorted(_drop.items()) if v)
              + " neuron(s) that are silent through the whole baseline "
              f"({PETH_BASE} s) and so have no z-score. They are in the Hz "
              "panel and in the MI unchanged.")
    fig.suptitle(f"ALIGNMENT 3 / the SPIKES around the push/pull  —  st_init, "
                 f"{PETH_BIN*1000:.0f} ms bins, solid = push, dashed = pull\n"
                 f"grey = KIN_RANGE (the rows), dotted outline = WIN_RANGE (the "
                 f"columns); the modulation must be locked to 0, which is what "
                 f"says st_init and 16 share the zero", fontsize=11)
    plt.tight_layout(rect=(0, 0, 1, 0.92))
    if SAVE_FIGS:
        _savefig_retry(fig, config.fig_dir(OUT_DIR_C) / f"{RUN_NAME}_align3_peth.png",
                       dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 3. the 2D landscape, per neuron ───────────────────────────────────────
# 3. cell --> compute the MI landscape for every neuron, and the nulls.
# Identical to 12's cell 3 in everything but where the data comes from: the same
# permutations are reused across cells so the null surface is smooth.
def _landscape2d_one(kin_stack, C_m, k, seed, n_shuf=0):
    """kin_stack (n_events, n_kin, d), C_m (n_events, n_neu)
    -> (mi, null_mean, null_std)"""
    n_kin_, n_neu_ = kin_stack.shape[1], C_m.shape[1]
    mi = np.full((n_kin_, n_neu_), np.nan)
    SH = np.full((n_shuf, n_kin_, n_neu_), np.nan) if n_shuf else None
    rs = np.random.default_rng(seed)
    perms = [rs.permutation(kin_stack.shape[0]) for _ in range(n_shuf)]
    ok_neu = [wi for wi in range(n_neu_) if np.unique(C_m[:, wi]).size >= 2]
    for ki in range(n_kin_):
        X = kin_stack[:, ki, :]
        if not np.isfinite(X).all(): continue
        for wi in ok_neu:
            c = C_m[:, wi]
            mi[ki, wi] = _ksg_md(X, c, k, seed + ki*97 + wi)
            for s, p in enumerate(perms):
                SH[s, ki, wi] = _ksg_md(X[p], c, k, seed + ki*97 + wi)
    if not n_shuf:
        return mi, np.full_like(mi, np.nan), np.full_like(mi, np.nan)
    return mi, np.nanmean(SH, axis=0), np.nanstd(SH, axis=0, ddof=1)


# which neurons to run, in a stable order (region, then the cache's own order)
_all_uids = [u for u in INIT if (INIT[u]["mouse"], INIT[u]["sess"]) in set(USE_KEYS)]
uids, run_idx, regions_run = [], [], []
for reg in ("CP", "CB"):
    sel = [u for u in _all_uids if INIT[u]["region"] == reg]
    if N_TEST is not None: sel = sel[:N_TEST]
    for u in sel:
        run_idx.append(len(uids)); uids.append(u); regions_run.append(reg)
regions_run = np.array(regions_run)
print(f"\nneurons to run: {len(run_idx)} "
      + ", ".join(f"{r} {int((regions_run == r).sum())}" for r in ("CP", "CB")))

_cells = len(run_idx) * n_kin * n_wins * (1 + (N_SHUF_2D if DO_NULL else 0))
print(f"ETA: ~{_cells:,} KSG calls; at ~2500 calls/s parallel that is "
      f"~{_cells/2500/60:.0f} min")

payload, keep = [], []
_t_pay = time.time()
for pi, u in enumerate(uids):
    n   = INIT[u]
    K   = KIN[(n["mouse"], n["sess"])]
    cd  = np.asarray(n["code"], float).ravel()
    msk = (np.isfinite(K).all(axis=(1, 2)) & np.isin(cd, list(PP_KEEP))
           & EV_OK[(n["mouse"], n["sess"])])
    if msk.sum() < MIN_TRIALS: payload.append(None); continue
    C = np.asarray(mi_rates.windowed(_spikes(u), win_starts, WIN_WIDTH,
                                     NEURAL_VAR), float)[msk]
    payload.append((K[msk], C)); keep.append(pi)
TIMES["spike counting + masks"] = time.time() - _t_pay
print(f"{len(keep)}/{len(run_idx)} neurons usable with MIN_TRIALS={MIN_TRIALS} "
      f"({TIMES['spike counting + masks']:.0f}s)")
_ses_ok = {(INIT[uids[pi]]["mouse"], INIT[uids[pi]]["sess"]) for pi in keep}
print(f"  {len(_ses_ok)}/{len(USE_KEYS)} sessions survive; events per neuron "
      f"{min(payload[pi][0].shape[0] for pi in keep) if keep else 0}.."
      f"{max(payload[pi][0].shape[0] for pi in keep) if keep else 0}")
if not keep:
    raise ValueError(f"no neuron has >= MIN_TRIALS ({MIN_TRIALS}) events of codes "
                     f"{PP_KEEP}. A single code is ~half a session's trials, so "
                     f"lower MIN_TRIALS or pool the codes.")

# ── cache: this is the expensive step, so never repeat it for free ───────────
CACHE_KEY = dict(features=list(FEATURES_PP), paw=PAW, PP_KEEP=sorted(PP_KEEP),
                 KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP,
                 KIN_RANGE=(None if KIN_RANGE is None else list(KIN_RANGE)),
                 n_kin=n_kin, n_wins=n_wins,
                 N_TEST=N_TEST, N_SHUF_2D=(N_SHUF_2D if DO_NULL else 0),
                 WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP, WIN_RANGE=list(WIN_RANGE),
                 K_KSG=K_KSG, NEURAL_VAR=NEURAL_VAR, MIN_TRIALS=MIN_TRIALS,
                 n_neurons=len(run_idx), align="pushpull", zero=ALIGN, clean=CLEAN_ON,
                 exclude_toblerone=EXCLUDE_TOBLERONE, v="v1")
_CSTEM = (f"{RUN_NAME}_{'-'.join(FEATURES_PP)}_{NW_TAG}_{PP_TAG}_kw{KIN_WIN}_ks{KIN_STEP}"
          f"_n{'all' if N_TEST is None else N_TEST}_s{N_SHUF_2D}{MICE_TAG}"
          f"{ALIGN_TAG}{CLEAN_TAG}")
CACHE_PATH = config.var_path(f"{_CSTEM}_2Dcache.pkl")

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
                                 keep=keep, uids=list(uids),
                                 regions=regions_run), f, protocol=4)
            print(f"cache -> {CACHE_PATH} "
                  f"({CACHE_PATH.stat().st_size/1e6:.0f} MB)")

EXC2 = MI2 - NU2 if DO_NULL else MI2


# %% ── 4. mean 2D landscape per region ───────────────────────────────────────
# 4. cell --> average over all neurons in each region, and the CB - CP difference.
ASPECT = "equal" if EQUAL_AXES else "auto"
fdir   = config.fig_dir(OUT_DIR_C)
extent = [win_centers[0], win_centers[-1], kin_centers[0], kin_centers[-1]]
lab    = (f"KSG k={K_KSG}, {NEURAL_VAR}, joint [{', '.join(FEATURES_PP)}], "
          f"{PAW} paw, {'+'.join(PP_CODE[c] for c in sorted(PP_KEEP))}"
          f"{', no Toblerone' if EXCLUDE_TOBLERONE else ''}"
          f", t=0 at {ALIGN}{', cleaned' if CLEAN_ON else ''}"
          f"{' — EXCESS' if DO_NULL else ''}")


def _draw2d(ax, M, title, vmin=None, vmax=None):
    im = ax.imshow(M, origin="lower", aspect=ASPECT, extent=extent, cmap=CMAP,
                   vmin=vmin, vmax=vmax)
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(d, d, color="w", lw=1.0, ls="--", alpha=0.7)     # simultaneity
    ax.axvline(0, color="w", lw=0.8, ls=":", alpha=0.7)
    ax.axhline(0, color="w", lw=0.8, ls=":", alpha=0.7)
    ax.set(xlabel=f"NEURAL window centre re {ZERO_TXT} (s)",
           ylabel=f"KINEMATIC window centre re {ZERO_TXT} (s)", title=title)
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
if SAVE_FIGS:
    _savefig_retry(fig, fdir / f"{RUN_NAME}_mean2D.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

if all(MEAN2[r] is not None for r in ("CP", "CB")):
    D = MEAN2["CB"] - MEAN2["CP"]
    lim = np.nanpercentile(np.abs(D), 99)
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(D, origin="lower", aspect=ASPECT, extent=extent, cmap="RdBu_r",
                   norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim))
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(d, d, color="k", lw=1.0, ls="--", alpha=0.6)
    ax.axvline(0, color="k", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set(xlabel=f"NEURAL window centre re {ZERO_TXT} (s)",
           ylabel=f"KINEMATIC window centre re {ZERO_TXT} (s)",
           title=f"CB - CP  ({lab})")
    plt.colorbar(im, ax=ax, label="Δ mean MI (bits)")
    plt.tight_layout()
    if SAVE_FIGS:
        _savefig_retry(fig, fdir / f"{RUN_NAME}_diff2D.png", dpi=150,
                       bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 4b-7. the analysis: setup ─────────────────────────────────────────────
# 4b-7. cell --> the analysis and figures of the maps come from mi2d_figures.py,
# one function per step -- the same functions as for the reach alignment, so
# the two are analysed by the same code. `ctx` holds this run's inputs; the
# labels name this alignment's zero (ZERO_TXT) instead of the reach endpoint.
from types import SimpleNamespace
import importlib, mi2d_figures
importlib.reload(mi2d_figures)          # pick up edits to the module

ctx = SimpleNamespace(
    EXC2=EXC2, keep=keep, uids=uids, run_idx=run_idx, regions_run=regions_run,
    win_centers=win_centers, kin_centers=kin_centers, n_kin=n_kin, n_wins=n_wins,
    WIN_STEP=WIN_STEP, KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP, CLR=CLR, lab=lab,
    fdir=fdir, RUN_NAME=RUN_NAME, SAVE_FIGS=SAVE_FIGS, SHOW_FIGS=SHOW_FIGS,
    extent=extent, ASPECT=ASPECT, CMAP=CMAP, zero_txt=ZERO_TXT, TIMES=TIMES)


# %% ── 4b. h(lag) of the mean landscape, per mouse ──────────────────────────
# z-scored h(lag) per mouse, robust peak and median lag with paired CB-vs-CP
# insets, the full two-row figure and the presentation figure (PNG + PDF).
HLAG = mi2d_figures.hlag_per_mouse(ctx, LAG_BIN=LAG_BIN, PER_MOUSE_MIN_N=PER_MOUSE_MIN_N)


# %% ── 5. DIAGONAL vs NEURAL TIME ───────────────────────────────────────────
# per neuron: does the map vary along the diagonal (lag) or along neural time?
DVN = mi2d_figures.diagonal_vs_neural(ctx, LAG_BIN=LAG_BIN, CV_GAP=CV_GAP)


# %% ── 5b. individual 2D landscapes ─────────────────────────────────────────
# Slow (one figure per example neuron, N_EXAMPLES per region): skip this cell if
# you do not need them -- nothing below depends on it.
mi2d_figures.individual_landscapes(ctx, DVN, OUT_DIR_C, N_EXAMPLES=N_EXAMPLES,
                                   EXAMPLE_MODE=EXAMPLE_MODE)


# %% ── 6. DIAGONAL vs NEURAL TIME: the figure ───────────────────────────────
mi2d_figures.diagonal_vs_neural_figure(ctx, DVN, FIG_SUBSET=FIG_SUBSET)


# %% ── 7. per-mouse mean landscapes, timing, results ───────────────────────
mi2d_figures.per_mouse_maps(ctx, PER_MOUSE_MIN_N=PER_MOUSE_MIN_N,
                            PER_MOUSE_FILES=PER_MOUSE_FILES,
                            save_pdf_all_mice=SAVE_PDF_ALL_MICE)
mi2d_figures.timing_summary(ctx, N_SHUF_2D=N_SHUF_2D, n_population=len(INIT))
if SAVE_VARS:
    mi2d_figures.save_results(
        config.var_path(f"{RUN_NAME}_{'-'.join(FEATURES_PP)}_{NW_TAG}_{PP_TAG}"
                        f"{MICE_TAG}{ALIGN_TAG}{CLEAN_TAG}_2D.pkl"),
        ctx, DVN, MI2, NU2, NSD2, DO_NULL,
        PP_KEEP=sorted(PP_KEEP), zero=ALIGN, clean=CLEAN_ON,
        features=list(FEATURES_PP), KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP,
        KIN_RANGE=KIN_RANGE, WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP,
        WIN_RANGE=WIN_RANGE, K=K_KSG, N_TEST=N_TEST, N_SHUF_2D=N_SHUF_2D,
        NEURAL_VAR=NEURAL_VAR)

# %%
