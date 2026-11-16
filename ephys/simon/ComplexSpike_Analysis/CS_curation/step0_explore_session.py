# -*- coding: utf-8 -*-
"""
step0_explore_session.py — CS / SS candidate exploration for ONE session.

Exploration only: enumerate and visualise CS units, their candidate SS partners
(CS→SS pause) and CS–CS duplicate candidates, with evidence scores. No decisions.

Spike times are in seconds on the imec clock (spike_times / sample_rate from params.py).
All unit and pair metrics use ONE analysis window, the same for every unit (ANALYSIS_WINDOW):
  "bhv_win_selected": the lab's selected good-behaviour + stable-ephys window, bhv.win_selected in
      behavior_fundamentals.mat (written by behavior/paw_and_body_trajectories/step3_bhv_struct.m:
      start = trial-init time of the first selected trial, end = behavior_duration end; Harp clock).
      It is mapped to the imec clock with the linear ephys→Harp fit recovered from the paired
      st_raw_all / st_all fields of neurons_session.mat (s01_NeuronStructAndPlot_v5.m).
  "full": the whole recording.
Overview figures still show the whole recording, with the analysis window shaded.

GOLDEN RULE: the Dropbox tree is READ-ONLY. This script never writes there
(see safe_out_dir) and opens every data file read-only.
"""

# %% 1. PARAMETERS + PATHS
from pathlib import Path
import json
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.widgets import Button, RadioButtons, Slider
from scipy.stats import poisson, mannwhitneyu
from scipy.signal import find_peaks
import scipy.io as sio
import h5py

# --- roots (local copy preferred, Dropbox original as read-only fallback)
LOCAL_ROOT = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior")
DROPBOX_ROOT = Path(r"D:\Learning Lab Dropbox")
FALLBACK_ROOT = DROPBOX_ROOT / r"Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_data\20230801_ChocolateGroup"
CODE_DIR = Path(r"C:\Users\simon\Desktop\Code\EPhys\ComplexSpike_Analysis")

# --- session
MOUSE = "1_CoteDor"
SESSION_R = "R5"      # behaviour name
SESSION_S = "S5"      # ephys name
IMEC = "imec0"        # CB probe (right-pawed mouse)
SORTER_REL = (r"1_CoteDor\04082023_CoteDor_Cer_S5_g0\04082023_CoteDor_Cer_S5_g0_imec0"
              r"\ibl_sorter_results_driftAdapt40")
SORTER_TYPE = "ibl_sorter"
RUN_TAG = "v8"                           # v2: tier-2 SS + duplicate evidence; v3: per-unit stable spans + template
                                         #   overlays; v4: unified script (CS-spike sharing flag, <2-CS / no-SS
                                         #   guards); v5: tier 2 without label gate (SS2_GROUPS = None);
                                         #   v6 (Simon 2026-11-07): pause measured AFTER the CS event (4–12 ms),
                                         #   tier-1 FR floor, CS event decomposition table + event-sharing flag;
                                         #   v7 (Simon 2026-11-07): pause window 3–8 ms (4–12 ms averaged in the
                                         #   post-pause rebound and dropped the clearest pairs).
                                         #   v8 (Simon 2026-11-16): ADDS literature metrics (latency test, SS|CS vs
                                         #   SS|SS, rate-scaled duration, candidate_lit), Lussac paper/code rules,
                                         #   spatial check, CS unit checks, manual-only and screen CS. Every v7
                                         #   column, `candidate`, `rank` and table keeps its meaning (see section 1b).
                                         #   Earlier runs stay untouched (each tag has its own folder)

# --- overrides (optional, environment variables; same output naming for single and batch runs):
#     STEP0_RUN_TAG="<tag>"           → output folder <mouse>_<R#>/step0_explore_bhvspan_<tag> (default RUN_TAG above)
#     STEP0_SESSION="<mouse>:<S#>"    → session taken from the sorter-comparison session table (one C4 sorter folder
#                                       per session). Without it, the hard-coded session above is used.
#     STEP0_SS1_FLOOR="20" | "10"     → tier-1 FR floor toggle (see SS1_FLOOR below)
import os as _os
RUN_TAG = _os.environ.get("STEP0_RUN_TAG", RUN_TAG)

# --- TOGGLE: tier-1 (extra-good) SS FR floor, Hz. "20" = same floor as tier 2 (Simon 2026-11-06, default);
#     "10" = keeps 10–20 Hz extra-good units such as CoteDor R5 938 (14 Hz; pauses 0.06–0.07 after CS 585/734 in v7)
#     and Ferrero R6 344 (13 Hz). Switch here or with env STEP0_SS1_FLOOR. A non-default floor writes to its own
#     folder (suffix "_fr10"), so the two settings never overwrite each other.
SS1_FLOOR = "20"
SS1_FLOOR = _os.environ.get("STEP0_SS1_FLOOR", SS1_FLOOR)
_SS1_FLOOR_HZ = {"20": 20.0, "10": 10.0}
assert SS1_FLOOR in _SS1_FLOOR_HZ, f"STEP0_SS1_FLOOR must be one of {list(_SS1_FLOOR_HZ)}, got {SS1_FLOOR!r}"
_FLOOR_SUFFIX = "" if SS1_FLOOR == "20" else f"_fr{SS1_FLOOR}"
if _os.environ.get("STEP0_SESSION"):
    _m, _s = _os.environ["STEP0_SESSION"].split(":")
    _tab = pd.read_csv(LOCAL_ROOT / "CS_curation" / "_sorter_comparison" / "session_table.csv")
    _row = _tab[(_tab["mouse"] == _m) & (_tab["S"] == _s)].iloc[0]
    MOUSE, SESSION_S, SESSION_R, IMEC = _row["mouse"], _row["S"], _row["R"], _row["imec"]
    SORTER_REL = str(Path(_row["sorter_dir"]).relative_to(FALLBACK_ROOT))
    SORTER_TYPE = {"ibl": "ibl_sorter", "ks4": "kilosort4"}[_row["sorter"]]
    print(f"STEP0_SESSION override: {MOUSE} {SESSION_S}/{SESSION_R} {IMEC} {SORTER_TYPE} | {SORTER_REL}")

# Candidate sorter-folder locations, tried in order. Add layouts of the local copy here.
SORTER_CANDIDATES = [
    LOCAL_ROOT / SORTER_REL.lstrip("\\/"),
    LOCAL_ROOT / "ephys_data" / "20230801_ChocolateGroup" / SORTER_REL.lstrip("\\/"),
    FALLBACK_ROOT / SORTER_REL.lstrip("\\/"),
]

# --- analysis window: "bhv_win_selected" (lab's stable window) or "full"
ANALYSIS_WINDOW = "bhv_win_selected"
EDGE_PAD = 0.1                           # s; CCG targets are taken from the window ± this pad
BHV_REL = Path("20230801_ChocolateGroup") / MOUSE / SESSION_R / "behavior_fundamentals.mat"
BHV_CANDIDATES = [
    LOCAL_ROOT / "mat_files" / BHV_REL,
    LOCAL_ROOT / BHV_REL,
    DROPBOX_ROOT / r"Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_and_behavior\mat_files" / BHV_REL,
]
# neurons_session.mat (s01 output) holds the paired ephys/Harp spike times used for the clock map
NEURONS_SESSION_CANDIDATES = [c / "neurons_session.mat" for c in SORTER_CANDIDATES]

# --- per-unit stable spans (bhv_win_selected only). win_selected is the REFERENCE: it is always kept, and each
#     CS/SS unit's window is extended outward, inside the behaviour period (first to last behavioural event,
#     bhv.all_trials), bin by bin while its firing rate and amplitude stay close to its own reference values.
#     A CS×SS (or CS–CS) pair uses the overlap of its two units' spans. False → win_selected for every unit.
STABLE_SPAN = True
STAB_BIN = 60.0                          # s; bins step outward from the reference window edges
STAB_FR_TOL = 1.5                        # bin FR within ×/÷ this of the unit's reference FR (CS ~1 Hz: ~60 spikes/bin)
STAB_AMP_TOL = 1.3                       # bin median amplitude within ×/÷ this of the reference median
STAB_MAX_GAP = 1                         # isolated unstable bins tolerated (only if a stable bin follows)
STAB_MIN_BIN = 30.0                      # s; shorter partial bins at the behaviour edges are not evaluated

_OUT_TAG = {"full": "step0_explore",
            "bhv_win_selected": "step0_explore_bhvspan" if STABLE_SPAN else "step0_explore_bhvwin"}[ANALYSIS_WINDOW]
OUT_DIR = (LOCAL_ROOT / "CS_curation" / f"{MOUSE}_{SESSION_R}"
           / ((f"{_OUT_TAG}_{RUN_TAG}" if RUN_TAG else _OUT_TAG) + _FLOOR_SUFFIX))

# --- population
REGION = "CB"
CS_LABEL = "PkC_cs"
# SS partner candidates: tier 1 = extra-good non-CS units; tier 2 = good/mua units that are not extra-good
SS2_GROUPS = None                        # None = NO label gate (v5, Simon 2026-11-05): any cluster_group label,
                                         #   incl. noise/unsorted. cluster_group means different things per sorter
                                         #   (KS4: manual Phy edits on top of KSLabel; IBL: automatic), so gating on
                                         #   it biased the KS4 vs IBL comparison. v2-v4 used ("good", "mua").
SS2_GROUP_TXT = "any label" if SS2_GROUPS is None else "/".join(SS2_GROUPS)
SS2_FR_MIN = 20.0                        # Hz in the analysis window (tier 2; SS fire tens of Hz)
SS1_FR_MIN = _SS1_FLOOR_HZ[SS1_FLOOR]    # toggle above (default 20 = SS2_FR_MIN). Extra-good units below
                                         #   it get role "SS_lowfr": never SS candidates, but their pause scores are
                                         #   kept for audit in tables/cs_ss_pairs_lowfr_excluded.csv

# --- correlograms (seconds)
COARSE_BIN, COARSE_WIN = 1e-3, 50e-3     # pauses
FINE_BIN, FINE_WIN = 0.1e-3, 10e-3       # duplicates, collision dips, zero-lag peaks
ACG_FINE_BIN, ACG_FINE_WIN = 0.25e-3, 20e-3
ACG_COARSE_BIN, ACG_COARSE_WIN = 1e-3, 150e-3

# --- CS→SS pause metrics (seconds, windows are [start, stop))
# v6: the pause is measured AFTER the CS event. 0–4 ms still holds the CS event itself (main spike + spikelets);
#   a unit that absorbed CS pieces has extra spikes there and looked like a weak pauser in v1–v5 (R5 raw layer:
#   CS 585 → 938 pauses 0.30 at 4–12 ms but not at 0–8 ms). The v5 window is kept as audit columns (*_0_8).
# v7 (Simon 2026-11-07): 3–8 ms. v6 used 4–12 ms, which reached into the post-pause rebound (1.3–1.8× baseline at
#   8–10 ms) and lost the clearest pairs (Toblerone R1 379→349 0.07→0.71, R5 871→841 0.12→0.76). 3–8 ms skips most of
#   the CS event and ends before the rebound (what-if: _sorter_comparison/step0_v6_vs_v5/window_whatif/).
PAUSE_WIN = (3e-3, 8e-3)
PAUSE_WIN_EXPLORE = (3e-3, 20e-3)
MIRROR_WIN = (-8e-3, -3e-3)
PAUSE_WIN_V5 = (0.0, 8e-3)               # audit only: pause_ratio_0_8 / mirror_ratio_0_8 / asym_0_8
MIRROR_WIN_V5 = (-8e-3, 0.0)
BASELINE_WIN = (-50e-3, -10e-3)          # coarse-CCG baseline (same CCG it normalises)
FINE_BASELINE_WIN = (-10e-3, -2e-3)      # fine-CCG baseline (same fine CCG)
ZERO_WIN = (-0.5e-3, 0.5e-3)
DIP_HALF = 0.5e-3                        # collision dip: [-DIP_HALF,0) and [0,DIP_HALF) in fine CCG
RECOVERY_FRACS = (0.5, 0.8)
LUSSAC_FRAC = 0.4                        # v7 `lussac_pass`: Lussac-STYLE medians on the v7 3–8 ms window (not Lussac)

# %% 1b. v8 LITERATURE BLOCK (Simon 2026-11-16). ADDS columns/flags; nothing here changes `candidate`, `rank`,
#     duplicate evidence, event decomposition or stable spans. All thresholds are PROVISIONAL.
# --- extra CS sources (both False → the v7 CS population, i.e. C4 PkC_cs only)
INCLUDE_MANUAL_CS = True                 # units flagged `y` in cluster_CS.tsv / cluster_potential_CS.tsv (any C4 label
                                         #   or none, extra-good or not, CB) enter as CS. Changes the 2026-11-03 decision
                                         #   "manual flag = diagnostic only" (Simon 2026-11-16).
INCLUDE_SCREEN_CS = True                 # C4-UNLABELLED extra-good CB units that pass the CS screen enter as CS
SCREEN_NON_EG = False                    # also screen non-extra-good CB units (no C4 label at all); later, default off
if _os.environ.get("STEP0_V8_EXTRA_CS") == "0":   # v7-reproduction check: C4 CS only (both sources off)
    INCLUDE_MANUAL_CS = INCLUDE_SCREEN_CS = False
    print("STEP0_V8_EXTRA_CS=0: manual and screen CS sources OFF (v7 CS population)")
#   Option A (Simon 2026-11-16): a manual / screen unit that is a v7 SS candidate (role SS_cand / SS_cand2) STAYS an SS
#   candidate and is NOT added as CS, so every v7 pair stays identical → tables/manual_kept_as_ss.csv.
# --- CS unit checks / screen (flags; the screen ENTRY rule is cs_fr_ok AND isi_5_35_ok)
CS_FR_RANGE = (0.2, 5.0)                 # Hz. Lussac paper v1 §4.2 (Llobet, Wyngaard & Barbour 2022, bioRxiv
                                         #   10.1101/2022.02.08.479192v1); Shin et al. 2007 EJN: 0.2–4.7 Hz in vivo
CS_ISI_5_35 = (5e-3, 35e-3)              # Lussac paper v1 §4.2: "< 5% of ISIs between 5 and 35 ms" (share of ISIs)
CS_ISI_5_35_MAX = 0.05
CS_ISI_10_35 = (10e-3, 35e-3)            # Lussac CODE v1 (tag v1, commit af2b9ef, params_cerebellar_cortex.json:
                                         #   "ISI range [10, 35], max 0.05", pipeline.categorize_units: ISIs in range ÷
                                         #   number of spikes, whole recording). Column only.
CS_MIN_IOI = 50e-3                       # minimum inter-CS interval in awake mice ~50 ms (Negrello et al. 2019 PLoS
                                         #   Comput Biol, "Quasiperiodic rhythms of the inferior olive") → isi_lt50_frac
CS_CV_REF, CS_CV2_REF = 0.74, 0.84       # awake mouse CS, near-Poisson (iScience 2022, Table 1); reference values only
PRESENCE_BIN = 60.0                      # s; confound flag partial_presence: span < PRESENCE_MIN_SPAN × reference window
PRESENCE_MIN_SPAN, PRESENCE_MIN_BINS = 0.5, 0.8   #   or < this share of PRESENCE_BIN bins in the span hold a spike
ISO_ISI1_MAX = 0.01                      # confound flag poor_isolation: frac_isi_lt1ms > this, or sorter score below
ISO_CONTAM_MAX = 20.0                    #   KS4 cluster_ContamPct.tsv (%)
ISO_NOISE_PRED_MAX = 0.5                 #   IBL cluster_predictions.tsv noise_predict
TMPL_LATE_MIN = 0.2                      # template late peak: |value| > this × |trough|, >= 0.5 ms after the trough
# --- rate-relative pause (Zhou et al. 2014 eLife, 2015 J Neurophysiol; Shin & De Schutter 2006 J Neurophysiol;
#     P-sort, Sedaghat-Nejad et al. 2022 J Neurophysiol). Same pair window and spike sets as the v7 metrics.
LAT_DEAD = PAUSE_WIN[0]                  # dead zone on BOTH sides of the CS (skips CS pieces, sorter blind spot and the
                                         #   0–2 ms ephaptic dip of neighbouring PCs, Han et al. 2020 Nat Neurosci)
REFR_WIN = 25e-3                         # ± lag range of the SS autocorrelogram used for SS|SS (1 ms bins, = COARSE_BIN)
LIT_LAT_P_MAX = 0.05                     # latency test: one-sided Mann-Whitney p (L_post > L_pre)
LIT_CLIFF_MIN = 0.1                      # and Cliff's delta (v7 data: candidates median 0.23, pairs > 500 µm ~0)
LIT_REFR_RATIO_MAX = 0.7                 # SS|CS ÷ SS|SS in PAUSE_WIN ("pause beyond refractoriness")
LIT_DUR_ISI_K = 1.0                      # dur_rate_ok: pause_dur50_ms > K × SS median ISI. COLUMN ONLY (Rule 1, Simon
                                         #   2026-11-16: v7 candidates 31% pass at K=1, pairs > 500 µm 1%)
#   candidate_lit = lat_p < LIT_LAT_P_MAX & lat_cliff >= LIT_CLIFF_MIN & refr_ratio < LIT_REFR_RATIO_MAX & ~collision_dip
# --- Lussac reference rules (columns only; `lussac_pass` above is v7's own 3–8 ms variant and is unchanged)
LUSSAC_P1_WIN, LUSSAC_P1_BASE, LUSSAC_P1_FRAC = (0.0, 8e-3), (-20e-3, 0.0), 0.5   # paper v1 §5: SS rate 0–8 ms after
                                         #   the CS < 50% of the rate in the 20 ms before → lussac_paper_v1_pass
LUSSAC_CM_HALF, LUSSAC_CM_WIN, LUSSAC_CM_FRAC = 12e-3, (0.0, 8e-3), 0.4   # code MAIN (Lussac 2, commit 2a0dce1,
                                         #   modules/find_purkinje_cells.py): SpikeInterface correlograms window_ms=25 =
                                         #   FULL width → ±12 ms (rounded down to whole 1 ms bins), raw counts, SS − CS;
                                         #   pause: median[0, 8) < 0.4 × median(12 negative-lag bins); REJECT if
                                         #   np.median(mirror bins < thr) is truthy (median of booleans, i.e. >= 4 of the
                                         #   8 bins at −8…0 ms below thr; reproduced as is) → lussac_code_main_pass
LUSSAC_CM_SS_MIN, LUSSAC_CM_CS_RANGE = 40.0, (0.5, 3.0)   # code main FR gates → lussac_code_main_fr_ok (column only)
# --- spatial check (flags only; Calame et al. 2023 Nat Neurosci: partner on adjacent channels, CS may be recorded on
#     the dendrite ~100–200 µm from the soma; Han et al. 2020: ephaptic suppression <= 50 µm)
FAR_UM = 200.0                           # far_flag
FP_UM = 500.0                            # pairs farther apart cannot be one cell → empirical false-positive rate
STRONG_PAUSE = 0.3                       # "strongest pairs" for the depth-offset sign: candidate & pause_ratio < this

# --- null: jitter of the CS train
JITTER_WIN = PAUSE_WIN                   # v6+: same window as the pause ratio (v1–v5: 1–8 ms)
JITTER_HALF = 25e-3
N_SURR = 500
RNG_SEED = 0

# --- stability
BLOCK_LEN = 600.0                        # s

# --- CS–CS duplicates
COINC_DELTA = 0.3e-3                     # ± window around the PEAK lag
CS_MIN_ISI = 10e-3                       # union ISIs shorter than this are "too short for a CS"
DUP_BASELINE_ABS = 5e-3                  # fine-CCG bins with |lag| > this are the duplicate baseline

# --- permissive candidate thresholds (exploration only, not decisions)
CAND_PAUSE_RATIO_MAX = 0.7
CAND_JITTER_P_MAX = 0.05
CAND_ASYM_MIN = 0.2                      # pause must be one-sided (a symmetric dip also gives small jitter p)
DUP_COINC_MIN = 0.2                      # duplicate candidate flag uses the coincident fraction only
DUP_PEAK_P_MAX = 1e-3                    # Poisson peak test: reported as a column, not used for the flag
CS_SHARE_FLAG = 0.2                      # CS→SS pair: share of CS spikes with an SS spike in ZERO_WIN, minus the
                                         #   Poisson chance for that SS rate, above this → cs_share_flag. FLAG ONLY:
                                         #   it does not change `candidate`; the "SS" may carry CS spikelets or be a
                                         #   split piece of the CS (2026-11-05 sorter comparison, zero-lag check)

# --- CS EVENT DECOMPOSITION (v6, section 5b; from spike times only, nothing is merged or excluded)
EVENT_WIN = (-0.5e-3, 4e-3)              # s after a CS spike: a spike of another cluster here is part of the CS event?
EVENT_NEIGH_UM = 250.0                   # clusters within this distance of the CS peak channel are tested
EVENT_MIN_EXCESS = 0.05                  # listed if (CS events with a spike of the cluster) − Poisson chance >= this
EVENT_MIN_N = 30                         #   and at least this many CS events carry a spike of the cluster
EVENT_ZERO_LAG = 0.3e-3                  # |median lag| <= this → "zero-lag copy", else "earlier/late piece"
EVENT_LOCKED_MOST = 0.5                  # >= this share of the cluster's own spikes are CS-locked → "mostly CS pieces"
CS_EVENT_FLAG = 0.2                      # CS→SS pair flag: SS spikes in EVENT_WIN on > this share of CS events
                                         #   (beyond chance). FLAG ONLY, candidate unchanged (pause is at 3–8 ms)

# --- duplicate EVIDENCE (section 6b). Nothing is merged: each line votes +1 (same cell), -1 (distinct
#     cells) or 0 (uninformative / not testable). Thresholds are provisional, tuned on nothing yet.
LAG_SPREAD_WIN = 1e-3                    # ± window around the peak lag used to measure lag jitter
DETREND_BIN = 60.0                       # s; amplitude/depth medians removed per bin before correlating
UNION_ISI_RANGE = (4e-3, 25e-3)          # union ISIs too short for two CS of one cell, too long for spikelets
DOUBLET_ISI = 5e-3                       # within-unit ISIs below this = multiple detections of one CS
MIN_SPIKES_COND = 50                     # minimum spikes for a conditional (subset-triggered) pause ratio
EV_COINC = (0.5, 0.2)                    # coinc_frac_peak: >= first → +1, < second → -1
EV_LAG_SD = (0.15e-3, 0.5e-3)            # lag jitter (robust SD): <= first → +1, >= second → -1
EV_UNION_NEW = (0.2, 0.5)                # new union ISIs in UNION_ISI_RANGE / number expected for independent
                                         #   trains: <= first → +1 (one cell), >= second → -1 (two cells)
UNION_MIN_EXPECTED = 5                   # fewer expected new ISIs → union line not testable
EV_COND_PR = (0.7, 0.9)                  # other-unit-only spikes on this unit's SS partner: <= → +1, >= → -1
EV_FOOTPRINT = (0.8, 0.3)                # template ptp-footprint cosine (only when footprints overlap)
EV_COVAR = 0.3                           # |detrended amp corr| or depth corr >= this → +1 (never -1)
FOOTPRINT_FRAC = 0.1                     # channel belongs to a footprint if ptp > this × peak ptp
FOOTPRINT_MIN_OVERLAP = 6                # fewer shared footprint channels → footprint "not testable"
EV_CONVERGE_MIN = 3                      # status "converges_*": >= this many lines one way and none the other way
# lag_sd and covar describe only the coincident spikes: they vote only if coinc_frac_peak >= EV_COINC[1]

# --- display
OPEN_QT_BROWSER = False                  # matplotlib/Qt browser (section 8); the HTML page (8b) is always written
TOP_N = 8
N_NEIGH_CH = 4
CB_COLOR = np.array([154, 123, 193]) / 255
BAR_COLOR = (0.2, 0.2, 0.2)
FR_FILL = (0.8, 0.8, 0.9)
DOT_COLOR = (0.2, 0.2, 0.2, 0.1)
SHADE = (0.88, 0.88, 0.88)


def safe_out_dir(p: Path) -> Path:
    """Create an output folder, refusing anything on Dropbox or inside the code folder."""
    p = Path(p).resolve()
    for forbidden in (DROPBOX_ROOT, CODE_DIR):
        try:
            p.relative_to(forbidden.resolve())
        except ValueError:
            continue
        raise RuntimeError(f"Refusing to write under {forbidden}: {p}")
    p.mkdir(parents=True, exist_ok=True)
    return p


def find_sorter_dir(candidates):
    for c in candidates:
        if (c / "spike_times.npy").exists():
            print(f"sorter folder: {c}")
            return c
    tried = "\n  ".join(str(c) for c in candidates)
    raise FileNotFoundError(f"No sorter folder with spike_times.npy found. Tried:\n  {tried}")


SORTER_DIR = find_sorter_dir(SORTER_CANDIDATES)
OUT_DIR = safe_out_dir(OUT_DIR)
FIG_DIR = safe_out_dir(OUT_DIR / "figs")
TAB_DIR = safe_out_dir(OUT_DIR / "tables")
print(f"output folder: {OUT_DIR}")


# %% 2. FUNCTIONS (+ synthetic self-test)

def set_lab_style():
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "font.size": 10,
    })


def _pair_diffs(ref, tgt, win, exclude_self=False, chunk=200_000):
    """All tgt - ref differences in [-win, win), without a loop over spikes (chunked over ref)."""
    out = []
    for s in range(0, len(ref), chunk):
        r = ref[s:s + chunk]
        lo = np.searchsorted(tgt, r - win, "left")
        hi = np.searchsorted(tgt, r + win, "left")
        n = hi - lo
        tot = int(n.sum())
        if tot == 0:
            continue
        ri = np.repeat(np.arange(len(r)), n)
        j = lo[ri] + (np.arange(tot) - np.repeat(np.cumsum(n) - n, n))
        if exclude_self:
            keep = j != (ri + s)
            ri, j = ri[keep], j[keep]
        out.append(tgt[j] - r[ri])
    return np.concatenate(out) if out else np.empty(0)


def ccg(ref, tgt, bin_s, win_s, acg=False):
    """Cross-correlogram of tgt relative to ref, as tgt rate (sp/s).

    Bins are [k*bin, (k+1)*bin) for k = -nb .. nb-1; returns (rate, left_edges, counts).
    """
    nb = int(round(win_s / bin_s))
    d = _pair_diffs(ref, tgt, nb * bin_s, exclude_self=acg)
    idx = np.floor(d / bin_s + 1e-9).astype(np.int64) + nb
    idx = idx[(idx >= 0) & (idx < 2 * nb)]
    counts = np.bincount(idx, minlength=2 * nb).astype(float)
    left = (np.arange(2 * nb) - nb) * bin_s
    rate = counts / (max(len(ref), 1) * bin_s)
    return rate, left, counts


def _win_mask(left, bin_s, win):
    """Bins fully inside [win0, win1)."""
    return (left >= win[0] - 1e-12) & (left + bin_s <= win[1] + 1e-12)


def pause_metrics(rate_c, left_c, rate_f, left_f):
    """CS→SS pause scores from ONE coarse and ONE fine CCG of the same pair."""
    m = {}
    base_mask = _win_mask(left_c, COARSE_BIN, BASELINE_WIN)
    base = rate_c[base_mask].mean()
    base_med = np.median(rate_c[base_mask])
    pw = _win_mask(left_c, COARSE_BIN, PAUSE_WIN)
    mw = _win_mask(left_c, COARSE_BIN, MIRROR_WIN)
    m["baseline_rate"] = base
    m["pause_ratio"] = rate_c[pw].mean() / base if base > 0 else np.nan
    m["mirror_ratio"] = rate_c[mw].mean() / base if base > 0 else np.nan
    m["asym"] = m["mirror_ratio"] - m["pause_ratio"]
    p8, m8 = _win_mask(left_c, COARSE_BIN, PAUSE_WIN_V5), _win_mask(left_c, COARSE_BIN, MIRROR_WIN_V5)
    m["pause_ratio_0_8"] = rate_c[p8].mean() / base if base > 0 else np.nan
    m["mirror_ratio_0_8"] = rate_c[m8].mean() / base if base > 0 else np.nan
    m["asym_0_8"] = m["mirror_ratio_0_8"] - m["pause_ratio_0_8"]
    pe = _win_mask(left_c, COARSE_BIN, PAUSE_WIN_EXPLORE)
    m["pause_ratio_explore"] = rate_c[pe].mean() / base if base > 0 else np.nan
    # recovery: first bin after the pause minimum (within PAUSE_WIN) whose rate is back above frac * baseline.
    # Searching from the minimum, not from lag 0, so a 1–2 ms CS→SS latency does not read as a 1 ms pause.
    pos = left_c >= -1e-12
    r_pos, l_pos = rate_c[pos], left_c[pos]
    i_min = int(np.argmin(np.where(_win_mask(l_pos, COARSE_BIN, (0.0, PAUSE_WIN[1])), r_pos, np.inf)))
    m["pause_min_lag_ms"] = l_pos[i_min] * 1e3
    for f in RECOVERY_FRACS:
        ok = np.nonzero(r_pos[i_min:] >= f * base)[0]
        m[f"pause_dur{int(f * 100)}_ms"] = l_pos[i_min + ok[0]] * 1e3 if len(ok) else np.nan
    # Lussac-style reference criterion (medians, one-sided)
    pos_med = np.median(rate_c[pw])
    neg_med = np.median(rate_c[mw])
    m["lussac_pass"] = bool(pos_med < LUSSAC_FRAC * base_med and not neg_med < LUSSAC_FRAC * base_med)
    # fine-CCG scores: baseline from the same fine CCG
    fb = rate_f[_win_mask(left_f, FINE_BIN, FINE_BASELINE_WIN)].mean()
    m["fine_baseline_rate"] = fb
    m["zero_peak_ratio"] = rate_f[_win_mask(left_f, FINE_BIN, ZERO_WIN)].mean() / fb if fb > 0 else np.nan
    dn = rate_f[_win_mask(left_f, FINE_BIN, (-DIP_HALF, 0.0))].mean()
    dp = rate_f[_win_mask(left_f, FINE_BIN, (0.0, DIP_HALF))].mean()
    late = rate_c[_win_mask(left_c, COARSE_BIN, (1e-3, PAUSE_WIN[1]))].mean()
    m["collision_dip"] = bool(fb > 0 and dn < 0.5 * fb and dp < 0.5 * fb and late >= 0.8 * base)
    return m


def count_in_window(ref, tgt, win):
    """Number of (ref, tgt) pairs with tgt - ref in [win0, win1)."""
    return int((np.searchsorted(tgt, ref + win[1], "left") - np.searchsorted(tgt, ref + win[0], "left")).sum())


def jitter_pause_test(cs, ss, win=JITTER_WIN, half=JITTER_HALF, n_surr=N_SURR, seed=RNG_SEED):
    """Jitter null: shift each CS by U(-half, half); statistic = SS count in win after CS.

    The window starts at 1 ms so that a ~1 ms sorter collision dip cannot drive it.
    One-sided (pause = deficit): p = (1 + #{surr <= obs}) / (1 + n_surr), floor 1/(1+n_surr).
    jitter_deficit = 1 - obs / surrogate mean is the effect size (0 = no pause, 1 = complete pause).
    """
    rng = np.random.default_rng(seed)
    obs = count_in_window(cs, ss, win)
    jit = cs[None, :] + rng.uniform(-half, half, (n_surr, len(cs)))
    flat = jit.ravel()
    c = np.searchsorted(ss, flat + win[1], "left") - np.searchsorted(ss, flat + win[0], "left")
    surr = c.reshape(n_surr, len(cs)).sum(1)
    p = (1 + np.sum(surr <= obs)) / (1 + n_surr)
    z = (obs - surr.mean()) / surr.std() if surr.std() > 0 else np.nan
    deficit = 1 - obs / surr.mean() if surr.mean() > 0 else np.nan
    return {"jitter_obs": obs, "jitter_mean": surr.mean(), "jitter_p": p, "jitter_z": z, "jitter_deficit": deficit}


def cs_ss_pair(cs, ss, with_jitter=True):
    rc, lc, _ = ccg(cs, ss, COARSE_BIN, COARSE_WIN)
    rf, lf, _ = ccg(cs, ss, FINE_BIN, FINE_WIN)
    m = pause_metrics(rc, lc, rf, lf)
    if with_jitter:
        m.update(jitter_pause_test(cs, ss))
    return m, (rc, lc), (rf, lf)


def cs_spike_sharing(cs, ss, t_start, t_end):
    """Do the CS spikes reappear in this SS unit at zero lag?  (cs, ss sorted; ss may carry an edge pad)

    cs_share_frac   = fraction of CS spikes with >= 1 SS spike within ZERO_WIN
    cs_share_null   = Poisson chance of that for the SS rate in [t_start, t_end)
    cs_share_excess = cs_share_frac - cs_share_null   (share of CS spikes also in the SS unit)
    ss_from_cs_frac = fraction of the SS unit's spikes within ZERO_WIN of a CS spike
    """
    ss_w = ss[(ss >= t_start) & (ss < t_end)]
    if len(cs) == 0 or len(ss_w) == 0 or t_end <= t_start:
        return dict(cs_share_frac=np.nan, cs_share_null=np.nan, cs_share_excess=np.nan, ss_from_cs_frac=np.nan)
    hit = np.searchsorted(ss, cs + ZERO_WIN[1], "right") > np.searchsorted(ss, cs + ZERO_WIN[0], "left")
    back = np.searchsorted(cs, ss_w - ZERO_WIN[0], "right") > np.searchsorted(cs, ss_w - ZERO_WIN[1], "left")
    null = 1 - np.exp(-len(ss_w) / (t_end - t_start) * (ZERO_WIN[1] - ZERO_WIN[0]))
    return dict(cs_share_frac=float(hit.mean()), cs_share_null=float(null),
                cs_share_excess=float(hit.mean() - null), ss_from_cs_frac=float(back.mean()))


def cs_event_capture(cs, k, t_start, t_end, win=EVENT_WIN):
    """Does cluster k fire inside the CS event (win after each CS spike) more often than chance?  (v6)

    cs sorted, inside [t_start, t_end); k sorted, may carry an edge pad.
    frac   = share of CS events with >= 1 spike of k in win;  null = Poisson chance for k's rate in the window
    excess = frac - null;  lags = first k spike in win per hit event (s)
    k_locked_frac = share of k's own spikes (in the window) that fall in win after some CS spike
    """
    k_w = k[(k >= t_start) & (k < t_end)]
    if len(cs) == 0 or len(k_w) == 0 or t_end <= t_start:
        return dict(n_hit=0, frac=np.nan, null=np.nan, excess=np.nan, k_locked_frac=np.nan), np.empty(0)
    i0, i1 = np.searchsorted(k, cs + win[0], "left"), np.searchsorted(k, cs + win[1], "left")
    hit = i1 > i0
    null = 1 - np.exp(-len(k_w) / (t_end - t_start) * (win[1] - win[0]))
    lags = k[i0[hit]] - cs[hit]
    locked = np.searchsorted(cs, k_w - win[0], "right") > np.searchsorted(cs, k_w - win[1], "right")
    return dict(n_hit=int(hit.sum()), frac=float(hit.mean()), null=float(null), excess=float(hit.mean() - null),
                k_locked_frac=float(locked.mean())), lags


def block_pause_ratios(cs, ss, t_start, t_end, block=BLOCK_LEN, return_ccgs=False):
    """Pause ratio per time block; CS refs and SS targets taken from the same block."""
    out, ccgs = [], []
    for t0 in np.arange(t_start, t_end, block):
        t1 = min(t0 + block, t_end)
        c = cs[(cs >= t0) & (cs < t1)]
        s = ss[(ss >= t0 - COARSE_WIN) & (ss < t1 + COARSE_WIN)]
        if len(c) < 20:
            out.append(np.nan)
            ccgs.append(None)
            continue
        rc, lc, _ = ccg(c, s, COARSE_BIN, COARSE_WIN)
        b = rc[_win_mask(lc, COARSE_BIN, BASELINE_WIN)].mean()
        out.append(rc[_win_mask(lc, COARSE_BIN, PAUSE_WIN)].mean() / b if b > 0 else np.nan)
        ccgs.append(rc)
    return (np.array(out), ccgs) if return_ccgs else np.array(out)


# --- v8 literature metrics (section 1b) -------------------------------------------------------------------------------
def latency_test(cs, ss, dead=LAT_DEAD):
    """Shin & De Schutter-style latency test. cs sorted (pair window), ss sorted (pair window ± pad).

    L_post = first SS spike >= dead after each CS, minus the CS; L_pre = CS minus the last SS spike >= dead before it.
    For an SS unrelated to the CS both are the same recurrence-time distribution, so L_pre is the null. The dead zone
    is applied on both sides (symmetric estimator). Only CS with both an L_post and an L_pre are used.
    One-sided Mann-Whitney (L_post > L_pre). NOTE: L_post and L_pre come from the same CS, so the two samples are not
    fully independent (accepted, Simon 2026-11-16). Cliff's delta = P(post > pre) − P(post < pre)."""
    nan = dict(lat_n=0, lat_post_med_ms=np.nan, lat_pre_med_ms=np.nan, lat_ratio=np.nan, lat_p=np.nan, lat_cliff=np.nan)
    if len(cs) < 5 or len(ss) < 5:
        return nan, np.empty(0), np.empty(0)
    i = np.searchsorted(ss, cs + dead, "left")
    j = np.searchsorted(ss, cs - dead, "right") - 1
    ok = (i < len(ss)) & (j >= 0)
    if ok.sum() < 5:
        return nan, np.empty(0), np.empty(0)
    post, pre = ss[i[ok]] - cs[ok], cs[ok] - ss[j[ok]]
    u = mannwhitneyu(post, pre, alternative="greater")
    cliff = 2 * u.statistic / (len(post) * len(pre)) - 1
    mp, mq = float(np.median(post)), float(np.median(pre))
    return dict(lat_n=int(ok.sum()), lat_post_med_ms=mp * 1e3, lat_pre_med_ms=mq * 1e3,
                lat_ratio=mp / mq if mq > 0 else np.nan, lat_p=float(u.pvalue), lat_cliff=float(cliff)), post, pre


def refr_metrics(rate_c, left_c, acg_rate, acg_left):
    """P-sort comparison: SS|CS (CS-triggered SS rate, v7 coarse CCG) vs SS|SS (SS autocorrelogram), same 1 ms bins.
    Both are conditional rates (sp/s per reference spike), so dividing each by the mean SS rate cancels in the ratio.
    refr_ratio = mean SS|CS ÷ mean SS|SS in PAUSE_WIN ("pause beyond refractoriness"; ~1 = only refractoriness).
    below_refr_ms = how long, from PAUSE_WIN[0], SS|CS stays continuously below SS|SS."""
    pw_c, pw_a = _win_mask(left_c, COARSE_BIN, PAUSE_WIN), _win_mask(acg_left, COARSE_BIN, PAUSE_WIN)
    a = acg_rate[pw_a].mean()
    out = dict(refr_ratio=rate_c[pw_c].mean() / a if a > 0 else np.nan)
    lag = np.round(acg_left / COARSE_BIN).astype(int)
    lut = dict(zip(lag, acg_rate))
    k0, n = int(round(PAUSE_WIN[0] / COARSE_BIN)), 0
    for k in range(k0, int(round(REFR_WIN / COARSE_BIN))):
        i = int(np.searchsorted(left_c, k * COARSE_BIN - 1e-12))
        if k not in lut or i >= len(rate_c) or not rate_c[i] < lut[k]:
            break
        n += 1
    out["below_refr_ms"] = n * COARSE_BIN * 1e3
    return out


def lussac_rules(rate_c, left_c):
    """Lussac reference rules on the v7 coarse CCG (1 ms bins, lag 0 at a bin edge, same as SpikeInterface).
    paper v1 §5: mean SS rate [0, 8) ms < 0.5 × mean rate [−20, 0) ms.
    code main (find_purkinje_cells.py, commit 2a0dce1): bins cropped to ±12 ms; pause if median[0, 8) < 0.4 × median of
    the 12 negative-lag bins; rejected if np.median(mirror < thr) is truthy (median over booleans, reproduced as is)."""
    p = _win_mask(left_c, COARSE_BIN, LUSSAC_P1_WIN)
    b = _win_mask(left_c, COARSE_BIN, LUSSAC_P1_BASE)
    pb = rate_c[b].mean()
    paper = bool(pb > 0 and rate_c[p].mean() < LUSSAC_P1_FRAC * pb)
    crop = _win_mask(left_c, COARSE_BIN, (-LUSSAC_CM_HALF, LUSSAC_CM_HALF))
    cc, lc = rate_c[crop], left_c[crop]
    mask = _win_mask(lc, COARSE_BIN, LUSSAC_CM_WIN)
    base = np.median(cc[lc < -1e-12])
    thr = base * LUSSAC_CM_FRAC
    code = bool(np.median(cc[mask]) < thr) and not bool(np.median(cc[mask[::-1]] < thr))
    return dict(lussac_paper_v1_pass=paper, lussac_code_main_pass=code)


def spike_train_checks(t, t0, t1):
    """CS unit checks on one spike train inside [t0, t1) (the unit's own stable span). Flags only."""
    t = t[(t >= t0) & (t < t1)]
    dur = max(t1 - t0, 1e-9)
    isi = np.diff(t)
    fr = len(t) / dur
    f = lambda lo, hi: float(np.mean((isi >= lo) & (isi < hi))) if len(isi) else np.nan
    cv = float(isi.std() / isi.mean()) if len(isi) > 1 and isi.mean() > 0 else np.nan
    cv2 = float(np.mean(2 * np.abs(np.diff(isi)) / (isi[1:] + isi[:-1]))) if len(isi) > 2 else np.nan
    edges = np.arange(t0, t1 + 1e-9, PRESENCE_BIN)
    occ = float(np.mean(np.histogram(t, edges)[0] > 0)) if len(edges) > 1 else np.nan
    i535 = f(*CS_ISI_5_35)
    return dict(scr_fr=fr, cs_fr_ok=bool(CS_FR_RANGE[0] <= fr <= CS_FR_RANGE[1]),
                isi_5_35_frac=i535, isi_5_35_ok=bool(i535 < CS_ISI_5_35_MAX) if np.isfinite(i535) else False,
                cs_isi_10_35_frac=float(np.sum((isi >= CS_ISI_10_35[0]) & (isi < CS_ISI_10_35[1])) / len(t)) if len(t) else np.nan,
                isi_lt50_frac=f(0.0, CS_MIN_IOI), cv=cv, cv2=cv2, presence_bins_frac=occ)


def template_descriptors(t, sr):
    """Columns only. t = (n_samples, n_channels) un-whitened template. Duration trough→following peak on the peak
    channel, footprint channels (ptp > FOOTPRINT_FRAC × peak ptp), late peaks >= 0.5 ms after the trough with
    |value| > TMPL_LATE_MIN × |trough|. Sorter templates span only ~±1.3 ms, so most spikelets are outside them."""
    ptp = t.max(0) - t.min(0)
    pk = int(np.argmax(ptp))
    w = t[:, pk]
    itr = int(np.argmin(w))
    ipk = itr + int(np.argmax(w[itr:]))
    late = np.abs(w[itr + int(round(0.5e-3 * sr)):])
    n_late = len(find_peaks(late, height=TMPL_LATE_MIN * abs(w[itr]))[0]) if len(late) > 2 else 0
    return dict(tmpl_dur_ms=(ipk - itr) / sr * 1e3, tmpl_n_fp_ch=int(np.sum(ptp > FOOTPRINT_FRAC * ptp.max())),
                tmpl_n_late_peaks=int(n_late), tmpl_post_trough_ms=(len(w) - itr) / sr * 1e3)


def lit_pair(cs, ss, rate_c, left_c, acg=None):
    """All v8 pair metrics that need only the two trains + the v7 coarse CCG. acg = (rate, left) of the SS
    autocorrelogram (1 ms bins, ±REFR_WIN) on the same SS spike set; computed here if not given."""
    m, post, pre = latency_test(cs, ss)
    if acg is None:
        acg = ccg(ss, ss, COARSE_BIN, REFR_WIN, acg=True)[:2]
    m.update(refr_metrics(rate_c, left_c, acg[0], acg[1]))
    m.update(lussac_rules(rate_c, left_c))
    return m, post, pre


def lit_candidate(m):
    """Rule 1 (Simon 2026-11-16). Works on a dict or a DataFrame."""
    if isinstance(m, dict):
        return bool(m["lat_p"] < LIT_LAT_P_MAX and m["lat_cliff"] >= LIT_CLIFF_MIN
                    and m["refr_ratio"] < LIT_REFR_RATIO_MAX and not bool(m["collision_dip"]))
    return ((m["lat_p"] < LIT_LAT_P_MAX) & (m["lat_cliff"] >= LIT_CLIFF_MIN) & (m["refr_ratio"] < LIT_REFR_RATIO_MAX)
            & ~m["collision_dip"].astype(bool))


def extra_cs_roles(u, include_manual=INCLUDE_MANUAL_CS, include_screen=INCLUDE_SCREEN_CS):
    """v8 CS sources on a unit table with columns role, manual_CS, cs_screen (bool), region, c4_label.
    Returns (new CS ids, ids kept as SS under option A). Only CB units; never C4 PkC_cs (those are CS already)."""
    want = pd.Series(False, index=u.index)
    if include_manual:
        want |= u["manual_CS"].astype(bool)
    if include_screen:
        want |= u["cs_screen"].astype(bool)
    want &= (u["region"] == REGION) & (u["role"] != "CS")
    kept = u.index[want & u["role"].isin(["SS_cand", "SS_cand2"])].tolist()      # option A
    new = u.index[want & ~u["role"].isin(["SS_cand", "SS_cand2"])].tolist()
    return new, kept


def coincident_mask(ref, tgt, lag, delta=COINC_DELTA):
    """Bool per ref spike: is there a tgt spike in [ref+lag-delta, ref+lag+delta]?"""
    lo = np.searchsorted(tgt, ref + lag - delta, "left")
    hi = np.searchsorted(tgt, ref + lag + delta, "right")
    return hi > lo


def cs_cs_pair(a, b):
    """Duplicate scores for two CS trains (b relative to a), all from one fine CCG."""
    m = {}
    rf, lf, cf = ccg(a, b, FINE_BIN, FINE_WIN)
    centres = lf + FINE_BIN / 2
    k = int(np.argmax(cf))
    lag = centres[k]
    bm = np.abs(centres) > DUP_BASELINE_ABS
    bmu = cf[bm].mean()
    # counts per fine bin are tiny for ~1 Hz trains, so use a Poisson test instead of a z-score
    lam = (cf[bm].sum() + 0.5) / bm.sum()
    m["peak_lag_ms"] = lag * 1e3
    m["peak_count"] = int(cf[k])
    m["base_count_per_bin"] = lam
    m["peak_p"] = min(1.0, poisson.sf(cf[k] - 1, lam) * len(cf))
    half = bmu + 0.5 * (cf[k] - bmu)
    i0 = k
    while i0 > 0 and cf[i0 - 1] >= half:
        i0 -= 1
    i1 = k
    while i1 < len(cf) - 1 and cf[i1 + 1] >= half:
        i1 += 1
    m["peak_fwhm_ms"] = (i1 - i0 + 1) * FINE_BIN * 1e3
    # coincidence relative to the smaller unit, at the peak lag and at 0 for contrast
    small, big, sl = (a, b, lag) if len(a) <= len(b) else (b, a, -lag)
    cm = coincident_mask(small, big, sl)
    m["coinc_frac_peak"] = cm.mean()
    m["coinc_frac_zero"] = coincident_mask(small, big, 0.0).mean()
    # union refractoriness: drop the big unit's partner of each coincident pair, then merge
    partner = coincident_mask(big, small, -sl)
    big_kept = big[~partner]
    union = np.sort(np.concatenate([small, big_kept]))
    n_short_union = int(np.sum(np.diff(union) < CS_MIN_ISI))
    n_short_own = int(np.sum(np.diff(small) < CS_MIN_ISI) + np.sum(np.diff(big_kept) < CS_MIN_ISI))
    m["union_new_short_isi"] = n_short_union - n_short_own
    m["union_new_short_frac"] = m["union_new_short_isi"] / max(len(union), 1)
    return m, (rf, lf)


def poisson_dead(rate, t_end, dead, rng):
    n = int(rate * t_end * 1.5) + 100
    isi = dead + rng.exponential(1 / rate - dead, n)
    t = np.cumsum(isi)
    return t[t < t_end]


def remove_after(ss, ref, w0, w1):
    """SS train with every spike in [ref + w0, ref + w1) removed (imposed pause / dip)."""
    lo = np.searchsorted(ss, ref + w0, "left")
    hi = np.searchsorted(ss, ref + w1, "left")
    bad = np.zeros(len(ss) + 1, int)
    np.add.at(bad, lo, 1)
    np.add.at(bad, hi, -1)
    return ss[np.cumsum(bad)[:-1] == 0]


# --- duplicate evidence on plain arrays (used by section 6b and by the synthetic evidence self-test)
EV_LINES = ("coinc", "lag_sd", "union", "ss_shared", "ss_cond", "footprint", "covar")


def nearest_partner(ref, tgt, lag, half):
    """For each ref spike: index of the tgt spike nearest ref+lag and its offset (NaN if none within half)."""
    j = np.clip(np.searchsorted(tgt, ref + lag), 1, len(tgt) - 1)
    j = np.where(np.abs(tgt[j] - ref - lag) < np.abs(tgt[j - 1] - ref - lag), j, j - 1)
    off = tgt[j] - ref - lag
    return j, np.where(np.abs(off) <= half, off, np.nan)


def detrend(x, t, t0):
    """x minus its per-DETREND_BIN median (slow drift affects both units of a pair)."""
    b = np.floor((t - t0) / DETREND_BIN).astype(int)
    return x - pd.Series(x).groupby(b).transform("median").to_numpy()


def detrended_corr(x, y, t, t0):
    if x is None or y is None or len(x) < 20:
        return np.nan
    dx, dy = detrend(x, t, t0), detrend(y, t, t0)
    return float(np.corrcoef(dx, dy)[0, 1]) if dx.std() > 0 and dy.std() > 0 else np.nan


def pause_after(refs, ss):
    """Pause ratio of SS train ss after the given reference spikes (same windows as pause_metrics)."""
    if len(refs) < MIN_SPIKES_COND:
        return np.nan
    rc, lc, _ = ccg(refs, ss, COARSE_BIN, COARSE_WIN)
    b = rc[_win_mask(lc, COARSE_BIN, BASELINE_WIN)].mean()
    return rc[_win_mask(lc, COARSE_BIN, PAUSE_WIN)].mean() / b if b > 0 else np.nan


def vote(v, plus, minus, higher_is_same=True):
    if v is None or not np.isfinite(v):
        return 0
    if higher_is_same:
        return 1 if v >= plus else (-1 if v < minus else 0)
    return 1 if v <= plus else (-1 if v >= minus else 0)


def dup_evidence_core(ta, tb, lag, coinc_frac_peak, amp_a, amp_b, depth_a, depth_b, ptp_a, ptp_b,
                      ss_a, ss_b, ss_trains, t0, t1, ss_spans=None):
    """Evidence lines for one CS pair, from plain arrays.

    ta, tb: sorted spike times (s) of units a and b in the analysis window; lag: time of b relative to a at
    the fine-CCG peak (s); amp_* / depth_*: per-spike values aligned with ta / tb (depth may be None);
    ptp_*: per-channel template peak-to-peak; ss_a / ss_b: sets of candidate SS partner ids;
    ss_trains: SS id -> spike train (window ± pad); [t0, t1): analysis window; ss_spans (optional): SS id ->
    (start, end) of that SS unit's own window: reference spikes outside it are not used for its pause.
    Returns (evidence dict, aux arrays for plotting).
    """
    def _pa(refs, s):
        if ss_spans is not None:
            refs = refs[(refs >= ss_spans[s][0]) & (refs < ss_spans[s][1])]
        return pause_after(refs, ss_trains[s])

    e = {"coinc_frac_peak": coinc_frac_peak}
    j, off = nearest_partner(ta, tb, lag, LAG_SPREAD_WIN)
    ok = np.isfinite(off)
    e["n_near_peak"] = int(ok.sum())
    e["lag_sd_ms"] = (1.4826 * np.median(np.abs(off[ok] - np.median(off[ok]))) * 1e3) if ok.sum() >= 20 else np.nan
    co = ok & (np.abs(np.nan_to_num(off, nan=np.inf)) <= COINC_DELTA)
    ia, ib = np.nonzero(co)[0], j[co]
    a_only = ta[~coincident_mask(ta, tb, lag)]
    b_only = tb[~coincident_mask(tb, ta, -lag)]
    e.update(n_coinc=len(ia), n_a_only=len(a_only), n_b_only=len(b_only))
    # amplitude / depth co-variation of the coincident spikes
    e["amp_corr_detr"] = detrended_corr(amp_a[ia], amp_b[ib], ta[ia], t0)
    e["depth_corr_detr"] = (detrended_corr(depth_a[ia], depth_b[ib], ta[ia], t0)
                            if depth_a is not None and depth_b is not None else np.nan)
    # slow rate co-variation (a split / hand-over between clusters gives a negative value)
    edges = np.arange(t0, t1 + DETREND_BIN, DETREND_BIN)
    ca, cb = np.histogram(ta, edges)[0], np.histogram(tb, edges)[0]
    e["rate_corr_60s"] = float(np.corrcoef(ca, cb)[0, 1]) if ca.std() > 0 and cb.std() > 0 else np.nan

    # union refractoriness in a range that excludes spikelet doublets
    def n_in(x):
        d = np.diff(x)
        return int(np.sum((d >= UNION_ISI_RANGE[0]) & (d < UNION_ISI_RANGE[1])))
    union = np.sort(np.concatenate([ta, b_only]))
    e["union_new_isi"] = n_in(union) - n_in(ta) - n_in(b_only)
    e["union_new_isi_frac"] = e["union_new_isi"] / max(min(len(ta), len(tb)), 1)
    # expected for two independent trains: pairs (a, b_only) within ±UNION_ISI_RANGE
    expected = len(ta) * len(b_only) * 2 * (UNION_ISI_RANGE[1] - UNION_ISI_RANGE[0]) / (t1 - t0)
    e["union_expected"] = expected
    e["union_new_ratio"] = e["union_new_isi"] / expected if expected >= UNION_MIN_EXPECTED else np.nan
    # within-unit doublets (several detections of one CS)
    e["doublet_frac_a"] = float(np.mean(np.diff(ta) < DOUBLET_ISI)) if len(ta) > 1 else np.nan
    e["doublet_frac_b"] = float(np.mean(np.diff(tb) < DOUBLET_ISI)) if len(tb) > 1 else np.nan
    # template footprints (only comparable where both footprints have signal)
    both = (ptp_a > FOOTPRINT_FRAC * ptp_a.max()) & (ptp_b > FOOTPRINT_FRAC * ptp_b.max())
    e["footprint_overlap_ch"] = int(both.sum())
    e["footprint_cos"] = (float(ptp_a @ ptp_b / np.linalg.norm(ptp_a) / np.linalg.norm(ptp_b))
                          if both.sum() >= FOOTPRINT_MIN_OVERLAP else np.nan)
    # SS co-silencing. ss_cond: pause of one unit's partner after the OTHER unit's non-coincident spikes, only
    # for partners that are not also the other unit's own partner (otherwise it repeats ss_shared)
    sa, sb = set(ss_a), set(ss_b)
    e["ss_a"], e["ss_b"], e["ss_shared"] = sorted(sa), sorted(sb), sorted(sa & sb)
    cond = []
    for owner, part, other_part, other_only in (("a", sa, sb, b_only), ("b", sb, sa, a_only)):
        for s in sorted(part - other_part):
            cond.append((owner, s, _pa(other_only, s)))
    finite = [c for c in cond if np.isfinite(c[2])]
    e["ss_cond_detail"] = ";".join(f"{o}:{s}={v:.2f}" for o, s, v in cond)
    e["ss_cond_pr_min"] = min(c[2] for c in finite) if finite else np.nan
    e["ss_coinc_pr"] = ";".join(f"{s}={_pa(ta[ia], s):.2f}" for s in sorted(sa | sb))
    # votes (lag_sd and covar describe only the coincident spikes: they vote only if coincidence is substantial)
    v = {}
    substantial = coinc_frac_peak >= EV_COINC[1]
    v["coinc"] = vote(coinc_frac_peak, *EV_COINC)
    v["lag_sd"] = (vote(e["lag_sd_ms"] * 1e-3 if np.isfinite(e["lag_sd_ms"]) else np.nan, *EV_LAG_SD,
                        higher_is_same=False) if substantial else 0)
    v["union"] = vote(e["union_new_ratio"], *EV_UNION_NEW, higher_is_same=False)
    v["ss_shared"] = 1 if sa & sb else (-1 if sa and sb else 0)
    v["ss_cond"] = vote(e["ss_cond_pr_min"], *EV_COND_PR, higher_is_same=False)
    v["footprint"] = vote(e["footprint_cos"], *EV_FOOTPRINT)
    cov = np.nanmax([abs(e["amp_corr_detr"]) if np.isfinite(e["amp_corr_detr"]) else np.nan,
                     e["depth_corr_detr"] if np.isfinite(e["depth_corr_detr"]) else np.nan] + [-np.inf])
    v["covar"] = 1 if (substantial and cov >= EV_COVAR) else 0
    e.update({f"ev_{k}": v[k] for k in EV_LINES})
    e["n_support"] = sum(v[k] > 0 for k in EV_LINES)
    e["n_against"] = sum(v[k] < 0 for k in EV_LINES)
    e["support"] = ",".join(k for k in EV_LINES if v[k] > 0)
    e["against"] = ",".join(k for k in EV_LINES if v[k] < 0)
    aux = dict(off=off[ok], ia=ia, ib=ib, a_only=a_only, b_only=b_only, edges=edges, ca=ca, cb=cb,
               union=union, both=both, ta=ta, tb=tb)
    return e, aux


def classify_pair(e, lag_ms, dup_candidate):
    """Descriptive type(s) of a CS pair from its evidence, plus the tests on the path (for the decision view)."""
    lsd = e["lag_sd_ms"] * 1e-3 if np.isfinite(e["lag_sd_ms"]) else np.nan
    t = {}
    t["candidate"] = bool(dup_candidate)
    t["tight"] = bool(np.isfinite(lsd) and lsd <= EV_LAG_SD[0] and e["coinc_frac_peak"] >= EV_COINC[0])
    t["zero_lag"] = abs(lag_ms) <= COINC_DELTA * 1e3
    t["split"] = bool(EV_COINC[1] <= e["coinc_frac_peak"] < 0.9 and (
        (np.isfinite(e["rate_corr_60s"]) and e["rate_corr_60s"] <= -0.3)
        or (np.isfinite(e["amp_corr_detr"]) and e["amp_corr_detr"] <= -EV_COVAR)))
    t["sync"] = bool((np.isfinite(lsd) and lsd >= EV_LAG_SD[1]) or e["ev_ss_shared"] < 0)
    types = []
    if t["tight"]:
        types.append("same_event_redetection" if t["zero_lag"] else "fixed_lag_component")
    if t["split"]:
        types.append("partial_split")
    if t["sync"]:
        types.append("synchrony_like")
    return ("+".join(types) if types else "unclear"), t


def evidence_status(e):
    """Summary of the votes. Not a merge decision: merging is always a manual call."""
    s, a = e["n_support"], e["n_against"]
    if s >= EV_CONVERGE_MIN and a == 0:
        return "converges_same_cell"
    if a >= EV_CONVERGE_MIN and s == 0:
        return "converges_distinct"
    if s + a < 2:
        return "insufficient"
    return "leans_same_cell" if s > a else ("leans_distinct" if a > s else "conflicting")


def run_evidence_self_test():
    """Synthetic CS pairs with known ground truth; checks every evidence line, the type and the status."""
    rng = np.random.default_rng(7)
    T = 1800.0
    ss0 = np.sort(rng.uniform(0, T, int(60 * T)))
    n_ch = 60
    ch = np.arange(n_ch)

    def ptp(center):
        return np.exp(-0.5 * ((ch - center) / 4.0) ** 2)

    def pause(ref):                             # SS silenced 1–12 ms after each reference spike
        return remove_after(ss0, ref, 1e-3, 12e-3)

    def run(name, ta, tb, amp_a, amp_b, ca, cb, ss_a, ss_b, ss_trains):
        o = np.argsort(tb, kind="stable")
        tb, amp_b = tb[o], amp_b[o]
        m, _ = cs_cs_pair(ta, tb)
        e, _ = dup_evidence_core(ta, tb, m["peak_lag_ms"] * 1e-3, m["coinc_frac_peak"], amp_a, amp_b, None, None,
                                 ptp(ca), ptp(cb), ss_a, ss_b, ss_trains, 0.0, T)
        typ, _ = classify_pair(e, m["peak_lag_ms"], m["coinc_frac_peak"] > DUP_COINC_MIN)
        st = evidence_status(e)
        votes = {k: e[f"ev_{k}"] for k in EV_LINES}
        print(f"  {name}: lag {m['peak_lag_ms']:.2f} ms coinc {m['coinc_frac_peak']:.2f} | votes {votes} | "
              f"{typ} | {st}")
        return votes, typ, st

    print("--- evidence self-test")
    cell = poisson_dead(1.0, T, 0.03, rng)
    L = 20 + 3 * rng.standard_normal(len(cell))           # per-event amplitude shared by both detections
    ss1 = pause(cell)

    # (1) same-event re-detection: B = 95% of A at +0.05 ms, both pause SS1, same footprint
    k = rng.random(len(cell)) < 0.95
    v, typ, st = run("same-event", cell, cell[k] + 0.05e-3 + rng.normal(0, 0.02e-3, k.sum()),
                     L + rng.standard_normal(len(cell)), 0.6 * L[k] + rng.standard_normal(k.sum()),
                     30, 30, {1}, {1}, {1: ss1})
    # union is not testable here: B is a subset of A, so merging adds (almost) no spikes
    assert v == dict(coinc=1, lag_sd=1, union=0, ss_shared=1, ss_cond=0, footprint=1, covar=1), v
    assert typ == "same_event_redetection" and st == "converges_same_cell"

    # (2) spikelet lock: B = 80% of A at +1.5 ms, 20% of B detected twice (+1.5 ms doublet)
    k = rng.random(len(cell)) < 0.8
    tb = cell[k] + 1.5e-3 + rng.normal(0, 0.03e-3, k.sum())
    d = rng.random(len(tb)) < 0.2
    tb2 = np.concatenate([tb, tb[d] + 1.5e-3])
    v, typ, st = run("spikelet", cell, tb2, L + rng.standard_normal(len(cell)), rng.normal(10, 2, len(tb2)),
                     30, 31, {1}, {1}, {1: ss1})
    assert v["coinc"] == 1 and v["lag_sd"] == 1 and v["union"] == 1 and v["ss_shared"] == 1, v
    assert typ == "fixed_lag_component" and st == "converges_same_cell"

    # (3) partial split: one cell's events detected by A (85%) and B (20% → 80% over the session);
    #     each event's energy is split between the two templates → anti-correlated amplitudes
    c3 = poisson_dead(1.2, T, 0.03, rng)
    inA = rng.random(len(c3)) < 0.85
    inB = rng.random(len(c3)) < 0.2 + 0.6 * c3 / T
    share = rng.uniform(0.2, 0.8, len(c3))
    L3 = 20 + 2 * rng.standard_normal(len(c3))
    v, typ, st = run("partial split", c3[inA], c3[inB] + rng.normal(0, 0.02e-3, inB.sum()),
                     (L3 * (1 - share))[inA], (L3 * share)[inB], 30, 30, {1}, {1}, {1: pause(c3)})
    assert "partial_split" in typ and v["coinc"] == 1 and v["covar"] == 1 and st == "converges_same_cell", (v, typ)

    # (4) synchronous distinct cells: independent trains, 8% of B's events synchronous with A (±2 ms),
    #     each cell pauses its own SS, overlapping but different footprints
    a4 = poisson_dead(1.0, T, 0.03, rng)
    b4 = poisson_dead(1.0, T, 0.03, rng)
    sync = rng.random(len(b4)) < 0.08
    b4[sync] = rng.choice(a4, sync.sum()) + rng.normal(0, 2e-3, sync.sum())
    b4 = np.sort(b4)
    v, typ, st = run("synchronous", a4, b4, rng.normal(20, 3, len(a4)), rng.normal(20, 3, len(b4)),
                     30, 38, {1}, {2}, {1: pause(a4), 2: pause(b4)})
    # ss_cond: v1–v5 (pause 0–8 ms) gave -1. With the v6 window (4–12 ms) the 8% of B spikes synchronous with A
    # (±2 ms) fall fully inside SS1's silence, so B-only on SS1 reads ~0.90–0.93 (edge of EV_COND_PR): -1 or 0
    assert v["coinc"] == -1 and v["union"] == -1 and v["ss_shared"] == -1 and v["ss_cond"] <= 0, v
    assert v["lag_sd"] == 0 and v["covar"] == 0 and st == "converges_distinct" and "synchrony_like" in typ

    # (5) co-silencing decides: B = 60% of A, only B is linked to SS1, footprints too far apart to compare.
    #     A's spikes that B missed must still silence SS1 → ss_cond +1
    k = rng.random(len(cell)) < 0.6
    v, typ, st = run("co-silencing", cell, cell[k] + 0.05e-3 + rng.normal(0, 0.02e-3, k.sum()),
                     L + rng.standard_normal(len(cell)), 0.6 * L[k] + rng.standard_normal(k.sum()),
                     10, 50, set(), {1}, {1: ss1})
    assert v["ss_cond"] == 1 and v["ss_shared"] == 0 and v["footprint"] == 0 and v["union"] == 0, v
    assert st == "converges_same_cell", st
    print("EVIDENCE SELF-TEST PASSED")


def run_self_test():
    rng = np.random.default_rng(1)
    T = 3600.0
    cs = poisson_dead(1.0, T, 0.02, rng)
    ss0 = np.sort(rng.uniform(0, T, int(60 * T)))

    print("--- self-test (a): imposed one-sided 15 ms pause")
    ss_a = remove_after(ss0, cs, 0.0, 15e-3)
    ma, _, _ = cs_ss_pair(cs, ss_a)
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in ma.items()})
    assert ma["pause_ratio"] < 0.05 and ma["mirror_ratio"] > 0.85 and ma["asym"] > 0.8
    assert 14 <= ma["pause_dur50_ms"] <= 16 and ma["jitter_p"] < 0.01 and ma["jitter_deficit"] > 0.95
    assert ma["lussac_pass"] and not ma["collision_dip"]

    print("--- self-test (b): symmetric ~1 ms collision dip")
    ss_b = remove_after(ss0, cs, -0.5e-3, 0.5e-3)
    mb, _, _ = cs_ss_pair(cs, ss_b)
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in mb.items()})
    assert mb["collision_dip"] and mb["pause_ratio"] > 0.85 and abs(mb["asym"]) < 0.1
    assert not mb["lussac_pass"]
    assert mb["jitter_p"] > CAND_JITTER_P_MAX and abs(mb["jitter_deficit"]) < 0.05, "dip leaks into jitter window"

    print("--- self-test (null): independent SS")
    mn, _, _ = cs_ss_pair(cs, ss0)
    print({k: round(mn[k], 3) for k in ("pause_ratio", "asym", "jitter_p", "jitter_z", "jitter_deficit")})
    assert 0.85 < mn["pause_ratio"] < 1.15 and not mn["lussac_pass"] and abs(mn["jitter_deficit"]) < 0.05

    print("--- self-test (c): duplicate CS train at +0.5 ms")
    keep = rng.random(len(cs)) < 0.95
    dup = np.sort(cs[keep] + 0.5e-3 + rng.normal(0, 0.03e-3, keep.sum()))
    mc, _ = cs_cs_pair(cs, dup)
    print({k: round(v, 3) for k, v in mc.items()})
    assert abs(mc["peak_lag_ms"] - 0.5) <= 0.15 and mc["coinc_frac_peak"] > 0.9
    assert mc["coinc_frac_zero"] < 0.05 and mc["union_new_short_isi"] == 0 and mc["peak_p"] < 1e-6

    print("--- self-test (c-control): independent CS train")
    other = poisson_dead(1.0, T, 0.02, rng)
    mo, _ = cs_cs_pair(cs, other)
    print({k: round(v, 3) for k, v in mo.items()})
    assert mo["coinc_frac_peak"] < 0.05 and mo["union_new_short_isi"] > 10 and mo["peak_p"] > DUP_PEAK_P_MAX

    print("--- self-test (d): CS spike sharing — SS unit carrying 50% of CS spikes vs independent SS")
    half = cs[rng.random(len(cs)) < 0.5] + 0.1e-3          # half the CS spikes, 0.1 ms later
    ss_d = np.sort(np.concatenate([ss0, half]))
    sd, sn = cs_spike_sharing(cs, ss_d, 0.0, T), cs_spike_sharing(cs, ss0, 0.0, T)
    print({k: round(v, 3) for k, v in sd.items()}, "| null:", {k: round(v, 3) for k, v in sn.items()})
    assert 0.4 < sd["cs_share_excess"] < 0.6 and sd["cs_share_excess"] > CS_SHARE_FLAG
    assert abs(sn["cs_share_excess"]) < 0.03 and sn["cs_share_excess"] < CS_SHARE_FLAG

    print("--- self-test (e, v6): SS paused 0–15 ms that also carries a CS piece at +1.8 ms on 60% of CS events")
    piece = cs[rng.random(len(cs)) < 0.6] + 1.8e-3
    ss_e = np.sort(np.concatenate([ss_a, piece]))
    me, _, _ = cs_ss_pair(cs, ss_e)
    ee, lag_e = cs_event_capture(cs, ss_e, 0.0, T)
    print({k: round(me[k], 3) for k in ("pause_ratio", "pause_ratio_0_8", "asym", "jitter_p")},
          {k: round(v, 3) for k, v in ee.items()}, f"median lag {np.median(lag_e) * 1e3:.2f} ms")
    assert me["pause_ratio"] < 0.05 and me["pause_ratio_0_8"] > 0.9, "post-event window must ignore the CS piece"
    # chance uses the unit's overall rate; a unit that pauses after the CS has less real chance, so the excess is
    # conservative (here ~0.38 for a 60% piece): it can miss pieces, it cannot invent them
    assert 0.55 < ee["frac"] < 0.65 and ee["excess"] > 0.3 and abs(np.median(lag_e) - 1.8e-3) < 0.2e-3
    en, _ = cs_event_capture(cs, ss0, 0.0, T)
    assert abs(en["excess"]) < 0.03
    # --- v8 literature metrics
    la, _, _ = lit_pair(cs, ss_a, *ccg(cs, ss_a, COARSE_BIN, COARSE_WIN)[:2])
    print("--- self-test v8 (a): 60 Hz SS, 15 ms pause", {k: round(v, 3) if isinstance(v, float) else v for k, v in la.items()})
    assert la["lussac_paper_v1_pass"] and la["lussac_code_main_pass"] and lit_candidate({**ma, **la})
    lb, _, _ = lit_pair(cs, ss_b, *ccg(cs, ss_b, COARSE_BIN, COARSE_WIN)[:2])
    assert not lb["lussac_code_main_pass"] and not lb["lussac_paper_v1_pass"] and not lit_candidate({**mb, **lb})
    ln, _, _ = lit_pair(cs, ss0, *ccg(cs, ss0, COARSE_BIN, COARSE_WIN)[:2])
    assert not ln["lussac_paper_v1_pass"] and not ln["lussac_code_main_pass"] and abs(ln["lat_cliff"]) < 0.05
    assert not lit_candidate({**mn, **ln})

    print("--- self-test v8 (f): slow SS (15 Hz) with a long true pause (0–25 ms)")
    ss_f = remove_after(poisson_dead(15.0, T, 3e-3, rng), cs, 0.0, 25e-3)
    mf, (rcf, lcf), _ = cs_ss_pair(cs, ss_f)
    lf, _, _ = lit_pair(cs, ss_f, rcf, lcf)
    print({k: round(v, 3) if isinstance(v, float) else v for k, v in {**lf, "pause_ratio": mf["pause_ratio"]}.items()})
    assert lf["lat_p"] < 1e-6 and lf["lat_cliff"] > 0.2 and lf["below_refr_ms"] >= 20 and lf["refr_ratio"] < 0.1 and lit_candidate({**mf, **lf})

    print("--- self-test v8 (g): fast regular SS (100 Hz, ISI 10 ± 0.7 ms); 'CS' = 1 Hz subset of its own spikes + 0.1 ms")
    ss_g = np.cumsum(rng.normal(10e-3, 0.7e-3, int(T * 100)))
    ss_g = ss_g[ss_g < T]
    cs_g = np.sort(rng.choice(ss_g, int(T), replace=False)) + 0.1e-3
    mg, (rcg, lcg), _ = cs_ss_pair(cs_g, ss_g)
    lg, _, _ = lit_pair(cs_g, ss_g, rcg, lcg)
    print({k: round(v, 3) if isinstance(v, float) else v for k, v in {**lg, "pause_ratio": mg["pause_ratio"]}.items()})
    assert mg["pause_ratio"] < 0.3, "the trap: a deep 3–8 ms 'pause' that is only refractoriness"
    assert lg["refr_ratio"] > 0.7 and lg["lat_cliff"] < LIT_CLIFF_MIN and not lit_candidate({**mg, **lg})
    print(f"  (Lussac paper v1 rule on this refractoriness-only dip: {lg['lussac_paper_v1_pass']})")

    print("--- self-test v8 (h): neighbour SS with only a 0–2 ms ephaptic dip")
    ss_h = remove_after(ss0, cs, 0.0, 2e-3)
    mh, (rch, lch), _ = cs_ss_pair(cs, ss_h)
    lh, _, _ = lit_pair(cs, ss_h, rch, lch)
    print({k: round(v, 3) if isinstance(v, float) else v for k, v in {**lh, "pause_ratio": mh["pause_ratio"]}.items()})
    assert mh["pause_ratio"] > 0.85 and abs(lh["lat_cliff"]) < 0.05 and not lit_candidate({**mh, **lh})

    print("--- self-test v8 (i): role assignment (manual-only CS, screen CS, option A)")
    u = pd.DataFrame(dict(role=["CS", "SS_cand", "SS_lowfr", "", "SS_cand2", "SS_lowfr", ""],
                          manual_CS=[True, True, True, True, False, False, False],
                          cs_screen=[False, False, False, False, False, True, False],
                          region=["CB", "CB", "CB", "CB", "CB", "CB", "DCN"],
                          c4_label=["PkC_cs", "PkC_ss", np.nan, np.nan, np.nan, np.nan, np.nan]), index=[1, 2, 3, 4, 5, 6, 7])
    u.loc[7, "manual_CS"] = True
    new, kept = extra_cs_roles(u, True, True)
    assert new == [3, 4, 6] and kept == [2], (new, kept)
    assert extra_cs_roles(u, False, False) == ([], [])
    assert extra_cs_roles(u, True, False)[0] == [3, 4] and extra_cs_roles(u, False, True)[0] == [6]
    print("SELF-TEST PASSED")


set_lab_style()
_t = time.time()
run_self_test()
run_evidence_self_test()
print(f"self-test took {time.time() - _t:.1f} s")


# %% 3. LOAD + INVENTORY + UNIT TABLE
_files = {p.name.lower(): p for p in SORTER_DIR.iterdir() if p.is_file()}


def sorter_file(name, required=True):
    p = _files.get(name.lower())
    if p is None and required:
        raise FileNotFoundError(f"Required file missing in {SORTER_DIR}: {name}")
    return p


INVENTORY = [
    ("params.py", True), ("spike_times.npy", True), ("spike_clusters.npy", True),
    ("spike_templates.npy", True), ("amplitudes.npy", True), ("templates.npy", True),
    ("whitening_mat_inv.npy", True), ("channel_positions.npy", True), ("channel_map.npy", True),
    ("cluster_info.tsv", True), ("cluster_group.tsv", True), ("cluster_Extragood.tsv", True),
    ("cluster_predicted_cell_type.tsv", True), ("cluster_confidence_ratio.tsv", True),
    ("cluster_pred_probability.tsv", True), ("cluster_model_votes.tsv", True),
    ("cluster_channel_region_map.txt", True),
    ("cluster_CS.tsv", False), ("cluster_potential_CS.tsv", False), ("similar_templates.npy", False),
    ("drift.times.npy", False), ("drift.um.npy", False), ("drift_depths.um.npy", False),
]
print(f"--- inventory of {SORTER_DIR}")
for name, req in INVENTORY:
    p = _files.get(name.lower())
    status = f"found   {p.name}  ({p.stat().st_size / 1e6:.1f} MB)" if p else ("MISSING (required)" if req else "absent (optional)")
    print(f"  {name:34s} {status}")
missing = [n for n, r in INVENTORY if r and n.lower() not in _files]
if missing:
    raise FileNotFoundError(f"Required files missing: {missing}")

# params.py: parse key = value lines (never hard-code the sample rate)
PARAMS = {}
for line in sorter_file("params.py").read_text().splitlines():
    if "=" in line:
        k, v = line.split("=", 1)
        try:
            PARAMS[k.strip()] = eval(v.strip(), {}, {})
        except Exception:
            PARAMS[k.strip()] = v.strip()
SR = float(PARAMS["sample_rate"])

_t = time.time()
spike_samples = np.load(sorter_file("spike_times.npy"), mmap_mode="r")
spike_clusters = np.load(sorter_file("spike_clusters.npy"), mmap_mode="r")
spike_templates = np.load(sorter_file("spike_templates.npy"), mmap_mode="r")
amplitudes = np.load(sorter_file("amplitudes.npy"), mmap_mode="r")
templates_w = np.load(sorter_file("templates.npy"), mmap_mode="r")
winv = np.load(sorter_file("whitening_mat_inv.npy"))
chan_pos = np.load(sorter_file("channel_positions.npy"))
chan_map = np.load(sorter_file("channel_map.npy")).ravel()

spike_t = np.asarray(spike_samples, dtype=np.float64) / SR      # s, imec clock
clu = np.asarray(spike_clusters).astype(np.int64)
assert np.all(np.diff(spike_t) >= 0), "spike_times not sorted"

# per-spike depth (duplicate evidence only): KS4 spike_positions.npy (y column) or IBL alf/spikes.depths.npy
SPIKE_DEPTH, _depth_src = None, "none"
for _p, _col in ((sorter_file("spike_positions.npy", False), 1), (SORTER_DIR / "alf" / "spikes.depths.npy", None)):
    if _p is not None and Path(_p).exists():
        _d = np.load(_p, mmap_mode="r")
        if len(_d) == len(spike_t):
            SPIKE_DEPTH, _depth_src = (_d[:, _col] if _col is not None else _d), Path(_p).name
            break
        print(f"  {Path(_p).name} has {len(_d)} rows, spike_times {len(spike_t)}: not used")
print(f"  per-spike depth source: {_depth_src}")
T_REC = spike_t[-1]
print(f"--- loaded in {time.time() - _t:.1f} s")
print(f"  sample_rate {SR:.0f} Hz (params.py) | n spikes {len(spike_t):,} | n clusters {len(np.unique(clu))} "
      f"| duration ≈ {T_REC:.1f} s ({T_REC / 60:.1f} min, last spike) | templates {templates_w.shape}")


def first_existing(candidates, what):
    for c in candidates:
        if c.exists():
            return c
    tried = "\n  ".join(str(c) for c in candidates)
    raise FileNotFoundError(f"No {what} found. Tried:\n  {tried}")


def ephys_to_harp_fit(path, n_check=3):
    """Linear ephys(imec)→Harp map from neurons_session.mat (v7.3): paired st_raw_all / st_all.

    s01 maps every spike with one linear function, so any unit gives the coefficients; a few units
    spread over the struct are fitted to check that they agree. Opened read-only.
    """
    fits = []
    with h5py.File(path, "r") as h:
        raw, harp = h["neurons"]["st_raw_all"], h["neurons"]["st_all"]
        for i in np.linspace(0, raw.shape[0] - 1, n_check).astype(int):
            x, y = h[raw[i, 0]][()].ravel(), h[harp[i, 0]][()].ravel()
            if len(x) < 2 or len(x) != len(y):
                continue
            p = np.polyfit(x, y, 1)
            fits.append((p[0], p[1], float(np.abs(y - np.polyval(p, x)).max())))
    fits = np.array(fits)
    assert len(fits), "no usable unit for the clock fit"
    assert np.ptp(fits[:, 0]) < 1e-9 and np.ptp(fits[:, 1]) < 1e-4 and fits[:, 2].max() < 1e-4, \
        f"clock fits disagree across units: {fits}"
    return float(fits[0, 0]), float(fits[0, 1]), float(fits[:, 2].max())


def read_win_selected(path):
    """bhv.win_selected [select_start_time, time_sess_end] (Harp s) and trials_vec. Read-only."""
    try:
        b = sio.loadmat(path, squeeze_me=True, struct_as_record=False, variable_names=["bhv"])["bhv"]
        w = b.win_selected
        return float(w.select_start_time), float(w.time_sess_end), np.atleast_1d(w.trials_vec)
    except NotImplementedError:  # v7.3
        with h5py.File(path, "r") as h:
            w = h["bhv"]["win_selected"]
            return (float(w["select_start_time"][()].ravel()[0]), float(w["time_sess_end"][()].ravel()[0]),
                    w["trials_vec"][()].ravel())


BHV_EVENT_FIELDS = ("initiation_times_all", "inva_push_times_all", "reach_times")


def read_bhv_span(path):
    """First and last behavioural event (Harp s) over bhv.all_trials (every trial, not only win_selected)."""
    vals = []
    try:
        b = sio.loadmat(path, squeeze_me=True, struct_as_record=False, variable_names=["bhv"])["bhv"]
        for f in BHV_EVENT_FIELDS:
            if hasattr(b.all_trials, f):
                vals.append(np.atleast_1d(np.asarray(getattr(b.all_trials, f), float)).ravel())
    except NotImplementedError:  # v7.3
        with h5py.File(path, "r") as h:
            g = h["bhv"]["all_trials"]
            vals = [g[f][()].ravel().astype(float) for f in BHV_EVENT_FIELDS if f in g]
    v = np.concatenate(vals) if vals else np.array([])
    v = v[np.isfinite(v)]
    assert len(v), f"no behavioural event times in bhv.all_trials {BHV_EVENT_FIELDS}"
    return float(v.min()), float(v.max())


SYNC = {}
if ANALYSIS_WINDOW == "full":
    T0, T1 = 0.0, T_REC
    B0, B1 = 0.0, T_REC
elif ANALYSIS_WINDOW == "bhv_win_selected":
    _bhv_path = first_existing(BHV_CANDIDATES, "behavior_fundamentals.mat")
    _ns_path = first_existing(NEURONS_SESSION_CANDIDATES, "neurons_session.mat")
    h0, h1, _trials = read_win_selected(_bhv_path)
    a, b, maxres = ephys_to_harp_fit(_ns_path)
    T0, T1 = (h0 - b) / a, (h1 - b) / a
    hb0, hb1 = read_bhv_span(_bhv_path)
    B0, B1 = max(0.0, min((hb0 - b) / a, T0)), min(T_REC, max((hb1 - b) / a, T1))
    SYNC = dict(bhv_file=str(_bhv_path), neurons_session_file=str(_ns_path), harp_slope=a, harp_intercept=b,
                fit_max_residual_s=maxres, win_harp=[h0, h1], trials=[int(_trials.min()), int(_trials.max())],
                bhv_span_harp=[hb0, hb1], bhv_span_imec=[B0, B1])
    print(f"  behaviour period (bhv.all_trials first/last event): imec [{B0:.1f}, {B1:.1f}] s")
    print(f"  clock map: harp = {a:.12f} * imec + {b:.6f} (max residual {maxres:.1e} s) from {_ns_path.name}")
    print(f"  bhv.win_selected: trials {int(_trials.min())}–{int(_trials.max())}, Harp [{h0:.3f}, {h1:.3f}] s")
    assert 0 <= T0 < T1 <= T_REC + 1, f"window [{T0:.1f}, {T1:.1f}] s is outside the recording (0–{T_REC:.1f} s)"
else:
    raise ValueError(f"unknown ANALYSIS_WINDOW {ANALYSIS_WINDOW!r}")
T_WIN = T1 - T0
print(f"  ANALYSIS WINDOW ({ANALYSIS_WINDOW}): imec [{T0:.1f}, {T1:.1f}] s = {T_WIN:.1f} s "
      f"({T_WIN / T_REC:.0%} of the recording)")

# per-cluster spike times (stable sort keeps time order within a cluster)
#   ST_FULL / SPK_IDX_FULL: whole recording (overview figures only)
#   ST / SPK_IDX:           analysis window (every metric); re-cut to each unit's stable span in section 3b
#   ST_PAD:                 analysis window ± EDGE_PAD (CCG targets, so that refs near the edges see full lags)
_order = np.argsort(clu, kind="stable")
_ids, _starts, _counts = np.unique(clu[_order], return_index=True, return_counts=True)
ST_FULL = {int(c): spike_t[_order[s:s + n]] for c, s, n in zip(_ids, _starts, _counts)}
SPK_IDX_FULL = {int(c): _order[s:s + n] for c, s, n in zip(_ids, _starts, _counts)}
_in = {c: (t >= T0) & (t < T1) for c, t in ST_FULL.items()}
ST = {c: t[_in[c]] for c, t in ST_FULL.items()}
SPK_IDX = {c: i[_in[c]] for c, i in SPK_IDX_FULL.items()}
ST_PAD = {c: t[(t >= T0 - EDGE_PAD) & (t < T1 + EDGE_PAD)] for c, t in ST_FULL.items()}

# cluster templates: spike-count-weighted mean of the templates each cluster uses (window spikes), un-whitened
n_tmpl = templates_w.shape[0]
_win_spk = (spike_t >= T0) & (spike_t < T1)
_tmpl = np.asarray(spike_templates).astype(np.int64)
_w = np.bincount(clu[_win_spk] * n_tmpl + _tmpl[_win_spk], minlength=(clu.max() + 1) * n_tmpl).reshape(-1, n_tmpl)
_w_full = np.bincount(clu * n_tmpl + _tmpl, minlength=(clu.max() + 1) * n_tmpl).reshape(-1, n_tmpl)
_no_win = [int(c) for c in _ids if _w[c].sum() == 0]
for c in _no_win:  # silent in the window: template from the whole recording (unit is excluded below anyway)
    _w[c] = _w_full[c]
print(f"  {len(_no_win)} clusters have no spikes in the analysis window")
del _tmpl
print(f"  templates have {templates_w.shape[2]} channels; channel_positions {chan_pos.shape}")


def cluster_template(c):
    """Un-whitened (n_samples, n_channels) mean template of cluster c."""
    w = _w[c]
    used = np.nonzero(w)[0]
    t = np.tensordot(w[used] / w[used].sum(), np.asarray(templates_w[used]), axes=1)
    return t @ winv


TEMPLATE = {}
rows = []
for c in _ids:
    c = int(c)
    t = cluster_template(c)
    ptp = t.max(0) - t.min(0)
    pk = int(np.argmax(ptp))
    TEMPLATE[c] = t
    st = ST[c]
    isi = np.diff(st)
    rows.append(dict(phyID=c, peak_ch=pk, x=chan_pos[pk, 0], y=chan_pos[pk, 1],
                     n_spikes=len(st), fr=len(st) / T_WIN,
                     frac_isi_lt1ms=float(np.mean(isi < 1e-3)) if len(isi) else np.nan,
                     n_templates_used=int(np.count_nonzero(_w[c])),
                     n_spikes_full=len(ST_FULL[c]), fr_full=len(ST_FULL[c]) / T_REC,
                     **template_descriptors(t, SR)))                                   # v8 columns
units = pd.DataFrame(rows).set_index("phyID")


def read_tsv(name, col=None):
    p = sorter_file(name, required=False)
    if p is None:
        return None
    df = pd.read_csv(p, sep="\t")
    df = df.rename(columns={df.columns[0]: "phyID"}).set_index("phyID")
    return df if col is None else df.iloc[:, 0].rename(col)


units = units.join(read_tsv("cluster_group.tsv", "group"))
eg = read_tsv("cluster_Extragood.tsv", "extragood")
units["extragood"] = units.index.isin(eg[eg.astype(str).str.lower() == "y"].index)
units = units.join(read_tsv("cluster_predicted_cell_type.tsv", "c4_label"))
units = units.join(read_tsv("cluster_confidence_ratio.tsv", "c4_conf"))
units = units.join(read_tsv("cluster_pred_probability.tsv", "c4_prob"))
units = units.join(read_tsv("cluster_model_votes.tsv", "c4_votes"))
_v = units["c4_votes"].dropna().astype(str).str.split("/", expand=True).astype(float)
units["c4_vote_ratio"] = (_v[0] / _v[1]).reindex(units.index)
reg = pd.read_csv(sorter_file("cluster_channel_region_map.txt"), sep="\t").set_index("cluster_id")
units["region"] = reg["brain_region"].reindex(units.index)
units["region_ch"] = reg["channel"].reindex(units.index)
man_file = "cluster_CS.tsv" if sorter_file("cluster_CS.tsv", False) else "cluster_potential_CS.tsv"
man = read_tsv(man_file, "manual_CS")
units["manual_CS"] = units.index.isin(man[man.astype(str).str.lower() == "y"].index) if man is not None else False
print(f"  manual CS flags from: {man_file if man is not None else 'none'}")
MANUAL_FILE = man_file if man is not None else "none"
# v8 confound inputs: the sorter's own isolation scores (KS4 ContamPct %, IBL noise_predict); columns only
_cp = read_tsv("cluster_ContamPct.tsv", "sorter_contam_pct")
units["sorter_contam_pct"] = pd.to_numeric(_cp, errors="coerce").reindex(units.index) if _cp is not None else np.nan
_pr = read_tsv("cluster_predictions.tsv")
units["ibl_noise_pred"] = (pd.to_numeric(_pr["noise_predict"], errors="coerce").reindex(units.index)
                           if _pr is not None and "noise_predict" in _pr.columns else np.nan)

# sanity checks against phy's cluster_info.tsv
ci = pd.read_csv(sorter_file("cluster_info.tsv"), sep="\t").set_index("cluster_id")
_both = units.index.intersection(ci.index)
print(f"  peak_ch == cluster_info.ch for {np.mean(units.loc[_both, 'peak_ch'] == ci.loc[_both, 'ch']):.0%} of clusters; "
      f"region-map channel == peak_ch for {np.mean(units['region_ch'] == units['peak_ch']):.0%}")

# population and roles
pop = units["extragood"] & (units["region"] == REGION)
units["role"] = ""
units.loc[pop & (units["c4_label"] == CS_LABEL), "role"] = "CS"
units.loc[pop & (units["c4_label"] != CS_LABEL), "role"] = "SS_cand"
units.loc[(units["role"] == "SS_cand") & (units["fr"] < SS1_FR_MIN), "role"] = "SS_lowfr"     # v6 tier-1 FR floor
LOWFR_IDS = units.index[units["role"] == "SS_lowfr"].tolist()
_gate = True if SS2_GROUPS is None else units["group"].isin(SS2_GROUPS)
pop2 = ((units["region"] == REGION) & ~units["extragood"] & _gate
        & (units["fr"] >= SS2_FR_MIN) & (units["c4_label"] != CS_LABEL))
units.loc[pop2, "role"] = "SS_cand2"
CS_IDS = units.index[units["role"] == "CS"].tolist()
SS_IDS = units.index[units["role"].isin(["SS_cand", "SS_cand2"])].tolist()
if _os.environ.get("STEP0_TEST_NO_SS"):          # test hook only: run the session as if no SS candidate existed
    SS_IDS = []
    print("STEP0_TEST_NO_SS: SS candidate list forced empty (test run)")
print(f"--- tier-1 FR floor (>= {SS1_FR_MIN} Hz): {len(LOWFR_IDS)} extra-good non-CS units excluded "
      f"(role SS_lowfr): " + ", ".join(f"{c} ({units.at[c, 'fr']:.1f} Hz)" for c in LOWFR_IDS))
print(f"--- tier-2 SS candidates ({SS2_GROUP_TXT}, not extragood, {REGION}, FR >= {SS2_FR_MIN} Hz, "
      f"not {CS_LABEL}): {int(pop2.sum())}")
print("  non-extragood units labelled", CS_LABEL, "(not used as CS or SS):",
      units.index[~units["extragood"] & (units["c4_label"] == CS_LABEL)].tolist())

print(f"--- population: extragood & region == {REGION}: {pop.sum()} units "
      f"(extragood total {units['extragood'].sum()})")
print("  extragood units per region:", units.loc[units["extragood"], "region"].value_counts(dropna=False).to_dict())
print("  C4 label within population:")
print(units.loc[pop, "c4_label"].fillna("unlabelled").value_counts().to_string())
print(f"  n CS = {len(CS_IDS)} | n SS candidates = {len(SS_IDS)}")
print("--- crosstab C4 label vs manual CS flag (population; diagnostic only)")
print(pd.crosstab(units.loc[pop, "c4_label"].fillna("unlabelled"), units.loc[pop, "manual_CS"], margins=True))
print("  manual CS but C4 != PkC_cs:", units.index[pop & units["manual_CS"] & (units["c4_label"] != CS_LABEL)].tolist())
print("  C4 PkC_cs but not manual CS:", units.index[pop & ~units["manual_CS"] & (units["c4_label"] == CS_LABEL)].tolist())
print("  manual CS outside population (not extragood or not CB):",
      units.index[~pop & units["manual_CS"]].tolist())

cols = ["role", "group", "extragood", "c4_label", "c4_conf", "c4_prob", "c4_votes", "region", "manual_CS",
        "peak_ch", "x", "y", "n_spikes", "fr", "frac_isi_lt1ms", "n_templates_used"]
units = units[cols + [c for c in units.columns if c not in cols]]
units.to_csv(TAB_DIR / "unit_table.csv")
with pd.option_context("display.width", 200, "display.max_columns", 30, "display.max_rows", 200):
    print("--- population unit table")
    print(units.loc[pop, cols].sort_values(["role", "y"]).round(3).to_string())
print(f"saved {TAB_DIR / 'unit_table.csv'}")


# %% 3b. PER-UNIT STABLE SPANS (reference = analysis window; extension only inside the behaviour period)
#   ST_REF / ST_PAD_REF / SPK_IDX_REF: the reference window [T0, T1) for every unit (kept for comparison columns)
#   SPAN[c] = (start, end): unit c's own window. ST / ST_PAD / SPK_IDX are re-cut to it for CS and SS units.
#   Pairs use pair_win(a, b), the overlap of the two spans, and cut both trains to it with win_st / win_idx.
ST_REF, ST_PAD_REF, SPK_IDX_REF = dict(ST), dict(ST_PAD), dict(SPK_IDX)   # copies: ST etc. are re-cut below
STAB = {}


def win_idx(c, w):
    """Indices into ST_FULL[c] / SPK_IDX_FULL[c] of the spikes in [w0, w1)."""
    t = ST_FULL[c]
    return slice(int(np.searchsorted(t, w[0], "left")), int(np.searchsorted(t, w[1], "left")))


def win_st(c, w, pad=0.0):
    t = ST_FULL[c]
    return t[int(np.searchsorted(t, w[0] - pad, "left")):int(np.searchsorted(t, w[1] + pad, "left"))]


def stab_bins():
    """Evaluation grid: STAB_BIN bins stepping left from T0 and right from T1, plus the reference bins."""
    left = [(max(x - STAB_BIN, 0.0), x) for x in np.arange(T0, 0, -STAB_BIN)]
    inside = [(x, min(x + STAB_BIN, T1)) for x in np.arange(T0, T1, STAB_BIN)]
    right = [(x, min(x + STAB_BIN, T_REC)) for x in np.arange(T1, T_REC, STAB_BIN)]
    return left, inside, right


def stable_span(c):
    """(span, per-bin table) for unit c. Bins are scored against the unit's reference FR and median amplitude."""
    t = ST_FULL[c]
    amp = np.asarray(amplitudes[SPK_IDX_FULL[c]], float).ravel()
    ref = (t >= T0) & (t < T1)
    fr_ref, amp_ref = ref.sum() / T_WIN, (np.median(amp[ref]) if ref.any() else np.nan)
    left, inside, right = stab_bins()
    rows = []

    def score(x0, x1, side):
        s = win_idx(c, (x0, x1))
        n = s.stop - s.start
        fr_r = (n / (x1 - x0)) / fr_ref if fr_ref > 0 else np.nan
        amp_r = np.median(amp[s]) / amp_ref if n and amp_ref > 0 else np.nan
        in_bhv = (x0 >= B0 - 1e-9) and (x1 <= B1 + 1e-9) and (x1 - x0) >= STAB_MIN_BIN
        ok = bool(np.isfinite(fr_r) and np.isfinite(amp_r) and abs(np.log(fr_r)) <= np.log(STAB_FR_TOL)
                  and abs(np.log(amp_r)) <= np.log(STAB_AMP_TOL))
        rows.append(dict(t0=x0, t1=x1, side=side, fr_ratio=fr_r, amp_ratio=amp_r, in_bhv=in_bhv, ok=ok))
        return in_bhv, ok

    def walk(bins, side):
        edge, gap = None, 0
        stop = False
        for x0, x1 in bins:
            in_bhv, ok = score(x0, x1, side)
            if stop or not in_bhv:
                stop = True
                continue
            if ok:
                edge, gap = (x0 if side == "left" else x1), 0
            else:
                gap += 1
                stop = gap > STAB_MAX_GAP
        return edge

    for x0, x1 in inside:
        score(x0, x1, "ref")
    e0 = walk(left, "left") if STABLE_SPAN else None
    e1 = walk(right, "right") if STABLE_SPAN else None
    span = (T0 if e0 is None else e0, T1 if e1 is None else e1)
    return span, pd.DataFrame(rows).sort_values("t0").reset_index(drop=True)


# v8: the extra CS sources need spans too: manual units (CB, any label, extra-good or not) and the screen pool
#   (C4-unlabelled extra-good CB units; + non-extra-good CB units if SCREEN_NON_EG). Spans depend only on the unit.
_cb = units["region"] == REGION
SCREEN_POOL = units.index[_cb & units["c4_label"].isna() & (units["extragood"] | SCREEN_NON_EG)].tolist()
CALIB_POOL = units.index[_cb & units["extragood"] & units["c4_label"].notna()].tolist()   # screen as columns only
MANUAL_IDS = units.index[_cb & units["manual_CS"].astype(bool)].tolist()
_V7_SPAN_UNITS = CS_IDS + SS_IDS + LOWFR_IDS
_SPAN_UNITS = list(dict.fromkeys(_V7_SPAN_UNITS + MANUAL_IDS + SCREEN_POOL + CALIB_POOL))
SPAN = {c: (T0, T1) for c in ST_FULL}


def _apply_span(c, span, tab):
    SPAN[c], STAB[c] = span, tab
    if STABLE_SPAN and ANALYSIS_WINDOW == "bhv_win_selected":
        ST[c] = win_st(c, SPAN[c])
        SPK_IDX[c] = SPK_IDX_FULL[c][win_idx(c, SPAN[c])]
        ST_PAD[c] = win_st(c, SPAN[c], EDGE_PAD)
    else:
        SPAN[c] = (T0, T1)


for c in _V7_SPAN_UNITS:                         # exactly the v7 set: spans of every other unit stay (T0, T1)
    _apply_span(c, *stable_span(c))
_XSPAN = {c: stable_span(c) for c in _SPAN_UNITS if c not in STAB}   # extra units: kept aside, applied only if CS
_span_of = lambda c: SPAN[c] if c in STAB else (_XSPAN[c][0] if STABLE_SPAN and ANALYSIS_WINDOW == "bhv_win_selected"
                                                else (T0, T1))

# --- v8 CS unit checks / screen (own stable span), for CS, manual units, the screen pool and the calibration pool
_chk = pd.DataFrame({c: spike_train_checks(ST_FULL[c], *_span_of(c)) for c in _SPAN_UNITS}).T
for _k in _chk.columns:
    units[_k] = _chk[_k].reindex(units.index)
for _k in ("cs_fr_ok", "isi_5_35_ok"):
    units[_k] = units[_k].astype("boolean")
units["screen_pass"] = (units["cs_fr_ok"] & units["isi_5_35_ok"]).fillna(False).astype(bool)   # entry rule
units["cs_screen"] = units.index.isin(SCREEN_POOL) & units["screen_pass"]
units["c4_soft_available"] = units["c4_label"].notna()        # C4 files hold rows (top class only) for labelled units
units["span_frac"] = [((_span_of(c)[1] - _span_of(c)[0]) / T_WIN) if c in _SPAN_UNITS else np.nan for c in units.index]
units["partial_presence"] = (units["span_frac"] < PRESENCE_MIN_SPAN) | (units["presence_bins_frac"] < PRESENCE_MIN_BINS)
units["poor_isolation"] = ((units["frac_isi_lt1ms"] > ISO_ISI1_MAX) | (units["sorter_contam_pct"] > ISO_CONTAM_MAX)
                           | (units["ibl_noise_pred"] > ISO_NOISE_PRED_MAX))

# --- v8 extra CS (manual / screen). Option A: v7 SS candidates stay SS. Units leaving the low-FR audit are reported.
_v7_cs = set(CS_IDS)
NEW_CS, KEPT_AS_SS = extra_cs_roles(units)
_left_lowfr = [c for c in NEW_CS if units.at[c, "role"] == "SS_lowfr"]
units["cs_c4"] = units.index.isin(_v7_cs)
units["cs_manual"] = units["manual_CS"].astype(bool) & _cb
units.loc[NEW_CS, "role"] = "CS"
for c in NEW_CS:
    if c in _XSPAN:
        _apply_span(c, *_XSPAN[c])
CS_IDS = units.index[units["role"] == "CS"].tolist()
LOWFR_IDS = units.index[units["role"] == "SS_lowfr"].tolist()
assert not set(CS_IDS) & set(SS_IDS), "a CS unit is also in the SS candidate pool"


def _source(r):
    parts = [n for n, k in (("C4", "cs_c4"), ("manual", "cs_manual"), ("screen", "cs_screen")) if bool(r[k])]
    if r.name not in CS_IDS:
        return ""
    return "C4" if parts == ["C4"] else "+".join(parts) if parts[0] == "C4" else (parts[0] + "_only" if len(parts) == 1
                                                                                  else "+".join(parts))


units["cs_source"] = units.apply(_source, axis=1)
units["manual_kept_as_ss"] = units.index.isin(KEPT_AS_SS)
print(f"--- v8 CS sources (INCLUDE_MANUAL_CS={INCLUDE_MANUAL_CS}, INCLUDE_SCREEN_CS={INCLUDE_SCREEN_CS}; manual file: "
      f"{MANUAL_FILE}): {len(CS_IDS)} CS = " + str(units.loc[CS_IDS, "cs_source"].value_counts().to_dict()))
print(f"  added as CS: {NEW_CS} | kept as SS (option A): {KEPT_AS_SS} | left the low-FR SS audit: {_left_lowfr}")
units["span_start"] = [SPAN[c][0] for c in units.index]
units["span_end"] = [SPAN[c][1] for c in units.index]
units["span_s"] = units["span_end"] - units["span_start"]
units["n_spikes_span"] = [len(ST[c]) for c in units.index]
units.to_csv(TAB_DIR / "unit_table.csv")
pd.concat([d.assign(unit=c) for c, d in STAB.items()]).to_csv(TAB_DIR / "unit_stability_bins.csv", index=False)


def pair_win(a, b):
    return max(SPAN[a][0], SPAN[b][0]), min(SPAN[a][1], SPAN[b][1])


def span_cut(ref, c):
    """Reference spikes inside unit c's span (a pause is only measured where c itself is stable)."""
    return ref[(ref >= SPAN[c][0]) & (ref < SPAN[c][1])]


_sp = units.loc[CS_IDS + SS_IDS]
print(f"--- stable spans ({'extended' if STABLE_SPAN else 'off'}; FR ×/÷{STAB_FR_TOL}, amp ×/÷{STAB_AMP_TOL}, "
      f"{STAB_BIN:.0f} s bins, behaviour [{B0:.0f}, {B1:.0f}] s, reference [{T0:.0f}, {T1:.0f}] s)")
for _role, _ids in (("CS", CS_IDS), ("SS", SS_IDS)):
    _s = _sp.loc[_ids, "span_s"]
    print(f"  {_role}: span / reference length median {np.median(_s / T_WIN):.2f} (min {(_s / T_WIN).min():.2f}, "
          f"max {(_s / T_WIN).max():.2f}); extended before T0: {int((_sp.loc[_ids, 'span_start'] < T0 - 1).sum())}, "
          f"after T1: {int((_sp.loc[_ids, 'span_end'] > T1 + 1).sum())} of {len(_ids)}")


# %% 4. CS OVERVIEW (one saved figure per CS; figures are not shown)
drift_t = drift_um = drift_depths = None
if sorter_file("drift.um.npy", False) is not None:
    drift_t = np.load(sorter_file("drift.times.npy"))
    drift_um = np.load(sorter_file("drift.um.npy"))
    drift_depths = np.load(sorter_file("drift_depths.um.npy")).ravel()


def bar_ccg(ax, rate, left, bin_s, color=BAR_COLOR):
    ax.bar(left * 1e3, rate, width=bin_s * 1e3, align="edge", color=color, alpha=0.9, linewidth=0)
    ax.set_xlim(left[0] * 1e3, (left[-1] + bin_s) * 1e3)


def shade_windows(ax):
    ax.axvspan(PAUSE_WIN[0] * 1e3, PAUSE_WIN[1] * 1e3, color=SHADE, zorder=0, linewidth=0)
    ax.axvspan(BASELINE_WIN[0] * 1e3, BASELINE_WIN[1] * 1e3, color=SHADE, alpha=0.5, zorder=0, linewidth=0)


def label_ccg(ax, xlabel=True, ylabel=True):
    if xlabel:
        ax.set_xlabel("time lag (ms)")
    if ylabel:
        ax.set_ylabel("spike rate (sp/s)")


def plot_template(ax, c, color=CB_COLOR, n_ch=2 * N_NEIGH_CH + 1):
    """Un-whitened template on the n_ch channels closest to the peak, laid out by probe geometry."""
    t = TEMPLATE[c]
    pk = int(units.at[c, "peak_ch"])
    d = np.hypot(*(chan_pos - chan_pos[pk]).T)
    chs = np.argsort(d)[:n_ch]
    amp = max(np.ptp(t[:, pk]), 1e-9)
    tt = np.arange(t.shape[0]) / SR * 1e3
    tspan = tt[-1] * 1.2
    for ch in chs:
        x0 = (chan_pos[ch, 0] - chan_pos[pk, 0]) / 32 * tspan
        y0 = (chan_pos[ch, 1] - chan_pos[pk, 1]) / 20 * amp * 0.8
        ax.plot(tt + x0, t[:, ch] + y0, color=color if ch == pk else (0.5, 0.5, 0.5), lw=1.2)
    ax.set_title(f"template, peak ch {pk}", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def plot_probe(ax, focus, cands=(), selected=None):
    ax.scatter(chan_pos[:, 0], chan_pos[:, 1], s=4, color=(0.85, 0.85, 0.85), marker="s")
    jit = {c: (i % 5 - 2) * 3 for i, c in enumerate(units.index)}
    for c in SS_IDS:
        ax.scatter(units.at[c, "x"] + jit[c], units.at[c, "y"], s=12, color=(0.6, 0.6, 0.6))
    for c in CS_IDS:
        ax.scatter(units.at[c, "x"] + jit[c], units.at[c, "y"], s=12, color=CB_COLOR, alpha=0.5)
    for c in cands:
        ax.scatter(units.at[c, "x"] + jit[c], units.at[c, "y"], s=40, facecolor="none", edgecolor=BAR_COLOR, lw=1.2)
    if selected is not None:
        ax.scatter(units.at[selected, "x"] + jit[selected], units.at[selected, "y"], s=90, facecolor="none",
                   edgecolor="k", lw=2)
    ax.scatter(units.at[focus, "x"] + jit[focus], units.at[focus, "y"], s=90, color=CB_COLOR, edgecolor="k", lw=1)
    y0 = units.at[focus, "y"]
    ax.set_ylim(y0 - 300, y0 + 300)
    ax.set_xlim(chan_pos[:, 0].min() - 20, chan_pos[:, 0].max() + 20)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_title("probe map", fontsize=10)


def plot_cs_overview(c):
    st = ST[c]
    fig = Figure(figsize=(15, 8))
    gs = fig.add_gridspec(3, 4, hspace=0.6, wspace=0.35)
    plot_template(fig.add_subplot(gs[0, 0]), c)
    ax = fig.add_subplot(gs[0, 1])
    r, l, _ = ccg(st, st, ACG_FINE_BIN, ACG_FINE_WIN, acg=True)
    ax.axvspan(-5, 5, color=SHADE, zorder=0, linewidth=0)
    bar_ccg(ax, r, l, ACG_FINE_BIN)
    ax.set_title("Autocorrelogram (fine); grey = ±5 ms", fontsize=10)
    label_ccg(ax)
    ax = fig.add_subplot(gs[0, 2])
    r, l, _ = ccg(st, st, ACG_COARSE_BIN, ACG_COARSE_WIN, acg=True)
    bar_ccg(ax, r, l, ACG_COARSE_BIN)
    ax.set_title("Autocorrelogram (coarse)", fontsize=10)
    label_ccg(ax)
    ax = fig.add_subplot(gs[0, 3])
    edges = np.logspace(-3, np.log10(5), 60)
    ax.hist(np.diff(st), bins=edges, color=BAR_COLOR, alpha=0.9)
    ax.set_xscale("log")
    ax.set_xlabel("ISI (s)")
    ax.set_ylabel("count")
    ax.set_title(f"ISI (window); {np.mean(np.diff(st) < CS_MIN_ISI):.2%} < {CS_MIN_ISI * 1e3:.0f} ms", fontsize=10)
    # amplitude across the whole recording (analysis window shaded)
    ax = fig.add_subplot(gs[1, 0:3])
    ax.axvspan(*SPAN[c], color=SHADE, alpha=0.3, zorder=0, linewidth=0)
    ax.axvspan(T0, T1, color=SHADE, alpha=0.5, zorder=0, linewidth=0)
    ax.scatter(ST_FULL[c], np.asarray(amplitudes[SPK_IDX_FULL[c]]), s=2, color=DOT_COLOR, rasterized=True)
    ax.set_ylabel("amplitude (a.u.)")
    ax.set_xlim(0, T_REC)
    ax.set_title(f"amplitude across the session; dark = reference window ({ANALYSIS_WINDOW}), "
                 f"light = this unit's stable span [{SPAN[c][0]:.0f}, {SPAN[c][1]:.0f}] s", fontsize=10)
    # FR + drift
    ax = fig.add_subplot(gs[2, 0:3])
    ax.axvspan(*SPAN[c], color=SHADE, alpha=0.3, zorder=0, linewidth=0)
    ax.axvspan(T0, T1, color=SHADE, alpha=0.5, zorder=0, linewidth=0)
    e = np.arange(0, T_REC + 30, 30)
    fr = np.histogram(ST_FULL[c], e)[0] / 30
    ax.fill_between(e[:-1], fr, step="post", color=FR_FILL)
    ax.set_ylabel("FR (sp/s), 30 s bins")
    ax.set_xlabel("time (s, imec clock)")
    ax.set_xlim(0, T_REC)
    if drift_um is not None:
        k = int(np.argmin(np.abs(drift_depths - units.at[c, "y"])))
        ax2 = ax.twinx()
        ax2.plot(drift_t, drift_um[:, k], color=CB_COLOR, lw=1)
        ax2.set_ylabel(f"drift (µm) @ {drift_depths[k]:.0f} µm", color=CB_COLOR)
        ax2.spines["right"].set_visible(True)
    plot_probe(fig.add_subplot(gs[1:, 3]), c)
    u = units.loc[c]
    fig.suptitle(f"CS phyID = {c} | C4 {u.c4_label} conf {u.c4_conf:.1f} votes {u.c4_votes} | "
                 f"span n = {len(st)} (reference {u.n_spikes}), FR = {u.fr:.2f} sp/s | manual CS = {u.manual_CS}",
                 fontsize=11)
    fig.savefig(FIG_DIR / f"CS{c}_overview.png", dpi=110)


_t = time.time()
for c in CS_IDS:
    plot_cs_overview(c)
print(f"saved {len(CS_IDS)} CS overview figures to {FIG_DIR} ({time.time() - _t:.0f} s)")


# %% 5. CS→SS CANDIDATES (every CS × every SS candidate; same CCG for every score)
def peak_dist(a, b):
    return float(np.hypot(units.at[a, "x"] - units.at[b, "x"], units.at[a, "y"] - units.at[b, "y"]))


CCG_COARSE, CCG_FINE = {}, {}
rows = []
_t = time.time()
PAIR_WIN = {}
ACG_REFR, LAT_HIST = {}, {}                      # v8: SS autocorrelogram per (SS, pair window); latency histograms
LAT_EDGES = np.arange(0, 100.1e-3, 2e-3)


def _acg_refr(ss, w, s_):
    """SS|SS on the same SS spike set as the pair's CCG targets (pair window ± pad), cached per window."""
    k = (ss, round(w[0], 6), round(w[1], 6))
    if k not in ACG_REFR:
        ACG_REFR[k] = ccg(s_, s_, COARSE_BIN, REFR_WIN, acg=True)[:2]
    return ACG_REFR[k]


def lit_row(cs, ss, c_, s_, w, rc, lc, m):
    """v8 pair columns (section 1b). Same pair window and spike sets as the v7 metrics."""
    out, post, pre = lit_pair(c_, s_, rc, lc, _acg_refr(ss, w, s_))
    LAT_HIST[(cs, ss)] = (np.histogram(post, LAT_EDGES)[0], np.histogram(pre, LAT_EDGES)[0])
    ssw = s_[(s_ >= w[0]) & (s_ < w[1])]
    isi = float(np.median(np.diff(ssw))) * 1e3 if len(ssw) > 2 else np.nan
    win_s = max(w[1] - w[0], 1e-9)
    out.update(ss_isi_med_ms=isi,
               pause_isi_units=m["pause_dur50_ms"] / isi if isi and np.isfinite(isi) else np.nan,
               lpost_isi_units=out["lat_post_med_ms"] / isi if isi and np.isfinite(isi) else np.nan,
               dy_um=float(units.at[cs, "y"] - units.at[ss, "y"]), dx_um=float(units.at[cs, "x"] - units.at[ss, "x"]),
               lussac_code_main_fr_ok=bool(len(ssw) / win_s >= LUSSAC_CM_SS_MIN
                                           and LUSSAC_CM_CS_RANGE[0] <= len(c_) / win_s <= LUSSAC_CM_CS_RANGE[1]),
               cs_source=units.at[cs, "cs_source"], cs_c4=bool(units.at[cs, "cs_c4"]),
               cs_manual=bool(units.at[cs, "cs_manual"]), cs_screen=bool(units.at[cs, "cs_screen"]))
    out["dur_rate_ok"] = bool(np.isfinite(out["pause_isi_units"]) and out["pause_isi_units"] > LIT_DUR_ISI_K)
    return out


def cs_ss_row(cs, ss, keep_ccg=True):
    # pair window = overlap of the two units' stable spans; CS refs cut to it, SS targets to it ± pad
    w = PAIR_WIN[(cs, ss)] = pair_win(cs, ss)
    c_, s_ = win_st(cs, w), win_st(ss, w, EDGE_PAD)
    m, (rc, lc), (rf, lf) = cs_ss_pair(c_, s_)
    if keep_ccg:
        CCG_COARSE[(cs, ss)], CCG_FINE[(cs, ss)] = rc, rf
    bl = block_pause_ratios(c_, s_, w[0], w[1])
    _lab = units.at[ss, "c4_label"]
    m.update(cs=cs, ss=ss, ss_tier=2 if units.at[ss, "role"] == "SS_cand2" else 1,
             ss_group=units.at[ss, "group"],
             ss_c4=_lab if isinstance(_lab, str) else "unlabelled", ss_fr=units.at[ss, "fr"],
             dist_um=peak_dist(cs, ss), n_cs=len(c_), n_ss=int(np.sum((s_ >= w[0]) & (s_ < w[1]))),
             win_start=w[0], win_end=w[1], win_s=w[1] - w[0],
             block_pr_min=np.nanmin(bl), block_pr_median=np.nanmedian(bl), block_pr_max=np.nanmax(bl),
             block_pause_ratios=";".join(f"{v:.2f}" for v in bl))
    m.update(cs_spike_sharing(c_, s_, w[0], w[1]))
    ev, _ = cs_event_capture(c_, s_, w[0], w[1])
    m.update(cs_event_frac=ev["frac"], cs_event_null=ev["null"], cs_event_excess=ev["excess"],
             ss_cs_locked_frac=ev["k_locked_frac"])
    m.update(lit_row(cs, ss, c_, s_, w, rc, lc, m))                             # v8
    # same metrics in the reference window (comparison only; jitter only where a candidate is still possible)
    mr = cs_ss_pair(ST_REF[cs], ST_PAD_REF[ss], with_jitter=False)[0]
    maybe = (mr["pause_ratio"] < CAND_PAUSE_RATIO_MAX) and (mr["asym"] > CAND_ASYM_MIN)
    jr = jitter_pause_test(ST_REF[cs], ST_PAD_REF[ss]) if maybe and len(ST_REF[cs]) else {"jitter_p": np.nan}
    m.update(pause_ratio_ref=mr["pause_ratio"], asym_ref=mr["asym"], jitter_p_ref=jr["jitter_p"],
             collision_dip_ref=mr["collision_dip"], n_cs_ref=len(ST_REF[cs]))
    return m


for cs in CS_IDS:
    for ss in SS_IDS:
        rows.append(cs_ss_row(cs, ss))
# v6 audit: CS × extra-good units below the tier-1 FR floor (never candidates; same metrics, nothing else uses them)
lowfr_pairs = pd.DataFrame([cs_ss_row(cs, ss, keep_ccg=False) for cs in CS_IDS for ss in LOWFR_IDS])
# lag axes straight from ccg (identical bins), so they exist even when no CS×SS pair was computed
LAGS_COARSE = ccg(np.empty(0), np.empty(0), COARSE_BIN, COARSE_WIN)[1]
LAGS_FINE = ccg(np.empty(0), np.empty(0), FINE_BIN, FINE_WIN)[1]
print(f"{len(rows)} CS×SS pairs in {time.time() - _t:.0f} s")

if rows:
    pairs = pd.DataFrame(rows)
else:   # no CS or no SS candidate: empty table with the same columns (metric keys from a dummy call)
    _d = np.arange(50.0)
    _keys = (list(cs_ss_pair(_d, _d + 0.5)[0].keys()) + list(cs_spike_sharing(_d, _d + 0.5, 0.0, 50.0).keys())
             + ["cs_event_frac", "cs_event_null", "cs_event_excess", "ss_cs_locked_frac"]
             + ["cs", "ss", "ss_tier", "ss_group", "ss_c4", "ss_fr", "dist_um", "n_cs", "n_ss", "win_start", "win_end",
                "win_s", "block_pr_min", "block_pr_median", "block_pr_max", "block_pause_ratios", "pause_ratio_ref",
                "asym_ref", "jitter_p_ref", "collision_dip_ref", "n_cs_ref"]
             + list(lit_pair(_d, _d + 0.5, *ccg(_d, _d + 0.5, COARSE_BIN, COARSE_WIN)[:2])[0].keys())
             + ["ss_isi_med_ms", "pause_isi_units", "lpost_isi_units", "dy_um", "dx_um", "lussac_code_main_fr_ok",
                "cs_source", "cs_c4", "cs_manual", "cs_screen", "dur_rate_ok"])
    _bool = ("lussac_pass", "collision_dip", "collision_dip_ref", "lussac_paper_v1_pass", "lussac_code_main_pass",
             "lussac_code_main_fr_ok", "cs_c4", "cs_manual", "cs_screen", "dur_rate_ok")
    pairs = pd.DataFrame({k: pd.Series(dtype=object if k in ("ss_group", "ss_c4", "block_pause_ratios", "cs_source") else
                                       bool if k in _bool else float)
                          for k in _keys})
    print("  (no CS×SS pairs: CS or SS candidate list is empty)")
pairs["candidate"] = ((pairs["pause_ratio"] < CAND_PAUSE_RATIO_MAX) & (pairs["jitter_p"] < CAND_JITTER_P_MAX)
                      & (pairs["asym"] > CAND_ASYM_MIN) & ~pairs["collision_dip"])
pairs["candidate_ref"] = ((pairs["pause_ratio_ref"] < CAND_PAUSE_RATIO_MAX) & (pairs["jitter_p_ref"] < CAND_JITTER_P_MAX)
                          & (pairs["asym_ref"] > CAND_ASYM_MIN) & ~pairs["collision_dip_ref"].astype(bool))
pairs["candidate_lit"] = lit_candidate(pairs).astype(bool)                  # v8 Rule 1; does not touch candidate/rank
pairs["far_flag"] = pairs["dist_um"] > FAR_UM                               # v8 flag only
pairs["cs_share_flag"] = pairs["cs_share_excess"] > CS_SHARE_FLAG          # flag only, candidate unchanged
pairs["cs_event_flag"] = pairs["cs_event_excess"] > CS_EVENT_FLAG          # v6, flag only, candidate unchanged
if len(lowfr_pairs):
    lowfr_pairs["would_be_candidate"] = ((lowfr_pairs["pause_ratio"] < CAND_PAUSE_RATIO_MAX)
                                         & (lowfr_pairs["jitter_p"] < CAND_JITTER_P_MAX)
                                         & (lowfr_pairs["asym"] > CAND_ASYM_MIN) & ~lowfr_pairs["collision_dip"])
    _lw = lowfr_pairs[lowfr_pairs["would_be_candidate"]]
    print(f"--- tier-1 FR floor: {len(_lw)} CS×SS pairs would have been candidates without the floor")
    if len(_lw):
        print(_lw[["cs", "ss", "ss_fr", "pause_ratio", "pause_ratio_0_8", "jitter_p", "asym", "dist_um",
                   "cs_event_excess"]].round(3).to_string(index=False))
_sh = pairs[pairs["cs_share_flag"]]
print(f"--- CS-spike sharing: {len(_sh)} CS×SS pairs where > {CS_SHARE_FLAG:.0%} of CS spikes (beyond chance) also "
      f"appear in the SS unit at zero lag ({int(_sh['candidate'].sum())} of them are candidates)")
if len(_sh):
    print(_sh[["cs", "ss", "ss_tier", "ss_group", "ss_c4", "dist_um", "cs_share_excess", "ss_from_cs_frac",
               "pause_ratio", "candidate"]].round(3).to_string(index=False))
_chg = pairs[pairs["candidate"] != pairs["candidate_ref"]]
print(f"--- candidates: {int(pairs['candidate'].sum())} with stable spans vs {int(pairs['candidate_ref'].sum())} in the "
      f"reference window; changed: {len(_chg)}")
if len(_chg):
    print(_chg[["cs", "ss", "candidate", "candidate_ref", "pause_ratio", "pause_ratio_ref", "jitter_p", "jitter_p_ref",
                "asym", "asym_ref", "n_cs", "n_cs_ref", "win_start", "win_end"]].round(3).to_string(index=False))
pairs = pairs.sort_values(["cs", "candidate", "jitter_z", "pause_ratio"], ascending=[True, False, True, True])
pairs["rank"] = pairs.groupby("cs").cumcount() + 1
front = ["cs", "ss", "rank", "candidate", "ss_tier", "ss_group", "ss_c4", "pause_ratio", "mirror_ratio", "asym",
         "pause_dur50_ms",
         "pause_dur80_ms", "jitter_p", "jitter_z", "jitter_deficit", "lussac_pass", "collision_dip", "zero_peak_ratio",
         "cs_share_excess", "ss_from_cs_frac", "cs_share_flag", "cs_event_excess", "ss_cs_locked_frac", "cs_event_flag",
         "pause_ratio_0_8", "asym_0_8",
         "block_pr_min", "block_pr_median", "block_pr_max", "dist_um"]
pairs = pairs[front + [c for c in pairs.columns if c not in front]].reset_index(drop=True)

best = pairs[pairs["rank"] == 1].set_index("cs")
best["n_candidates"] = pairs.groupby("cs")["candidate"].sum()
best["n_lussac"] = pairs.groupby("cs")["lussac_pass"].sum()
with pd.option_context("display.width", 220, "display.max_rows", 200):
    print("--- best-ranked SS per CS")
    print(best[["ss", "ss_tier", "ss_c4", "candidate", "n_candidates", "n_lussac", "pause_ratio", "asym", "pause_dur50_ms",
                "jitter_z", "jitter_p", "jitter_deficit", "block_pr_min", "block_pr_max", "zero_peak_ratio",
                "dist_um"]].round(3).to_string())
_flat = pairs["pause_ratio"] > 0.9
if _flat.any():
    print(f"jitter null check: of {_flat.sum()} pairs with pause_ratio > 0.9, "
          f"{(pairs.loc[_flat, 'jitter_p'] < CAND_JITTER_P_MAX).mean():.0%} have jitter_p < {CAND_JITTER_P_MAX}")


# %% 5b. CS EVENT DECOMPOSITION (v6; spike times only; nothing merged or excluded)
# Which clusters carry a piece of each CS event (EVENT_WIN after the CS spike), at which lag, and how much of their
# own firing is CS-locked. Lag ~0 = zero-lag copy (re-detection of the main spike, see section 6); lag > 0.3 ms =
# later piece (spikelet / late component). "mostly CS pieces" = most of the cluster's spikes sit in CS events.
_t = time.time()
_rank1 = pairs[pairs["rank"] == 1].set_index("cs")["ss"].to_dict() if len(pairs) else {}
_cand = pairs[pairs["candidate"]].groupby("cs")["ss"].apply(set).to_dict() if len(pairs) else {}
dec_rows = []
for cs in CS_IDS:
    for k in units.index:
        if k == cs or k not in ST_FULL or peak_dist(cs, k) > EVENT_NEIGH_UM:
            continue
        w = pair_win(cs, k)
        c_, k_ = win_st(cs, w), win_st(k, w, EDGE_PAD)
        if len(c_) < EVENT_MIN_N:
            continue
        ev, lags = cs_event_capture(c_, k_, w[0], w[1])
        if ev["n_hit"] < EVENT_MIN_N or not ev["excess"] >= EVENT_MIN_EXCESS:
            continue
        lag_med = float(np.median(lags))
        piece = ("zero-lag copy" if abs(lag_med) <= EVENT_ZERO_LAG else
                 "late piece" if lag_med > 0 else "earlier piece")
        _lab = units.at[k, "c4_label"]
        dec_rows.append(dict(
            cs=cs, unit=k, role=units.at[k, "role"] or "other", group=units.at[k, "group"],
            c4=_lab if isinstance(_lab, str) else "unlabelled", fr=units.at[k, "fr"], dist_um=peak_dist(cs, k),
            n_cs=len(c_), n_hit=ev["n_hit"], frac_cs_events=ev["frac"], chance=ev["null"], excess=ev["excess"],
            unit_locked_frac=ev["k_locked_frac"], lag_med_ms=lag_med * 1e3,
            lag_q25_ms=float(np.quantile(lags, .25)) * 1e3, lag_q75_ms=float(np.quantile(lags, .75)) * 1e3,
            piece=piece, mostly_cs_pieces=ev["k_locked_frac"] >= EVENT_LOCKED_MOST,
            strong=(ev["excess"] > CS_EVENT_FLAG) or (ev["k_locked_frac"] >= EVENT_LOCKED_MOST),
            is_cs_unit=k in CS_IDS, is_ss_candidate=k in _cand.get(cs, set()), is_rank1_ss=_rank1.get(cs) == k,
            win_start=w[0], win_end=w[1]))
_dec_cols = ["cs", "unit", "role", "group", "c4", "fr", "dist_um", "n_cs", "n_hit", "frac_cs_events", "chance", "excess",
             "unit_locked_frac", "lag_med_ms", "lag_q25_ms", "lag_q75_ms", "piece", "mostly_cs_pieces", "strong", "is_cs_unit",
             "is_ss_candidate", "is_rank1_ss", "win_start", "win_end"]
decomp = (pd.DataFrame(dec_rows, columns=_dec_cols).sort_values(["cs", "excess"], ascending=[True, False])
          .reset_index(drop=True))
_n_piece = decomp.groupby("cs").size().reindex(CS_IDS, fill_value=0)
print(f"--- CS event decomposition ({time.time() - _t:.0f} s): {len(decomp)} CS×cluster pieces "
      f"(excess >= {EVENT_MIN_EXCESS}, >= {EVENT_MIN_N} events, <= {EVENT_NEIGH_UM:.0f} µm); "
      f"CS with >= 1 piece: {int((_n_piece > 0).sum())} / {len(CS_IDS)}; strong (excess > {CS_EVENT_FLAG} or >= "
      f"{EVENT_LOCKED_MOST:.0%} of the cluster CS-locked): {int(decomp['strong'].sum())}")
_ssp = decomp[decomp["is_ss_candidate"]]
print(f"  SS candidates that also carry pieces of their CS's event: {len(_ssp)} "
      f"(rank-1 SS: {int(_ssp['is_rank1_ss'].sum())})")
with pd.option_context("display.width", 220, "display.max_rows", 300):
    if len(decomp):
        print(decomp.loc[decomp["strong"] | decomp["is_ss_candidate"],
                         ["cs", "unit", "role", "dist_um", "excess", "unit_locked_frac", "lag_med_ms", "piece",
                          "is_ss_candidate", "is_rank1_ss"]].round(3).to_string(index=False))


# %% 6. CS–CS DUPLICATE CANDIDATES
sim_t = np.load(sorter_file("similar_templates.npy")) if sorter_file("similar_templates.npy", False) else None
dom_tmpl = {c: int(np.argmax(_w[c])) for c in CS_IDS}


def template_corr(a, b, radius=60.0):
    pa, pb = chan_pos[int(units.at[a, "peak_ch"])], chan_pos[int(units.at[b, "peak_ch"])]
    chs = np.nonzero((np.hypot(*(chan_pos - pa).T) <= radius) | (np.hypot(*(chan_pos - pb).T) <= radius))[0]
    return float(np.corrcoef(TEMPLATE[a][:, chs].ravel(), TEMPLATE[b][:, chs].ravel())[0, 1])


rows, CCG_DUP = [], {}
_t = time.time()
for i, a in enumerate(CS_IDS):
    for b in CS_IDS[i + 1:]:
        w = pair_win(a, b)
        ta_, tb_ = win_st(a, w), win_st(b, w)
        m, (rf, lf) = cs_cs_pair(ta_, tb_)
        CCG_DUP[(a, b)] = rf
        mr = cs_cs_pair(ST_REF[a], ST_REF[b])[0]
        m.update(cs_a=a, cs_b=b, n_a=len(ta_), n_b=len(tb_), dist_um=peak_dist(a, b),
                 win_start=w[0], win_end=w[1], coinc_frac_peak_ref=mr["coinc_frac_peak"],
                 peak_lag_ms_ref=mr["peak_lag_ms"],
                 template_corr=template_corr(a, b),
                 similar_templates=float(sim_t[dom_tmpl[a], dom_tmpl[b]]) if sim_t is not None else np.nan)
        rows.append(m)
if rows:
    dups = pd.DataFrame(rows)
else:   # < 2 CS: empty table with the same columns (keys of cs_cs_pair taken from a dummy call)
    _keys = list(cs_cs_pair(np.arange(50.0), np.arange(50.0) + 0.5)[0].keys())
    dups = pd.DataFrame(columns=_keys + ["cs_a", "cs_b", "n_a", "n_b", "dist_um", "win_start", "win_end",
                                         "coinc_frac_peak_ref", "peak_lag_ms_ref", "template_corr",
                                         "similar_templates"]).astype(float)
dups["dup_candidate"] = dups["coinc_frac_peak"] > DUP_COINC_MIN
dups["sharp_peak_only"] = ~dups["dup_candidate"] & (dups["peak_p"] < DUP_PEAK_P_MAX)  # e.g. shared artefacts
front = ["cs_a", "cs_b", "dup_candidate", "sharp_peak_only", "coinc_frac_peak", "coinc_frac_zero", "peak_lag_ms", "peak_fwhm_ms",
         "peak_p", "peak_count", "union_new_short_isi", "union_new_short_frac", "template_corr",
         "similar_templates", "dist_um"]
dups = dups[front + [c for c in dups.columns if c not in front]].sort_values("coinc_frac_peak", ascending=False)
dups = dups.reset_index(drop=True)
print(f"{len(dups)} CS–CS pairs in {time.time() - _t:.0f} s; {dups['dup_candidate'].sum()} duplicate candidates; "
      f"{dups['sharp_peak_only'].sum()} sharp-peak-only pairs (low coincidence, listed in the table, not flagged)")
with pd.option_context("display.width", 220, "display.max_rows", 200):
    print(dups[dups["dup_candidate"]][front].round(3).to_string())


# %% 6b. DUPLICATE EVIDENCE (no merging: independent evidence lines + a descriptive type per pair)
# Evidence lines (each votes +1 same cell / -1 distinct cells / 0 uninformative or not testable):
#   coinc      fraction of the smaller unit's spikes coincident at the peak lag
#   lag_sd     jitter of the coincident lag (robust SD): duplicates are sub-0.1 ms, CF synchrony is ~ms
#   union      new ISIs in UNION_ISI_RANGE after merging the trains (two cells create them; spikelets are < 4 ms)
#   ss_shared  the two units pause the same SS candidate (one-to-one: same SS → same cell)
#   ss_cond    co-silencing: spikes of one unit that the OTHER missed still pause the other's SS partner
#   footprint  template ptp-footprint cosine, only when the two footprints overlap (else not testable)
#   covar      detrended amplitude or depth co-variation of coincident spikes (one event → co-vary)
# Descriptive types (several may apply): same_event_redetection, fixed_lag_component, partial_split,
# synchrony_like. They describe the evidence pattern; they are not decisions.

SS_OF = pairs[pairs["candidate"]].groupby("cs")["ss"].apply(set).to_dict()


def _per_spike(arr, c, w=None):
    """Per-spike values of unit c, aligned with win_st(c, w) (or ST[c] when w is None)."""
    if arr is None:
        return None
    idx = SPK_IDX[c] if w is None else SPK_IDX_FULL[c][win_idx(c, w)]
    return np.asarray(arr[idx], float).ravel()


def dup_evidence(a, b, lag_ms, coinc_frac_peak):
    """Session wrapper around dup_evidence_core (lag = time of b relative to a at the fine-CCG peak).
    Both trains are cut to the pair window; SS-partner pauses only use reference spikes inside that SS's span."""
    w = pair_win(a, b)
    return dup_evidence_core(win_st(a, w), win_st(b, w), lag_ms * 1e-3, coinc_frac_peak,
                             _per_spike(amplitudes, a, w), _per_spike(amplitudes, b, w),
                             _per_spike(SPIKE_DEPTH, a, w), _per_spike(SPIKE_DEPTH, b, w),
                             np.ptp(TEMPLATE[a], 0), np.ptp(TEMPLATE[b], 0),
                             SS_OF.get(a, set()), SS_OF.get(b, set()), ST_PAD, w[0], w[1], ss_spans=SPAN)


ev_rows = []
for _, r in dups[dups["dup_candidate"] | dups["sharp_peak_only"]].iterrows():
    a, b = int(r.cs_a), int(r.cs_b)
    e = dict(cs_a=a, cs_b=b, dup_candidate=bool(r.dup_candidate),
             peak_lag_ms=r.peak_lag_ms, peak_fwhm_ms=r.peak_fwhm_ms, dist_um=r.dist_um, n_a=r.n_a, n_b=r.n_b,
             win_start=r.win_start, win_end=r.win_end, coinc_frac_peak_ref=r.coinc_frac_peak_ref)
    e.update(dup_evidence(a, b, r.peak_lag_ms, r.coinc_frac_peak)[0])
    e["type_hint"] = classify_pair(e, r.peak_lag_ms, r.dup_candidate)[0]
    e["evidence_status"] = evidence_status(e)
    ev_rows.append(e)
dup_ev = pd.DataFrame(ev_rows)

# duplicate families: connected components of dup-candidate pairs (descriptive grouping, not a merge)
_parent = {c: c for c in CS_IDS}
def _find(c):
    while _parent[c] != c:
        c = _parent[c]
    return c
for _, r in (dup_ev[dup_ev["dup_candidate"]].iterrows() if len(dup_ev) else []):
    _parent[_find(int(r.cs_a))] = _find(int(r.cs_b))
_fam = {}
for c in CS_IDS:
    _fam.setdefault(_find(c), []).append(c)
FAMILIES = {min(m): sorted(m) for m in _fam.values() if len(m) > 1}
# v8: families among C4 CS only (= v7's definition; written to cs_dup_families.json so downstream counts keep their
#     meaning); FAMILIES (all CS incl. manual / screen) → cs_dup_families_all.json
_C4 = set(units.index[units["cs_c4"]])
_p4 = {c: c for c in CS_IDS if c in _C4}


def _find4(c):
    while _p4[c] != c:
        c = _p4[c]
    return c


for _, r in (dup_ev[dup_ev["dup_candidate"]].iterrows() if len(dup_ev) else []):
    if int(r.cs_a) in _p4 and int(r.cs_b) in _p4:
        _p4[_find4(int(r.cs_a))] = _find4(int(r.cs_b))
_f4 = {}
for c in _p4:
    _f4.setdefault(_find4(c), []).append(c)
FAMILIES_C4 = {min(m): sorted(m) for m in _f4.values() if len(m) > 1}
if len(dup_ev):
    dup_ev["family"] = [min(_fam[_find(int(a))]) if len(_fam[_find(int(a))]) > 1 else np.nan for a in dup_ev.cs_a]
    front_ev = ["family", "cs_a", "cs_b", "dup_candidate", "type_hint", "n_support", "n_against", "support", "against",
                "coinc_frac_peak", "peak_lag_ms", "lag_sd_ms", "union_new_isi_frac", "ss_a", "ss_b", "ss_shared",
                "ss_cond_pr_min", "footprint_cos", "amp_corr_detr", "depth_corr_detr", "rate_corr_60s", "dist_um"]
    dup_ev = dup_ev[front_ev + [c for c in dup_ev.columns if c not in front_ev]]
    with pd.option_context("display.width", 260, "display.max_rows", 200, "display.max_colwidth", 40):
        print("--- duplicate evidence (dup candidates + sharp-peak-only pairs); NOTHING IS MERGED")
        print(dup_ev[front_ev].round(3).to_string(index=False))
print("--- duplicate families (connected dup candidates):", FAMILIES if FAMILIES else "none")


# %% 7. CONFLICT VIEW (no resolution)
cand = pairs[pairs["candidate"]]
ss_of = cand.groupby("cs")["ss"].apply(lambda s: sorted(s)).to_dict()
cs_of = cand.groupby("ss")["cs"].apply(lambda s: sorted(s)).to_dict()
conf = []
for ss, css in cs_of.items():
    if len(css) > 1:
        conf.append(dict(type="SS_claimed_by_multiple_CS", unit=ss, members=css, detail=""))
for cs, sss in ss_of.items():
    if len(sss) > 1:
        conf.append(dict(type="CS_with_multiple_SS", unit=cs, members=sss, detail=""))
for _, r in dups[dups["dup_candidate"]].iterrows():
    sa, sb = set(ss_of.get(r.cs_a, [])), set(ss_of.get(r.cs_b, []))
    if sa & sb:
        tag = "same_SS"
    elif sa and sb:
        tag = "different_SS"
    elif sa or sb:
        tag = "one_SS"
    else:
        tag = "no_SS"
    conf.append(dict(type=f"CS_duplicate_candidate:{tag}", unit=int(r.cs_a), members=[int(r.cs_b)],
                     detail=f"SS(a)={sorted(sa)} SS(b)={sorted(sb)} coinc={r.coinc_frac_peak:.2f} "
                            f"lag={r.peak_lag_ms:.2f}ms tcorr={r.template_corr:.2f} d={r.dist_um:.0f}um"))
conflicts = pd.DataFrame(conf, columns=["type", "unit", "members", "detail"])
with pd.option_context("display.width", 220, "display.max_rows", 300, "display.max_colwidth", 120):
    print("--- conflicts (permissive candidates; input to the future adjudication step)")
    print(conflicts.to_string(index=False) if len(conflicts) else "none")


# %% 7b. v8 REPORTS (flags and comparisons; nothing here changes a v7 column)
cs_tab = units.loc[CS_IDS, ["cs_source", "cs_c4", "cs_manual", "cs_screen", "c4_label", "extragood", "fr", "scr_fr",
                            "cs_fr_ok", "isi_5_35_frac", "isi_5_35_ok", "cs_isi_10_35_frac", "isi_lt50_frac", "cv", "cv2",
                            "presence_bins_frac", "span_frac", "partial_presence", "poor_isolation", "frac_isi_lt1ms",
                            "tmpl_dur_ms", "tmpl_n_fp_ch", "tmpl_n_late_peaks"]].copy()
cs_tab["n_candidate"] = pairs.groupby("cs")["candidate"].sum().reindex(CS_IDS).fillna(0).astype(int)
cs_tab["n_candidate_lit"] = pairs.groupby("cs")["candidate_lit"].sum().reindex(CS_IDS).fillna(0).astype(int)
cs_tab["paired"] = (cs_tab["n_candidate"] > 0) | (cs_tab["n_candidate_lit"] > 0)
# duplicate / piece of a C4 CS? (family with a C4 CS, or a strong event-decomposition row with a C4 CS, either way)
_fam_of = {c: m for m in FAMILIES.values() for c in m}
_strong = decomp[decomp["strong"]] if len(decomp) else decomp


def _dup_of_c4(c):
    hits = {x for x in _fam_of.get(c, []) if x != c and x in _C4}
    if len(_strong):
        hits |= set(_strong.loc[(_strong.unit == c) & _strong.cs.isin(_C4), "cs"])
        hits |= set(_strong.loc[(_strong.cs == c) & _strong.unit.isin(_C4), "unit"])
    return sorted(int(h) for h in hits)


cs_tab["dup_or_piece_of_c4"] = [";".join(map(str, _dup_of_c4(c))) if not units.at[c, "cs_c4"] else "" for c in CS_IDS]
cs_tab["recovered_cs"] = ~cs_tab["cs_c4"] & cs_tab["paired"] & (cs_tab["dup_or_piece_of_c4"] == "")
cs_tab["fails"] = [";".join([k for k in ("cs_fr_ok", "isi_5_35_ok") if not bool(cs_tab.at[c, k])]
                            + (["isi_lt50>0"] if cs_tab.at[c, "isi_lt50_frac"] > 0 else [])) for c in CS_IDS]
# screen block: the pool (C4-unlabelled) + calibration on labelled extra-good CB units (columns only)
_scr_cols = ["c4_label", "extragood", "manual_CS", "role", "cs_source", "scr_fr", "cs_fr_ok", "isi_5_35_frac", "isi_5_35_ok",
             "cs_isi_10_35_frac", "isi_lt50_frac", "cv", "cv2", "screen_pass", "cs_screen", "partial_presence",
             "poor_isolation", "tmpl_dur_ms", "tmpl_n_fp_ch", "tmpl_n_late_peaks", "c4_soft_available"]
screen_tab = units.loc[list(dict.fromkeys(SCREEN_POOL + CALIB_POOL)), _scr_cols].copy()
screen_tab["pool"] = np.where(screen_tab.index.isin(SCREEN_POOL), "screen (C4-unlabelled)", "calibration (C4-labelled)")
calib = (screen_tab[screen_tab.pool.str.startswith("calib")].assign(c4=lambda d: d.c4_label.fillna("?"))
         .groupby("c4")["screen_pass"].agg(n="size", pass_share="mean"))
_agree = units.loc[sorted(set(SCREEN_POOL) | set(MANUAL_IDS))]
agree_tab = (pd.crosstab(_agree["manual_CS"].astype(bool).rename("manual flag"),
                         _agree["screen_pass"].astype(bool).rename("screen pass")) if len(_agree) else pd.DataFrame())
# manual / screen units kept as SS (option A): FR, C4 label, tier, best pause (as SS)
kept_rows = []
for c in KEPT_AS_SS:
    q = pairs[pairs.ss == c].sort_values("pause_ratio")
    b = q.iloc[0] if len(q) else None
    kept_rows.append(dict(unit=c, role=units.at[c, "role"], fr=units.at[c, "fr"], c4_label=units.at[c, "c4_label"],
                          extragood=bool(units.at[c, "extragood"]), manual=bool(units.at[c, "manual_CS"]),
                          screen_pass=bool(units.at[c, "screen_pass"]),
                          best_cs=None if b is None else int(b.cs), best_pause_ratio=None if b is None else b.pause_ratio,
                          best_candidate=None if b is None else bool(b.candidate),
                          best_candidate_lit=None if b is None else bool(b.candidate_lit),
                          n_candidate_as_ss=int(q.candidate.sum()) if len(q) else 0))
kept_tab = pd.DataFrame(kept_rows, columns=["unit", "role", "fr", "c4_label", "extragood", "manual", "screen_pass", "best_cs",
                                            "best_pause_ratio", "best_candidate", "best_candidate_lit", "n_candidate_as_ss"])
# spatial check: distances of all CS×SS pairs (chance) vs rank 1 / candidate / candidate_lit; depth-offset sign
_grp = {"all pairs (chance)": pairs, "rank 1": pairs[pairs["rank"] == 1], "candidate": pairs[pairs.candidate],
        "candidate_lit": pairs[pairs.candidate_lit], "candidate & C4 CS": pairs[pairs.candidate & pairs.cs_c4]}
dist_rows = []
for g, d in _grp.items():
    x = d["dist_um"]
    strong = d[d.candidate & (d.pause_ratio < STRONG_PAUSE)] if g != "all pairs (chance)" else d.iloc[:0]
    dist_rows.append(dict(group=g, n=len(d), med_um=x.median(), q25_um=x.quantile(.25), q75_um=x.quantile(.75),
                          le50=(x <= 50).mean(), le100=(x <= 100).mean(), le200=(x <= 200).mean(), gt500=(x > FP_UM).mean(),
                          dy_pos_share=(d["dy_um"] > 0).mean() if len(d) else np.nan,
                          n_strong=len(strong), strong_dy_pos_share=(strong["dy_um"] > 0).mean() if len(strong) else np.nan))
_far = pairs[pairs.dist_um > FP_UM]
for k in ("candidate", "candidate_lit", "lussac_pass", "lussac_paper_v1_pass", "lussac_code_main_pass"):
    dist_rows.append(dict(group=f"> {FP_UM:.0f} µm passing {k} (false-positive rate)", n=len(_far),
                          gt500=float(_far[k].mean()) if len(_far) else np.nan))
dist_tab = pd.DataFrame(dist_rows)
with pd.option_context("display.width", 240, "display.max_rows", 300, "display.max_columns", 40):
    print("=" * 30, "v8 REPORTS", "=" * 30)
    print("--- CS by source")
    print(cs_tab.round(3).to_string())
    _fail = cs_tab[cs_tab.fails != ""]
    print(f"--- CS unit-check failures: {len(_fail)} of {len(cs_tab)} CS; of these paired: {int(_fail.paired.sum())} | "
          f"CS passing every check that are paired: {int(cs_tab[cs_tab.fails == ''].paired.sum())} of "
          f"{int((cs_tab.fails == '').sum())}")
    print(f"--- SCREEN (pool: {len(SCREEN_POOL)} C4-unlabelled {'' if SCREEN_NON_EG else 'extra-good '}CB units; entry = "
          f"FR {CS_FR_RANGE[0]}–{CS_FR_RANGE[1]} Hz AND ISIs {CS_ISI_5_35[0] * 1e3:.0f}–{CS_ISI_5_35[1] * 1e3:.0f} ms < "
          f"{CS_ISI_5_35_MAX:.0%}): pass {int(units.loc[SCREEN_POOL, 'screen_pass'].sum())} | "
          f"fr_ok {int(units.loc[SCREEN_POOL, 'cs_fr_ok'].fillna(False).sum())} | "
          f"isi_ok {int(units.loc[SCREEN_POOL, 'isi_5_35_ok'].fillna(False).sum())}")
    _sp_ = screen_tab[screen_tab.pool.str.startswith("screen") & (screen_tab.screen_pass | screen_tab.manual_CS.astype(bool)
                                                                   | screen_tab.cs_fr_ok.fillna(False).astype(bool))]
    if len(_sp_):
        print(_sp_.round(3).to_string())
    _scs = cs_tab[cs_tab.cs_screen]
    print(f"  screen CS: {len(_scs)} | paired (candidate / candidate_lit): {int((_scs.n_candidate > 0).sum())} / "
          f"{int((_scs.n_candidate_lit > 0).sum())} | dup/piece of a C4 CS: {int((_scs.dup_or_piece_of_c4 != '').sum())} | "
          f"also manual: {int(_scs.cs_manual.sum())}")
    print("--- screen calibration on C4-labelled extra-good CB units (share passing; PkC_cs = sensitivity, others = false alarms)")
    print(calib.round(3).to_string() if len(calib) else "  none")
    print("--- manual flag vs screen pass (screen pool ∪ manual units)")
    print(agree_tab.to_string() if len(agree_tab) else "  none")
    _man = cs_tab[cs_tab.cs_manual & ~cs_tab.cs_c4]
    print(f"--- manual-only CS (file {MANUAL_FILE}): {len(_man)} added | paired (candidate / candidate_lit): "
          f"{int((_man.n_candidate > 0).sum())} / {int((_man.n_candidate_lit > 0).sum())} | dup/piece of a C4 CS: "
          f"{int((_man.dup_or_piece_of_c4 != '').sum())}")
    print(f"--- manual flag, kept as SS (option A): {len(kept_tab)}")
    if len(kept_tab):
        print(kept_tab.round(3).to_string(index=False))
    print(f"--- recovered CS (manual or screen, paired, not a duplicate/piece of a C4 CS): "
          f"{cs_tab.index[cs_tab.recovered_cs].tolist()}")
    print("--- spatial check (dy_um > 0: CS peak channel farther from the probe tip than the SS's)")
    print(dist_tab.round(3).to_string(index=False))
    _c, _f = pairs[pairs.candidate], pairs[pairs.dist_um > FP_UM]
    print("--- rate-relative pause: quantiles 5/25/50/75/95 % (candidates vs pairs > 500 µm)")
    for _n, _d in (("candidate", _c), (f"> {FP_UM:.0f} um", _f)):
        if len(_d):
            print(f"  {_n} (n {len(_d)}): " + " | ".join(
                f"{k} " + "/".join(f"{v:.2f}" for v in _d[k].quantile([.05, .25, .5, .75, .95]))
                for k in ("pause_isi_units", "lat_ratio", "lat_cliff", "refr_ratio")))
    print(f"  candidate vs candidate_lit: both {int((pairs.candidate & pairs.candidate_lit).sum())} | candidate only "
          f"{int((pairs.candidate & ~pairs.candidate_lit).sum())} | candidate_lit only "
          f"{int((~pairs.candidate & pairs.candidate_lit).sum())}")


# %% 8. INTERACTIVE BROWSER (re-runnable; nothing is merged anywhere)
# One CS at a time (Prev / Next, or type a phyID in the box). Two views (radio "view"):
#   pairing     funnel of the CS→SS candidate criteria (tier 1 / tier 2 stacked), a pass/fail grid of every
#               criterion for the top-ranked SS, and their CCGs. Click a grid row or a CCG → CS×SS detail figure.
#   duplicates  vote matrix of every CS–CS pair of this CS (evidence lines × partner CS), the decision path of the
#               selected pair, its fine CCG and its co-silencing CCGs. Click a matrix column to select a pair;
#               "pair detail" opens the full evidence figure. The path always ends in manual review.
# Headless (Agg backend): snapshots of both views and of every duplicate-candidate pair are saved instead.
try:
    get_ipython().run_line_magic("matplotlib", "qt")  # noqa: F821
except Exception:
    pass

from matplotlib.colors import ListedColormap
from matplotlib.widgets import TextBox

RANK_KEYS = {"rank": ("rank", True), "jitter z": ("jitter_z", True), "pause ratio": ("pause_ratio", True),
             "distance": ("dist_um", True)}
VOTE_CMAP = ListedColormap([(0.93, 0.6, 0.6), (0.92, 0.92, 0.92), (0.6, 0.83, 0.6)])
TIER_COLORS = {1: CB_COLOR, 2: (0.55, 0.55, 0.55)}
# candidate criteria, in the order of the funnel (all four make pairs["candidate"])
CAND_STEPS = [
    (f"pause ratio < {CAND_PAUSE_RATIO_MAX}", lambda d: d["pause_ratio"] < CAND_PAUSE_RATIO_MAX),
    (f"asym > {CAND_ASYM_MIN}", lambda d: d["asym"] > CAND_ASYM_MIN),
    (f"jitter p < {CAND_JITTER_P_MAX}", lambda d: d["jitter_p"] < CAND_JITTER_P_MAX),
    ("no collision dip", lambda d: ~d["collision_dip"].astype(bool)),
]
# grid columns: (header, value column, pass test, is a candidate criterion); info columns are shown, not used
GRID_COLS = [
    ("pause", "pause_ratio", lambda v: v < CAND_PAUSE_RATIO_MAX, True),
    ("asym", "asym", lambda v: v > CAND_ASYM_MIN, True),
    ("jitter p", "jitter_p", lambda v: v < CAND_JITTER_P_MAX, True),
    ("no dip", "collision_dip", lambda v: not bool(v), True),
    ("blocks\nmax", "block_pr_max", lambda v: v < CAND_PAUSE_RATIO_MAX, False),
    ("dur50\nms", "pause_dur50_ms", lambda v: v >= 3, False),
    ("zero\npeak", "zero_peak_ratio", lambda v: v < 2, False),
    ("Lussac", "lussac_pass", lambda v: bool(v), False),
]
EV_LABELS = {"coinc": "coincidence", "lag_sd": "lag jitter", "union": "union ISIs", "ss_shared": "shared SS",
             "ss_cond": "co-silencing", "footprint": "footprint", "covar": "amp/depth covar"}


def _fmt(v, nd=2):
    if isinstance(v, (bool, np.bool_)):
        return "y" if v else "n"
    try:
        return "–" if not np.isfinite(v) else f"{v:.{nd}f}"
    except TypeError:
        return str(v)


def ev_value_text(r, k):
    """Short value shown in a vote-matrix cell."""
    if k == "coinc":
        return f"{r.coinc_frac_peak:.2f}"
    if k == "lag_sd":
        return f"{_fmt(r.lag_sd_ms)} ms"
    if k == "union":
        return f"{_fmt(r.union_new_ratio)}×"
    if k == "ss_shared":
        return ",".join(map(str, r.ss_shared)) if len(r.ss_shared) else (
            "diff" if len(r.ss_a) and len(r.ss_b) else "–")
    if k == "ss_cond":
        return _fmt(r.ss_cond_pr_min)
    if k == "footprint":
        return _fmt(r.footprint_cos)
    if k == "covar":
        return f"{_fmt(r.amp_corr_detr)}/{_fmt(r.depth_corr_detr)}"
    return ""


def pairs_of_cs(cs):
    """Rows of dup_ev that involve cs, oriented so that cs is unit a (lag = partner relative to cs)."""
    if not len(dup_ev):
        return []
    out = []
    for _, r in dup_ev[(dup_ev.cs_a == cs) | (dup_ev.cs_b == cs)].iterrows():
        out.append((int(r.cs_b) if r.cs_a == cs else int(r.cs_a), r))
    return sorted(out, key=lambda x: -x[1].coinc_frac_peak)


def norm_ccg(refs, ss):
    """CCG of ss after refs, divided by its own baseline (same windows as pause_metrics)."""
    if len(refs) < MIN_SPIKES_COND:
        return None
    rc, lc, _ = ccg(refs, ss, COARSE_BIN, COARSE_WIN)
    b = rc[_win_mask(lc, COARSE_BIN, BASELINE_WIN)].mean()
    return rc / b if b > 0 else None


def plot_cosilencing(ax, r, aux, max_ss=2):
    """SS partner CCGs after coincident events, a-only and b-only spikes (the co-silencing evidence)."""
    ta = aux["ta"]
    sss = sorted(set(r.ss_a) | set(r.ss_b))[:max_ss]
    if not sss:
        ax.text(0.5, 0.5, "no SS partner candidate for either unit:\nco-silencing not testable",
                ha="center", va="center", fontsize=9, transform=ax.transAxes)
        ax.set_axis_off()
        return
    styles = {"coincident": ("k", "-"), f"{int(r.cs_a)}-only": (CB_COLOR, "-"), f"{int(r.cs_b)}-only": ((0.2, 0.5, 0.8), "-")}
    refs = {"coincident": ta[aux["ia"]], f"{int(r.cs_a)}-only": aux["a_only"], f"{int(r.cs_b)}-only": aux["b_only"]}
    shade_windows(ax)
    for i, s in enumerate(sss):
        for lab, ref in refs.items():
            y = norm_ccg(span_cut(ref, s), ST_PAD[s])
            if y is None:
                continue
            c, _ = styles[lab]
            ax.step(LAGS_COARSE * 1e3, y, where="post", color=c, lw=1.2, ls="-" if i == 0 else ":",
                    label=f"SS {s} after {lab} (n={len(ref)}): pr {pause_after(span_cut(ref, s), ST_PAD[s]):.2f}")
    ax.axhline(1, color=(0.6, 0.6, 0.6), lw=0.8, ls="--")
    ax.set_xlabel("time lag (ms)")
    ax.set_ylabel("SS rate / baseline")
    ax.legend(fontsize=6, frameon=False, loc="lower right")
    ax.set_title("co-silencing: SS partner after each spike subset", fontsize=9)


def draw_decision_path(ax, r):
    """The type-hint logic of classify_pair as a tree, with the path of this pair highlighted."""
    ax.clear()
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    _, t = classify_pair(r, r.peak_lag_ms, r.dup_candidate)
    on, off = dict(fc=(0.93, 0.89, 0.97), ec=CB_COLOR, lw=1.8), dict(fc="white", ec=(0.75, 0.75, 0.75), lw=1)
    lsd = _fmt(r.lag_sd_ms)

    def node(x, y, txt, active, fs=7.5):
        ax.text(x, y, txt, ha="center", va="center", fontsize=fs, color="k" if active else (0.55, 0.55, 0.55),
                bbox=dict(boxstyle="round,pad=0.35", **(on if active else off)))

    def edge(p, q, lab, active):
        ax.annotate("", xy=q, xytext=p, arrowprops=dict(arrowstyle="->", lw=1.8 if active else 0.8,
                                                        color=CB_COLOR if active else (0.75, 0.75, 0.75)))
        ax.text((p[0] + q[0]) / 2, (p[1] + q[1]) / 2, lab, fontsize=7, ha="center",
                color=CB_COLOR if active else (0.6, 0.6, 0.6))

    node(0.5, 0.95, f"pair {int(r.cs_a)} / {int(r.cs_b)}: dup candidate = {t['candidate']} "
                    f"(coinc {r.coinc_frac_peak:.2f} at {r.peak_lag_ms:+.2f} ms)", True, fs=8)
    node(0.5, 0.80, f"tight lock?\nlag SD {lsd} ms ≤ {EV_LAG_SD[0] * 1e3:.2f}  and  coinc ≥ {EV_COINC[0]}", True)
    edge((0.5, 0.75), (0.27, 0.62), "yes", t["tight"])
    edge((0.5, 0.75), (0.73, 0.62), "no", not t["tight"])
    node(0.27, 0.59, f"|lag| {abs(r.peak_lag_ms):.2f} ≤ {COINC_DELTA * 1e3:.1f} ms?", t["tight"])
    node(0.73, 0.59, f"lag SD ≥ {EV_LAG_SD[1] * 1e3:.1f} ms\nor different SS partners?", not t["tight"])
    edge((0.27, 0.55), (0.13, 0.44), "yes", t["tight"] and t["zero_lag"])
    edge((0.27, 0.55), (0.40, 0.44), "no", t["tight"] and not t["zero_lag"])
    edge((0.73, 0.55), (0.62, 0.44), "yes", not t["tight"] and t["sync"])
    edge((0.73, 0.55), (0.87, 0.44), "no", not t["tight"] and not t["sync"])
    node(0.13, 0.40, "same-event\nre-detection", t["tight"] and t["zero_lag"])
    node(0.40, 0.40, "fixed-lag component\n(spikelet / other phase?)", t["tight"] and not t["zero_lag"])
    node(0.62, 0.40, "synchrony-like\n(distinct cells?)", not t["tight"] and t["sync"])
    node(0.87, 0.40, "unclear", not t["tight"] and not t["sync"])
    node(0.5, 0.27, f"partial-split check: {EV_COINC[1]} ≤ coinc < 0.9 and (rate corr ≤ −0.3 or amp corr ≤ −{EV_COVAR})"
                    f"  →  {'YES' if t['split'] else 'no'}  (rate {_fmt(r.rate_corr_60s)}, amp {_fmt(r.amp_corr_detr)})",
         t["split"], fs=7)
    sup = r.support if isinstance(r.support, str) and r.support else "–"
    ag = r.against if isinstance(r.against, str) and r.against else "–"
    node(0.5, 0.15, f"votes   +{r.n_support}: {sup}    −{r.n_against}: {ag}", True, fs=7.5)
    node(0.5, 0.04, f"evidence status: {r.evidence_status}   →   MANUAL REVIEW (no automatic merge)", True, fs=8)


def open_detail(cs, ss):
    r = pairs[(pairs.cs == cs) & (pairs.ss == ss)].iloc[0]
    fig = plt.figure(figsize=(15, 8))
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.3, bottom=0.12)
    ax = fig.add_subplot(gs[0, 0])
    shade_windows(ax)
    bar_ccg(ax, CCG_COARSE[(cs, ss)], LAGS_COARSE, COARSE_BIN)
    ax.set_title("Cross-correlogram: SS with respect to CS", fontsize=10)
    label_ccg(ax)
    ax = fig.add_subplot(gs[0, 1])
    bar_ccg(ax, CCG_FINE[(cs, ss)], LAGS_FINE, FINE_BIN)
    ax.set_title(f"fine CCG; zero-peak ratio {r.zero_peak_ratio:.2f}, dip {r.collision_dip}", fontsize=10)
    label_ccg(ax)
    ax = fig.add_subplot(gs[0, 2])
    ra, la, _ = ccg(ST[ss], ST[ss], ACG_FINE_BIN, ACG_FINE_WIN, acg=True)
    bar_ccg(ax, ra, la, ACG_FINE_BIN)
    ax.set_title(f"Autocorrelogram SS {ss}", fontsize=10)
    label_ccg(ax)
    plot_template(fig.add_subplot(gs[1, 0]), ss, color=BAR_COLOR)
    bl, bccgs = block_pause_ratios(win_st(cs, PAIR_WIN[(cs, ss)]), ST_PAD[ss], *PAIR_WIN[(cs, ss)], return_ccgs=True)
    ax = fig.add_subplot(gs[1, 1])
    ax.plot((PAIR_WIN[(cs, ss)][0] + np.arange(len(bl)) * BLOCK_LEN) / 60, bl, "o-", color=CB_COLOR)
    ax.axhline(1, color=(0.6, 0.6, 0.6), lw=1, ls="--")
    ax.axhline(CAND_PAUSE_RATIO_MAX, color=(0.6, 0.6, 0.6), lw=1, ls=":")
    ax.set_xlabel("block start (min)")
    ax.set_ylabel("pause ratio")
    ax.set_title(f"per-block pause ratio ({BLOCK_LEN / 60:.0f} min blocks)", fontsize=10)
    ax = fig.add_subplot(gs[1, 2])
    shade_windows(ax)
    cmap = plt.get_cmap("viridis")
    for k, rc in enumerate(bccgs):
        if rc is not None:
            ax.step(LAGS_COARSE * 1e3, rc, where="post", color=cmap(k / max(len(bccgs) - 1, 1)), lw=1,
                    label=f"{(PAIR_WIN[(cs, ss)][0] + k * BLOCK_LEN) / 60:.0f}′")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.set_title("per-block CCGs", fontsize=10)
    label_ccg(ax)
    checks = "  ".join(f"{h.replace(chr(10), ' ')} {'✓' if f(r[c]) else '✗'}" for h, c, f, _ in GRID_COLS)
    fig.suptitle(f"SS phyID = {ss} (tier {r.ss_tier}, {r.ss_group}, C4 {r.ss_c4}) | CS phyID = {cs} | "
                 f"d = {r.dist_um:.0f} µm | candidate = {r.candidate}\n{checks}", fontsize=10)
    bsave = Button(fig.add_axes([0.88, 0.02, 0.09, 0.05]), "Save PNG")
    bsave.on_clicked(lambda _e: (fig.savefig(FIG_DIR / f"CS{cs}_SS{ss}_pair.png", dpi=120),
                                 print(f"saved {FIG_DIR / f'CS{cs}_SS{ss}_pair.png'}")))
    fig._keep = [bsave]
    if plt.get_backend().lower() != "agg":
        fig.show()
    return fig


def open_dup_detail(a, b):
    """Full evidence figure for CS pair (a, b): every evidence line with the data behind it."""
    rr = [r for p, r in pairs_of_cs(a) if p == b]
    if not rr:
        print(f"no evidence row for pair {a}/{b}")
        return None
    r = rr[0]
    a, b = int(r.cs_a), int(r.cs_b)
    e, aux = dup_evidence(a, b, r.peak_lag_ms, r.coinc_frac_peak)
    ta, tb = aux["ta"], aux["tb"]
    fig = plt.figure(figsize=(17, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.55, wspace=0.3, top=0.9, bottom=0.07)
    # 1. fine CCG
    ax = fig.add_subplot(gs[0, 0])
    rf, lf, _ = ccg(ta, tb, FINE_BIN, FINE_WIN)
    ax.axvspan(r.peak_lag_ms - COINC_DELTA * 1e3, r.peak_lag_ms + COINC_DELTA * 1e3, color=SHADE, zorder=0)
    bar_ccg(ax, rf, lf, FINE_BIN)
    ax.set_title(f"fine CCG {b} rel. {a}; peak {r.peak_lag_ms:+.2f} ms, coinc {r.coinc_frac_peak:.2f}", fontsize=9)
    label_ccg(ax)
    # 2. lag jitter
    ax = fig.add_subplot(gs[0, 1])
    ax.hist(aux["off"] * 1e3, bins=np.linspace(-LAG_SPREAD_WIN * 1e3, LAG_SPREAD_WIN * 1e3, 41), color=BAR_COLOR)
    ax.set_xlabel("offset from peak lag (ms)")
    ax.set_ylabel("count")
    ax.set_title(f"lag jitter: robust SD {_fmt(e['lag_sd_ms'], 3)} ms (vote {e['ev_lag_sd']:+d})", fontsize=9)
    # 3. rates over time
    ax = fig.add_subplot(gs[0, 2])
    ed = aux["edges"]
    co_t = ta[aux["ia"]]
    for y, lab, c in ((aux["ca"], str(a), CB_COLOR), (aux["cb"], str(b), (0.2, 0.5, 0.8)),
                      (np.histogram(co_t, ed)[0], "coincident", "k")):
        ax.step(ed[:-1] / 60, y / DETREND_BIN, where="post", color=c, label=lab)
    ax.set_xlabel("time (min, imec)")
    ax.set_ylabel("rate (sp/s)")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(f"rates, {DETREND_BIN:.0f} s bins: corr {_fmt(e['rate_corr_60s'])}", fontsize=9)
    # 4–5. co-silencing
    plot_cosilencing(fig.add_subplot(gs[1, 0:2]), r, aux)
    # 6. footprints
    ax = fig.add_subplot(gs[1, 2])
    pa, pb = np.ptp(TEMPLATE[a], 0), np.ptp(TEMPLATE[b], 0)
    order = np.argsort(chan_pos[:, 1])
    ax.plot(pa[order] / pa.max(), chan_pos[order, 1], color=CB_COLOR, lw=1, label=str(a))
    ax.plot(pb[order] / pb.max(), chan_pos[order, 1], color=(0.2, 0.5, 0.8), lw=1, label=str(b))
    yb = chan_pos[aux["both"], 1]
    if len(yb):
        ax.axhspan(yb.min(), yb.max(), color=SHADE, zorder=0)
    y0 = (units.at[a, "y"] + units.at[b, "y"]) / 2
    ax.set_ylim(y0 - 250, y0 + 250)
    ax.set_xlabel("template ptp / own max")
    ax.set_ylabel("probe y (µm)")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(f"footprints: overlap {e['footprint_overlap_ch']} ch, cos {_fmt(e['footprint_cos'])}", fontsize=9)
    # 7–8. amplitude / depth co-variation of coincident spikes
    for k, (arr, lab) in enumerate(((amplitudes, "amplitude"), (SPIKE_DEPTH, "depth"))):
        ax = fig.add_subplot(gs[2, k])
        if arr is None or len(aux["ia"]) < 20:
            ax.text(0.5, 0.5, f"{lab}: not available", ha="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue
        xa, xb = _per_spike(arr, a, pair_win(a, b))[aux["ia"]], _per_spike(arr, b, pair_win(a, b))[aux["ib"]]
        ax.scatter(detrend(xa, co_t, pair_win(a, b)[0]), detrend(xb, co_t, pair_win(a, b)[0]), s=3, color=DOT_COLOR, rasterized=True)
        ax.set_xlabel(f"{lab} {a} (detrended)")
        ax.set_ylabel(f"{lab} {b} (detrended)")
        val = e["amp_corr_detr"] if k == 0 else e["depth_corr_detr"]
        ax.set_title(f"{lab} of coincident spikes: r = {_fmt(val)}", fontsize=9)
    # 9. vote card
    ax = fig.add_subplot(gs[2, 2])
    ax.set_axis_off()
    lines = [f"{EV_LABELS[k]:16s} {e[f'ev_{k}']:+d}   {ev_value_text(pd.Series({**r.to_dict(), **e}), k)}"
             for k in EV_LINES]
    ax.text(0, 1, "\n".join(lines + ["", f"union: {e['union_new_isi']} new / {e['union_expected']:.1f} expected",
                                     f"doublets: {_fmt(e['doublet_frac_a'], 3)} / {_fmt(e['doublet_frac_b'], 3)}",
                                     f"ss_cond detail: {e['ss_cond_detail'] or '–'}",
                                     f"SS after coincident: {e['ss_coinc_pr'] or '–'}",
                                     "", f"type: {r.type_hint}", f"status: {r.evidence_status}",
                                     "→ manual review, nothing merged"]),
            va="top", family="monospace", fontsize=8, transform=ax.transAxes)
    fam = r.family if np.isfinite(r.family) else None
    fig.suptitle(f"CS pair {a} / {b}  |  d = {r.dist_um:.0f} µm  |  n = {len(ta)} / {len(tb)} "
                 f"(coincident {e['n_coinc']}, {a}-only {e['n_a_only']}, {b}-only {e['n_b_only']})"
                 f"{f'  |  family {int(fam)}: {FAMILIES.get(int(fam))}' if fam else ''}", fontsize=11)
    bsave = Button(fig.add_axes([0.88, 0.01, 0.09, 0.04]), "Save PNG")
    bsave.on_clicked(lambda _e: (fig.savefig(FIG_DIR / f"CS{a}_CS{b}_evidence.png", dpi=120),
                                 print(f"saved {FIG_DIR / f'CS{a}_CS{b}_evidence.png'}")))
    fig._keep = [bsave]
    if plt.get_backend().lower() != "agg":
        fig.show()
    return fig


def open_browser(start_cs=None):
    B = {"i": CS_IDS.index(start_cs) if start_cs in CS_IDS else 0, "rank": "rank", "view": "pairing",
         "bin": COARSE_BIN, "win": COARSE_WIN, "top": [], "dup_rows": [], "sel": 0, "details": []}
    fig = plt.figure(figsize=(19, 10.5))
    # left: the CS itself
    ax_acg = fig.add_axes([0.04, 0.72, 0.24, 0.19])
    ax_tmp = fig.add_axes([0.04, 0.40, 0.11, 0.26])
    ax_probe = fig.add_axes([0.19, 0.40, 0.09, 0.26])
    ax_info = fig.add_axes([0.03, 0.12, 0.27, 0.24])
    # pairing view
    ax_funnel = fig.add_axes([0.40, 0.63, 0.22, 0.28])
    ax_grid = fig.add_axes([0.67, 0.63, 0.31, 0.28])
    ccg_axes = [fig.add_axes([0.36 + (k % 4) * 0.16, 0.36 if k < 4 else 0.12, 0.13, 0.17]) for k in range(TOP_N)]
    pairing_axes = [ax_funnel, ax_grid] + ccg_axes
    # duplicates view
    ax_votes = fig.add_axes([0.40, 0.63, 0.36, 0.28])
    ax_fam = fig.add_axes([0.79, 0.63, 0.20, 0.28])
    ax_tree = fig.add_axes([0.34, 0.08, 0.38, 0.50])
    ax_dccg = fig.add_axes([0.76, 0.39, 0.22, 0.15])
    ax_cosil = fig.add_axes([0.76, 0.12, 0.22, 0.19])
    dup_axes = [ax_votes, ax_fam, ax_tree, ax_dccg, ax_cosil]

    def draw_left(cs):
        for ax in (ax_acg, ax_tmp, ax_probe, ax_info):
            ax.clear()
        ra, la, _ = ccg(ST[cs], ST[cs], ACG_FINE_BIN, ACG_FINE_WIN, acg=True)
        ax_acg.axvspan(-5, 5, color=SHADE, zorder=0, linewidth=0)
        bar_ccg(ax_acg, ra, la, ACG_FINE_BIN, color=CB_COLOR)
        ax_acg.set_title(f"Autocorrelogram CS {cs} (window)", fontsize=10)
        label_ccg(ax_acg)
        plot_template(ax_tmp, cs)
        plot_probe(ax_probe, cs, cands=B["top"] if B["view"] == "pairing" else [p for p, _ in B["dup_rows"]])
        u = units.loc[cs]
        cand_rows = pairs[(pairs.cs == cs) & pairs.candidate]
        fam = [f for f, m in FAMILIES.items() if cs in m]
        conf_cs = conflicts[conflicts.apply(lambda q: q.unit == cs or cs in q.members, axis=1)] if len(conflicts) else []
        txt = [f"CS {cs}   C4 {u.c4_label} conf {u.c4_conf:.1f} ({u.c4_votes})   manual CS {u.manual_CS}",
               f"window: n = {u.n_spikes}, FR = {u.fr:.2f} sp/s   (full: {u.fr_full:.2f})",
               f"peak ch {int(u.peak_ch)}, y = {u.y:.0f} µm, ISI<1ms {u.frac_isi_lt1ms:.4f}",
               "", f"SS partner candidates ({len(cand_rows)}):"]
        txt += [f"  SS {int(q.ss)} tier {q.ss_tier} {q.ss_group}: pr {q.pause_ratio:.2f}, dur50 {_fmt(q.pause_dur50_ms, 0)} ms,"
                f" d {q.dist_um:.0f} µm" for _, q in cand_rows.iterrows()] or ["  none"]
        txt += ["", f"duplicate family: {FAMILIES[fam[0]] if fam else 'none'}  (not merged)",
                f"conflicts involving this CS: {len(conf_cs)}"]
        ax_info.text(0, 1, "\n".join(txt), va="top", fontsize=8.5, family="monospace", transform=ax_info.transAxes)
        ax_info.set_axis_off()

    def draw_pairing(cs):
        for ax in pairing_axes:
            ax.clear()
        key, asc = RANK_KEYS[B["rank"]]
        sub = pairs[pairs.cs == cs].sort_values(key, ascending=asc)
        B["top"] = sub["ss"].head(TOP_N).tolist()
        # funnel
        m = np.ones(len(sub), bool)
        steps = [(f"all SS candidates", m.copy())]
        for lab, f in CAND_STEPS:
            m &= f(sub).to_numpy()
            steps.append(("+ " + lab, m.copy()))
        claimed = sub["ss"].map(lambda s: len(cs_of.get(s, []))).to_numpy()
        steps.append(("+ one-to-one", m & (claimed <= 1)))
        tiers = sub["ss_tier"].to_numpy()
        for k, (lab, mk) in enumerate(steps):
            left = 0
            for tier in (1, 2):
                n = int((mk & (tiers == tier)).sum())
                ax_funnel.barh(k, n, left=left, color=TIER_COLORS[tier], label=f"tier {tier}" if k == 0 else None)
                left += n
            ax_funnel.text(left + 0.5, k, str(left), va="center", fontsize=8)
        ax_funnel.set_yticks(range(len(steps)))
        ax_funnel.set_yticklabels([s for s, _ in steps], fontsize=8)
        ax_funnel.invert_yaxis()
        ax_funnel.set_xlabel("number of SS units")
        ax_funnel.legend(fontsize=7, frameon=False, loc="lower right")
        ax_funnel.set_title("pairing funnel (cumulative criteria)", fontsize=10)
        # pass/fail grid for the ranked SS
        top = sub.head(TOP_N)
        code = np.array([[1 if f(q[c]) else -1 for _, c, f, _ in GRID_COLS] for _, q in top.iterrows()], float)
        if len(top):
            ax_grid.imshow(code, cmap=VOTE_CMAP, vmin=-1.5, vmax=1.5, aspect="auto")
            for i, (_, q) in enumerate(top.iterrows()):
                for j, (_, c, _, crit) in enumerate(GRID_COLS):
                    nd = 0 if c == "pause_dur50_ms" else (3 if c == "jitter_p" else 2)
                    ax_grid.text(j, i, _fmt(q[c], nd), ha="center", va="center",
                                 fontsize=7, fontweight="bold" if crit else "normal")
            ax_grid.set_yticks(range(len(top)))
            ax_grid.set_yticklabels([f"SS {int(q.ss)} t{q.ss_tier}{' ★' if q.candidate else ''}"
                                     for _, q in top.iterrows()], fontsize=8)
            ax_grid.set_xticks(range(len(GRID_COLS)))
            ax_grid.set_xticklabels([h for h, *_ in GRID_COLS], fontsize=7)
            for j, (*_, crit) in enumerate(GRID_COLS):
                if not crit:
                    ax_grid.get_xticklabels()[j].set_color((0.5, 0.5, 0.5))
            ax_grid.axvline(3.5, color="k", lw=1.5)
        ax_grid.set_title(f"criteria per SS (ranked by {B['rank']}; ★ = candidate; grey = info only; click a row)",
                          fontsize=9)
        # CCGs
        for k, ax in enumerate(ccg_axes):
            if k >= len(B["top"]):
                ax.set_visible(False)
                continue
            ss = B["top"][k]
            q = sub[sub.ss == ss].iloc[0]
            rc, lc, _ = ccg(win_st(cs, PAIR_WIN[(cs, ss)]), ST_PAD[ss], B["bin"], B["win"])
            shade_windows(ax)
            bar_ccg(ax, rc, lc, B["bin"])
            ax.set_title(f"SS {ss} t{q.ss_tier} ({q.ss_c4}) {q.dist_um:.0f} µm\npr {q.pause_ratio:.2f} "
                         f"asym {q.asym:.2f} p {q.jitter_p:.3f}", fontsize=8,
                         color=CB_COLOR if q.candidate else "k", fontweight="bold" if q.candidate else "normal")
            label_ccg(ax, xlabel=k >= 4, ylabel=k % 4 == 0)

    def draw_duplicates(cs):
        for ax in dup_axes:
            ax.clear()
        rows = pairs_of_cs(cs)
        B["dup_rows"] = rows
        if not rows:
            ax_votes.text(0.5, 0.5, "no duplicate-candidate or sharp-peak pair involves this CS",
                          ha="center", va="center", transform=ax_votes.transAxes)
            for ax in dup_axes:
                ax.set_axis_off()
            return
        B["sel"] = min(B["sel"], len(rows) - 1)
        M = np.array([[r[f"ev_{k}"] for _, r in rows] for k in EV_LINES], float)
        ax_votes.set_axis_on()
        ax_votes.imshow(M, cmap=VOTE_CMAP, vmin=-1.5, vmax=1.5, aspect="auto")
        for j, (_, r) in enumerate(rows):
            for i, k in enumerate(EV_LINES):
                ax_votes.text(j, i, ev_value_text(r, k), ha="center", va="center", fontsize=7)
        ax_votes.add_patch(plt.Rectangle((B["sel"] - 0.5, -0.5), 1, len(EV_LINES), fill=False, ec="k", lw=2.5))
        ax_votes.set_yticks(range(len(EV_LINES)))
        ax_votes.set_yticklabels([EV_LABELS[k] for k in EV_LINES], fontsize=8)
        ax_votes.set_xticks(range(len(rows)))
        ax_votes.set_xticklabels([f"{p}{' *' if r.dup_candidate else ''}\n{r.evidence_status.replace('_', ' ')}"
                                  f"\n+{r.n_support} / −{r.n_against}" for p, r in rows], fontsize=7)
        ax_votes.set_title(f"duplicate evidence of CS {cs} vs partner CS (* = dup candidate; green +1 same cell, "
                           f"red −1 distinct, grey untested; click a column)", fontsize=9)
        p, r = rows[B["sel"]]
        # selected-pair summary
        ax_fam.set_axis_off()
        fam = [m for m in FAMILIES.values() if cs in m]
        ax_fam.text(0, 1, "\n".join([
            f"selected: {int(r.cs_a)} / {int(r.cs_b)}",
            f"type: {r.type_hint}", f"status: {r.evidence_status}",
            f"lag {r.peak_lag_ms:+.2f} ms, FWHM {r.peak_fwhm_ms:.2f} ms",
            f"d = {r.dist_um:.0f} µm", f"SS a: {r.ss_a}  SS b: {r.ss_b}",
            f"family: {fam[0] if fam else 'none'}", "", "nothing is merged:", "merging is a manual decision",
            "", "button 'pair detail' → full figure"]),
            va="top", fontsize=8.5, family="monospace", transform=ax_fam.transAxes)
        draw_decision_path(ax_tree, r)
        rf, lf, _ = ccg(*(win_st(int(r[k]), pair_win(int(r.cs_a), int(r.cs_b))) for k in ("cs_a", "cs_b")), FINE_BIN, FINE_WIN)
        ax_dccg.axvspan(r.peak_lag_ms - COINC_DELTA * 1e3, r.peak_lag_ms + COINC_DELTA * 1e3, color=SHADE, zorder=0)
        bar_ccg(ax_dccg, rf, lf, FINE_BIN)
        ax_dccg.set_title(f"fine CCG {int(r.cs_b)} rel. {int(r.cs_a)}", fontsize=9)
        label_ccg(ax_dccg)
        _, aux = dup_evidence(int(r.cs_a), int(r.cs_b), r.peak_lag_ms, r.coinc_frac_peak)
        plot_cosilencing(ax_cosil, r, aux)

    def draw():
        cs = CS_IDS[B["i"]]
        pv = B["view"] == "pairing"
        for ax in pairing_axes:
            ax.set_visible(pv)
        for ax in dup_axes:
            ax.set_visible(not pv)
        if pv:
            draw_pairing(cs)
        else:
            draw_duplicates(cs)
        draw_left(cs)
        u = units.loc[cs]
        fig.suptitle(f"CS phyID = {cs}  ({B['i'] + 1}/{len(CS_IDS)})  |  view: {B['view']}  |  "
                     f"{MOUSE} {SESSION_R} {IMEC} ({SORTER_TYPE}), window {ANALYSIS_WINDOW} "
                     f"[{T0:.0f}, {T1:.0f}] s", fontsize=11)
        fig.canvas.draw_idle()

    def go(i):
        B["i"] = i % len(CS_IDS)
        B["sel"] = 0
        draw()

    def on_click(ev):
        if ev.button != 1 or ev.inaxes is None or not ev.inaxes.get_visible():
            return
        cs = CS_IDS[B["i"]]
        if ev.inaxes in ccg_axes:
            k = ccg_axes.index(ev.inaxes)
            if k < len(B["top"]):
                B["details"].append(open_detail(cs, B["top"][k]))
        elif ev.inaxes is ax_grid and ev.ydata is not None:
            k = int(round(ev.ydata))
            if 0 <= k < len(B["top"]):
                B["details"].append(open_detail(cs, B["top"][k]))
        elif ev.inaxes is ax_votes and ev.xdata is not None and B["dup_rows"]:
            k = int(round(ev.xdata))
            if 0 <= k < len(B["dup_rows"]):
                B["sel"] = k
                draw()

    def jump(text):
        try:
            c = int(text)
        except ValueError:
            return
        if c in CS_IDS:
            go(CS_IDS.index(c))
        else:
            print(f"{c} is not a CS in this session; CS phyIDs: {CS_IDS}")

    def pair_detail(_e):
        if B["view"] == "duplicates" and B["dup_rows"]:
            p, r = B["dup_rows"][B["sel"]]
            B["details"].append(open_dup_detail(int(r.cs_a), int(r.cs_b)))

    def save_view(_e=None):
        out = safe_out_dir(FIG_DIR / "browser")
        p = out / f"CS{CS_IDS[B['i']]}_{B['view']}.png"
        fig.savefig(p, dpi=100)
        print(f"saved {p}")

    widgets = []
    b_prev = Button(fig.add_axes([0.03, 0.03, 0.05, 0.045]), "< Prev")
    b_next = Button(fig.add_axes([0.085, 0.03, 0.05, 0.045]), "Next >")
    b_prev.on_clicked(lambda _e: go(B["i"] - 1))
    b_next.on_clicked(lambda _e: go(B["i"] + 1))
    tb = TextBox(fig.add_axes([0.18, 0.03, 0.05, 0.045]), "phyID ")
    tb.on_submit(jump)
    r_view = RadioButtons(fig.add_axes([0.25, 0.005, 0.07, 0.09]), ["pairing", "duplicates"], active=0)
    r_view.on_clicked(lambda lab: (B.update(view=lab, sel=0), draw()))
    r_rank = RadioButtons(fig.add_axes([0.33, 0.005, 0.07, 0.09]), list(RANK_KEYS), active=0)
    r_rank.on_clicked(lambda lab: (B.update(rank=lab), draw()))
    s_bin = Slider(fig.add_axes([0.47, 0.055, 0.18, 0.02]), "bin (ms)", 0.5, 2.0, valinit=COARSE_BIN * 1e3,
                   valstep=0.25)
    s_win = Slider(fig.add_axes([0.47, 0.02, 0.18, 0.02]), "window ± (ms)", 25, 150, valinit=COARSE_WIN * 1e3,
                   valstep=5)
    s_bin.on_changed(lambda v: (B.update(bin=v * 1e-3), draw()))
    s_win.on_changed(lambda v: (B.update(win=v * 1e-3), draw()))
    b_pair = Button(fig.add_axes([0.80, 0.03, 0.08, 0.045]), "pair detail")
    b_pair.on_clicked(pair_detail)
    b_save = Button(fig.add_axes([0.90, 0.03, 0.08, 0.045]), "save view")
    b_save.on_clicked(save_view)
    widgets += [b_prev, b_next, tb, r_view, r_rank, s_bin, s_win, b_pair, b_save]
    fig.canvas.mpl_connect("button_press_event", on_click)
    B.update(fig=fig, widgets=widgets, draw=draw, go=go, save_view=save_view, on_click=on_click)
    draw()
    if plt.get_backend().lower() != "agg":
        fig.show()
    return B


if OPEN_QT_BROWSER:
    BROWSER = open_browser()

if OPEN_QT_BROWSER and plt.get_backend().lower() == "agg":
    # headless run: snapshots of both views for CS with a partner or a duplicate pair, and every dup-candidate pair
    _snap = sorted(set(pairs[pairs.candidate].cs) | set(dup_ev.cs_a) | set(dup_ev.cs_b)) if len(dup_ev) else []
    for _c in _snap:
        BROWSER["go"](CS_IDS.index(_c))
        for _v in ("pairing", "duplicates"):
            BROWSER.update(view=_v, sel=0)
            BROWSER["draw"]()
            BROWSER["save_view"]()
    for _, _r in (dup_ev[dup_ev.dup_candidate].iterrows() if len(dup_ev) else []):
        _f = open_dup_detail(int(_r.cs_a), int(_r.cs_b))
        if _f is not None:
            _f.savefig(safe_out_dir(FIG_DIR / "browser") / f"CS{int(_r.cs_a)}_CS{int(_r.cs_b)}_evidence.png", dpi=100)
            plt.close(_f)
    print(f"headless: saved browser snapshots to {FIG_DIR / 'browser'}")


# %% 8b. HTML BROWSER EXPORT
# Writes OUT_DIR/browser_data.js (every number the page shows, precomputed here) and OUT_DIR/browser.html (a copy
# of browser_template.html next to this script). Open browser.html in any web browser; it works offline. To change
# the page, edit browser_template.html and copy it over browser.html again: the data does not need regenerating.
HTML_CCG_WIN = 100e-3                    # exported coarse CCGs: 1 ms bins over ± this (page crops / rebins)
HTML_MAX_POINTS = 2000                   # amplitude-over-time points per CS; scatter points per evidence panel
TMPL_OVERLAY_NCH = 4                     # template overlay of a CS pair: each unit's top-ptp channels (union)


def _r(x, nd=3):
    """JSON-safe: numpy → python, arrays → lists, NaN/inf → None, floats rounded to nd decimals."""
    if isinstance(x, dict):
        return {str(k): _r(v, nd) for k, v in x.items()}
    if isinstance(x, (list, tuple, set)):
        return [_r(v, nd) for v in (sorted(x) if isinstance(x, set) else x)]
    if isinstance(x, np.ndarray):
        return _r(x.tolist(), nd)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return round(float(x), nd) if np.isfinite(x) else None
    if x is None or isinstance(x, str):
        return x
    try:
        if pd.isna(x):
            return None
    except (TypeError, ValueError):
        pass
    return str(x)


def template_export(c):
    t = TEMPLATE[c]
    pk = int(units.at[c, "peak_ch"])
    chs = np.argsort(np.hypot(*(chan_pos - chan_pos[pk]).T))[:2 * N_NEIGH_CH + 1]
    return dict(peak=pk, ch=chs, dx=chan_pos[chs, 0] - chan_pos[pk, 0], dy=chan_pos[chs, 1] - chan_pos[pk, 1],
                w=np.round(t[:, chs].T, 3))


# sample of the spike time inside a template (sorters align templates on a fixed sample): the most common
# position of the largest deflection on the peak channel, over the CS and SS units
T_ZERO_SAMPLE = int(pd.Series([int(np.argmax(np.abs(TEMPLATE[c][:, int(units.at[c, "peak_ch"])])))
                               for c in CS_IDS + SS_IDS]).mode()[0])
print(f"template time zero: sample {T_ZERO_SAMPLE} of {TEMPLATE[CS_IDS[0]].shape[0]}")


def template_overlay(a, b, lag_s, offsets_s):
    """Templates of a (at its spike time) and b on the union of both units' top channels, plus the lag of b
    relative to a (fine-CCG peak) and its spread over coincident events, for the overlay panel."""
    pa, pb = np.ptp(TEMPLATE[a], 0), np.ptp(TEMPLATE[b], 0)
    chs = np.unique(np.r_[np.argsort(pa)[::-1][:TMPL_OVERLAY_NCH], np.argsort(pb)[::-1][:TMPL_OVERLAY_NCH]])
    chs = chs[np.argsort(-chan_pos[chs, 1])]                       # top of the probe first
    lags = lag_s + np.asarray(offsets_s)
    q = np.percentile(lags, [10, 25, 50, 75, 90]) * 1e3 if len(lags) >= 5 else None
    return dict(ch=chs, x=chan_pos[chs, 0], y=chan_pos[chs, 1], a=np.round(TEMPLATE[a][:, chs].T, 2),
                b=np.round(TEMPLATE[b][:, chs].T, 2), peak_a=int(np.argmax(pa)), peak_b=int(np.argmax(pb)),
                t_zero=T_ZERO_SAMPLE, sr=SR, lag_ms=lag_s * 1e3, lag_q_ms=q, n_lags=len(lags))


def _subsample(n, k, seed=0):
    return np.arange(n) if n <= k else np.sort(np.random.default_rng(seed).choice(n, k, replace=False))


_t = time.time()
H = dict(meta=dict(mouse=MOUSE, session_r=SESSION_R, session_s=SESSION_S, imec=IMEC, sorter_type=SORTER_TYPE,
                   sorter_dir=str(SORTER_DIR), window=ANALYSIS_WINDOW, t0=T0, t1=T1, t_rec=T_REC, sample_rate=SR,
                   b0=B0, b1=B1, stable_span=STABLE_SPAN,
                   run_tag=RUN_TAG, created=time.strftime("%Y-%m-%d %H:%M")))
H["params"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in globals().items()
               if k.isupper() and isinstance(v, (int, float, str, tuple)) and not k.startswith("_")}
H["chan_pos"] = chan_pos
ucols = ["role", "group", "extragood", "c4_label", "c4_conf", "c4_votes", "region", "manual_CS", "peak_ch", "x", "y",
         "n_spikes", "fr", "fr_full", "frac_isi_lt1ms", "span_start", "span_end", "n_spikes_span",
         # v8
         "cs_source", "cs_c4", "cs_manual", "cs_screen", "scr_fr", "cs_fr_ok", "isi_5_35_frac", "isi_5_35_ok",
         "cs_isi_10_35_frac", "isi_lt50_frac", "cv", "cv2", "presence_bins_frac", "partial_presence", "poor_isolation",
         "tmpl_dur_ms", "tmpl_n_fp_ch", "tmpl_n_late_peaks", "manual_kept_as_ss"]
H["units"] = {int(c): dict(units.loc[c, ucols].to_dict(), template=template_export(c)) for c in CS_IDS + SS_IDS}
H["cs_ids"], H["ss_ids"] = CS_IDS, SS_IDS
H["stab"] = {int(c): _r(d[["t0", "t1", "fr_ratio", "amp_ratio", "in_bhv", "ok"]].to_dict("list"), 3)
             for c, d in STAB.items()}

# CS overview data
H["cs"] = {}
_isi_edges = np.logspace(-3, np.log10(5), 61)
_e30 = np.arange(0, T_REC + 30, 30)
for c in CS_IDS:
    st = ST[c]
    acg_f = ccg(st, st, ACG_FINE_BIN, ACG_FINE_WIN, acg=True)[0]
    acg_c = ccg(st, st, ACG_COARSE_BIN, ACG_COARSE_WIN, acg=True)[0]
    k = _subsample(len(ST_FULL[c]), HTML_MAX_POINTS, seed=c)
    d = dict(acg_fine=acg_f, acg_coarse=acg_c, isi_edges=_isi_edges, isi=np.histogram(np.diff(st), _isi_edges)[0],
             isi_short=float(np.mean(np.diff(st) < CS_MIN_ISI)) if len(st) > 1 else None,
             fr_t=_e30[:-1], fr=np.histogram(ST_FULL[c], _e30)[0] / 30,
             amp_t=ST_FULL[c][k], amp=np.asarray(amplitudes[SPK_IDX_FULL[c][k]], float).ravel())
    if drift_um is not None:
        kk = int(np.argmin(np.abs(drift_depths - units.at[c, "y"])))
        step = max(1, len(drift_t) // 1500)
        d.update(drift_t=drift_t[::step], drift=drift_um[::step, kk], drift_depth=drift_depths[kk])
    H["cs"][c] = _r(d, 3)

# SS fine ACGs (detail panel)
H["ss_acg"] = {s: _r(ccg(ST[s], ST[s], ACG_FINE_BIN, ACG_FINE_WIN, acg=True)[0], 2) for s in SS_IDS}

# CS×SS pairs: every pair gets a ±HTML_CCG_WIN coarse CCG; fine CCG + per-block CCGs for candidates and top ranks
pcols = ["cs", "ss", "rank", "candidate", "ss_tier", "ss_group", "ss_c4", "pause_ratio", "mirror_ratio", "asym",
         "pause_min_lag_ms", "pause_dur50_ms", "pause_dur80_ms", "jitter_p", "jitter_z", "jitter_deficit",
         "lussac_pass", "collision_dip", "zero_peak_ratio", "cs_share_frac", "cs_share_null", "cs_share_excess",
         "ss_from_cs_frac", "cs_share_flag", "cs_event_excess", "ss_cs_locked_frac", "cs_event_flag",
         "pause_ratio_0_8", "mirror_ratio_0_8", "asym_0_8", "block_pr_min", "block_pr_median", "block_pr_max",
         "block_pause_ratios", "dist_um", "baseline_rate", "ss_fr", "n_cs", "n_ss", "win_start", "win_end",
         "candidate_ref", "pause_ratio_ref", "asym_ref", "jitter_p_ref", "n_cs_ref",
         # v8
         "candidate_lit", "lat_n", "lat_post_med_ms", "lat_pre_med_ms", "lat_ratio", "lat_p", "lat_cliff", "refr_ratio",
         "below_refr_ms", "ss_isi_med_ms", "pause_isi_units", "lpost_isi_units", "dur_rate_ok", "lussac_paper_v1_pass",
         "lussac_code_main_pass", "lussac_code_main_fr_ok", "dy_um", "dx_um", "far_flag", "cs_source"]
_detail = set(map(tuple, pairs.loc[pairs.candidate | (pairs["rank"] <= 2 * TOP_N), ["cs", "ss"]].to_numpy()))
H["pairs"] = {}
for _, q in pairs.iterrows():
    cs, ss = int(q.cs), int(q.ss)
    d = {k: q[k] for k in pcols}
    d["ccg"] = np.round(ccg(win_st(cs, PAIR_WIN[(cs, ss)]), ST_PAD[ss], COARSE_BIN, HTML_CCG_WIN)[0], 1)
    if (cs, ss) in _detail:
        d["fine"] = np.round(CCG_FINE[(cs, ss)], 1)
        bl, bc = block_pause_ratios(win_st(cs, PAIR_WIN[(cs, ss)]), ST_PAD[ss], *PAIR_WIN[(cs, ss)], return_ccgs=True)
        d["blocks"] = [None if b is None else np.round(b, 1) for b in bc]
        d["block_t0"] = PAIR_WIN[(cs, ss)][0] + np.arange(len(bc)) * BLOCK_LEN
        # v8: latency histograms (2 ms bins, 0–100 ms) and SS|SS (SS autocorrelogram, 1 ms bins, ±REFR_WIN)
        if (cs, ss) in LAT_HIST:
            d["lat_post_hist"], d["lat_pre_hist"] = LAT_HIST[(cs, ss)]
        _w = PAIR_WIN[(cs, ss)]
        _k = (ss, round(_w[0], 6), round(_w[1], 6))
        if _k in ACG_REFR:
            d["acg_refr"] = np.round(ACG_REFR[_k][0], 1)
    H["pairs"][f"{cs}_{ss}"] = _r(d, 3)
H["cs_of"] = {int(s): [int(c) for c in v] for s, v in cs_of.items()}

# CS–CS duplicate evidence, with the data behind every line
H["dups"] = []
for _, r in (dup_ev.iterrows() if len(dup_ev) else []):
    a, b = int(r.cs_a), int(r.cs_b)
    e, aux = dup_evidence(a, b, r.peak_lag_ms, r.coinc_frac_peak)
    _, tests = classify_pair(r, r.peak_lag_ms, r.dup_candidate)
    ta = aux["ta"]
    co_t = ta[aux["ia"]]
    d = dict(r.to_dict(), tests=tests, n_coinc=e["n_coinc"], n_a_only=e["n_a_only"], n_b_only=e["n_b_only"],
             union_expected=e["union_expected"])
    d["fine"] = np.round(ccg(ta, aux["tb"], FINE_BIN, FINE_WIN)[0], 1)
    d["lag_hist"] = np.histogram(aux["off"] * 1e3, np.linspace(-LAG_SPREAD_WIN * 1e3, LAG_SPREAD_WIN * 1e3, 41))[0]
    d["rates"] = dict(t=aux["edges"][:-1], a=aux["ca"] / DETREND_BIN, b=aux["cb"] / DETREND_BIN,
                      coinc=np.histogram(co_t, aux["edges"])[0] / DETREND_BIN)
    subsets = {"coincident": co_t, f"{a}-only": aux["a_only"], f"{b}-only": aux["b_only"]}
    d["cosil"] = []
    for s in sorted(set(r.ss_a) | set(r.ss_b))[:2]:
        for lab, ref in subsets.items():
            y = norm_ccg(span_cut(ref, s), ST_PAD[s])
            d["cosil"].append(dict(ss=s, subset=lab, n=len(ref), pr=pause_after(span_cut(ref, s), ST_PAD[s]),
                                   y=None if y is None else np.round(y, 3)))
    pa, pb = np.ptp(TEMPLATE[a], 0), np.ptp(TEMPLATE[b], 0)
    o = np.argsort(chan_pos[:, 1])
    yb = chan_pos[aux["both"], 1]
    d["footprint"] = dict(y=chan_pos[o, 1], x=chan_pos[o, 0], a=pa[o] / pa.max(), b=pb[o] / pb.max(),
                          overlap=[float(yb.min()), float(yb.max())] if len(yb) else None)
    for key, arr in (("amp", amplitudes), ("depth", SPIKE_DEPTH)):
        if arr is None or len(aux["ia"]) < 20:
            d[key] = None
            continue
        xa, xb = _per_spike(arr, a, pair_win(a, b))[aux["ia"]], _per_spike(arr, b, pair_win(a, b))[aux["ib"]]
        k = _subsample(len(xa), HTML_MAX_POINTS)
        d[key] = dict(x=detrend(xa, co_t, pair_win(a, b)[0])[k], y=detrend(xb, co_t, pair_win(a, b)[0])[k])
    d["ss_cond_detail"], d["ss_coinc_pr"] = e["ss_cond_detail"], e["ss_coinc_pr"]
    d["overlay"] = template_overlay(a, b, r.peak_lag_ms * 1e-3, aux["off"])
    H["dups"].append(_r(d, 3))
H["families"] = {int(k): v for k, v in FAMILIES.items()}
H["decomp"] = _r(decomp.drop(columns=["win_start", "win_end"]).to_dict("records"))
H["lowfr"] = _r(lowfr_pairs[[c for c in ["cs", "ss", "ss_fr", "pause_ratio", "pause_ratio_0_8", "jitter_p", "asym",
                                         "dist_um", "cs_event_excess", "would_be_candidate"]
                             if c in lowfr_pairs.columns]].to_dict("records")) if len(lowfr_pairs) else []
H["conflicts"] = _r(conflicts.to_dict("records"))
# v8 tables (CS units by source, screen + calibration, kept-as-SS, spatial check, agreement)
H["v8"] = dict(
    lat_edges_ms=LAT_EDGES * 1e3, refr_win_ms=REFR_WIN * 1e3, manual_file=MANUAL_FILE,
    include_manual=INCLUDE_MANUAL_CS, include_screen=INCLUDE_SCREEN_CS, screen_non_eg=SCREEN_NON_EG,
    cs=cs_tab.reset_index().rename(columns={"index": "phyID", "phyID": "phyID"}).to_dict("records"),
    screen=screen_tab.reset_index().rename(columns={"index": "phyID"}).to_dict("records"),
    calib=calib.reset_index().to_dict("records") if len(calib) else [],
    agree=({f"manual={bool(i)}|screen={bool(j)}": int(agree_tab.loc[i, j]) for i in agree_tab.index for j in agree_tab.columns}
           if len(agree_tab) else {}),
    kept=kept_tab.to_dict("records"), dist=dist_tab.to_dict("records"),
    families_c4={int(k): v for k, v in FAMILIES_C4.items()})
H["v8"] = _r(H["v8"])

_HR = _r({k: v for k, v in H.items() if k != "params"})
_HR["params"] = _r(H["params"], 12)          # full precision: bins and windows are fractions of a ms
_js = "window.CS_DATA = " + json.dumps(_HR, separators=(",", ":"), allow_nan=False) + ";\n"
(OUT_DIR / "browser_data.js").write_text(_js, encoding="utf-8")
_tpl = CODE_DIR / "CS_curation" / "browser_template.html"
if _tpl.exists():
    (OUT_DIR / "browser.html").write_text(_tpl.read_text(encoding="utf-8"), encoding="utf-8")
    print(f"HTML browser: {OUT_DIR / 'browser.html'}  (data {len(_js) / 1e6:.1f} MB, {time.time() - _t:.0f} s)")
else:
    print(f"browser_data.js written ({len(_js) / 1e6:.1f} MB); template missing: {_tpl}")


# %% 9. SAVE + SUMMARY
pairs.to_csv(TAB_DIR / "cs_ss_pairs.csv", index=False)
lowfr_pairs.to_csv(TAB_DIR / "cs_ss_pairs_lowfr_excluded.csv", index=False)
decomp.to_csv(TAB_DIR / "cs_event_decomposition.csv", index=False)
dups.to_csv(TAB_DIR / "cs_cs_pairs.csv", index=False)
dup_ev.to_csv(TAB_DIR / "cs_cs_evidence.csv", index=False)
(TAB_DIR / "cs_dup_families.json").write_text(json.dumps({str(k): v for k, v in FAMILIES_C4.items()}, indent=2))  # C4 only (= v7)
(TAB_DIR / "cs_dup_families_all.json").write_text(json.dumps({str(k): v for k, v in FAMILIES.items()}, indent=2))  # v8: all CS
cs_tab.to_csv(TAB_DIR / "cs_units_v8.csv")                                     # v8 reports
screen_tab.to_csv(TAB_DIR / "cs_screen.csv")
kept_tab.to_csv(TAB_DIR / "manual_kept_as_ss.csv", index=False)
dist_tab.to_csv(TAB_DIR / "distance_summary.csv", index=False)
units.to_csv(TAB_DIR / "unit_table.csv")                                       # with the v8 columns
conflicts.assign(members=conflicts["members"].astype(str)).to_csv(TAB_DIR / "conflicts.csv", index=False)
try:
    pairs.to_parquet(TAB_DIR / "cs_ss_pairs.parquet", index=False)
    dups.to_parquet(TAB_DIR / "cs_cs_pairs.parquet", index=False)
except Exception as e:
    print(f"parquet skipped: {e}")
params_used = {k: (list(v) if isinstance(v, tuple) else v) for k, v in globals().items()
               if k.isupper() and isinstance(v, (int, float, str, tuple)) and not k.startswith("_")}
params_used.update(sorter_dir=str(SORTER_DIR), sorter_type=SORTER_TYPE, sample_rate=SR, t_rec_s=T_REC,
                   window_imec_s=[T0, T1], window_len_s=T_WIN, behaviour_imec_s=[B0, B1], sync=SYNC)
(OUT_DIR / "params_used.json").write_text(json.dumps(params_used, indent=2, default=str))

_ctype = conflicts["type"].value_counts().to_dict() if len(conflicts) else {}
print("=" * 70)
print(f"SUMMARY  {MOUSE} {SESSION_R}/{SESSION_S} {IMEC} ({SORTER_TYPE})  |  sorter: {SORTER_DIR.name}")
print(f"  reference window ({ANALYSIS_WINDOW}): imec [{T0:.1f}, {T1:.1f}] s = {T_WIN / 60:.1f} min | behaviour "
      f"[{B0:.0f}, {B1:.0f}] s | stable spans {'ON' if STABLE_SPAN else 'off'}: CS span median "
      f"{units.loc[CS_IDS, 'span_s'].median() / 60:.1f} min")
print(f"  CS×SS candidates: {int(pairs['candidate'].sum())} (reference window: {int(pairs['candidate_ref'].sum())}, "
      f"changed {int((pairs['candidate'] != pairs['candidate_ref']).sum())})")
print(f"  population (extragood & {REGION}): {int(pop.sum())} | CS: {len(CS_IDS)} | SS candidates: {len(SS_IDS)}")
print(f"  SS candidates by tier: 1 (extragood) {int((units['role'] == 'SS_cand').sum())} | "
      f"2 ({SS2_GROUP_TXT}, >= {SS2_FR_MIN} Hz) {int((units['role'] == 'SS_cand2').sum())}")
print(f"  CS with >=1 permissive SS candidate: {int((best['n_candidates'] > 0).sum())} / {len(CS_IDS)} "
      f"(tier-1 partner: {pairs[pairs.candidate & (pairs.ss_tier == 1)].cs.nunique()}, "
      f"tier-2 only: {len(set(pairs[pairs.candidate].cs) - set(pairs[pairs.candidate & (pairs.ss_tier == 1)].cs))})")
print(f"  CS with >=1 Lussac-pass SS:          {int((best['n_lussac'] > 0).sum())} / {len(CS_IDS)}")
print(f"  SS units with >=1 candidate CS:      {len(cs_of)} / {len(SS_IDS)}")
print(f"  SS claimed by >1 CS: {_ctype.get('SS_claimed_by_multiple_CS', 0)} | "
      f"CS with >1 SS: {_ctype.get('CS_with_multiple_SS', 0)}")
print(f"  CS–CS duplicate candidates: {int(dups['dup_candidate'].sum())} "
      + str({k.split(':')[1]: v for k, v in _ctype.items() if k.startswith('CS_duplicate')})
      + f" | sharp-peak-only (not flagged): {int(dups['sharp_peak_only'].sum())}")
print(f"  duplicate families (not merged): {FAMILIES if FAMILIES else 'none'}")
_c4p = pairs[pairs.cs_c4]
print(f"  v8 | C4 CS only (= v7 population): CS with >=1 candidate {_c4p[_c4p.candidate].cs.nunique()} / "
      f"{int(units['cs_c4'].sum())}, candidate_lit {_c4p[_c4p.candidate_lit].cs.nunique()}; pairs candidate "
      f"{int(_c4p.candidate.sum())}, candidate_lit {int(_c4p.candidate_lit.sum())}, both {int((_c4p.candidate & _c4p.candidate_lit).sum())}")
print(f"  v8 | all CS {len(CS_IDS)} " + str(units.loc[CS_IDS, 'cs_source'].value_counts().to_dict())
      + f": paired (candidate | candidate_lit) {int(cs_tab.paired.sum())}; recovered CS {cs_tab.index[cs_tab.recovered_cs].tolist()}")
print(f"  v8 | manual file {MANUAL_FILE}; kept as SS (option A) {KEPT_AS_SS}; screen pool {len(SCREEN_POOL)}, "
      f"screen CS {int(cs_tab.cs_screen.sum())}; left the low-FR audit {_left_lowfr}")
print(f"  v8 | false-positive rate > {FP_UM:.0f} µm: candidate {_far.candidate.mean():.3f}, candidate_lit "
      f"{_far.candidate_lit.mean():.3f}, Lussac paper v1 {_far.lussac_paper_v1_pass.mean():.3f}, code main "
      f"{_far.lussac_code_main_pass.mean():.3f} (n {len(_far)})" if len(_far) else "  v8 | no pairs > 500 µm")
print(f"  v8 | duplicate families C4 only (cs_dup_families.json): {FAMILIES_C4 if FAMILIES_C4 else 'none'}")
print(f"  outputs: {OUT_DIR}")
print("=" * 70)
