# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 03 -- MAIN ANALYSIS: point-process GLM with smooth causal kernels (gctools/kernels.py).
#   y = RAW spike counts, 2.5 ms bins; x = RAW kinematic position KIN_VAR (default med-lat; 2.5 ms = camera frames).
#   kin -> neuron (x->y), Poisson:  log lam = b0 + own_y(H) [+ event] [+ nuisance] + cross(x)
#       own history H chosen PER NEURON by trial-wise CV among H_LIST (never < 150 ms), or
#       fixed at 150 ms (HIST_MODES); built from causal bases (2 single lags + 6 log bumps).
#   neuron -> kin (y->x), linear:   x = c + own_x(H_kin) [+ event] [+ nuisance] + cross(y)
#       own history H_kin chosen per session by trial-wise CV among H_KIN_LIST, fixed (150 ms)
#       or the longest (KIN_HIST_MODES). With KIN_SOURCE = "long" (REACH_LONG, scripts/00b) the
#       kinematic history before -375 ms is real data, so a long H_kin costs no rows and the
#       fitted rows start at -375 ms (with "short": -225 ms, and H_kin > 150 cuts rows).
#   cross(.) = the other signal's past through N_CROSS_BASES causal raised-cosine bumps.
#   Significance: trial shuffle within spout side (+ adaptive extension), BH-FDR per region
#   x direction pooled over sessions. chi2 / F p secondary.
#   Effect size (primary): EXCESS over the shuffle = observed statistic minus the mean of the
#       neuron's own shuffle null, in bits/spike (kin->neuron) and % partial R2 (neuron->kin):
#       the TRIAL-SPECIFIC part, the same thing significance tests (~0 without an effect, also
#       under the shared-drive null). Secondary: CROSS-VALIDATED gain of full over reduced
#       (held-out) = total predictive gain, INCLUDING the shared reach shape of the spout side.
#   Lags: kernel beta(tau), peak and centroid, trial-bootstrap CIs (significant neurons).
#   Adequacy: flags neurons whose reduced model still depends on OLDER own history.
# Wording: "predictive precedence", never "drives"/"causes".
#
# ── HOW TO RUN ───────────────────────────────────────────────────────────────
#   simulations first (DATA_SOURCE = "sim"), then the real data:
#   DATA_SOURCE = "real"   SESSIONS = "all"   N_SHUFFLE = 1000          <<<
#   cached per session: an interrupted run resumes; CACHE_MODE = "load" re-draws figures
#   without refitting.
import os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

_here = globals().get("__vsc_ipynb_file__") or globals().get("__file__")
_starts = ([Path(_here).resolve().parent] if _here else []) + [Path.cwd()]
_cands = [q for s in _starts for p in [s, *s.parents]
          for q in (p, p / "04_population_linear_dependence_py",
                    p / "ephys" / "04_population_linear_dependence_py")]
ROOT = next((p for p in _cands if (p / "gctools").is_dir()), None)
if ROOT is None:
    raise RuntimeError(f"gctools/ not found from {_starts}")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
os.environ["PYTHONPATH"] = os.pathsep.join([str(ROOT), os.environ.get("PYTHONPATH", "")])
import importlib
from gctools import (config, io, preprocess as pp, models, nulls, stats, sim, cache, plotting,
                     runs, kernels)
for _m in (config, io, pp, models, nulls, stats, sim, cache, plotting, runs, kernels):
    importlib.reload(_m)

# ── data ─────────────────────────────────────────────────────────────────────
DATA_SOURCE  = "real"               # "sim" | "real"                        <<< real run: "real"
SESSIONS     = "all"             # real only: "first" | "all" | list of (mouse, sess)  <<< "all"
MIN_STRATUM  = 20                  # skip a session if any spout side has fewer trials
SIM_SCENARIO = "y_to_x" # "null_independent" | "x_to_y" | "y_to_x" | "null_shared_drive"
                                   # (rate follows the spout-side MEAN reach: task-locked, no
                                   # trial-specific coupling -> must NOT be detected)
SIM_PARAMS   = dict(rate_hz="match",   # target rates = the REAL neurons' rates
                    n_neurons=20,      # subset per seed (None = all) -- 03 is slower than 01/02
                    n_drivers=5,       # y_to_x: leaders per simulated session
                    coupling=0.3,      # x_to_y: log-rate gain; y_to_x: kin SD units per leader
                    lag_ms=50.0)       # mode of the true coupling kernel
SIM_SEEDS    = tuple(range(3))    # each seed = one simulated session (5 x 20 neurons)
SESSION      = config.FIRST_PASS[0]   # template session for simulations
REACH_CAT    = config.REACH_CAT
MIN_RATE_HZ  = 1.0                 # per-NEURON inclusion
KIN_VAR      = "y"                 # kinematic variable: "y" (med-lat) | "x" | "z" (position, mm).
                                   # Not "y" -> run name gets _kin<var>, cache keyed by it.
KIN_SOURCE   = "long"              # "long" = REACH_LONG (real history before -375 ms; run 00b
                                   # first) | "short" = the 201-frame reach_px (as 01/02)

# ── model ────────────────────────────────────────────────────────────────────
BIN_MS        = 2.5
HIST_MODES    = ("cv",)            # NEURON own history: "cv" = chosen per neuron among H_LIST by
                                   # CV; "fixed" = 150 ms; "max" = max(H_LIST). Several, e.g.
                                   # ("fixed", "cv", "max"), run side by side and flag
                                   # HISTORY-SENSITIVE neurons (significance differs).
H_LIST        = (150, 200, 300)    # neuron own-history candidates (ms); never shorter than the
                                   # cross range (01: too short -> fake y->x)
KIN_HIST_MODES = ("max",)          # KINEMATICS own history: "max" (max(H_KIN_LIST); default --
                                   # sims: CV always picked the longest; 900 ms clears the
                                   # adequacy flag) | "cv" (per session, held-out SSE) | "fixed" (H_KIN)
H_KIN_LIST    = (450, 600, 900)    # kinematics own-history candidates (ms); > 150 needs
                                   # KIN_SOURCE = "long"
H_KIN         = 150                # the "fixed" kinematics own history (ms)
FIT_FROM_MS   = -225               # fitted rows start here (ms re max reach); None = as early as every
                                   # term exists (-375 with KIN_SOURCE "long", -225 with "short").
                                   # -225 = the rows of the short-kinematics runs (diagnostic)
SCALE_BASES   = True               # own-history log bumps grow with H (~2 per doubling:
                                   # 150 -> 6, 300 -> 8, 450 -> 9); False = 6 for every H
N_CROSS_BASES = 6                  # raised-cosine bumps over CROSS_LAG_RANGE_MS
BASIS_SPACING = "log"              # "log" | "uniform"
EVENT_KERNEL  = False              # True: 4 causal bumps active only AFTER max reach (t = 0)
TASK_MEAN_NUISANCE = False         # True: + leave-one-out same-spout mean of the target at the
                                   # same time (a time-varying intercept; never the trial's own)
RIDGE_ALPHA   = 0.1                # Poisson, on standardised coefficients (see 02)
N_FOLDS       = 5                  # trial-wise CV folds, stratified by spout side
ADEQUACY_MS   = 150                # adequacy test: own history over (H, H + ADEQUACY_MS]
ADEQ_P        = 0.01               # flag a neuron if the adequacy p < ADEQ_P

# ── inference ────────────────────────────────────────────────────────────────
N_SHUFFLE    = 1000                 #                                      <<< real run: 1000
SHUFFLE_SEED = 0
FDR_Q        = 0.05
ADAPTIVE_SHUFFLE = True            # as in 01/02 (Besag & Clifford; same shuffle sequence)
ADAPTIVE_CAP = "auto"
EXTRA_BATCH  = 1000
SEQ_H        = 10
N_BOOT       = 200                  # trial bootstrap of the lag kernel (FDR-significant only)  <<< real run: 200

N_JOBS       = 16
CACHE_MODE   = "load"              # "load" | "rerun"
SAVE_CACHE   = True
SHOW_FIGS    = True
SAVE_FIGS    = True
SAVE_VARS    = True
OUT_SUBDIR   = "03_KERNELS"

KIN_NAME = config.KIN_NAME[KIN_VAR]
LAG_TAG  = config.lag_tag(config.CROSS_LAG_RANGE_MS, config.TRIAL_WINDOW_MS)
REGIONS_PLOT = ["CP", "CB"] if DATA_SOURCE == "real" else ["SIM"]
SESS_LIST = runs.session_list(DATA_SOURCE, SESSIONS, MIN_STRATUM, REACH_CAT, SESSION)
SETTING_COLS = ["hist_mode", "kin_hist_mode"]
SETTINGS = [(hm, km) for hm in HIST_MODES for km in KIN_HIST_MODES]
SPECS = {v: kernels.make_spec(bin_ms=BIN_MS, h_list=H_LIST, hist_mode=v[0], h_fixed=150,
                              h_kin=H_KIN, h_kin_list=H_KIN_LIST, kin_hist_mode=v[1],
                              scale_bases=SCALE_BASES, fit_from_ms=FIT_FROM_MS,
                              n_cross=N_CROSS_BASES, cross_range=config.CROSS_LAG_RANGE_MS,
                              spacing=BASIS_SPACING, event=EVENT_KERNEL,
                              nuisance=TASK_MEAN_NUISANCE, ridge=RIDGE_ALPHA, n_folds=N_FOLDS,
                              adequacy_ms=ADEQUACY_MS, seed=SHUFFLE_SEED) for v in SETTINGS}
PAD_MS = max(kernels.pad_needed_ms(s) for s in SPECS.values())
if KIN_SOURCE == "short" and any(max(sp["h_kin_list"]) > 150 for sp in SPECS.values()):
    print("WARNING: KIN_SOURCE = 'short' with kinematic history > 150 ms -> fitted rows start later "
          "(no kinematics before -375 ms). Use KIN_SOURCE = 'long'.")
_HL = {"cv": f"CV {list(H_LIST)}", "fixed": "150 ms", "max": f"max {max(H_LIST)} ms"}
_KL = {"cv": f"CV {list(H_KIN_LIST)}", "fixed": f"{H_KIN} ms", "max": f"max {max(H_KIN_LIST)} ms"}
def _label(v):
    return f"neuron history {_HL[v[0]]} | kin history {_KL[v[1]]}"

# RUN_NAME: <data>_kernels [_hist<modes>] [_kinLong] [_kinH<modes>] [_ownBases6] [_bases..] [_event] [_nuisance]
_tok = lambda m, fixed, top: {"cv": "CV", "fixed": f"{fixed}", "max": f"max{top}"}[m]
_extra = (("_hist" + "+".join(_tok(h, 150, max(H_LIST)) for h in HIST_MODES)
           if HIST_MODES != ("cv",) else "")
          + ("_kinLong" if KIN_SOURCE == "long" else "")
          + ("_kinH" + "+".join(_tok(h, H_KIN, max(H_KIN_LIST)) for h in KIN_HIST_MODES)
             if KIN_HIST_MODES != ("fixed",) or H_KIN != 150 else "")
          + ("" if SCALE_BASES else "_ownBases6")
          + (f"_from{FIT_FROM_MS:g}" if FIT_FROM_MS is not None else "")
          + (f"_bases{BASIS_SPACING.capitalize()}{N_CROSS_BASES}"
             if (BASIS_SPACING, N_CROSS_BASES) != ("log", 6) else "")
          + ("_event" if EVENT_KERNEL else "") + ("_nuisance" if TASK_MEAN_NUISANCE else "")
          + (f"_kin{KIN_VAR}" if KIN_VAR != config.KIN_VAR else ""))
RUN_NAME = f"{runs.data_name(DATA_SOURCE, SESS_LIST, SIM_SCENARIO, SIM_PARAMS, SIM_SEEDS)}_kernels{_extra}"
FIG_DIR  = config.fig_dir(LAG_TAG, OUT_SUBDIR, RUN_NAME) if SAVE_FIGS else None
config.warn_long(config.OUT_DIR / LAG_TAG / OUT_SUBDIR / RUN_NAME / "8_effect_strength_fixed.png")
def _fig(base, *v):
    if len(v) == 1 and isinstance(v[0], tuple):
        v = v[0]
    return base if len(SETTINGS) == 1 or not v else f"{base}_{v[0]}-{v[1]}"

RUN_INFO = "\n".join(
    [f"run: {RUN_NAME}", "script: 03_ppglm_kernels.py (Poisson GLM with smooth causal kernels; MAIN)"]
    + runs.common_info(DATA_SOURCE, SESS_LIST, REACH_CAT, MIN_STRATUM, MIN_RATE_HZ, SIM_SCENARIO,
                       SIM_PARAMS, SIM_SEEDS, SESSION, N_SHUFFLE, FDR_Q, kin_var=KIN_VAR,
                       kin_source=KIN_SOURCE)
    + [f"bin {BIN_MS:g} ms raw counts, raw {KIN_NAME} ({KIN_VAR}); neuron own history modes "
       f"{list(HIST_MODES)} (candidates {list(H_LIST)} ms); kinematics own history modes "
       f"{list(KIN_HIST_MODES)} (candidates {list(H_KIN_LIST)} ms, fixed {H_KIN} ms); "
       f"own-history log bumps {'scaled with H' if SCALE_BASES else '6 for every H'}",
       f"cross bases: {N_CROSS_BASES} {BASIS_SPACING} over {config.CROSS_LAG_RANGE_MS} ms; "
       f"event kernel {EVENT_KERNEL}; task-mean nuisance {TASK_MEAN_NUISANCE}; ridge {RIDGE_ALPHA:g}",
       f"CV {N_FOLDS} trial folds (stratified by spout); adequacy test over (H, H+{ADEQUACY_MS}] ms, "
       f"flag p < {ADEQ_P}; lag bootstrap {N_BOOT} (significant neurons)",
       f"adaptive shuffles: {ADAPTIVE_SHUFFLE} (cap {ADAPTIVE_CAP}, batch {EXTRA_BATCH}, h {SEQ_H})"])
runs.write_run_info(RUN_INFO, FIG_DIR, LAG_TAG, OUT_SUBDIR, RUN_NAME, SAVE_FIGS, SAVE_VARS)


# %% ── 1. data: real sessions, or one simulated session per seed ────────────
# simulations must start early enough for the longest own history + adequacy window
SESSIONS_DATA, NT_ALL = runs.load_data(DATA_SOURCE, SESS_LIST, REACH_CAT, SESSION, SIM_SCENARIO,
                                       dict(SIM_PARAMS, hist_pad_ms=PAD_MS + 100), SIM_SEEDS,
                                       MIN_RATE_HZ, kin_var=KIN_VAR, kin_source=KIN_SOURCE)


# %% ── 2. bin (raw counts, raw kinematics) ──────────────────────────────────
BINNED = {tag: pp.bin_session(S, BIN_MS, hist_pad_ms=PAD_MS, kin_smooth="raw", fr_kernel=None)
          for tag, S in SESSIONS_DATA}
_B0 = next(iter(BINNED.values()))
_D0 = kernels.build_design(_B0["Y"][0], _B0["X"], _B0["in_window"], _B0["cond"], _B0["t"],
                           SPECS[SETTINGS[0]])
print(f"{len(_B0['t'])} bins of {BIN_MS:g} ms; fitted rows per trial {int(_D0['rows'].sum())} "
      f"({_B0['t'][_D0['rows']][0]*1000:.1f}..{_B0['t'][_D0['rows']][-1]*1000:.1f} ms); "
      f"cross bases centred at {np.round(SPECS[SETTINGS[0]]['cross_centers']).astype(int).tolist()} ms")


# %% ── 3. fit every neuron (both directions, shuffle null) -- cached per session
def _fit_session(tag, S, v):
    B = BINNED[tag]
    perms = nulls.make_perms(B["cond"], N_SHUFFLE, seed=SHUFFLE_SEED)
    idx = np.flatnonzero(S["_keep"])
    res = Parallel(n_jobs=N_JOBS, backend="loky")(
        delayed(kernels.kernel_neuron)(B["Y"][i], B["X"], B["in_window"], B["cond"], B["t"],
                                       SPECS[v], perms) for i in idx)
    for r in res:
        for d in ("xy", "yx"):
            r[d]["null"] = r[d]["null"].astype(np.float32)
    return dict(idx=idx, res=res)

FITS, FIT_PARAMS = {}, {}
t0 = time.time()
for si, (tag, S) in enumerate(SESSIONS_DATA):
    stamp = io.mi_stamp() if S["source"] == "real" else dict(SIM_PARAMS, scenario=SIM_SCENARIO, tag=tag)
    if KIN_SOURCE == "long":
        stamp = dict(stamp, reach_long_mtime=os.path.getmtime(cache._long(io.reach_long_path())))
    for v in SETTINGS:
        sp = {k: val for k, val in SPECS[v].items() if not k.startswith(("B_", "L_"))}
        params = dict(stamp=stamp, session=S["key"], reach_cat=REACH_CAT,
                      uids=[n["uid"] for n in S["neurons"]], keep=S["_keep"], spec=sp,
                      n_shuffle=N_SHUFFLE, shuffle_seed=SHUFFLE_SEED, window=config.TRIAL_WINDOW_MS,
                      shuffle_within=config.SHUFFLE_WITHIN, kin_var=KIN_VAR, kin_source=KIN_SOURCE,
                      pad=PAD_MS)
        FIT_PARAMS[(tag, v)] = params
        FITS[(tag, v)] = cache.cached(f"03ker_{tag}_{v[0]}-{v[1]}", params, lambda: _fit_session(tag, S, v),
                                      mode=CACHE_MODE, save=SAVE_CACHE)
    print(f"  session {si+1}/{len(SESSIONS_DATA)} ({tag}) done, {time.time()-t0:.0f}s elapsed", flush=True)
print(f"fits ready in {time.time()-t0:.0f}s")


# %% ── 4. per-neuron table, adaptive shuffles, BH-FDR, summaries ───────────
rows, KERN = [], {}
for tag, S in SESSIONS_DATA:
    for v in SETTINGS:
        F = FITS[(tag, v)]
        for i, r in zip(F["idx"], F["res"]):
            xy, yx = r["xy"], r["yx"]
            row = dict(runs.base_row(tag, S, i), hist_mode=v[0], kin_hist_mode=v[1], H=r["H"],
                       H_kin=r["H_kin"], n_spikes=r["n_spikes"], rows_per_trial=r["rows_per_trial"],
                       pS_xy=xy["p_shuf"], z_xy=xy["z"], score_xy=xy["score"], lr_xy=xy["lr"],
                       pChi2_xy=xy["p_chi2"], cv_gain_xy=xy["cv_gain_bits"],
                       adequacy_p_xy=xy["adequacy_p"], dispersion=xy["dispersion"],
                       converged=xy["converged"],
                       pS_yx=yx["p_shuf"], z_yx=yx["z"], gc_yx=yx["gc"], pF_yx=yx["p_F"],
                       cv_gain_yx=100.0 * yx["cv_gain_r2"], adequacy_p_yx=yx["adequacy_p"],
                       # excess over the shuffle (observed and null = the SAME statistic)
                       bpsx_xy=models.excess_bits_per_spike(xy["score"], xy["null"], r["n_spikes"]),
                       # the same BEFORE subtracting the shuffle mean (figure 3b), and that mean
                       bps_xy=xy["score"] / (2 * max(r["n_spikes"], 1.0) * np.log(2)),
                       bps0_xy=float(np.mean(xy["null"])) / (2 * max(r["n_spikes"], 1.0) * np.log(2))
                               if np.size(xy["null"]) else np.nan,
                       pr2_yx=100.0 * (1 - np.exp(-yx["gc"])),
                       pr20_yx=100.0 * float(np.mean(1 - np.exp(-np.asarray(yx["null"], float))))
                               if np.size(yx["null"]) else np.nan,
                       pr2x_yx=models.excess_partial_r2(yx["gc"], yx["null"]))
            for d in ("xy", "yx"):
                row.update({f"lag_peak_{d}": r[d]["peak_ms"], f"lag_sign_{d}": r[d]["peak_sign"],
                            f"lag_centroid_{d}": r[d]["centroid_ms"]})
                KERN[(tag, v, S["neurons"][i]["uid"], d)] = dict(tau_ms=r[d]["tau_ms"], beta=r[d]["beta"])
            rows.append(row)
RES = pd.DataFrame(rows)
for d in ("xy", "yx"):
    RES[f"adequacy_flag_{d}"] = RES[f"adequacy_p_{d}"] < ADEQ_P

if ADAPTIVE_SHUFFLE:
    RES = runs.adaptive_extend(
        RES, SETTING_COLS, {"xy": "score_xy", "yx": "gc_yx"}, N_SHUFFLE, FDR_Q, ADAPTIVE_CAP,
        SEQ_H, EXTRA_BATCH, SHUFFLE_SEED,
        cond_of=lambda tag, v: BINNED[tag]["cond"],
        null_of=lambda row, d: kernels.KernelNull(
            BINNED[row.tag]["Y"][row.neuron_idx], BINNED[row.tag]["X"], BINNED[row.tag]["in_window"],
            BINNED[row.tag]["cond"], BINNED[row.tag]["t"], SPECS[(row.hist_mode, row.kin_hist_mode)],
            row.H if d == "xy" else row.H_kin, d),
        cache_params=lambda tag, v: FIT_PARAMS[(tag, tuple(v))], n_jobs=N_JOBS, cache_name="03ker",
        cache_mode=CACHE_MODE, save_cache=SAVE_CACHE, label=_label)

RES = runs.fdr_and_classes(RES, SETTING_COLS, FDR_Q)
runs.print_direction_summary(RES, SETTING_COLS, FDR_Q, len(SESSIONS_DATA), _label)
FLOOR = runs.print_floor_summary(RES, SETTING_COLS, N_SHUFFLE, FDR_Q, _label)
EFFECT_COLS = {"xy": ("bpsx_xy", "excess over shuffle (bits/spike)"),
               "yx": ("pr2x_yx", f"excess over shuffle, partial R² of {KIN_NAME} (%)")}
CVGAIN_COLS = {"xy": ("cv_gain_xy", "CV gain, held-out (bits/spike)"),
               "yx": ("cv_gain_yx", f"CV gain, held-out {KIN_NAME} variance (%)")}
print("\n(secondary) CV GAIN = total predictive gain, INCLUDING the shared reach shape:")
runs.print_effect_summary(RES, SETTING_COLS, REGIONS_PLOT, CVGAIN_COLS, _label)
runs.print_effect_summary(RES, SETTING_COLS, REGIONS_PLOT, EFFECT_COLS, _label)
print("\n--- same, WITHOUT neurons > 100 Hz (Poisson effect sizes are mis-scaled there) ---")
runs.print_effect_summary(RES[RES.rate_hz <= 100], SETTING_COLS, REGIONS_PLOT, EFFECT_COLS, _label)

print(f"\n=== OWN HISTORY: chosen H per region, and ADEQUACY flags (p < {ADEQ_P}) ===")
for v in runs.settings_of(RES, SETTING_COLS):
    R = runs.select(RES, SETTING_COLS, v)
    print(f"[{_label(v)}]")
    print("  neuron model, chosen H (ms), fraction of neurons:")
    print("  " + (pd.crosstab(R.region, R.H, normalize="index").round(2).to_string().replace("\n", "\n  ")))
    hk = R.groupby("tag").H_kin.first()
    print(f"  kinematics model, chosen H_kin (ms), number of sessions: {hk.value_counts().sort_index().to_dict()}"
          f"  | fitted rows per trial: {sorted(R.rows_per_trial.unique().tolist())}")
    for d in ("xy", "yx"):
        fl = R.groupby("region")[f"adequacy_flag_{d}"].mean().round(2).to_dict()
        print(f"  adequacy flag, {plotting.DIR_LABEL[d]}: {fl}"
              + ("   (kinematics model: same for every neuron of a session)" if d == "yx" else ""))
    ok = R[~R.adequacy_flag_xy]
    if {"CP", "CB"} <= set(R.region) and len(ok) < len(R):
        print("  direction classes WITHOUT neurons flagged in the neuron model:")
        print("  " + (pd.crosstab(ok.region, ok.direction, normalize="index").round(2)
                      .to_string().replace("\n", "\n  ")))

if len(SETTINGS) > 1:                           # history-sensitive neurons
    _n = RES.groupby("uid")[["sig_xy", "sig_yx"]].nunique()
    flip = (_n > 1).any(axis=1)
    RES = RES.merge(flip.rename("history_sensitive").reset_index(), on="uid")
    print(f"\nhistory-sensitive neurons (significance differs between the {len(SETTINGS)} history "
          f"settings): {int(flip.sum())}/{len(flip)}  "
          f"(kin -> neuron: {int((_n.sig_xy > 1).sum())}, neuron -> kin: {int((_n.sig_yx > 1).sum())})")

print(f"\n=== {kernels.KERNEL_LAG_NOTE} ===")
runs.print_lag_summary(RES, SETTING_COLS, REGIONS_PLOT, ("centroid", "peak"), "KERNEL LAG", _label,
                       unreliable=("peak",))

if DATA_SOURCE == "sim":
    DET = runs.sim_truth_check(RES, SETTING_COLS, _label, lag_keys=("centroid", "peak"),
                               true_lag=SIM_PARAMS["lag_ms"])
    print("\nSECONDARY tests (p < 0.05), same neurons: chi2 LR (kin->neuron), F (neuron->kin) vs shuffle"
          "\n  (expected False = false-positive rate; the chi2/F rows are expected to be miscalibrated)")
    for v in runs.settings_of(RES, SETTING_COLS):
        R = runs.select(RES, SETTING_COLS, v)
        for d, e, pc, nm in (("xy", "exp_xy", "pChi2_xy", "chi2 LR"), ("yx", "exp_yx", "pF_yx", "F    ")):
            for exp in (False, True):
                g = R[R[e] == exp]
                if len(g):
                    print(f"  [{_label(v)}] {plotting.DIR_LABEL[d]} expected {str(exp):5s} n={len(g):4d}: "
                          f"{nm} {np.mean(g[pc] < 0.05):.3f}   shuffle {np.mean(g[f'pS_{d}'] < 0.05):.3f}")


# %% ── 4b. lag bootstrap for FDR-significant neurons (cached per session) ──
def _boot_session(tag, v, todo):
    B = BINNED[tag]
    out = Parallel(n_jobs=N_JOBS, backend="loky")(
        delayed(kernels.kernel_lag_bootstrap)(B["Y"][ni], B["X"], B["in_window"], B["cond"], B["t"],
                                              SPECS[v], H, d, N_BOOT, SHUFFLE_SEED)
        for ni, H, d, _ in todo)
    return {(u, d): o for (_, _, d, u), o in zip(todo, out)}

BOOT = {}
for tag, _ in SESSIONS_DATA:
    for v in SETTINGS:
        R = runs.select(RES[RES.tag == tag], SETTING_COLS, v)
        todo = [(int(row.neuron_idx), int(row.H if d == "xy" else row.H_kin), d, row.uid)
                for row in R.itertuples() for d in ("xy", "yx") if getattr(row, f"sig_{d}")]
        if not todo:
            continue
        BOOT.update({(tag, v) + k: o for k, o in cache.cached(
            f"03boot_{tag}_{v[0]}-{v[1]}", dict(fit=FIT_PARAMS[(tag, v)], n_boot=N_BOOT, todo=todo),
            lambda: _boot_session(tag, v, todo), mode=CACHE_MODE, save=SAVE_CACHE).items()})
for key in ("peak", "centroid"):
    for d in ("xy", "yx"):
        RES[f"lag_{key}_lo_{d}"] = [BOOT.get((r.tag, (r.hist_mode, r.kin_hist_mode), r.uid, d), {}).get(f"{key}_ci", (np.nan,))[0]
                                    for r in RES.itertuples()]
        RES[f"lag_{key}_hi_{d}"] = [BOOT.get((r.tag, (r.hist_mode, r.kin_hist_mode), r.uid, d), {}).get(f"{key}_ci", (np.nan, np.nan))[1]
                                    for r in RES.itertuples()]
print(f"lag bootstrap: {len(BOOT)} neuron x direction kernels with 95% CIs")

if SAVE_VARS:
    runs.save_tables(dict(NT=NT_ALL, kernels=KERN, boot=BOOT, lag_note=kernels.KERNEL_LAG_NOTE,
                          params=dict(DATA_SOURCE=DATA_SOURCE, SESS_LIST=SESS_LIST,
                                      SIM_SCENARIO=SIM_SCENARIO, SIM_PARAMS=SIM_PARAMS,
                                      SIM_SEEDS=SIM_SEEDS, KIN_VAR=KIN_VAR, KIN_SOURCE=KIN_SOURCE,
                                      HIST_MODES=HIST_MODES, H_LIST=H_LIST, KIN_HIST_MODES=KIN_HIST_MODES,
                                      H_KIN_LIST=H_KIN_LIST, H_KIN=H_KIN, SCALE_BASES=SCALE_BASES,
                                      N_CROSS_BASES=N_CROSS_BASES,
                                      BASIS_SPACING=BASIS_SPACING, EVENT_KERNEL=EVENT_KERNEL,
                                      TASK_MEAN_NUISANCE=TASK_MEAN_NUISANCE, RIDGE_ALPHA=RIDGE_ALPHA,
                                      N_FOLDS=N_FOLDS, N_SHUFFLE=N_SHUFFLE, FDR_Q=FDR_Q,
                                      N_BOOT=N_BOOT, LAG_TAG=LAG_TAG)),
                     RES, LAG_TAG, OUT_SUBDIR, RUN_NAME)


# %% ── 5a/5b. direction classes (CP vs CB) and direction map per region ─────
runs.plot_direction_figures(RES, SETTING_COLS, REGIONS_PLOT, N_SHUFFLE, FDR_Q, FIG_DIR, _fig,
                            SAVE_FIGS, SHOW_FIGS, f"03 KERNELS  |  {RUN_NAME}", _label)


# %% ── 5c. CV gain per neuron: FDR-significant (top) vs all (bottom) ────────
# Held-out gain of the full over the reduced model: ~<= 0 by chance (0 = dotted line).
# Symmetric-log axis (negative values exist). Compare regions WITHIN a panel only.
for v in SETTINGS:
    R = runs.select(RES, SETTING_COLS, v)
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), squeeze=False)
    for c_, d in enumerate(("xy", "yx")):
        col, lab = EFFECT_COLS[d]
        va = R[col].dropna().to_numpy()
        s = float(np.median(np.abs(va))) if va.size else 1.0
        s = s if s > 0 else 1e-6
        edges = s * np.sinh(np.linspace(np.arcsinh(min(va.min() if va.size else -s, -s) / s),
                                        np.arcsinh(max(va.max() if va.size else s, s) / s), 41))
        for r_, (nm, only_sig) in enumerate((("FDR-significant", True), ("all neurons", False))):
            ax = axs[r_, c_]
            for reg in REGIONS_PLOT:
                g = R[R.region == reg]
                g = g[g[f"sig_{d}"]] if only_sig else g
                x_ = g[col].dropna().to_numpy()
                if not x_.size:
                    continue
                c = config.CLR.get(reg, (0.4, 0.4, 0.4)); w = np.full(x_.size, 1 / x_.size)
                ax.hist(x_, bins=edges, weights=w, color=c, alpha=0.5,
                        label=f"{reg} (n={x_.size}, median {np.median(x_):.2g})")
                ax.hist(x_, bins=edges, weights=w, color=c, histtype="step", lw=1.2)
                ax.axvline(np.median(x_), color=c, ls="--", lw=1.4)
            ax.axvline(0, color="k", ls=":", lw=0.9)
            ax.set_xscale("symlog", linthresh=s)
            ax.set(xlim=(edges[0], edges[-1]), xlabel=lab, ylabel="fraction of neurons",
                   title=f"{plotting.DIR_LABEL[d]} — {nm}")
            ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"03 CV gain (held-out; SECONDARY: includes the shared reach shape of the spout side, "
                 f"so > 0 even without trial-specific coupling)  |  {_label(v)}", fontsize=9)
    plotting.finish(fig, FIG_DIR, _fig("3_cv_gain", v), SAVE_FIGS, SHOW_FIGS)


# %% ── 5c2. effect BEFORE subtracting the shuffle mean (log scale) ──────────
# Same layout as 3 (top: FDR-significant, bottom: all). In-sample effect, always > 0, so a
# log axis works: kin->neuron bits/spike = score / (2 n_spikes ln2); neuron->kin partial R2.
# CAUTION: NOT chance-corrected -- it includes the chance level (extra parameters, finite
# trials), which differs between neurons (rate, n spikes). Grey outline = each neuron's own
# chance level (mean of its shuffle null). The primary effect size is the EXCESS (figure 8).
RAW_COLS = {"xy": ("bps_xy", "bps0_xy", "in-sample gain from kinematics (bits/spike)"),
            "yx": ("pr2_yx", "pr20_yx", f"in-sample partial R² of {KIN_NAME} (%)")}
if "bps_xy" not in RES:                       # RES from before these columns existed: add them from FITS
    _raw = {}
    for (tag, v), F in FITS.items():
        S_ = dict(SESSIONS_DATA)[tag]
        for i, r in zip(F["idx"], F["res"]):
            k = 2 * max(r["n_spikes"], 1.0) * np.log(2)
            _raw[(tag, v, S_["neurons"][i]["uid"])] = dict(
                bps_xy=r["xy"]["score"] / k, bps0_xy=float(np.mean(r["xy"]["null"])) / k,
                pr2_yx=100.0 * (1 - np.exp(-r["yx"]["gc"])),
                pr20_yx=100.0 * float(np.mean(1 - np.exp(-np.asarray(r["yx"]["null"], float)))))
    _keys = [(r.tag, (r.hist_mode, r.kin_hist_mode), r.uid) for r in RES.itertuples()]
    for c in ("bps_xy", "bps0_xy", "pr2_yx", "pr20_yx"):
        RES[c] = [_raw[k][c] for k in _keys]
for v in SETTINGS:
    R = runs.select(RES, SETTING_COLS, v)
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), squeeze=False)
    for c_, d in enumerate(("xy", "yx")):
        col, col0, lab = RAW_COLS[d]
        allv = np.concatenate([R[col].to_numpy(float), R[col0].to_numpy(float)])
        allv = allv[np.isfinite(allv) & (allv > 0)]
        if not allv.size:
            continue
        edges = np.logspace(np.log10(np.percentile(allv, 0.5)), np.log10(allv.max()), 41)
        for r_, (nm, only_sig) in enumerate((("FDR-significant", True), ("all neurons", False))):
            ax = axs[r_, c_]
            for reg in REGIONS_PLOT:
                g = R[R.region == reg]
                g = g[g[f"sig_{d}"]] if only_sig else g
                x_ = g[col].to_numpy(float); x0 = g[col0].to_numpy(float)
                x_, x0 = x_[np.isfinite(x_) & (x_ > 0)], x0[np.isfinite(x0) & (x0 > 0)]
                if not x_.size:
                    continue
                c = config.CLR.get(reg, (0.4, 0.4, 0.4))
                ax.hist(x_, bins=edges, weights=np.full(x_.size, 1 / x_.size), color=c, alpha=0.5,
                        label=f"{reg} (n={x_.size}, median {np.median(x_):.2g})")
                ax.axvline(np.median(x_), color=c, ls="--", lw=1.4)
                if x0.size:
                    ax.hist(x0, bins=edges, weights=np.full(x0.size, 1 / x0.size), histtype="step",
                            color=c, lw=1.0, ls=":", label=f"{reg} chance level (shuffle mean), median {np.median(x0):.2g}")
            ax.set_xscale("log")
            ax.set(xlabel=lab, ylabel="fraction of neurons", title=f"{plotting.DIR_LABEL[d]} — {nm}")
            ax.legend(fontsize=6.5, frameon=False)
    fig.suptitle(f"03 effect BEFORE subtracting the shuffle mean (log scale)  |  {_label(v)}\n"
                 "CAUTION: not chance-corrected -- includes each neuron's chance level (dotted), which varies with "
                 "rate and spike count; use the EXCESS (figure 8) to compare strengths", fontsize=8.5, color="tab:red")
    plotting.finish(fig, FIG_DIR, _fig("3b_effect_raw_log", v), SAVE_FIGS, SHOW_FIGS)


# %% ── 5d. observed vs the trial-shuffle null, per region ────────────────────
for v in SETTINGS:
    PAIRS = {}
    for tag, S in SESSIONS_DATA:
        F = FITS[(tag, v)]
        for i, r in zip(F["idx"], F["res"]):
            reg, ns = S["neurons"][i]["region"], max(r["n_spikes"], 1.0)
            k = 2 * ns * np.log(2)                     # score -> bits/spike (as in 02)
            PAIRS.setdefault((reg, "xy"), []).append(
                (r["xy"]["score"] / k, np.asarray(r["xy"]["null"], float) / k))
            PAIRS.setdefault((reg, "yx"), []).append(
                (1 - np.exp(-r["yx"]["gc"]), 1 - np.exp(-np.asarray(r["yx"]["null"], float))))
    runs.plot_null_vs_observed(PAIRS, REGIONS_PLOT, {"xy": "in-sample gain from kinematics (bits/spike)",
                                                     "yx": "in-sample partial R² from spikes"},
                               FIG_DIR, _fig("3d_null_vs_observed", v), SAVE_FIGS, SHOW_FIGS,
                               f"03 KERNELS  |  {_label(v)}")


# %% ── 5e. lag distributions (FDR-significant), peak and centroid ──────────
for v in SETTINGS:
    R = runs.select(RES, SETTING_COLS, v)
    fig, axs = plt.subplots(2, 2, figsize=(11, 7))
    bins = np.linspace(0, config.CROSS_LAG_RANGE_MS[1], 16)
    for c_, d in enumerate(("xy", "yx")):
        for r_, key in enumerate(("centroid", "peak")):
            plotting.lag_hist(axs[r_, c_], R, d, key, REGIONS_PLOT, bins)
            axs[r_, c_].set_title(f"{plotting.DIR_LABEL[d]} — {key}"
                                  + ("  (UNRELIABLE, see note)" if key == "peak" else "  (use this)"),
                                  fontsize=9, color="tab:red" if key == "peak" else "k")
    fig.suptitle(f"03 kernel lags of FDR-significant neurons  |  {_label(v)}\n{kernels.KERNEL_LAG_NOTE}",
                 fontsize=8)
    plotting.finish(fig, FIG_DIR, _fig("4_lags", v), SAVE_FIGS, SHOW_FIGS)


# %% ── 5e2. lag of each FDR-significant neuron with its 95% trial-bootstrap CI ─
# One row per neuron (sorted by lag, regions stacked); bar = 95% CI from N_BOOT trial
# resamples (refitting the full model). Wide bars = the lag of that neuron is not pinned down.
for v in SETTINGS:
    R = runs.select(RES, SETTING_COLS, v)
    fig, axs = plt.subplots(2, 2, figsize=(12, 9), squeeze=False)
    for r_, d in enumerate(("xy", "yx")):
        for c_, key in enumerate(("centroid", "peak")):
            ax, y0, yt, yl = axs[r_, c_], 0, [], []
            for reg in REGIONS_PLOT:
                g = R[(R.region == reg) & R[f"sig_{d}"]].dropna(subset=[f"lag_{key}_{d}"])
                g = g.sort_values(f"lag_{key}_{d}")
                if not len(g):
                    continue
                c = config.CLR.get(reg, (0.4, 0.4, 0.4))
                m, lo, hi = (g[f"lag_{key}_{d}"].to_numpy(), g[f"lag_{key}_lo_{d}"].to_numpy(),
                             g[f"lag_{key}_hi_{d}"].to_numpy())
                yy = y0 + np.arange(len(g))
                ok = np.isfinite(lo) & np.isfinite(hi)
                w = np.median(hi[ok] - lo[ok]) if ok.any() else np.nan
                ax.hlines(yy[ok], lo[ok], hi[ok], color=c, lw=0.8, alpha=0.6)
                ax.plot(m, yy, "o", ms=2.5, color=c,
                        label=f"{reg} n={len(g)}, median lag {np.median(m):.0f} ms, median CI width {w:.0f} ms")
                yt.append(y0 + len(g) / 2); yl.append(reg)
                y0 += len(g) + max(3, len(g) // 10)
            if DATA_SOURCE == "sim":
                ax.axvline(SIM_PARAMS["lag_ms"], color="k", ls="--", lw=1, label="true lag")
            ax.set_yticks(yt, yl)
            ax.set(xlim=(0, config.CROSS_LAG_RANGE_MS[1]), xlabel=f"{key} lag (ms)",
                   title=f"{plotting.DIR_LABEL[d]} — {key} (95% CI, {N_BOOT} trial bootstraps)"
                         + ("\nUNRELIABLE, see note" if key == "peak" else ""))
            if key == "peak":
                ax.title.set_color("tab:red")
            ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle(f"03 lag per FDR-significant neuron with bootstrap CI  |  {_label(v)}\n{kernels.KERNEL_LAG_NOTE}",
                 fontsize=8)
    plotting.finish(fig, FIG_DIR, _fig("4b_lag_CIs", v), SAVE_FIGS, SHOW_FIGS)


# %% ── 5e3. centroid lag distribution, COUNTS of FDR-significant neurons ─────
# One panel per direction, CP and CB in the same panel: in each 10 ms bin the regions'
# bars sit side by side. Each bar is split into neurons significant in THIS direction only
# (solid) and in BOTH directions (hatched). Same count scale in both panels.
import textwrap
for v in SETTINGS:
    R = runs.select(RES, SETTING_COLS, v)
    bw = 10.0
    edges = np.arange(0, config.CROSS_LAG_RANGE_MS[1] + bw, bw)
    nreg = len(REGIONS_PLOT)
    w = 0.85 * bw / nreg
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
    for ax, d in zip(axs, ("xy", "yx")):
        for j, reg in enumerate(REGIONS_PLOT):
            g = R[(R.region == reg) & R[f"sig_{d}"]]
            both = (g.sig_xy & g.sig_yx).to_numpy()
            lag = g[f"lag_centroid_{d}"].to_numpy(float)
            n_only, _ = np.histogram(lag[~both & np.isfinite(lag)], edges)
            n_both, _ = np.histogram(lag[both & np.isfinite(lag)], edges)
            x = edges[:-1] + 0.075 * bw + (j + 0.5) * w
            c = config.CLR.get(reg, (0.4, 0.4, 0.4))
            n_tot = int((R.region == reg).sum())
            ax.bar(x, n_only, width=w, color=c, edgecolor="w", lw=0.5,
                   label=f"{reg}: {len(g)}/{n_tot} significant — {d} only (n={int(n_only.sum())})")
            ax.bar(x, n_both, bottom=n_only, width=w, color=c, alpha=0.45, hatch="///",
                   edgecolor="w", lw=0.5, label=f"{reg}: both directions (n={int(n_both.sum())})")
            if np.isfinite(lag).any():
                ax.axvline(np.nanmedian(lag), color=c, ls="--", lw=1.5,
                           label=f"{reg} median {np.nanmedian(lag):.0f} ms")
        ax.set_xticks(edges[::2])
        ax.set(xlim=(edges[0], edges[-1]), xlabel="centroid lag (ms)", title=plotting.DIR_LABEL[d],
               ylabel="# significant neurons" if d == "xy" else None)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"03 centroid lag of FDR-significant neurons (counts)  |  {_label(v)}\n"
                 + textwrap.fill(kernels.KERNEL_LAG_NOTE, 150), fontsize=8)
    plotting.finish(fig, FIG_DIR, _fig("4c_centroid_lag_counts", v), SAVE_FIGS, SHOW_FIGS)


# %% ── 5f. mean kernel per region (each neuron scaled to max |beta| = 1) ───
for v in SETTINGS:
    R = runs.select(RES, SETTING_COLS, v)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, d in zip(axs, ("xy", "yx")):
        for reg in REGIONS_PLOT:
            c = config.CLR.get(reg, (0.4, 0.4, 0.4))
            for sel, ls, lab in (((R.region == reg) & R[f"sig_{d}"], "-", "significant"),
                                 ((R.region == reg), "--", "all")):
                K = [KERN[(r.tag, v, r.uid, d)]["beta"] for r in R[sel].itertuples()]
                if not K:
                    continue
                K = np.array(K); K = K / np.maximum(np.abs(K).max(1, keepdims=True), 1e-12)
                tau = KERN[(R[sel].tag.iloc[0], v, R[sel].uid.iloc[0], d)]["tau_ms"]
                m = K.mean(0); se = K.std(0, ddof=1) / np.sqrt(len(K)) if len(K) > 1 else 0 * m
                ax.plot(tau, m, color=c, ls=ls, lw=1.8 if ls == "-" else 1, label=f"{reg} {lab} (n={len(K)})")
                if ls == "-":
                    ax.fill_between(tau, m - se, m + se, color=c, alpha=0.2)
        ax.axhline(0, color="k", lw=0.4)
        ax.set(xlabel="lag (ms)", ylabel="mean normalised β(τ) (± SEM)", title=plotting.DIR_LABEL[d])
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"03 mean kernel per region  |  {_label(v)}  (biphasic = the neuron follows a CHANGE, "
                 f"e.g. velocity, not a delay)", fontsize=9)
    plotting.finish(fig, FIG_DIR, _fig("5_mean_kernel", v), SAVE_FIGS, SHOW_FIGS)


# %% ── 5g. own-history choice and adequacy flags, per region ───────────────
for v in SETTINGS:
    R = runs.select(RES, SETTING_COLS, v)
    fig, axs = plt.subplots(1, 3, figsize=(15, 3.8))
    tab = pd.crosstab(R.region, R.H, normalize="index").reindex(REGIONS_PLOT)
    tab.T.plot(kind="bar", ax=axs[0], color=[config.CLR.get(r, (0.4, 0.4, 0.4)) for r in tab.index], rot=0)
    axs[0].set(xlabel="chosen own history H (ms)", ylabel="fraction of neurons",
               title=f"neuron model ({_HL[v[0]]})")
    hk = R.groupby("tag").H_kin.first().value_counts().reindex(sorted(set(SPECS[v]["h_kin_list"])), fill_value=0)
    axs[2].bar([str(h) for h in hk.index], hk.values, color="0.5")
    axs[2].set(xlabel="chosen kinematic own history H_kin (ms)", ylabel="# sessions",
               title=f"kinematics model ({_KL[v[1]]})")
    fl = pd.DataFrame({plotting.DIR_LABEL[d]: R.groupby("region")[f"adequacy_flag_{d}"].mean()
                       for d in ("xy", "yx")}).reindex(REGIONS_PLOT)
    fl.T.plot(kind="bar", ax=axs[1], color=[config.CLR.get(r, (0.4, 0.4, 0.4)) for r in fl.index], rot=0)
    axs[1].set(ylabel=f"fraction flagged (p < {ADEQ_P})", ylim=(0, 1.05),
               title="adequacy: older own history still predicts?")
    fig.suptitle(f"03 own history  |  {_label(v)}", fontsize=9)
    plotting.finish(fig, FIG_DIR, _fig("6_own_history", v), SAVE_FIGS, SHOW_FIGS)


# %% ── 5h. effect strength and fraction of neurons with an effect (pi1) ────
runs.plot_effect_figures(RES, SETTING_COLS, REGIONS_PLOT, EFFECT_COLS, FIG_DIR, _fig,
                         SAVE_FIGS, SHOW_FIGS, f"03 KERNELS  |  {RUN_NAME}", _label, rate_flag_hz=100)


# %% ── 5h2. fraction of neurons with an effect: pi1 + session CI, CP vs CB ───
# Grey band = pi1 of a simulated NULL run of this script (no neuron has an effect): real pi1
# should be read against it. PI1_NULL_RUN = run name of a saved null simulation (e.g. the
# shared-drive null), or None for no band.
PI1_NULL_RUN = "sim-sharedDrive_true50ms_c1_5seeds_kernels_kinLong_kinHCV"   # e.g. "sim-sharedDrive_true50ms_c1_5seeds_kernels_kinLong_kinHCV"
_null = None
if PI1_NULL_RUN:
    _p = config.var_path(LAG_TAG, OUT_SUBDIR, f"{PI1_NULL_RUN}_neurons.csv")
    _null = pd.read_csv(cache._long(_p)) if os.path.exists(cache._long(_p)) else None
    print(f"pi1 null reference: {_p.name}" if _null is not None else f"pi1 null reference NOT found: {_p}")
runs.plot_pi1(RES, SETTING_COLS, REGIONS_PLOT, FIG_DIR, _fig, SAVE_FIGS, SHOW_FIGS,
              f"03 KERNELS  |  {RUN_NAME}", _label, null_res=_null, rate_flag_hz=100)


# %% ── 5i. simulations only: detection and where neurons land vs the truth ─
if DATA_SOURCE == "sim":
    runs.plot_detection(DET, SIM_SCENARIO, FIG_DIR, _fig("7_detection"), SAVE_FIGS, SHOW_FIGS)
    runs.plot_sim_truth(RES, SETTING_COLS, ("peak", "centroid"), SIM_PARAMS["lag_ms"],
                        FIG_DIR, _fig, SAVE_FIGS, SHOW_FIGS, _label)

# %%
