# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 02 -- Point-process GLM with discrete lag windows (sanity check between 01 and 03).
#   y = RAW spike counts in 2.5 ms bins (no smoothing)
#   x = RAW med-lat; at 2.5 ms zero-order hold is the identity (bins = camera frames)
#   Every predictor = mean of the signal's PAST lags in a window (lags >= 1 bin).
#   kin -> neuron (x->y): Poisson GLM   y ~ 1 + own-history windows [+ kinematic windows]
#   neuron -> kin (y->x): Gaussian LM   x ~ 1 + own-history windows [+ spike-count windows]
#   Own history reaches 150 ms (= the cross lags): the lesson from 01, where too short an
#   own history produced spurious neuron->kin.
#   Significance: trial shuffle within LCR (primary; score statistic for the Poisson
#   equation), BH-FDR per region x direction pooled over sessions. chi2 / F p secondary.
#   Effect sizes: bits/spike (x->y), partial R^2 (y->x). Lags: coarse windows (see note).
# Wording: "predictive precedence", never "drives"/"causes".
#
# ── HOW TO RUN THE FULL DATASET (CP vs CB) ───────────────────────────────────
#   DATA_SOURCE = "real"      # <<<
#   SESSIONS    = "all"       # <<< every session passing MIN_STRATUM (22 of 29 with cat 1)
#   N_SHUFFLE   = 1000        # <<<
#   cached per session, so an interrupted run resumes.
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
from gctools import config, io, preprocess as pp, models, nulls, stats, sim, cache, plotting, runs
from gctools import design as ds
for _m in (config, io, pp, models, nulls, stats, sim, cache, plotting, ds, runs):
    importlib.reload(_m)

# ── data ─────────────────────────────────────────────────────────────────────
DATA_SOURCE  = "sim"               # "sim" | "real"                        <<< full run: "real"
SESSIONS     = "first"             # real data only: "first" | "all" | list of (mouse, sess)
                                   #                                       <<< full run: "all"
MIN_STRATUM  = 20                  # skip a session if any LCR stratum has fewer trials
SIM_SCENARIO = "y_to_x"  # "null_independent" | "x_to_y" | "y_to_x"
SIM_PARAMS   = dict(rate_hz="match",   # target rates = the REAL neurons' window rates,
                                       # or a tuple of Hz values (e.g. (100.0, 150.0))
                    n_neurons=None,    # subset per seed (None = all)
                    n_drivers=1,       # y_to_x: neurons that LEAD the kinematics in EACH
                                       # simulated session (seed); e.g. 10 = 10 leaders in one
                                       # session. Each adds `coupling` as if alone (summed).
                    coupling=0.3,     # x_to_y: log-rate gain; y_to_x: kin SD units
                    lag_ms=50.0)       # mode of the true coupling kernel
SIM_SEEDS    = (0,3)
SESSION      = config.FIRST_PASS[0]   # template session for simulations
REACH_CAT    = config.REACH_CAT
MIN_RATE_HZ  = 1.0                 # per-NEURON inclusion (io.neuron_table), never per region

# ── design ───────────────────────────────────────────────────────────────────
BIN_MS           = 2.5             # = camera frame: kinematics enter unresampled
HIST_WINDOWS_MS  = [(0, 2.5), (2.5, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 150)]
                                   # own history, lags (lo, hi] ms -- reaches 150 ms, as the
                                   # cross lags do (both equations use the same windows)
CROSS_WINDOWS_MS = [(0, 10), (10, 20), (20, 40), (40, 80), (80, 150)]
                                   # tested cross lags, (lo, hi] ms (PLAN list cut at 150 ms)
RIDGE_ALPHA      = 0.1             # L2 on the (standardised) Poisson coefficients, identical in
                                   # full, reduced and every shuffle. Needed: for sparse neurons
                                   # the shortest own-history window is quasi-separated (a spike
                                   # is never followed by one 2.5 ms later), so its coefficient
                                   # runs to -inf without it. 0.1 changed cross-term t-values by
                                   # <= 0.1 and deviance by < 0.03 on the M3 sims.

# ── inference ────────────────────────────────────────────────────────────────
N_SHUFFLE    = 200                 #                                       <<< full run: 1000
SHUFFLE_SEED = 20
FDR_Q        = 0.05
ADAPTIVE_SHUFFLE = True            # as in 01: neurons that beat ALL N_SHUFFLE shuffles and are
                                   # not yet FDR-significant get more (same sequence continued)
                                   # until SEQ_H shuffles reach them or the cap (Besag & Clifford)
ADAPTIVE_CAP = "auto"              # "auto" = ceil(n_neurons_region / FDR_Q), or an integer
EXTRA_BATCH  = 1000
SEQ_H        = 10

N_JOBS       = 8
CACHE_MODE   = "load"              # "load" | "rerun"
SAVE_CACHE   = True
SHOW_FIGS    = True
SAVE_FIGS    = True
SAVE_VARS    = True
OUT_SUBDIR   = "02_GLM"            # figures + tables of this script live under <LAG_TAG>/02_GLM/

LAG_TAG  = config.lag_tag(config.CROSS_LAG_RANGE_MS, config.TRIAL_WINDOW_MS)
REGIONS_PLOT = ["CP", "CB"] if DATA_SOURCE == "real" else ["SIM"]
SESS_LIST = runs.session_list(DATA_SOURCE, SESSIONS, MIN_STRATUM, REACH_CAT, SESSION)
SETTING_COLS = []                  # one setting per run
CROSS_CENTERS = np.array([(a + b) / 2 for a, b in CROSS_WINDOWS_MS])
WIN_LABELS = [f"{a:g}-{b:g}" for a, b in CROSS_WINDOWS_MS]

# RUN_NAME: <data>_poissonWindows [+ _bin<ms>ms if not 2.5] [+ _ridge<alpha> if not 0.1]
#           [+ _rates<Hz>-<Hz>Hz for sims with explicit rates]; details in run_info.txt
_rates = SIM_PARAMS.get("rate_hz")
RUN_NAME = (f"{runs.data_name(DATA_SOURCE, SESS_LIST, SIM_SCENARIO, SIM_PARAMS, SIM_SEEDS)}"
            f"_poissonWindows" + (f"_bin{BIN_MS:g}ms" if BIN_MS != 2.5 else "")
            + (f"_ridge{RIDGE_ALPHA:g}" if RIDGE_ALPHA != 0.1 else "")
            + (f"_rates{'-'.join(f'{r:g}' for r in _rates)}Hz"
               if DATA_SOURCE == "sim" and not isinstance(_rates, str) else ""))
FIG_DIR  = config.fig_dir(LAG_TAG, OUT_SUBDIR, RUN_NAME) if SAVE_FIGS else None
config.warn_long(config.OUT_DIR / LAG_TAG / OUT_SUBDIR / RUN_NAME / "5_lag_windows_PREVIEW.png")
_fig = lambda base, *v: base      # single setting: plain figure names

RUN_INFO = "\n".join(
    [f"run: {RUN_NAME}", "script: 02_ppglm_discrete_lags.py (Poisson GLM on raw 2.5 ms counts, "
                         "discrete lag windows)"]
    + runs.common_info(DATA_SOURCE, SESS_LIST, REACH_CAT, MIN_STRATUM, MIN_RATE_HZ, SIM_SCENARIO,
                       SIM_PARAMS, SIM_SEEDS, SESSION, N_SHUFFLE, FDR_Q)
    + [f"bin: {BIN_MS:g} ms raw spike counts; raw kinematics (no smoothing)",
       f"own-history windows (ms): {HIST_WINDOWS_MS}",
       f"cross (tested) windows (ms): {CROSS_WINDOWS_MS}",
       f"ridge (Poisson equation): {RIDGE_ALPHA:g}",
       "shuffle statistic: Poisson score statistic (kin->neuron), ln RSS ratio (neuron->kin)",
       f"adaptive shuffles: {ADAPTIVE_SHUFFLE} (cap {ADAPTIVE_CAP}, batch {EXTRA_BATCH}, h {SEQ_H})"])
runs.write_run_info(RUN_INFO, FIG_DIR, LAG_TAG, OUT_SUBDIR, RUN_NAME, SAVE_FIGS, SAVE_VARS)


# %% ── 1. data: real sessions, or one simulated session per seed ────────────
SESSIONS_DATA, NT_ALL = runs.load_data(DATA_SOURCE, SESS_LIST, REACH_CAT, SESSION, SIM_SCENARIO,
                                       SIM_PARAMS, SIM_SEEDS, MIN_RATE_HZ)


# %% ── 2. bin (raw counts, raw kinematics) and build the lag windows ─────────
HIST_PAD_MS = max(b for _, b in HIST_WINDOWS_MS + CROSS_WINDOWS_MS)
BINNED = {tag: pp.bin_session(S, BIN_MS, hist_pad_ms=HIST_PAD_MS, kin_smooth="raw", fr_kernel=None)
          for tag, S in SESSIONS_DATA}
HIST_SETS  = ds.ms_windows_to_lags(HIST_WINDOWS_MS, BIN_MS)
CROSS_SETS = ds.ms_windows_to_lags(CROSS_WINDOWS_MS, BIN_MS)
B0 = next(iter(BINNED.values()))
_W = ds.lagged_block(B0["X"][:1], ds.window_means, HIST_SETS + CROSS_SETS)[0]   # kin windows, trial 0
_rows = B0["in_window"] & np.isfinite(_W).all(axis=1)
print(f"{len(B0['t'])} bins of {BIN_MS:g} ms; fitted rows per trial ~{_rows.sum()} "
      f"({B0['t'][_rows][0]*1000:.1f}..{B0['t'][_rows][-1]*1000:.1f} ms); "
      f"own-history windows {len(HIST_SETS)}, cross windows {len(CROSS_SETS)}")


# %% ── 3. fit every neuron (both directions, shuffle null) -- cached per session
def _fit_session(tag, S):
    B = BINNED[tag]
    perms = nulls.make_perms(B["cond"], N_SHUFFLE, seed=SHUFFLE_SEED)
    idx = np.flatnonzero(S["_keep"])
    res = Parallel(n_jobs=N_JOBS, backend="loky")(
        delayed(models.glm_neuron)(B["Y"][i], B["X"], B["in_window"], HIST_SETS, CROSS_SETS,
                                   CROSS_CENTERS, perms, RIDGE_ALPHA)
        for i in idx)
    for r in res:
        for d in ("xy", "yx"):
            r[d]["null"] = r[d]["null"].astype(np.float32)
    return dict(idx=idx, res=res)

FITS, FIT_PARAMS = {}, {}
t0 = time.time()
for si, (tag, S) in enumerate(SESSIONS_DATA):
    stamp = io.mi_stamp() if S["source"] == "real" else dict(SIM_PARAMS, scenario=SIM_SCENARIO, tag=tag)
    params = dict(stamp=stamp, session=S["key"], reach_cat=REACH_CAT,
                  uids=[n["uid"] for n in S["neurons"]], keep=S["_keep"], bin_ms=BIN_MS,
                  hist=HIST_WINDOWS_MS, cross=CROSS_WINDOWS_MS, ridge=RIDGE_ALPHA,
                  n_shuffle=N_SHUFFLE, shuffle_seed=SHUFFLE_SEED, window=config.TRIAL_WINDOW_MS,
                  shuffle_within=config.SHUFFLE_WITHIN, kin_var=config.KIN_VAR)
    FIT_PARAMS[tag] = params
    FITS[tag] = cache.cached(f"02glm_{tag}", params, lambda: _fit_session(tag, S),
                             mode=CACHE_MODE, save=SAVE_CACHE)
    print(f"  session {si+1}/{len(SESSIONS_DATA)} ({tag}) done, {time.time()-t0:.0f}s elapsed", flush=True)
print(f"fits ready in {time.time()-t0:.0f}s")


# %% ── 4. per-neuron table, BH-FDR, direction class, summaries ────────────────
rows, PROFILES = [], {}
for tag, S in SESSIONS_DATA:
    F = FITS[tag]
    for i, r in zip(F["idx"], F["res"]):
        row = runs.base_row(tag, S, i)
        xy, yx = r["xy"], r["yx"]
        row.update(n_spikes=r["n_spikes"], rows_per_trial=r["rows_per_trial"],
                   pS_xy=xy["p_shuf"], z_xy=xy["z"], lr_xy=xy["lr"], pChi2_xy=xy["p_chi2"],
                   score_xy=xy["score"],
                   bits_per_spike_xy=xy["bits_per_spike"], dev_expl_red_xy=xy["dev_expl_red"],
                   dev_expl_full_xy=xy["dev_expl_full"], dispersion=xy["dispersion"],
                   converged=xy["converged"],
                   pS_yx=yx["p_shuf"], z_yx=yx["z"], gc_yx=yx["gc"], pF_yx=yx["p_F"],
                   r2_red_yx=yx["r2_red"], r2_full_yx=yx["r2_full"], partial_r2_yx=yx["partial_r2"],
                   bpsx_xy=models.excess_bits_per_spike(xy["lr"], xy["null"], r["n_spikes"]),
                   pr2x_yx=models.excess_partial_r2(yx["gc"], yx["null"]))
        for d in ("xy", "yx"):
            L = r[d]["lag"]
            row.update({f"lag_best_{d}": L["best_ms"], f"lag_sign_{d}": L["best_sign"],
                        f"lag_centroid_{d}": L["centroid_ms"]})
            row.update({f"t_{d}_{w}": t for w, t in zip(WIN_LABELS, L["t"])})
            PROFILES[(tag, S["neurons"][i]["uid"], d)] = L
        rows.append(row)
RES = pd.DataFrame(rows)

# ---- adaptive shuffles for neurons at the p floor (before FDR) ---------------
if ADAPTIVE_SHUFFLE:
    RES = runs.adaptive_extend(
        RES, SETTING_COLS, {"xy": "score_xy", "yx": "gc_yx"}, N_SHUFFLE, FDR_Q, ADAPTIVE_CAP,
        SEQ_H, EXTRA_BATCH, SHUFFLE_SEED,
        cond_of=lambda tag, v: BINNED[tag]["cond"],
        null_of=lambda row, d: models.GlmNull(BINNED[row.tag]["Y"][row.neuron_idx], BINNED[row.tag]["X"],
                                              BINNED[row.tag]["in_window"], HIST_SETS, CROSS_SETS,
                                              RIDGE_ALPHA, d),
        cache_params=lambda tag, v: FIT_PARAMS[tag], n_jobs=N_JOBS, cache_name="02glm",
        cache_mode=CACHE_MODE, save_cache=SAVE_CACHE)

RES = runs.fdr_and_classes(RES, SETTING_COLS, FDR_Q)
runs.print_direction_summary(RES, SETTING_COLS, FDR_Q, len(SESSIONS_DATA))
FLOOR = runs.print_floor_summary(RES, SETTING_COLS, N_SHUFFLE, FDR_Q)
EFFECT_COLS = {"xy": ("bpsx_xy", "excess bits/spike from kinematics"),
               "yx": ("pr2x_yx", "excess partial R² of med-lat (%)")}
runs.print_effect_summary(RES, SETTING_COLS, REGIONS_PLOT, EFFECT_COLS)

print(f"\n=== {models.GLM_LAG_NOTE} ===")
runs.print_lag_summary(RES, SETTING_COLS, REGIONS_PLOT, ("best", "centroid"), "LAG WINDOW")

print("\n=== EFFECT SIZE (median [IQR]; significant neurons of that direction) ===")
for reg in REGIONS_PLOT:
    g = RES[RES.region == reg]
    for d, col, unit in (("xy", "bits_per_spike_xy", "bits/spike"), ("yx", "partial_r2_yx", "partial R2")):
        v = g[g[f"sig_{d}"]][col]
        if len(v):
            print(f"  {reg} {plotting.DIR_LABEL[d]}: {unit} median {v.median():.4f} "
                  f"[{v.quantile(.25):.4f}-{v.quantile(.75):.4f}] n={len(v)}")
disp = RES.groupby("region")["dispersion"].describe(percentiles=[.1, .5, .9])[["10%", "50%", "90%"]]
print(f"\nPoisson dispersion (Pearson chi2/df, reduced model; 1 = Poisson, >1 overdispersed):\n"
      f"{disp.round(2).to_string()}\nnot converged: {int((~RES.converged).sum())}")

if DATA_SOURCE == "sim":
    DET = runs.sim_truth_check(RES, SETTING_COLS, lag_keys=("best", "centroid"),
                               true_lag=SIM_PARAMS["lag_ms"])
    print("  chi2 LR false-positive rate (x->y, expected False): "
          f"{np.mean(RES.loc[~RES.exp_xy, 'pChi2_xy'] < 0.05):.3f}   "
          f"F-test (y->x, expected False): {np.mean(RES.loc[~RES.exp_yx, 'pF_yx'] < 0.05):.3f}")

if SAVE_VARS:
    runs.save_tables(dict(NT=NT_ALL, profiles=PROFILES, lag_note=models.GLM_LAG_NOTE, params=dict(
        DATA_SOURCE=DATA_SOURCE, SESS_LIST=SESS_LIST, SIM_SCENARIO=SIM_SCENARIO,
        SIM_PARAMS=SIM_PARAMS, SIM_SEEDS=SIM_SEEDS, BIN_MS=BIN_MS, HIST_WINDOWS_MS=HIST_WINDOWS_MS,
        CROSS_WINDOWS_MS=CROSS_WINDOWS_MS, RIDGE_ALPHA=RIDGE_ALPHA, N_SHUFFLE=N_SHUFFLE,
        FDR_Q=FDR_Q, MIN_RATE_HZ=MIN_RATE_HZ, MIN_STRATUM=MIN_STRATUM, LAG_TAG=LAG_TAG)),
        RES, LAG_TAG, OUT_SUBDIR, RUN_NAME)


# %% ── 5a/5b. direction classes (CP vs CB) and direction map per region ─────
runs.plot_direction_figures(RES, SETTING_COLS, REGIONS_PLOT, N_SHUFFLE, FDR_Q, FIG_DIR, _fig,
                            SAVE_FIGS, SHOW_FIGS, f"02 GLM  |  {RUN_NAME}")


# %% ── 5c. effect sizes per region: FDR-significant (top row) and all neurons (bottom)
# Columns = direction; compare regions WITHIN a panel only (never across columns).
fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), squeeze=False)
for c_, (d, col, lab) in enumerate((("xy", "bits_per_spike_xy", "gain from kinematics (bits/spike)"),
                                    ("yx", "partial_r2_yx", "partial R² from spikes"))):
    for r_, (sel_name, only_sig) in enumerate((("FDR-significant", True), ("all neurons", False))):
        ax = axs[r_, c_]
        data, labels, cols = [], [], []
        for reg in REGIONS_PLOT:
            g = RES[RES.region == reg]
            if only_sig:
                g = g[g[f"sig_{d}"]]
            v = g[col].dropna().to_numpy()
            data.append(np.clip(v, 1e-6, None)); labels.append(f"{reg}\n(n={len(v)})")
            cols.append(config.CLR.get(reg, (0.4, 0.4, 0.4)))
        if all(len(v) == 0 for v in data):
            ax.text(0.5, 0.5, "no neurons", ha="center", va="center", transform=ax.transAxes, color="0.4")
        else:
            bp = ax.boxplot([v if len(v) else [np.nan] for v in data], patch_artist=True,
                            showfliers=False, widths=0.6)
            for patch, c in zip(bp["boxes"], cols):
                patch.set_facecolor(c); patch.set_alpha(0.85 if only_sig else 0.45)
            for j, v in enumerate(data):
                ax.scatter(np.full(v.size, j + 1) + np.random.default_rng(j).uniform(-0.15, 0.15, v.size),
                           v, s=4, color="k", alpha=0.3)
            ax.set_yscale("log")
        ax.set_xticks(range(1, len(labels) + 1), labels, fontsize=8)
        ax.set(xlim=(0.5, len(labels) + 0.5), ylabel=lab, title=f"{plotting.DIR_LABEL[d]} — {sel_name}")
fig.suptitle("effect size per neuron: top = FDR-significant, bottom = all neurons "
             "(compare regions WITHIN a panel only)", fontsize=9)
plotting.finish(fig, FIG_DIR, "3_effect_sizes", SAVE_FIGS, SHOW_FIGS)

# same content as histograms (pick whichever reads better): log-spaced bins, each region
# normalised to its own n (fraction of neurons), dashed line = median
fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), squeeze=False)
for c_, (d, col, lab) in enumerate((("xy", "bits_per_spike_xy", "gain from kinematics (bits/spike)"),
                                    ("yx", "partial_r2_yx", "partial R² from spikes"))):
    vals_all = np.clip(RES[col].dropna().to_numpy(), 1e-6, None)
    edges = (np.logspace(np.log10(vals_all.min()), np.log10(vals_all.max()), 31)
             if vals_all.size else np.logspace(-6, 0, 31))                    # shared per column
    for r_, (sel_name, only_sig) in enumerate((("FDR-significant", True), ("all neurons", False))):
        ax = axs[r_, c_]
        any_data = False
        for reg in REGIONS_PLOT:
            g = RES[RES.region == reg]
            if only_sig:
                g = g[g[f"sig_{d}"]]
            v = np.clip(g[col].dropna().to_numpy(), 1e-6, None)
            if not v.size:
                continue
            any_data = True
            c = config.CLR.get(reg, (0.4, 0.4, 0.4))
            ax.hist(v, bins=edges, weights=np.full(v.size, 1 / v.size), color=c, alpha=0.5,
                    label=f"{reg} (n={v.size}, median {np.median(v):.2g})")
            ax.hist(v, bins=edges, weights=np.full(v.size, 1 / v.size), color=c, histtype="step", lw=1.2)
            ax.axvline(np.median(v), color=c, ls="--", lw=1.4)
        if not any_data:
            ax.text(0.5, 0.5, "no neurons", ha="center", va="center", transform=ax.transAxes, color="0.4")
        else:
            ax.legend(fontsize=7, frameon=False)
        ax.set_xscale("log")
        ax.set(xlim=(edges[0], edges[-1]), xlabel=lab, ylabel="fraction of neurons",
               title=f"{plotting.DIR_LABEL[d]} — {sel_name}")
fig.suptitle("effect size per neuron (histograms): top = FDR-significant, bottom = all neurons "
             "(compare regions WITHIN a panel only)", fontsize=9)
plotting.finish(fig, FIG_DIR, "3b_effect_sizes_hist", SAVE_FIGS, SHOW_FIGS)

# SHUFFLE-CORRECTED version: effect minus the mean of the neuron's own shuffles
# (0 = chance level; can be negative). Symmetric-log x axis: linear within +-s around 0,
# log beyond; s = median |excess| over all neurons (the size of the null scatter), so null
# neurons spread around 0 and real effects separate on a log-like scale.
fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), squeeze=False)
for c_, (d, col, lab) in enumerate((("xy", "bpsx_xy", "excess gain from kinematics (bits/spike)"),
                                    ("yx", "pr2x_yx", "excess partial R² from spikes (%)"))):
    va = RES[col].dropna().to_numpy()
    s = float(np.median(np.abs(va))) if va.size else 1.0
    s = s if s > 0 else 1e-6
    lo_, hi_ = (va.min(), va.max()) if va.size else (-1, 1)
    edges = s * np.sinh(np.linspace(np.arcsinh(min(lo_, -s) / s), np.arcsinh(max(hi_, s) / s), 41))
    for r_, (sel_name, only_sig) in enumerate((("FDR-significant", True), ("all neurons", False))):
        ax = axs[r_, c_]
        any_data = False
        for reg in REGIONS_PLOT:
            g = RES[RES.region == reg]
            if only_sig:
                g = g[g[f"sig_{d}"]]
            v = g[col].dropna().to_numpy()
            if not v.size:
                continue
            any_data = True
            c = config.CLR.get(reg, (0.4, 0.4, 0.4))
            w = np.full(v.size, 1 / v.size)
            ax.hist(v, bins=edges, weights=w, color=c, alpha=0.5,
                    label=f"{reg} (n={v.size}, median {np.median(v):.2g})")
            ax.hist(v, bins=edges, weights=w, color=c, histtype="step", lw=1.2)
            ax.axvline(np.median(v), color=c, ls="--", lw=1.4)
        ax.axvline(0, color="k", lw=0.9, ls=":")
        if not any_data:
            ax.text(0.5, 0.5, "no neurons", ha="center", va="center", transform=ax.transAxes, color="0.4")
        else:
            ax.legend(fontsize=7, frameon=False)
        ax.set_xscale("symlog", linthresh=s)
        k0 = int(np.floor(np.log10(s))) + 1            # decades beyond the linear zone only
        dec = 10.0 ** np.arange(k0, int(np.ceil(np.log10(max(abs(edges[0]), abs(edges[-1]))))) + 1)
        ticks = [0.0] + [t for t in dec if t <= edges[-1]] + [-t for t in dec if -t >= edges[0]]
        ax.set_xticks(sorted(ticks))
        ax.set(xlim=(edges[0], edges[-1]), xlabel=lab, ylabel="fraction of neurons",
               title=f"{plotting.DIR_LABEL[d]} — {sel_name}")
fig.suptitle("SHUFFLE-CORRECTED effect size (effect − mean of own shuffles; dotted 0 = chance): "
             "top = FDR-significant, bottom = all neurons  (symlog axis)", fontsize=9)
plotting.finish(fig, FIG_DIR, "3c_effect_sizes_excess_hist", SAVE_FIGS, SHOW_FIGS)

# observed vs the trial-shuffle null (= chance level), same units as 3 / 3b. Rows = region.
# x->y shuffles are score statistics; shown as bits/spike via score / (2 n_spikes ln 2)
# (score ~ LR for large samples). y->x shuffles: partial R^2 = 1 - exp(-ln RSS ratio).
PAIRS = {}
for tag, S in SESSIONS_DATA:
    F = FITS[tag]
    for i, r in zip(F["idx"], F["res"]):
        reg, ns = S["neurons"][i]["region"], max(r["n_spikes"], 1.0)
        PAIRS.setdefault((reg, "xy"), []).append(
            (r["xy"]["bits_per_spike"], np.asarray(r["xy"]["null"], float) / (2 * ns * np.log(2))))
        PAIRS.setdefault((reg, "yx"), []).append(
            (r["yx"]["partial_r2"], 1 - np.exp(-np.asarray(r["yx"]["null"], float))))
runs.plot_null_vs_observed(PAIRS, REGIONS_PLOT,
                           {"xy": "gain from kinematics (bits/spike)", "yx": "partial R² from spikes"},
                           FIG_DIR, "3d_null_vs_observed", SAVE_FIGS, SHOW_FIGS, f"02 GLM  |  {RUN_NAME}")


# %% ── 5d. Poisson dispersion of the reduced model, per region ───────────────
fig, ax = plt.subplots(figsize=(6, 3.6))
bins = np.linspace(0.5, max(2.0, RES.dispersion.quantile(0.99)), 40)
for reg in REGIONS_PLOT:
    v = RES[RES.region == reg].dispersion
    ax.hist(v, bins=bins, alpha=0.6, color=config.CLR.get(reg, (0.4, 0.4, 0.4)),
            label=f"{reg} median {v.median():.2f}")
ax.axvline(1, color="k", lw=0.8, ls="--")
ax.set(xlabel="Pearson χ²/df (reduced Poisson model)", ylabel="# neurons",
       title="dispersion: 1 = Poisson; >1 over-, <1 under-dispersed")
ax.legend(fontsize=7, frameon=False)
plotting.finish(fig, FIG_DIR, "4_dispersion", SAVE_FIGS, SHOW_FIGS)


# %% ── 5e. lag windows (coarse; conditional on own history) ──────────────────
fig, axs = plt.subplots(2, 2, figsize=(11, 7))
for c_, d in enumerate(("xy", "yx")):
    ax = axs[0, c_]                          # which window carries the largest |t|
    wpos = np.arange(len(WIN_LABELS))
    for k, reg in enumerate(REGIONS_PLOT):
        g = RES[(RES.region == reg) & RES[f"sig_{d}"]]
        cnt = np.array([np.sum(np.isclose(g[f"lag_best_{d}"], c)) for c in CROSS_CENTERS])
        ax.bar(wpos + (k - (len(REGIONS_PLOT) - 1) / 2) * 0.35, cnt / max(len(g), 1), width=0.35,
               color=config.CLR.get(reg, (0.4, 0.4, 0.4)), label=f"{reg} (n={len(g)})")
    ax.set_xticks(wpos, [f"{w} ms" for w in WIN_LABELS], fontsize=7)
    ax.set(ylabel="fraction of significant neurons", title=f"strongest window, {plotting.DIR_LABEL[d]}")
    ax.legend(fontsize=7, frameon=False)
    ax = axs[1, c_]                          # mean |t|-normalised profile
    for reg in REGIONS_PLOT:
        g = RES[(RES.region == reg) & RES[f"sig_{d}"]]
        if not len(g):
            continue
        T = g[[f"t_{d}_{w}" for w in WIN_LABELS]].to_numpy(float)
        T = T / np.maximum(np.abs(T).max(axis=1, keepdims=True), 1e-12)
        m, se = T.mean(0), T.std(0, ddof=1) / np.sqrt(len(T)) if len(T) > 1 else 0 * T.mean(0)
        c = config.CLR.get(reg, (0.4, 0.4, 0.4))
        ax.errorbar(CROSS_CENTERS, m, yerr=se, color=c, marker="o", lw=1.5, capsize=3,
                    label=f"{reg} (n={len(g)})")
    ax.axhline(0, color="k", lw=0.4)
    ax.set(xlabel="lag window centre (ms)", ylabel="mean t / max|t|  (± SEM)",
           title=f"mean lag profile, {plotting.DIR_LABEL[d]}")
    ax.legend(fontsize=7, frameon=False)
fig.suptitle("LAG WINDOWS — coarse, conditional on own history (grain of salt); FDR-significant "
             "neurons only\nsims, true lag 50 ms: kin→neuron strongest 40-80 ms (centroid 65); "
             "neuron→kin strongest 20-40 ms (centroid 37; 27 with 10 leaders/session) — "
             "UNDERESTIMATE", fontsize=9)
plotting.finish(fig, FIG_DIR, "5_lag_windows_PREVIEW", SAVE_FIGS, SHOW_FIGS)


# %% ── 5f. simulations only: detection vs ground truth ───────────────────────
if DATA_SOURCE == "sim":
    runs.plot_detection(DET, SIM_SCENARIO, FIG_DIR, "6_detection", SAVE_FIGS, SHOW_FIGS)
    runs.plot_sim_truth(RES, SETTING_COLS, ("best", "centroid"), SIM_PARAMS["lag_ms"],
                        FIG_DIR, _fig, SAVE_FIGS, SHOW_FIGS)

# %% ── 5g. effect STRENGTH and FRACTION of neurons with an effect ────────────
# As in 01 (5h). Effect: excess bits/spike (kin->neuron; mis-scaled above 100 Hz,
# marked x) and excess partial R^2 of med-lat (neuron->kin). Compare within a row only.
runs.plot_effect_figures(RES, SETTING_COLS, REGIONS_PLOT, EFFECT_COLS, FIG_DIR, _fig,
                         SAVE_FIGS, SHOW_FIGS, f"02 GLM  |  {RUN_NAME}", rate_flag_hz=100)

# %%
