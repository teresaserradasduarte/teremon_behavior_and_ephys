# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 01 -- Bivariate Gaussian VAR on continuous signals (baseline, not final inference:
# a Gaussian model on smoothed spike counts is misspecified; the GLMs of 02/03 are).
#   y = spike train -> causal gamma rate (FR_KERNEL) -> BIN_MS bins
#   x = med-lat, zero-order hold into the same bins; "raw" or "matched" (same kernel)
#   Per neuron and direction: GC = ln(RSS_red / RSS_full), trial-shuffle p within LCR
#   (primary), F-test p (secondary), null-z; BH-FDR per region POOLED over sessions.
#   Direction class per neuron: kin→neuron / neuron→kin / both / none.
#   Lag: PREVIEW estimate from a smooth causal cross-kernel (see models.LAG_NOTE):
#   use the CENTROID (sims, true 50 ms: 52 / 43 ms); the PEAK is unreliable.
# Wording: "predictive precedence", never "drives"/"causes".
#
# ── HOW TO RUN THE FULL DATASET (CP vs CB) ───────────────────────────────────
#   DATA_SOURCE = "real"      # <<<
#   SESSIONS    = "all"       # <<< every session passing MIN_STRATUM (22 of 29 with cat 1)
#   N_SHUFFLE   = 1000        # <<< FDR over hundreds of neurons needs fine p-values
#   ~1 h on 8 cores for one (smoothing x order) setting; cached per session, so an
#   interrupted run resumes where it stopped. Outputs: see the paths printed below.
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
DATA_SOURCE  = "real"               # "sim" | "real"                        <<< full run: "real"
SESSIONS     = "all"             # real data only: "first" (FIRST_PASS[0]) | "all" |
                                   # list of (mouse, sess), e.g. [("5_FerreroRocher", "R1")]
                                   #                                       <<< full run: "all"
MIN_STRATUM  = 20                  # skip a session if any LCR stratum has fewer trials
                                   # (the within-side shuffle needs partners)
SIM_SCENARIO = "null_independent"  # "null_independent" | "x_to_y" | "y_to_x"
SIM_PARAMS   = dict(rate_hz="match",   # target rates = the REAL neurons' window rates
                    n_neurons=None,    # subset per seed (None = all)
                    n_drivers=1,       # y_to_x: neurons that LEAD the kinematics in EACH
                                       # simulated session (seed); e.g. 10 = 10 leaders in one
                                       # session. Each adds `coupling` as if alone (summed).
                    coupling=1.0,     # x_to_y: log-rate gain; y_to_x: kin SD units
                    lag_ms=50.0)       # mode of the true coupling kernel
SIM_SEEDS    = (0,)                # each seed = one simulated session
SESSION      = config.FIRST_PASS[0]   # template session for simulations
REACH_CAT    = config.REACH_CAT
MIN_RATE_HZ  = 1.0                 # per-NEURON inclusion (io.neuron_table), never per region

# ── signals ──────────────────────────────────────────────────────────────────
BIN_MS       = 10.0
FR_KERNEL_SHAPE    = 2.0           # causal gamma rate kernel; mean delay = shape*scale
FR_KERNEL_SCALE_MS = 10.0
KIN_SMOOTHS  = ("matched",)        # "matched": x smoothed with the same kernel as y, so
                                   # RELATIVE timing (and the lag readout) is preserved.
                                   # "raw": only y carries the ~20 ms kernel delay -> x->y
                                   # lags read longer, y->x shorter. ("raw", "matched")
                                   # runs both and adds the comparison figure.
REMOVE_ENSEMBLE_MEAN = False       # subtract leave-one-out LCR mean from y and x

# ── model / inference ────────────────────────────────────────────────────────
P_MAX        = int(config.CROSS_LAG_RANGE_MS[1] / BIN_MS)   # 15 lags = 150 ms
ORDERS       = ("pmax",)           # "pmax" = p = P_MAX always (default, validated)
                                   # "bic"  = p by BIC. M2 sims: BIC picks p~8, too short
                                   # for the kinematics' own history -> spurious y->x (42%
                                   # of neurons under x_to_y). ("bic", "pmax") runs both
                                   # and flags ORDER-SENSITIVE neurons.
N_SHUFFLE    = 1000                 #                                       <<< full run: 1000
SHUFFLE_SEED = 0
FDR_Q        = 0.05
ADAPTIVE_SHUFFLE = True            # neurons that beat ALL N_SHUFFLE shuffles get more (same
                                   # shuffle sequence continued) until SEQ_H shuffles reach
                                   # them or the cap; p stays purely permutation-based
                                   # (Besag & Clifford). Needed so BH-FDR is not limited by
                                   # the p floor 1/(N_SHUFFLE+1) when few neurons carry an effect.
ADAPTIVE_CAP = "auto"              # "auto" = ceil(n_neurons_region / FDR_Q): beyond it more
                                   # shuffles cannot change the FDR decision. Or an integer.
EXTRA_BATCH  = 1000                # shuffles added per step
SEQ_H        = 10                  # stop once this many shuffles reach the neuron (p = h / L)

# ── lag readout (PREVIEW, biased -- see models.LAG_NOTE) ────────────────────
LAG_READOUT  = True
N_BASES      = 6                   # raised-cosine bases over CROSS_LAG_RANGE_MS
BASIS_SPACING = "log"              # "log" | "uniform" spacing of the bases.
                                   # Measured on sims (true lag 50 ms, matched kin, 6 bases),
                                   # centroid read: log 52 (x->y) / 43 (y->x) ms;
                                   # uniform 65 / 47 ms. "log" = finer at short lags.
N_BOOT       = 100                 # trial bootstrap for lag CIs

N_JOBS       = 8
CACHE_MODE   = "load"              # "load" | "rerun"
SAVE_CACHE   = True
SHOW_FIGS    = True
SAVE_FIGS    = True
SAVE_VARS    = True
OUT_SUBDIR   = "01_VAR"            # figures + tables of this script live under <LAG_TAG>/01_VAR/
N_EXAMPLE_KERNELS = 4              # example kernels per region x direction

LAG_TAG  = config.lag_tag(config.CROSS_LAG_RANGE_MS, config.TRIAL_WINDOW_MS)
FR_KERNEL = (FR_KERNEL_SHAPE, FR_KERNEL_SCALE_MS)
REGIONS_PLOT = ["CP", "CB"] if DATA_SOURCE == "real" else ["SIM"]

SESS_LIST = runs.session_list(DATA_SOURCE, SESSIONS, MIN_STRATUM, REACH_CAT, SESSION)
SETTINGS = [(ks, o) for ks in KIN_SMOOTHS for o in ORDERS]
SETTING_COLS = ["kin_smooth", "order"]
_LAGS_NAME = {"pmax": f"{P_MAX}lags", "bic": "BIClags"}
def _label(v):                    # human label of one setting (kin_smooth, order)
    return f"{v[0]} kin, {_LAGS_NAME[v[1]]}"

# ── RUN_NAME: readable, but short enough for Windows/Dropbox path limits ─────
#   real_<session or N sessions>  |  sim-<scenario>_true<lag>ms_c<coupling>_<N>seeds
#   _kin<Matched|Raw|RawMatched>   kinematic smoothing
#   _lags<15|BIC|BIC+15>           VAR order: fixed P_MAX lags, BIC, or both
#   only when not the default:  _bin<ms>ms (default 10), _rateG<shape>x<scale>ms (rate
#   kernel, default gamma 2x10 ms), _bases<Log|Uniform><N> (default Log6), _meanRemoved
# Everything (incl. defaults) is spelled out in run_info.txt next to the figures.
_kin  = "kin" + "".join(k.capitalize() for k in KIN_SMOOTHS)
_lags = "lags" + "+".join({"pmax": f"{P_MAX}", "bic": "BIC"}[o] for o in ORDERS)
_extra = ((f"_rateG{FR_KERNEL_SHAPE:g}x{FR_KERNEL_SCALE_MS:g}ms" if FR_KERNEL != (2.0, 10.0) else "")
          + (f"_bases{BASIS_SPACING.capitalize()}{N_BASES}"
             if (BASIS_SPACING, N_BASES) != ("log", 6) else "")
          + ("_meanRemoved" if REMOVE_ENSEMBLE_MEAN else ""))
RUN_NAME = (f"{runs.data_name(DATA_SOURCE, SESS_LIST, SIM_SCENARIO, SIM_PARAMS, SIM_SEEDS)}"
            f"_{_kin}_{_lags}" + (f"_bin{BIN_MS:g}ms" if BIN_MS != 10.0 else "") + _extra)
FIG_DIR  = config.fig_dir(LAG_TAG, OUT_SUBDIR, RUN_NAME) if SAVE_FIGS else None
config.warn_long(config.OUT_DIR / LAG_TAG / OUT_SUBDIR / RUN_NAME / "7_mean_kernel_matched_BIClags.png")

def _fig(base, *v):               # figure file name; setting suffix only if >1 setting ran
    if len(v) == 1 and isinstance(v[0], tuple):
        v = v[0]
    if len(SETTINGS) == 1 or not v:
        return base
    ks, o = (list(v) + [None, None])[:2]
    return "_".join([base] + ([ks] if ks else []) + ([_LAGS_NAME[o]] if o else []))

RUN_INFO = "\n".join(
    [f"run: {RUN_NAME}", "script: 01_var_continuous.py (bivariate Gaussian VAR; didactic baseline)"]
    + runs.common_info(DATA_SOURCE, SESS_LIST, REACH_CAT, MIN_STRATUM, MIN_RATE_HZ, SIM_SCENARIO,
                       SIM_PARAMS, SIM_SEEDS, SESSION, N_SHUFFLE, FDR_Q)
    + [f"bin: {BIN_MS:g} ms;  neural rate kernel: causal gamma shape {FR_KERNEL_SHAPE:g} x scale "
       f"{FR_KERNEL_SCALE_MS:g} ms (mean delay {FR_KERNEL_SHAPE*FR_KERNEL_SCALE_MS:g} ms)",
       f"kinematic smoothing: {list(KIN_SMOOTHS)};  VAR order: {list(ORDERS)} "
       f"(pmax = {P_MAX} lags = {P_MAX*BIN_MS:g} ms)",
       f"lag readout: {'on' if LAG_READOUT else 'off'}, {N_BASES} {BASIS_SPACING} raised-cosine bases, "
       f"{N_BOOT} trial bootstraps (PREVIEW, biased; use the centroid)",
       f"ensemble mean removed: {REMOVE_ENSEMBLE_MEAN}",
       f"adaptive shuffles: {ADAPTIVE_SHUFFLE} (cap {ADAPTIVE_CAP}, batch {EXTRA_BATCH}, h {SEQ_H})"])
runs.write_run_info(RUN_INFO, FIG_DIR, LAG_TAG, OUT_SUBDIR, RUN_NAME, SAVE_FIGS, SAVE_VARS)

# %% ── 1. data: real sessions, or one simulated session per seed ────────────
SESSIONS_DATA, NT_ALL = runs.load_data(DATA_SOURCE, SESS_LIST, REACH_CAT, SESSION, SIM_SCENARIO,
                                       SIM_PARAMS, SIM_SEEDS, MIN_RATE_HZ)

# %% ── 2. bin once per session x kinematic smoothing ─────────────────────────
HIST_PAD_MS = P_MAX * BIN_MS
BINNED = {(tag, ks): pp.bin_session(S, BIN_MS, hist_pad_ms=HIST_PAD_MS, kin_smooth=ks,
                                    fr_kernel=FR_KERNEL, remove_mean=REMOVE_ENSEMBLE_MEAN)
          for tag, S in SESSIONS_DATA for ks in KIN_SMOOTHS}
B0 = next(iter(BINNED.values()))
_, _, _rows = models.var_lags(B0["Y"][0], B0["X"], B0["in_window"], P_MAX)
BASIS, BASIS_LAGS, BASIS_C = ds.raised_cosine_basis(N_BASES, config.CROSS_LAG_RANGE_MS, BIN_MS,
                                                    BASIS_SPACING)
print(f"rate kernel gamma({FR_KERNEL_SHAPE:g}, {FR_KERNEL_SCALE_MS:g} ms): mean delay "
      f"{B0['fr_kernel_mean_delay_ms']:.0f} ms | fitted rows per trial {_rows.sum()} "
      f"({B0['t'][_rows][0]*1000:.0f}..{B0['t'][_rows][-1]*1000:.0f} ms)")
print(f"lag readout: {N_BASES} {BASIS_SPACING} raised-cosine bases, centres {np.round(BASIS_C).astype(int).tolist()} ms")


# %% ── 3. fit every neuron (both directions, shuffle null, lag readout) -- cached
def _fit_session(tag, S, ks, order):
    B = BINNED[(tag, ks)]
    perms = nulls.make_perms(B["cond"], N_SHUFFLE, seed=SHUFFLE_SEED)
    idx = np.flatnonzero(S["_keep"])
    p_fixed = P_MAX if order == "pmax" else None
    lr = (dict(B=BASIS, lags=BASIS_LAGS, bin_ms=BIN_MS, n_boot=N_BOOT, seed=SHUFFLE_SEED)
          if LAG_READOUT else None)
    res = Parallel(n_jobs=N_JOBS, backend="loky")(
        delayed(models.var_neuron)(B["Y"][i], B["X"], B["in_window"], P_MAX, perms, p_fixed, lr)
        for i in idx)
    for r in res:                                   # keep nulls compact
        for d in ("xy", "yx"):
            r[d]["null"] = r[d]["null"].astype(np.float32)
    return dict(idx=idx, res=res)

FITS, FIT_PARAMS = {}, {}
t0 = time.time()
for si, (tag, S) in enumerate(SESSIONS_DATA):
    uids = [n["uid"] for n in S["neurons"]]
    stamp = io.mi_stamp() if S["source"] == "real" else dict(SIM_PARAMS, scenario=SIM_SCENARIO, tag=tag)
    for ks in KIN_SMOOTHS:
        for order in ORDERS:
            params = dict(stamp=stamp, session=S["key"], reach_cat=REACH_CAT, uids=uids,
                          keep=S["_keep"], bin_ms=BIN_MS, fr_kernel=FR_KERNEL, kin_smooth=ks,
                          order=order, erm=REMOVE_ENSEMBLE_MEAN, p_max=P_MAX,
                          n_shuffle=N_SHUFFLE, shuffle_seed=SHUFFLE_SEED,
                          window=config.TRIAL_WINDOW_MS, shuffle_within=config.SHUFFLE_WITHIN,
                          kin_var=config.KIN_VAR,
                          lag=(N_BASES, BASIS_SPACING, N_BOOT, config.CROSS_LAG_RANGE_MS)
                              if LAG_READOUT else None)
            FIT_PARAMS[(tag, ks, order)] = params
            FITS[(tag, ks, order)] = cache.cached(f"01var_{tag}_{ks}_{order}", params,
                                                  lambda: _fit_session(tag, S, ks, order),
                                                  mode=CACHE_MODE, save=SAVE_CACHE)
    print(f"  session {si+1}/{len(SESSIONS_DATA)} ({tag}) done, {time.time()-t0:.0f}s elapsed", flush=True)
print(f"fits ready in {time.time()-t0:.0f}s")


# %% ── 4. per-neuron table, BH-FDR (region pooled over sessions), direction class
def _table(tag, S, ks, order):
    F = FITS[(tag, ks, order)]
    rows = []
    for i, r in zip(F["idx"], F["res"]):
        row = dict(runs.base_row(tag, S, i), kin_smooth=ks, order=order, p=r["p"])
        for d in ("xy", "yx"):
            row.update({f"pS_{d}": r[d]["p_shuf"], f"z_{d}": r[d]["z"], f"gc_{d}": r[d]["gc"],
                        f"pF_{d}": r[d]["p_F"],
                        f"pr2x_{d}": models.excess_partial_r2(r[d]["gc"], r[d]["null"])})
            if "lag" in r[d]:
                L = r[d]["lag"]
                row.update({f"lag_peak_{d}": L["peak_ms"], f"lag_peak_lo_{d}": L["peak_ci"][0],
                            f"lag_peak_hi_{d}": L["peak_ci"][1], f"lag_sign_{d}": L["peak_sign"],
                            f"lag_centroid_{d}": L["centroid_ms"],
                            f"lag_centroid_lo_{d}": L["centroid_ci"][0],
                            f"lag_centroid_hi_{d}": L["centroid_ci"][1]})
        rows.append(row)
    return pd.DataFrame(rows)

RES = pd.concat([_table(tag, S, ks, o) for tag, S in SESSIONS_DATA for ks in KIN_SMOOTHS
                 for o in ORDERS], ignore_index=True)

# ---- adaptive shuffles for neurons at the p floor (before FDR) ---------------
if ADAPTIVE_SHUFFLE:
    RES = runs.adaptive_extend(
        RES, SETTING_COLS, {"xy": "gc_xy", "yx": "gc_yx"}, N_SHUFFLE, FDR_Q, ADAPTIVE_CAP,
        SEQ_H, EXTRA_BATCH, SHUFFLE_SEED,
        cond_of=lambda tag, v: BINNED[(tag, v[0])]["cond"],
        null_of=lambda row, d: models.VarNull(BINNED[(row.tag, row.kin_smooth)]["Y"][row.neuron_idx],
                                              BINNED[(row.tag, row.kin_smooth)]["X"],
                                              BINNED[(row.tag, row.kin_smooth)]["in_window"],
                                              P_MAX, int(row.p), d),
        cache_params=lambda tag, v: FIT_PARAMS[(tag, *v)], n_jobs=N_JOBS, cache_name="01var",
        cache_mode=CACHE_MODE, save_cache=SAVE_CACHE, label=_label)

RES = runs.fdr_and_classes(RES, SETTING_COLS, FDR_Q)
runs.print_direction_summary(RES, SETTING_COLS, FDR_Q, len(SESSIONS_DATA), _label)
FLOOR = runs.print_floor_summary(RES, SETTING_COLS, N_SHUFFLE, FDR_Q, _label)
EFFECT_COLS = {"xy": ("pr2x_xy", "excess partial R² of the neuron's rate (%)"),
               "yx": ("pr2x_yx", "excess partial R² of med-lat (%)")}
runs.print_effect_summary(RES, SETTING_COLS, REGIONS_PLOT, EFFECT_COLS, _label)

if LAG_READOUT:
    print(f"\n=== {models.LAG_NOTE} ===")
    runs.print_lag_summary(RES, SETTING_COLS, REGIONS_PLOT, ("centroid", "peak"), "PREVIEW LAG",
                           _label, unreliable=("peak",))

if len(ORDERS) > 1:                              # order-sensitive neurons
    _b = RES[RES.order == "bic"].set_index(["uid", "kin_smooth"])
    _m = RES[RES.order == "pmax"].set_index(["uid", "kin_smooth"]).loc[_b.index]
    flip = (_b.sig_xy != _m.sig_xy) | (_b.sig_yx != _m.sig_yx)
    RES = RES.merge(flip.rename("order_sensitive").reset_index(), on=["uid", "kin_smooth"])
    _fl = flip[flip].reset_index().merge(RES[["uid", "region"]].drop_duplicates(), on="uid")
    print(f"\norder-sensitive neurons (FDR significance differs bic vs pmax): "
          f"{int(flip.sum())}/{len(flip)}  by region: {_fl.region.value_counts().to_dict()}")

if DATA_SOURCE == "sim":
    DET = runs.sim_truth_check(RES, SETTING_COLS, _label,
                               lag_keys=("centroid", "peak") if LAG_READOUT else (),
                               true_lag=SIM_PARAMS["lag_ms"])

kern = {(tag, ks, o, S["neurons"][i]["uid"]): {d: r[d].get("lag") for d in ("xy", "yx")}
        for (tag, ks, o), F in FITS.items() for S in [dict(SESSIONS_DATA)[tag]]
        for i, r in zip(F["idx"], F["res"])} if LAG_READOUT else None

if SAVE_VARS:
    runs.save_tables(dict(NT=NT_ALL, kernels=kern, lag_note=models.LAG_NOTE, params=dict(
        DATA_SOURCE=DATA_SOURCE, SESS_LIST=SESS_LIST, SIM_SCENARIO=SIM_SCENARIO,
        SIM_PARAMS=SIM_PARAMS, SIM_SEEDS=SIM_SEEDS, BIN_MS=BIN_MS, FR_KERNEL=FR_KERNEL,
        KIN_SMOOTHS=KIN_SMOOTHS, ORDERS=ORDERS, REMOVE_ENSEMBLE_MEAN=REMOVE_ENSEMBLE_MEAN,
        P_MAX=P_MAX, N_SHUFFLE=N_SHUFFLE, FDR_Q=FDR_Q, MIN_RATE_HZ=MIN_RATE_HZ,
        MIN_STRATUM=MIN_STRATUM, N_BASES=N_BASES, BASIS_SPACING=BASIS_SPACING, N_BOOT=N_BOOT,
        LAG_TAG=LAG_TAG)), RES, LAG_TAG, OUT_SUBDIR, RUN_NAME)

# %% ── 5a/5b. direction classes (CP vs CB) and direction map per region ─────
runs.plot_direction_figures(RES, SETTING_COLS, REGIONS_PLOT, N_SHUFFLE, FDR_Q, FIG_DIR, _fig,
                            SAVE_FIGS, SHOW_FIGS, f"01 VAR  |  {RUN_NAME}", _label)


# %% ── 5c. PREVIEW lag distributions, CP vs CB (grain of salt) ──────────────
_SIM_CENTROID = {"log": "52 (x→y) / 43 (y→x; 42 with 10 leaders)", "uniform": "65 (x→y) / 47 (y→x)"}  # 6 bases, matched
if LAG_READOUT:
    lag_bins = np.arange(0, config.CROSS_LAG_RANGE_MS[1] + BIN_MS, BIN_MS) + BIN_MS / 2
    for ks, o in SETTINGS:
        R = RES[(RES.kin_smooth == ks) & (RES.order == o)]
        fig, axs = plt.subplots(2, 2, figsize=(11, 7))
        for c_, d in enumerate(("xy", "yx")):
            for r_, key in enumerate(("centroid", "peak")):
                plotting.lag_hist(axs[r_, c_], R, d, key, REGIONS_PLOT,
                                  lag_bins if key == "peak" else np.linspace(0, config.CROSS_LAG_RANGE_MS[1], 16))
                if key == "peak":
                    axs[r_, c_].set_title("peak: UNRELIABLE in simulations (lands at 10-30 ms)", fontsize=8)
                else:
                    axs[r_, c_].set_title(f"centroid ({BASIS_SPACING} bases): sims with true 50 ms read "
                                          f"{_SIM_CENTROID.get(BASIS_SPACING, 'not measured')}", fontsize=8)
        fig.suptitle("PREVIEW LAGS — biased, grain of salt (Gaussian VAR on smoothed rates; "
                     f"final lags from 03)\n{ks} kin, order {o}, FDR-significant neurons only", fontsize=9)
        plotting.finish(fig, FIG_DIR, _fig("3_lags_PREVIEW", ks, o), SAVE_FIGS, SHOW_FIGS)


# %% ── 5d. example cross-kernels behind the lag numbers ──────────────────────
if LAG_READOUT:
    for ks, o in SETTINGS:
        R = RES[(RES.kin_smooth == ks) & (RES.order == o)]
        fig, axs = plt.subplots(len(REGIONS_PLOT), 2, figsize=(10, 3.2 * len(REGIONS_PLOT)), squeeze=False)
        for r_, reg in enumerate(REGIONS_PLOT):
            for c_, d in enumerate(("xy", "yx")):
                ax = axs[r_, c_]
                g = R[(R.region == reg) & R[f"sig_{d}"]].nsmallest(N_EXAMPLE_KERNELS, f"q_{d}")
                for j, (_, row) in enumerate(g.iterrows()):
                    L = kern[(row.tag, ks, o, row.uid)][d]
                    c = plt.cm.tab10(j)
                    ax.plot(L["tau_ms"], L["beta"], color=c, lw=1.3,
                            label=f"{row.uid.split('|', 1)[1]}  centroid {L['centroid_ms']:.0f} ms")
                    if L["beta_lo"] is not None:
                        ax.fill_between(L["tau_ms"], L["beta_lo"], L["beta_hi"], color=c, alpha=0.15)
                ax.axhline(0, color="k", lw=0.4)
                ax.set(xlabel="lag (ms)", ylabel="cross-kernel β(τ)",
                       title=f"{reg}: {plotting.DIR_LABEL[d]} (top {N_EXAMPLE_KERNELS} by q)")
                ax.legend(fontsize=6, frameon=False)
        fig.suptitle(f"PREVIEW cross-kernels with bootstrap 95% bands, {ks} kin, order {o}", fontsize=9)
        plotting.finish(fig, FIG_DIR, _fig("4_example_kernels", ks, o), SAVE_FIGS, SHOW_FIGS)


# %% ── 5e. only when both settings ran: raw vs matched, bic vs pmax ──────────
def _paired(ka, oa, kb, ob):
    A = RES[(RES.kin_smooth == ka) & (RES.order == oa)].set_index("uid")
    return A, RES[(RES.kin_smooth == kb) & (RES.order == ob)].set_index("uid").loc[A.index]

if set(KIN_SMOOTHS) >= {"raw", "matched"}:
    for o in ORDERS:
        A, M = _paired("raw", o, "matched", o)
        fig, axs = plt.subplots(1, 2, figsize=(10, 4.4))
        for ax, d in zip(axs, ("xy", "yx")):
            plotting.z_compare(ax, A[f"z_{d}"], M[f"z_{d}"], A, "null-z, raw kinematics",
                               "null-z, matched kinematics", plotting.DIR_LABEL[d])
        plotting.finish(fig, FIG_DIR, _fig("5_raw_vs_matched", "", o), SAVE_FIGS, SHOW_FIGS)
if set(ORDERS) >= {"bic", "pmax"}:
    for ks in KIN_SMOOTHS:
        Bb, Pm = _paired(ks, "bic", ks, "pmax")
        fig, axs = plt.subplots(1, 2, figsize=(10, 4.4))
        for ax, d in zip(axs, ("xy", "yx")):
            plotting.z_compare(ax, Bb[f"z_{d}"], Pm[f"z_{d}"], Bb, "null-z, BIC order",
                               f"null-z, p = {P_MAX}", plotting.DIR_LABEL[d])
        plotting.finish(fig, FIG_DIR, _fig("5_bic_vs_15lags", ks), SAVE_FIGS, SHOW_FIGS)


# %% ── 5f. simulations only: detection vs ground truth ───────────────────────
if DATA_SOURCE == "sim":
    runs.plot_detection(DET, SIM_SCENARIO, FIG_DIR, _fig("6_detection"), SAVE_FIGS, SHOW_FIGS)
    runs.plot_sim_truth(RES, SETTING_COLS, ("centroid", "peak") if LAG_READOUT else (),
                        SIM_PARAMS["lag_ms"], FIG_DIR, _fig, SAVE_FIGS, SHOW_FIGS, _label)

# %% ── 5g. mean cross-kernel over neurons, per region x direction (PREVIEW) ───
# Each neuron's kernel beta(tau) is divided by its own max |beta| (sign kept) before
# averaging: raw amplitudes are in different units per neuron (rate vs mm) and a
# few strong cells would otherwise dominate. Solid = FDR-significant neurons in that
# direction (mean ± SEM); dashed = all kept neurons. The centroid of |mean kernel| is
# printed and marked. Same caveat as every lag here (models.LAG_NOTE).
if LAG_READOUT:
    def _norm_kernels(R, d):
        K = np.array([kern[(row.tag, row.kin_smooth, row.order, row.uid)][d]["beta"]
                      for row in R.itertuples()]) if len(R) else np.zeros((0, len(BASIS_LAGS)))
        mx = np.max(np.abs(K), axis=1, keepdims=True) if len(K) else 1.0
        return K / np.where(mx > 0, mx, 1.0)

    tau = BASIS_LAGS * BIN_MS
    for ks, o in SETTINGS:
        R = RES[(RES.kin_smooth == ks) & (RES.order == o)]
        fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
        for ax, d in zip(axs, ("xy", "yx")):
            for reg in REGIONS_PLOT:
                c = config.CLR.get(reg, (0.4, 0.4, 0.4))
                for sel, ls, lab in [((R.region == reg) & R[f"sig_{d}"], "-", "significant"),
                                     ((R.region == reg), "--", "all kept")]:
                    K = _norm_kernels(R[sel], d)
                    if len(K) == 0:
                        continue
                    m = K.mean(0); se = K.std(0, ddof=1) / np.sqrt(len(K)) if len(K) > 1 else 0 * m
                    ax.plot(tau, m, color=c, ls=ls, lw=1.8 if ls == "-" else 1.0,
                            label=f"{reg} {lab} (n={len(K)})")
                    if ls == "-":
                        ax.fill_between(tau, m - se, m + se, color=c, alpha=0.2)
                        a = np.abs(m); cen = np.sum(tau * a) / a.sum() if a.sum() else np.nan
                        ax.axvline(cen, color=c, lw=1, ls=":")
                        print(f"  MEAN KERNEL [{ks}, {o}] {plotting.DIR_LABEL[d]}, {reg}: "
                              f"n={len(K)} significant, centroid of |mean| = {cen:.0f} ms (PREVIEW)")
            ax.axhline(0, color="k", lw=0.4)
            ax.set(xlabel="lag (ms)", ylabel="mean normalised β(τ)  (± SEM)", title=plotting.DIR_LABEL[d])
            ax.legend(fontsize=6, frameon=False)
        fig.suptitle(f"PREVIEW mean cross-kernel over neurons (each scaled to max |β| = 1), {ks} kin, "
                     f"order {o}, {BASIS_SPACING} bases — grain of salt; dotted = centroid", fontsize=9)
        plotting.finish(fig, FIG_DIR, _fig("7_mean_kernel", ks, o), SAVE_FIGS, SHOW_FIGS)

# %% ── 5h. effect STRENGTH and FRACTION of neurons with an effect ────────────
# p-values say how SURE we are, not how STRONG: they mix effect size with the amount of
# data and saturate at the shuffle floor. A = p-value histograms + pi1 (estimated
# fraction of neurons with an effect, incl. ones too weak for FDR); B = excess partial R^2
# (effect minus its shuffle mean; 0 = chance); C = strength vs evidence. Compare regions
# within a direction (row) only.
runs.plot_effect_figures(RES, SETTING_COLS, REGIONS_PLOT, EFFECT_COLS, FIG_DIR, _fig,
                         SAVE_FIGS, SHOW_FIGS, f"01 VAR  |  {RUN_NAME}", _label)

# observed partial R^2 vs the trial-shuffle null (= chance level), per region (rows) x
# direction (columns); shuffles converted with partial R^2 = 1 - exp(-GC). Intuition only.
for ks, o in SETTINGS:
    PAIRS = {}
    for tag, S in SESSIONS_DATA:
        F = FITS[(tag, ks, o)]
        for i, r in zip(F["idx"], F["res"]):
            reg = S["neurons"][i]["region"]
            for d in ("xy", "yx"):
                PAIRS.setdefault((reg, d), []).append(
                    (1 - np.exp(-r[d]["gc"]), 1 - np.exp(-np.asarray(r[d]["null"], float))))
    runs.plot_null_vs_observed(PAIRS, REGIONS_PLOT,
                               {"xy": "partial R² of the neuron rate", "yx": "partial R² of med-lat"},
                               FIG_DIR, _fig("8b_null_vs_observed", ks, o), SAVE_FIGS, SHOW_FIGS,
                               f"01 VAR  |  {RUN_NAME}  |  {_label((ks, o))}")

# %%
