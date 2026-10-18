# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 03b -- EXPLAINER figures for script 03 (how the model works; not results).
#   M3.0a  pi1 with a neuron- vs session-bootstrap CI           (needs 01/02 neuron CSVs)
#   M3.0b  what adaptive shuffles change                         (needs a 02 y_to_x sim CSV)
#   M3.1a  every basis on one lag axis (cross, own history, event)
#   M3.1b  one real trial unrolled into the design columns
#   M3.1c  the optional task-locked nuisance mean vs single trials
#   M3.2   script 03 on one real neuron per region: CV over H, fit, kernel + bootstrap band,
#          shuffle null, adequacy test
# Each cell runs on its own after cell 0 (cells 3-5 share the session loaded in cell 2).
# Saved under <OUT_DIR>/<LAG_TAG>/03_KERNELS/explainers/.
import os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
import importlib
from gctools import config, io, preprocess as pp, design as ds, kernels, nulls, stats, plotting, cache
from gctools.models import poisson_irls
for _m in (config, io, pp, ds, kernels, nulls, stats, plotting, cache):
    importlib.reload(_m)

SESSION    = config.FIRST_PASS[0]          # session used for M3.1b/c and M3.2
BIN_MS     = 2.5
KIN_SOURCE = "long"                       # as script 03 (REACH_LONG history; run 00b first)
SPEC       = kernels.make_spec(bin_ms=BIN_MS, h_kin_list=(150, 300, 450), kin_hist_mode="cv")  # script 03 defaults
N_SHUFFLE  = 500                           # M3.2 shuffle null
N_BOOT     = 100                           # M3.2 lag bootstrap
NEURON_UIDS = None                         # M3.1b / M3.2: None = pick from REF_CSV, or
                                           # {"CB": uid, "CP": uid} to choose yourself
SHOW_FIGS  = True
SAVE_FIGS  = True

V = config.VAR_DIR / config.LAG_TAG
# neuron CSVs of YOUR earlier runs (None = skip the figure that needs it)
REF_CSV   = V / "02_GLM" / "real_22sessions_poissonWindows_neurons.csv"     # picks example neurons
REAL_CSVS = {"REAL 01 (VAR)": V / "01_VAR" / "real_22sessions_kinMatched_lags15_neurons.csv",
             "REAL 02 (GLM)": V / "02_GLM" / "real_22sessions_poissonWindows_neurons.csv"}
SIM_CSV   = None   # M3.0a/b: a 02 neurons CSV from a y_to_x SIMULATION run with adaptive
                   # shuffles, e.g. V / "02_GLM" / "<sim run name>_neurons.csv"

OUT = config.fig_dir(config.LAG_TAG, "03_KERNELS", "explainers") if SAVE_FIGS else None
_exists = lambda f: f is not None and os.path.exists(cache._long(f))
_read = lambda f: pd.read_csv(cache._long(f))
LCRC = {1: "tab:blue", 2: "tab:orange", 3: "tab:red"}


# %% ── M3.0a. pi1: neuron vs session bootstrap ───────────────────────────────
sets = []
if _exists(SIM_CSV):
    R = _read(SIM_CSV)
    sets.append((f"SIMULATION ({R.tag.nunique()} sessions,\n{int(R.exp_yx.sum())} true leaders / {len(R)})",
                 R, ["SIM"], {"xy": R.exp_xy.mean(), "yx": R.exp_yx.mean()}))
for lab, f in REAL_CSVS.items():
    if _exists(f):
        R = _read(f)
        if "kin_smooth" in R:
            R = R[(R.kin_smooth == "matched") & (R.order == "pmax")]
        sets.append((f"{lab}, {R.tag.nunique()} sessions", R, ["CP", "CB"], None))
if not sets:
    print("M3.0a skipped: none of SIM_CSV / REAL_CSVS exists")
else:
    fig, axs = plt.subplots(1, len(sets), figsize=(5.2 * len(sets), 4.3), squeeze=False)
    for ax, (lab, R, regs, truth) in zip(axs[0], sets):
        x = 0; ticks, tl = [], []
        for d in ("xy", "yx"):
            for reg in regs:
                g = R[R.region == reg]
                pi, lo1, hi1 = stats.storey_pi1(g[f"pS_{d}"], n_boot=1000)
                _, lo2, hi2 = stats.storey_pi1(g[f"pS_{d}"], n_boot=1000, groups=g["tag"])
                c = config.CLR.get(reg, (0.4, 0.4, 0.4))
                ax.errorbar(x - 0.12, pi, yerr=[[pi - lo1], [hi1 - pi]], fmt="o", color=c, mfc="w", capsize=4,
                            label="neuron bootstrap" if x == 0 else None)
                if np.isfinite(lo2):
                    ax.errorbar(x + 0.12, pi, yerr=[[pi - lo2], [hi2 - pi]], fmt="o", color=c, capsize=4, lw=2,
                                label="session bootstrap (used)" if x == 0 else None)
                if truth is not None:
                    ax.plot([x - 0.35, x + 0.35], [truth[d]] * 2, color="tab:green", lw=2,
                            label="true fraction" if x == 0 else None)
                ticks.append(x); tl.append(f"{reg}\n{plotting.DIR_LABEL[d]}")
                x += 1
        ax.set_xticks(ticks, tl, fontsize=8)
        ax.set(ylim=(-0.02, 1.02), ylabel="π₁ = estimated fraction of neurons with an effect", title=lab)
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("π₁ and its 95% CI: resampling NEURONS (open) vs whole SESSIONS (filled). Neurons of a session "
                 "share kinematics and shuffles,\nso their p-values are correlated -> the session bootstrap is the "
                 "honest one (wider when sessions differ)", fontsize=9)
    plotting.finish(fig, OUT, "M3.0a_pi1_session_bootstrap", SAVE_FIGS, SHOW_FIGS)


# %% ── M3.0b. adaptive shuffles: before / after (one simulated session) ──────
if not (_exists(SIM_CSV) and "n_shuf_yx" in _read(SIM_CSV)):
    print("M3.0b skipped: set SIM_CSV to a 02 y_to_x simulation run with ADAPTIVE_SHUFFLE = True")
else:
    A = _read(SIM_CSV)
    A = A[A.tag == A[A.driver].tag.iloc[0]].reset_index(drop=True)        # one session with a leader
    N0 = int(A.n_shuf_yx.min())
    before = A.copy()
    for d in ("xy", "yx"):
        ext = before[f"n_shuf_{d}"] > N0
        before.loc[ext, f"pS_{d}"] = 1 / (N0 + 1)                         # what it was without extension
        before[f"sig_{d}"] = stats.bh_fdr(before[f"pS_{d}"].to_numpy(), 0.05)[0]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (lab, T) in zip(axs, ((f"BEFORE: {N0} shuffles only", before), ("AFTER: adaptive shuffles", A))):
        lp = -np.log10(T["pS_yx"]); leader = T["driver"].astype(bool).to_numpy(); sig = T["sig_yx"].astype(bool).to_numpy()
        ax.scatter(np.arange(len(T))[~leader], lp[~leader], s=12, color="0.6", label="null neurons")
        ax.scatter(np.arange(len(T))[leader], lp[leader], s=80, color="tab:red", marker="*",
                   label=f"true leader(s) ({int(sig[leader].sum())}/{leader.sum()} FDR-significant)")
        thr = 0.05 / len(T)
        ax.axhline(-np.log10(1 / (N0 + 1)), color="k", ls=":", lw=1, label=f"floor with {N0} shuffles: p = 1/{N0+1}")
        ax.axhline(-np.log10(thr), color="tab:blue", ls="--", lw=1, label=f"BH threshold for a lone neuron: q/m = {thr:.1e}")
        ax.set(xlabel="neuron", ylabel="-log10 shuffle p (neuron → kin)", title=lab,
               ylim=(0, max(4, float(lp.max()) + 0.5)))
        ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.suptitle(f"Adaptive shuffles (ONE simulated session, {len(A)} neurons): with {N0} shuffles a leader cannot go "
                 f"below p = 1/{N0+1}, which can be above what BH needs\nfor a lone effect. Extending only those "
                 "neurons' shuffles (same sequence) resolves it.", fontsize=9)
    plotting.finish(fig, OUT, "M3.0b_adaptive_shuffles", SAVE_FIGS, SHOW_FIGS)


# %% ── M3.1a. every basis on one lag axis ────────────────────────────────────
CROSS_B, CROSS_L = SPEC["B_cross"], SPEC["L_cross"]
EV_B, EV_L, _ = ds.raised_cosine_basis(4, (0, 125), BIN_MS, "log")
HB = {H: kernels.own_bases(H, SPEC) for H in sorted(set(SPEC["h_list"]) | set(SPEC["h_kin_list"]))}
fig, axs = plt.subplots(1, 3, figsize=(16, 4))
ax = axs[0]
for k in range(CROSS_B.shape[1]):
    ax.plot(CROSS_L * BIN_MS, CROSS_B[:, k], lw=1.6)
ax.set(title=f"CROSS bases (both directions):\nother signal's past, {SPEC['n_cross']} {SPEC['spacing']} bumps "
             f"over (0, {SPEC['cross_range'][1]}] ms", xlabel="lag into the past (ms)", ylabel="weight", xlim=(0, 310))
ax = axs[1]
for (H, (B_, L_)), c in zip(HB.items(), ("tab:blue", "tab:orange", "tab:green", "tab:purple")):
    for k in range(B_.shape[1]):
        ax.plot(L_ * BIN_MS, B_[:, k], color=c, lw=1.2 if k >= 2 else 2.5,
                ls="-" if k >= 2 else ":", label=f"H = {H} ms: {B_.shape[1] - 2} log bumps" if k == 2 else None)
ax.annotate("lags 1 and 2 (2.5, 5 ms):\none column each\n(refractoriness)", (5, 1), xytext=(60, 0.75),
            fontsize=8, arrowprops=dict(arrowstyle="->"))
ax.set(title=f"OWN-HISTORY bases, H = {' / '.join(map(str, HB))} ms (~2 more bumps per doubling)\n"
             f"neuron {list(SPEC['h_list'])}: CV per neuron; kinematics {list(SPEC['h_kin_list'])}: CV per session",
       xlabel="lag into the past (ms)", xlim=(0, 460))
ax.legend(fontsize=7, frameon=False)
ax = axs[2]
for k in range(EV_B.shape[1]):
    ax.plot(EV_L * BIN_MS, EV_B[:, k], lw=1.6)
ax.axvline(0, color="k", ls=":", lw=1)
ax.set(title="EVENT kernel (optional, EVENT_KERNEL = True):\nactive only AFTER max reach (t = 0), 4 bumps over (0, 125] ms",
       xlabel="time after the event (ms)", xlim=(-20, 140))
fig.suptitle("Script 03 building blocks: each smooth bump becomes ONE regression column = the signal's past "
             "weighted by that bump (all strictly causal: lag >= 1 bin)", fontsize=9)
plotting.finish(fig, OUT, "M3.1a_bases", SAVE_FIGS, SHOW_FIGS)


# %% ── 2. load the session and pick example neurons (for M3.1b/c and M3.2) ──
S = io.load_session(*SESSION, verbose=False, kin_source=KIN_SOURCE)
UIDS = [n["uid"] for n in S["neurons"]]
TAG = f"{SESSION[0].split('_', 1)[-1]}-{SESSION[1]}"
RATE = io.neuron_table(S)["fr_win_hz"].to_numpy()
if NEURON_UIDS is None:
    NEURON_UIDS = {}
    R0 = _read(REF_CSV) if _exists(REF_CSV) else None
    for reg in ("CB", "CP"):
        g = None if R0 is None else R0[(R0.tag == TAG) & (R0.region == reg) & R0.sig_xy
                                       & (R0.rate_hz > 3) & (R0.rate_hz < 80)]
        if g is not None and len(g):
            NEURON_UIDS[reg] = g.sort_values("z_xy", ascending=False).uid.iloc[0]
        else:                                       # no reference run: a mid-rate neuron
            cand = [i for i, n in enumerate(S["neurons"]) if n["region"] == reg and 3 < RATE[i] < 80]
            if cand:
                NEURON_UIDS[reg] = UIDS[cand[len(cand) // 2]]
print(f"{TAG}: example neurons {NEURON_UIDS}"
      + ("" if _exists(REF_CSV) else "  (REF_CSV missing: picked by rate, not by significance)"))
B = pp.bin_session(S, BIN_MS, hist_pad_ms=kernels.pad_needed_ms(SPEC), kin_smooth="raw", fr_kernel=None)


# %% ── M3.1b. one real trial unrolled into the design ────────────────────────
uid = NEURON_UIDS.get("CB") or next(iter(NEURON_UIDS.values()))
ni = UIDS.index(uid)
y, x, t = B["Y"][ni].astype(float), B["X"], B["t"] * 1000
tr = int(np.argmax(np.ptp(np.nan_to_num(x), axis=1)))           # trial with the largest med-lat range
Hb, Hl = HB[150]
own = ds.basis_features(y[tr], Hb, Hl)
crs = ds.basis_features(x[tr], CROSS_B, CROSS_L)
ev = ds.event_regressors(B["t"], 0.0, EV_B, EV_L)
nui = ds.nuisance_mean(x, B["cond"])[tr]
now = 40.0
fig, axs = plt.subplots(5, 1, figsize=(11, 11), sharex=True,
                        gridspec_kw=dict(height_ratios=[0.6, 1, 1.3, 1.3, 1]))
sp = np.asarray(S["neurons"][ni]["spikes"][tr]) * 1000
axs[0].vlines(sp[(sp > t[0]) & (sp < t[-1])], 0, 1, color="k", lw=0.8)
axs[0].set(yticks=[], title=f"one real trial: {uid}, trial {int(S['trials'].trial_idx[tr])} "
                             f"(spout {int(S['trials'].LCR[tr])})")
axs[1].plot(t, x[tr], color="tab:blue")
axs[1].plot(t, nui, color="tab:blue", ls="--", lw=1, alpha=0.7, label="nuisance: LOO mean of same spout")
axs[1].set(ylabel="med-lat (mm)"); axs[1].legend(fontsize=7, frameon=False)
for k in range(own.shape[1]):
    axs[2].plot(t, own[:, k] / (np.nanmax(own[:, k]) or 1) + k, lw=0.9)
axs[2].set(ylabel="own-history\ncolumns (scaled)", yticks=[])
for k in range(crs.shape[1]):
    v = crs[:, k]; axs[3].plot(t, (v - np.nanmean(v)) / (np.nanstd(v) or 1) * 0.3 + k, lw=1)
axs[3].set(ylabel="cross columns\n(past med-lat)", yticks=[])
for k in range(ev.shape[1]):
    axs[4].plot(t, ev[:, k] + 1.1 * k, lw=1.2)
axs[4].set(ylabel="event columns", yticks=[], xlabel="time re max reach (ms)")
for ax in axs:
    ax.axvline(now, color="tab:red", lw=1.2)
    ax.axvspan(now - 150, now, color="tab:red", alpha=0.07)
    ax.axvspan(config.TRIAL_WINDOW_MS[0], config.TRIAL_WINDOW_MS[0] + 150, color="0.85", alpha=0.5, zorder=-2)
    ax.axvline(0, color="k", ls=":", lw=0.8)
axs[0].text(now + 3, 0.5, "'now' (one fitted row)", color="tab:red", fontsize=8, va="center")
axs[0].text(now - 148, 0.5, "past it may use (150 ms)", color="tab:red", fontsize=8, va="center")
axs[1].text(config.TRIAL_WINDOW_MS[0] + 3, np.nanmax(x[tr]),
            f"kinematics start\n({config.TRIAL_WINDOW_MS[0]} ms); fitted rows\nstart at "
            f"{config.TRIAL_WINDOW_MS[0] + 150} ms", fontsize=7, va="top", color="0.3")
fig.suptitle("One trial unrolled into script 03's design: each line is one regression column over time. At the "
             "red 'now', every column only summarises\nthe shaded past (lags >= 1 bin); the event columns switch on "
             "only after max reach; NaN where the past is not available (grey).", fontsize=9)
plotting.finish(fig, OUT, "M3.1b_one_trial_design", SAVE_FIGS, SHOW_FIGS)


# %% ── M3.1c. nuisance mean vs single trials ─────────────────────────────────
fig, axs = plt.subplots(2, 3, figsize=(15, 7), sharex=True)
NU = ds.nuisance_mean(x, B["cond"])
rate = pp.causal_filter(y, pp.gamma_kernel(2.0, 10.0, BIN_MS)[0], pre="zero") / (BIN_MS / 1000)
NUr = ds.nuisance_mean(rate, B["cond"])
inw = B["in_window"]
for c_, lcr in enumerate((1, 2, 3)):
    idx = np.flatnonzero(S["trials"].LCR.to_numpy() == lcr)
    if not idx.size:
        continue
    c = LCRC[lcr]
    for i in idx[:12]:
        axs[0, c_].plot(t[inw], x[i, inw], color=c, lw=0.6, alpha=0.4)
        axs[1, c_].plot(t[inw], rate[i, inw], color=c, lw=0.5, alpha=0.3)
    axs[0, c_].plot(t[inw], NU[idx[0], inw], color="k", lw=2, label="nuisance for trial 1\n(mean of the other trials)")
    axs[1, c_].plot(t[inw], NUr[idx[0], inw], color="k", lw=2)
    axs[0, c_].set(title=f"spout {['left', 'centre', 'right'][c_]} (n={idx.size})")
    axs[1, c_].set(xlabel="time re max reach (ms)")
    for ax in axs[:, c_]:
        ax.axvline(0, color="k", ls=":", lw=0.8)
axs[0, 0].set_ylabel("med-lat (mm)"); axs[1, 0].set_ylabel(f"rate (Hz), {uid.split('|')[-1]}")
axs[0, 0].legend(fontsize=7, frameon=False)
fig.suptitle("Optional task-locked NUISANCE regressor (TASK_MEAN_NUISANCE): for each trial, the mean of the OTHER "
             "trials of the same spout at the same time (black).\nIt soaks up what is common to all trials of a side, "
             "leaving the trial-by-trial deviations (thin lines) for the coupling terms. Never uses the trial's own values.",
             fontsize=9)
plotting.finish(fig, OUT, "M3.1c_nuisance_mean", SAVE_FIGS, SHOW_FIGS)


# %% ── M3.2. script 03 on one real neuron per region ─────────────────────────
# Slow-ish: N_SHUFFLE shuffles + N_BOOT bootstrap refits per neuron (~1-3 min each).
perms = nulls.make_perms(B["cond"], N_SHUFFLE, seed=0)
lcr = S["trials"]["LCR"].to_numpy()
for reg, uid in NEURON_UIDS.items():
    i = UIDS.index(uid)
    y, x, t = B["Y"][i], B["X"], B["t"]
    r = kernels.kernel_neuron(y, x, B["in_window"], B["cond"], t, SPEC, perms)
    H = r["H"]
    bt = {d: kernels.kernel_lag_bootstrap(y, x, B["in_window"], B["cond"], t, SPEC, H if d == "xy" else r["H_kin"], d, n_boot=N_BOOT)
          for d in ("xy", "yx")}
    D = kernels.build_design(y, x, B["in_window"], B["cond"], t, SPEC)
    rows, n_tr = D["rows"], y.shape[0]
    one = np.ones((D["yt"].size, 1))                  # reduced / full fits, for display only
    Xr = np.hstack([one, kernels._z(D[f"own_y{H}"])[0]]); Xc = kernels._z(D["cross_x"])[0]
    mu_r = np.exp(Xr @ poisson_irls(Xr, D["yt"], SPEC["ridge"])["b"])
    mu_f = np.exp(np.hstack([Xr, Xc]) @ poisson_irls(np.hstack([Xr, Xc]), D["yt"], SPEC["ridge"])["b"])
    Zr = np.column_stack([np.ones(D["xt"].size), D[f"own_x{r['H_kin']}"]]); Zf = np.column_stack([Zr, D["cross_y"]])
    xr = Zr @ np.linalg.lstsq(Zr, D["xt"], rcond=None)[0]; xf = Zf @ np.linalg.lstsq(Zf, D["xt"], rcond=None)[0]
    tt = t[rows] * 1000
    sec = BIN_MS / 1000

    fig, axs = plt.subplots(2, 4, figsize=(19, 8.5))
    # row 1: kin -> neuron
    ax = axs[0, 0]
    Hs = list(r["cv_ll_H"]); ll = np.array([r["cv_ll_H"][h] for h in Hs])
    ax.plot(Hs, ll - ll.max(), "o-", color="k"); ax.axvline(H, color="tab:red", ls="--", label=f"chosen H = {H} ms")
    ax.set(xlabel="own-history length H (ms)", ylabel="held-out log-lik − best",
           title=f"(a) CV picks the own history\n(reduced model, {SPEC['n_folds']} trial folds)")
    ax.legend(fontsize=8, frameon=False)
    ax = axs[0, 1]
    for c in (1, 2, 3):
        m = (lcr == c)
        if not m.any():
            continue
        ax.plot(tt, D["yt"].reshape(n_tr, -1)[m].mean(0) / sec, color=LCRC[c], lw=1, alpha=0.5)
        ax.plot(tt, mu_r.reshape(n_tr, -1)[m].mean(0) / sec, color=LCRC[c], ls="--", lw=1.2)
        ax.plot(tt, mu_f.reshape(n_tr, -1)[m].mean(0) / sec, color=LCRC[c], lw=2)
    ax.plot([], [], color="0.4", lw=1, alpha=0.5, label="observed")
    ax.plot([], [], color="0.4", ls="--", label="reduced (own history)")
    ax.plot([], [], color="0.4", lw=2, label="full (+ kinematics)")
    ax.axvline(0, color="k", ls=":", lw=0.8)
    ax.set(xlabel="time re max reach (ms)", ylabel="rate (Hz), mean per spout",
           title=f"(b) fit, per spout side (colours)\nCV gain = {r['xy']['cv_gain_bits']:.4f} bits/spike")
    ax.legend(fontsize=7, frameon=False)
    # kernels
    for ax, d in ((axs[0, 2], "xy"), (axs[1, 2], "yx")):
        R = r[d]; tau = R["tau_ms"]
        ax.fill_between(tau, bt[d]["beta_lo"], bt[d]["beta_hi"], color=config.CLR.get(reg, (0.4, 0.4, 0.4)),
                        alpha=0.25, label=f"95% trial bootstrap ({N_BOOT})")
        ax.plot(tau, R["beta"], color="k", lw=2, label="fitted kernel β(τ)")
        ax.axhline(0, color="k", lw=0.5)
        ax.axvline(R["peak_ms"], color="tab:red", ls="--", lw=1,
                   label=f"peak {R['peak_ms']:.0f} ms [{bt[d]['peak_ci'][0]:.0f}, {bt[d]['peak_ci'][1]:.0f}]")
        ax.axvline(R["centroid_ms"], color="tab:blue", ls=":", lw=1.4,
                   label=f"centroid {R['centroid_ms']:.0f} ms [{bt[d]['centroid_ci'][0]:.0f}, {bt[d]['centroid_ci'][1]:.0f}]")
        ax.set(xlabel="lag (ms)", ylabel="β(τ)  (" + ("per mm" if d == "xy" else "mm per spike") + ")",
               title=f"(c) {plotting.DIR_LABEL[d]} kernel (PREVIEW lag)")
        ax.legend(fontsize=7, frameon=False)
    # shuffle nulls
    for ax, d, stat in ((axs[0, 3], "xy", "Poisson score statistic"), (axs[1, 3], "yx", "ln RSS_red / RSS_full")):
        R = r[d]; obs = R["score"] if d == "xy" else R["gc"]
        ax.hist(R["null"], bins=40, color="0.6", label=f"{len(R['null'])} within-side shuffles")
        ax.axvline(obs, color="tab:red", lw=2, label=f"observed  (p = {R['p_shuf']:.3g})")
        ax.set(xlabel=stat, ylabel="# shuffles", title="(d) significance: observed vs its own shuffle null")
        ax.legend(fontsize=7, frameon=False)
        if obs > 3 * np.max(R["null"]) and np.min(R["null"]) > 0:
            ax.set_xscale("log")
    # row 2: adequacy text + one-trial prediction
    ax = axs[1, 0]
    ax.axis("off")
    ax.text(0.0, 0.95, "(a) adequacy: was own history long enough?", fontsize=11, va="top", weight="bold")
    for yy, lab, Hh, p_ in ((0.75, "kin → neuron", H, r["xy"]["adequacy_p"]),
                            (0.40, "neuron → kin", r["H_kin"], r["yx"]["adequacy_p"])):
        what = "spike" if lab.startswith("kin") else "med-lat"
        ax.text(0.0, yy, f"{lab}: own history {Hh} ms\n  add {what} history ({Hh}, {Hh + SPEC['adequacy_ms']}] ms"
                         f" -> p = {p_:.3g}\n  {'FLAG: older history still predicts' if p_ < 0.01 else 'OK'}",
                fontsize=9, va="top", family="monospace")
    ax = axs[1, 1]
    trk = int(np.argmax(np.ptp(D["xt"].reshape(n_tr, -1), axis=1)))
    ax.plot(tt, D["xt"].reshape(n_tr, -1)[trk], color="k", lw=2, label="observed med-lat")
    ax.plot(tt, xr.reshape(n_tr, -1)[trk], color="0.5", ls="--", label="reduced: own past only")
    ax.plot(tt, xf.reshape(n_tr, -1)[trk], color="tab:red", lw=1, label="full: + neuron's past")
    ax.set(xlabel="time re max reach (ms)", ylabel="med-lat (mm)",
           title=f"(b) one trial, one-step-ahead predictions\nCV gain = {100*r['yx']['cv_gain_r2']:.3f} % of held-out variance")
    ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"Script 03 on a real {reg} neuron: {uid}  ({n_tr} trials, {RATE[i]:.1f} Hz)  —  "
                 "top: kinematics → neuron (Poisson), bottom: neuron → kinematics (linear)", fontsize=10)
    plotting.finish(fig, OUT, f"M3.2_{reg}_neuron_{uid.split('|')[-1]}", SAVE_FIGS, SHOW_FIGS)
    print(f"{reg} {uid}: H {H} ms | p_shuf kin->neuron {r['xy']['p_shuf']:.3g}, neuron->kin {r['yx']['p_shuf']:.3g}"
          f" | adequacy p {r['xy']['adequacy_p']:.2g} / {r['yx']['adequacy_p']:.2g}")

# %%
