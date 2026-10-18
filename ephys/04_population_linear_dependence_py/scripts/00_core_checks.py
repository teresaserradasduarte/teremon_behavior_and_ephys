# %% ── 0. imports & parameters ────────────────────────────────────────────────
# Milestone-1 EVIDENCE for the gctools core (the pass/fail asserts are in tests/):
#   1. causal kernels: gamma smoothing kernel + raised-cosine lag bases
#   2. raw vs "matched" (causally smoothed) med-lat: the smoothing DELAY made visible
#   3. rasters + med-lat by spout side, one real CP and one real CB neuron, and sims
#   4. spike-triggered past kinematics on REAL kinematics: null vs x_to_y
#   5. session table (trials per LCR stratum, neurons per region)
# t=0 is the MAX-REACH frame everywhere (not reach onset).
import sys, time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# project root = the folder holding gctools/. Found from this file's path (python run:
# __file__; VS Code interactive window: __vsc_ipynb_file__), else from the cwd -- up
# the tree, or down into 04_population_linear_dependence_py (e.g. cwd = repo root).
_here = globals().get("__vsc_ipynb_file__") or globals().get("__file__")
_starts = ([Path(_here).resolve().parent] if _here else []) + [Path.cwd()]
_cands = [q for s in _starts for p in [s, *s.parents]
          for q in (p, p / "04_population_linear_dependence_py",
                    p / "ephys" / "04_population_linear_dependence_py")]
ROOT = next((p for p in _cands if (p / "gctools").is_dir()), None)
if ROOT is None:
    raise RuntimeError(f"gctools/ not found from {_starts}; open the interactive window "
                       f"from the 04_population_linear_dependence_py folder")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import importlib
from gctools import config, io, preprocess as pp, design as ds, sim, cache
for _m in (config, io, pp, ds, sim, cache):
    importlib.reload(_m)

# ── data ─────────────────────────────────────────────────────────────────────
SESSION     = config.FIRST_PASS[0]     # ("5_FerreroRocher", "R1")
REACH_CAT   = config.REACH_CAT         # (1,)
LAG_TAG     = config.LAG_TAG           # "lag150_win-375to125" -- first folder level

# ── kernels ──────────────────────────────────────────────────────────────────
FR_KERNEL   = (2.0, 10.0)   # causal gamma (shape, scale_ms) for rates / "matched" kin;
                            # mean delay = shape*scale = 20 ms
N_BASES     = 6             # raised-cosine cross bases over CROSS_LAG_RANGE_MS
BIN_MS      = 2.5

# ── simulation (for panels 3-4) ──────────────────────────────────────────────
SIM_RATES   = "match"       # target rates drawn from the REAL neurons of SESSION
SIM_LAG_MS  = 50.0
SIM_COUPLING = 1.0
N_SHUF      = 30            # within-LCR shuffles for the STA band

SHOW_FIGS   = True
SAVE_FIGS   = True
OUT_SUBDIR  = "00_core_checks"
LCR_CLR     = {1: "tab:blue", 2: "tab:orange", 3: "tab:red"}
LCR_NAME    = {1: "left", 2: "centre", 3: "right"}

def _finish(fig, name):
    for ax in fig.axes:
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    if SAVE_FIGS:
        cache.savefig_retry(fig, config.fig_dir(LAG_TAG, OUT_SUBDIR) / f"{name}.png")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

if "SESS" not in globals():
    SESS = io.load_session(*SESSION, reach_cat=REACH_CAT)
print(f"figures -> {config.fig_dir(LAG_TAG, OUT_SUBDIR)}")


# %% ── 1. causal kernels ─────────────────────────────────────────────────────
k, mean_delay = pp.gamma_kernel(*FR_KERNEL, BIN_MS)
fig, axs = plt.subplots(1, 3, figsize=(12, 3))
tk = np.arange(k.size) * BIN_MS
axs[0].bar(tk, k, width=BIN_MS, align="edge", color="0.5")
axs[0].axvline(0, color="k", ls=":", lw=0.8)
axs[0].axvline(mean_delay, color="tab:red", ls="--", lw=1, label=f"mean delay {mean_delay:.0f} ms")
axs[0].set(xlabel="lag (ms, past ->)", ylabel="weight",
           title=f"causal gamma kernel a={FR_KERNEL[0]:g}, th={FR_KERNEL[1]:g} ms")
axs[0].legend(fontsize=7)
for ax, spacing in zip(axs[1:], ("log", "uniform")):
    B, lags, c = ds.raised_cosine_basis(N_BASES, config.CROSS_LAG_RANGE_MS, BIN_MS, spacing)
    ax.plot(lags * BIN_MS, B)
    ax.axvline(0, color="k", ls=":", lw=0.8)
    ax.set(xlabel="lag (ms)", xlim=(0, config.CROSS_LAG_RANGE_MS[1] + 5),
           title=f"raised-cosine, {spacing}  (first lag {lags[0]*BIN_MS:g} ms)")
_finish(fig, "1_kernels")


# %% ── 2. raw vs matched med-lat (smoothing delay) ───────────────────────────
t_ms = SESS["kin"]["t"] * 1000
X = SESS["kin"]["X"]
Xm = pp.causal_filter(X, k, pre="hold")
lag_xc = []
for i in range(len(X)):                             # delay by cross-correlation of velocity
    a, b = np.diff(X[i]), np.diff(Xm[i])
    xc = [np.dot(a[:a.size - L], b[L:]) for L in range(40)]
    lag_xc.append(np.argmax(xc) * BIN_MS)
fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))
i = int(np.argmax(np.ptp(X, axis=1)))
axs[0].plot(t_ms, X[i], color="k", lw=1.2, label="raw")
axs[0].plot(t_ms, Xm[i], color="tab:red", lw=1.2, label="matched (causal gamma)")
axs[0].axvline(0, color="k", ls=":", lw=0.8)
axs[0].set(xlabel="time re max reach (ms)", ylabel="med-lat (mm)",
           title=f"trial {SESS['trials']['trial_idx'][i]}")
axs[0].legend(fontsize=7)
axs[1].hist(lag_xc, bins=np.arange(0, 60, BIN_MS), color="0.5")
axs[1].axvline(mean_delay, color="tab:red", ls="--", label=f"kernel mean {mean_delay:.0f} ms")
axs[1].set(xlabel="delay matched vs raw (ms, xcorr of velocity peak)", ylabel="# trials")
axs[1].legend(fontsize=7)
_finish(fig, "2_raw_vs_matched")
print(f"matched-vs-raw delay: median {np.median(lag_xc):.1f} ms (kernel mean {mean_delay:.0f} ms, "
      f"xcorr peak ~ kernel MODE {(FR_KERNEL[0]-1)*FR_KERNEL[1]:.0f} ms)")


# %% ── 3. rasters + med-lat by spout side (real and simulated) ───────────────
NT = io.neuron_table(SESS)
print(NT.groupby("region")[["fr_win_hz", "frac_empty", "cv2"]].describe(percentiles=[.5]).round(2).T)

t0 = time.time()
SIM1 = sim.simulate_session(SESS, "x_to_y", rate_hz=SIM_RATES, coupling=SIM_COUPLING,
                            lag_ms=SIM_LAG_MS, seed=0)
SIM0 = sim.simulate_session(SESS, "null_independent", rate_hz=SIM_RATES, seed=0)
print(f"simulated {len(SIM1['neurons'])} neurons x {len(SIM1['trials'])} trials x 2 scenarios "
      f"in {time.time()-t0:.1f}s")

def _raster(ax, spikes, lcr, title):
    order = np.argsort(lcr, kind="stable")
    for row, tr in enumerate(order):
        s = spikes[tr]; s = s[(s >= -0.5) & (s < 0.2)] * 1000
        ax.plot(s, np.full(s.size, row), "|", ms=1.5, color=LCR_CLR[int(lcr[tr])])
    for w in config.TRIAL_WINDOW_MS: ax.axvline(w, color="0.5", ls="--", lw=0.8)
    ax.axvline(0, color="k", ls=":", lw=0.8)
    ax.set(xlabel="time re max reach (ms)", ylabel="trial (sorted by LCR)", title=title,
           xlim=(-500, 200))

lcr = SESS["trials"]["LCR"].to_numpy()
pick = {r: NT[NT.region == r].sort_values("fr_win_hz").index[len(NT[NT.region == r]) // 2]
        for r in config.REGIONS}                      # median-rate neuron per region
j_sim = int(np.argmin(np.abs(SIM1["truth"]["rate_hz"] - 20)))
fig, axs = plt.subplots(1, 4, figsize=(16, 3.6))
for ax, (r, ii) in zip(axs, pick.items()):
    _raster(ax, SESS["neurons"][ii]["spikes"], lcr,
            f"{SESS['neurons'][ii]['uid']}  {NT.fr_win_hz[ii]:.1f} Hz")
_raster(axs[2], SIM1["neurons"][j_sim]["spikes"], SIM1["trials"]["LCR"].to_numpy(),
        f"SIM x_to_y L={SIM_LAG_MS:g} ms  {SIM1['truth']['rate_hz'][j_sim]:.1f} Hz")
for c in (1, 2, 3):
    m = lcr == c
    mu, sd = X[m].mean(0), X[m].std(0)
    axs[3].plot(t_ms, mu, color=LCR_CLR[c], label=f"{LCR_NAME[c]} (n={m.sum()})")
    axs[3].fill_between(t_ms, mu - sd, mu + sd, color=LCR_CLR[c], alpha=0.15)
axs[3].axvline(0, color="k", ls=":", lw=0.8)
axs[3].set(xlabel="time re max reach (ms)", ylabel="med-lat (mm)", title="med-lat by spout, mean ± SD")
axs[3].legend(fontsize=7)
_finish(fig, "3_rasters_medlat")


# %% ── 4. spike-triggered PAST med-lat on REAL kinematics: null vs x_to_y ─────
# Real med-lat is smooth and task-locked, so the STA is smeared and has a large
# task-locked component. Shown raw and after removing the leave-one-out LCR mean.
# This is illustration only -- the directional tests are the GLMs of M2+.
LAGS = np.arange(1, int(config.CROSS_LAG_RANGE_MS[1] / BIN_MS) + 1)
rng = np.random.default_rng(0)
fig, axs = plt.subplots(1, 2, figsize=(11, 3.4), sharey=False)
for ax, resid in zip(axs, (False, True)):
    for S_, lab, clr in ((SIM0, "null_independent", "0.4"), (SIM1, "x_to_y", "tab:red")):
        B = pp.bin_session(S_, BIN_MS, hist_pad_ms=config.CROSS_LAG_RANGE_MS[1] + 10)
        rc = B["cond"] if resid else None
        sta = np.mean([sim.sta_kin(B["Y"][n], B["X"], B["in_window"], LAGS, resid_cond=rc)
                       for n in range(B["Y"].shape[0])], axis=0)
        ax.plot(LAGS * BIN_MS, sta, color=clr, lw=1.5, label=lab)
        if lab == "null_independent":
            shuf = np.array([np.mean([sim.sta_kin(B["Y"][n], B["X"], B["in_window"], LAGS,
                                                  pair=sim.perm_within(B["cond"], rng),
                                                  resid_cond=rc)
                                      for n in range(B["Y"].shape[0])], axis=0)
                             for _ in range(N_SHUF)])
            lo, hi = np.percentile(shuf, [0.5, 99.5], axis=0)
            ax.fill_between(LAGS * BIN_MS, lo, hi, color="0.8", label="within-LCR shuffle 99%")
    ax.axvline(SIM_LAG_MS, color="tab:red", ls="--", lw=0.8, label=f"true mode {SIM_LAG_MS:g} ms")
    ax.axhline(0, color="k", lw=0.4)
    ax.set(xlabel="lag of med-lat before the spike (ms)", ylabel="STA − mean (mm)",
           title="raw med-lat" if not resid else "med-lat − LOO LCR mean")
    ax.legend(fontsize=7)
_finish(fig, "4_sta_real_kin")


# %% ── 5. session table ─────────────────────────────────────────────────────
LS = io.list_sessions(reach_cat=REACH_CAT)
print(LS.sort_values("min_stratum", ascending=False).to_string(index=False))
if SAVE_FIGS:
    cat_tag = "".join(map(str, REACH_CAT))
    LS.to_csv(config.fig_dir(LAG_TAG, OUT_SUBDIR) / f"sessions_cat{cat_tag}.csv", index=False)

# %%
