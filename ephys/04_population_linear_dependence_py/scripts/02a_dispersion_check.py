# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 02a -- Poisson DISPERSION of every real neuron, as seen by script 02's model.
# Dispersion = Pearson chi2 / df of the REDUCED Poisson GLM (intercept + own-history
# windows, raw 2.5 ms counts) -- exactly the model and rows 02 uses. It measures how
# variable the counts are AFTER the neuron's own recent history is accounted for:
#     1   = Poisson-like (the model's assumption)
#     < 1 = UNDER-dispersed (more regular than Poisson: e.g. refractoriness inside a bin)
#     > 1 = OVER-dispersed (bursty / extra variability the model misses)
# Why it matters: significance in 02 comes from the trial shuffle, which stays valid
# either way (validated at 100-150 Hz). What dispersion affects is the Poisson-based
# effect size (bits/spike) and the chi2 p -- and whether 03 needs 1 ms bins / another
# count model. No shuffles here, so this is quick.
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
from gctools import config, io, preprocess as pp, models, cache, plotting, runs
from gctools import design as ds
for _m in (config, io, pp, models, cache, plotting, ds, runs):
    importlib.reload(_m)

SESSIONS     = "all"               # "first" | "all" | list of (mouse, sess)
MIN_STRATUM  = 20
REACH_CAT    = config.REACH_CAT
MIN_RATE_HZ  = 1.0
# must match 02_ppglm_discrete_lags.py
BIN_MS           = 2.5
HIST_WINDOWS_MS  = [(0, 2.5), (2.5, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 150)]
CROSS_WINDOWS_MS = [(0, 10), (10, 20), (20, 40), (40, 80), (80, 150)]
RIDGE_ALPHA      = 0.1
BANDS        = (0.8, 1.2)          # dispersion bands for the summary: under / ~Poisson / over

N_JOBS       = 8
SHOW_FIGS    = True
SAVE_FIGS    = True
SAVE_VARS    = True
OUT_SUBDIR   = "02_GLM"
LAG_TAG  = config.lag_tag(config.CROSS_LAG_RANGE_MS, config.TRIAL_WINDOW_MS)
SESS_LIST = runs.session_list("real", SESSIONS, MIN_STRATUM, REACH_CAT, config.FIRST_PASS[0])
RUN_NAME = f"dispersion_{runs.data_name('real', SESS_LIST, None, None, None)}"
FIG_DIR  = config.fig_dir(LAG_TAG, OUT_SUBDIR, RUN_NAME) if SAVE_FIGS else None
print(f"RUN_NAME {RUN_NAME}\nfigures  {FIG_DIR}")


# %% ── 1. data ───────────────────────────────────────────────────────────────
SESSIONS_DATA, NT_ALL = runs.load_data("real", SESS_LIST, REACH_CAT, None, None, None, None, MIN_RATE_HZ)


# %% ── 2. reduced Poisson model per neuron (no shuffles) ──────────────────────
HIST_SETS  = ds.ms_windows_to_lags(HIST_WINDOWS_MS, BIN_MS)
CROSS_SETS = ds.ms_windows_to_lags(CROSS_WINDOWS_MS, BIN_MS)
CENTERS = np.array([(a + b) / 2 for a, b in CROSS_WINDOWS_MS])
rows, t0 = [], time.time()
for tag, S in SESSIONS_DATA:
    B = pp.bin_session(S, BIN_MS, hist_pad_ms=150, kin_smooth="raw", fr_kernel=None)
    idx = np.flatnonzero(S["_keep"])
    res = Parallel(n_jobs=N_JOBS, backend="loky")(
        delayed(models.glm_neuron)(B["Y"][i], B["X"], B["in_window"], HIST_SETS, CROSS_SETS,
                                   CENTERS, None, RIDGE_ALPHA) for i in idx)
    for i, r in zip(idx, res):
        Y = B["Y"][i][:, B["in_window"]]
        rows.append(dict(runs.base_row(tag, S, i), dispersion=r["xy"]["dispersion"],
                         converged=r["xy"]["converged"], n_spikes=r["n_spikes"],
                         frac_bins_2plus=float(np.mean(Y >= 2)),
                         fano_raw=float(Y.var() / max(Y.mean(), 1e-12))))
    print(f"  {tag:18s} {len(idx):3d} neurons  ({time.time()-t0:.0f}s)", flush=True)
DISP = pd.DataFrame(rows)
DISP["band"] = pd.cut(DISP.dispersion, [0, BANDS[0], BANDS[1], np.inf],
                      labels=[f"under (<{BANDS[0]})", f"~Poisson ({BANDS[0]}-{BANDS[1]})", f"over (>{BANDS[1]})"])


# %% ── 3. summary ────────────────────────────────────────────────────────────
print("\nPearson chi2/df of the reduced model (1 = Poisson):")
print(DISP.groupby("region").dispersion.describe(percentiles=[.1, .25, .5, .75, .9]).round(2).to_string())
print("\nfraction of neurons per band:")
print((pd.crosstab(DISP.region, DISP.band, normalize="index") * 100).round(1).to_string())
DISP["rate_bin"] = pd.cut(DISP.rate_hz, [0, 5, 10, 20, 50, 100, np.inf],
                          labels=["<5", "5-10", "10-20", "20-50", "50-100", ">100"])
print("\nmedian dispersion by firing rate (Hz) and region (n in brackets):")
g = DISP.groupby(["rate_bin", "region"], observed=True).dispersion
print(pd.concat([g.median().round(2).rename("median"), g.size().rename("n")], axis=1).unstack("region").to_string())
print(f"\nnot converged: {int((~DISP.converged).sum())} of {len(DISP)}")
if SAVE_VARS:
    p = config.var_path(LAG_TAG, OUT_SUBDIR, f"{RUN_NAME}.csv")
    DISP.to_csv(cache._long(p), index=False); print(f"saved -> {p}")


# %% ── 4. figure ─────────────────────────────────────────────────────────────
REG = ["CP", "CB"]
fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), gridspec_kw=dict(width_ratios=[1.1, 1.3, 0.9]))
ax = axs[0]                                         # A: distribution per region
bins = np.linspace(0, max(2.0, DISP.dispersion.quantile(0.995)), 50)
ax.axvspan(0, BANDS[0], color="tab:blue", alpha=0.07)
ax.axvspan(BANDS[1], bins[-1], color="tab:red", alpha=0.07)
for reg in REG:
    v = DISP[DISP.region == reg].dispersion
    ax.hist(v, bins=bins, alpha=0.6, color=config.CLR[reg], label=f"{reg} n={len(v)}, median {v.median():.2f}")
ax.axvline(1, color="k", lw=1, ls="--")
ax.text(BANDS[0] / 2, ax.get_ylim()[1] * 0.95, "under-\ndispersed", ha="center", va="top", fontsize=7, color="tab:blue")
ax.text((BANDS[1] + bins[-1]) / 2, ax.get_ylim()[1] * 0.95, "over-\ndispersed", ha="center", va="top", fontsize=7, color="tab:red")
ax.set(xlabel="dispersion (Pearson χ²/df, reduced Poisson model)", ylabel="# neurons",
       title="A. how Poisson-like are the neurons?")
ax.legend(fontsize=7, frameon=False, loc="center right")

ax = axs[1]                                         # B: dispersion vs firing rate
for reg in REG:
    g = DISP[DISP.region == reg]
    ax.scatter(g.rate_hz, g.dispersion, s=9, alpha=0.5, color=config.CLR[reg], label=reg)
    rb = g.groupby("rate_bin", observed=True)
    ax.plot(rb.rate_hz.median(), rb.dispersion.median(), color=config.CLR[reg], lw=2, marker="o",
            markeredgecolor="k", label=f"{reg} median per rate bin")
ax.axhline(1, color="k", lw=1, ls="--")
ax.axhspan(BANDS[0], BANDS[1], color="0.9", zorder=-1)
ax.axhline(0.71, color="0.4", lw=0.8, ls=":")
ax.text(1.05, 0.72, "simulated 100-150 Hz neurons: 0.71", fontsize=7, color="0.4", va="bottom")
YMAX = 1.6
n_hi = int((DISP.dispersion > YMAX).sum())
ax.text(0.99, 0.99, f"{n_hi} neuron(s) above {YMAX} not shown", transform=ax.transAxes,
        ha="right", va="top", fontsize=7, color="0.3")
ax.set(xscale="log", ylim=(0.3, YMAX), xlabel="firing rate in the analysis window (Hz)",
       ylabel="dispersion", title="B. dispersion vs firing rate (grey band = ~Poisson)")
ax.legend(fontsize=7, frameon=False, loc="lower left")

ax = axs[2]                                       # C: bands per region
tab = pd.crosstab(DISP.region, DISP.band, normalize="index").reindex(REG)
bottom = np.zeros(len(REG))
for band, c in zip(tab.columns, ("tab:blue", "0.75", "tab:red")):
    ax.bar(range(len(REG)), tab[band], bottom=bottom, color=c, width=0.6, label=band)
    for i, (b0, f) in enumerate(zip(bottom, tab[band])):
        if f > 0.04:
            ax.text(i, b0 + f / 2, f"{100*f:.0f}%", ha="center", va="center", fontsize=8)
    bottom += tab[band].to_numpy()
ax.set_xticks(range(len(REG)), [f"{r}\n(n={int((DISP.region == r).sum())})" for r in REG])
ax.set(ylabel="fraction of neurons", title="C. share per band", ylim=(0, 1.02))
ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
fig.suptitle(f"Poisson dispersion of real neurons in script 02's model (2.5 ms counts, own history "
             f"to 150 ms)  |  {len(SESS_LIST)} sessions", fontsize=9)
plotting.finish(fig, FIG_DIR, "dispersion", SAVE_FIGS, SHOW_FIGS)

# %%
