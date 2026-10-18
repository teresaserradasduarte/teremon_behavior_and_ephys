# %% ── 0. imports & constants ─────────────────────────────────────────────────
# Inspect the THREE neural-variable modes (see mi_rates.py): example traces per
# neuron and the population distribution of the windowed value.
# NOTHING is cached here. A smoothed rate trace is a property of the spike train
# and the smoothing params alone; the sliding window belongs to the MI step. The
# window-independent traces would cost 1-5 GB per mode to store, and re-deriving
# them costs ~2 min per mode, so 07 just derives them when it needs them.
#   "count" : spike count in the window            (integer, project default)
#   "isi"   : 1/ISI + Gaussian sigma 20 ms         (matches the MATLAB pipeline)
#   "conv"  : spike train * Gaussian sigma 70 ms   (window-averaged, Hz)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import importlib
import mi_io, config, mi_rates
importlib.reload(config); importlib.reload(mi_rates)

if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

CLR = {"CP": np.array([106,176,134])/255, "CB": np.array([154,123,193])/255}

# ── window parameters ─────────────────────────────────────────────────────────
# Most of this script is WINDOW-INDEPENDENT: the isi/conv example traces come
# from the smoothed rate itself, so changing the window below does not change
# them and there is no reason to re-run 06 when you re-window. The only parts
# that do depend on the window are the 'count' example column (a count is only
# defined given a window) and the cell-3 distributions. Match 04 so those two
# reflect what the MI step will see. Caching lives in 07, not here.
WIN_WIDTH  = 0.100        # sliding window width (s)
WIN_STEP   = 0.025        # step size (s)
WIN_RANGE  = (-1.0, 1.0)  # landscape range (s), re reach

# ── rate parameters (defaults come from mi_rates / the MATLAB pipeline) ───────
RATE_BIN   = mi_rates.RATE_BIN     # 2 ms fine grid
ISI_SIGMA  = mi_rates.ISI_SIGMA    # 20 ms  (MATLAB sig_pk)
CONV_SIGMA = mi_rates.CONV_SIGMA   # 70 ms

MODES      = ("count", "isi", "conv")
MODE_CLR   = {"count": "#444444", "isi": "#d95f02", "conv": "#1f78b4"}

N_PER_GROUP      = 3      # how many example neurons per density group, per region
                          # (groups = sparse / median / dense -> 3 x N_PER_GROUP x 2
                          # regions figures). They are taken as consecutive neurons
                          # centred on that group's quantile, so they sit in the same
                          # firing-rate band.
N_EXAMPLE_TRIALS = 3      # single trials drawn in the example figures
SHOW_FIGS  = True
SAVE_FIGS  = True

OUT_SUBDIR = "rate_modes"     # figures folder under OUT_DIR

win_starts  = np.arange(WIN_RANGE[0], WIN_RANGE[1] - WIN_WIDTH + WIN_STEP/2, WIN_STEP)
win_centers = win_starts + WIN_WIDTH / 2
n_wins      = len(win_starts)
print(f"{len(NEURONS)} neurons  |  windows: {n_wins}  width={WIN_WIDTH*1000:.0f} ms  "
      f"step={WIN_STEP*1000:.0f} ms  range={WIN_RANGE}")
print(f"rate grid {RATE_BIN*1000:.0f} ms  |  isi sigma {ISI_SIGMA*1000:.0f} ms  "
      f"|  conv sigma {CONV_SIGMA*1000:.0f} ms")


# %% ── 1. pick example neurons (sparse -> dense, per region) ──────────────────
# Mean count in the analysis range, used only to choose illustrative examples.
mean_ct = np.array([
    np.mean([((s >= WIN_RANGE[0]) & (s < WIN_RANGE[1])).sum() for s in n["spikes"]])
    if n["n_trials"] else 0.0 for n in NEURONS])

GROUPS = (("sparse", 0.10), ("median", 0.50), ("dense", 0.90))

def _pick(order, frac, k):
    """k neurons centred on the `frac` quantile of the rate-sorted order."""
    c  = int(np.clip(round(frac * (len(order) - 1)), 0, len(order) - 1))
    lo = int(np.clip(c - k // 2, 0, max(len(order) - k, 0)))
    return order[lo:lo + k]

examples = []
for reg in ("CP", "CB"):
    idx = np.array([i for i, n in enumerate(NEURONS)
                    if n["region"] == reg and n["n_trials"] >= 20])
    if idx.size == 0:
        continue
    order = idx[np.argsort(mean_ct[idx])]
    for gname, frac in GROUPS:
        for j, i in enumerate(_pick(order, frac, N_PER_GROUP)):
            tag = f"{reg} {gname}" + (f" #{j+1}" if N_PER_GROUP > 1 else "")
            examples.append((tag, int(i)))

print("\nexample neurons:")
for lab, i in examples:
    print(f"  {lab:14s}  {NEURONS[i]['uid']:32s}  mean count in range = {mean_ct[i]:.1f}")


# %% ── 2. example figures: raster + the three modes ───────────────────────────
ex_dir = config.fig_dir(OUT_SUBDIR, "examples")
rng_ex = np.random.default_rng(0)

for lab, i in examples:
    n = NEURONS[i]
    spikes = n["spikes"]
    n_tr = len(spikes)
    if n_tr == 0:
        continue
    show_tr = rng_ex.choice(n_tr, size=min(N_EXAMPLE_TRIALS, n_tr), replace=False)

    fig = plt.figure(figsize=(16, 5.6))
    gs  = gridspec.GridSpec(2, 4, width_ratios=[1.1, 1, 1, 1], figure=fig)

    # raster spanning both rows
    axr = fig.add_subplot(gs[:, 0])
    for ti, spk in enumerate(spikes):
        spk = np.asarray(spk)
        spk = spk[(spk >= WIN_RANGE[0]) & (spk < WIN_RANGE[1])]
        axr.scatter(spk, np.full(spk.size, ti), s=1.2,
                    c=[CLR[n["region"]]], linewidths=0)
    for ti in show_tr:
        axr.axhline(ti, color="k", lw=0.4, alpha=0.35)
    axr.axvline(0, color="k", lw=0.8, ls="--")
    axr.set(xlabel="time re reach (s)", ylabel="trial", title="raster",
            xlim=WIN_RANGE, ylim=(-1, n_tr))
    axr.spines["top"].set_visible(False); axr.spines["right"].set_visible(False)

    for mi_, mode in enumerate(MODES):
        if mode == "count":
            # plot what the estimator actually consumes: the integer count per
            # sliding window (raw 2 ms bins would just show 1/bin = 500 Hz spikes)
            t    = win_centers
            R    = np.asarray(mi_rates.windowed(spikes, win_starts, WIN_WIDTH,
                                                "count"), float)
            unit = f"spikes / {WIN_WIDTH*1000:.0f} ms win"
        else:
            t, R = mi_rates.rate_traces(spikes, mode, WIN_RANGE[0], WIN_RANGE[1],
                                        bin_width=RATE_BIN, isi_sigma=ISI_SIGMA,
                                        conv_sigma=CONV_SIGMA)
            unit = "Hz"

        # top row: mean +- SEM across trials
        ax = fig.add_subplot(gs[0, mi_ + 1])
        mn  = R.mean(axis=0)
        sem = R.std(axis=0) / np.sqrt(max(R.shape[0], 1))
        ax.plot(t, mn, color=MODE_CLR[mode], lw=1.6)
        ax.fill_between(t, mn - sem, mn + sem, color=MODE_CLR[mode], alpha=0.25)
        ax.axvline(0, color="k", lw=0.8, ls="--")
        ax.set(title=f"{mi_rates.MODE_LABEL[mode]}\nmean +- SEM", ylabel=unit,
               xlim=WIN_RANGE)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

        # bottom row: a few single trials
        ax = fig.add_subplot(gs[1, mi_ + 1])
        for k, ti in enumerate(show_tr):
            ax.plot(t, R[ti], lw=1.0, alpha=0.85, color=MODE_CLR[mode],
                    ls=["-", "--", ":"][k % 3], label=f"trial {ti}")
        ax.axvline(0, color="k", lw=0.8, ls="--")
        ax.set(xlabel="time re reach (s)", ylabel=unit,
               title="single trials", xlim=WIN_RANGE)
        ax.legend(fontsize=6, loc="upper right")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    fig.suptitle(f"{lab}  |  {n['uid']}  ({n_tr} trials)", fontsize=11)
    plt.tight_layout()
    if SAVE_FIGS:
        fig.savefig(ex_dir / f"{n['uid'].replace('|','_')}.png", dpi=130,
                    bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

print(f"example figures -> {ex_dir}")


# %% ── 3. windowed-variable distributions across the population ──────────────
# Sanity check that the three modes live on sensible, comparable scales.
sub = rng_ex.choice(len(NEURONS), size=min(300, len(NEURONS)), replace=False)
fig, axes = plt.subplots(1, 3, figsize=(15, 3.6))
for ax, mode in zip(axes, MODES):
    for reg in ("CP", "CB"):
        vals = []
        for i in sub:
            if NEURONS[i]["region"] != reg or NEURONS[i]["n_trials"] < 20:
                continue
            W = mi_rates.windowed(NEURONS[i]["spikes"], win_starts, WIN_WIDTH, mode,
                                  bin_width=RATE_BIN, isi_sigma=ISI_SIGMA,
                                  conv_sigma=CONV_SIGMA)
            vals.append(np.asarray(W, float).ravel())
        if vals:
            v = np.concatenate(vals)
            ax.hist(v, bins=60, histtype="step", lw=1.6, density=True,
                    color=CLR[reg], label=reg)
    ax.set(title=mi_rates.MODE_LABEL[mode],
           xlabel=f"windowed value ({mi_rates.MODE_UNIT[mode]})", ylabel="density")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.suptitle(f"windowed neural variable, {len(sub)}-neuron sample "
             f"(win {WIN_WIDTH*1000:.0f} ms)", fontsize=11)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir(OUT_SUBDIR) / "windowed_value_distributions.png",
                dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

# %%
