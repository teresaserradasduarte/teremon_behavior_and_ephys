# %% ── 0. imports & constants ─────────────────────────────────────────────────
import numpy as np
import matplotlib.pyplot as plt
import pickle, time, importlib
import mi_io, config, gpu_mi
importlib.reload(config)

if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

CLR     = {"CP": np.array([106,176,134])/255, "CB": np.array([154,123,193])/255}
rng     = np.random.default_rng(0)

print(f"gpu_mi device: {gpu_mi.DEVICE}")

# ── analysis parameters — change here only ────────────────────────────────────
WIN_RANGE  = (-1.0, 1.0)  # (s) re reach onset — both axes of the 2D landscape
GRID_STEP  = 0.050        # (s) spacing of the start/stop grid
MIN_WIN    = 0.050        # (s) minimum window width (stop - start); skip narrower
BIN_SIZE   = 0.0025       # (s) spike-bin resolution for fast cumsum counting
K_KSG      = 3            # k for KSG and Ross estimators
N_BINS_KIN = 16           # kin bins for plug-in estimator
COORD      = 1            # kinematic coord: 0=x, 1=y med-lat, 2=z
FRAMES     = (147, 167)   # endpoint window, MATLAB 1-indexed
CAT_KEEP   = {1, 2}       # reach categories to include
MIN_TRIALS = 20           # min valid trials to compute MI

N_SHUF_2D   = 20          # shuffle-null size per (neuron, window) cell, for Excess_MI
CHUNK_SIZE  = 40          # windows batched together per GPU call (KSG/Ross)

N_TEST    = None             # neurons per region for test; set None for all
SHOW_FIGS = True          # set False for full run (only saves)
SAVE_FIGS = True

CMAP      = "viridis"     # colormap for 2D landscapes

MASK_NONSIG = True       # if True, cells where |Excess_MI| isn't SIG_Z null-SDs from
                          # zero become NaN (rendered white) in the plotted landscapes
SIG_Z       = 1.96        # significance threshold, in units of the per-cell null std
                          # (z = mi_excess / mi_null_std; ~1.96 -> two-sided p<0.05)
# %%
MEAN_VMIN = None          # colorbar lower bound for the MEAN landscape plots only.
                          # None -> 0 (current default)
MEAN_VMAX = None          # colorbar upper bound for the MEAN landscape plots only.
                          # None -> auto (98th percentile of positive values, shared
                          # across CP/CB so the two panels stay comparable). Set an
                          # explicit float (e.g. 0.05) to fix the contrast yourself.


# %% ── 1. select neurons + precompute kin / mask / z-scored kin ───────────────
uids, regions, kin_list, kin_z_list, mask_list = [], [], [], [], []
for n in NEURONS:
    kin  = mi_io.endpoint_kin(n, REACH_PX, coord=COORD, frames=FRAMES)
    cat  = n["labels"]["idx_reach_cat"]
    mask = np.isfinite(kin) & np.isin(cat, list(CAT_KEEP))
    uids.append(n["uid"]); regions.append(n["region"])
    kin_m = kin[mask]
    kin_list.append(kin)
    kin_z_list.append(gpu_mi.zscore(kin_m) if kin_m.size > 1 else kin_m)  # z-scored once per neuron
    mask_list.append(mask)

cp_idx = [i for i, r in enumerate(regions) if r == "CP"]
cb_idx = [i for i, r in enumerate(regions) if r == "CB"]
if N_TEST is not None:
    run_idx = cp_idx[:N_TEST] + cb_idx[:N_TEST]
    print(f"TEST MODE: {len(run_idx)} neurons ({N_TEST} CP + {N_TEST} CB)")
else:
    run_idx = list(range(len(NEURONS)))
    print(f"FULL RUN: {len(run_idx)} neurons")


# %% ── 2. define grid + precompute cumsum bin edges ───────────────────────────
grid   = np.arange(WIN_RANGE[0], WIN_RANGE[1] + GRID_STEP/2, GRID_STEP)
n_grid = len(grid)

# valid (i_start, j_stop) pairs: stop > start + MIN_WIN
valid_pairs = [(i, j) for i in range(n_grid) for j in range(n_grid)
               if grid[j] - grid[i] >= MIN_WIN]
n_pairs = len(valid_pairs)
print(f"Grid: {n_grid} × {n_grid}  |  valid (start,stop) pairs: {n_pairs}  |  "
      f"N_SHUF_2D={N_SHUF_2D}  chunk_size={CHUNK_SIZE}")

# bin edges for cumsum (extend slightly beyond WIN_RANGE to catch all spikes)
cum_edges = np.arange(WIN_RANGE[0] - BIN_SIZE, WIN_RANGE[1] + 2*BIN_SIZE, BIN_SIZE)
# grid time → bin index lookup (precomputed once)
grid_bins = np.searchsorted(cum_edges, grid)   # shape (n_grid,)


# %% ── per-neuron cumsum builder ──────────────────────────────────────────────
def make_cumcount(spikes, cum_edges):
    """Returns (n_trials, n_bins+1) cumulative spike count matrix."""
    n_bins = len(cum_edges) - 1
    C = np.zeros((len(spikes), n_bins + 1), dtype=np.int16)
    for ti, spk in enumerate(spikes):
        if len(spk) == 0: continue
        C[ti, 1:] = np.cumsum(np.histogram(spk, bins=cum_edges)[0])
    return C

def window_counts(cumC, mask, i_start, j_stop):
    """Spike count in (grid[i_start], grid[j_stop]) for all valid trials."""
    return (cumC[mask, grid_bins[j_stop]] - cumC[mask, grid_bins[i_start]]).astype(float)

def all_window_counts(cumC, mask, pairs):
    """(n_pairs, n_trials) raw spike counts, one row per valid (start,stop) pair."""
    return np.stack([window_counts(cumC, mask, i_s, j_s) for (i_s, j_s) in pairs])


# %% ── significance masking (non-significant cells -> NaN -> rendered white) ───
def sig_mask(mi_excess, mi_null_std, z_thresh=SIG_Z, apply=MASK_NONSIG):
    """
    Returns a copy of mi_excess with cells where |mi_excess / mi_null_std| < z_thresh
    set to NaN (so _imshow_2d renders them white — "not significantly far from the
    shuffle null"). Pass apply=False (or leave MASK_NONSIG=False) to return mi_excess
    unchanged, e.g. for the raw/unmasked landscape.
    """
    if not apply:
        return mi_excess
    z = mi_excess / np.maximum(mi_null_std, 1e-12)
    out = mi_excess.copy()
    out[np.abs(z) < z_thresh] = np.nan
    return out


# %% ── shared plotting helpers (plot Excess_MI as the primary landscape) ──────
def _imshow_2d(ax, mi_2d, title, vmax=None, vmin=0, cmap=CMAP, cbar_label="Excess MI (bits)"):
    """Plot one neuron's or mean 2D Excess_MI matrix. NaN cells (invalid grid region,
    and — when MASK_NONSIG is on — non-significant cells) render as white."""
    finite = mi_2d[np.isfinite(mi_2d)]
    vm = vmax if vmax is not None else (np.nanpercentile(finite[finite > 0], 98)
                                        if np.any(finite > 0) else 1)
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("white")
    im = ax.imshow(mi_2d, origin="lower", aspect="equal", cmap=cmap_obj,
                   extent=[grid[0], grid[-1], grid[0], grid[-1]],
                   vmin=vmin, vmax=vm)
    ax.axhline(0, color="w", lw=0.8, ls="--"); ax.axvline(0, color="w", lw=0.8, ls="--")
    ax.set(xlabel="stop time (s)", ylabel="start time (s)", title=title)
    return im

def _plot_mean_2d(mi_stack, regions_run, method_label, save_dir, vmin=MEAN_VMIN, vmax=MEAN_VMAX):
    """Mean 2D Excess_MI landscape per region. vmin/vmax fix the contrast explicitly;
    leave as None (MEAN_VMIN/MEAN_VMAX) for the auto 98th-percentile behaviour, shared
    across the CP/CB panels so the two stay comparable."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    finite = mi_stack[np.isfinite(mi_stack)]
    vmn = vmin if vmin is not None else 0
    vmx = vmax if vmax is not None else (np.nanpercentile(finite[finite > 0], 98)
                                         if np.any(finite > 0) else 1)
    for ax, reg in zip(axes, ["CP", "CB"]):
        idx = [i for i, r in enumerate(regions_run) if r == reg]
        if not idx:
            ax.set_visible(False); continue
        mean_2d = np.nanmean(mi_stack[idx], axis=0)
        im = _imshow_2d(ax, mean_2d, f"{reg}  (n={len(idx)})\n{method_label}", vmin=vmn, vmax=vmx)
        plt.colorbar(im, ax=ax, label="Excess MI (bits)")
    plt.suptitle(f"Mean 2D Excess_MI landscape  |  {method_label}  (vmin={vmn:.3g}, vmax={vmx:.3g})",
                 fontsize=11)
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(save_dir / "mean_landscape_2D.png", dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

def _plot_lead_lag_2d(mi_stack, regions_run, method_label, save_dir):
    """Distribution of peak-window center time per region (peak of Excess_MI)."""
    fig, ax = plt.subplots(figsize=(6, 4))
    for reg in ["CP", "CB"]:
        idx = [i for i, r in enumerate(regions_run) if r == reg]
        peaks = []
        for pi in idx:
            m = mi_stack[pi]
            if np.any(np.isfinite(m) & (m > 0)):
                i_s, j_s = np.unravel_index(np.nanargmax(m), m.shape)
                peaks.append((grid[i_s] + grid[j_s]) / 2)   # window centre
        if peaks:
            ax.hist(peaks, bins=20, color=CLR[reg], alpha=0.6, density=False,
                    label=f"{reg} (n={len(peaks)})")
    ax.axvline(0, color="k", lw=0.8, ls=":", label="reach onset")
    ax.set(xlabel="centre of best window re reach (s)", ylabel="# neurons",
           title=f"{method_label}  —  lead / lag (peak Excess_MI)")
    ax.legend()
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(save_dir / "lead_lag_dist.png", dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

def _save_individuals(mi_stack, run_idx, regions_run, method_label, save_dir):
    """One Excess_MI heatmap per neuron → individual subfolder."""
    for pi, i in enumerate(run_idx):
        fig, ax = plt.subplots(figsize=(5, 4.5))
        row = mi_stack[pi]
        finite = row[np.isfinite(row)]
        vmax = np.nanpercentile(finite[finite > 0], 98) if np.any(finite > 0) else 1
        im = _imshow_2d(ax, row, f"{method_label}\n{uids[i]}", vmax=vmax)
        plt.colorbar(im, ax=ax, label="Excess MI (bits)")
        plt.tight_layout()
        if SAVE_FIGS:
            fig.savefig(save_dir / f"{uids[i].replace('|','_')}.png",
                        dpi=100, bbox_inches="tight")
        if SHOW_FIGS: plt.show()
        plt.close(fig)


"""
 ________  ___  ________   ________   ___  ________   ________
|\   __  \|\  \|\   ___  \|\   ___  \|\  \|\   ___  \|\   ____\
\ \  \|\ /\ \  \ \  \\ \  \ \  \\ \  \ \  \ \  \\ \  \ \  \___|
 \ \   __  \ \  \ \  \\ \  \ \  \\ \  \ \  \ \  \\ \  \ \  \  ___
  \ \  \|\  \ \  \ \  \\ \  \ \  \\ \  \ \  \ \  \\ \  \ \  \|\  \
   \ \_______\ \__\ \__\\ \__\ \__\\ \__\ \__\ \__\\ \__\ \_______\
    \|_______|\|__|\|__| \|__|\|__| \|__|\|__|\|__| \|__|\|_______|

""" 

# %% ── METHOD 1: Binning (plug-in) — CPU, no z-scoring (scale-invariant), shuffle-null ──
# Neural variable: spike COUNT (not FR as in MATLAB). Estimator: plug-in.
# Purpose: rough comparison with old MATLAB results. Known positive bias — Excess_MI
# (shuffle-null-subtracted) still helps here since it removes the entropy-ceiling
# confound even for a biased estimator.
print("\n=== METHOD: binning (plug-in) ===")
shape = (len(run_idx), n_grid, n_grid)
mi_bin_obs   = np.full(shape, np.nan, dtype=np.float32)
mi_bin_exc   = np.full(shape, np.nan, dtype=np.float32)
mi_bin_nmean = np.full(shape, np.nan, dtype=np.float32)
mi_bin_nstd  = np.full(shape, np.nan, dtype=np.float32)
regions_run = [regions[i] for i in run_idx]
t0 = time.time()

for pi, i in enumerate(run_idx):
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_m = kin_list[i][mask]           # raw kin, Binning uses quantile bins (scale-invariant)
    if mask.sum() < MIN_TRIALS: continue

    cumC = make_cumcount(n["spikes"], cum_edges)

    for (i_s, j_s) in valid_pairs:
        c_w = window_counts(cumC, mask, i_s, j_s)
        if np.unique(c_w).size < 2: continue
        mi_obs, mi_exc, nmean, nstd = gpu_mi.excess_mi(
            gpu_mi.bin_mi_batch, kin_m, c_w, N_SHUF_2D, rng, n_bins=N_BINS_KIN)
        mi_bin_obs[pi, i_s, j_s]   = mi_obs
        mi_bin_exc[pi, i_s, j_s]   = mi_exc
        mi_bin_nmean[pi, i_s, j_s] = nmean
        mi_bin_nstd[pi, i_s, j_s]  = nstd

    if (pi+1) % max(1, len(run_idx)//10) == 0:
        el = time.time()-t0; rate = (pi+1)/el
        print(f"  {pi+1}/{len(run_idx)}  {el:.0f}s  ~{(len(run_idx)-pi-1)/rate:.0f}s left",
              flush=True)

print(f"Done in {time.time()-t0:.1f}s")

res_bin = dict(method="bin_2D", grid=grid, WIN_RANGE=WIN_RANGE, GRID_STEP=GRID_STEP,
               mi=mi_bin_obs, mi_obs=mi_bin_obs, mi_excess=mi_bin_exc,
               mi_null_mean=mi_bin_nmean, mi_null_std=mi_bin_nstd, N_SHUF_2D=N_SHUF_2D,
               uids=[uids[i] for i in run_idx], regions=regions_run,
               COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP))
with open(config.var_path("landscape2D_bin.pkl"), "wb") as f:
    pickle.dump(res_bin, f, protocol=4)
print(f"Saved → {config.var_path('landscape2D_bin.pkl')}")


# %% ── METHOD 1: Binning — plots (Excess_MI) ─────────────────────────────────
mi_bin_plot = sig_mask(mi_bin_exc, mi_bin_nstd)
_save_individuals(mi_bin_plot, run_idx, regions_run, "Binning", config.fig_dir("landscapes","bin_2D","landscapes_individual"))

# %%
_plot_mean_2d(mi_bin_plot, regions_run, "Binning", config.fig_dir("landscapes","bin_2D"), vmin=0, vmax=0.12)
_plot_lead_lag_2d(mi_bin_plot, regions_run, "Binning", config.fig_dir("landscapes","bin_2D"))
print("Binning plots done.")


"""
 ___  __    ________  ________
|\  \|\  \ |\   ____\|\   ____\
\ \  \/  /|\ \  \___|\ \  \___|
 \ \   ___  \ \_____  \ \  \  ___
  \ \  \\ \  \|____|\  \ \  \|\  \
   \ \__\\ \__\____\_\  \ \_______\
    \|__| \|__|\_________\|_______|
              \|_________|
"""
# %% ── METHOD 2: KSG k=3 — GPU, z-scored both axes, chunked Excess_MI ─────────
print(f"\n=== METHOD: KSG k={K_KSG} (GPU, z-scored, Excess_MI) ===")
mi_ksg_obs   = np.full(shape, np.nan, dtype=np.float32)
mi_ksg_exc   = np.full(shape, np.nan, dtype=np.float32)
mi_ksg_nmean = np.full(shape, np.nan, dtype=np.float32)
mi_ksg_nstd  = np.full(shape, np.nan, dtype=np.float32)
t0 = time.time()

for pi, i in enumerate(run_idx):
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_z = kin_z_list[i]                       # z-scored once per neuron
    if mask.sum() < MIN_TRIALS: continue

    cumC = make_cumcount(n["spikes"], cum_edges)
    count_all_raw = all_window_counts(cumC, mask, valid_pairs)    # (n_pairs, n_trials)
    valid_win = np.array([np.unique(row).size >= 2 for row in count_all_raw])
    if not np.any(valid_win): continue

    count_all_z = gpu_mi.zscore(count_all_raw[valid_win], axis=1)  # z-score PER window
    kin_stack = gpu_mi.build_kin_stack(kin_z, N_SHUF_2D, rng)      # (1+N_SHUF, n_trials), same for all windows

    mi_obs, mi_exc, nmean, nstd = gpu_mi.excess_mi_windows(
        gpu_mi.ksg_mi_batch, kin_stack, count_all_z, k=K_KSG, chunk_size=CHUNK_SIZE)

    valid_idx = np.where(valid_win)[0]
    for local_k, pair_idx in enumerate(valid_idx):
        i_s, j_s = valid_pairs[pair_idx]
        mi_ksg_obs[pi, i_s, j_s]   = mi_obs[local_k]
        mi_ksg_exc[pi, i_s, j_s]   = mi_exc[local_k]
        mi_ksg_nmean[pi, i_s, j_s] = nmean[local_k]
        mi_ksg_nstd[pi, i_s, j_s]  = nstd[local_k]

    if (pi+1) % max(1, len(run_idx)//10) == 0:
        el = time.time()-t0; rate = (pi+1)/el
        print(f"  {pi+1}/{len(run_idx)}  {el:.0f}s  ~{(len(run_idx)-pi-1)/rate:.0f}s left",
              flush=True)

print(f"Done in {time.time()-t0:.1f}s")

res_ksg = dict(method="ksg_2D", grid=grid, WIN_RANGE=WIN_RANGE, GRID_STEP=GRID_STEP,
               mi=mi_ksg_obs, mi_obs=mi_ksg_obs, mi_excess=mi_ksg_exc,
               mi_null_mean=mi_ksg_nmean, mi_null_std=mi_ksg_nstd, N_SHUF_2D=N_SHUF_2D,
               uids=[uids[i] for i in run_idx], regions=regions_run,
               COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), K=K_KSG)
with open(config.var_path("landscape2D_ksg.pkl"), "wb") as f:
    pickle.dump(res_ksg, f, protocol=4)
print(f"Saved → {config.var_path('landscape2D_ksg.pkl')}")


# %% ── METHOD 2: KSG — plots (Excess_MI) ─────────────────────────────────────
mi_ksg_plot = sig_mask(mi_ksg_exc, mi_ksg_nstd)
_save_individuals(mi_ksg_plot, run_idx, regions_run, f"KSG k={K_KSG}", config.fig_dir("landscapes","ksg_2D","landscapes_individual"))

# %%
_plot_mean_2d(mi_ksg_exc, regions_run, f"KSG k={K_KSG}", config.fig_dir("landscapes","ksg_2D"), vmin=0, vmax=0.1)
_plot_lead_lag_2d(mi_ksg_plot, regions_run, f"KSG k={K_KSG}", config.fig_dir("landscapes","ksg_2D"))
print("KSG plots done.")


"""
 ________  ________  ________   ________
|\   __  \|\   __  \|\   ____\ |\   ____\
\ \  \|\  \ \  \|\  \ \  \___|_\ \  \___|_
 \ \   _  _\ \  \\\  \ \_____  \\ \_____  \
  \ \  \\  \\ \  \\\  \|____|\  \\|____|\  \
   \ \__\\ _\\ \_______\____\_\  \ ____\_\  \
    \|__|\|__|\|_______|\_________\\_________\
                       \|_________\|_________|
"""
# %% ── METHOD 3: Ross — GPU, kin z-scored (count raw int), CP only, chunked Excess_MI ──
print(f"\n=== METHOD: Ross k={K_KSG} (CP only, GPU, Excess_MI) ===")
cp_run = [(pi, run_idx[pi]) for pi, i in enumerate(run_idx) if regions[i] == "CP"]
mi_ross_obs   = np.full(shape, np.nan, dtype=np.float32)
mi_ross_exc   = np.full(shape, np.nan, dtype=np.float32)
mi_ross_nmean = np.full(shape, np.nan, dtype=np.float32)
mi_ross_nstd  = np.full(shape, np.nan, dtype=np.float32)
t0 = time.time()

for count_cp, (pi, i) in enumerate(cp_run):
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_z = kin_z_list[i]                       # z-scored kin; count stays raw int (class label)
    if mask.sum() < MIN_TRIALS: continue

    cumC = make_cumcount(n["spikes"], cum_edges)
    count_all_raw = all_window_counts(cumC, mask, valid_pairs)    # (n_pairs, n_trials)
    valid_win = np.array([np.unique(row).size >= 2 for row in count_all_raw])
    if not np.any(valid_win): continue

    count_all_int = count_all_raw[valid_win].astype(int)          # raw int class labels, NOT z-scored
    kin_stack = gpu_mi.build_kin_stack(kin_z, N_SHUF_2D, rng)

    mi_obs, mi_exc, nmean, nstd = gpu_mi.excess_mi_windows(
        gpu_mi.ross_mi_batch, kin_stack, count_all_int, k=K_KSG, chunk_size=CHUNK_SIZE)

    valid_idx = np.where(valid_win)[0]
    for local_k, pair_idx in enumerate(valid_idx):
        i_s, j_s = valid_pairs[pair_idx]
        mi_ross_obs[pi, i_s, j_s]   = mi_obs[local_k]
        mi_ross_exc[pi, i_s, j_s]   = mi_exc[local_k]
        mi_ross_nmean[pi, i_s, j_s] = nmean[local_k]
        mi_ross_nstd[pi, i_s, j_s]  = nstd[local_k]

    if (count_cp+1) % max(1, len(cp_run)//10) == 0:
        el = time.time()-t0; rate = (count_cp+1)/el
        print(f"  {count_cp+1}/{len(cp_run)} CP  {el:.0f}s  ~{(len(cp_run)-count_cp-1)/rate:.0f}s left",
              flush=True)

print(f"Done in {time.time()-t0:.1f}s")

res_ross = dict(method="ross_2D", grid=grid, WIN_RANGE=WIN_RANGE, GRID_STEP=GRID_STEP,
                mi=mi_ross_obs, mi_obs=mi_ross_obs, mi_excess=mi_ross_exc,
                mi_null_mean=mi_ross_nmean, mi_null_std=mi_ross_nstd, N_SHUF_2D=N_SHUF_2D,
                uids=[uids[i] for i in run_idx], regions=regions_run,
                COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), K=K_KSG)
with open(config.var_path("landscape2D_ross.pkl"), "wb") as f:
    pickle.dump(res_ross, f, protocol=4)
print(f"Saved → {config.var_path('landscape2D_ross.pkl')}")


# %% ── METHOD 3: Ross — plots + KSG vs Ross Excess_MI difference (CP) ────────
cp_pi = [pi for pi, i in cp_run]
mi_ross_plot = sig_mask(mi_ross_exc, mi_ross_nstd)
#_save_individuals(mi_ross_plot[cp_pi], [run_idx[pi] for pi in cp_pi],
#                  ["CP"]*len(cp_pi), f"Ross k={K_KSG}",
#                  config.fig_dir("landscapes","ross_2D_cp","landscapes_individual"))
# %%
_plot_mean_2d(mi_ross_plot[cp_pi], ["CP"]*len(cp_pi), f"Ross k={K_KSG} (CP)", config.fig_dir("landscapes","ross_2D_cp"), vmin=0, vmax=0.15)
_plot_lead_lag_2d(mi_ross_plot[cp_pi], ["CP"]*len(cp_pi), f"Ross k={K_KSG} (CP)", config.fig_dir("landscapes","ross_2D_cp"))

# KSG - Ross Excess_MI difference for CP (uses the same significance-masked arrays)
diff = mi_ksg_plot[cp_pi] - mi_ross_plot[cp_pi]
mean_diff = np.nanmean(diff, axis=0)
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
finite_diff = mean_diff[np.isfinite(mean_diff)]
vmax = np.nanpercentile(np.abs(finite_diff), 98) if finite_diff.size else 1
for ax, (data, title) in zip(axes, [
        (np.nanmean(mi_ksg_plot[cp_pi], axis=0), f"KSG k={K_KSG} (CP mean Excess_MI)"),
        (np.nanmean(mi_ross_plot[cp_pi], axis=0), f"Ross k={K_KSG} (CP mean Excess_MI)"),
        (mean_diff, "KSG − Ross (CP mean Excess_MI)")]):
    if "−" in title:
        im = ax.imshow(data, origin="lower", aspect="equal", cmap="RdBu_r",
                       extent=[grid[0],grid[-1],grid[0],grid[-1]],
                       vmin=-vmax, vmax=vmax)
    else:
        im = _imshow_2d(ax, data, title)
    ax.set_title(title); plt.colorbar(im, ax=ax, label="Excess MI (bits)")
plt.suptitle(f"KSG vs Ross Excess_MI cross-check (CP, n={len(cp_pi)})", fontsize=11)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir("landscapes","ross_2D_cp") / "diff_ksg_ross_CP.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)
print("Ross plots done.")

# %%
