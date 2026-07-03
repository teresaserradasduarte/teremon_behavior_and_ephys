# %% ── 0. imports & constants ─────────────────────────────────────────────────
import numpy as np
import matplotlib.pyplot as plt
import pickle, time, importlib
import mi_io, config
importlib.reload(config)
from sklearn.feature_selection import mutual_info_regression, mutual_info_classif

if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

CLR     = {"CP": np.array([106,176,134])/255, "CB": np.array([154,123,193])/255}
NAT2BIT = 1 / np.log(2)
rng     = np.random.default_rng(0)

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

N_TEST    = None             # neurons per region for test; set None for all
SHOW_FIGS = True          # set False for full run (only saves)
SAVE_FIGS = True

CMAP      = "viridis"     # colormap for 2D landscapes


# %% ── 1. select neurons + precompute kin / mask ──────────────────────────────
uids, regions, kin_list, mask_list = [], [], [], []
for n in NEURONS:
    kin  = mi_io.endpoint_kin(n, REACH_PX, coord=COORD, frames=FRAMES)
    cat  = n["labels"]["idx_reach_cat"]
    mask = np.isfinite(kin) & np.isin(cat, list(CAT_KEEP))
    uids.append(n["uid"]); regions.append(n["region"])
    kin_list.append(kin); mask_list.append(mask)

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
print(f"Grid: {n_grid} × {n_grid}  |  valid (start,stop) pairs: {n_pairs}")

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


# %% ── shared plotting helpers ────────────────────────────────────────────────
def _imshow_2d(ax, mi_2d, title, vmax=None):
    """Plot one neuron's or mean 2D MI matrix."""
    vm = vmax or np.nanpercentile(mi_2d[mi_2d > 0], 98) if np.any(mi_2d > 0) else 1
    im = ax.imshow(mi_2d, origin="lower", aspect="equal", cmap=CMAP,
                   extent=[grid[0], grid[-1], grid[0], grid[-1]],
                   vmin=0, vmax=vm)
    ax.axhline(0, color="w", lw=0.8, ls="--"); ax.axvline(0, color="w", lw=0.8, ls="--")
    ax.set(xlabel="stop time (s)", ylabel="start time (s)", title=title)
    return im

def _plot_mean_2d(mi_stack, regions_run, method_label, save_dir):
    """Mean 2D landscape per region."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    vmax = np.nanpercentile(mi_stack[mi_stack > 0], 98) if np.any(mi_stack > 0) else 1
    for ax, reg in zip(axes, ["CP", "CB"]):
        idx = [i for i, r in enumerate(regions_run) if r == reg]
        if not idx:
            ax.set_visible(False); continue
        mean_2d = np.nanmean(mi_stack[idx], axis=0)
        im = _imshow_2d(ax, mean_2d, f"{reg}  (n={len(idx)})\n{method_label}", vmax=vmax)
        plt.colorbar(im, ax=ax, label="MI (bits)")
    plt.suptitle(f"Mean 2D MI landscape  |  {method_label}", fontsize=11)
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(save_dir / "mean_landscape_2D.png", dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

def _plot_lead_lag_2d(mi_stack, regions_run, method_label, save_dir):
    """Distribution of peak-window center time per region."""
    fig, ax = plt.subplots(figsize=(6, 4))
    for reg in ["CP", "CB"]:
        idx = [i for i, r in enumerate(regions_run) if r == reg]
        peaks = []
        for pi in idx:
            m = mi_stack[pi]
            if np.any(m > 0):
                i_s, j_s = np.unravel_index(np.nanargmax(m), m.shape)
                peaks.append((grid[i_s] + grid[j_s]) / 2)   # window centre
        if peaks:
            ax.hist(peaks, bins=20, color=CLR[reg], alpha=0.6, density=False,
                    label=f"{reg} (n={len(peaks)})")
    ax.axvline(0, color="k", lw=0.8, ls=":", label="reach onset")
    ax.set(xlabel="centre of best window re reach (s)", ylabel="# neurons",
           title=f"{method_label}  —  lead / lag")
    ax.legend()
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(save_dir / "lead_lag_dist.png", dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

def _save_individuals(mi_stack, run_idx, regions_run, method_label, save_dir):
    """One heatmap per neuron → individual subfolder."""
    for pi, i in enumerate(run_idx):
        fig, ax = plt.subplots(figsize=(5, 4.5))
        vmax = np.nanpercentile(mi_stack[pi][mi_stack[pi] > 0], 98) \
               if np.any(mi_stack[pi] > 0) else 1
        im = _imshow_2d(ax, mi_stack[pi],
                        f"{method_label}\n{uids[i]}", vmax=vmax)
        plt.colorbar(im, ax=ax, label="MI (bits)")
        plt.tight_layout()
        if SAVE_FIGS:
            fig.savefig(save_dir / f"{uids[i].replace('|','_')}.png",
                        dpi=100, bbox_inches="tight")
        if SHOW_FIGS: plt.show()
        plt.close(fig)


# %% ── METHOD 1: Binning (plug-in) — compute 2D landscape ────────────────────
# Neural variable: spike COUNT (not FR as in MATLAB). Estimator: plug-in.
# Purpose: rough comparison with old MATLAB results. Known positive bias.
print("\n=== METHOD: binning (plug-in) ===")
mi_bin = np.full((len(run_idx), n_grid, n_grid), np.nan, dtype=np.float32)
regions_run = [regions[i] for i in run_idx]
t0 = time.time()

for pi, i in enumerate(run_idx):
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_m = kin_list[i][mask]
    if mask.sum() < MIN_TRIALS: continue

    cumC = make_cumcount(n["spikes"], cum_edges)

    for (i_s, j_s) in valid_pairs:
        c_w = window_counts(cumC, mask, i_s, j_s)
        if np.unique(c_w).size < 2: continue
        edges = np.unique(np.quantile(kin_m, np.linspace(0, 1, N_BINS_KIN + 1)))
        kb = np.clip(np.digitize(kin_m, edges[1:-1]), 0, len(edges) - 2)
        cv = np.unique(c_w); cm_ = {v: k for k, v in enumerate(cv)}
        ci = np.array([cm_[v] for v in c_w])
        J  = np.zeros((len(cv), len(edges)-1)); np.add.at(J, (ci, kb), 1.); P = J / J.sum()
        pc = P.sum(1, keepdims=True); pk_ = P.sum(0, keepdims=True); nz = P > 0
        mi_bin[pi, i_s, j_s] = float(np.sum(P[nz]*np.log2(P[nz]/((pc@pk_+1e-300)[nz]))))

    if (pi+1) % max(1, len(run_idx)//10) == 0:
        el = time.time()-t0; rate = (pi+1)/el
        print(f"  {pi+1}/{len(run_idx)}  {el:.0f}s  ~{(len(run_idx)-pi-1)/rate:.0f}s left",
              flush=True)

print(f"Done in {time.time()-t0:.1f}s")

res_bin = dict(method="bin_2D", grid=grid, WIN_RANGE=WIN_RANGE, GRID_STEP=GRID_STEP,
               mi=mi_bin, uids=[uids[i] for i in run_idx], regions=regions_run,
               COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP))
with open(config.var_path("landscape2D_bin.pkl"), "wb") as f:
    pickle.dump(res_bin, f, protocol=4)
print(f"Saved → {config.var_path('landscape2D_bin.pkl')}")


# %% ── METHOD 1: Binning — plots ─────────────────────────────────────────────
_save_individuals(mi_bin, run_idx, regions_run, "Binning", config.fig_dir("landscapes","bin_2D","landscapes_individual"))
_plot_mean_2d(mi_bin, regions_run, "Binning", config.fig_dir("landscapes","bin_2D"))
_plot_lead_lag_2d(mi_bin, regions_run, "Binning", config.fig_dir("landscapes","bin_2D"))
print("Binning plots done.")


# %% ── METHOD 2: KSG k=3 — compute 2D landscape ──────────────────────────────
print(f"\n=== METHOD: KSG k={K_KSG} ===")
mi_ksg = np.full((len(run_idx), n_grid, n_grid), np.nan, dtype=np.float32)
t0 = time.time()

for pi, i in enumerate(run_idx):
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_m = kin_list[i][mask]
    if mask.sum() < MIN_TRIALS: continue

    cumC = make_cumcount(n["spikes"], cum_edges)
    X    = kin_m.reshape(-1, 1)

    for (i_s, j_s) in valid_pairs:
        c_w = window_counts(cumC, mask, i_s, j_s)
        if np.unique(c_w).size < 2: continue
        mi_ksg[pi, i_s, j_s] = mutual_info_regression(
            X, c_w + rng.normal(0, 1e-9, c_w.size),
            n_neighbors=K_KSG, random_state=0)[0] * NAT2BIT

    if (pi+1) % max(1, len(run_idx)//10) == 0:
        el = time.time()-t0; rate = (pi+1)/el
        print(f"  {pi+1}/{len(run_idx)}  {el:.0f}s  ~{(len(run_idx)-pi-1)/rate:.0f}s left",
              flush=True)

print(f"Done in {time.time()-t0:.1f}s")

res_ksg = dict(method="ksg_2D", grid=grid, WIN_RANGE=WIN_RANGE, GRID_STEP=GRID_STEP,
               mi=mi_ksg, uids=[uids[i] for i in run_idx], regions=regions_run,
               COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), K=K_KSG)
with open(config.var_path("landscape2D_ksg.pkl"), "wb") as f:
    pickle.dump(res_ksg, f, protocol=4)
print(f"Saved → {config.var_path('landscape2D_ksg.pkl')}")


# %% ── METHOD 2: KSG — plots ─────────────────────────────────────────────────
_save_individuals(mi_ksg, run_idx, regions_run, f"KSG k={K_KSG}", config.fig_dir("landscapes","ksg_2D","landscapes_individual"))
_plot_mean_2d(mi_ksg, regions_run, f"KSG k={K_KSG}", config.fig_dir("landscapes","ksg_2D"))
_plot_lead_lag_2d(mi_ksg, regions_run, f"KSG k={K_KSG}", config.fig_dir("landscapes","ksg_2D"))
print("KSG plots done.")


# %% ── METHOD 3: Ross — CP only, cross-check ─────────────────────────────────
print(f"\n=== METHOD: Ross k={K_KSG} (CP only) ===")
cp_run = [(pi, run_idx[pi]) for pi, i in enumerate(run_idx) if regions[i] == "CP"]
mi_ross = np.full((len(run_idx), n_grid, n_grid), np.nan, dtype=np.float32)
t0 = time.time()

for count_cp, (pi, i) in enumerate(cp_run):
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_m = kin_list[i][mask]
    if mask.sum() < MIN_TRIALS: continue

    cumC = make_cumcount(n["spikes"], cum_edges)
    X    = kin_m.reshape(-1, 1)

    for (i_s, j_s) in valid_pairs:
        c_w = window_counts(cumC, mask, i_s, j_s)
        if np.unique(c_w).size < 2: continue
        mi_ross[pi, i_s, j_s] = mutual_info_classif(
            X, c_w.astype(int), n_neighbors=K_KSG, random_state=0)[0] * NAT2BIT

    if (count_cp+1) % max(1, len(cp_run)//10) == 0:
        el = time.time()-t0; rate = (count_cp+1)/el
        print(f"  {count_cp+1}/{len(cp_run)} CP  {el:.0f}s  ~{(len(cp_run)-count_cp-1)/rate:.0f}s left",
              flush=True)

print(f"Done in {time.time()-t0:.1f}s")

res_ross = dict(method="ross_2D", grid=grid, WIN_RANGE=WIN_RANGE, GRID_STEP=GRID_STEP,
                mi=mi_ross, uids=[uids[i] for i in run_idx], regions=regions_run,
                COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), K=K_KSG)
with open(config.var_path("landscape2D_ross.pkl"), "wb") as f:
    pickle.dump(res_ross, f, protocol=4)
print(f"Saved → {config.var_path('landscape2D_ross.pkl')}")


# %% ── METHOD 3: Ross — plots + KSG vs Ross difference (CP) ──────────────────
cp_pi = [pi for pi, i in cp_run]
_save_individuals(mi_ross[cp_pi], [run_idx[pi] for pi in cp_pi],
                  ["CP"]*len(cp_pi), f"Ross k={K_KSG}",
                  config.fig_dir("landscapes","ross_2D_cp","landscapes_individual"))
_plot_mean_2d(mi_ross[cp_pi], ["CP"]*len(cp_pi), f"Ross k={K_KSG} (CP)", config.fig_dir("landscapes","ross_2D_cp"))
_plot_lead_lag_2d(mi_ross[cp_pi], ["CP"]*len(cp_pi), f"Ross k={K_KSG} (CP)", config.fig_dir("landscapes","ross_2D_cp"))

# KSG - Ross difference for CP
diff = mi_ksg[cp_pi] - mi_ross[cp_pi]
mean_diff = np.nanmean(diff, axis=0)
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
vmax = np.nanpercentile(np.abs(mean_diff[np.isfinite(mean_diff)]), 98) if np.any(np.isfinite(mean_diff)) else 1
for ax, (data, title) in zip(axes, [
        (np.nanmean(mi_ksg[cp_pi], axis=0), f"KSG k={K_KSG} (CP mean)"),
        (np.nanmean(mi_ross[cp_pi], axis=0), f"Ross k={K_KSG} (CP mean)"),
        (mean_diff, "KSG − Ross (CP mean)")]):
    if "−" in title:
        im = ax.imshow(data, origin="lower", aspect="equal", cmap="RdBu_r",
                       extent=[grid[0],grid[-1],grid[0],grid[-1]],
                       vmin=-vmax, vmax=vmax)
    else:
        im = _imshow_2d(ax, data, title)
    ax.set_title(title); plt.colorbar(im, ax=ax, label="MI (bits)")
plt.suptitle(f"KSG vs Ross cross-check (CP, n={len(cp_pi)})", fontsize=11)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir("landscapes","ross_2D_cp") / "diff_ksg_ross_CP.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)
print("Ross plots done.")

# %%
