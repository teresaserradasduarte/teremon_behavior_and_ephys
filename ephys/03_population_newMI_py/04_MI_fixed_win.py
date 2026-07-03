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
WIN_WIDTH   = 0.100       # sliding window width (s)
WIN_STEP    = 0.025       # step size (s)  →  temporal resolution of the landscape
WIN_RANGE   = (-1.0, 1.0) # landscape range (s), re reach onset
K_KSG       = 3           # k for KSG estimator
N_BINS_KIN  = 16          # kin bins for plug-in estimator
COORD       = 1           # kinematic coord: 0=x, 1=y med-lat, 2=z
FRAMES      = (147, 167)  # endpoint window, MATLAB 1-indexed
CAT_KEEP    = {1, 2}      # reach categories to include
MIN_TRIALS  = 20          # min valid trials to attempt MI

N_TEST           = None      # neurons per region for test run; set None for all neurons
N_SHUF_LEADLAG   = 50     # shuffles per neuron for lead-lag significance (peak window only)
SMOOTH_WIN       = 3      # moving-average kernel (# windows) for peak detection on individual plots
PEAK_Z           = 3.0    # peak must exceed mean + PEAK_Z*sigma of the curve
PEAK_RANGE       = None   # e.g. (-0.3, 0.3) to restrict lead-lag peak search/median to that window; None = full WIN_RANGE
SHOW_FIGS        = True   # set False for full run (only saves, faster)
SAVE_FIGS        = True


# %% ── 1. select neurons to process ──────────────────────────────────────────
# endpoint kin + category mask — computed once and stored in lists
uids, regions, kin_list, mask_list = [], [], [], []
for n in NEURONS:
    kin  = mi_io.endpoint_kin(n, REACH_PX, coord=COORD, frames=FRAMES)
    cat  = n["labels"]["idx_reach_cat"]
    mask = np.isfinite(kin) & np.isin(cat, list(CAT_KEEP))
    uids.append(n["uid"]); regions.append(n["region"])
    kin_list.append(kin); mask_list.append(mask)

# subset for test run
cp_idx = [i for i, r in enumerate(regions) if r == "CP"]
cb_idx = [i for i, r in enumerate(regions) if r == "CB"]
if N_TEST is not None:
    run_idx = cp_idx[:N_TEST] + cb_idx[:N_TEST]
    print(f"TEST MODE: {len(run_idx)} neurons ({N_TEST} CP + {N_TEST} CB)")
else:
    run_idx = list(range(len(NEURONS)))
    print(f"FULL RUN: {len(run_idx)} neurons")


# %% ── 2. define sliding windows ─────────────────────────────────────────────
win_starts  = np.arange(WIN_RANGE[0], WIN_RANGE[1] - WIN_WIDTH + WIN_STEP/2, WIN_STEP)
win_centers = win_starts + WIN_WIDTH / 2
n_wins      = len(win_starts)
print(f"Windows: {n_wins} steps  |  width={WIN_WIDTH*1000:.0f} ms  "
      f"step={WIN_STEP*1000:.0f} ms  range={WIN_RANGE}")

# %% ── helpers ────────────────────────────────────────────────────────────────
"""
 ___  ___  _______   ___       ________  _______   ________  ________
|\  \|\  \|\  ___ \ |\  \     |\   __  \|\  ___ \ |\   __  \|\   ____\
\ \  \\\  \ \   __/|\ \  \    \ \  \|\  \ \   __/|\ \  \|\  \ \  \___|_
 \ \   __  \ \  \_|/_\ \  \    \ \   ____\ \  \_|/_\ \   _  _\ \_____  \
  \ \  \ \  \ \  \_|\ \ \  \____\ \  \___|\ \  \_|\ \ \  \\  \\|____|\  \
   \ \__\ \__\ \_______\ \_______\ \__\    \ \_______\ \__\\ _\ ____\_\  \
    \|__|\|__|\|_______|\|_______|\|__|     \|_______|\|__|\|__|\_________\
                                                               \|_________|


"""

def counts_landscape(spikes, win_starts, WIN_WIDTH):
    """Returns (n_trials, n_wins) int array via sorted searchsorted."""
    hi = win_starts + WIN_WIDTH
    C  = np.zeros((len(spikes), len(win_starts)), dtype=np.int16)
    for ti, spk in enumerate(spikes):
        if len(spk) == 0: continue
        s = np.sort(spk)
        C[ti] = np.searchsorted(s, hi) - np.searchsorted(s, win_starts)
    return C

def _ksg(kin_m, c_w, k):
    X = kin_m.reshape(-1, 1)
    return mutual_info_regression(X, c_w + rng.normal(0, 1e-9, c_w.size),
                                  n_neighbors=k, random_state=0)[0] * NAT2BIT

def _ross(kin_m, c_w, k):
    X = kin_m.reshape(-1, 1)
    return mutual_info_classif(X, c_w.astype(int),
                               n_neighbors=k, random_state=0)[0] * NAT2BIT

def _bin(kin_m, c_w, nb):
    edges = np.unique(np.quantile(kin_m, np.linspace(0, 1, nb + 1)))
    kb = np.clip(np.digitize(kin_m, edges[1:-1]), 0, len(edges) - 2)
    cv = np.unique(c_w); cm = {v: i for i, v in enumerate(cv)}
    ci = np.array([cm[v] for v in c_w])
    J  = np.zeros((len(cv), len(edges)-1)); np.add.at(J, (ci, kb), 1.0); P = J / J.sum()
    pc = P.sum(1, keepdims=True); pk = P.sum(0, keepdims=True); nz = P > 0
    return float(np.sum(P[nz] * np.log2(P[nz] / (pc @ pk + 1e-300)[nz])))

def _smooth(arr, win):
    """Centered moving average, NaN-aware, edge-safe (partial windows at boundaries)."""
    if win <= 1: return arr.copy()
    n = len(arr)
    out = np.full(n, np.nan)
    half = win // 2
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        seg = arr[lo:hi]
        if np.any(~np.isnan(seg)):
            out[i] = np.nanmean(seg)
    return out

def _plot_individual(uid, reg, mi_arr, win_centers, method_label, save_dir, sig=None, entropy=None):
    fig, ax = plt.subplots(figsize=(7, 2.8))
    ax.plot(win_centers, mi_arr, color=CLR[reg], lw=0.8, alpha=0.35, label="raw")

    mi_sm = _smooth(mi_arr, SMOOTH_WIN)
    ax.plot(win_centers, mi_sm, color=CLR[reg], lw=1.6, label=f"smoothed (win={SMOOTH_WIN})")

    peaked = False
    if np.any(~np.isnan(mi_sm)):
        pk_val = np.nanmax(mi_sm)
        pk_t   = win_centers[np.nanargmax(mi_sm)]
        mn_c   = np.nanmean(mi_sm)
        sd_c   = np.nanstd(mi_sm)
        thresh = mn_c + PEAK_Z * sd_c
        peaked = pk_val > thresh
        ax.axvline(pk_t, color=CLR[reg], lw=0.8, ls="--", alpha=0.7)
        ax.axhline(thresh, color="gray", lw=0.8, ls="--", alpha=0.7, label=f"μ+{PEAK_Z:.0f}σ (smoothed)")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.axhline(0, color="k", lw=0.4)
    sig_str  = ("  * sig" if sig else "  ns") if sig is not None else ""
    ent_str  = f"  H={entropy:.2f} bits" if (entropy is not None and np.isfinite(entropy)) else ""
    peak_str = f"  peak > μ+{PEAK_Z:.0f}σ" if peaked else f"  peak ≤ μ+{PEAK_Z:.0f}σ"
    ax.set(xlabel="window center reach (s)", ylabel="MI (bits)",
           title=f"{method_label}  |  {uid}{sig_str}{ent_str}{peak_str}")
    ax.legend(fontsize=7, loc="upper right")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    if SAVE_FIGS:
        fig.savefig(save_dir / f"{uid.replace('|','_')}.png", dpi=100, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

def _plot_mean_landscape(mi_mat, regions_run, win_centers, method_label, save_dir):
    """Mean ± SEM landscape per region. Adds a dashed vertical line at the peak of the mean MI curve."""
    fig, ax = plt.subplots(figsize=(8, 4))
    for reg in ["CP", "CB"]:
        idx = [i for i, r in enumerate(regions_run) if r == reg]
        if not idx: continue
        m   = mi_mat[idx]
        ok  = ~np.all(np.isnan(m), axis=1)
        m   = m[ok]
        mn  = np.nanmean(m, axis=0)
        sem = np.nanstd(m, axis=0) / np.sqrt(np.maximum(np.sum(~np.isnan(m), axis=0), 1))
        ax.plot(win_centers, mn, color=CLR[reg], lw=2, label=f"{reg} (n={ok.sum()})")
        ax.fill_between(win_centers, mn-sem, mn+sem, color=CLR[reg], alpha=0.2)

        # peak of the mean MI curve
        pk = win_centers[np.nanargmax(mn)]
        ax.axvline(pk, color=CLR[reg], lw=1.5, ls="--", alpha=0.85,
                   label=f"{reg} peak = {pk*1000:.0f} ms")

    ax.axvline(0, color="k", lw=0.8, ls=":", label="reach endpoint")
    ax.axhline(0, color="k", lw=0.4)
    ax.set(xlabel="window center reach (s)", ylabel="mean MI ± SEM (bits)",
           title=f"{method_label}  —  mean landscape")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(save_dir / "mean_landscape.png", dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

def _plot_lead_lag(mi_mat, regions_run, win_centers, method_label, save_dir, sig=None, peak_range=None):
    """Lead-lag histogram of peak MI time per neuron.
    sig : 1-D bool array (same length as mi_mat rows).
          Pass sig=None to include all neurons; pass a bool array to keep only significant ones.
    peak_range : optional (lo, hi) in seconds. Restricts BOTH the search window used to find
                 each neuron's peak AND the displayed/median x-range to [lo, hi].
                 Pass None to use the full landscape range (original behaviour)."""
    sig_label = (f"significant only  (p<0.05, n={N_SHUF_LEADLAG} shuffle)"
                 if sig is not None else "all neurons (no significance filter)")
    if peak_range is not None:
        lo, hi = peak_range
        win_mask = (win_centers >= lo) & (win_centers <= hi)
    else:
        lo, hi = win_centers[0], win_centers[-1]
        win_mask = np.ones(len(win_centers), dtype=bool)
    bin_edges = np.linspace(lo, hi, 21)   # fixed, shared edges — same for both regions, no auto-range surprises

    fig, ax = plt.subplots(figsize=(6, 4))
    for reg in ["CP", "CB"]:
        idx = [i for i, r in enumerate(regions_run) if r == reg]
        if not idx: continue
        m = mi_mat[idx]
        s = sig[idx] if sig is not None else np.ones(len(idx), dtype=bool)
        peaks = []
        for row, keep in zip(m, s):
            if not keep or not np.any(~np.isnan(row)): continue
            row_sm = _smooth(row, SMOOTH_WIN)
            row_sm_masked = np.where(win_mask, row_sm, np.nan)
            if np.any(~np.isnan(row_sm_masked)) and np.nanmax(row_sm_masked) > 0:
                peaks.append(win_centers[np.nanargmax(row_sm_masked)])
        if peaks:
            ax.hist(peaks, bins=bin_edges, color=CLR[reg], alpha=0.6, edgecolor="white",
                    label=f"{reg} (n={len(peaks)})", density=False)
            med = np.median(peaks)
            ax.axvline(med, color=CLR[reg], lw=1.5, ls="--", alpha=0.9,
                       label=f"{reg} median = {med*1000:.0f} ms")
    ax.axvline(0, color="k", lw=0.8, ls=":", label="reach endpoint")
    range_str = f"  |  peak search restricted to [{peak_range[0]*1000:.0f}, {peak_range[1]*1000:.0f}] ms" if peak_range else ""
    ax.set(xlabel="time of peak MI reach (s)", ylabel="# neurons",
           title=f"{method_label}  —  lead / lag\n{sig_label}{range_str}")
    if peak_range is not None:
        ax.set_xlim(peak_range)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend()
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(save_dir / "lead_lag_dist.png", dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


"""
 _____ ______   ___          ________  ___  ________
|\   _ \  _   \|\  \        |\   __  \|\  \|\   ___  \
\ \  \\\__\ \  \ \  \       \ \  \|\ /\ \  \ \  \\ \  \
 \ \  \\|__| \  \ \  \       \ \   __  \ \  \ \  \\ \  \
  \ \  \    \ \  \ \  \       \ \  \|\  \ \  \ \  \\ \  \
   \ \__\    \ \__\ \__\       \ \_______\ \__\ \__\\ \__\
    \|__|     \|__|\|__|        \|_______|\|__|\|__| \|__|

"""
# %% ── METHOD 1: BINNING (plug-in) — compute + save ──────────────────────────
# Purpose: rough comparison with old MATLAB results; known positive bias.
print("\n=== METHOD: binning (plug-in) ===")
mi_bin = np.full((len(run_idx), n_wins), np.nan)
regions_run = [regions[i] for i in run_idx]
t0 = time.time()

for pi, i in enumerate(run_idx):
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_m = kin_list[i][mask]
    if mask.sum() < MIN_TRIALS: continue

    C   = counts_landscape(n["spikes"], win_starts, WIN_WIDTH)
    C_m = C[mask]

    for wi in range(n_wins):
        c_w = C_m[:, wi]
        if np.unique(c_w).size < 2: continue
        mi_bin[pi, wi] = _bin(kin_m, c_w, N_BINS_KIN)

    if (pi + 1) % 5 == 0:
        elapsed = time.time() - t0
        rate = (pi+1) / elapsed
        print(f"  {pi+1}/{len(run_idx)}  ({elapsed:.0f}s elapsed, "
              f"~{(len(run_idx)-pi-1)/rate:.0f}s remaining)", flush=True)

print(f"Done in {time.time()-t0:.1f}s")

res_bin = dict(method="bin", win_centers=win_centers, win_starts=win_starts,
               WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP, mi=mi_bin,
               uids=[uids[i] for i in run_idx], regions=regions_run,
               COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP))
with open(config.var_path("landscape_bin.pkl"), "wb") as f:
    pickle.dump(res_bin, f, protocol=4)
print(f"Saved → {config.var_path('landscape_bin.pkl')}")


# %% ── METHOD 1: BINNING — lead-lag significance (shuffle at peak window) ─────
# Shuffle kin at each neuron's peak window; only peaks > p95 null enter lead-lag.
print(f"Lead-lag significance: {N_SHUF_LEADLAG} shuffles per neuron (binning)")
sig_bin = np.zeros(len(run_idx), dtype=bool)
t0 = time.time()

for pi, i in enumerate(run_idx):
    mi_arr = mi_bin[pi]
    if np.all(np.isnan(mi_arr)): continue
    mi_sm   = _smooth(mi_arr, SMOOTH_WIN)          # peak window from smoothed curve
    if np.all(np.isnan(mi_sm)): continue
    pk_wi   = int(np.nanargmax(mi_sm))
    real_mi = mi_arr[pk_wi]                         # significance uses the actual (raw) MI at that window
    if real_mi <= 0: continue

    n    = NEURONS[i]
    mask = mask_list[i]
    if mask.sum() < MIN_TRIALS: continue
    kin_m = kin_list[i][mask]

    c_pk = mi_io.counts(n, (float(win_starts[pk_wi]),
                             float(win_starts[pk_wi] + WIN_WIDTH)))[mask]
    if np.unique(c_pk).size < 2: continue

    null = np.array([_bin(rng.permutation(kin_m), c_pk, N_BINS_KIN)
                     for _ in range(N_SHUF_LEADLAG)])
    sig_bin[pi] = real_mi > np.percentile(null, 95)

print(f"  {sig_bin.sum()}/{len(run_idx)} significant  ({time.time()-t0:.1f}s)")


# %% ── METHOD 1: BINNING — plots individual landscapes ─────────────────────────
ind_dir = config.fig_dir("landscapes", "bin", "landscapes_individual")
for pi, i in enumerate(run_idx):
    _plot_individual(uids[i], regions[i], mi_bin[pi], win_centers, "Binning", ind_dir, sig=sig_bin[pi])
# %% ── METHOD 1: BINNING — plots mean / lead-lag ───────────────────────────────


_plot_mean_landscape(mi_bin, regions_run, win_centers, "Binning",
                     config.fig_dir("landscapes", "bin"))
_plot_lead_lag(mi_bin, regions_run, win_centers, "Binning",
               config.fig_dir("landscapes", "bin"), sig=sig_bin, peak_range=PEAK_RANGE)
print("Binning plots done.")


"""
 _____ ______   ___          ___  __    ________  ________
|\   _ \  _   \|\  \        |\  \|\  \ |\   ____\|\   ____\
\ \  \\\__\ \  \ \  \       \ \  \/  /|\ \  \___|\ \  \___|
 \ \  \\|__| \  \ \  \       \ \   ___  \ \_____  \ \  \  ___
  \ \  \    \ \  \ \  \       \ \  \\ \  \|____|\  \ \  \|\  \
   \ \__\    \ \__\ \__\       \ \__\\ \__\____\_\  \ \_______\
    \|__|     \|__|\|__|        \|__| \|__|\_________\|_______|
                                          \|_________|


"""
# %% ── METHOD 2: KSG k=3 — compute + save ────────────────────────────────────
print(f"\n=== METHOD: KSG  k={K_KSG} ===")
mi_ksg = np.full((len(run_idx), n_wins), np.nan)
t0 = time.time()

for pi, i in enumerate(run_idx):
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_m = kin_list[i][mask]
    if mask.sum() < MIN_TRIALS: continue

    C   = counts_landscape(n["spikes"], win_starts, WIN_WIDTH)
    C_m = C[mask]

    for wi in range(n_wins):
        c_w = C_m[:, wi]
        if np.unique(c_w).size < 2: continue
        mi_ksg[pi, wi] = _ksg(kin_m, c_w, K_KSG)

    if (pi + 1) % 5 == 0:
        elapsed = time.time() - t0
        rate = (pi+1) / elapsed
        print(f"  {pi+1}/{len(run_idx)}  ({elapsed:.0f}s elapsed, "
              f"~{(len(run_idx)-pi-1)/rate:.0f}s remaining)", flush=True)

print(f"Done in {time.time()-t0:.1f}s")

res_ksg = dict(method="ksg", win_centers=win_centers, win_starts=win_starts,
               WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP, mi=mi_ksg,
               uids=[uids[i] for i in run_idx], regions=regions_run,
               COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), K=K_KSG)
with open(config.var_path("landscape_ksg.pkl"), "wb") as f:
    pickle.dump(res_ksg, f, protocol=4)
print(f"Saved → {config.var_path('landscape_ksg.pkl')}")


# %% ── METHOD 2: KSG — lead-lag significance (shuffle at peak window) ─────────
print(f"Lead-lag significance: {N_SHUF_LEADLAG} shuffles per neuron (KSG)")
sig_ksg = np.zeros(len(run_idx), dtype=bool)
t0 = time.time()

for pi, i in enumerate(run_idx):
    mi_arr = mi_ksg[pi]
    if np.all(np.isnan(mi_arr)): continue
    mi_sm   = _smooth(mi_arr, SMOOTH_WIN)          # peak window from smoothed curve
    if np.all(np.isnan(mi_sm)): continue
    pk_wi   = int(np.nanargmax(mi_sm))
    real_mi = mi_arr[pk_wi]                         # significance uses the actual (raw) MI at that window
    if real_mi <= 0: continue

    n    = NEURONS[i]
    mask = mask_list[i]
    if mask.sum() < MIN_TRIALS: continue
    kin_m = kin_list[i][mask]

    c_pk = mi_io.counts(n, (float(win_starts[pk_wi]),
                             float(win_starts[pk_wi] + WIN_WIDTH)))[mask]
    if np.unique(c_pk).size < 2: continue

    null = np.array([_ksg(rng.permutation(kin_m), c_pk, K_KSG)
                     for _ in range(N_SHUF_LEADLAG)])
    sig_ksg[pi] = real_mi > np.percentile(null, 95)

print(f"  {sig_ksg.sum()}/{len(run_idx)} significant  ({time.time()-t0:.1f}s)")


# %% ── METHOD 2: KSG — MI landscape entropy (selectivity filter) ─────────────
# Entropy of the MI-vs-time curve (treating the positive MI as a mass distribution).
# Low entropy = peaked landscape (clear preferred time).
# High entropy = flat / noisy landscape (no clear peak).
# Inspect the distribution below, then set ENT_CUTOFF_KSG and re-run the plots cell.
ENT_CUTOFF_KSG = None   # e.g. 3.5 to keep only peaked neurons; None = no entropy filter

ent_ksg = np.full(len(run_idx), np.nan)
for pi in range(len(run_idx)):
    mi_arr = mi_ksg[pi]
    if np.all(np.isnan(mi_arr)): continue
    mi_pos = np.clip(mi_arr, 0, None)
    total  = mi_pos.sum()
    if total <= 0: continue
    p = mi_pos / total
    p = p[p > 0]
    ent_ksg[pi] = -np.sum(p * np.log2(p))

fig, ax = plt.subplots(figsize=(6, 4))
for reg in ["CP", "CB"]:
    idx = [i for i, r in enumerate(regions_run) if r == reg]
    e   = ent_ksg[idx]; e = e[~np.isnan(e)]
    ax.hist(e, bins=30, color=CLR[reg], alpha=0.6, label=f"{reg} (n={len(e)})", density=True)
    med = np.median(e)
    ax.axvline(med, color=CLR[reg], lw=1.5, ls="--", label=f"{reg} median = {med:.2f} bits")
if ENT_CUTOFF_KSG is not None:
    ax.axvline(ENT_CUTOFF_KSG, color="k", lw=1.5, ls="-", label=f"cutoff = {ENT_CUTOFF_KSG:.2f}")
ax.set(xlabel="MI landscape entropy (bits)", ylabel="# neurons",
       title=f"KSG k={K_KSG}  —  MI curve entropy per neuron\n"
             f"low = peaked  |  high = flat/noisy  |  max = {np.log2(n_wins):.1f} bits ({n_wins} windows)",
             xlim=(3, 6.5))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir("landscapes", "ksg") / "entropy_dist.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

print(f"{'Region':6s}  {'n':>5s}  {'p10':>5s}  {'p25':>5s}  {'p50':>5s}  {'p75':>5s}  {'p90':>5s}  bits")
for reg in ["CP", "CB"]:
    idx = [i for i, r in enumerate(regions_run) if r == reg]
    e   = ent_ksg[idx]; e = e[~np.isnan(e)]
    qs  = np.percentile(e, [10, 25, 50, 75, 90])
    print(f"{reg:6s}  {len(e):5d}  {qs[0]:5.2f}  {qs[1]:5.2f}  {qs[2]:5.2f}  {qs[3]:5.2f}  {qs[4]:5.2f}")
print(f"→ set ENT_CUTOFF_KSG above, then re-run this cell and the plots cell")


# %% ── METHOD 2: KSG — plots ─────────────────────────────────────────────────
ind_dir = config.fig_dir("landscapes", "ksg", "landscapes_individual")
for pi, i in enumerate(run_idx):
    _plot_individual(uids[i], regions[i], mi_ksg[pi], win_centers, f"KSG k={K_KSG}", ind_dir,
                     sig=sig_ksg[pi], entropy=ent_ksg[pi])

# %% ── METHOD 2: KSG — plots mean / lead-lag ──────────────────────────────────
_plot_mean_landscape(mi_ksg, regions_run, win_centers, f"KSG k={K_KSG}",
                     config.fig_dir("landscapes", "ksg"))
# combine significance filter + optional entropy cutoff
sel_ksg = sig_ksg & (ent_ksg <= ENT_CUTOFF_KSG) if ENT_CUTOFF_KSG is not None else sig_ksg
_plot_lead_lag(mi_ksg, regions_run, win_centers, f"KSG k={K_KSG}",
               config.fig_dir("landscapes", "ksg"), sig=sel_ksg, peak_range=PEAK_RANGE)
print("KSG plots done.")


"""
 _____ ______   ___          ________  ________  ________   ________
|\   _ \  _   \|\  \        |\   __  \|\   __  \|\   ____\ |\   ____\
\ \  \\\__\ \  \ \  \       \ \  \|\  \ \  \|\  \ \  \___|_\ \  \___|_
 \ \  \\|__| \  \ \  \       \ \   _  _\ \  \\\  \ \_____  \\ \_____  \
  \ \  \    \ \  \ \  \       \ \  \\  \\ \  \\\  \|____|\  \\|____|\  \
   \ \__\    \ \__\ \__\       \ \__\\ _\\ \_______\____\_\  \ ____\_\  \
    \|__|     \|__|\|__|        \|__|\|__|\|_______|\_________\\_________\
                                                   \|_________\|_________|


"""
# %% ── METHOD 3: ROSS — CP only, cross-check against KSG ─────────────────────
print("\n=== METHOD: Ross (CP only) ===")
cp_run_idx   = [pi for pi, i in enumerate(run_idx) if regions[i] == "CP"]
mi_ross      = np.full((len(run_idx), n_wins), np.nan)   # same shape, CB rows stay nan
t0 = time.time()

for count_cp, pi in enumerate(cp_run_idx):
    i     = run_idx[pi]
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_m = kin_list[i][mask]
    if mask.sum() < MIN_TRIALS: continue

    C   = counts_landscape(n["spikes"], win_starts, WIN_WIDTH)
    C_m = C[mask]

    for wi in range(n_wins):
        c_w = C_m[:, wi]
        if np.unique(c_w).size < 2: continue
        mi_ross[pi, wi] = _ross(kin_m, c_w, K_KSG)

    if (count_cp + 1) % 5 == 0:
        elapsed = time.time() - t0
        rate = (count_cp+1) / elapsed
        print(f"  {count_cp+1}/{len(cp_run_idx)} CP neurons  ({elapsed:.0f}s elapsed, "
              f"~{(len(cp_run_idx)-count_cp-1)/rate:.0f}s remaining)", flush=True)

print(f"Done in {time.time()-t0:.1f}s")

res_ross = dict(method="ross", win_centers=win_centers, win_starts=win_starts,
                WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP, mi=mi_ross,
                uids=[uids[i] for i in run_idx], regions=regions_run,
                COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), K=K_KSG)
with open(config.var_path("landscape_ross.pkl"), "wb") as f:
    pickle.dump(res_ross, f, protocol=4)
print(f"Saved → {config.var_path('landscape_ross.pkl')}")


# %% ── METHOD 3: ROSS — lead-lag significance (CP only) ──────────────────────
print(f"Lead-lag significance: {N_SHUF_LEADLAG} shuffles per neuron (Ross, CP only)")
sig_ross = np.zeros(len(run_idx), dtype=bool)
t0 = time.time()

for count_cp, pi in enumerate(cp_run_idx):
    mi_arr = mi_ross[pi]
    if np.all(np.isnan(mi_arr)): continue
    mi_sm   = _smooth(mi_arr, SMOOTH_WIN)          # peak window from smoothed curve
    if np.all(np.isnan(mi_sm)): continue
    pk_wi   = int(np.nanargmax(mi_sm))
    real_mi = mi_arr[pk_wi]                         # significance uses the actual (raw) MI at that window
    if real_mi <= 0: continue

    i    = run_idx[pi]
    n    = NEURONS[i]
    mask = mask_list[i]
    if mask.sum() < MIN_TRIALS: continue
    kin_m = kin_list[i][mask]

    c_pk = mi_io.counts(n, (float(win_starts[pk_wi]),
                             float(win_starts[pk_wi] + WIN_WIDTH)))[mask]
    if np.unique(c_pk).size < 2: continue

    null = np.array([_ross(rng.permutation(kin_m), c_pk, K_KSG)
                     for _ in range(N_SHUF_LEADLAG)])
    sig_ross[pi] = real_mi > np.percentile(null, 95)

print(f"  {sig_ross[cp_run_idx].sum()}/{len(cp_run_idx)} CP significant  ({time.time()-t0:.1f}s)")


# %% ── METHOD 3: ROSS — plots + KSG vs Ross difference (CP) ──────────────────
ind_dir = config.fig_dir("ross_cp", "landscapes_individual")
for pi in cp_run_idx:
    _plot_individual(uids[run_idx[pi]], "CP", mi_ross[pi], win_centers,
                     f"Ross k={K_KSG} (CP)", ind_dir, sig=sig_ross[pi])

_plot_mean_landscape(mi_ross[cp_run_idx], ["CP"]*len(cp_run_idx),
                     win_centers, f"Ross k={K_KSG} (CP only)",
                     config.fig_dir("landscapes", "ross_cp"))
_plot_lead_lag(mi_ross[cp_run_idx], ["CP"]*len(cp_run_idx),
               win_centers, f"Ross k={K_KSG} (CP only)",
               config.fig_dir("landscapes", "ross_cp"), sig=sig_ross[cp_run_idx], peak_range=PEAK_RANGE)

# difference landscape: KSG - Ross (CP only)
diff = mi_ksg[cp_run_idx] - mi_ross[cp_run_idx]
fig, ax = plt.subplots(figsize=(8, 4))
mn  = np.nanmean(diff, axis=0)
sem = np.nanstd(diff, axis=0) / np.sqrt(np.maximum(np.sum(~np.isnan(diff), axis=0), 1))
ax.plot(win_centers, mn, color=CLR["CP"], lw=2, label="KSG − Ross")
ax.fill_between(win_centers, mn-sem, mn+sem, color=CLR["CP"], alpha=0.2)
ax.axhline(0, color="k", lw=0.8, ls="--")
ax.axvline(0, color="k", lw=0.8, ls=":")
ax.set(xlabel="window center re reach (s)", ylabel="MI difference (bits)",
       title=f"KSG vs Ross cross-check (CP, n={len(cp_run_idx)})\n"
             f"positive = KSG > Ross; expected ~0 for qualitative agreement")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir("landscapes", "ross_cp") / "diff_ksg_ross_CP.png",
                dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)
print("Ross plots done.")
