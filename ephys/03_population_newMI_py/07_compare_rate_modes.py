# %% ── 0. imports & constants ─────────────────────────────────────────────────
# Compare the THREE neural-variable modes (count / isi / conv, see mi_rates.py)
# using the SAME fixed-window pipeline and parameters as 04_MI_fixed_win.py.
# Methods: Binning (plug-in) and KSG only -- Ross needs a discrete variable and so
# cannot be applied to the continuous "isi"/"conv" rates.
# Outputs: a few example landscapes, the mean landscape, and the peak (lead/lag)
# distribution -- each with the three modes overlaid. No per-neuron figures.
import numpy as np
import matplotlib.pyplot as plt
import pickle, time, importlib
import mi_io, config, mi_rates
importlib.reload(config); importlib.reload(mi_rates)
from sklearn.feature_selection import mutual_info_regression
from scipy.signal import find_peaks
from joblib import Parallel, delayed

if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

CLR      = {"CP": np.array([106,176,134])/255, "CB": np.array([154,123,193])/255}
MODE_CLR = {"count": "#444444", "isi": "#d95f02", "conv": "#1f78b4"}
NAT2BIT  = 1 / np.log(2)
rng      = np.random.default_rng(0)

# ── analysis parameters — KEEP IN SYNC WITH 04_MI_fixed_win.py ────────────────
WIN_WIDTH   = 0.100       # sliding window width (s)
WIN_STEP    = 0.025       # step size (s)
WIN_RANGE   = (-1.0, 1.0) # landscape range (s), re reach
K_KSG       = 3           # k for KSG estimator
N_BINS_KIN  = 16          # kin bins for the plug-in estimator
N_BINS_NEURAL = 16        # plug-in bins for a CONTINUOUS neural variable
COORD       = 1           # 0=x, 1=y med-lat, 2=z
FRAMES      = (147, 167)  # endpoint window, MATLAB 1-indexed
CAT_KEEP    = {1, 2}      # reach categories to include
MIN_TRIALS  = 20

N_SHUF_LEADLAG    = 50    # shuffles per neuron for peak significance
SMOOTH_WIN        = 2     # moving-average kernel (# windows) for peak detection
PEAK_PROM_RATIO   = 0.5   # local check  (see 04 for the rationale)
PEAK_RUNNERUP_MAX = 0.7   # global check (see 04 for the rationale)
PEAK_RANGE        = (-0.3, 0.3)  # e.g. (-0.3, 0.3); None = full WIN_RANGE
RATE_BIN   = mi_rates.RATE_BIN     # 2 ms fine grid for the rate modes
ISI_SIGMA  = mi_rates.ISI_SIGMA    # 20 ms Gaussian on 1/ISI (MATLAB sig_pk)
CONV_SIGMA = 0.03 #mi_rates.CONV_SIGMA   # 70 ms Gaussian on the spike train

DO_NULL     = True        # per-window shuffle null for every mode x method (drawn as
                          # a dashed line in the matching mode colour). 6 combinations,
                          # so this is the slowest part: ~33 min for the full population
                          # at N_JOBS=-1, ~56 min at 8 (S=10).
N_SHUF_NULL = 10          # shuffles per WINDOW for the null
N_JOBS      = -1          # cores for the null (-1 = all); see 04 for the scaling numbers

MODES   = ("count", "isi", "conv")
METHODS = ("bin", "ksg")

N_TEST      = None          # neurons PER REGION. Set None for the full run --
                          # 3 modes x 2 methods over ~2900 neurons takes a while,
                          # so test small first and check the timing printout.
N_EXAMPLES  = 4           # example neurons drawn in the example figure
SHOW_FIGS   = True
SAVE_FIGS   = True
SAVE_VARS   = True
CACHE_RATES = True        # cache the windowed isi/conv variable, keyed by the
                          # window config in the filename. First run at a given
                          # window builds it for the FULL population (~1 min per
                          # mode) and saves; later runs at that window just load
                          # it. Set False for a quick throwaway test run.

OUT_SUBDIR  = "rate_mode_MI_littleConv"   # figures folder under OUT_DIR
RUN_NAME    = "fixedwin100_littleConv"    # tags saved pickles; keep in sync with 04 / 06

win_starts  = np.arange(WIN_RANGE[0], WIN_RANGE[1] - WIN_WIDTH + WIN_STEP/2, WIN_STEP)
win_centers = win_starts + WIN_WIDTH / 2
n_wins      = len(win_starts)
print(f"windows: {n_wins}  width={WIN_WIDTH*1000:.0f} ms  step={WIN_STEP*1000:.0f} ms")


# %% ── 1. select neurons + endpoint kin ───────────────────────────────────────
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
regions_run = [regions[i] for i in run_idx]


# %% ── 2. helpers (same definitions as 04, kept local per project style) ──────
def _ksg(kin_m, c_w, k):
    X = kin_m.reshape(-1, 1)
    return mutual_info_regression(X, c_w + rng.normal(0, 1e-9, c_w.size),
                                  n_neighbors=k, random_state=0)[0] * NAT2BIT

def _quantile_codes(v, nb):
    edges = np.unique(np.quantile(v, np.linspace(0, 1, nb + 1)))
    if len(edges) < 3:
        return np.zeros(len(v), dtype=int), 1
    return np.clip(np.digitize(v, edges[1:-1]), 0, len(edges) - 2), len(edges) - 1

def _bin(kin_m, c_w, nb, discrete):
    """Plug-in MI. kin always quantile-binned; the neural variable keeps its exact
    integer levels when discrete ("count"), else is quantile-binned."""
    kb, n_k = _quantile_codes(kin_m, nb)
    c_w = np.asarray(c_w)
    if discrete:
        cv = np.unique(c_w); cm = {v: i for i, v in enumerate(cv)}
        ci = np.array([cm[v] for v in c_w]); n_c = len(cv)
    else:
        ci, n_c = _quantile_codes(c_w, N_BINS_NEURAL)
    J  = np.zeros((n_c, n_k)); np.add.at(J, (ci, kb), 1.0); P = J / J.sum()
    pc = P.sum(1, keepdims=True); pk = P.sum(0, keepdims=True); nz = P > 0
    return float(np.sum(P[nz] * np.log2(P[nz] / (pc @ pk + 1e-300)[nz])))

def _smooth(arr, win):
    if win <= 1: return arr.copy()
    n = len(arr); out = np.full(n, np.nan); half = win // 2
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        seg = arr[lo:hi]
        if np.any(~np.isnan(seg)):
            out[i] = np.nanmean(seg)
    return out

def _find_peak(arr, prom_ratio=PEAK_PROM_RATIO, runnerup_max=PEAK_RUNNERUP_MAX):
    """Tallest INTERIOR local max that is both locally clean and globally unique.
    Returns (idx, prominence_ratio, height) or None. Identical to 04."""
    valid = ~np.isnan(arr)
    if valid.sum() < 3:
        return None
    a = np.where(valid, arr, -np.inf)
    idx_peaks, props = find_peaks(a, prominence=0, distance=max(1, SMOOTH_WIN))
    if len(idx_peaks) == 0:
        return None
    proms = props["prominences"]; heights = a[idx_peaks]
    order = np.argsort(heights)[::-1]; top = order[0]
    top_h, top_prom = heights[top], proms[top]
    runner_h = heights[order[1]] if len(order) > 1 else 0.0
    if top_h <= 0 or top_prom / top_h < prom_ratio or runner_h > runnerup_max * top_h:
        return None
    return int(idx_peaks[top]), float(top_prom / top_h), float(top_h)

def _null_one(kin_m, C_m, method, n_shuf, seed, k, nb, discrete):
    """Per-window shuffle null for ONE neuron -> (n_wins,): the MEAN MI over
    `n_shuf` trial-permutations of kin at every window. Same permutations reused
    across windows (unbiased per window, smoother curve). Top-level so joblib can
    ship it to worker processes."""
    r    = np.random.default_rng(seed)
    shuf = [r.permutation(kin_m) for _ in range(n_shuf)]
    out  = np.full(C_m.shape[1], np.nan)
    for wi in range(C_m.shape[1]):
        c_w = C_m[:, wi]
        if np.unique(c_w).size < 2: continue
        out[wi] = np.mean([_bin(ks, c_w, nb, discrete) if method == "bin"
                           else _ksg(ks, c_w, k) for ks in shuf])
    return out


def _mi_call(method, kin_m, c_w, discrete):
    return (_bin(kin_m, c_w, N_BINS_KIN, discrete) if method == "bin"
            else _ksg(kin_m, c_w, K_KSG))


# %% ── 3. neural variable per mode (load cache, else build it and save) ──────
# "count" is always derived (one searchsorted). For isi/conv: if a cache for THIS
# window config exists and still matches the smoothing params, load it; otherwise
# build it for the FULL population and save, so every later run at this window is
# instant. Building the full population even under N_TEST is deliberate -- a cache
# covering only the test subset would be silently useless to a later full run.
all_uids = [n["uid"] for n in NEURONS]
WIN_VAR = {}

for mode in MODES:
    t0 = time.time()

    if mode == "count":
        WIN_VAR[mode] = {i: np.asarray(mi_rates.windowed(
            NEURONS[i]["spikes"], win_starts, WIN_WIDTH, "count"), np.float32)
            for i in run_idx}
        print(f"count  derived ({len(run_idx)} neurons, {time.time()-t0:.0f}s)")
        continue

    path = config.var_path(mi_rates.cache_name(mode, WIN_WIDTH, WIN_STEP, WIN_RANGE,
                                       RATE_BIN, ISI_SIGMA, CONV_SIGMA))
    if path.exists():
        with open(path, "rb") as f:
            res = pickle.load(f)
        if mi_rates.cache_matches(res, all_uids, win_starts, WIN_WIDTH,
                                  RATE_BIN, ISI_SIGMA, CONV_SIGMA):
            WIN_VAR[mode] = {i: res["X"][i] for i in run_idx}
            del res
            print(f"{mode:6s} loaded  {path.name}  ({time.time()-t0:.0f}s)")
            continue
        del res
        print(f"{mode:6s} cache exists but window/smoothing differ - rebuilding")

    build = range(len(NEURONS)) if CACHE_RATES else run_idx
    if CACHE_RATES:
        print(f"{mode:6s} building cache for all {len(NEURONS)} neurons "
              f"(~1 min; set CACHE_RATES=False to skip)...", flush=True)
    X = {i: np.asarray(mi_rates.windowed(
            NEURONS[i]["spikes"], win_starts, WIN_WIDTH, mode, bin_width=RATE_BIN,
            isi_sigma=ISI_SIGMA, conv_sigma=CONV_SIGMA), np.float32) for i in build}
    WIN_VAR[mode] = {i: X[i] for i in run_idx}

    if CACHE_RATES:
        res = dict(mode=mode, X=[X[i] for i in range(len(NEURONS))], uids=all_uids,
                   win_starts=win_starts, WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP,
                   WIN_RANGE=WIN_RANGE, RATE_BIN=RATE_BIN,
                   ISI_SIGMA=ISI_SIGMA, CONV_SIGMA=CONV_SIGMA)
        with open(path, "wb") as f:
            pickle.dump(res, f, protocol=4)
        del res
        print(f"{mode:6s} built + saved {path.name} "
              f"({path.stat().st_size/1e6:.0f} MB, {time.time()-t0:.0f}s)")
    else:
        print(f"{mode:6s} derived, not cached ({len(run_idx)} neurons, "
              f"{time.time()-t0:.0f}s)")
    del X


# %% ── 4. MI landscapes: 3 modes x 2 methods ─────────────────────────────────
MI = {(m, meth): np.full((len(run_idx), n_wins), np.nan)
      for m in MODES for meth in METHODS}
t_all = time.time()

for mode in MODES:
    discrete = mi_rates.is_discrete(mode)
    for method in METHODS:
        t0 = time.time()
        for pi, i in enumerate(run_idx):
            mask  = mask_list[i]
            if mask.sum() < MIN_TRIALS: continue
            kin_m = kin_list[i][mask]
            C_m   = np.asarray(WIN_VAR[mode][i], float)[mask]
            for wi in range(n_wins):
                c_w = C_m[:, wi]
                if np.unique(c_w).size < 2: continue
                MI[(mode, method)][pi, wi] = _mi_call(method, kin_m, c_w, discrete)
            if (pi + 1) % max(1, len(run_idx) // 5) == 0:
                el = time.time() - t0
                print(f"  {mode}/{method}: {pi+1}/{len(run_idx)}  {el:.0f}s "
                      f"(~{el/(pi+1)*(len(run_idx)-pi-1):.0f}s left)", flush=True)
        print(f"{mode:6s} {method:4s} done in {time.time()-t0:.0f}s")

print(f"\nall landscapes done in {time.time()-t_all:.0f}s")


# %% ── 4b. per-window shuffle null for every mode x method ───────────────────
# Null MI at EVERY window (not only at the peak), so the landscapes can be read
# against chance. Parallelised across neurons; this is the slowest cell here.
NULL = {}
if DO_NULL:
    t_all = time.time()
    for mode in MODES:
        discrete = mi_rates.is_discrete(mode)
        payload = []
        for pi, i in enumerate(run_idx):
            mask = mask_list[i]
            if mask.sum() < MIN_TRIALS:
                payload.append(None); continue
            payload.append((kin_list[i][mask],
                            np.asarray(WIN_VAR[mode][i], float)[mask]))
        keep = [pi for pi, it in enumerate(payload) if it is not None]
        for method in METHODS:
            t0 = time.time()
            res_n = Parallel(n_jobs=N_JOBS, backend="loky", verbose=0)(
                delayed(_null_one)(payload[pi][0], payload[pi][1], method,
                                   N_SHUF_NULL, 1000 + pi, K_KSG, N_BINS_KIN, discrete)
                for pi in keep)
            arr = np.full((len(run_idx), n_wins), np.nan)
            for pi, v in zip(keep, res_n):
                arr[pi] = v
            NULL[(mode, method)] = arr
            print(f"  null {mode:6s} {method:4s}: {time.time()-t0:.0f}s", flush=True)
    print(f"all nulls done in {time.time()-t_all:.0f}s "
          f"(S={N_SHUF_NULL}, {N_JOBS} jobs)")
else:
    print("DO_NULL is False - no per-window null computed")


# %% ── 5. peak significance (shuffle at each neuron's peak window) ───────────
SIG = {}
for mode in MODES:
    discrete = mi_rates.is_discrete(mode)
    for method in METHODS:
        t0 = time.time()
        sig = np.zeros(len(run_idx), dtype=bool)
        for pi, i in enumerate(run_idx):
            mi_arr = MI[(mode, method)][pi]
            if np.all(np.isnan(mi_arr)): continue
            pk = _find_peak(_smooth(mi_arr, SMOOTH_WIN))
            if pk is None: continue
            pk_wi   = pk[0]
            real_mi = mi_arr[pk_wi]
            if real_mi <= 0: continue
            mask = mask_list[i]
            if mask.sum() < MIN_TRIALS: continue
            kin_m = kin_list[i][mask]
            c_pk  = np.asarray(WIN_VAR[mode][i], float)[mask][:, pk_wi]
            if np.unique(c_pk).size < 2: continue
            null = np.array([_mi_call(method, rng.permutation(kin_m), c_pk, discrete)
                             for _ in range(N_SHUF_LEADLAG)])
            sig[pi] = real_mi > np.percentile(null, 95)
        SIG[(mode, method)] = sig
        print(f"{mode:6s} {method:4s}: {sig.sum()}/{len(run_idx)} significant "
              f"({time.time()-t0:.0f}s)")

if SAVE_VARS:
    out = dict(modes=list(MODES), methods=list(METHODS),
               mi={f"{m}|{meth}": MI[(m, meth)] for m in MODES for meth in METHODS},
               sig={f"{m}|{meth}": SIG[(m, meth)] for m in MODES for meth in METHODS},
               null={f"{m}|{meth}": NULL[(m, meth)] for m in MODES for meth in METHODS
                     if (m, meth) in NULL}, N_SHUF_NULL=N_SHUF_NULL,
               win_starts=win_starts, win_centers=win_centers,
               WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP,
               uids=[uids[i] for i in run_idx], regions=regions_run,
               COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), K=K_KSG,
               N_SHUF_LEADLAG=N_SHUF_LEADLAG, N_TEST=N_TEST)
    p = config.var_path(f"{RUN_NAME}_ratemodes_MI.pkl")
    with open(p, "wb") as f:
        pickle.dump(out, f, protocol=4)
    print(f"saved -> {p}")


# %% ── 6. example landscapes (a few neurons, three modes overlaid) ───────────
fig_dir = config.fig_dir(OUT_SUBDIR)

# prefer neurons that are significant in as many mode/method combos as possible
score = np.sum([SIG[(m, meth)].astype(int) for m in MODES for meth in METHODS], axis=0)
ex_pi = list(np.argsort(score)[::-1][:N_EXAMPLES])

for pi in ex_pi:
    fig, axes = plt.subplots(1, 2, figsize=(13, 3.4), sharex=True)
    for ax, method in zip(axes, METHODS):
        for mode in MODES:
            row = MI[(mode, method)][pi]
            ax.plot(win_centers, _smooth(row, SMOOTH_WIN), lw=1.7,
                    color=MODE_CLR[mode],
                    label=f"{mode}{' *' if SIG[(mode, method)][pi] else ''}")
            ax.plot(win_centers, row, lw=0.7, alpha=0.3, color=MODE_CLR[mode])
            if (mode, method) in NULL:
                ax.plot(win_centers, NULL[(mode, method)][pi], lw=1.0, ls="--",
                        color=MODE_CLR[mode], alpha=0.55)
        ax.axvline(0, color="k", lw=0.8, ls=":")
        ax.axhline(0, color="k", lw=0.4)
        ax.set(xlabel="window center re reach (s)", ylabel="MI (bits)",
               title=f"{method.upper()}")
        ax.legend(fontsize=7)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.suptitle(f"{uids[run_idx[pi]]}  ({regions_run[pi]})   * = significant peak",
                 fontsize=10)
    plt.tight_layout()
    if SAVE_FIGS:
        fig.savefig(fig_dir / f"example_{uids[run_idx[pi]].replace('|','_')}.png",
                    dpi=130, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 7. mean landscape: rows = method, cols = region, lines = mode ─────────
fig, axes = plt.subplots(len(METHODS), 2, figsize=(13, 7), sharex=True)
axes = np.atleast_2d(axes)
for r, method in enumerate(METHODS):
    for c, reg in enumerate(("CP", "CB")):
        ax  = axes[r, c]
        idx = [i for i, rr in enumerate(regions_run) if rr == reg]
        for mode in MODES:
            if not idx: continue
            m  = MI[(mode, method)][idx]
            ok = ~np.all(np.isnan(m), axis=1)
            if ok.sum() == 0: continue
            m  = m[ok]
            mn  = np.nanmean(m, axis=0)
            sem = np.nanstd(m, axis=0) / np.sqrt(
                np.maximum(np.sum(~np.isnan(m), axis=0), 1))
            ax.plot(win_centers, mn, lw=2, color=MODE_CLR[mode],
                    label=f"{mode} (n={ok.sum()})")
            ax.fill_between(win_centers, mn - sem, mn + sem,
                            color=MODE_CLR[mode], alpha=0.18)
            if (mode, method) in NULL:
                nm   = NULL[(mode, method)][idx][ok]
                nmn  = np.nanmean(nm, axis=0)
                nsem = np.nanstd(nm, axis=0) / np.sqrt(
                    np.maximum(np.sum(~np.isnan(nm), axis=0), 1))
                ax.plot(win_centers, nmn, lw=1.2, ls="--", color=MODE_CLR[mode],
                        alpha=0.5,
                        label="shuffle null (dashed)" if mode == MODES[0] else None)
                ax.fill_between(win_centers, nmn - nsem, nmn + nsem,
                                color=MODE_CLR[mode], alpha=0.06)
            pk = win_centers[np.nanargmax(mn)]
            ax.axvline(pk, color=MODE_CLR[mode], lw=1.1, ls="--", alpha=0.8)
        ax.axvline(0, color="k", lw=0.8, ls=":")
        ax.axhline(0, color="k", lw=0.4)
        ax.set(title=f"{method.upper()}  -  {reg}", ylabel="mean MI +- SEM (bits)")
        if r == len(METHODS) - 1:
            ax.set_xlabel("window center re reach (s)")
        if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=8)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.suptitle(f"Mean MI landscape by neural-variable mode  "
             f"(win {WIN_WIDTH*1000:.0f} ms, step {WIN_STEP*1000:.0f} ms)", fontsize=11)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(fig_dir / "mean_landscape_by_mode.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 8. peak (lead/lag) distributions: rows = method, cols = region ────────
lo, hi = PEAK_RANGE if PEAK_RANGE is not None else (win_centers[0], win_centers[-1])
win_mask  = (win_centers >= lo) & (win_centers <= hi)
bin_edges = np.linspace(lo, hi, 21)

fig, axes = plt.subplots(len(METHODS), 2, figsize=(13, 7), sharex=True)
axes = np.atleast_2d(axes)
peak_summary = []
for r, method in enumerate(METHODS):
    for c, reg in enumerate(("CP", "CB")):
        ax  = axes[r, c]
        idx = [i for i, rr in enumerate(regions_run) if rr == reg]
        for mode in MODES:
            peaks = []
            for pi in idx:
                if not SIG[(mode, method)][pi]: continue
                row = MI[(mode, method)][pi]
                if not np.any(~np.isnan(row)): continue
                rs = np.where(win_mask, _smooth(row, SMOOTH_WIN), np.nan)
                pk = _find_peak(rs)
                if pk is not None:
                    peaks.append(win_centers[pk[0]])
            if peaks:
                ax.hist(peaks, bins=bin_edges, histtype="step", lw=1.8,
                        color=MODE_CLR[mode], label=f"{mode} (n={len(peaks)})")
                med = np.median(peaks)
                ax.axvline(med, color=MODE_CLR[mode], lw=1.4, ls="--", alpha=0.9)
                peak_summary.append((method, reg, mode, len(peaks), med))
            else:
                peak_summary.append((method, reg, mode, 0, np.nan))
        ax.axvline(0, color="k", lw=0.8, ls=":")
        ax.set(title=f"{method.upper()}  -  {reg}", ylabel="# neurons")
        if r == len(METHODS) - 1:
            ax.set_xlabel("time of peak MI re reach (s)")
        if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=8)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.suptitle("Peak-MI (lead/lag) distribution by neural-variable mode  "
             "(significant peaks only; dashed = median)", fontsize=11)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(fig_dir / "peak_distribution_by_mode.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

print(f"\n{'method':7s} {'region':7s} {'mode':6s} {'n_sig_peaks':>12s} {'median_ms':>10s}")
for method, reg, mode, npk, med in peak_summary:
    med_s = f"{med*1000:.0f}" if np.isfinite(med) else "-"
    print(f"{method:7s} {reg:7s} {mode:6s} {npk:12d} {med_s:>10s}")
print(f"\nfigures -> {fig_dir}")

# %%
