# %% ── 0. imports & constants ─────────────────────────────────────────────────
import numpy as np
import matplotlib.pyplot as plt
import pickle, time, importlib
import mi_io, config, mi_rates
importlib.reload(config)
from sklearn.feature_selection import mutual_info_regression, mutual_info_classif
from scipy.signal import find_peaks
from joblib import Parallel, delayed

if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

CLR     = {"CP": np.array([106,176,134])/255, "CB": np.array([154,123,193])/255}
NAT2BIT = 1 / np.log(2)
rng     = np.random.default_rng(0)

# ── analysis parameters — change here only ────────────────────────────────────
WIN_WIDTH   = 0.1       # sliding window width (s)
WIN_STEP    = 0.025       # step size (s)  →  temporal resolution of the landscape
WIN_RANGE   = (-1.0, 1.0) # landscape range (s), re reach onset
K_KSG       = 3           # k for KSG estimator
N_BINS_KIN  = 16          # kin bins for plug-in estimator
COORD       = 1           # kinematic coord: 0=x, 1=y med-lat, 2=z
FRAMES      = (147, 167)  # endpoint window, MATLAB 1-indexed
CAT_KEEP    = {1, 2}      # reach categories to include
MIN_TRIALS  = 20          # min valid trials to attempt MI

NEURAL_VAR  = "count"     # NEURAL VARIABLE fed to the MI estimators. One of:
                          #   "count" -> spike COUNT in the window (project default,
                          #              integer; the only mode Ross can use)
                          #   "isi"   -> 1/ISI instantaneous FR, lightly Gaussian
                          #              smoothed per trial (matches the MATLAB
                          #              pipeline: 2 ms bins, sigma 20 ms)
                          #   "conv"  -> spike train convolved with a Gaussian
                          #              (sigma 70 ms), window-averaged
                          # "isi"/"conv" are CONTINUOUS (mean Hz over the window):
                          # the plug-in estimator quantile-bins them (N_BINS_NEURAL)
                          # and Ross is disabled. See mi_rates.py.
N_BINS_NEURAL = 16        # plug-in bins for the neural variable when it is continuous
                          # (ignored for "count", which uses its exact integer levels)

N_TEST           = None      # neurons per region for test run; set None for all neurons
N_SHUF_LEADLAG   = 50     # shuffles per neuron for lead-lag significance (peak window only)
SMOOTH_WIN       = 3      # moving-average kernel (# windows) for peak detection on individual plots.
                          # A wider WIN_WIDTH (75ms) gives a less jagged raw MI curve than 50ms
                          # did, so this can go back down to 3 (~75ms) — narrowing WIN_WIDTH
                          # while ALSO smoothing more heavily was fighting itself (a narrower
                          # count window makes raw counts noisier; heavier smoothing then
                          # flattens real peaks too, for the wrong reason).
PEAK_PROM_RATIO  = 0.5    # a local max counts as a real "peak" only if its prominence
                          # (height above the deepest valley separating it from any
                          # taller neighbor) is >= this fraction of its own height. This
                          # is a purely LOCAL check ("does this bump drop to near-baseline
                          # on both sides") — it says nothing about whether OTHER bumps
                          # elsewhere in the curve are comparably tall, so a curve with
                          # several isolated, similar-height bumps can have every single
                          # one individually score ratio~1.0.
PEAK_RUNNERUP_MAX = 0.7   # the GLOBAL check that fixes the above: the second-tallest
                          # interior peak's height must not exceed this fraction of the
                          # tallest peak's height, or the "peak" is considered not
                          # prominent overall (several comparable bumps -> no significant
                          # peak), even if it individually passed PEAK_PROM_RATIO. Lighter
                          # touch than the earlier pure top-vs-runner-up "dominance" gate
                          # (which needed runner-up <= 0.65x top and rejected almost every
                          # neuron population-wide) — this only rejects clearly-tied cases.
                          # Boundary points (peak sitting at the very edge of the window —
                          # an edge artifact, not a real temporal peak) are excluded
                          # automatically since find_peaks only reports interior local
                          # maxima. A minimum peak-to-peak distance (tied to SMOOTH_WIN,
                          # see _find_peak) stops a single-bin noise wobble from being
                          # counted as a separate peak.
N_SHUF_NULL      = 10     # shuffles per WINDOW for the per-bin shuffle null (the grey
                          # dashed reference line). Cost is ~N_SHUF_NULL x the landscape
                          # itself, so it is parallelised across neurons (N_JOBS).
NULL_METHODS     = ("bin", "ksg")   # Ross is excluded on purpose: it is already only a
                          # cross-check for KSG, and its per-window null would cost ~9 h
                          # even parallelised (its landscape alone is 1.9 s/neuron).
N_JOBS           = 8      # cores for the null. Benchmarked on 120 neurons (enough to keep
                          # every worker fed; KSG, S=10, extrapolated to the full
                          # population): 8 jobs 9.3 min, 12 -> 8.2, 16 -> 7.2, 24 -> 5.5.
                          # Scaling does NOT plateau after 8 -- an earlier benchmark that
                          # suggested it did was run on only 16 neurons, too few to fill
                          # 24 workers. 8 keeps the machine responsive while 04 runs;
                          # set -1 for the fastest unattended run (07 uses -1).
PEAK_RANGE       = None  # e.g. (-0.3, 0.3) to restrict lead-lag peak search/median to that window; None = full WIN_RANGE
SHOW_FIGS        = True   # set False for full run (only saves, faster)
SAVE_FIGS        = True
N_EXAMPLES       = 40     # individual landscapes drawn PER REGION (CP and CB).
                          # None = every neuron, one PNG each (slow: ~2900 files
                          # per method). Only the individual plots are affected;
                          # the MI, the mean landscape and the lead-lag figure
                          # always use all neurons.
EXAMPLE_MODE     = "spread"   # which ones, when N_EXAMPLES is not None:
                          # "spread" = evenly spaced over the neurons ranked by
                          #            their peak MI, so strong AND weak cells
                          #            are shown, not just the best
                          # "top"    = the N with the highest peak MI
                          # "random" = N at random (fixed seed)
                          # The choice is made ONCE, by the first method that is
                          # plotted, and reused by the others, so bin / KSG /
                          # Ross show the SAME neurons and can be compared.

RUN_TAG          = "v2"   # free suffix for the run, "" for none. "v2" = computed after the
                          # empty-trial fix in mi_io (2026-10-12), so these never land on
                          # top of the August results, which are named without it.

# The names are BUILT from the parameters, never typed: the width, the reach
# categories (only when they are not the default {1, 2}) and RUN_TAG. Changing
# WIN_WIDTH or CAT_KEEP therefore changes the folder and the pickle names with
# it, and a run can no longer be saved under another width's name.
#   WIN_WIDTH 0.1, CAT_KEEP {1, 2}  ->  fixed_win_100ms_v2 / fixedwin100_v2
#   WIN_WIDTH 0.1, CAT_KEEP {1}     ->  fixed_win_100ms_cat1_v2 / fixedwin100_cat1_v2
# Figure subfolders are "{RUN_NAME}_bin" / "_ksg" / "_ross_cp" under OUT_SUBDIR,
# and the pickles are "{RUN_NAME}_bin.pkl" etc.
_W_MS      = int(round(WIN_WIDTH * 1000))
_CAT_TAG   = "" if set(CAT_KEEP) == {1, 2} else "_cat" + "".join(str(c) for c in sorted(CAT_KEEP))
_RUN_TAG   = f"_{RUN_TAG}" if RUN_TAG else ""
OUT_SUBDIR = f"fixed_win_{_W_MS}ms{_CAT_TAG}{_RUN_TAG}"
RUN_NAME   = f"fixedwin{_W_MS}{_CAT_TAG}{_RUN_TAG}"
print(f"run: WIN_WIDTH {_W_MS} ms, cat {sorted(CAT_KEEP)}  ->  figures in "
      f"{OUT_SUBDIR}/, pickles {RUN_NAME}_*.pkl")


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

def neural_landscape(neuron, win_starts, WIN_WIDTH):
    """(n_trials, n_wins) neural variable for every sliding window, per NEURAL_VAR.
    "count" is bit-identical to counts_landscape above; "isi"/"conv" return the mean
    firing rate (Hz) in each window (see mi_rates.py)."""
    return mi_rates.windowed(neuron["spikes"], win_starts, WIN_WIDTH, NEURAL_VAR)

def neural_window(neuron, win):
    """Same neural variable for a SINGLE window (lo, hi) -> (n_trials,).
    Used by the shuffle-significance cells at each neuron's peak window.
    Deliberately takes the column out of the FULL landscape rather than
    recomputing on a short grid: for the smoothed rate modes a short grid can
    differ by a sample at the edges, and the significance test must be run on
    exactly the value the landscape peaked at."""
    lo, _ = win
    wi = int(np.argmin(np.abs(win_starts - lo)))
    return neural_landscape(neuron, win_starts, WIN_WIDTH)[:, wi]

def _ksg(kin_m, c_w, k):
    X = kin_m.reshape(-1, 1)
    return mutual_info_regression(X, c_w + rng.normal(0, 1e-9, c_w.size),
                                  n_neighbors=k, random_state=0)[0] * NAT2BIT

def _ross(kin_m, c_w, k):
    X = kin_m.reshape(-1, 1)
    return mutual_info_classif(X, c_w.astype(int),
                               n_neighbors=k, random_state=0)[0] * NAT2BIT

def _quantile_codes(v, nb):
    """Quantile-bin a continuous vector into <=nb levels -> (codes, n_levels)."""
    edges = np.unique(np.quantile(v, np.linspace(0, 1, nb + 1)))
    if len(edges) < 3:
        return np.zeros(len(v), dtype=int), 1
    return np.clip(np.digitize(v, edges[1:-1]), 0, len(edges) - 2), len(edges) - 1

def _bin(kin_m, c_w, nb, discrete=None):
    """Plug-in (binning) MI. kin is always quantile-binned into nb levels.
    The neural variable keeps its exact integer levels when discrete (NEURAL_VAR=
    "count" -- identical to the previous behaviour), otherwise it is quantile-binned
    into N_BINS_NEURAL levels (continuous "isi"/"conv" rates)."""
    if discrete is None:
        discrete = mi_rates.is_discrete(NEURAL_VAR)
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

def _find_peak(arr, prom_ratio=PEAK_PROM_RATIO, runnerup_max=PEAK_RUNNERUP_MAX):
    """Tallest INTERIOR local maximum in arr (NaN-safe), if it passes BOTH:
      - LOCAL check: its own prominence / height >= prom_ratio (drops to near-baseline
        on both sides, not riding on a shoulder of a taller nearby feature).
      - GLOBAL check: the second-tallest interior peak is no more than runnerup_max of
        its height (nothing else in the curve comes close) — this is what actually
        catches "several comparable, individually-clean bumps" landscapes, which the
        local check alone cannot (every isolated bump can score prom_ratio~1.0
        regardless of how many similar-height rivals exist elsewhere).
    Returns (idx, prominence_ratio, height) or None if no interior peak exists or
    either check fails ("not a significant peak").
    find_peaks only reports points with lower neighbors on BOTH sides, so a value
    sitting at the very edge of the window (e.g. a monotonic trend truncated by the
    window boundary) can never be returned here — that's an edge artifact, not a
    real temporal peak. A minimum peak-to-peak distance (tied to SMOOTH_WIN) stops a
    single-bin noise wobble from being counted as a separate, competing peak."""
    valid = ~np.isnan(arr)
    if valid.sum() < 3:
        return None
    a = np.where(valid, arr, -np.inf)
    idx_peaks, props = find_peaks(a, prominence=0, distance=max(1, SMOOTH_WIN))
    if len(idx_peaks) == 0:
        return None
    proms   = props["prominences"]
    heights = a[idx_peaks]
    order   = np.argsort(heights)[::-1]         # tallest first
    top     = order[0]
    top_h, top_prom = heights[top], proms[top]
    runner_h = heights[order[1]] if len(order) > 1 else 0.0
    if top_h <= 0:
        return None
    if top_prom / top_h < prom_ratio:           # fails LOCAL cleanliness
        return None
    if runner_h > runnerup_max * top_h:         # fails GLOBAL uniqueness
        return None
    return int(idx_peaks[top]), float(top_prom / top_h), float(top_h)

def _null_one(kin_m, C_m, method, n_shuf, seed, k, nb, discrete):
    """Per-window shuffle null for ONE neuron -> (n_wins,), the MEAN MI over
    `n_shuf` trial-permutations of kin at every window. The same permutations are
    reused across windows on purpose: each window's null stays unbiased and the
    resulting curve is smoother window-to-window than redrawing every time.
    Top-level (not a closure) so joblib can ship it to worker processes."""
    r    = np.random.default_rng(seed)
    shuf = [r.permutation(kin_m) for _ in range(n_shuf)]
    out  = np.full(C_m.shape[1], np.nan)
    for wi in range(C_m.shape[1]):
        c_w = C_m[:, wi]
        if np.unique(c_w).size < 2: continue
        out[wi] = np.mean([_bin(ks, c_w, nb, discrete) if method == "bin"
                           else _ksg(ks, c_w, k) for ks in shuf])
    return out


def _null_landscape(method, run_idx, tag):
    """Per-window null for every neuron in run_idx -> (n_run, n_wins)."""
    discrete = mi_rates.is_discrete(NEURAL_VAR)
    payload  = []
    for pi, i in enumerate(run_idx):
        mask = mask_list[i]
        if mask.sum() < MIN_TRIALS:
            payload.append(None); continue
        C_m = neural_landscape(NEURONS[i], win_starts, WIN_WIDTH)[mask]
        payload.append((kin_list[i][mask], np.asarray(C_m, float)))
    t0 = time.time()
    res = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
        delayed(_null_one)(kin_m, C_m, method, N_SHUF_NULL, 1000 + pi,
                           K_KSG, N_BINS_KIN, discrete)
        for pi, item in enumerate(payload) if item is not None
        for kin_m, C_m in [item])
    out = np.full((len(run_idx), n_wins), np.nan)
    for pi, arr in zip([pi for pi, it in enumerate(payload) if it is not None], res):
        out[pi] = arr
    print(f"{tag}: null done in {time.time()-t0:.0f}s "
          f"(S={N_SHUF_NULL}, {N_JOBS} jobs)")
    return out


def _example_set(mi_mat):
    """-> set of positions in run_idx whose individual landscape is drawn:
    N_EXAMPLES per region, chosen by EXAMPLE_MODE from the peak MI of mi_mat
    (rows aligned with run_idx). Chosen once and cached in EX_SEL, so every
    method draws the same neurons; delete EX_SEL (or restart) to choose again."""
    global EX_SEL
    if N_EXAMPLES is None:
        return set(range(len(run_idx)))
    if "EX_SEL" in globals() and EX_SEL.get("key") == (N_EXAMPLES, EXAMPLE_MODE,
                                                      len(run_idx)):
        return EX_SEL["sel"]
    rng, sel = np.random.default_rng(0), set()
    for reg in ("CP", "CB"):
        idx = [pi for pi, i in enumerate(run_idx)
               if regions[i] == reg and np.isfinite(mi_mat[pi]).any()]
        if not idx: continue
        pk    = np.array([np.nanmax(mi_mat[pi]) for pi in idx])
        order = [idx[j] for j in np.argsort(-pk)]          # strongest first
        k = min(N_EXAMPLES, len(order))
        if EXAMPLE_MODE == "top":      pick = order[:k]
        elif EXAMPLE_MODE == "random": pick = list(rng.choice(order, k, replace=False))
        elif EXAMPLE_MODE == "spread": pick = [order[int(round(t))]
                                               for t in np.linspace(0, len(order)-1, k)]
        else: raise ValueError(f"EXAMPLE_MODE must be spread|top|random, "
                               f"got {EXAMPLE_MODE!r}")
        sel.update(int(p) for p in pick)
        print(f"  individual landscapes, {reg}: {len(set(pick))} of {len(idx)} "
              f"({EXAMPLE_MODE})")
    EX_SEL = dict(key=(N_EXAMPLES, EXAMPLE_MODE, len(run_idx)), sel=sel)
    return sel


def _plot_individual(uid, reg, mi_arr, win_centers, method_label, save_dir,
                     sig=None, entropy=None, null_arr=None):
    fig, ax = plt.subplots(figsize=(7, 2.8))
    ax.plot(win_centers, mi_arr, color=CLR[reg], lw=0.8, alpha=0.35, label="raw")

    mi_sm = _smooth(mi_arr, SMOOTH_WIN)
    ax.plot(win_centers, mi_sm, color=CLR[reg], lw=1.6, label=f"smoothed (win={SMOOTH_WIN})")

    pk = _find_peak(mi_sm)
    if pk is not None:
        pk_idx, prom_ratio, height = pk
        pk_t = win_centers[pk_idx]
        ax.axvline(pk_t, color=CLR[reg], lw=0.8, ls="--", alpha=0.7,
                   label=f"peak (prom/h={prom_ratio:.2f})")
    if null_arr is not None and np.any(~np.isnan(null_arr)):
        ax.plot(win_centers, null_arr, color="0.45", lw=1.2, ls="--", alpha=0.9,
                label=f"shuffle null (S={N_SHUF_NULL})")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.axhline(0, color="k", lw=0.4)
    sig_str  = ("  * sig" if sig else "  ns") if sig is not None else ""
    ent_str  = f"  H={entropy:.2f} bits" if (entropy is not None and np.isfinite(entropy)) else ""
    peak_str = f"  peak prom/h={prom_ratio:.2f}" if pk is not None else "  no dominant peak"
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

def _plot_mean_landscape(mi_mat, regions_run, win_centers, method_label, save_dir,
                         null_mat=None):
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

        if null_mat is not None:
            nm   = null_mat[idx][ok]
            nmn  = np.nanmean(nm, axis=0)
            nsem = np.nanstd(nm, axis=0) / np.sqrt(
                np.maximum(np.sum(~np.isnan(nm), axis=0), 1))
            ax.plot(win_centers, nmn, color=CLR[reg], lw=1.4, ls="--", alpha=0.55,
                    label=f"{reg} shuffle null")
            ax.fill_between(win_centers, nmn - nsem, nmn + nsem,
                            color=CLR[reg], alpha=0.08)

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
            pk = _find_peak(row_sm_masked)          # dominant, interior local max only
            if pk is not None:
                peaks.append(win_centers[pk[0]])
        if peaks:
            ax.hist(peaks, bins=bin_edges, color=CLR[reg], alpha=0.6,
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

    C   = neural_landscape(n, win_starts, WIN_WIDTH)
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
               COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP),
               NEURAL_VAR=NEURAL_VAR)
with open(config.var_path(f"{RUN_NAME}_bin.pkl"), "wb") as f:
    pickle.dump(res_bin, f, protocol=4)
print(f"Saved → {config.var_path(f'{RUN_NAME}_bin.pkl')}")


# %% ── METHOD 1: BINNING — lead-lag significance (shuffle at peak window) ─────
# Shuffle kin at each neuron's peak window; only peaks > p95 null enter lead-lag.
print(f"Lead-lag significance: {N_SHUF_LEADLAG} shuffles per neuron (binning)")
sig_bin = np.zeros(len(run_idx), dtype=bool)
t0 = time.time()

for pi, i in enumerate(run_idx):
    mi_arr = mi_bin[pi]
    if np.all(np.isnan(mi_arr)): continue
    mi_sm = _smooth(mi_arr, SMOOTH_WIN)             # peak window from smoothed curve
    pk = _find_peak(mi_sm)                          # dominant, interior local max only
    if pk is None: continue                          # no real peak (flat/multi-modal/edge) -> skip
    pk_wi   = pk[0]
    real_mi = mi_arr[pk_wi]                         # significance uses the actual (raw) MI at that window
    if real_mi <= 0: continue

    n    = NEURONS[i]
    mask = mask_list[i]
    if mask.sum() < MIN_TRIALS: continue
    kin_m = kin_list[i][mask]

    c_pk = neural_window(n, (float(win_starts[pk_wi]),
                             float(win_starts[pk_wi] + WIN_WIDTH)))[mask]
    if np.unique(c_pk).size < 2: continue

    null = np.array([_bin(rng.permutation(kin_m), c_pk, N_BINS_KIN)
                     for _ in range(N_SHUF_LEADLAG)])
    sig_bin[pi] = real_mi > np.percentile(null, 95)

print(f"  {sig_bin.sum()}/{len(run_idx)} significant  ({time.time()-t0:.1f}s)")

# keep the significance flags in the pickle so downstream scripts (08) can reuse
# them instead of recomputing the shuffle test.
res_bin["sig"] = sig_bin
res_bin["N_SHUF_LEADLAG"] = N_SHUF_LEADLAG
with open(config.var_path(f"{RUN_NAME}_bin.pkl"), "wb") as f:
    pickle.dump(res_bin, f, protocol=4)
print(f"Re-saved with sig -> {config.var_path(f'{RUN_NAME}_bin.pkl')}")


# %% ── METHOD 1: BINNING — per-window shuffle null ───────────────────────────
# Null MI at EVERY window (not just the peak): kin is permuted N_SHUF_NULL times
# and the mean MI is taken per window. Drawn as the grey dashed reference line on
# the individual plots and as mean +- SEM on the mean landscape.
null_bin = np.full((len(run_idx), n_wins), np.nan)
if "bin" in NULL_METHODS:
    null_bin = _null_landscape("bin", run_idx, "Binning")
    res_bin["null"] = null_bin
    res_bin["N_SHUF_NULL"] = N_SHUF_NULL
    with open(config.var_path(f"{RUN_NAME}_bin.pkl"), "wb") as f:
        pickle.dump(res_bin, f, protocol=4)
    print(f"Re-saved with null -> {config.var_path(f'{RUN_NAME}_bin.pkl')}")
else:
    print("Binning null skipped (not in NULL_METHODS)")


# %% ── METHOD 1: BINNING — plots individual landscapes ─────────────────────────
ind_dir = config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_bin", "landscapes_individual")
_ex = _example_set(mi_bin)
for pi, i in enumerate(run_idx):
    if pi not in _ex: continue
    _plot_individual(uids[i], regions[i], mi_bin[pi], win_centers, "Binning", ind_dir,
                     sig=sig_bin[pi], null_arr=null_bin[pi])
# %% ── METHOD 1: BINNING — plots mean / lead-lag ───────────────────────────────

_plot_mean_landscape(mi_bin, regions_run, win_centers, "Binning",
                     config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_bin"), null_mat=null_bin)
_plot_lead_lag(mi_bin, regions_run, win_centers, "Binning",
               config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_bin"), sig=sig_bin, peak_range=PEAK_RANGE)
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

    C   = neural_landscape(n, win_starts, WIN_WIDTH)
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
               COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), K=K_KSG,
               NEURAL_VAR=NEURAL_VAR)
with open(config.var_path(f"{RUN_NAME}_ksg.pkl"), "wb") as f:
    pickle.dump(res_ksg, f, protocol=4)
print(f"Saved → {config.var_path(f'{RUN_NAME}_ksg.pkl')}")


# %% ── METHOD 2: KSG — lead-lag significance (shuffle at peak window) ─────────
print(f"Lead-lag significance: {N_SHUF_LEADLAG} shuffles per neuron (KSG)")
sig_ksg = np.zeros(len(run_idx), dtype=bool)
t0 = time.time()

for pi, i in enumerate(run_idx):
    mi_arr = mi_ksg[pi]
    if np.all(np.isnan(mi_arr)): continue
    mi_sm = _smooth(mi_arr, SMOOTH_WIN)             # peak window from smoothed curve
    pk = _find_peak(mi_sm)                          # dominant, interior local max only
    if pk is None: continue                          # no real peak (flat/multi-modal/edge) -> skip
    pk_wi   = pk[0]
    real_mi = mi_arr[pk_wi]                         # significance uses the actual (raw) MI at that window
    if real_mi <= 0: continue

    n    = NEURONS[i]
    mask = mask_list[i]
    if mask.sum() < MIN_TRIALS: continue
    kin_m = kin_list[i][mask]

    c_pk = neural_window(n, (float(win_starts[pk_wi]),
                             float(win_starts[pk_wi] + WIN_WIDTH)))[mask]
    if np.unique(c_pk).size < 2: continue

    null = np.array([_ksg(rng.permutation(kin_m), c_pk, K_KSG)
                     for _ in range(N_SHUF_LEADLAG)])
    sig_ksg[pi] = real_mi > np.percentile(null, 95)

print(f"  {sig_ksg.sum()}/{len(run_idx)} significant  ({time.time()-t0:.1f}s)")

# keep the significance flags in the pickle so downstream scripts (08) can reuse
# them instead of recomputing the shuffle test.
res_ksg["sig"] = sig_ksg
res_ksg["N_SHUF_LEADLAG"] = N_SHUF_LEADLAG
with open(config.var_path(f"{RUN_NAME}_ksg.pkl"), "wb") as f:
    pickle.dump(res_ksg, f, protocol=4)
print(f"Re-saved with sig -> {config.var_path(f'{RUN_NAME}_ksg.pkl')}")


# %% ── METHOD 2: KSG — per-window shuffle null ───────────────────────────────
null_ksg = np.full((len(run_idx), n_wins), np.nan)
if "ksg" in NULL_METHODS:
    null_ksg = _null_landscape("ksg", run_idx, f"KSG k={K_KSG}")
    res_ksg["null"] = null_ksg
    res_ksg["N_SHUF_NULL"] = N_SHUF_NULL
    with open(config.var_path(f"{RUN_NAME}_ksg.pkl"), "wb") as f:
        pickle.dump(res_ksg, f, protocol=4)
    print(f"Re-saved with null -> {config.var_path(f'{RUN_NAME}_ksg.pkl')}")
else:
    print("KSG null skipped (not in NULL_METHODS)")


# %% ── METHOD 2: KSG — MI landscape entropy (diagnostic only, NOT used for selection) ──
# Entropy of the MI-vs-time curve (treating the positive MI as a mass distribution).
# Low entropy = peaked landscape (clear preferred time).
# High entropy = flat / noisy landscape (no clear peak).
# Kept for visualization/diagnostics only — selection of significant peaks for the
# lead-lag histogram uses the same strategy as Binning/Ross: sig_ksg alone, which
# already requires a dominant interior peak (_find_peak, PEAK_PROM_RATIO) AND that its
# real MI clears the shuffle null. Entropy didn't separate "clean peak" vs "noisy
# multi-modal" well enough to use as a cutoff — no entropy-based filtering here,
# this cell is purely descriptive.

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
ax.set(xlabel="MI landscape entropy (bits)", ylabel="# neurons",
       title=f"KSG k={K_KSG}  —  MI curve entropy per neuron\n"
             f"low = peaked  |  high = flat/noisy  |  max = {np.log2(n_wins):.1f} bits ({n_wins} windows)",
             xlim=(3, 6.5))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend()
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_ksg") / "entropy_dist.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

print(f"{'Region':6s}  {'n':>5s}  {'p10':>5s}  {'p25':>5s}  {'p50':>5s}  {'p75':>5s}  {'p90':>5s}  bits")
for reg in ["CP", "CB"]:
    idx = [i for i, r in enumerate(regions_run) if r == reg]
    e   = ent_ksg[idx]; e = e[~np.isnan(e)]
    qs  = np.percentile(e, [10, 25, 50, 75, 90])
    print(f"{reg:6s}  {len(e):5d}  {qs[0]:5.2f}  {qs[1]:5.2f}  {qs[2]:5.2f}  {qs[3]:5.2f}  {qs[4]:5.2f}")


# %% ── METHOD 2: KSG — plots ─────────────────────────────────────────────────
ind_dir = config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_ksg", "landscapes_individual")
_ex = _example_set(mi_ksg)
for pi, i in enumerate(run_idx):
    if pi not in _ex: continue
    _plot_individual(uids[i], regions[i], mi_ksg[pi], win_centers, f"KSG k={K_KSG}", ind_dir,
                     sig=sig_ksg[pi], entropy=ent_ksg[pi], null_arr=null_ksg[pi])

# %% ── METHOD 2: KSG — plots mean / lead-lag ──────────────────────────────────
_plot_mean_landscape(mi_ksg, regions_run, win_centers, f"KSG k={K_KSG}",
                     config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_ksg"), null_mat=null_ksg)
_plot_lead_lag(mi_ksg, regions_run, win_centers, f"KSG k={K_KSG}",
               config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_ksg"), sig=sig_ksg, peak_range=PEAK_RANGE)
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
if not mi_rates.is_discrete(NEURAL_VAR):
    raise RuntimeError(
        f"Ross needs a DISCRETE neural variable, but NEURAL_VAR={NEURAL_VAR!r} is a "
        f"continuous rate. Run the Binning/KSG cells for this mode, or set "
        f"NEURAL_VAR='count' to use Ross.")
cp_run_idx   = [pi for pi, i in enumerate(run_idx) if regions[i] == "CP"]
mi_ross      = np.full((len(run_idx), n_wins), np.nan)   # same shape, CB rows stay nan
t0 = time.time()

for count_cp, pi in enumerate(cp_run_idx):
    i     = run_idx[pi]
    n     = NEURONS[i]
    mask  = mask_list[i]
    kin_m = kin_list[i][mask]
    if mask.sum() < MIN_TRIALS: continue

    C   = neural_landscape(n, win_starts, WIN_WIDTH)
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
                COORD=COORD, FRAMES=FRAMES, CAT_KEEP=list(CAT_KEEP), K=K_KSG,
               NEURAL_VAR=NEURAL_VAR)
with open(config.var_path(f"{RUN_NAME}_ross.pkl"), "wb") as f:
    pickle.dump(res_ross, f, protocol=4)
print(f"Saved → {config.var_path(f'{RUN_NAME}_ross.pkl')}")


# %% ── METHOD 3: ROSS — lead-lag significance (CP only) ──────────────────────
print(f"Lead-lag significance: {N_SHUF_LEADLAG} shuffles per neuron (Ross, CP only)")
sig_ross = np.zeros(len(run_idx), dtype=bool)
t0 = time.time()

for count_cp, pi in enumerate(cp_run_idx):
    mi_arr = mi_ross[pi]
    if np.all(np.isnan(mi_arr)): continue
    mi_sm = _smooth(mi_arr, SMOOTH_WIN)             # peak window from smoothed curve
    pk = _find_peak(mi_sm)                          # dominant, interior local max only
    if pk is None: continue                          # no real peak (flat/multi-modal/edge) -> skip
    pk_wi   = pk[0]
    real_mi = mi_arr[pk_wi]                         # significance uses the actual (raw) MI at that window
    if real_mi <= 0: continue

    i    = run_idx[pi]
    n    = NEURONS[i]
    mask = mask_list[i]
    if mask.sum() < MIN_TRIALS: continue
    kin_m = kin_list[i][mask]

    c_pk = neural_window(n, (float(win_starts[pk_wi]),
                             float(win_starts[pk_wi] + WIN_WIDTH)))[mask]
    if np.unique(c_pk).size < 2: continue

    null = np.array([_ross(rng.permutation(kin_m), c_pk, K_KSG)
                     for _ in range(N_SHUF_LEADLAG)])
    sig_ross[pi] = real_mi > np.percentile(null, 95)

print(f"  {sig_ross[cp_run_idx].sum()}/{len(cp_run_idx)} CP significant  ({time.time()-t0:.1f}s)")

# keep the significance flags in the pickle so downstream scripts (08) can reuse
# them instead of recomputing the shuffle test.
res_ross["sig"] = sig_ross
res_ross["N_SHUF_LEADLAG"] = N_SHUF_LEADLAG
with open(config.var_path(f"{RUN_NAME}_ross.pkl"), "wb") as f:
    pickle.dump(res_ross, f, protocol=4)
print(f"Re-saved with sig -> {config.var_path(f'{RUN_NAME}_ross.pkl')}")


# %% ── METHOD 3: ROSS — plots + KSG vs Ross difference (CP) ──────────────────
ind_dir = config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_ross_cp", "landscapes_individual")
_ex = _example_set(mi_ross)
for pi in cp_run_idx:
    if pi not in _ex: continue
    _plot_individual(uids[run_idx[pi]], "CP", mi_ross[pi], win_centers,
                     f"Ross k={K_KSG} (CP)", ind_dir, sig=sig_ross[pi])
# %%
_plot_mean_landscape(mi_ross[cp_run_idx], ["CP"]*len(cp_run_idx),
                     win_centers, f"Ross k={K_KSG} (CP only)",
                     config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_ross_cp"))
_plot_lead_lag(mi_ross[cp_run_idx], ["CP"]*len(cp_run_idx),
               win_centers, f"Ross k={K_KSG} (CP only)",
               config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_ross_cp"), sig=sig_ross[cp_run_idx], peak_range=PEAK_RANGE)

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
    fig.savefig(config.fig_dir(OUT_SUBDIR, f"{RUN_NAME}_ross_cp") / "diff_ksg_ross_CP.png",
                dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)
print("Ross plots done.")

# %%
