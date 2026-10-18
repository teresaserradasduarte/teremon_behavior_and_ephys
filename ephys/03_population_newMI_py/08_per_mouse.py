# %% ── 0. imports & constants ─────────────────────────────────────────────────
# Per-MOUSE version of the 04 fixed-window summary: mean landscape + lead/lag,
# KSG only, spike COUNT by default. No per-neuron figures.
#
# The MI landscape is NOT recomputed if 04 already saved it: this script loads
# {RUN_NAME}_ksg.pkl and only computes what is missing. Keep the window/estimator
# parameters below identical to the 04 run that produced that pickle -- they are
# checked against its metadata and a mismatch forces a recompute.
import numpy as np
import matplotlib.pyplot as plt
import pickle, time, importlib
import mi_io, config, mi_rates
importlib.reload(config); importlib.reload(mi_rates)
from sklearn.feature_selection import mutual_info_regression
from scipy.signal import find_peaks

if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

CLR     = {"CP": np.array([106,176,134])/255, "CB": np.array([154,123,193])/255}
NAT2BIT = 1 / np.log(2)
rng     = np.random.default_rng(0)

# ── must match the 04 run whose pickle you want to reuse ──────────────────────
# (defaults point at the 100 ms run: RUN_NAME=fixedwin100, WIN_WIDTH=0.100)
WIN_WIDTH   = 0.100       # sliding window width (s)
WIN_STEP    = 0.025       # step size (s)
WIN_RANGE   = (-1.0, 1.0) # landscape range (s), re reach
K_KSG       = 3           # k for KSG
COORD       = 1           # 0=x, 1=y med-lat, 2=z
FRAMES      = (147, 167)  # endpoint window, MATLAB 1-indexed
CAT_KEEP    = {1, 2}      # reach categories
MIN_TRIALS  = 20

NEURAL_VAR  = "count"     # "count" (default) | "isi" | "conv"  -- see mi_rates.py.
                          # Only KSG is used here, so any of the three is valid
                          # (unlike Ross, which needs a discrete variable).

# ── peak / significance settings (same meaning as in 04) ──────────────────────
N_SHUF_LEADLAG    = 50
SMOOTH_WIN        = 2
PEAK_PROM_RATIO   = 0.5
PEAK_RUNNERUP_MAX = 0.7
PEAK_RANGE        = None  # e.g. (-0.3, 0.3); None = full WIN_RANGE

MIN_NEURONS_MOUSE = 20    # skip a mouse/region with fewer than this many neurons --
                          # a "mean landscape" over a handful of cells is noise.

COMPUTE_IF_MISSING = True # if no matching {RUN_NAME}_ksg.pkl exists (or its parameters
                          # differ), compute the landscape + significance HERE.
                          # ~10 min for the full population with KSG. Set False to stop
                          # with a message instead, so a long run never starts by accident.
SAVE_COMPUTED      = True # write what was computed to {RUN_NAME}_ksg.pkl so the next run
                          # just loads it. NOTE: 08 does not compute the per-window null
                          # (that is 04's job, ~12 min) -- a pickle written here therefore
                          # has no "null" key and the mean plots come without the dashed
                          # null until 04 has been run for this parameter set.
SHOW_FIGS   = True
SAVE_FIGS   = True

RUN_NAME    = "fixedwin100"     # which 04 pickle to reuse ({RUN_NAME}_ksg.pkl)
OUT_SUBDIR  = "per_mouse_100ms" # figures folder under OUT_DIR

win_starts  = np.arange(WIN_RANGE[0], WIN_RANGE[1] - WIN_WIDTH + WIN_STEP/2, WIN_STEP)
win_centers = win_starts + WIN_WIDTH / 2
n_wins      = len(win_starts)
print(f"windows: {n_wins}  width={WIN_WIDTH*1000:.0f} ms  step={WIN_STEP*1000:.0f} ms  "
      f"neural variable: {NEURAL_VAR}")


# %% ── 1. helpers (same definitions as 04, kept local per project style) ──────
def _ksg(kin_m, c_w, k):
    X = kin_m.reshape(-1, 1)
    return mutual_info_regression(X, c_w + rng.normal(0, 1e-9, c_w.size),
                                  n_neighbors=k, random_state=0)[0] * NAT2BIT

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
    """Tallest INTERIOR local max that is both locally clean and globally unique."""
    valid = ~np.isnan(arr)
    if valid.sum() < 3: return None
    a = np.where(valid, arr, -np.inf)
    idx_peaks, props = find_peaks(a, prominence=0, distance=max(1, SMOOTH_WIN))
    if len(idx_peaks) == 0: return None
    proms = props["prominences"]; heights = a[idx_peaks]
    order = np.argsort(heights)[::-1]; top = order[0]
    top_h, top_prom = heights[top], proms[top]
    runner_h = heights[order[1]] if len(order) > 1 else 0.0
    if top_h <= 0 or top_prom / top_h < prom_ratio or runner_h > runnerup_max * top_h:
        return None
    return int(idx_peaks[top]), float(top_prom / top_h), float(top_h)


# %% ── 2. load the 04 landscape, or compute what is missing ──────────────────
pkl = config.var_path(f"{RUN_NAME}_ksg.pkl")
res, loaded = None, False
if pkl.exists():
    with open(pkl, "rb") as f:
        res = pickle.load(f)
    same = (np.isclose(res.get("WIN_WIDTH", np.nan), WIN_WIDTH)
            and len(res.get("win_centers", [])) == n_wins
            and np.allclose(res["win_centers"], win_centers)
            and res.get("K", K_KSG) == K_KSG
            and res.get("COORD", COORD) == COORD
            and tuple(res.get("FRAMES", FRAMES)) == tuple(FRAMES)
            and set(res.get("CAT_KEEP", CAT_KEEP)) == set(CAT_KEEP)
            # pickles written before NEURAL_VAR was recorded are assumed "count"
            and res.get("NEURAL_VAR", "count") == NEURAL_VAR)
    if same:
        loaded = True
        # NB these are presence checks (is the key in the file?), not results.
        print(f"loaded {pkl.name}: mi{res['mi'].shape}"
              f" | per-window null stored: {'yes' if 'null' in res else 'no'}"
              f" | significance flags stored: "
              f"{'yes' if 'sig' in res else 'no - will be recomputed here'}")
    else:
        print(f"{pkl.name} exists but its parameters differ from this script - recomputing")
        res = None
else:
    print(f"{pkl.name} not found - computing the landscape here")

uids    = [n["uid"] for n in NEURONS]
regions = [n["region"] for n in NEURONS]
mice    = [u.split("|")[1] for u in uids]          # uid = REG|mouse|sess|phyID

kin_list, mask_list = [], []
for n in NEURONS:
    kin  = mi_io.endpoint_kin(n, REACH_PX, coord=COORD, frames=FRAMES)
    cat  = n["labels"]["idx_reach_cat"]
    kin_list.append(kin)
    mask_list.append(np.isfinite(kin) & np.isin(cat, list(CAT_KEEP)))

if not loaded and not COMPUTE_IF_MISSING:
    raise RuntimeError(
        f"No matching landscape for these parameters and COMPUTE_IF_MISSING is False.\n"
        f"  expected: {pkl}\n"
        f"  either run 04_MI_fixed_win.py with RUN_NAME={RUN_NAME!r}, WIN_WIDTH="
        f"{WIN_WIDTH}, NEURAL_VAR={NEURAL_VAR!r},\n"
        f"  or set COMPUTE_IF_MISSING = True to compute it here (~10 min).")

if loaded:
    assert list(res["uids"]) == uids, "pickle neuron list differs from the cache"
    mi_ksg   = res["mi"]
    null_ksg = res.get("null")
else:
    print(f"computing the KSG landscape here for {len(NEURONS)} neurons...")
    mi_ksg, null_ksg = np.full((len(NEURONS), n_wins), np.nan), None
    t0 = time.time()
    for i, n in enumerate(NEURONS):
        mask = mask_list[i]
        if mask.sum() < MIN_TRIALS: continue
        kin_m = kin_list[i][mask]
        C_m = np.asarray(mi_rates.windowed(n["spikes"], win_starts, WIN_WIDTH,
                                           NEURAL_VAR), float)[mask]
        for wi in range(n_wins):
            c_w = C_m[:, wi]
            if np.unique(c_w).size < 2: continue
            mi_ksg[i, wi] = _ksg(kin_m, c_w, K_KSG)
        if (i + 1) % 250 == 0:
            el = time.time() - t0
            print(f"  {i+1}/{len(NEURONS)}  {el:.0f}s "
                  f"(~{el/(i+1)*(len(NEURONS)-i-1):.0f}s left)", flush=True)
    print(f"landscape computed in {time.time()-t0:.0f}s")


# %% ── 3. peak significance (reuse from the pickle when it is there) ─────────
if loaded and "sig" in res:
    sig_ksg = np.asarray(res["sig"], bool)
    print(f"significance reused from the pickle: {sig_ksg.sum()}/{len(sig_ksg)}")
else:
    print(f"computing peak significance ({N_SHUF_LEADLAG} shuffles/neuron)...")
    sig_ksg = np.zeros(len(NEURONS), dtype=bool)
    t0 = time.time()
    for i, n in enumerate(NEURONS):
        mi_arr = mi_ksg[i]
        if np.all(np.isnan(mi_arr)): continue
        pk = _find_peak(_smooth(mi_arr, SMOOTH_WIN))
        if pk is None: continue
        pk_wi = pk[0]; real_mi = mi_arr[pk_wi]
        if real_mi <= 0: continue
        mask = mask_list[i]
        if mask.sum() < MIN_TRIALS: continue
        kin_m = kin_list[i][mask]
        c_pk = np.asarray(mi_rates.windowed(n["spikes"], win_starts, WIN_WIDTH,
                                            NEURAL_VAR), float)[mask][:, pk_wi]
        if np.unique(c_pk).size < 2: continue
        null = np.array([_ksg(rng.permutation(kin_m), c_pk, K_KSG)
                         for _ in range(N_SHUF_LEADLAG)])
        sig_ksg[i] = real_mi > np.percentile(null, 95)
        if (i + 1) % 250 == 0:
            print(f"  {i+1}/{len(NEURONS)}  {time.time()-t0:.0f}s", flush=True)
    print(f"  {sig_ksg.sum()}/{len(NEURONS)} significant ({time.time()-t0:.0f}s)")

if not loaded and SAVE_COMPUTED:
    out = dict(method="ksg", win_centers=win_centers, win_starts=win_starts,
               WIN_WIDTH=WIN_WIDTH, WIN_STEP=WIN_STEP, mi=mi_ksg, sig=sig_ksg,
               uids=uids, regions=regions, COORD=COORD, FRAMES=FRAMES,
               CAT_KEEP=list(CAT_KEEP), K=K_KSG, NEURAL_VAR=NEURAL_VAR,
               N_SHUF_LEADLAG=N_SHUF_LEADLAG)
    with open(pkl, "wb") as f:
        pickle.dump(out, f, protocol=4)
    print(f"saved computed landscape + sig -> {pkl.name} "
          f"(no 'null' key: run 04 for the per-window null)")

MOUSE_LIST = sorted(set(mice))
print("\nneurons per mouse:")
for ms in MOUSE_LIST:
    line = "  ".join(f"{reg}={sum(1 for j in range(len(uids)) if mice[j]==ms and regions[j]==reg)}"
                     for reg in ("CP", "CB"))
    print(f"  {ms:18s} {line}")


# %% ── 4. plotting helpers (per-mouse variants of the 04 figures) ────────────
def _draw_mean_landscape(ax, sel, title, fontsize=7, xlabel=True, ylabel=True):
    """Draw the mean +- SEM MI landscape per region into `ax` (with the shuffle
    null dashed/transparent when available). Split out from the figure-making
    wrapper so the same drawing can go into a panel of a grid."""
    for reg in ("CP", "CB"):
        idx = [i for i in sel if regions[i] == reg]
        if len(idx) < MIN_NEURONS_MOUSE: continue
        m  = mi_ksg[idx]
        ok = ~np.all(np.isnan(m), axis=1)
        m  = m[ok]
        if m.shape[0] == 0: continue
        mn  = np.nanmean(m, axis=0)
        sem = np.nanstd(m, axis=0) / np.sqrt(np.maximum(np.sum(~np.isnan(m), axis=0), 1))
        ax.plot(win_centers, mn, color=CLR[reg], lw=2, label=f"{reg} (n={ok.sum()})")
        ax.fill_between(win_centers, mn - sem, mn + sem, color=CLR[reg], alpha=0.2)
        if null_ksg is not None:
            nm   = null_ksg[idx][ok]
            nmn  = np.nanmean(nm, axis=0)
            nsem = np.nanstd(nm, axis=0) / np.sqrt(
                np.maximum(np.sum(~np.isnan(nm), axis=0), 1))
            ax.plot(win_centers, nmn, color=CLR[reg], lw=1.4, ls="--", alpha=0.55,
                    label=f"{reg} shuffle null")
            ax.fill_between(win_centers, nmn - nsem, nmn + nsem,
                            color=CLR[reg], alpha=0.08)
        pk = win_centers[np.nanargmax(mn)]
        ax.axvline(pk, color=CLR[reg], lw=1.5, ls="--", alpha=0.85,
                   label=f"{reg} peak = {pk*1000:.0f} ms")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.axhline(0, color="k", lw=0.4)
    ax.set_title(title)
    if xlabel: ax.set_xlabel("window center re reach (s)")
    if ylabel: ax.set_ylabel("mean MI ± SEM (bits)")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=fontsize)


def _plot_mean_landscape(sel, title, save_path):
    fig, ax = plt.subplots(figsize=(8, 4))
    _draw_mean_landscape(ax, sel, title)
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


def _peaks_for(idx):
    """Peak times of the significant, genuinely-peaked neurons in idx."""
    lo, hi = PEAK_RANGE if PEAK_RANGE is not None else (win_centers[0], win_centers[-1])
    win_mask = (win_centers >= lo) & (win_centers <= hi)
    out = []
    for i in idx:
        if not sig_ksg[i]: continue
        row = mi_ksg[i]
        if not np.any(~np.isnan(row)): continue
        pk = _find_peak(np.where(win_mask, _smooth(row, SMOOTH_WIN), np.nan))
        if pk is not None:
            out.append(win_centers[pk[0]])
    return out


def _draw_lead_lag(ax, sel, title, fontsize=8, xlabel=True, ylabel=True):
    """Draw the significant-peak-time distribution per region into `ax`."""
    lo, hi = PEAK_RANGE if PEAK_RANGE is not None else (win_centers[0], win_centers[-1])
    bin_edges = np.linspace(lo, hi, 21)
    for reg in ("CP", "CB"):
        peaks = _peaks_for([i for i in sel if regions[i] == reg])
        if not peaks: continue
        ax.hist(peaks, bins=bin_edges, color=CLR[reg], alpha=0.6, edgecolor="white",
                label=f"{reg} (n={len(peaks)})")
        med = np.median(peaks)
        ax.axvline(med, color=CLR[reg], lw=1.5, ls="--", alpha=0.9,
                   label=f"{reg} median = {med*1000:.0f} ms")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.set_title(title)
    if xlabel: ax.set_xlabel("time of peak MI re reach (s)")
    if ylabel: ax.set_ylabel("# neurons")
    if PEAK_RANGE is not None: ax.set_xlim(PEAK_RANGE)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if ax.get_legend_handles_labels()[0]: ax.legend(fontsize=fontsize)


def _plot_lead_lag(sel, title, save_path):
    fig, ax = plt.subplots(figsize=(6, 4))
    _draw_lead_lag(ax, sel, title)
    plt.tight_layout()
    if SAVE_FIGS: fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 5. one mean landscape + one lead/lag figure per mouse ─────────────────
fdir = config.fig_dir(OUT_SUBDIR)
lab  = f"KSG k={K_KSG}, {NEURAL_VAR}, win {WIN_WIDTH*1000:.0f} ms"

# pooled reference (all mice together), then one figure per mouse
_plot_mean_landscape(range(len(NEURONS)), f"ALL MICE  —  {lab}",
                     fdir / "mean_landscape_ALL.png")
_plot_lead_lag(range(len(NEURONS)), f"ALL MICE  —  {lab}",
               fdir / "lead_lag_ALL.png")

for ms in MOUSE_LIST:
    sel = [i for i in range(len(NEURONS)) if mice[i] == ms]
    _plot_mean_landscape(sel, f"{ms}  —  {lab}",
                         fdir / f"mean_landscape_{ms}.png")
    _plot_lead_lag(sel, f"{ms}  —  {lab}",
                   fdir / f"lead_lag_{ms}.png")
print(f"per-mouse figures -> {fdir}")


# %% ── 6. across-mouse summary: is the CP/CB pattern consistent per animal? ──
fig, axes = plt.subplots(1, 2, figsize=(13, 4), sharey=True)
for ax, reg in zip(axes, ("CP", "CB")):
    for k, ms in enumerate(MOUSE_LIST):
        idx = [i for i in range(len(NEURONS)) if mice[i] == ms and regions[i] == reg]
        if len(idx) < MIN_NEURONS_MOUSE: continue
        m  = mi_ksg[idx]
        ok = ~np.all(np.isnan(m), axis=1)
        if ok.sum() == 0: continue
        mn = np.nanmean(m[ok], axis=0)
        ax.plot(win_centers, mn, lw=1.6, alpha=0.85,
                color=plt.cm.tab10(k), label=f"{ms} (n={ok.sum()})")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.set(xlabel="window center re reach (s)", title=f"{reg}  —  {lab}")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(fontsize=7)
axes[0].set_ylabel("mean MI (bits)")
plt.suptitle("Mean MI landscape per mouse", fontsize=11)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(fdir / "mean_landscape_by_mouse.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

# numeric summary
print(f"\n{'mouse':18s} {'reg':4s} {'n':>5s} {'n_sig':>6s} {'peak_ms':>8s} {'median_ms':>10s}")
for ms in MOUSE_LIST + ["ALL"]:
    for reg in ("CP", "CB"):
        idx = [i for i in range(len(NEURONS))
               if regions[i] == reg and (ms == "ALL" or mice[i] == ms)]
        if not idx: continue
        m  = mi_ksg[idx]; ok = ~np.all(np.isnan(m), axis=1)
        pk = (f"{win_centers[np.nanargmax(np.nanmean(m[ok], axis=0))]*1000:.0f}"
              if ok.sum() else "-")
        peaks = _peaks_for(idx)
        med = f"{np.median(peaks)*1000:.0f}" if peaks else "-"
        print(f"{ms:18s} {reg:4s} {ok.sum():5d} {len(peaks):6d} {pk:>8s} {med:>10s}")


# %% ── 7. six-panel overviews: one panel per mouse + the pooled mean ─────────
# Same two figures as above, but with every mouse side by side and the pooled
# "ALL MICE" panel last, so per-animal consistency is readable at a glance.
PANELS = [(ms, [i for i in range(len(NEURONS)) if mice[i] == ms]) for ms in MOUSE_LIST]
PANELS.append(("ALL MICE", list(range(len(NEURONS)))))
NCOL = 3
NROW = int(np.ceil(len(PANELS) / NCOL))

# --- mean landscape grid (shared y: the panels are meant to be compared) ---
fig, axes = plt.subplots(NROW, NCOL, figsize=(6*NCOL, 3.6*NROW),
                         sharex=True, sharey=True)
axes = np.atleast_1d(axes).ravel()
for k, (name, sel) in enumerate(PANELS):
    _draw_mean_landscape(axes[k], sel, name, fontsize=6,
                         xlabel=(k // NCOL == NROW - 1), ylabel=(k % NCOL == 0))
for ax in axes[len(PANELS):]:
    ax.set_visible(False)
plt.suptitle(f"Mean MI landscape per mouse  —  {lab}", fontsize=12)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(fdir / "panels_mean_landscape.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

# --- lead/lag grid (y NOT shared: neuron counts differ a lot between mice,
#     and the pooled panel would otherwise flatten every individual one) ---
fig, axes = plt.subplots(NROW, NCOL, figsize=(6*NCOL, 3.6*NROW), sharex=True)
axes = np.atleast_1d(axes).ravel()
for k, (name, sel) in enumerate(PANELS):
    _draw_lead_lag(axes[k], sel, name, fontsize=7,
                   xlabel=(k // NCOL == NROW - 1), ylabel=(k % NCOL == 0))
for ax in axes[len(PANELS):]:
    ax.set_visible(False)
plt.suptitle(f"Peak-MI (lead/lag) distribution per mouse  —  {lab}", fontsize=12)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(fdir / "panels_lead_lag.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)
print(f"panel figures -> {fdir}")

# %%
