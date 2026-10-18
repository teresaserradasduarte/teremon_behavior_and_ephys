# %% ── 0. parameters ─────────────────────────────────────────────────────────
# STANDALONE lead/lag figure. Reads ONLY the pickle written by 09 cell 7, so it
# needs no spike times, no kinematics and no recompute -- open it, set the labels,
# run it. Nothing here touches the neural "count" pipeline.
#
# Everything you would want to restyle lives in this cell.
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import mannwhitneyu
import config

# ---- which run to plot ------------------------------------------------------
PKL_NAME   = "multiD100_x-y-z-speed-dist_water_ksg_multiD.pkl"
OUT_SUBDIR = "multiD_kin_x-y-z-speed-dist_water"   # figures land here
FIG_NAME   = "lead_lag_joint_probab"                 # extension added per format

# ---- text: title, axes, legend ----------------------------------------------
# Set any of these to None to fall back to the automatic version.
TITLE       = None
SUBTITLE    = None          # small grey line under the title; None = the run label
XLABEL      = "Time of peak MI relative to reach (ms)"
YLABEL      = None          # None -> matches HIST_NORM, e.g. "Proportion of neurons"
REG_LABEL   = {"CP": "Basal Ganglia", "CB": "Cerebellum"}
SHOW_N      = True          # append "(n=123)" to each legend entry
SHOW_MEDIAN = True          # dashed vertical line + legend entry per region
ZERO_LABEL  = "Reach endpoint"      # "" or None to drop the line

# ---- histogram --------------------------------------------------------------
# "probability" -> each region's bars sum to 1. This is the one that makes CP and
#                  CB comparable when they have different neuron counts.
# "density"     -> area sums to 1 (units are 1/time). NOT the same thing.
# "count"       -> raw neuron counts, as in 09.
HIST_NORM  = "probability"
N_BINS     = 30
PEAK_RANGE = None           # (lo, hi) in SECONDS to restrict peaks; None = all windows
TIME_UNIT  = "ms"           # "ms" or "s" -- axis scaling only
XLIM       = (-750, 750)    # (lo, hi) in TIME_UNIT, e.g. (-500, 500); None = auto
HIST_ALPHA = 0.5
HIST_STYLE = "bars"         # "bars" or "step"

# ---- peak detection (keep identical to 09 unless you mean to change it) -----
SMOOTH_WIN        = 3
PEAK_PROM_RATIO   = 0.5
PEAK_RUNNERUP_MAX = 0.7

# ---- output -----------------------------------------------------------------
FIGSIZE  = (6, 4)
DPI      = 300
SAVE_PNG = True
SAVE_PDF = True
SHOW_FIG = True
# Keep PDF text as real text (not outlines) so it stays editable in Illustrator.
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"]  = 42

CLR = {"CP": np.array([106, 176, 134]) / 255,
       "CB": np.array([154, 123, 193]) / 255}

#CLR = {"CP": np.array([215, 159, 177]) / 255,
#       "CB": np.array([220, 200, 158]) / 255}
REGIONS = ("CP", "CB")


# %% ── 1. load ───────────────────────────────────────────────────────────────
PKL = config.var_path(PKL_NAME)
if not os.path.exists(PKL):
    raise FileNotFoundError(
        f"{PKL}\nRun 09 with the matching FEATURES first (it writes this on cell 7).")
with open(PKL, "rb") as f:
    D = pickle.load(f)

MI          = D["mi"]                 # (n_neurons, n_windows)
SIG         = D["sig"]                # bool, per neuron
REGION_OF   = list(D["regions"])
win_centers = np.asarray(D["win_centers"], float)
FEATURES    = list(D["features"])
K_KSG       = D["K"]
NEURAL_VAR  = D["NEURAL_VAR"]

RUN_LAB = f"KSG k={K_KSG}, {NEURAL_VAR}, joint [{', '.join(FEATURES)}]"
print(f"loaded {PKL_NAME}")
print(f"  {MI.shape[0]} neurons x {MI.shape[1]} windows | {RUN_LAB}")
print(f"  significant: {int(SIG.sum())}")


# %% ── 2. helpers (verbatim from 09, so these are the same peaks) ────────────
def _smooth(arr, win):
    if win <= 1: return arr.copy()
    n = len(arr); out = np.full(n, np.nan); half = win // 2
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        seg = arr[lo:hi]
        if np.any(~np.isnan(seg)): out[i] = np.nanmean(seg)
    return out


def _find_peak(arr, prom_ratio=PEAK_PROM_RATIO, runnerup_max=PEAK_RUNNERUP_MAX,
               smooth_win=SMOOTH_WIN):
    valid = ~np.isnan(arr)
    if valid.sum() < 3: return None
    a = np.where(valid, arr, -np.inf)
    idx_peaks, props = find_peaks(a, prominence=0, distance=max(1, smooth_win))
    if len(idx_peaks) == 0: return None
    proms = props["prominences"]; heights = a[idx_peaks]
    order = np.argsort(heights)[::-1]; top = order[0]
    top_h, top_prom = heights[top], proms[top]
    runner_h = heights[order[1]] if len(order) > 1 else 0.0
    if top_h <= 0 or top_prom / top_h < prom_ratio or runner_h > runnerup_max * top_h:
        return None
    return int(idx_peaks[top]), float(top_prom / top_h), float(top_h)


# %% ── 3. peak times per region ──────────────────────────────────────────────
lo, hi = PEAK_RANGE if PEAK_RANGE is not None else (win_centers[0], win_centers[-1])
win_mask = (win_centers >= lo) & (win_centers <= hi)

PEAKS = {}
for reg in REGIONS:
    v = []
    for pi, r in enumerate(REGION_OF):
        if r != reg or not SIG[pi]: continue
        row = MI[pi]
        if not np.any(~np.isnan(row)): continue
        pk = _find_peak(np.where(win_mask, _smooth(row, SMOOTH_WIN), np.nan))
        if pk is not None: v.append(win_centers[pk[0]])
    PEAKS[reg] = np.asarray(v, float)

print(f"\npeak times in [{lo:.3f}, {hi:.3f}] s")
for reg in REGIONS:
    p = PEAKS[reg]
    if p.size == 0:
        print(f"  {reg}: no peaks"); continue
    q1, q3 = np.percentile(p, [25, 75])
    print(f"  {reg}: n={p.size:4d}  median={np.median(p)*1000:+7.1f} ms  "
          f"IQR=[{q1*1000:+.1f}, {q3*1000:+.1f}]  mean={p.mean()*1000:+.1f}")
if all(PEAKS[r].size > 1 for r in REGIONS):
    u, pv = mannwhitneyu(PEAKS["CP"], PEAKS["CB"])
    d = (np.median(PEAKS["CB"]) - np.median(PEAKS["CP"])) * 1000
    print(f"  CB - CP = {d:+.1f} ms   Mann-Whitney p = {pv:.2g}")


# %% ── 4. the figure ─────────────────────────────────────────────────────────
SC    = 1000.0 if TIME_UNIT == "ms" else 1.0
edges = np.linspace(lo, hi, N_BINS + 1) * SC

_auto_y = {"probability": "Proportion of neurons",
           "density":     "Probability density",
           "count":       "# neurons"}[HIST_NORM]

fig, ax = plt.subplots(figsize=FIGSIZE)
for reg in REGIONS:
    p = PEAKS[reg] * SC
    if p.size == 0: continue
    # "probability" = fraction of THIS region's neurons per bin, so the two
    # regions stay comparable even with very different n.
    w = np.full(p.size, 1.0 / p.size) if HIST_NORM == "probability" else None
    lbl = REG_LABEL.get(reg, reg) + (f" (n={p.size})" if SHOW_N else "")
    if HIST_STYLE == "step":
        ax.hist(p, bins=edges, weights=w, density=(HIST_NORM == "density"),
                histtype="step", lw=2, color=CLR[reg], label=lbl)
    else:
        ax.hist(p, bins=edges, weights=w, density=(HIST_NORM == "density"),
                color=CLR[reg], alpha=HIST_ALPHA, edgecolor="white", label=lbl)
    if SHOW_MEDIAN:
        med = np.median(p)
        ax.axvline(med, color=CLR[reg], lw=1.5, ls="--", alpha=0.9,
                   label=f"{reg} median = {med:.0f} {TIME_UNIT}")

if ZERO_LABEL:
    ax.axvline(0, color="k", lw=0.8, ls=":", label=ZERO_LABEL)
ax.set(xlabel=(XLABEL if XLABEL is not None
               else f"time of peak MI re reach endpoint ({TIME_UNIT})"),
       ylabel=YLABEL if YLABEL is not None else _auto_y,
       )
_sub = RUN_LAB if SUBTITLE is None else SUBTITLE
# Pad the title to leave room for the subtitle line, otherwise the two overlap.
ax.set_title(TITLE if TITLE is not None else f"{RUN_LAB}  -  lead / lag",
             pad=16 if _sub else 6)
if _sub:
    ax.text(0.5, 1.005, _sub, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=7, color="0.35")
if XLIM is not None: ax.set_xlim(XLIM)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.legend(fontsize=8, frameon=False)
plt.tight_layout()

fdir  = config.fig_dir(OUT_SUBDIR)
saved = []
if SAVE_PNG:
    q = fdir / f"{FIG_NAME}.png"
    fig.savefig(q, dpi=DPI, bbox_inches="tight"); saved.append(q)
if SAVE_PDF:
    q = fdir / f"{FIG_NAME}.pdf"
    fig.savefig(q, bbox_inches="tight"); saved.append(q)
for q in saved: print(f"saved -> {q}")
if SHOW_FIG: plt.show()
plt.close(fig)

# %%
