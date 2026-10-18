# %% ── 0. imports & parameters ────────────────────────────────────────────────
# Per-neuron MI curves from script 04, REPLOTTED with the raw data beside them.
#
# 04 computes MI( endpoint y ; spike count ) for every sliding window and saves
# one curve per neuron. It never recomputes here: this script LOADS those
# pickles. What it adds is the data behind the curve, on the same row:
#
#   col 1  the MI curve      exactly as 04 draws it (raw + smoothed + peak +
#                            shuffle null), one line per method if several
#   col 2  the raster        trials sorted by endpoint y, the scatter window
#                            shaded -- so the curve can be read against spikes
#   col 3  the SCATTER       endpoint y (px) vs spike count in one window,
#                            which is the single number the MI is computed from
#
# Pick neurons either as N representatives spread over the population, or by
# name: PICK = [("1_CoteDor", "R2", 322), ...]. In the interactive window,
#   show("1_CoteDor", "R2", 322)
# draws one neuron on demand.
import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

# The VS Code interactive window does not always start in this folder.
_HERE = (os.path.dirname(os.path.abspath(__file__))
         if "__file__" in globals() else os.getcwd())
if not os.path.isfile(os.path.join(_HERE, "mi_io.py")):
    _sub = os.path.join("ephys", "03_population_newMI_py")
    _p, _found = os.getcwd(), None
    while _found is None:
        for _c in (_p, os.path.join(_p, _sub)):
            if os.path.isfile(os.path.join(_c, "mi_io.py")): _found = _c; break
        if _found is None:
            if _p == os.path.dirname(_p):
                raise ModuleNotFoundError(
                    f"mi_io.py not found from the current directory ({os.getcwd()}). "
                    "Open this script's folder in VS Code.")
            _p = os.path.dirname(_p)
    _HERE = _found
if _HERE not in sys.path: sys.path.insert(0, _HERE)
os.chdir(_HERE)
import mi_io, mi_rates, config

# ── which saved run to read ──────────────────────────────────────────────────
MI_RUN     = "fixedwin100_v2"       # 04's RUN_NAME: fixedwin50 | fixedwin70 | fixedwin100
MI_METHODS = ("ksg",)           # curves to draw: any of "ksg", "bin", "ross"
                                # (04 saves one pickle per method; bin and ksg
                                #  also carry the shuffle null)

# ── the panels ───────────────────────────────────────────────────────────────
PANELS      = ("mi", "raster", "scatter")    # order of the columns
SCATTER_WIN = "peak"            # which window the counts come from:
                                #   "peak"  -> that neuron's peak of the smoothed
                                #              curve of MI_METHODS[0]
                                #   (lo,hi) -> a fixed window in seconds, e.g.
                                #              (0.0, 0.07)
RASTER_WIN  = (-0.5, 1.0)       # (s) span of the raster column
RASTER_SQUARE = True            # draw the raster box square (trials are not a
                                # time axis, so a wide box just stretches it)
PANEL_W = {"mi": 5.4, "raster": 4.0, "scatter": 4.4}   # relative column widths
FIG_H   = 4.2                   # row height (inches)

# ── which neurons ────────────────────────────────────────────────────────────
PICK        = None              # [(mouse, sess, phyID), ...] to name them, or
                                # None to take N representatives
N_EXAMPLES  = 40                # how many representatives
EXAMPLE_MODE = "spread"         # "spread" = evenly over the peak-MI ranking
                                # "top"    = the highest peak MI
                                # "random" = N at random
ACROSS_MICE = True              # rotate over mice/sessions while picking, so the
                                # examples are not all from one recording
REGIONS     = ("CP", "CB")

# ── peak detection (kept identical to 04, so the marked peak is the same) ────
SMOOTH_WIN        = 3
PEAK_PROM_RATIO   = 0.5
PEAK_RUNNERUP_MAX = 0.7

# ── output ───────────────────────────────────────────────────────────────────
OUT_SUBDIR = "MI_curve_examples"
SAVE_FIGS  = True
SAVE_PDF   = False              # vector copies as well (slower, bigger)
SHOW_FIGS  = True
DPI        = 130
CLR = {"CP": np.array([106, 176, 134]) / 255, "CB": np.array([154, 123, 193]) / 255}
MCLR = {"ksg": None, "bin": "0.35", "ross": "0.65"}   # None = the region colour
plt.rcParams["pdf.fonttype"] = 42


# %% ── 1. load 04's curves and the spike/kinematics cache ───────────────────
RES = {}
for m in MI_METHODS:
    p = config.var_path(f"{MI_RUN}_{m}.pkl")
    if not p.is_file():
        raise FileNotFoundError(f"{p}\nRun 04 with RUN_NAME={MI_RUN!r} first, or "
                                f"point MI_RUN at a run that exists.")
    with open(p, "rb") as f: RES[m] = pickle.load(f)
    print(f"loaded {p.name}: mi {RES[m]['mi'].shape}, "
          f"null {'yes' if RES[m].get('null') is not None else 'no'}")

R0 = RES[MI_METHODS[0]]
UIDS    = list(R0["uids"])
REGION  = list(R0["regions"])
WC      = np.asarray(R0["win_centers"], float)
WS      = np.asarray(R0["win_starts"], float)
WW      = float(R0["WIN_WIDTH"])
COORD   = int(R0["COORD"])
FRAMES  = tuple(R0["FRAMES"])
CAT_KEEP = set(R0["CAT_KEEP"])
for m, r in RES.items():                      # the curves must describe one run
    if list(r["uids"]) != UIDS or not np.allclose(r["win_centers"], WC):
        raise ValueError(f"{m}: different neurons/windows than {MI_METHODS[0]}")
print(f"{len(UIDS)} neurons | {WC.size} windows, {WW*1000:.0f} ms wide, "
      f"{WC[0]*1000:.0f}..{WC[-1]*1000:.0f} ms | kin = endpoint coord {COORD} "
      f"(px), frames {FRAMES}, cat {sorted(CAT_KEEP)}")

if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())
BY_UID = {n["uid"]: i for i, n in enumerate(NEURONS)}
ROW    = {u: j for j, u in enumerate(UIDS)}          # uid -> row of the mi matrix
missing = [u for u in UIDS if u not in BY_UID]
if missing: print(f"!! {len(missing)} uids in the pickle are not in the cache")


# %% ── 2. helpers (peak finder copied from 04 so the marked peak matches) ───
def _smooth(arr, win=SMOOTH_WIN):
    if win <= 1: return arr.copy()
    n = len(arr); out = np.full(n, np.nan); half = win // 2
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        seg = arr[lo:hi]
        if np.any(~np.isnan(seg)): out[i] = np.nanmean(seg)
    return out


def _find_peak(arr, prom_ratio=PEAK_PROM_RATIO, runnerup_max=PEAK_RUNNERUP_MAX):
    from scipy.signal import find_peaks
    valid = ~np.isnan(arr)
    if valid.sum() < 3: return None
    a = np.where(valid, arr, -np.inf)
    idx, props = find_peaks(a, prominence=0, distance=max(1, SMOOTH_WIN))
    if len(idx) == 0: return None
    proms, heights = props["prominences"], a[idx]
    order = np.argsort(heights)[::-1]; top = order[0]
    th, tp = heights[top], proms[top]
    runner = heights[order[1]] if len(order) > 1 else 0.0
    if th <= 0 or tp / th < prom_ratio or runner > runnerup_max * th: return None
    return int(idx[top]), float(tp / th), float(th)


def trial_data(uid):
    """-> (endpoint y px, counts (n_trials, n_wins), spike lists) for kept trials.
    Rebuilt with the parameters stored IN THE PICKLE, so it is the same data the
    saved MI was computed from."""
    n   = NEURONS[BY_UID[uid]]
    kin = mi_io.endpoint_kin(n, REACH_PX, coord=COORD, frames=FRAMES)
    cat = n["labels"]["idx_reach_cat"]
    msk = np.isfinite(kin) & np.isin(cat, list(CAT_KEEP))
    C   = np.asarray(mi_rates.windowed(n["spikes"], WS, WW, "count"), float)[msk]
    spk = [np.asarray(s, float) for s, keep in zip(n["spikes"], msk) if keep]
    return kin[msk], C, spk


def _win_index(uid, mi_sm):
    """which window the scatter uses -> (index, label)"""
    if SCATTER_WIN == "peak":
        pk = _find_peak(mi_sm)
        if pk is not None: return pk[0], "peak of the smoothed curve"
        return int(np.nanargmax(mi_sm)), "max (no dominant peak)"
    lo, hi = SCATTER_WIN
    return int(np.argmin(np.abs(WC - (lo + hi) / 2))), "fixed window"


# %% ── 3. choose the neurons ─────────────────────────────────────────────────
def _peak_mi(j, m=None):
    v = RES[m or MI_METHODS[0]]["mi"][j]
    return np.nanmax(v) if np.any(np.isfinite(v)) else np.nan

if PICK:
    SEL = []
    for mouse, sess, phy in PICK:
        hits = [u for u in UIDS if u.split("|")[1] == mouse
                and u.split("|")[2] == sess and u.split("|")[3] == str(phy)]
        if not hits:
            raise KeyError(f"no neuron {mouse}|{sess}|{phy} in {MI_RUN}; "
                           f"check the phyID (uids look like {UIDS[0]!r})")
        SEL += hits
else:
    rng = np.random.default_rng(0)
    SEL = []
    for reg in REGIONS:
        rows = [j for j, r in enumerate(REGION)
                if r == reg and np.isfinite(_peak_mi(j)) and UIDS[j] in BY_UID]
        rows.sort(key=lambda j: -_peak_mi(j))
        k = min(N_EXAMPLES, len(rows))
        if EXAMPLE_MODE == "top":      take = rows[:k]
        elif EXAMPLE_MODE == "random": take = list(rng.choice(rows, k, replace=False))
        else:                          take = [rows[int(round(t))]
                                               for t in np.linspace(0, len(rows) - 1, k)]
        if ACROSS_MICE:
            # reorder so consecutive examples come from different mice/sessions
            by = {}
            for j in take: by.setdefault(UIDS[j].split("|")[1], []).append(j)
            mice, out = sorted(by), []
            while any(by.values()):
                for mm in mice:
                    if by[mm]: out.append(by[mm].pop(0))
            take = out
        SEL += [UIDS[j] for j in take]
print(f"\n{len(SEL)} neurons selected"
      + (f" (PICK)" if PICK else f" ({EXAMPLE_MODE}, {N_EXAMPLES}/region)"))


# %% ── 4. the figure ─────────────────────────────────────────────────────────
fdir = config.fig_dir(OUT_SUBDIR)

def plot_uid(uid, save=True):
    """One row of panels for one neuron. Returns the figure."""
    j   = ROW[uid]
    reg = REGION[j]
    col = CLR[reg]
    kin, C, spk = trial_data(uid)
    mi0 = RES[MI_METHODS[0]]["mi"][j]
    sm0 = _smooth(mi0)
    wi, wlab = _win_index(uid, sm0)
    cnt = C[:, wi]
    w_lo, w_hi = WS[wi], WS[wi] + WW

    _w = [PANEL_W.get(q, 5.0) for q in PANELS]
    fig, axes = plt.subplots(1, len(PANELS), figsize=(sum(_w), FIG_H),
                             gridspec_kw=dict(width_ratios=_w), squeeze=False)
    axes = axes[0]
    for ax, what in zip(axes, PANELS):
        if what == "mi":
            for m in MI_METHODS:
                c = MCLR.get(m) or col
                v = RES[m]["mi"][j]
                ax.plot(WC, v, color=c, lw=0.8, alpha=0.35)
                ax.plot(WC, _smooth(v), color=c, lw=1.8, label=m)
                nul = RES[m].get("null")
                if nul is not None and np.any(np.isfinite(nul[j])):
                    ax.plot(WC, nul[j], color=c, lw=1.0, ls="--", alpha=0.6,
                            label=f"{m} null")
            pk = _find_peak(sm0)
            if pk is not None:
                ax.axvline(WC[pk[0]], color=col, lw=1.2, ls="--", alpha=0.8,
                           label=f"peak {WC[pk[0]]*1000:+.0f} ms")
            ax.axvspan(w_lo, w_hi, color="0.85", zorder=0, lw=0)
            ax.axvline(0, color="k", lw=0.8, ls=":")
            ax.axhline(0, color="k", lw=0.4)
            ax.set(xlabel="window centre re reach (s)", ylabel="MI (bits)",
                   title=f"MI curve ({', '.join(MI_METHODS)})")
            ax.legend(fontsize=6.5)
        elif what == "raster":
            order = np.argsort(kin)
            for row, ti in enumerate(order):
                s = spk[ti]
                s = s[(s >= RASTER_WIN[0]) & (s < RASTER_WIN[1])]
                ax.plot(s, np.full(s.size, row), "|", ms=2, color=col, mew=0.5)
            ax.axvspan(w_lo, w_hi, color="0.85", zorder=0, lw=0)
            ax.axvline(0, color="k", lw=0.8, ls=":")
            ax.set(xlabel="time re reach (s)", ylabel="trial (sorted by endpoint y)",
                   xlim=RASTER_WIN, ylim=(0, len(order)),
                   title=f"raster, {len(order)} trials")
            if RASTER_SQUARE: ax.set_box_aspect(1)
        elif what == "scatter":
            # the trials themselves, nothing else: KSG bins nothing and fits
            # nothing. It counts k nearest neighbours in the joint space, so it
            # responds to any dependence between these two numbers -- including
            # one with no trend in the mean -- and a binned curve drawn on top
            # would suggest a quantity it never computes.
            ax.plot(kin, cnt, "o", ms=4, alpha=0.4, color=col, mew=0)
            _mi = RES[MI_METHODS[0]]["mi"][j][wi]
            ax.set(xlabel=f"endpoint y (px), frames {FRAMES}",
                   ylabel=f"spike count in [{w_lo*1000:.0f}, {w_hi*1000:.0f}] ms",
                   title=f"{wlab}: MI = {_mi:.3f} bits")
        ax.spines[["top", "right"]].set_visible(False)

    _pk = _find_peak(sm0)
    fig.suptitle(f"{uid}   |   {reg}   |   peak MI {np.nanmax(sm0):.3f} bits"
                 + (f" at {WC[_pk[0]]*1000:+.0f} ms" if _pk else "  (no dominant peak)")
                 + f"   |   {MI_RUN}, {WW*1000:.0f} ms windows", fontsize=10)
    plt.tight_layout()
    if save and SAVE_FIGS:
        stem = fdir / uid.replace("|", "_")
        fig.savefig(f"{stem}.png", dpi=DPI, bbox_inches="tight")
        if SAVE_PDF: fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    return fig


def show(mouse, sess, phyID):
    """Draw one neuron on demand:  show("1_CoteDor", "R2", 322)"""
    hits = [u for u in UIDS if u.split("|")[1] == mouse
            and u.split("|")[2] == sess and u.split("|")[3] == str(phyID)]
    if not hits:
        here = sorted({u.split("|")[3] for u in UIDS
                       if u.split("|")[1] == mouse and u.split("|")[2] == sess},
                      key=lambda x: int(x))
        raise KeyError(f"no neuron {mouse}|{sess}|{phyID}. phyIDs in that session: "
                       f"{here[:20]}{' ...' if len(here) > 20 else ''}")
    f = plot_uid(hits[0], save=False)
    plt.show()
    return f


for _u in SEL:
    _f = plot_uid(_u)
    if SHOW_FIGS: plt.show()
    plt.close(_f)
print(f"{len(SEL)} figures -> {fdir}")
print('single neuron on demand:  show("1_CoteDor", "R2", 322)')
