# %% ── 0. imports & parameters ────────────────────────────────────────────────
# WHY THIS SCRIPT EXISTS: it draws the OLD method, the one used before the KSG
# landscape, so the thesis can show what it did and where it breaks.
#
# The old recipe, per neuron:
#   1 take a kinematic scalar per trial   (mean ML position in the endpoint window)
#   2 take a neural scalar per trial      (mean instantaneous 1/ISI rate in a window)
#   3 BIN both into histograms, bin the JOINT, and read MI off the plug-in formula
#        I = SUM p(x,y) log2[ p(x,y) / (p(x) p(y)) ]
#
# Step 3 is the problem, and the figure is built to make it visible:
#   - with N_BINS_KIN x N_BINS_NEU cells and only ~200-300 trials, most joint bins
#     are EMPTY, and the few that are filled hold 1-2 trials
#   - every empty bin biases the plug-in estimate UPWARDS, because missing mass
#     looks like structure
#   - the bias grows with the number of bins and with how SPARSE the neuron is,
#     so a low-firing CP cell and a high-firing CB cell are not on a fair scale --
#     which is exactly the comparison the project is about
#   - the figure prints the plug-in MI on SHUFFLED trials, where the true MI is
#     0 by construction. Whatever it reports there is pure bias.
# The current method (KSG) is printed next to it on the same trials for contrast.
#
# NOTE on t=0: the kinematic axis is drawn relative to the START of the endpoint
# window (frame 147), the same convention as 09/12. The true maximum extension is
# nearer frame 151, so "0" here is the window start, not peak reach.
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D            # noqa: F401  (registers "3d")
from sklearn.feature_selection import mutual_info_regression
import mi_io, mi_rates, config

# ---- what to show -----------------------------------------------------------
REGIONS    = ("CP", "CB")
NEU_MODE   = "isi"          # "isi" = instantaneous 1/ISI rate (the OLD variable)
                            # "count" = spike count (what the current method uses)
KIN_AXIS   = 1              # 0=x, 1=y (medio-lateral, the old endpoint scalar), 2=z
KIN_UNITS  = "mm"           # "mm" (project standard) or "px" (what the old
                            # MATLAB figures showed). Only rescales the axis.
FRAMES     = (147, 167)     # endpoint window, MATLAB 1-indexed, as in 09/12
FPS        = 400.0
CAT_KEEP   = {1, 2}         # reach categories to include
MIN_TRIALS = 20

TRACE_RANGE = (-0.5, 1.0)   # (s) re reach, the span drawn in the heatmaps
NEU_WIN     = (0.0, 0.5)    # (s) re reach, the window the neural scalar comes from

# ---- the binning, i.e. the thing being criticised ---------------------------
N_BINS_KIN = 10             # bins for p(kin)
N_BINS_NEU = 10             # bins for p(neural); the joint gets KIN x NEU cells
N_SHUF     = 50             # shuffles for the "plug-in MI when the truth is 0"

# ---- which neurons ----------------------------------------------------------
# Picked per region by their mean firing in NEU_WIN: a SPARSE one and a DENSE one.
SPARSE_Q = 10               # percentile of mean rate -> the sparse example
DENSE_Q  = 90               # percentile of mean rate -> the dense example
PICK_UIDS = None            # or e.g. {"CP sparse": "CP|3_Toblerone|R3|322", ...}
                            # to pin exact neurons; overrides the percentiles

# ---- output -----------------------------------------------------------------
OUT_SUBDIR = "method_binning_explained"
RUN_NAME   = "binning"
SAVE_FIGS  = True
SAVE_PDF   = True           # vector copy for the thesis
SHOW_FIGS  = True
DPI        = 200
CMAP_KIN   = "bone"         # matches the old MATLAB figures
CMAP_NEU   = "bone"
CMAP_JOINT = "viridis"
CLR = {"CP": np.array([106, 176, 134]) / 255,
       "CB": np.array([154, 123, 193]) / 255}
plt.rcParams["pdf.fonttype"] = 42     # keep PDF text editable in Illustrator
plt.rcParams["ps.fonttype"]  = 42

NAT2BIT = 1.0 / np.log(2)
EP      = slice(FRAMES[0] - 1, FRAMES[1])
AXNAME  = {0: "x (AP)", 1: "y (ML)", 2: "z (DV)"}[KIN_AXIS]


# %% ── 1. data ───────────────────────────────────────────────────────────────
if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

PX2MM    = np.asarray(config.PX2MM, float)
REACH_MM = {k: v * PX2MM[None, :, None] for k, v in REACH_PX.items()}
REACH    = REACH_MM if KIN_UNITS == "mm" else REACH_PX
print(f"{len(NEURONS)} neurons | kinematics in {KIN_UNITS} | neural variable "
      f"{NEU_MODE!r} in window {NEU_WIN} s")

# Per-neuron trial vectors. The kinematic scalar is the mean over the endpoint
# window, the neural scalar the mean over NEU_WIN -- one number each per trial,
# which is all the old method ever used.
TRIALS = {}
for i, n in enumerate(NEURONS):
    key = (n["mouse"], n["sess"])
    rp  = REACH[key]                                   # (trials, xyz, frames)
    cat = n["labels"]["idx_reach_cat"]
    kin = rp[:, KIN_AXIS, EP].mean(axis=1)             # endpoint scalar
    msk = np.isfinite(kin) & np.isin(cat, list(CAT_KEEP))
    if msk.sum() < MIN_TRIALS: continue
    neu = np.asarray(mi_rates.windowed(n["spikes"], np.array([NEU_WIN[0]]),
                                       NEU_WIN[1] - NEU_WIN[0], NEU_MODE),
                     float).ravel()[msk]
    if np.unique(neu).size < 2: continue
    TRIALS[i] = dict(kin=kin[msk], neu=neu, mask=msk, key=key,
                     region=n["region"], uid=n["uid"], rate=float(np.mean(neu)))
print(f"{len(TRIALS)} neurons usable ({MIN_TRIALS}+ trials, non-constant neural "
      f"variable)")

# Pick the sparse and dense example in each region.
PICKS = {}
if PICK_UIDS:
    _by_uid = {v["uid"]: k for k, v in TRIALS.items()}
    for lab, u in PICK_UIDS.items():
        if u not in _by_uid: raise KeyError(f"{u} not among the usable neurons")
        PICKS[lab] = _by_uid[u]
else:
    for reg in REGIONS:
        idx  = [i for i, t in TRIALS.items() if t["region"] == reg]
        rate = np.array([TRIALS[i]["rate"] for i in idx])
        for lab, q in (("sparse", SPARSE_Q), ("dense", DENSE_Q)):
            PICKS[f"{reg} {lab}"] = idx[int(np.argmin(np.abs(
                rate - np.percentile(rate, q))))]
for lab, i in PICKS.items():
    t = TRIALS[i]
    print(f"  {lab:10s} {t['uid']:28s} n={t['kin'].size:4d} trials, "
          f"mean {'rate' if NEU_MODE=='isi' else 'count'} {t['rate']:.2f}")


# %% ── 2. the binned (plug-in) estimator, and what it does to noise ──────────
def mi_plugin(x, y, nbx=N_BINS_KIN, nby=N_BINS_NEU, edges=None):
    """The OLD estimator: histogram both variables, histogram the joint, and
    evaluate  I = SUM p log2[p / (px py)]  over the filled cells.
    -> (bits, joint counts, x edges, y edges). Empty cells contribute nothing to
    the sum, which is precisely why they bias it: the mass that should have been
    there is simply missing, and missing mass reads as structure."""
    if edges is None:
        H, xe, ye = np.histogram2d(x, y, bins=[nbx, nby])
    else:
        H, xe, ye = np.histogram2d(x, y, bins=edges)
    P  = H / H.sum()
    px = P.sum(axis=1, keepdims=True)
    py = P.sum(axis=0, keepdims=True)
    nz = P > 0
    return float(np.sum(P[nz] * np.log2(P[nz] / (px @ py)[nz]))), H, xe, ye


def mi_ksg(x, y, seed=0):
    """The CURRENT estimator at d=1, for contrast. No bins at all."""
    r = np.random.default_rng(seed)
    yj = y + 1e-10 * max(1.0, np.abs(y).mean()) * r.standard_normal(y.size)
    return float(mutual_info_regression(x.reshape(-1, 1), yj, n_neighbors=3,
                                        random_state=0)[0] * NAT2BIT)


def mi_plugin_null(x, y, n_shuf=N_SHUF, seed=0):
    """Plug-in MI after breaking the pairing. The true value is 0, so whatever
    comes back is the estimator's own bias at this bin count and trial count.
    Edges are held fixed at the real ones, so only the pairing changes."""
    _, _, xe, ye = mi_plugin(x, y)
    r = np.random.default_rng(seed)
    v = [mi_plugin(x, y[r.permutation(y.size)], edges=[xe, ye])[0]
         for _ in range(n_shuf)]
    return float(np.mean(v)), float(np.std(v, ddof=1))


# %% ── 3. one figure per example neuron ──────────────────────────────────────
fdir = config.fig_dir(OUT_SUBDIR)
NEU_LAB = ("instantaneous rate (sp/s)" if NEU_MODE == "isi"
           else "spike count")

SUMMARY = {}
for lab, i in PICKS.items():
    t    = TRIALS[i]
    reg  = t["region"]
    col  = CLR[reg]
    kin, neu = t["kin"], t["neu"]
    n_tr = kin.size

    # the numbers the figure is really about
    mi_b, H, xe, ye = mi_plugin(kin, neu)
    mi_b0, mi_b0_sd = mi_plugin_null(kin, neu)
    mi_k  = mi_ksg(kin, neu)
    empty = 100 * np.mean(H == 0)
    SUMMARY[lab] = dict(uid=t["uid"], n=n_tr, rate=t["rate"], mi_binned=mi_b,
                        mi_binned_shuffled=mi_b0, mi_ksg=mi_k, pct_empty=empty)

    # traces, with the trials ORDERED by the kinematic scalar so the heatmaps
    # show the gradient the MI is supposed to pick up
    order = np.argsort(kin)
    rp    = REACH[t["key"]][t["mask"]][:, KIN_AXIS, :][order]     # (trials, frames)
    t_kin = (np.arange(rp.shape[1]) - EP.start) / FPS
    t_neu, R = mi_rates.rate_traces(NEURONS[i]["spikes"], NEU_MODE,
                                    TRACE_RANGE[0], TRACE_RANGE[1])
    R = R[t["mask"]][order]

    fig = plt.figure(figsize=(19, 6.4))
    # margins set here rather than with tight_layout: a 3-D axes is not
    # compatible with it and matplotlib warns and lays the figure out wrongly
    gs  = fig.add_gridspec(2, 4, width_ratios=[1.25, 1.25, 1.0, 1.5],
                           height_ratios=[3, 1.25], wspace=0.33, hspace=0.38,
                           left=0.045, right=0.955, top=0.84, bottom=0.095)

    # --- kinematic heatmap + the per-trial scalar underneath ----------------
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(rp, aspect="auto", origin="lower", cmap=CMAP_KIN,
                   extent=[t_kin[0], t_kin[-1], 0, n_tr])
    ax.axvspan(0, (EP.stop - EP.start) / FPS, color="w", alpha=0.18, lw=0)
    ax.axvline(0, color="w", ls="--", lw=1.2, alpha=0.8)
    ax.set(xlabel="time re endpoint window (s)", ylabel="trial (sorted by kin)",
           title=f"paw {AXNAME} trajectory")
    plt.colorbar(im, ax=ax, label=f"{AXNAME} ({KIN_UNITS})", fraction=0.046, pad=0.02)

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(np.arange(n_tr), kin[order], ".", ms=3, color=col)
    ax.set(xlabel="trial (sorted)", ylabel=f"endpoint\n{AXNAME} ({KIN_UNITS})",
           xlim=(0, n_tr))
    ax.spines[["top", "right"]].set_visible(False)

    # --- neural heatmap + the per-trial scalar underneath -------------------
    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(R, aspect="auto", origin="lower", cmap=CMAP_NEU,
                   extent=[t_neu[0], t_neu[-1], 0, n_tr],
                   vmax=np.nanpercentile(R, 99) if np.isfinite(R).any() else None)
    ax.axvspan(NEU_WIN[0], NEU_WIN[1], color="w", alpha=0.18, lw=0)
    ax.axvline(0, color="w", ls="--", lw=1.2, alpha=0.8)
    ax.set(xlabel="time re reach endpoint (s)", ylabel="trial (sorted by kin)",
           title=f"neural: {NEU_LAB}")
    plt.colorbar(im, ax=ax, label=NEU_LAB, fraction=0.046, pad=0.02)

    ax = fig.add_subplot(gs[1, 1])
    ax.plot(np.arange(n_tr), neu[order], ".", ms=3, color=col)
    ax.set(xlabel="trial (sorted by kin)", ylabel=f"window\n{NEU_LAB}",
           xlim=(0, n_tr))
    ax.spines[["top", "right"]].set_visible(False)

    # --- the two marginal histograms ---------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    ax.hist(kin, bins=xe, color=col, alpha=0.75, edgecolor="white")
    ax.set(xlabel=f"endpoint {AXNAME} ({KIN_UNITS})", ylabel="# trials",
           title=f"p(kin), {N_BINS_KIN} bins")
    ax.spines[["top", "right"]].set_visible(False)

    ax = fig.add_subplot(gs[1, 2])
    ax.hist(neu, bins=ye, color=col, alpha=0.75, edgecolor="white")
    ax.set(xlabel=NEU_LAB, ylabel="# trials", title=f"p(neural), {N_BINS_NEU} bins")
    ax.spines[["top", "right"]].set_visible(False)

    # --- the joint, in 3D, which is where the trials run out ---------------
    ax = fig.add_subplot(gs[:, 3], projection="3d")
    xc = (xe[:-1] + xe[1:]) / 2
    yc = (ye[:-1] + ye[1:]) / 2
    XX, YY = np.meshgrid(xc, yc, indexing="ij")
    dx = (xe[1] - xe[0]) * 0.85
    dy = (ye[1] - ye[0]) * 0.85
    filled = H > 0
    if filled.any():
        cols = plt.get_cmap(CMAP_JOINT)(H[filled] / H.max())
        ax.bar3d(XX[filled] - dx / 2, YY[filled] - dy / 2, np.zeros(filled.sum()),
                 dx, dy, H[filled], color=cols, shade=True, zsort="average")
    # the empty cells are the whole point: mark them on the floor
    if (~filled).any():
        ax.scatter(XX[~filled], YY[~filled], np.zeros((~filled).sum()),
                   marker="x", s=14, color="0.55", depthshade=False)
    ax.set(xlabel=f"endpoint {AXNAME} ({KIN_UNITS})", ylabel=NEU_LAB,
           zlabel="# trials",
           title=f"joint p(kin, neural): {N_BINS_KIN}x{N_BINS_NEU} cells\n"
                 f"{empty:.0f}% EMPTY (grey x), {n_tr} trials for "
                 f"{N_BINS_KIN*N_BINS_NEU} cells")
    ax.view_init(elev=32, azim=-57)

    fig.suptitle(
        f"{lab}  |  {t['uid']}  |  n={n_tr} trials, cat {sorted(CAT_KEEP)}, "
        f"mean {'rate' if NEU_MODE=='isi' else 'count'} {t['rate']:.2f}"
        f"{' sp/s' if NEU_MODE=='isi' else ''}\n"
        f"BINNED (old) MI = {mi_b:.3f} bits   |   same estimator on SHUFFLED "
        f"trials = {mi_b0:.3f} ± {mi_b0_sd:.3f} bits, where the truth is 0   "
        f"-> bias {100*mi_b0/max(mi_b, 1e-9):.0f}% of the reported value"
        f"   |   KSG (current) = {mi_k:.3f} bits",
        fontsize=10)

    if SAVE_FIGS:
        stem = f"{RUN_NAME}_{lab.replace(' ', '_')}_{t['uid'].replace('|', '_')}"
        fig.savefig(fdir / f"{stem}.png", dpi=DPI, bbox_inches="tight")
        if SAVE_PDF:
            fig.savefig(fdir / f"{stem}.pdf", bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

print(f"\n{'example':12s} {'n':>5s} {'mean':>7s} {'binned':>8s} {'shuffled':>9s} "
      f"{'KSG':>7s} {'% empty':>8s}")
print("-" * 62)
for lab, s in SUMMARY.items():
    print(f"{lab:12s} {s['n']:5d} {s['rate']:7.2f} {s['mi_binned']:8.3f} "
          f"{s['mi_binned_shuffled']:9.3f} {s['mi_ksg']:7.3f} {s['pct_empty']:7.0f}%")
print("\nthe 'shuffled' column is the bias: the true MI there is 0 by construction.")
if SAVE_FIGS: print(f"figures -> {fdir}")
