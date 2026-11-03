# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# The TIME x TIME MI maps of 12, split by CELL TYPE (CB only).
#
# No MI is computed here. 12's cache already holds the map of EVERY neuron, so
# this script only chooses which neurons to look at:
#   per cell type   its mean map and its individual neuron maps (own folder)
#   across types    the mean maps on one colour scale, the per-type average of
#                   f(neural) and h(lag), and dR2 (diagonal vs neural time) per type
# No CP and no CB - CP comparison: this is about the cell types of one region.
#
# Cell types come from the `cell_type` field of the .mat (config.DATA_MAT),
# matched to 12's neurons by uid. mi_io does not keep that field, so it is read
# here directly; mi_io and the parsed cache are unchanged.
# In st_eg_neurons.mat only CB has the field:
#   PkC_cs 351, NaN 285, MFB 249, PkC_ss 132, GoC 54, MLI 54   (1125 neurons)
# NaN is called "unlabelled" here.
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy.stats import kruskal, mannwhitneyu
import pickle, time, importlib, collections
import config, mi_io
import mi2d_figures
importlib.reload(mi2d_figures)          # pick up edits to the module

# ── which cell types ─────────────────────────────────────────────────────────
CELL_TYPES = "each"    # "each"  = every type separately, one after the other
                       #           (incl. "unlabelled"), each with its own outputs
                       # a list  = only these, e.g. ["PkC_ss", "PkC_cs"];
                       #           a tuple inside the list POOLS types into one
                       #           set, e.g. ["MFB", ("PkC_ss", "PkC_cs")]
                       # None    = no split: all of the region as one set
TYPE_REGION = "CB"     # the region whose neurons are split by type
MIN_TYPE_N  = 5       # a type with fewer neurons (with a map) is skipped
DIFF_PAIRS  = [("PkC_ss", "PkC_cs")]   # cell 4b: difference maps A - B, one figure
                       # per pair; pairs whose sets were not run are skipped

# ── which run of 12 to read (must match 12's settings of that run) ───────────
FEATURES_2D  = ("x", "y", "z", "speed", "dist_water")
REACH_CAT    = (1, )
TRIAL_SELECT = None    # as in 12; None / {"succ": 1} = every reach of the category
KIN_WIN, KIN_STEP, KIN_RANGE = 21, 4, (-0.34, 0.10)
NEU_WIN_WIDTH, NEU_WIN_STEP, NEU_WIN_RANGE = 0.07, 0.025, (-0.5, 0.5)
N_SHUF_2D    = 20
MIN_TRIALS   = 20
SRC_RUN_NAME = "td2D_v2"               # 12's RUN_NAME
SRC_SUBDIR   = "MI_2D_kin_x_neuron_v2" # 12's OUT_SUBDIR

# ── figures ──────────────────────────────────────────────────────────────────
N_EXAMPLES   = 20      # individual maps per type. None = EVERY neuron of the type
                       # (PkC_cs alone is ~350 PNGs); 0 = none
EXAMPLE_MODE = "spread"  # "spread" over the dR2 order | "top" | "random" (as in 12)
LAG_BIN, CV_GAP = None, None   # as in 12 (None = the neural step / the overlap)
EQUAL_AXES = False
CMAP       = "viridis"

SHOW_FIGS = True
SAVE_FIGS = True
SAVE_VARS = True
OUT_SUBDIR = "MI2D_celltypes"   # kept short: the example-map paths come close to
                              # Windows' 260-character limit
RUN_NAME   = "td2D_ct"

FPS = 400.0                                   # camera frame rate (as in 09)
N_FRAMES = 201                                # frames per tracked reach


# %% ── 1. 12's maps and grid ─────────────────────────────────────────────────
# 1. cell --> load 12's cache, check it is the run asked for, rebuild the axes.
SEL_TAG = "" if not TRIAL_SELECT else "-".join(
    f"{k}{'+'.join(str(x) for x in np.atleast_1d(v))}"
    for k, v in sorted(TRIAL_SELECT.items()))
CAT_TAG = "".join(str(c) for c in sorted(REACH_CAT))
NW_TAG  = f"nw{int(round(NEU_WIN_WIDTH * 1000))}"
SRC = config.var_path(f"{SRC_RUN_NAME}_{'-'.join(FEATURES_2D)}_{NW_TAG}_kw{KIN_WIN}"
                      f"_ks{KIN_STEP}_nall_s{N_SHUF_2D}"
                      f"{'_' + SEL_TAG if SEL_TAG else ''}_cat{CAT_TAG}_2Dcache.pkl")
if not SRC.is_file():
    raise FileNotFoundError(f"12's cache not found: {SRC.name} -- run 12 with these "
                            f"settings first (N_TEST = None)")
t0 = time.time()
with open(SRC, "rb") as f:
    _c = pickle.load(f)
KEY = _c["key"]
_want = dict(features=list(FEATURES_2D), KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP,
             KIN_RANGE=list(KIN_RANGE), N_SHUF_2D=N_SHUF_2D, WIN_WIDTH=NEU_WIN_WIDTH,
             WIN_STEP=NEU_WIN_STEP, WIN_RANGE=list(NEU_WIN_RANGE),
             CAT_KEEP=sorted(REACH_CAT), MIN_TRIALS=MIN_TRIALS, N_TEST=None)
_bad = [k for k, v in _want.items() if KEY.get(k) != v]
if (KEY.get("select") or "") != SEL_TAG: _bad.append("select")
if _bad:
    raise ValueError(f"{SRC.name} was made with other settings: "
                     + ", ".join(f"{k}={KEY.get(k)!r}" for k in _bad))

EXC2        = _c["mi"] - _c["null"]           # excess MI, as in 12
KEEP_12     = list(_c["keep"])                # neurons with a map in 12
uids        = list(_c["uids"])                # one per row of the maps
regions_run = list(_c["regions"])
run_idx     = np.arange(len(uids))            # rows ARE the run (uids[run_idx[pi]])
K_KSG, NEURAL_VAR = KEY["K_KSG"], KEY["NEURAL_VAR"]
del _c

# the axes, built as 12 builds them; checked against the cache's sizes
win_starts  = np.arange(NEU_WIN_RANGE[0], NEU_WIN_RANGE[1] - NEU_WIN_WIDTH
                        + NEU_WIN_STEP / 2, NEU_WIN_STEP)
win_centers = win_starts + NEU_WIN_WIDTH / 2
kin_starts  = np.arange(0, N_FRAMES - KIN_WIN + 1, KIN_STEP)
kin_centers = (kin_starts + KIN_WIN / 2 - config.REACH_ENDPOINT_FRAME) / FPS
kin_centers = kin_centers[(kin_centers >= KIN_RANGE[0]) & (kin_centers <= KIN_RANGE[1])]
n_kin, n_wins = len(kin_centers), len(win_centers)
assert (n_kin, n_wins) == (KEY["n_kin"], KEY["n_wins"]) == EXC2.shape[1:], \
    f"grid {n_kin}x{n_wins} does not match the cache {EXC2.shape[1:]}"

OUT_DIR_12 = f"{SRC_SUBDIR}_{NW_TAG}_cat{CAT_TAG}" + (f"_{SEL_TAG}" if SEL_TAG else "")
OUT_RUN    = f"{OUT_SUBDIR}_{TYPE_REGION}_{OUT_DIR_12.replace(SRC_SUBDIR + '_', '')}"
print(f"12's maps: {SRC.name}  ({time.time() - t0:.0f}s)")
print(f"  {len(KEEP_12)}/{len(uids)} neurons with a map, grid {n_kin} kin x {n_wins} neural")
print(f"  outputs -> {OUT_RUN}/")


# %% ── 2. cell types from the .mat ───────────────────────────────────────────
# 2. cell --> cell_type per uid, then the neuron sets.
def _cell_types(path, region):
    """-> {uid: cell_type} for every neuron of `region`, or None when its struct
    has no cell_type field. Missing / NaN / empty -> "unlabelled"."""
    import h5py
    out = {}
    with h5py.File(str(path), "r") as f:
        g = f[config.REGION_STRUCT[region]]
        if "cell_type" not in g: return None
        for i in range(g["phyID"].shape[1]):
            mouse = mi_io._str_from_ref(f, g["mouse"][0, i])
            sess  = mi_io._str_from_ref(f, g["sess"][0, i])
            phy   = int(np.asarray(f[g["phyID"][0, i]], float).ravel()[0])
            ref   = g["cell_type"][0, i]
            try:
                ct = "" if mi_io._is_matlab_empty(f[ref]) else mi_io._str_from_ref(f, ref)
            except (ValueError, TypeError):
                ct = ""                                   # numeric NaN, not a string
            out[f"{region}|{mouse}|{sess}|{phy}"] = ("unlabelled"
                                                     if ct in ("", "NaN", "nan") else ct)
    return out

t0 = time.time()
CTYPE = _cell_types(config.DATA_MAT, TYPE_REGION)
if CTYPE is None:
    raise ValueError(f"{TYPE_REGION} has no cell_type in {config.DATA_MAT.name}")
print(f"cell_type read from {config.DATA_MAT.name} ({time.time() - t0:.0f}s)")

# every neuron of the region with a map must have a type, or the uids do not match
_ty = {pi: CTYPE.get(uids[pi]) for pi in KEEP_12 if regions_run[pi] == TYPE_REGION}
_miss = [uids[pi] for pi, t in _ty.items() if t is None]
if _miss:
    raise ValueError(f"{len(_miss)} {TYPE_REGION} uids of 12 not found in the .mat, "
                     f"e.g. {_miss[:3]}")
_cnt = collections.Counter(_ty.values())

# the sets: name -> the types it pools
if CELL_TYPES is None:
    SETS = {"all": tuple(_cnt)}
elif CELL_TYPES == "each":
    SETS = {t: (t,) for t, _ in _cnt.most_common()}
else:
    SETS = {(t if isinstance(t, str) else "+".join(t)):
            ((t,) if isinstance(t, str) else tuple(t)) for t in CELL_TYPES}
_unknown = {t for ts in SETS.values() for t in ts} - set(_cnt)
if _unknown:
    raise ValueError(f"unknown cell types {sorted(_unknown)}; there are {sorted(_cnt)}")

KEEP_SET = {}
for name, ts in SETS.items():
    sel = sorted(pi for pi, t in _ty.items() if t in ts)
    if len(sel) < MIN_TYPE_N:
        print(f"  skip {name}: {len(sel)} neurons < MIN_TYPE_N={MIN_TYPE_N}")
        continue
    KEEP_SET[name] = sel

# neurons per type and mouse
_mice = sorted({uids[pi].split("|")[1] for pi in _ty})
print(f"\n{TYPE_REGION} neurons with a map in 12, per type and mouse:")
print(f"  {'':12s}{'total':>7s}" + "".join(f"{m[:8]:>9s}" for m in _mice))
for name, kp in KEEP_SET.items():
    c = collections.Counter(uids[pi].split("|")[1] for pi in kp)
    print(f"  {name:12s}{len(kp):7d}" + "".join(f"{c.get(m, 0):9d}" for m in _mice))

# colours of the C4 classifier (Beau et al.), copied from npyx/c4/plots_functions.py
# (C4_COLORS); C4 has no "unlabelled", which is grey here. A pooled set ("A+B")
# or a type not in the table falls back to tab10.
C4_COLORS = {"PkC_ss": [28, 120, 181], "PkC_cs": [0, 0, 0], "MLI": [224, 85, 159],
             "MFB": [214, 37, 41], "GrC": [143, 103, 169], "GoC": [56, 174, 62],
             "unlabelled": [150, 150, 150]}
TYPE_CLR = {n: (np.array(C4_COLORS[n]) / 255 if n in C4_COLORS else plt.cm.tab10(i % 10))
            for i, n in enumerate(KEEP_SET)}


# %% ── 3. per type: mean map and individual maps ─────────────────────────────
# 3. cell --> for each set: its mean map (mean2D.png) and N_EXAMPLES individual
# maps (individual_landscapes/), in <OUT_RUN>/<type>/. The individual maps are
# 12's figure (map, f(neural), h(lag)), ordered by dR2 as in 12.
ASPECT = "equal" if EQUAL_AXES else "auto"
extent = [win_centers[0], win_centers[-1], kin_centers[0], kin_centers[-1]]
ZERO_TXT = "reach endpoint"
_XL = f"NEURAL window centre re {ZERO_TXT} (s)"
_YL = f"KINEMATIC window centre re {ZERO_TXT} (s)"
LAB = (f"KSG k={K_KSG}, {NEURAL_VAR}, joint [{', '.join(FEATURES_2D)}], "
       f"cat {sorted(REACH_CAT)}{', ' + SEL_TAG if SEL_TAG else ''} — EXCESS")


def _draw(ax, M, title, vmin, vmax):
    im = ax.imshow(M, origin="lower", aspect=ASPECT, extent=extent, cmap=CMAP,
                   vmin=vmin, vmax=vmax)
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    ax.plot(d, d, color="w", lw=1.0, ls="--", alpha=0.7)     # simultaneity
    ax.axvline(0, color="w", lw=0.8, ls=":", alpha=0.7)
    ax.axhline(0, color="w", lw=0.8, ls=":", alpha=0.7)
    ax.set(xlabel=_XL, ylabel=_YL, title=title)
    return im


def _ctx(name, keep, fdir):
    return SimpleNamespace(
        EXC2=EXC2, keep=keep, uids=uids, run_idx=run_idx, regions_run=regions_run,
        win_centers=win_centers, kin_centers=kin_centers, n_kin=n_kin, n_wins=n_wins,
        WIN_STEP=NEU_WIN_STEP, KIN_WIN=KIN_WIN, KIN_STEP=KIN_STEP,
        CLR={TYPE_REGION: TYPE_CLR[name]}, lab=f"{LAB} — {TYPE_REGION} {name}",
        fdir=fdir, RUN_NAME=f"{RUN_NAME}_{name.replace('+', '-')}",
        SAVE_FIGS=SAVE_FIGS, SHOW_FIGS=SHOW_FIGS, extent=extent, ASPECT=ASPECT,
        CMAP=CMAP, zero_txt=ZERO_TXT, TIMES={})


RES = {}
for name, keep in KEEP_SET.items():
    t0 = time.time()
    sub  = f"{OUT_RUN}/{name.replace('+', '-')}"
    fdir = config.fig_dir(sub) if SAVE_FIGS else None
    print(f"\n{'=' * 70}\n{TYPE_REGION} {name}: {len(keep)} neurons"
          f"{'' if fdir is None else f'  -> {sub}/'}\n{'=' * 70}")
    ctx = _ctx(name, keep, fdir)
    M = np.nanmean(EXC2[keep], axis=0)
    fig, ax = plt.subplots(figsize=(7, 5))
    im = _draw(ax, M, f"{TYPE_REGION} {name} (n={len(keep)})",
               min(0.0, np.nanpercentile(M, 1)), np.nanpercentile(M, 99))
    plt.colorbar(im, ax=ax, label="mean excess MI (bits)")
    fig.suptitle(LAB, fontsize=9)
    plt.tight_layout()
    if SAVE_FIGS:
        mi2d_figures._savefig_retry(fig, fdir / f"{ctx.RUN_NAME}_mean2D.png",
                                    dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)
    # per neuron: f(neural), h(lag), dR2 (12's cell 5) -- prints its table
    dvn = mi2d_figures.diagonal_vs_neural(ctx, LAG_BIN=LAG_BIN, CV_GAP=CV_GAP)
    if N_EXAMPLES != 0 and SAVE_FIGS:
        mi2d_figures.individual_landscapes(ctx, dvn, sub, N_EXAMPLES=N_EXAMPLES,
                                           EXAMPLE_MODE=EXAMPLE_MODE)
    RES[name] = dict(types=SETS[name], keep=keep, uids=[uids[pi] for pi in keep],
                     mean=M, dvn=dvn)
    print(f"    [{name}] {time.time() - t0:.0f}s")


# %% ── 4. across types: mean maps side by side ───────────────────────────────
# 4. cell --> every type's mean map side by side, twice:
#   shared  ONE colour scale, so a type with weaker MI also looks weaker
#           (mean2D_by_type.png)
#   own     each map scaled to itself (1st-99th percentile, own colour bar), so
#           the SHAPE of a weak type is visible too (mean2D_by_type_ownscale.png)
names = list(RES)
fdir_all = config.fig_dir(OUT_RUN) if SAVE_FIGS else None

def _lims(M):
    return min(0.0, np.nanpercentile(M, 1)), np.nanpercentile(M, 99)

if names:
    _shared = _lims(np.concatenate([RES[n]["mean"].ravel() for n in names]))
    nc = min(3, len(names)); nr = int(np.ceil(len(names) / nc))
    for scale, fname, ttl in (
            ("shared", "mean2D_by_type", "one colour scale"),
            ("own", "mean2D_by_type_ownscale", "each map on its own colour scale")):
        fig, axes = plt.subplots(nr, nc, figsize=(5.6 * nc + (1.2 if scale == "own" else 0),
                                                  4.3 * nr), squeeze=False)
        for ax in axes.flat[len(names):]: ax.set_visible(False)
        for i, (ax, n) in enumerate(zip(axes.flat, names)):
            vmin, vmax = _shared if scale == "shared" else _lims(RES[n]["mean"])
            im = _draw(ax, RES[n]["mean"], f"{n} (n={len(RES[n]['keep'])})", vmin, vmax)
            ax.title.set_color(TYPE_CLR[n])
            if i % nc: ax.set_ylabel("")                       # label the outer
            if i + nc < len(names): ax.set_xlabel("")          # edges only
            if scale == "own":
                plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03,
                             label="bits" if i % nc == nc - 1 or i == len(names) - 1 else None)
        fig.subplots_adjust(hspace=0.3, wspace=0.25 if scale == "own" else 0.15)
        if scale == "shared":
            fig.colorbar(im, ax=axes, label="mean excess MI (bits)", shrink=0.8)
        fig.suptitle(f"{TYPE_REGION} mean 2D landscape by cell type — {ttl}\n{LAB}",
                     fontsize=10)
        if SAVE_FIGS:
            mi2d_figures._savefig_retry(fig, fdir_all / f"{RUN_NAME}_{fname}.png",
                                        dpi=150, bbox_inches="tight")
        if SHOW_FIGS: plt.show()
        plt.close(fig)


# %% ── 4b. difference maps between types ─────────────────────────────────────
# 4b. cell --> for each pair in DIFF_PAIRS: map A, map B (one colour scale) and
# z(A) - z(B). For the difference only, each mean map is z-scored over all its
# cells, so a type with much less MI overall (PkC_cs) does not turn the whole
# difference one colour: red = relatively more MI in A there, blue = in B.
from matplotlib.colors import TwoSlopeNorm
for a, b in DIFF_PAIRS:
    if a not in RES or b not in RES:
        print(f"skip {a} - {b}: not both in this run ({list(RES)})"); continue
    Ma, Mb = RES[a]["mean"], RES[b]["mean"]
    vmin, vmax = _lims(np.r_[Ma.ravel(), Mb.ravel()])
    fig, axes = plt.subplots(1, 3, figsize=(18, 4.8))
    for ax, n, M in ((axes[0], a, Ma), (axes[1], b, Mb)):
        im = _draw(ax, M, f"{n} (n={len(RES[n]['keep'])})", vmin, vmax)
        ax.title.set_color(TYPE_CLR[n])
        plt.colorbar(im, ax=ax, label="mean excess MI (bits)")
    _zm = lambda M: (M - np.nanmean(M)) / np.nanstd(M)   # over all cells of the map
    D = _zm(Ma) - _zm(Mb)
    lim = np.nanpercentile(np.abs(D), 99)
    im = axes[2].imshow(D, origin="lower", aspect=ASPECT, extent=extent, cmap="RdBu_r",
                        norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim))
    d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
    axes[2].plot(d, d, color="k", lw=1.0, ls="--", alpha=0.6)
    axes[2].axvline(0, color="k", lw=0.8, ls=":"); axes[2].axhline(0, color="k", lw=0.8, ls=":")
    axes[2].set(xlabel=_XL, ylabel=_YL, title=f"z({a}) - z({b})   (shape only)")
    plt.colorbar(im, ax=axes[2], label="Δ z")
    fig.suptitle(f"{TYPE_REGION} {a} vs {b}\n{LAB}", fontsize=10)
    plt.tight_layout()
    if SAVE_FIGS:
        mi2d_figures._savefig_retry(fig, fdir_all / f"{RUN_NAME}_diff2D_{a}-minus-{b}.png",
                                    dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 5. across types: f(neural), h(lag) and dR2 ────────────────────────────
# 5. cell --> per type, the mean +- SEM over its neurons of
#   f(neural)  each neuron's map averaged down each column  (WHEN the spikes are)
#   h(lag)     ... averaged along each diagonal              (the LAG to the movement)
# raw (bits: also shows how much MI) and z-scored per type (shape only, so a
# weak type's timing can be compared with a strong one's). Right: dR2 per neuron
# (> 0 diagonal-like, < 0 neural-time-like), Kruskal-Wallis across types and
# each type against the rest.
def _mean_sem(P):
    P = np.asarray(P, float)
    n = np.isfinite(P).sum(axis=0)
    return np.nanmean(P, axis=0), np.nanstd(P, axis=0, ddof=1) / np.sqrt(np.maximum(n, 1))

def _z(m):
    sd = np.nanstd(m)
    return (m - np.nanmean(m)) / sd if sd > 0 else m * np.nan

if names:
    lag_axis = RES[names[0]]["dvn"]["lag_axis"]
    lag_full = RES[names[0]]["dvn"]["lag_full"]
    fig = plt.figure(figsize=(17, 8))
    gs  = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.8], wspace=0.3, hspace=0.35)
    ax_n, ax_nz = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    ax_l, ax_lz = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[:, 2])
    for n in names:
        dv, kp, c = RES[n]["dvn"], RES[n]["keep"], TYPE_CLR[n]
        for (ax, axz, x, prof) in ((ax_n, ax_nz, win_centers, "PROF_NEU"),
                                   (ax_l, ax_lz, lag_axis, "PROF_LAG")):
            m, s = _mean_sem([dv[prof][pi] for pi in kp])
            ax.fill_between(x, m - s, m + s, color=c, alpha=0.2, lw=0)
            ax.plot(x, m, color=c, lw=1.8, label=f"{n} (n={len(kp)})")
            axz.plot(x, _z(m), color=c, lw=1.8)
    for ax in (ax_l, ax_lz):
        mi2d_figures.shade_partial(ax, lag_axis, lag_full, scale=1.0)
    for ax, t in ((ax_n, "f(neural), mean ± SEM"), (ax_nz, "f(neural), z-scored per type"),
                  (ax_l, "h(lag), mean ± SEM"), (ax_lz, "h(lag), z-scored per type")):
        ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6, ls=":")
        ax.set_title(t, fontsize=10); ax.spines[["top", "right"]].set_visible(False)
    ax_n.set_ylabel("excess MI (bits)"); ax_nz.set_ylabel("z")
    ax_nz.set_xlabel(_XL)
    ax_lz.set_xlabel("lag = neural - kinematic (s)   < 0: spikes BEFORE the movement")
    ax_n.legend(fontsize=8, frameon=False)

    D = {n: np.array([RES[n]["dvn"]["DR2"][pi] for pi in RES[n]["keep"]]) for n in names}
    D = {n: d[np.isfinite(d)] for n, d in D.items()}
    rs = np.random.default_rng(0)
    for i, n in enumerate(names):
        ax_d.boxplot(D[n], positions=[i], widths=0.5, showfliers=False,
                     medianprops=dict(color="k"))
        ax_d.scatter(i + rs.uniform(-0.18, 0.18, D[n].size), D[n], s=5,
                     color=TYPE_CLR[n], alpha=0.5, lw=0)
    ax_d.axhline(0, color="k", lw=0.6)
    ax_d.set_xticks(range(len(names)), names, rotation=30, ha="right")
    ax_d.set(ylabel="dR2 = R2_lag - R2_neural   (> 0 diagonal-like)")
    ax_d.spines[["top", "right"]].set_visible(False)
    _ok = [D[n] for n in names if D[n].size >= 5]
    if len(_ok) >= 2:
        H, p_kw = kruskal(*_ok)
        ax_d.set_title(f"dR2 per neuron\nKruskal-Wallis p = {p_kw:.1e}", fontsize=10)
        print(f"\ndR2 across types: Kruskal-Wallis H = {H:.1f}, p = {p_kw:.1e}")
    print(f"{'type':12s} {'n':>5s} {'dR2 median [IQR]':>24s} {'% diagonal':>11s} "
          f"{'vs rest (MW p)':>15s}")
    for n in names:
        rest = np.concatenate([D[o] for o in names if o != n]) if len(names) > 1 else []
        p = (mannwhitneyu(D[n], rest).pvalue if len(rest) and D[n].size else np.nan)
        q = np.percentile(D[n], [25, 75]) if D[n].size else (np.nan, np.nan)
        print(f"{n:12s} {D[n].size:5d} {np.median(D[n]):+9.3f} [{q[0]:+.2f},{q[1]:+.2f}]"
              f"{100 * np.mean(D[n] > 0):10.0f}% {p:15.1e}")
    fig.suptitle(f"{TYPE_REGION} cell types compared — {LAB}", fontsize=10)
    if SAVE_FIGS:
        mi2d_figures._savefig_retry(fig, fdir_all / f"{RUN_NAME}_profiles_by_type.png",
                                    dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 6. save ───────────────────────────────────────────────────────────────
# 6. cell --> per set: its neurons, mean map and per-neuron f(neural), h(lag),
# dR2. The maps themselves stay in 12's cache (src), not copied here.
if SAVE_VARS and RES:
    p = config.var_path(f"{RUN_NAME}_{TYPE_REGION}_{'-'.join(FEATURES_2D)}_{NW_TAG}"
                        f"_ks{KIN_STEP}_cat{CAT_TAG}{'_' + SEL_TAG if SEL_TAG else ''}"
                        f"_celltypes.pkl")
    with open(p, "wb") as f:
        pickle.dump(dict(src=SRC.name, src_key=KEY, type_region=TYPE_REGION,
                         cell_type={uids[pi]: t for pi, t in _ty.items()},
                         sets={n: dict(types=r["types"], uids=r["uids"], mean=r["mean"],
                                       R2_neural=[r["dvn"]["R2N"][pi] for pi in r["keep"]],
                                       R2_lag=[r["dvn"]["R2L"][pi] for pi in r["keep"]],
                                       dR2=[r["dvn"]["DR2"][pi] for pi in r["keep"]],
                                       prof_neural=[r["dvn"]["PROF_NEU"][pi] for pi in r["keep"]],
                                       prof_lag=[r["dvn"]["PROF_LAG"][pi] for pi in r["keep"]])
                               for n, r in RES.items()},
                         lag_axis=lag_axis, kin_centers=kin_centers,
                         win_centers=win_centers), f, protocol=4)
    print(f"saved -> {p}")

# %%
