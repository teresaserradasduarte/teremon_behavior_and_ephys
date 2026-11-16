# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# CHOOSE THE MAPS TO LABEL BY HAND (for 27_label_gui.py, then 28's classifier)
#
# 23 (hot cluster) and 25 (templates) mislabel too many maps, so the shapes are
# labelled by hand: column, diagonal, tiling, broad, nothing, mixed/unsure.
# This script picks WHICH maps, and stores three views of each for the GUI:
#     excess MI   MI - null mean (bits)
#     z           (MI - null mean) / null SD
#     smoothed z  z, NaN-aware Gaussian (23's SMOOTH_SIG)
#
# First pass, 200 maps:
#   180 real, half CP half CB; in each region
#       RANDOM_FRAC  picked at random from all the region's maps (what a typical
#                    map looks like; an unbiased check of the class mix)
#       N_DISAGREE   where 23 and 25 give different classes
#       the rest     spread evenly over 25's classes
#   20 shuffle maps (23's), mixed in and not marked in the GUI: how often noise
#      is called a shape
# Repeat pass, 50 more: maps of the first pass shown again in a new order, to
# measure the agreement of the labeller with herself.
# 23's and 25's classes only SPREAD the picks; the GUI never shows them.
import os, sys, pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

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
                raise ModuleNotFoundError("mi_io.py not found from " + os.getcwd())
            _p = os.path.dirname(_p)
    _HERE = _found
if _HERE not in sys.path: sys.path.insert(0, _HERE)
os.chdir(_HERE)
import config

# ── which maps: 12's run, 23's shuffle maps and classes, 25's classes ───────
REACH_CAT     = (1, )
FEATURES_2D   = ("x", "y", "z", "speed", "dist_water")
KIN_WIN, KIN_STEP, KIN_RANGE = 21, 4, (-0.34, 0.10)
NEU_WIN_WIDTH, NEU_WIN_STEP, NEU_WIN_RANGE = 0.07, 0.025, (-0.5, 0.5)
SRC_12     = "td2D_v2_x-y-z-speed-dist_water_nw70_kw21_ks4_nall_s20_cat1_2Dcache.pkl"
SRC_23     = "shapes2D_x-y-z-speed-dist_water_nw70_ks4_cat1_nall_shufmap_cache.pkl"
SRC_23_RES = "shapes2D_x-y-z-speed-dist_water_nw70_ks4_cat1_nall_results.pkl"
SRC_25_RES = "tmpl2D_x-y-z-speed-dist_water_nw70_ks4_cat1_nall_results.pkl"

# ── the selection ────────────────────────────────────────────────────────────
N_PER_REGION  = 90     # real maps per region (CP, CB)
RANDOM_FRAC   = 0.10   # share of them picked at random from the whole region
N_DISAGREE    = 27     # maps where 23 and 25 disagree, per region
N_SHUF_PER_REGION = 10 # shuffle maps per region
N_REPEAT      = 50     # repeat pass, in addition to the first pass
N_REPEAT_SHUF = 5      # of which shuffle maps
MUST_INCLUDE  = ("Ferrero|R4|278", "Lindt|R6|996")   # "mouse|session|unit": the mouse
                       # part only has to be IN the mouse name (5_FerreroRocher)
MIN_FINITE    = 0.6    # maps defined in fewer cells are left out
SMOOTH_SIG    = (1.0, 0.5)   # (kin, neural) cells, as 23
SEED          = 7

SAVE_VARS = True
SAVE_FIGS = True
SHOW_FIGS = True
OUT_SUBDIR = "MI_2D_labels"
LABELS = ("column", "diagonal", "tiling", "2 tiles", "broad", "nothing", "mixed/unsure")
VIEWS  = ("excess MI", "z", "smoothed z")

# 23's and 25's classes on one coarse scale, only to find where they DISAGREE
COARSE = {"column": "column", "diagonal": "diagonal", "tiling": "tiling", "broad": "broad",
          "not informative": "nothing", "noise-sized": "nothing", "no cluster": "nothing",
          "scattered": "nothing", "ambiguous": "unsure", "too sparse": "too sparse"}


# %% ── 1. load the maps and the two classifications ───────────────────────────
# 1. cell --> 12's maps, 23's shuffle maps, 23's and 25's classes, by uid.
FPS = 400.0
win_starts  = np.arange(NEU_WIN_RANGE[0], NEU_WIN_RANGE[1] - NEU_WIN_WIDTH + NEU_WIN_STEP / 2,
                        NEU_WIN_STEP)
win_centers = win_starts + NEU_WIN_WIDTH / 2
kin_starts  = np.arange(0, config.REACH_N_FRAMES - KIN_WIN + 1, KIN_STEP)
kin_centers = (kin_starts + KIN_WIN / 2 - config.REACH_ENDPOINT_FRAME) / FPS
kin_centers = kin_centers[(kin_centers >= KIN_RANGE[0]) & (kin_centers <= KIN_RANGE[1])]
TAG = f"{'-'.join(FEATURES_2D)}_nw{int(round(NEU_WIN_WIDTH * 1000))}_cat{''.join(map(str, REACH_CAT))}"

with open(config.var_path(SRC_12), "rb") as fh: C12 = pickle.load(fh)
if C12["mi"].shape[1:] != (len(kin_centers), len(win_centers)):
    raise ValueError(f"12's maps {C12['mi'].shape[1:]} do not match the grid "
                     f"{len(kin_centers)} x {len(win_centers)}")
with open(config.var_path(SRC_23), "rb") as fh: C23 = pickle.load(fh)
if C23["key"].get("src") != SRC_12:
    raise ValueError(f"23's shuffle maps were not made from {SRC_12}")
with open(config.var_path(SRC_23_RES), "rb") as fh: R23 = pickle.load(fh)
with open(config.var_path(SRC_25_RES), "rb") as fh: R25 = pickle.load(fh)

uids = list(C12["uids"]); regions = np.asarray(C12["regions"])
_row23 = {u: i for i, u in enumerate(C23["key"]["uids"])}
_c23 = dict(zip(R23["uids"], R23["cls"])); _c25 = dict(zip(R25["uids"], R25["cls"]))
with np.errstate(all="ignore"):
    _fin = np.isfinite(C12["mi"] - C12["null"]).mean(axis=(1, 2))
POOL = np.array([u in _row23 and u in _c23 and u in _c25 and _fin[i] >= MIN_FINITE
                 and _c23[u] != "too sparse" and _c25[u] != "too sparse"
                 for i, u in enumerate(uids)])
C23_ = np.array([str(_c23.get(u, "")) for u in uids])
C25_ = np.array([str(_c25.get(u, "")) for u in uids])
DISAG = np.array([COARSE.get(a, a) != COARSE.get(b, b) for a, b in zip(C23_, C25_)])
for reg in ("CP", "CB"):
    m = POOL & (regions == reg)
    print(f"{reg}: {int(m.sum())} maps to choose from; 23 and 25 disagree on "
          f"{int((m & DISAG).sum())} ({100 * (m & DISAG).sum() / m.sum():.0f}%)")


# %% ── 2. the selection ──────────────────────────────────────────────────────
# 2. cell --> per region: MUST_INCLUDE, random, disagreements, spread over 25's
# classes; plus shuffle maps; then the order and the repeat pass.
rs = np.random.default_rng(SEED)


def _must(reg):
    out = []
    for e in MUST_INCLUDE:
        mo, se, un = e.split("|")
        hit = [i for i, u in enumerate(uids) if regions[i] == reg and POOL[i]
               and mo in u.split("|")[1] and u.split("|")[2] == se and u.split("|")[3] == un]
        out += hit
    return out


ITEMS = []                                    # dicts: idx (row of 12), kind, stratum
for reg in ("CP", "CB"):
    pool = set(np.flatnonzero(POOL & (regions == reg)).tolist())
    chosen = []
    def _take(cands, n, stratum):
        cands = [i for i in cands if i in pool and i not in chosen]
        pick = list(rs.choice(cands, min(n, len(cands)), replace=False)) if cands and n > 0 else []
        chosen.extend(pick)
        ITEMS.extend(dict(idx=int(i), kind="real", stratum=stratum) for i in pick)
        return len(pick)
    for i in _must(reg):
        print(f"  MUST_INCLUDE {uids[i]}  (23: {C23_[i]}, 25: {C25_[i]})")
        chosen.append(i); ITEMS.append(dict(idx=int(i), kind="real", stratum="must"))
    left = N_PER_REGION - len(chosen)
    _take(sorted(pool), int(round(RANDOM_FRAC * N_PER_REGION)), "random")
    _take([i for i in pool if DISAG[i]], N_DISAGREE, "disagree")
    cls25 = sorted({C25_[i] for i in pool})
    quota = (N_PER_REGION - len(chosen)) // len(cls25)
    for c in cls25:
        _take([i for i in pool if C25_[i] == c], quota, f"25:{c}")
    _take(sorted(pool), N_PER_REGION - len(chosen), "fill")
    # shuffle maps of this region (any neuron of the pool)
    for i in rs.choice(sorted(pool), N_SHUF_PER_REGION, replace=False):
        ITEMS.append(dict(idx=int(i), kind="shuffle", stratum="shuffle"))

for it in ITEMS:
    i = it["idx"]
    it.update(uid=uids[i], region=str(regions[i]), c23=C23_[i], c25=C25_[i],
              key=f"{it['kind']}:{uids[i]}")
assert len({it["key"] for it in ITEMS}) == len(ITEMS), "duplicate map in the selection"

ORDER = list(rs.permutation(len(ITEMS)))                       # first pass
_shuf = [j for j, it in enumerate(ITEMS) if it["kind"] == "shuffle"]
_real = [j for j, it in enumerate(ITEMS) if it["kind"] == "real"]
REPEAT = (list(rs.choice(_shuf, N_REPEAT_SHUF, replace=False))
          + list(rs.choice(_real, N_REPEAT - N_REPEAT_SHUF, replace=False)))
REPEAT = [int(j) for j in rs.permutation(REPEAT)]

print(f"\n{len(ITEMS)} maps in the first pass, {len(REPEAT)} in the repeat pass")
import collections
for reg in ("CP", "CB"):
    c = collections.Counter(it["stratum"] for it in ITEMS if it["region"] == reg)
    print(f"  {reg}: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
    c = collections.Counter(COARSE.get(it["c25"], it["c25"]) for it in ITEMS
                            if it["region"] == reg and it["kind"] == "real")
    print(f"      25's classes among the real maps: " + ", ".join(f"{k} {v}" for k, v in c.most_common()))


# %% ── 3. the three views of each map ────────────────────────────────────────
# 3. cell --> excess MI, z, smoothed z (from the shuffle map for shuffle items).
def _smooth(Zm):
    ok = np.isfinite(Zm)
    if not ok.any(): return Zm.copy()
    num = gaussian_filter(np.where(ok, Zm, 0.0), SMOOTH_SIG)
    den = gaussian_filter(ok.astype(float), SMOOTH_SIG)
    with np.errstate(all="ignore"):
        return np.where(ok, num / den, np.nan)


V = np.full((len(ITEMS), len(VIEWS), len(kin_centers), len(win_centers)), np.nan, np.float32)
for j, it in enumerate(ITEMS):
    i = it["idx"]
    M = C12["mi"][i] if it["kind"] == "real" else C23["sh"][_row23[it["uid"]]]
    with np.errstate(all="ignore"):
        exc = M - C12["null"][i]
        z = exc / np.where(C12["null_std"][i] > 0, C12["null_std"][i], np.nan)
    V[j] = (exc, z, _smooth(z))

# a few maps in each view, to check
fig, axes = plt.subplots(len(VIEWS), 6, figsize=(20, 8), squeeze=False)
ext = [win_centers[0], win_centers[-1], kin_centers[0], kin_centers[-1]]
for c, j in enumerate(ORDER[:6]):
    for r, vn in enumerate(VIEWS):
        M = V[j, r]; ax = axes[r, c]
        ax.imshow(M, origin="lower", aspect="auto", extent=ext, cmap="viridis",
                  vmin=np.nanpercentile(M, 1), vmax=np.nanpercentile(M, 99))
        ax.set_title(f"{vn}  |  {ITEMS[j]['kind']} {ITEMS[j]['uid']}", fontsize=7)
        ax.tick_params(labelsize=6)
fig.suptitle("the first 6 maps of the first pass, in the three views (the GUI hides the names)",
             fontsize=10)
plt.tight_layout()
if SAVE_FIGS: fig.savefig(config.fig_dir(OUT_SUBDIR) / f"label_set_{TAG}_check.png", dpi=130,
                          bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 4. save the label set ─────────────────────────────────────────────────
# 4. cell --> one file for the GUI. Overwriting it after labelling has begun
# would re-draw the maps, so it refuses if the file exists (delete it yourself
# to start again).
if SAVE_VARS:
    p = config.var_path(f"label_set_{TAG}.pkl")
    if p.is_file():
        print(f"{p.name} exists -- NOT overwritten (labels may refer to it). Delete it to redo.")
    else:
        with open(p, "wb") as fh:
            pickle.dump(dict(items=[{k: v for k, v in it.items() if k != "idx"} for it in ITEMS],
                             views=V, view_names=VIEWS, labels=LABELS, order=[int(j) for j in ORDER],
                             repeat=REPEAT, win_centers=win_centers, kin_centers=kin_centers,
                             tag=TAG,
                             params=dict(N_PER_REGION=N_PER_REGION, RANDOM_FRAC=RANDOM_FRAC,
                                         N_DISAGREE=N_DISAGREE, N_SHUF_PER_REGION=N_SHUF_PER_REGION,
                                         N_REPEAT=N_REPEAT, MUST_INCLUDE=MUST_INCLUDE,
                                         SMOOTH_SIG=SMOOTH_SIG, SEED=SEED, SRC_12=SRC_12,
                                         SRC_23=SRC_23, SRC_23_RES=SRC_23_RES,
                                         SRC_25_RES=SRC_25_RES)), fh, protocol=4)
        print(f"label set -> {p}  ({p.stat().st_size / 1e6:.1f} MB)")
        print(f"label it with:  python 27_label_gui.py")

# %%
