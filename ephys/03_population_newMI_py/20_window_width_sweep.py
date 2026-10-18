# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# HOW WIDE SHOULD THE SPIKE-COUNT WINDOW BE?
#
# The width has been 70 ms in 04 and 100 ms in 09/12, chosen by eye. This script
# measures it instead. Everything is held fixed (the kinematic FEATURES, the endpoint
# kinematic window, KSG, the same neurons, the same trials, the same window
# CENTRES) and only the WIDTH of the neural window changes:
#
#        excess MI( FEATURES at the endpoint ; spike count in [c - w/2, c + w/2] )   for each w
#
# FEATURES is one variable ("y", as in 04) or several, in which case the MI is the
# JOINT multivariate KSG over all of them, exactly as in 09/12. The features
# themselves are built by 09's own code (mm, Savitzky-Golay speed, water), so
# they cannot differ from the ones the real analyses use.
#
# A narrow window has little entropy to carry information (most CP trials are 0
# spikes); a wide one integrates past the timescale on which the activity
# relates to the movement and blurs timing. The excess MI therefore rises with
# the width and then flattens. THE CRITERION: the smallest width that reaches
# KNEE_FRAC of the plateau -- most of the information, at the best temporal
# resolution that still has it.
#
# Three things keep the choice honest:
#   * the criterion is the AMOUNT of information, not the lead/lag result the
#     project is about, so the width is not tuned towards the hypothesis;
#   * the knee is read on HALF of the sessions (fold A) and the other half
#     (fold B) is drawn next to it as a check that it is not noise;
#   * ONE width is to be used for both regions, so CP and CB are both shown and
#     the decision between their knees is left visible, not made silently.
# Panel E then shows what the width does to the TIMING (peak time of each region's mean curve),
# which is the thing that must not depend on it.
#
# The summary is the PEAK OF THE POPULATION-MEAN curve, not the mean of each
# neuron's own peak: a per-neuron peak is a maximum over noisy windows, and that
# maximum is larger for noisier (narrower) windows, which would bias the sweep.
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import re, time, pickle
from joblib import Parallel, delayed

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
import mi_io, config, mi_rates

# ── what is swept ────────────────────────────────────────────────────────────
WIDTHS      = (0.025, 0.050, 0.070, 0.100, 0.150, 0.200)   # (s) neural window widths
KNEE_FRAC   = 0.50     # the knee = smallest width reaching this fraction of the
                       # largest value across WIDTHS

# ── what is held fixed (as in 04) ────────────────────────────────────────────
FEATURES    = ("x", "y", "z", "speed", "dist_water")
                       # the kinematic variable(s), read in the endpoint window.
                       # any of: x, y, z, speed, accel, PC1, PC2, dist_water,
                       # speed_water. One name, e.g. ("y",), is the 1-D case of 04
                       # (and ~1.5x faster); several = joint MI as in 09/12.
                       # dist_water / speed_water need 10_behavior_water.py's pickle.
FRAMES      = (147, 167)        # endpoint window, MATLAB 1-indexed
CAT_KEEP    = {1,}   # reach categories
MIN_TRIALS  = 20
K_KSG       = 3
CENTER_STEP = 0.025    # (s) step between window CENTRES
CENTER_RANGE = (-0.5, 0.5)      # (s) re reach endpoint. The CENTRES are the same
                       # for every width (a window is [c - w/2, c + w/2]), so the
                       # widths are compared at identical times.
SUMMARY_RANGE = (-0.3, 0.3)     # (s) centres over which the secondary summary
                       # (mean excess MI) is taken, and inside which the peak
                       # of the mean curve is searched

# ── how much to run ──────────────────────────────────────────────────────────
N_TEST      = None     # neurons PER REGION, drawn at random (fixed seed) from
                       # those with enough trials, so every mouse is sampled.
                       # None = all. Cost = neurons x len(WIDTHS) x centres x
                       # (1 + N_SHUF) one-dimensional KSG calls; cell 2 prints
                       # the ETA.
N_SHUF      = 10       # trial shuffles per window for the null (excess = obs - null)
N_BOOT      = 300      # bootstrap resamples (over neurons) for the error bars
N_JOBS      = -1
SEED        = 0

CACHE_MODE  = "load"   # "load" = reuse a run whose parameters match; "rerun"
SAVE_VARS   = True
SAVE_FIGS   = True
SHOW_FIGS   = True
OUT_SUBDIR  = "window_width_sweep"
RUN_NAME    = "wsweep"

CLR = {"CP": np.array([106, 176, 134]) / 255,
       "CB": np.array([154, 123, 193]) / 255}
WIDTHS  = tuple(sorted(WIDTHS))
centers = np.arange(CENTER_RANGE[0], CENTER_RANGE[1] + CENTER_STEP / 2, CENTER_STEP)
n_c, n_w = centers.size, len(WIDTHS)
print(f"widths {[int(round(w*1000)) for w in WIDTHS]} ms  x  {n_c} centres "
      f"({centers[0]*1000:.0f}..{centers[-1]*1000:.0f} ms)")


# %% ── 1. data, estimator, neurons, folds ───────────────────────────────────
# 1. cell --> 09's cells 0-2 are exec'd (the trick 12 uses) with this script's
# FEATURES / FRAMES / CAT_KEEP, which gives the endpoint feature matrix per
# session, the trial masks and the KSG estimator from the code the real analyses
# use. Then the neurons are picked and the SESSIONS split into two folds.
_src9 = open(os.path.join(_HERE, "09_MI_multiD_kin.py"), encoding="utf-8").read()
_m3   = re.search(r"(?m)^# %% ── 3\.", _src9)
if _m3 is None: raise RuntimeError("could not find 09's cell 3 -- renumbered?")
_head = _src9[:_m3.start()]
for _nm, _val in (("FEATURES", tuple(FEATURES)), ("FRAMES", tuple(FRAMES)),
                  ("CAT_KEEP", set(CAT_KEEP)), ("N_TEST", None)):
    _head, _n = re.subn(rf"(?m)^{_nm}\s*=.*$", f"{_nm} = {_val!r}", _head, count=1)
    if not _n: raise RuntimeError(f"could not find {_nm} in 09 -- renamed?")
# 09 also assigns names this script owns; they are put back after the exec
_MINE = {k: globals()[k] for k in
         ("FEATURES", "FRAMES", "CAT_KEEP", "MIN_TRIALS", "K_KSG", "N_TEST",
          "N_JOBS", "CACHE_MODE", "SAVE_VARS", "SAVE_FIGS", "SHOW_FIGS",
          "OUT_SUBDIR", "RUN_NAME", "CLR")}
# 09's cells print a description of 09's OWN run (its windows, its output
# folder, its PCA, its neuron count). None of that describes this script, so
# it is not shown: everything printed below states what THIS run uses. If 09
# fails, its output is printed, so the error is not hidden.
import io as _io, contextlib as _ctx
_buf = _io.StringIO()
try:
    with _ctx.redirect_stdout(_buf):
        exec(compile(_head, "09_head", "exec"))
except BaseException:
    print(_buf.getvalue()); raise
globals().update(_MINE)
FEAT_TAG = "-".join(FEATURES)
N_FEAT   = len(FEATURES)
# 09 leaves, per neuron i: kin_list[i] (n_trials, n_features) and mask_list[i]
# (all features finite AND the reach in CAT_KEEP); _ksg_md is its joint KSG,
# which at one feature is numerically the estimator 04 uses.

rng = np.random.default_rng(SEED)
cand = {"CP": [], "CB": []}
KIN, MSK = {}, {}
for i, n in enumerate(NEURONS):
    if n["region"] not in cand: continue
    msk = mask_list[i]
    if msk.sum() < MIN_TRIALS: continue
    KIN[i], MSK[i] = kin_list[i][msk], msk
    cand[n["region"]].append(i)
run_idx = []
for reg in ("CP", "CB"):
    c = cand[reg]
    run_idx += sorted(c if N_TEST is None or N_TEST >= len(c)
                      else rng.choice(c, N_TEST, replace=False).tolist())
regions_run = np.array([NEURONS[i]["region"] for i in run_idx])
uids        = [NEURONS[i]["uid"] for i in run_idx]

# folds: sessions alternate A/B WITHIN each mouse, so both folds hold every mouse
_sess = sorted({(NEURONS[i]["mouse"], NEURONS[i]["sess"]) for i in run_idx})
FOLD, _cnt = {}, {}
for m, s in _sess:
    FOLD[(m, s)] = "AB"[_cnt.get(m, 0) % 2]; _cnt[m] = _cnt.get(m, 0) + 1
fold_run = np.array([FOLD[(NEURONS[i]["mouse"], NEURONS[i]["sess"])] for i in run_idx])
print(f"features: {list(FEATURES)} ({'KSG' if N_FEAT == 1 else 'joint KSG'}, "
      f"k={K_KSG}), endpoint frames {FRAMES}  |  neural variable: count  |  "
      f"reach cat {sorted(CAT_KEEP)}")
for reg in ("CP", "CB"):
    r = regions_run == reg
    print(f"{reg}: {int(r.sum())} neurons of {len(cand[reg])} usable  "
          f"(fold A {int((r & (fold_run == 'A')).sum())}, "
          f"fold B {int((r & (fold_run == 'B')).sum())})")


# %% ── 2. the sweep ─────────────────────────────────────────────────────────
# 2. cell --> for every neuron and width: the MI curve over the centres, and its
# shuffle null. The same trial permutations are used for every window and every
# width of a neuron, so the widths differ only in the window, not in the null's
# luck.
def _sweep_one(kin, spikes, widths, centers, k, n_shuf, seed):
    """-> (mi, null) each (n_widths, n_centres), plus p_empty (n_widths,) at the
    centre nearest 0."""
    rs    = np.random.default_rng(seed)
    perms = [rs.permutation(kin.shape[0]) for _ in range(n_shuf)]   # whole ROWS
    # are permuted, so the features stay together and only the pairing with
    # the spikes is broken
    mi  = np.full((len(widths), centers.size), np.nan)
    nu  = np.full_like(mi, np.nan)
    pe  = np.full(len(widths), np.nan)
    c0  = int(np.argmin(np.abs(centers)))
    for wi, w in enumerate(widths):
        C = np.asarray(mi_rates.windowed(spikes, centers - w / 2, w, "count"), float)
        pe[wi] = float(np.mean(C[:, c0] == 0))
        for ci in range(centers.size):
            c = C[:, ci]
            if np.unique(c).size < 2: continue
            mi[wi, ci] = _ksg_md(kin, c, k, seed + ci)
            nu[wi, ci] = np.mean([_ksg_md(kin[p], c, k, seed + ci) for p in perms])
    return mi, nu, pe


KEY = dict(widths=list(WIDTHS), centers=[float(centers[0]), float(centers[-1]),
                                         CENTER_STEP],
           features=list(FEATURES), FRAMES=list(FRAMES), CAT_KEEP=sorted(CAT_KEEP),
           MIN_TRIALS=MIN_TRIALS, K_KSG=K_KSG, N_SHUF=N_SHUF, N_TEST=N_TEST,
           SEED=SEED, uids=uids, v="v1")
# the file names carry the features, the categories and the neuron count, so
# sweeps with different settings are different files, never one overwriting
# the other
_STEM = (f"{RUN_NAME}_{FEAT_TAG}_cat{''.join(str(c) for c in sorted(CAT_KEEP))}"
         f"_n{'all' if N_TEST is None else N_TEST}")
_P = config.var_path(f"{_STEM}.pkl")
_c = None
if CACHE_MODE == "load" and _P.is_file():
    with open(_P, "rb") as fh: _c = pickle.load(fh)
    if _c.get("key") != KEY:
        print("saved sweep was run with other parameters -> recomputing"); _c = None
if _c is not None:
    MI, NU, PE = _c["mi"], _c["null"], _c["p_empty"]
    print(f"loaded {_P.name} -- no recompute")
else:
    _calls = len(run_idx) * n_w * n_c * (1 + N_SHUF)
    _rate = 4000 if N_FEAT == 1 else 2500        # joint KSG is slower than 1-D
    print(f"ETA: ~{_calls:,} KSG calls; at ~{_rate} calls/s ({N_FEAT}-D) that is "
          f"~{_calls/_rate/60:.0f} min")
    _t = time.time()
    res = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
        delayed(_sweep_one)(KIN[i], [s for s, m in zip(NEURONS[i]["spikes"], MSK[i]) if m],
                            WIDTHS, centers, K_KSG, N_SHUF, 100 + pi)
        for pi, i in enumerate(run_idx))
    MI = np.stack([r[0] for r in res]); NU = np.stack([r[1] for r in res])
    PE = np.stack([r[2] for r in res])
    print(f"sweep done in {(time.time()-_t)/60:.1f} min")
    if SAVE_VARS:
        with open(_P, "wb") as fh:
            pickle.dump(dict(key=KEY, mi=MI, null=NU, p_empty=PE, uids=uids,
                             regions=regions_run, fold=fold_run, widths=WIDTHS,
                             centers=centers), fh, protocol=4)
        print(f"saved -> {_P}")
EXC = MI - NU                                   # (neurons, widths, centres)


# %% ── 3. summaries and the knee ────────────────────────────────────────────
# 3. cell --> per region and fold: the population-mean curve at each width, its
# peak (primary) and its mean over SUMMARY_RANGE (secondary), with bootstrap
# error bars over neurons, and the knee.
_sr = (centers >= SUMMARY_RANGE[0]) & (centers <= SUMMARY_RANGE[1])


def _summ(E):
    """E (n, widths, centres) -> the peak of the population-mean curve (its
    height), the TIME of that peak, and the mean of the curve, all inside
    SUMMARY_RANGE, each (widths,). The peak is that of the MEAN curve -- the
    same quantity 04's mean-landscape figure marks -- not a per-neuron peak."""
    m  = np.nanmean(E, axis=0)[:, _sr]
    pk = np.nanargmax(m, axis=1)
    return m[np.arange(n_w), pk], centers[_sr][pk], np.nanmean(m, axis=1)


def _knee(v):
    v = np.asarray(v, float)
    return WIDTHS[int(np.flatnonzero(v >= KNEE_FRAC * np.nanmax(v))[0])]


S = {}
_rb = np.random.default_rng(SEED + 1)
for reg in ("CP", "CB"):
    for fold in ("A", "B", "all"):
        sel = (regions_run == reg) & ((fold_run == fold) if fold != "all" else True)
        E = EXC[sel]
        if E.shape[0] < 5: continue
        pk, tpk, mn = _summ(E)
        bs = [_summ(E[_rb.integers(0, E.shape[0], E.shape[0])]) for _ in range(N_BOOT)]
        S[(reg, fold)] = dict(
            n=E.shape[0], curve=np.nanmean(E, axis=0), peak=pk, t_peak=tpk, mean=mn,
            peak_ci=np.percentile([b[0] for b in bs], [2.5, 97.5], axis=0),
            t_ci=np.percentile([b[1] for b in bs], [2.5, 97.5], axis=0),
            mean_ci=np.percentile([b[2] for b in bs], [2.5, 97.5], axis=0),
            knee=_knee(pk), knee_boot=np.array([_knee(b[0]) for b in bs]),
            p_empty=np.nanmedian(PE[sel], axis=0))

print(f"\npeak of the population-mean excess MI (bits), by width   "
      f"[knee = smallest width >= {KNEE_FRAC:.0%} of the largest]")
print(f"{'':10s} " + " ".join(f"{int(round(w*1000)):>7d}" for w in WIDTHS) + "   knee")
for (reg, fold), d in S.items():
    kb = d["knee_boot"]
    vals, cnts = np.unique(kb, return_counts=True)
    dist = ", ".join(f"{int(round(v*1000))}:{100*c/kb.size:.0f}%"
                     for v, c in zip(vals, cnts))
    print(f"{reg} {fold:3s} {d['n']:4d}" + " "
          + " ".join(f"{v:7.4f}" for v in d["peak"])
          + f"   {int(round(d['knee']*1000))} ms   (bootstrap: {dist})")
for reg in ("CP", "CB"):
    if (reg, "all") in S and np.argmax(S[(reg, "all")]["peak"]) == n_w - 1:
        print(f"  NOTE {reg}: the largest value is at the WIDEST window tested, so "
              f"the plateau has not been reached -- the knee is a lower bound; "
              f"add wider WIDTHS.")
print(f"\nPEAK time of the population-mean curve (ms re reach endpoint), by width "
      f"[95% bootstrap CI]:")
for reg in ("CP", "CB"):
    if (reg, "all") in S:
        d = S[(reg, "all")]
        print(f"  {reg}: " + "  ".join(
            f"{int(round(w*1000))}ms {t*1000:+.0f} [{lo*1000:+.0f},{hi*1000:+.0f}]"
            for w, t, lo, hi in zip(WIDTHS, d["t_peak"], d["t_ci"][0], d["t_ci"][1])))
if all((r, "all") in S for r in ("CP", "CB")):
    print("  CB - CP: " + "  ".join(
        f"{int(round(w*1000))}ms {(b - a)*1000:+.0f}"
        for w, a, b in zip(WIDTHS, S[("CP", "all")]["t_peak"], S[("CB", "all")]["t_peak"]))
        + "   (< 0: CB earlier)")
print(f"\nmedian fraction of trials with 0 spikes in the window at the endpoint:")
for reg in ("CP", "CB"):
    if (reg, "all") in S:
        print(f"  {reg}: " + "  ".join(f"{int(round(w*1000))} ms {100*p:.0f}%"
                                       for w, p in zip(WIDTHS, S[(reg, 'all')]["p_empty"])))


# %% ── 4. the figure ────────────────────────────────────────────────────────
# A, B  the population-mean excess-MI curve at every width, CP and CB
# C     the criterion: peak of that curve vs width, fold A (solid, where the knee
#       is read) and fold B (dashed, the check), normalised to each curve's max
# D     the same in bits, all neurons, with the secondary summary
# E     what the width does to TIMING: peak time of the mean curve vs width
# F     why narrow windows fail: trials with 0 spikes vs width
_wms  = np.array(WIDTHS) * 1000
_cmap = plt.cm.viridis(np.linspace(0.05, 0.92, n_w))
fig, axes = plt.subplots(2, 3, figsize=(18, 9.6))

for ax, reg in zip(axes[0, :2], ("CP", "CB")):
    d = S.get((reg, "all"))
    if d is None: ax.set_visible(False); continue
    for wi, w in enumerate(WIDTHS):
        ax.plot(centers * 1000, d["curve"][wi], color=_cmap[wi], lw=1.8,
                label=f"{int(round(w*1000))} ms")
    ax.axvspan(SUMMARY_RANGE[0]*1000, SUMMARY_RANGE[1]*1000, color="0.94", zorder=0, lw=0)
    ax.axvline(0, color="k", lw=0.7, ls=":"); ax.axhline(0, color="k", lw=0.6)
    ax.set(xlabel="neural window CENTRE re reach endpoint (ms)",
           ylabel="population-mean excess MI (bits)",
           title=f"{'AB'[reg == 'CB']}  {reg} (n={d['n']}): mean curve at each width\n"
                 f"grey = SUMMARY_RANGE, where the peak is read")
    ax.title.set_fontsize(10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=7, title="width", title_fontsize=7)

ax = axes[0, 2]
for reg in ("CP", "CB"):
    for fold, ls, lw in (("A", "-", 2.2), ("B", "--", 1.4)):
        d = S.get((reg, fold))
        if d is None: continue
        mx = np.nanmax(d["peak"])
        ax.plot(_wms, d["peak"] / mx, color=CLR[reg], ls=ls, lw=lw, marker="o", ms=4,
                label=f"{reg} fold {fold} (n={d['n']})"
                      + (f": knee {int(round(d['knee']*1000))} ms" if fold == "A" else ""))
        if fold == "A":
            ax.fill_between(_wms, d["peak_ci"][0] / mx, d["peak_ci"][1] / mx,
                            color=CLR[reg], alpha=0.15, lw=0)
            ax.axvline(d["knee"] * 1000, color=CLR[reg], lw=1.2, ls=":")
ax.axhline(KNEE_FRAC, color="k", lw=0.8, ls=":")
ax.set(xlabel="neural window WIDTH (ms)", ylabel="peak excess MI / its maximum",
       title=f"C  THE CRITERION: smallest width reaching {KNEE_FRAC:.0%}\n"
             f"read on fold A (solid, 95% CI); fold B (dashed) is the check")
ax.title.set_fontsize(10)
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=7)

ax = axes[1, 0]
for reg in ("CP", "CB"):
    d = S.get((reg, "all"))
    if d is None: continue
    ax.plot(_wms, d["peak"], color=CLR[reg], lw=2, marker="o", ms=4,
            label=f"{reg} peak of mean curve")
    ax.fill_between(_wms, d["peak_ci"][0], d["peak_ci"][1], color=CLR[reg], alpha=0.15, lw=0)
    ax.plot(_wms, d["mean"], color=CLR[reg], lw=1.3, ls="--", marker="s", ms=3,
            label=f"{reg} mean over {SUMMARY_RANGE} s")
ax.axhline(0, color="k", lw=0.6)
ax.set(xlabel="neural window WIDTH (ms)", ylabel="excess MI (bits)",
       title="D  the same in bits, all neurons\n(dashed = the secondary summary)")
ax.title.set_fontsize(10)
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=7)

ax = axes[1, 1]
for reg in ("CP", "CB"):
    d = S.get((reg, "all"))
    if d is None: continue
    ax.plot(_wms, d["t_peak"] * 1000, color=CLR[reg], lw=2, marker="o", ms=4, label=reg)
    ax.fill_between(_wms, d["t_ci"][0] * 1000, d["t_ci"][1] * 1000,
                    color=CLR[reg], alpha=0.15, lw=0)
ax.axhline(0, color="k", lw=0.6, ls=":")
ax.set(xlabel="neural window WIDTH (ms)",
       ylabel="peak time of the mean curve (ms re reach endpoint)",
       title="E  what the width does to TIMING: time of the peak of the\n"
             "population-mean curve (95% CI). The CP/CB order should hold")
ax.title.set_fontsize(10)
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

ax = axes[1, 2]
for reg in ("CP", "CB"):
    d = S.get((reg, "all"))
    if d is None: continue
    ax.plot(_wms, 100 * d["p_empty"], color=CLR[reg], lw=2, marker="o", ms=4, label=reg)
ax.set(xlabel="neural window WIDTH (ms)", ylim=(0, 100),
       ylabel="trials with 0 spikes in the window (%, median neuron)",
       title="F  why narrow windows fail\n(window centred on the reach endpoint)")
ax.title.set_fontsize(10)
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

fig.suptitle(f"Neural window width sweep  —  KSG k={K_KSG}, endpoint "
             f"[{', '.join(FEATURES)}] (frames {FRAMES}), cat {sorted(CAT_KEEP)}, "
             f"{len(run_idx)} neurons, excess = MI - mean of {N_SHUF} shuffles",
             fontsize=12)
plt.tight_layout(rect=(0, 0, 1, 0.96))
if SAVE_FIGS:
    _f = config.fig_dir(OUT_SUBDIR) / f"{_STEM}.png"
    fig.savefig(_f, dpi=150, bbox_inches="tight"); print(f"figure -> {_f}")
if SHOW_FIGS: plt.show()
plt.close(fig)

# %%
