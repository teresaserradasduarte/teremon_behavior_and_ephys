# %% ── 0. parameters ─────────────────────────────────────────────────────────
# STEP 03 -- complex-spike-triggered PAUSE: an identity test that does not
# depend on C4.
#
# A Purkinje cell's simple spikes stop for ~10-20 ms after each of its own
# complex spikes. So if a unit pauses right after the complex spikes (CS) of a
# nearby unit, that unit is the simple-spike (SS) train of the same Purkinje
# cell -- whatever C4 said about it -- and the CS unit is confirmed as a CS.
#
# For every candidate CS unit (C4's best guess is PkC_cs, reported or not, or
# s04 called it PkC_cs) and every other CB unit of the same session whose peak
# channel is within MAX_CH_DIST channels:
#   CS-triggered histogram of the other unit (1 ms bins, +-50 ms)
#   baseline  = its rate 50-10 ms BEFORE the CS
#   ratio     = rate 3-8 ms after the CS / baseline            (0 = full pause)
#   asym      = (rate 8-3 ms BEFORE the CS / baseline) - ratio (one-sided pause ~ 1)
#   refr      = rate 3-8 ms after the CS / the unit's own ACG at 3-8 ms
#   share     = fraction of CS spikes with a spike of the unit within +-0.5 ms,
#               minus what its firing rate gives by chance
#   z         = Poisson z-score of the drop
# A pair PASSES when all the criteria below hold. Each one closes a known trap
# (windows and traps follow Simon's CS curation, ephys/simon/ComplexSpike_Analysis,
# step0 v7 -- re-implemented here; his code and outputs are not read or changed):
#   window 3-8 ms   0-2 ms is left out: a CS briefly suppresses the simple spikes
#                   of NEIGHBOURING Purkinje cells through its field (Han et al.
#                   2020) and the sorter misses spikes overlapping the CS; after
#                   8 ms a rebound starts.
#   SS_FR_MIN       slow units (often MFB) pass a ratio test by chance.
#   MIN_ASYM        a sorter collision dip is symmetric around the CS.
#   MAX_REFR        if pieces of the CS are sorted INTO the unit, the unit's own
#                   refractory period looks like a pause after the "CS".
#   MAX_SHARE       the same trap seen directly: the two units fire together.
# One CS per SS and one SS per CS: the strongest pause wins; the other passing
# pairs of a unit stay in the table (used = False) but are not used.
# Only extra-good units are in eg_neurons_CB.mat, so a CS whose SS partner is not
# extra-good stays unpaired here (Simon's tier 2 covers those).
# This step only measures and shows; step 04 uses the passing pairs.
import collections, time
import numpy as np
import matplotlib.pyplot as plt
import ct_config as C, ct_io as io, ct_features as F

MAX_CH_DIST = 25       # channels between the two peak channels (~250 um on NP1)
BIN, WIN    = 0.001, 0.05
BASE_WIN, PAUSE_WIN = (-0.05, -0.01), (0.003, 0.008)
# a pair passes when:
MIN_CS      = 100      # at least this many complex spikes
SS_FR_MIN   = 20.0     # Hz, the paused unit's rate before the CS (Simon: 20; 10 is his open alternative)
MAX_RATIO   = 0.30     # rate after / before   (Simon's "strong pause")
MIN_ASYM    = 0.20     # one-sided
MAX_REFR    = 0.70     # pause beyond the unit's own refractoriness
MAX_SHARE   = 0.20     # excess fraction of CS spikes shared with the unit
MAX_Z       = -5.0     # drop is this many Poisson SDs

N_GALLERY = 40         # pairs per gallery figure
SHOW_FIGS = True
SAVE_FIGS = True
SAVE_VARS = True       # writes cb_pairs.pkl


# %% ── 1. test every CS candidate against its neighbours ─────────────────────
FEAT = io.load_features()
C4 = {r["uid"]: r for r in io.read_csv(C.CB_C4_CSV)}
if not C4: raise FileNotFoundError(f"{C.CB_C4_CSV} not found -- run 02_CB_c4_review.py")
PK = {r["uid"]: r["peak_ch"] for r in FEAT["rows"]}
_i = {r["uid"]: i for i, r in enumerate(FEAT["rows"])}
_al = (np.arange(FEAT["arrays"]["acg"].shape[1]) + 0.5) * FEAT["meta"]["acg_bin"]
_aw = (_al >= PAUSE_WIN[0]) & (_al < PAUSE_WIN[1])
ACG_WIN_HZ = {u: float(FEAT["arrays"]["acg"][i][_aw].mean()) for u, i in _i.items()}   # ACG at 3-8 ms
is_cs = lambda u: C4[u]["c4_pred"] == "PkC_cs" or C4[u]["old_cell_type"] == "PkC_cs"

pairs, lags = [], None
t00 = time.time()
for mouse, sess, kind, path in io.sessions(kinds=("CB",)):
    neus = io.load_eg_neurons(path, regions=("CB",), with_wf=False)
    U = {C.uid("CB", mouse, sess, d["phyID"]): d["st"] for d in neus}
    U = {u: st for u, st in U.items() if u in C4}
    cs = [u for u in U if is_cs(u) and U[u].size >= MIN_CS]
    n0 = len(pairs)
    for c in cs:
        for u, st in U.items():
            if u == c or abs(PK[u] - PK[c]) > MAX_CH_DIST or st.size < 100: continue
            cnt, lags = F.ccg(U[c], st, BIN, WIN)
            s = F.pause_stats(cnt, lags, U[c].size, BIN, BASE_WIN, PAUSE_WIN)
            a = ACG_WIN_HZ.get(u, np.nan)
            near = (np.searchsorted(st, U[c] + 0.0005) - np.searchsorted(st, U[c] - 0.0005)) > 0
            s.update(cs=c, tgt=u, n_cs=int(U[c].size), ch_dist=int(abs(PK[u] - PK[c])),
                     refr=float(s["win_hz"] / a) if a > 0 else np.nan,
                     share=float(near.mean() - min(1.0, 0.001 * s["base_hz"])),
                     rate=(cnt / U[c].size / BIN).astype(np.float32))
            s["passed"] = bool(s["n_cs"] >= MIN_CS and s["base_hz"] >= SS_FR_MIN
                               and s["ratio"] <= MAX_RATIO and s["asym"] >= MIN_ASYM
                               and not s["refr"] > MAX_REFR and s["share"] <= MAX_SHARE
                               and s["z"] <= MAX_Z)
            pairs.append(s)
    print(f"{mouse} {sess}: {len(cs):3d} CS candidates, {len(pairs) - n0:4d} pairs tested, "
          f"{sum(p['passed'] for p in pairs[n0:]):3d} pass")
print(f"{len(pairs)} pairs, {sum(p['passed'] for p in pairs)} pass  ({time.time() - t00:.0f}s)")

PASS = [p for p in pairs if p["passed"]]
# one-to-one: go through the passing pairs from the strongest pause down and keep
# a pair only if neither of its units is already taken
best_cs_of, best_ss_of = {}, {}
for p in sorted(PASS, key=lambda p: p["ratio"]):
    p["used"] = p["tgt"] not in best_cs_of and p["cs"] not in best_ss_of
    if p["used"]: best_cs_of[p["tgt"]] = p["cs"]; best_ss_of[p["cs"]] = p["tgt"]
print(f"{len(best_cs_of)} units pause after a CS; {len(best_ss_of)} CS units have a pausing partner")
print("what the pausing units were called:")
for k, lab in (("old_cell_type", "s04"), ("c4_pred", "C4 best guess")):
    print(f"  {lab:14s}", dict(collections.Counter(C4[u][k] or "not_run" for u in best_cs_of)))
if SAVE_VARS:
    io.save_pkl(C.CB_PAIRS_PKL, dict(
        pairs=pairs, lags=lags, best_cs_of=best_cs_of, best_ss_of=best_ss_of,
        params=dict(MAX_CH_DIST=MAX_CH_DIST, BIN=BIN, WIN=WIN, BASE_WIN=BASE_WIN,
                    PAUSE_WIN=PAUSE_WIN, MIN_CS=MIN_CS, SS_FR_MIN=SS_FR_MIN,
                    MAX_RATIO=MAX_RATIO, MIN_ASYM=MIN_ASYM, MAX_REFR=MAX_REFR,
                    MAX_SHARE=MAX_SHARE, MAX_Z=MAX_Z)))
    print(f"saved -> {C.CB_PAIRS_PKL}")


# %% ── 2. figure: all tested pairs in the criteria plane ────────────────────
# Each dot = one (CS, neighbour) pair. A real CS-SS pair sits bottom right: low
# ratio, long pause. The cluster should be clearly apart from the cloud of
# unrelated pairs (ratio ~ 1, pause ~ 0); if it is not, move the thresholds.
fdir = C.fig_dir("03_CB_cs_ss_pairs") if SAVE_FIGS else None
ok = [p for p in pairs if p["base_hz"] >= SS_FR_MIN and np.isfinite(p["ratio"])]
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
ax = axes[0]
for passed, kw in ((False, dict(color="0.6", s=6, alpha=0.5)), (True, dict(color="#0b2a6f", s=14))):
    q = [p for p in ok if p["passed"] == passed]
    ax.scatter([p["asym"] for p in q], [p["ratio"] for p in q], lw=0,
               label=f"{'pass' if passed else 'fail'} ({len(q)})", **kw)
ax.axhline(MAX_RATIO, color="k", ls="--", lw=0.8); ax.axvline(MIN_ASYM, color="k", ls="--", lw=0.8)
ax.set(xlabel="asymmetry (before - after)", ylabel="rate after CS / before", ylim=(0, 2),
       xlim=(-1, 1.5), title=f"all pairs whose target fires >= {SS_FR_MIN:g} Hz")
ax.legend(frameon=False, fontsize=8)
ax = axes[1]
ax.hist([p["ratio"] for p in ok], np.linspace(0, 2, 81), color="0.5")
ax.axvline(MAX_RATIO, color="k", ls="--", lw=0.8); ax.set_yscale("log")
ax.set(xlabel="rate after CS / before", ylabel="pairs", title="ratio (log count)")
ax = axes[2]                                    # what the pausing units were called
labs = ["PkC_ss", "PkC_cs", "MLI", "GoC", "MFB", C.UNLAB]
co = collections.Counter(C4[u]["old_cell_type"] for u in best_cs_of)
cn = collections.Counter(C4[u]["c4_pred"] or C.UNLAB for u in best_cs_of)
x = np.arange(len(labs))
ax.bar(x - 0.2, [co.get(l, 0) for l in labs], 0.4, color=[C.COLORS[l] for l in labs], edgecolor="k", label="s04 label")
ax.bar(x + 0.2, [cn.get(l, 0) for l in labs], 0.4, color=[C.COLORS[l] for l in labs], alpha=0.45,
       edgecolor="k", label="C4 best guess (pale)")
ax.set_xticks(x, labs); ax.legend(frameon=False, fontsize=8)
ax.set(ylabel="units", title=f"the {len(best_cs_of)} units that pause after a CS")
for a in axes: a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "pairs_criteria.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 3. figures: galleries of CS-triggered histograms ──────────────────────
# (a) the passing pairs, strongest pause first; (b) the best FAILING pairs, i.e.
# the ones closest to the thresholds -- to see what the criteria leave out.
def _gallery(sel, title, fname):
    for page in range(0, len(sel), N_GALLERY):
        sub = sel[page:page + N_GALLERY]
        nc = 8; nr = int(np.ceil(len(sub) / nc))
        fig, axes = plt.subplots(nr, nc, figsize=(2.6 * nc, 2.1 * nr), squeeze=False)
        for ax in axes.ravel()[len(sub):]: ax.set_visible(False)
        for ax, p in zip(axes.ravel(), sub):
            t = C4[p["tgt"]]
            ax.bar(lags * 1e3, p["rate"], width=1.0, color=C.COLORS.get(t["old_cell_type"], "0.5"))
            ax.axvline(0, color="k", lw=0.6); ax.axhline(p["base_hz"], color="k", lw=0.5, ls=":")
            m, s, ph = p["tgt"].split("|")[1:]
            ax.set_title(f"{m[2:]} {s} cs{p['cs'].split('|')[-1]}>{ph}\n"
                         f"r={p['ratio']:.2f} as={p['asym']:.2f} rf={p['refr']:.2f}  s04:{t['old_cell_type']}",
                         fontsize=6.5)
            ax.tick_params(labelsize=6); ax.spines[["top", "right"]].set_visible(False)
        fig.suptitle(f"{title} — bar colour = s04 label of the paused unit; x: ms from CS, y: sp/s",
                     fontsize=10)
        fig.tight_layout()
        if SAVE_FIGS: fig.savefig(fdir / f"{fname}_{page // N_GALLERY + 1:02d}.png", dpi=130)
        if SHOW_FIGS: plt.show()
        plt.close(fig)

if lags is not None:
    _gallery(sorted(PASS, key=lambda p: p["ratio"]), "03 — PASSING CS→unit pairs", "gallery_pass")
    near = sorted([p for p in ok if not p["passed"]], key=lambda p: p["ratio"])[:N_GALLERY]
    _gallery(near, "03 — best FAILING pairs (lowest ratio)", "gallery_nearmiss")
# %%
