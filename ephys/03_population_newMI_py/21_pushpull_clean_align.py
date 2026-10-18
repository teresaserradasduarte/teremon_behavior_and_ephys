# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# CLEAN THE INVALID PUSH/PULLS, AND FIND AN ALIGNMENT VALID AND INVALID SHARE.
#
# Reads 16's pushpull_kin.pkl (both paws + load cell, +-3 s around EVERY event,
# valid and invalid) and writes pushpull_clean.pkl: per event, whether it is a
# real push/pull, and the time of its FORCE ONSET. Nothing is re-cut here.
#
# WHAT THE TIMESTAMPS ARE. Each event time is the moment the rig DETECTED it,
# which is also when the cue fires (valid: 50 ms sound + 200 ms bar vibration;
# invalid: blue light for 500 ms). Detection is not the same moment relative to
# the force for the two kinds:
#   valid    |force| must stay above its threshold for 50 ms -> detected ~70 ms
#            after the threshold is crossed
#   invalid  no hold -> detected at the crossing (~10 ms after it)
# and the thresholds themselves can differ per session (cell 2 reads them off
# the data). So "t = 0" means a different point of the force rise for valid and
# invalid events, and aligning both on it mixes two different phases.
#
# TWO ALIGNMENTS are therefore offered, and saved:
#   "cue"    the detection timestamp as recorded (what 16/18 use now). The
#            right zero for anything driven by the CUE.
#   "onset"  the FORCE ONSET: the moment the force starts to rise towards the
#            threshold. Defined the same way for every event, so valid and
#            invalid line up on the same phase of the movement.
#
# WHAT IS CLEANED. An invalid push/pull is sometimes detected while the mouse is
# REACHING: the reaching (dom) paw is up, and only the non-dominant paw is on
# the resting bar, pressing it. A real push/pull has both paws down on the bar.
# Cell 3 measures, AT THE MOMENT OF DETECTION, how HIGH each paw is above the
# height it has during the VALID push/pulls of the same session, and classifies:
#   "bar"         dom paw down on the bar                      -> real, kept
#   "reach-like"  dom paw up, non-dom paw down                 -> removed
#   "both off"    both paws up                                 -> removed
# Height only, and only at detection, on purpose. A first version used the
# largest x-z distance anywhere in [-200, +50] ms and flagged 30% of VALID pulls:
# a push/pull moves the paws FORWARD/BACK (x), and before many valid events the
# dom paw is still coming DOWN from a previous reach and lands on the bar just in
# time. Neither is a reach at the moment the rig detects the event.
#   "untracked"   too little paw tracking to judge            -> removed unless
#                                                               KEEP_UNTRACKED
# VALID events are classified the same way and FLAGGED, never removed: whatever
# the paws were doing, a valid event started a trial. The flag is saved
# (`reach_like`) so it can be used or ignored downstream, and its rate is also
# the most that the rule can be wrong by on real push/pulls.
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import pickle, time
import scipy.io as sio

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

# ── the force onset ──────────────────────────────────────────────────────────
REF_FORCE   = 1000.0   # (mV/V, signed: + for push, - for pull) the level whose
                       # crossing anchors the search. The same for EVERY event,
                       # whatever threshold the rig used in that session, so the
                       # onset is defined identically for valid and invalid. It
                       # must sit BELOW every threshold in use: cell 2 reads them
                       # off the data (2000-3000 for most sessions, down to ~1000
                       # in Toblerone); at 2000, ~6% of valid events never reach
                       # it and get no onset.
ANCHOR_WIN  = (-0.60, 0.0)      # (s re the timestamp) where that crossing is
                       # looked for; the LAST upward crossing in it is the anchor.
                       # It ends AT the timestamp: the rig can only detect an
                       # event after its force has crossed REF_FORCE (which is
                       # below every threshold), so a later crossing belongs to
                       # a second rise, not to this event.
BASE_WIN    = (-0.50, -0.30)    # (s re the anchor) the force before the rise
ONSET_FRAC  = 0.10     # onset = the force first exceeds base + this fraction of
                       # (REF_FORCE - base) on the way up to the anchor
MAX_RISE    = 0.60     # (s) a rise slower than this has no defined onset (NaN)

# ── the cleaning ─────────────────────────────────────────────────────────────
POSTURE_WIN = (-0.05, 0.05)     # (s re the timestamp) the moment of detection,
                       # where the paws are judged (median height over it)
H_UP        = 3.0      # (mm) a paw is "up" when it is this much HIGHER (z) than
                       # its height during the session's valid events. Measured
                       # at detection: valid events are below 1 mm in 90% of
                       # cases, while the invalid ones split off at 7-10 mm --
                       # the distribution is bimodal and 3 mm sits in the gap.
MIN_TRACKED = 0.5      # fraction of POSTURE_WIN frames a paw must be tracked in
KEEP_UNTRACKED = False # keep invalid events whose posture cannot be judged
REACH_TOL   = 0.50     # (s) an independent CHECK, not used to classify: is a
                       # reach detected by the behaviour pipeline within this
                       # of the event? Should be common for "reach-like" and rare
                       # for "bar".

# ── what is shown ────────────────────────────────────────────────────────────
PLOT_WIN    = (-1.0, 0.6)
EX_EVENTS   = 5        # single events per panel column, from different mice
SESSIONS    = None     # None = all, or [(mouse, sess), ...]
EXCLUDE_TOBLERONE = False   # Toblerone cannot push; set True to leave it out here
                       # too. Off by default: its events are worth seeing.

PP_PKL      = "pushpull_kin.pkl"    # from 16
OUT_PKL     = "pushpull_clean.pkl"
SAVE_VARS   = True
SAVE_FIGS   = True
SAVE_TABLES = True     # the cleaning tables as CSV, next to the figures
SHOW_FIGS   = True
OUT_SUBDIR  = "pushpull_clean"

CODE   = {0: "valid push", 1: "valid pull", 2: "invalid push", 3: "invalid pull"}
SIGN   = {0: +1, 1: -1, 2: +1, 3: -1}        # push = positive force, pull = negative
_push  = np.array([132, 160, 124]) / 256; _pull = np.array([26, 94, 99]) / 256
CLR    = {0: _push, 1: _pull, 2: _push * 0.55 + 0.45, 3: _pull * 0.55 + 0.45}
CLS_CLR = {"bar": "0.25", "reach-like": np.array([200, 60, 60]) / 255,
           "both off": np.array([230, 150, 30]) / 255, "untracked": "0.7"}


def _save(fig, name):
    if SAVE_FIGS:
        p = config.fig_dir(OUT_SUBDIR) / name
        fig.savefig(p, dpi=150, bbox_inches="tight"); print(f"figure -> {p}")
    if SHOW_FIGS: plt.show()
    plt.close(fig)


# %% ── 1. load ──────────────────────────────────────────────────────────────
# 1. cell --> 16's windows, and the reach times (only for the independent check)
with open(config.var_path(PP_PKL), "rb") as fh: _PK = pickle.load(fh)
PP, T = _PK["sessions"], np.asarray(_PK["t"], float)
FPS   = float(_PK["fps"])
if "dom" not in _PK["paws"] or "nondom" not in _PK["paws"]:
    raise KeyError("16 must have extracted BOTH paws (PAWS = ('dom', 'nondom'))")
if not _PK.get("loadcell", False):
    raise KeyError("16 must have extracted the load cell (LOADCELL = True)")
KEYS = sorted(PP) if SESSIONS is None else [tuple(k) for k in SESSIONS]
if EXCLUDE_TOBLERONE:
    KEYS = [k for k in KEYS if "toblerone" not in k[0].lower()]
REACH_T = {}
for k in KEYS:
    b = sio.loadmat(config.bhv_fundamentals_mat(*k), squeeze_me=True,
                    struct_as_record=False)["bhv"]
    REACH_T[k] = np.sort(np.asarray(b.reach_times_inVec, float).ravel())
print(f"{len(KEYS)} sessions, window {T[0]:+.1f}..{T[-1]:+.1f} s at {FPS:.0f} fps; "
      "events: " + ", ".join(f"{CODE[c]} {sum(int((np.asarray(PP[k]['code']) == c).sum()) for k in KEYS)}"
                             for c in CODE))


# %% ── 2. per event: threshold crossing, force onset, detection lag ─────────
# 2. cell --> everything in seconds RE THE TIMESTAMP (the cue), per event.
#   t_cross  last upward crossing of REF_FORCE before the timestamp (the anchor)
#   t_onset  when the force started rising towards it
#   f_det    the signed force AT the timestamp -- for invalid events, which are
#            detected at their crossing, this reads off the threshold the rig used
def _onset(fe):
    """fe: signed force on T -> (t_cross, t_onset), NaN where undefined."""
    w = np.flatnonzero((T >= ANCHOR_WIN[0]) & (T <= ANCHOR_WIN[1]))
    w = w[w > 0]
    up = w[(fe[w] >= REF_FORCE) & (fe[w - 1] < REF_FORCE)]
    if not up.size: return np.nan, np.nan
    ia = int(up[-1])
    bm = (T >= T[ia] + BASE_WIN[0]) & (T <= T[ia] + BASE_WIN[1])
    if not np.isfinite(fe[bm]).any(): return T[ia], np.nan
    base  = float(np.nanmedian(fe[bm]))
    level = base + ONSET_FRAC * (REF_FORCE - base)
    j = ia
    while j > 0 and T[j] > T[ia] - MAX_RISE and np.isfinite(fe[j - 1]) and fe[j - 1] > level:
        j -= 1
    if T[j] <= T[ia] - MAX_RISE or not np.isfinite(fe[j - 1]): return T[ia], np.nan
    return T[ia], T[j]


_z = int(np.argmin(np.abs(T)))
EV = {}
for k in KEYS:
    s, cd = PP[k], np.asarray(PP[k]["code"]).astype(int)
    F = np.asarray(s["force"], float)
    tc, to, fd = (np.full(cd.size, np.nan) for _ in range(3))
    for e in range(cd.size):
        if cd[e] not in SIGN or not np.isfinite(F[e]).any(): continue
        fe = SIGN[cd[e]] * F[e]
        tc[e], to[e] = _onset(fe)
        fd[e] = fe[_z]
    EV[k] = dict(code=cd, t_cross=tc, t_onset=to, f_det=fd)

print(f"\n{'':14s} {'n':>5s} {'onset found':>12s} {'crossing re cue (ms)':>22s} "
      f"{'onset re cue (ms)':>20s} {'onset re crossing (ms)':>24s}")
for c in CODE:
    cr = np.concatenate([EV[k]["t_cross"][EV[k]["code"] == c] for k in KEYS])
    on = np.concatenate([EV[k]["t_onset"][EV[k]["code"] == c] for k in KEYS])
    q = lambda a: (f"{np.nanmedian(a)*1000:+5.0f} [{np.nanpercentile(a,25)*1000:+.0f},"
                   f"{np.nanpercentile(a,75)*1000:+.0f}]") if np.isfinite(a).any() else "-"
    print(f"{CODE[c]:14s} {cr.size:5d} {100*np.mean(np.isfinite(on)):11.0f}% "
          f"{q(cr):>22s} {q(on):>20s} {q(on - cr):>24s}")
print(f"  (median [IQR]. 'crossing re cue' = how long before detection the force "
      f"passed {REF_FORCE:.0f}; valid events are detected later in their rise "
      f"than invalid ones, because valid detection waits for the 50 ms hold.)")
print("why an onset is missing:")
for c in CODE:
    nof = sum(int((~np.isfinite(np.asarray(PP[k]["force"], float)[EV[k]["code"] == c][:, _z])).sum())
              for k in KEYS)
    nox = sum(int((np.isfinite(np.asarray(PP[k]["force"], float)[EV[k]["code"] == c][:, _z])
                   & ~np.isfinite(EV[k]["t_cross"][EV[k]["code"] == c])).sum()) for k in KEYS)
    noo = sum(int((np.isfinite(EV[k]["t_cross"]) & ~np.isfinite(EV[k]["t_onset"])
                   & (EV[k]["code"] == c)).sum()) for k in KEYS)
    print(f"  {CODE[c]:13s} no force (16 could not place the event): {nof:4d} | "
          f"force never crosses {REF_FORCE:.0f}: {nox:4d} | rise too slow / no baseline: {noo:4d}")

# the threshold each session's rig used for INVALID events: invalid events are
# detected at their crossing, so the force at the timestamp sits on it
print(f"\nforce at the timestamp of INVALID events, per session (median, mV/V) -- "
      f"the invalid threshold in use:")
_thr = {}
for k in KEYS:
    for c in (2, 3):
        v = EV[k]["f_det"][EV[k]["code"] == c]
        if np.isfinite(v).sum() >= 5: _thr[(k, c)] = float(np.nanmedian(v))
_line = {}
for (k, c), v in _thr.items(): _line.setdefault(k, {})[c] = v
for k in KEYS:
    d = _line.get(k, {})
    print(f"  {k[0]+'|'+k[1]:22s} " + "  ".join(f"{CODE[c].split()[1]} {d[c]:6.0f}" if c in d
                                                else f"{CODE[c].split()[1]}      -" for c in (2, 3)))


# %% ── 3. CLEANING: where are the paws? ─────────────────────────────────────
# 3. cell --> for each paw, its "bar height" in a session is the median height
# (z) over that session's VALID events at detection. Each event is then judged
# by how much higher than that each paw is, median over POSTURE_WIN.
_pw = (T >= POSTURE_WIN[0]) & (T <= POSTURE_WIN[1])
for k in KEYS:
    s, ev = PP[k], EV[k]
    val = np.isin(ev["code"], (0, 1))
    for p in ("dom", "nondom"):
        Z   = np.asarray(s[f"kin_{p}"], float)[:, 2, :][:, _pw]        # height
        trk = np.isfinite(Z).mean(axis=1)
        with np.errstate(all="ignore"):
            h = np.nanmedian(Z, axis=1)
        h[trk < MIN_TRACKED] = np.nan
        ref = np.nanmedian(h[val])
        ev[f"h_{p}"]   = h - ref                     # mm above the bar height
        ev[f"ref_{p}"] = ref
    cls = np.full(ev["code"].size, "bar", dtype=object)
    hd, hn = ev["h_dom"], ev["h_nondom"]
    cls[(hd > H_UP) & (hn <= H_UP)] = "reach-like"
    cls[(hd > H_UP) & (hn > H_UP)] = "both off"
    cls[~np.isfinite(hd) | ~np.isfinite(hn)] = "untracked"
    ev["cls"] = cls
    # keep: valid always (reported, never removed); invalid only when on the bar
    ev["keep"] = val | (cls == "bar") | ((cls == "untracked") & KEEP_UNTRACKED)
    # the flag, for every event: valid ones carry it but are still kept
    ev["reach_like"] = np.isin(cls, ("reach-like", "both off"))
    # the independent check: a reach detected nearby?
    rt, tt = REACH_T[k], np.asarray(s["t"], float)
    ev["reach_near"] = np.array([np.any(np.abs(rt - t) <= REACH_TOL) for t in tt])

CLASSES = ("bar", "reach-like", "both off", "untracked")
print(f"\nheight of the paws at detection ({POSTURE_WIN} s re the timestamp; "
      f"up = > {H_UP:.0f} mm above the session's push/pull height)")
print(f"{'':14s} {'n':>5s} " + " ".join(f"{c:>11s}" for c in CLASSES)
      + f"   {'reach within ' + str(REACH_TOL) + ' s: bar / reach-like':>34s}")
for c in CODE:
    cl = np.concatenate([EV[k]["cls"][EV[k]["code"] == c] for k in KEYS])
    rn = np.concatenate([EV[k]["reach_near"][EV[k]["code"] == c] for k in KEYS])
    if not cl.size: continue
    frac = lambda m: f"{100*np.mean(rn[m]):.0f}%" if m.any() else "-"
    print(f"{CODE[c]:14s} {cl.size:5d} "
          + " ".join(f"{100*np.mean(cl == x):10.0f}%" for x in CLASSES)
          + f"   {frac(cl == 'bar'):>16s} / {frac(cl == 'reach-like'):<8s}")
print("  valid events are FLAGGED (reach_like), never removed: they started a trial "
      "whatever the paws did.")
print("  for them the 'reach-like' rate is also the most the rule can be wrong by on "
      "real push/pulls.")

print(f"\nper mouse, INVALID events kept / total:")
for m in sorted({k[0] for k in KEYS}):
    ks = [k for k in KEYS if k[0] == m]
    inv = np.concatenate([np.isin(EV[k]["code"], (2, 3)) for k in ks])
    kp  = np.concatenate([EV[k]["keep"] for k in ks])
    rl  = np.concatenate([EV[k]["cls"] == "reach-like" for k in ks])
    print(f"  {m:18s} {int((kp & inv).sum()):4d} / {int(inv.sum()):4d}   "
          f"(reach-like {int((rl & inv).sum())})")


# %% ── 3a. the same numbers as tables (CSV, in the figure folder) ───────────
# 3a. cell --> one row per event type (all mice), and one per mouse x event
# type. Counts AND percentages, so they can be quoted directly.
#   kept      = what is used downstream (valid: always; invalid: on the bar)
#   flagged   = reach-like + both off (for valid events: flagged but kept)
#   reach_near_* = % of that group with a reach detected within REACH_TOL --
#                  the independent check, not used to classify
import pandas as pd


def _row(sel_keys, c):
    cl = np.concatenate([EV[k]["cls"][EV[k]["code"] == c] for k in sel_keys])
    kp = np.concatenate([EV[k]["keep"][EV[k]["code"] == c] for k in sel_keys])
    rn = np.concatenate([EV[k]["reach_near"][EV[k]["code"] == c] for k in sel_keys])
    n  = int(cl.size)
    if not n: return None
    r  = dict(event=CODE[c], n=n)
    for x in CLASSES:
        r[f"n_{x.replace(' ', '_').replace('-', '_')}"] = int((cl == x).sum())
    fl = np.isin(cl, ("reach-like", "both off"))
    r.update(n_flagged=int(fl.sum()), n_kept=int(kp.sum()), n_excluded=int(n - kp.sum()))
    for x in CLASSES:
        r[f"pct_{x.replace(' ', '_').replace('-', '_')}"] = round(100 * float(np.mean(cl == x)), 1)
    r.update(pct_flagged=round(100 * float(fl.mean()), 1),
             pct_kept=round(100 * float(kp.mean()), 1),
             pct_excluded=round(100 * float(1 - kp.mean()), 1),
             pct_reach_near_bar=(round(100 * float(rn[cl == "bar"].mean()), 1)
                                 if (cl == "bar").any() else np.nan),
             pct_reach_near_flagged=(round(100 * float(rn[fl].mean()), 1)
                                     if fl.any() else np.nan))
    return r


TAB_ALL = pd.DataFrame([r for c in CODE if (r := _row(KEYS, c)) is not None])
TAB_MOUSE = pd.DataFrame([dict(mouse=m, **r) for m in sorted({k[0] for k in KEYS})
                          for c in CODE
                          if (r := _row([k for k in KEYS if k[0] == m], c)) is not None])
_note = (f"# 21_pushpull_clean_align.py | rule: dom paw > {H_UP} mm above its push/pull "
         f"height in {POSTURE_WIN} s re the timestamp; valid events are never excluded "
         f"(flagged only); reach_near = a reach within {REACH_TOL} s" + chr(10))
if SAVE_TABLES:
    for nm, tab in (("cleaning_summary.csv", TAB_ALL), ("cleaning_per_mouse.csv", TAB_MOUSE)):
        _p = config.fig_dir(OUT_SUBDIR) / nm
        with open(_p, "w", encoding="utf-8", newline="") as fh:
            fh.write(_note); tab.to_csv(fh, index=False)
        print(f"table -> {_p}")
print(TAB_ALL[["event", "n", "n_kept", "n_excluded", "pct_kept", "pct_excluded",
               "pct_flagged", "pct_reach_near_bar", "pct_reach_near_flagged"]]
      .to_string(index=False))


# %% ── 3b. FIGURE: the rule ─────────────────────────────────────────────────
# Each dot is an event: how high the dom paw (x) and the non-dom paw (y) were
# above their push/pull height at detection. The rule is the two dashed lines.
# Valid events should sit in the bottom-left corner; invalid events detected
# during a reach should sit bottom-RIGHT (dom up, non-dom on the bar).
fig, axes = plt.subplots(2, 4, figsize=(20, 9),
                         gridspec_kw=dict(height_ratios=[2.2, 1]))
_lim = (-3, 15)
_hb = np.arange(_lim[0], _lim[1] + 0.25, 0.5)
for col, c in enumerate(CODE):
    ax = axes[0, col]
    dd = np.concatenate([EV[k]["h_dom"][EV[k]["code"] == c] for k in KEYS])
    dn = np.concatenate([EV[k]["h_nondom"][EV[k]["code"] == c] for k in KEYS])
    cl = np.concatenate([EV[k]["cls"][EV[k]["code"] == c] for k in KEYS])
    for x in CLASSES[:3]:
        m = cl == x
        ax.plot(np.clip(dd[m], *_lim), np.clip(dn[m], *_lim), "o", ms=3,
                alpha=0.35, mew=0, color=CLS_CLR[x],
                label=f"{x} {100*m.mean():.0f}%")
    ax.axvline(H_UP, color="k", lw=1, ls="--"); ax.axhline(H_UP, color="k", lw=1, ls="--")
    ax.set(xlim=_lim, ylim=_lim, aspect="equal",
           title=f"{CODE[c]} (n={cl.size}, {100*np.mean(cl == 'untracked'):.0f}% untracked)")
    if col == 0: ax.set_ylabel("NON-DOM paw height above its push/pull height (mm)")
    ax.title.set_fontsize(10)
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=7, loc="upper right")
    # the same dom heights as a histogram: the gap the threshold sits in
    ax = axes[1, col]
    ax.hist(np.clip(dd[np.isfinite(dd)], *_lim), bins=_hb, color=CLR[c],
            weights=np.full(np.isfinite(dd).sum(), 1 / max(np.isfinite(dd).sum(), 1)))
    ax.axvline(H_UP, color="k", lw=1, ls="--")
    ax.set(xlim=_lim, xlabel="DOM paw height above its push/pull height (mm)",
           ylabel="proportion" if col == 0 else "")
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"CLEANING 1 / paw HEIGHT at detection ({POSTURE_WIN} s re the timestamp), "
             f"relative to each paw's height during that session's VALID events "
             f"(clipped to {_lim} mm)\nbottom-right = dom paw UP, non-dom on the bar "
             f"= detected while REACHING. Bottom row: the dom heights alone", fontsize=11)
plt.tight_layout(rect=(0, 0, 1, 0.9), h_pad=2.5)
_save(fig, "clean1_posture_rule.png")


# %% ── 3c. FIGURE: what the classes look like ───────────────────────────────
# Mean paw HEIGHT (z), paw FORWARD position (x) and force, for valid events,
# invalid events kept and invalid events removed. Both paws share an axis:
# dom paw solid, non-dom paw dashed. A removed "reach-like" event should show
# the dom paw UP and FORWARD while the non-dom stays down; a kept one should
# look like a valid push/pull.
from matplotlib.lines import Line2D
_sm = (T >= PLOT_WIN[0]) & (T <= PLOT_WIN[1]); TS = T[_sm]


def _mean_trace(c, cls, what):
    """mean +- SEM over the events of code c (and class cls, if given) of
    what = ("dom"|"nondom", axis index) or "force"."""
    A = []
    for k in KEYS:
        m = EV[k]["code"] == c
        if cls is not None: m &= EV[k]["cls"] == cls
        if not m.any(): continue
        A.append(np.asarray(PP[k]["force"], float)[m][:, _sm] if what == "force"
                 else np.asarray(PP[k][f"kin_{what[0]}"], float)[m][:, what[1], :][:, _sm])
    if not A: return None
    A = np.concatenate(A)
    return (np.nanmean(A, axis=0),
            np.nanstd(A, axis=0) / np.sqrt(np.maximum(np.isfinite(A).sum(0), 1)), A.shape[0])


fig, axes = plt.subplots(2, 3, figsize=(18, 8.6), sharex=True)
for r, (codes, nm) in enumerate((((0, 2), "PUSH"), ((1, 3), "PULL"))):
    groups = [(f"{CODE[codes[0]]}", codes[0], None, CLR[codes[0]]),
              (f"{CODE[codes[1]]}, kept", codes[1], "bar", "0.25"),
              (f"{CODE[codes[1]]}, reach-like (removed)", codes[1], "reach-like",
               CLS_CLR["reach-like"])]
    for col, (ax_i, lab) in enumerate(((2, "paw height z (mm)"), (0, "paw forward position x (mm)"),
                                       (None, "load cell (mV/V)"))):
        ax = axes[r, col]
        for g, c, x, clr in groups:
            if ax_i is None:
                res = _mean_trace(c, x, "force")
                if res is None: continue
                mn, se, n = res
                ax.plot(TS, mn, color=clr, lw=2, label=f"{g} (n={n})")
                ax.fill_between(TS, mn - se, mn + se, color=clr, alpha=0.18, lw=0)
                continue
            for paw, ls, lw in (("dom", "-", 2.0), ("nondom", "--", 1.5)):
                res = _mean_trace(c, x, (paw, ax_i))
                if res is None: continue
                mn, se, n = res
                ax.plot(TS, mn, color=clr, lw=lw, ls=ls,
                        label=f"{g} (n={n})" if paw == "dom" else None)
                ax.fill_between(TS, mn - se, mn + se, color=clr, alpha=0.15, lw=0)
        ax.axvspan(*POSTURE_WIN, color="0.92", zorder=0, lw=0)
        ax.axvline(0, color="orange", lw=1.4, ls="--")
        ax.set(ylabel=lab, title=f"{nm}: {lab.split(' (')[0]}")
        ax.title.set_fontsize(10)
        if r == 1: ax.set_xlabel("time from the timestamp / cue (s)")
        ax.spines[["top", "right"]].set_visible(False)
        h, l = ax.get_legend_handles_labels()
        if ax_i is not None:              # say which line style is which paw
            h += [Line2D([], [], color="k", lw=2, ls="-"),
                  Line2D([], [], color="k", lw=1.5, ls="--")]
            l += ["dom paw", "non-dom paw"]
        ax.legend(h, l, fontsize=7)
fig.suptitle("CLEANING 2 / valid events, invalid events kept and invalid events removed — "
             "dom paw solid, non-dom paw dashed; grey = POSTURE_WIN, where the rule looks\n"
             "a kept invalid event should look like a valid one; a reach-like one has the "
             "dom paw UP and FORWARD while the non-dom paw stays down", fontsize=11)
plt.tight_layout(rect=(0, 0, 1, 0.92))
_save(fig, "clean2_class_means.png")


# %% ── 3d. FIGURE: single events, kept and removed ─────────────────────────
# Round-robin over mice, so the examples are not all from one animal. Each panel:
# dom paw height (colour), non-dom paw height (dashed) on the left axis, force
# (grey) on the right axis. Onset (green) and REF_FORCE crossing (blue) marked.
def _pick(cond, n, seed=0):
    rng, pool = np.random.default_rng(seed), {}
    for k in KEYS:
        for e in np.flatnonzero(cond(EV[k])): pool.setdefault(k[0], []).append((k, int(e)))
    for m in pool: pool[m] = [pool[m][j] for j in rng.permutation(len(pool[m]))]
    out, r = [], 0
    while len(out) < n and any(len(v) > r for v in pool.values()):
        for m in sorted(pool):
            if len(pool[m]) > r and len(out) < n: out.append(pool[m][r])
        r += 1
    return out


COLS = [("invalid, KEPT (on the bar)", lambda ev: np.isin(ev["code"], (2, 3)) & (ev["cls"] == "bar")),
        ("invalid, REMOVED (reach-like)", lambda ev: np.isin(ev["code"], (2, 3)) & (ev["cls"] == "reach-like")),
        ("VALID, reach-like: FLAGGED but KEPT (they started a trial)",
         lambda ev: np.isin(ev["code"], (0, 1)) & (ev["cls"] == "reach-like"))]
EXS = [_pick(f, EX_EVENTS) for _, f in COLS]
fig, axes = plt.subplots(EX_EVENTS, len(COLS), figsize=(6.6 * len(COLS), 2.0 * EX_EVENTS),
                         squeeze=False, sharex=True)
for col, ((ttl, _), ex) in enumerate(zip(COLS, EXS)):
    for r in range(EX_EVENTS):
        ax = axes[r, col]
        if r >= len(ex): ax.axis("off"); continue
        k, e = ex[r]; c = int(EV[k]["code"][e])
        ax.plot(TS, np.asarray(PP[k]["kin_dom"], float)[e, 2, _sm], color=CLR[c], lw=1.7)
        ax.plot(TS, np.asarray(PP[k]["kin_nondom"], float)[e, 2, _sm], color=CLR[c],
                lw=1.1, ls="--")
        ax.set_ylabel(f"{k[0].split('_')[-1]}|{k[1]} #{e}\n{CODE[c]}\nz (mm)", fontsize=6.5)
        ax.tick_params(labelsize=6.5)
        a2 = ax.twinx()
        a2.plot(TS, np.asarray(PP[k]["force"], float)[e, _sm], color="0.5", lw=1)
        a2.axhline(0, color="0.8", lw=0.5); a2.tick_params(labelsize=6, colors="0.5")
        ax.axvspan(*POSTURE_WIN, color="0.93", zorder=0, lw=0)
        ax.axvline(0, color="orange", lw=1.3, ls="--")
        if np.isfinite(EV[k]["t_cross"][e]): ax.axvline(EV[k]["t_cross"][e], color="tab:blue", lw=1)
        if np.isfinite(EV[k]["t_onset"][e]): ax.axvline(EV[k]["t_onset"][e], color="tab:green", lw=1.2)
        ax.text(0.99, 0.95, f"up: dom {EV[k]['h_dom'][e]:+.1f} mm, non-dom "
                            f"{EV[k]['h_nondom'][e]:+.1f} mm"
                            + ("  | reach nearby" if EV[k]["reach_near"][e] else ""),
                transform=ax.transAxes, ha="right", va="top", fontsize=6.5, color="0.3")
        if r == 0: ax.set_title(ttl, fontsize=10)
        if r == EX_EVENTS - 1: ax.set_xlabel("time from the timestamp / cue (s)")
        ax.spines[["top"]].set_visible(False)
fig.suptitle("CLEANING 3 / single events — dom paw height (solid), non-dom height "
             "(dashed), force (grey, right axis); orange = timestamp/cue, blue = "
             f"crossing of {REF_FORCE:.0f}, green = force onset", fontsize=11)
plt.tight_layout(rect=(0, 0, 1, 0.96))
_save(fig, "clean3_single_events.png")


# %% ── 4. ALIGNMENT: cue vs force onset ─────────────────────────────────────
# 4. cell --> left: the force of the KEPT events aligned on the CUE (timestamp);
# valid and invalid rise at different times because detection differs.
# Middle: the same events aligned on the FORCE ONSET; if the onset is the shared
# phase, valid and invalid of the same direction now rise together.
# Right: where the onset falls relative to the cue, per code.
def _realign(A, shift):
    """A (n, T) on T; shift (n,) s -> A re-sampled so that t=0 is at `shift`."""
    out = np.full((A.shape[0], TS.size), np.nan)
    for i in range(A.shape[0]):
        if np.isfinite(shift[i]):
            out[i] = np.interp(TS + shift[i], T, A[i], left=np.nan, right=np.nan)
    return out


fig, axes = plt.subplots(2, 3, figsize=(18, 8.6))
for r, codes in enumerate(((0, 2), (1, 3))):
    for col, mode in enumerate(("cue", "onset")):
        ax = axes[r, col]
        for c in codes:
            A, S = [], []
            for k in KEYS:
                m = (EV[k]["code"] == c) & EV[k]["keep"] & np.isfinite(EV[k]["t_onset"])
                if not m.any(): continue
                A.append(np.asarray(PP[k]["force"], float)[m]); S.append(EV[k]["t_onset"][m])
            if not A: continue
            A, S = np.concatenate(A), np.concatenate(S)
            Y = A[:, _sm] if mode == "cue" else _realign(A, S)
            mn = np.nanmean(Y, axis=0)
            se = np.nanstd(Y, axis=0) / np.sqrt(np.maximum(np.isfinite(Y).sum(0), 1))
            ax.plot(TS, mn, color=CLR[c], lw=2.2, label=f"{CODE[c]} (n={Y.shape[0]})")
            ax.fill_between(TS, mn - se, mn + se, color=CLR[c], alpha=0.2, lw=0)
        ax.axvline(0, color="orange" if mode == "cue" else "tab:green", lw=1.5, ls="--")
        ax.axhline(0, color="k", lw=0.5)
        ax.set(xlabel=f"time from the {'CUE (timestamp)' if mode == 'cue' else 'FORCE ONSET'} (s)",
               ylabel="load cell (mV/V)",
               title=f"aligned on the {'cue' if mode == 'cue' else 'force onset'}")
        ax.title.set_fontsize(10)
        ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)
    ax = axes[r, 2]
    _b = np.arange(-600, 101, 10)
    for c in codes:
        on = np.concatenate([EV[k]["t_onset"][(EV[k]["code"] == c) & EV[k]["keep"]] for k in KEYS])
        on = on[np.isfinite(on)] * 1000
        if not on.size: continue
        ax.hist(on, bins=_b, weights=np.full(on.size, 1 / on.size), color=CLR[c], alpha=0.6,
                label=f"{CODE[c]}: median {np.median(on):+.0f} ms")
    ax.axvline(0, color="orange", lw=1.5, ls="--")
    ax.set(xlabel="force onset re the cue (ms)", ylabel="proportion of events",
           title="where the onset falls relative to the cue")
    ax.title.set_fontsize(10)
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)
fig.suptitle(f"ALIGNMENT / KEPT events only — onset = force first exceeds "
             f"{ONSET_FRAC:.0%} of the way from its baseline to {REF_FORCE:.0f} mV/V, "
             f"on the rise that crosses it\naligned on the onset, valid and invalid of "
             f"the same direction should rise together", fontsize=11)
plt.tight_layout(rect=(0, 0, 1, 0.92))
_save(fig, "align_cue_vs_onset.png")


# %% ── 5. save ──────────────────────────────────────────────────────────────
# 5. cell --> per session, aligned index-for-index with 16's events (and so with
# st_init): keep flag, class, and the onset as a SHIFT re the timestamp. To align
# anything on the onset, subtract `onset_shift` from its time axis.
if SAVE_VARS:
    out = dict(
        params=dict(REF_FORCE=REF_FORCE, ANCHOR_WIN=ANCHOR_WIN, BASE_WIN=BASE_WIN,
                    ONSET_FRAC=ONSET_FRAC, MAX_RISE=MAX_RISE, POSTURE_WIN=POSTURE_WIN,
                    H_UP=H_UP, MIN_TRACKED=MIN_TRACKED, KEEP_UNTRACKED=KEEP_UNTRACKED),
        source=PP_PKL, codes=CODE,
        sessions={k: dict(code=EV[k]["code"], keep=EV[k]["keep"], cls=EV[k]["cls"],
                          reach_like=EV[k]["reach_like"],
                          onset_shift=EV[k]["t_onset"], cross_shift=EV[k]["t_cross"],
                          h_dom=EV[k]["h_dom"], h_nondom=EV[k]["h_nondom"],
                          reach_near=EV[k]["reach_near"], t=np.asarray(PP[k]["t"]))
                  for k in KEYS})
    _p = config.var_path(OUT_PKL)
    with open(_p, "wb") as fh: pickle.dump(out, fh, protocol=4)
    print(f"saved -> {_p}")
    print("  keep = valid always + invalid on the bar; reach_like = the flag, for "
          "valid AND invalid")
    print("  sessions[(mouse, sess)] = dict(code, keep, cls, reach_like, onset_shift (s re the "
          "timestamp), cross_shift, h_dom, h_nondom (mm above the bar), "
          "reach_near, t)")

# %%
