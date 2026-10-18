# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 0. cell --> parameters
# EMPTY TRIALS ARE MIS-READ FROM THE .mat AS TWO SPIKES, AT 0 s AND AT 1 s.
#
# MATLAB is right: a trial in which the neuron fired nothing is
#        st_init{k}  =  0x1 empty double
# The problem is on the PYTHON side. A v7.3 .mat is an HDF5 file, and HDF5 has
# no zero-length dataset of this kind, so MATLAB stores an empty array as a
# small dataset holding its DIMENSIONS -- the two integers [0, 1], for "0 rows
# by 1 column" -- and flags it with the attribute MATLAB_empty = 1.
# Until 2026-10-12 mi_io._spikes_from_ref did not look at that flag. It read the
# dataset as numbers, so the SIZE of the empty array came back as its CONTENT
# (mi_io is fixed now; this script keeps a copy of the OLD reader, _read_naive,
# so the before/after can still be shown):
#
#        0x1 empty   ->   [0.0, 1.0]   ->   "a spike at 0 s and a spike at 1 s"
#
# This script shows that, measures it, and shows the correction. It changes
# nothing on disk except the figure it saves.
#
# THE CORRECTION turns those trials back into EMPTY trials. It does NOT remove
# them: the neuron keeps exactly the same number of trials, in the same order,
# so trial k still pairs with kinematics k and with label k. Only the two fake
# spike times go. Cell 2 checks that the value rule (a trial that is exactly
# [0.0, 1.0]) and the file's own MATLAB_empty flag pick out the same trials.
import os, sys
import numpy as np
import matplotlib.pyplot as plt

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
import mi_io, config

FIELD    = "st_init"        # "st_init" (push/pull) or "st_reach"
REGION   = "CP"             # the region to illustrate; CP is where it is common
SESSIONS = [("1_CoteDor", "R2"), ("2_Lindt", "R5"), ("3_Toblerone", "R1")]
                            # sessions to pool for the measurements. None = all,
                            # which takes ~2.5 min on the 9 GB file.
PLACEHOLDER = (0.0, 1.0)    # the pair that stands in for a silent trial

RASTER_N   = 60             # trials drawn in the raster
WIN        = (-0.05, 0.05)  # (s) the count window panel F uses -- a window that
                            # straddles 0, as every one of these analyses has
PETH_BIN   = 0.010
PETH_WIN   = (-1.5, 1.5)

CLR = {"CP": np.array([106, 176, 134]) / 255,
       "CB": np.array([154, 123, 193]) / 255}
BAD = np.array([200, 60, 60]) / 255          # the colour of the placeholder

SHOW_FIGS  = True
SAVE_FIGS  = True
OUT_SUBDIR = "diagnostics"
RUN_NAME   = "placeholder"


# %% ── 1. read the raw spike times ──────────────────────────────────────────
# 1. cell --> read FIELD straight from the .mat, UNCHANGED. mi_io's own readers
# are used, so this is exactly what every other script gets.
import h5py

KEYS = None if SESSIONS is None else {tuple(k) for k in SESSIONS}
REG  = config.REGION_STRUCT[REGION]


def _read_naive(f, ref):
    """mi_io._spikes_from_ref AS IT WAS before the fix: every dataset is read as
    numbers and the MATLAB_empty flag is ignored. Kept here only so this script
    can still show what the old reader returned."""
    out = []
    for tr in np.asarray(f[ref]).ravel():
        spk = np.asarray(f[tr], float).ravel()
        out.append(spk[np.isfinite(spk)])
    return out


def _is_ph(t):
    """is this trial the placeholder, and nothing else?"""
    return (t.size == 2 and t[0] == PLACEHOLDER[0] and t[1] == PLACEHOLDER[1])


RAW  = []           # [(uid, [spike arrays per trial])]  -- as the OLD reader did
FLAG = {}           # uid -> bool per trial: the file's own MATLAB_empty flag
with h5py.File(str(config.DATA_MAT), "r") as f:
    g = f[REG]
    for i in range(g[FIELD].shape[1]):
        m = mi_io._str_from_ref(f, g["mouse"][0, i])
        s = mi_io._str_from_ref(f, g["sess"][0, i])
        if KEYS is not None and (m, s) not in KEYS: continue
        pid = int(np.asarray(f[g["phyID"][0, i]], float).ravel()[0])
        uid = f"{REGION}|{m}|{s}|{pid}"
        RAW.append((uid, _read_naive(f, g[FIELD][0, i])))
        FLAG[uid] = np.array([bool(f[r].attrs.get("MATLAB_empty", 0))
                              for r in np.asarray(f[g[FIELD][0, i]]).ravel()])

n_tr = sum(len(tr) for _, tr in RAW)
PH   = {u: np.array([_is_ph(t) for t in tr]) for u, tr in RAW}
n_ph = sum(int(v.sum()) for v in PH.values())
print(f"{REGION} {FIELD}: {len(RAW)} neurons, {n_tr} trials, "
      f"{n_ph} placeholder trials ({100*n_ph/n_tr:.1f}%) in "
      f"{sum(1 for v in PH.values() if v.any())} neurons")
print(f"truly empty trials: {sum(1 for _, tr in RAW for t in tr if t.size == 0)}")


# %% ── 2. is the pair ever a real pair of spikes? ───────────────────────────
# 2. cell --> the correction must not delete data. A real trial could in
# principle contain exactly two spikes, or a spike exactly at 0.0 or 1.0. Count
# how often each of those happens, so the rule can be read as safe or not.
n_two   = sum(1 for _, tr in RAW for t in tr if t.size == 2)
n_zero  = sum(int(np.sum(t == 0.0)) for _, tr in RAW for t in tr)
n_one   = sum(int(np.sum(t == 1.0)) for _, tr in RAW for t in tr)
n_other = sum(1 for _, tr in RAW for t in tr
              if (np.any(t == 0.0) or np.any(t == 1.0)) and not _is_ph(t))
print(f"\ntrials with exactly 2 spikes      : {n_two}   "
      f"(of which the placeholder: {n_ph})")
print(f"spikes at EXACTLY t = 0.0         : {n_zero}")
print(f"spikes at EXACTLY t = 1.0         : {n_one}   "
      f"<- identical count: they come in pairs")
print(f"trials holding 0.0 or 1.0 but NOT the placeholder: {n_other}")
print("\nso the rule 'exactly [0.0, 1.0] and nothing else' picks the "
      f"{n_ph} mis-read trials and {n_other} real trials.")

# the file's own flag is the ground truth for "this trial was empty in MATLAB"
_fl = np.concatenate([FLAG[u] for u, _ in RAW])
_vl = np.concatenate([PH[u] for u, _ in RAW])
print(f"\nMATLAB_empty flag vs the value rule, over {_fl.size} trials:")
print(f"   both                    : {int((_fl & _vl).sum())}")
print(f"   flag only (rule misses) : {int((_fl & ~_vl).sum())}")
print(f"   rule only (false alarm) : {int((~_fl & _vl).sum())}")

# ── the corrected spike times ────────────────────────────────────────────────
# t[:0] is the same array with no elements: the trial stays in the list, at the
# same position, with zero spikes. Nothing is deleted and nothing is reordered.
CORR = [(u, [t[:0] if _is_ph(t) else t for t in tr]) for u, tr in RAW]
_same_n   = all(len(a) == len(b) for (_, a), (_, b) in zip(RAW, CORR))
_n_empty  = sum(1 for _, tr in CORR for t in tr if t.size == 0)
_untouched = all(np.array_equal(a, b) for (_, ta), (_, tb) in zip(RAW, CORR)
                 for a, b in zip(ta, tb) if not _is_ph(a))
print(f"\ntrials before: {n_tr}   trials after: "
      f"{sum(len(tr) for _, tr in CORR)}   "
      f"-> same count in every neuron: {_same_n}")
print(f"empty trials before: 0   after: {_n_empty}   (= the {n_ph} mis-read ones)")
print(f"every other trial identical before and after: {_untouched}")


# %% ── 3. the figure ────────────────────────────────────────────────────────
# 3. cell --> six panels: what it looks like in the data (A, B), what it does to
# the spike-time histogram (C, D), what it does to a spike COUNT (E), and what
# it does to the population PETH (F).
# "before" = the file as it is. "after" = the same thing with every exactly
# [0.0, 1.0] trial restored to empty. Nothing else differs between them.
col = CLR[REGION]

# the neuron for the raster: the one whose placeholder fraction is nearest 25%,
# so both kinds of trial are visible in one picture
_f = {u: v.mean() for u, v in PH.items() if v.any()}
EG = min(_f, key=lambda u: abs(_f[u] - 0.25))
eg_tr = dict(RAW)[EG]; eg_ph = PH[EG]
print(f"\nraster neuron: {EG}  "
      f"({100*eg_ph.mean():.0f}% of its {eg_ph.size} trials are the placeholder)")

fig = plt.figure(figsize=(22, 9.8))
gs  = fig.add_gridspec(2, 4, hspace=0.42, wspace=0.28)

# ── A: the raster ───────────────────────────────────────────────────────────
ax = fig.add_subplot(gs[0, 0])
# the placeholder trials are not spread evenly, so show the block of
# RASTER_N CONSECUTIVE trials that contains the most of them -- consecutive,
# because the point is what the trial list looks like, not a selection of it
_w  = min(RASTER_N, len(eg_tr))
_cs = np.r_[0, np.cumsum(eg_ph.astype(int))]
_i0 = int(np.argmax(_cs[_w:] - _cs[:-_w])) if len(eg_tr) > _w else 0
_sel = np.arange(_i0, _i0 + _w)
eg_corr = dict(CORR)[EG]
for r, k in enumerate(_sel):
    t = eg_tr[k]
    if eg_ph[k]:
        ax.plot(t, [r, r], "o", ms=5, color=BAD, mew=0, zorder=3)
    else:
        ax.plot(t, np.full(t.size, r), "|", ms=3.4, color=col, mew=0.8)
ax.axvline(0.0, color=BAD, lw=0.9, ls=":", alpha=0.8)
ax.axvline(1.0, color=BAD, lw=0.9, ls=":", alpha=0.8)
ax.set(xlim=(-1.6, 1.6), ylim=(-1.5, len(_sel) + 0.5),
       xlabel=f"time from the event (s), {FIELD}",
       ylabel=f"trial {_sel[0]}-{_sel[-1]}  ({len(_sel)} rows)",
       title=f"A  AS READ by python\n{EG}")
ax.title.set_fontsize(10)
ax.text(0.015, 0.015, f"red = {int(eg_ph[_sel].sum())} of these {len(_sel)} trials "
                      f"read as [0.0, 1.0]:\nthe SIZE of an empty array, "
                      f"not two spikes",
        transform=ax.transAxes, va="bottom", ha="left", fontsize=8.5, color=BAD,
        bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
ax.spines[["top", "right"]].set_visible(False)

# ── A2: the SAME rows after the correction ─────────────────────────────────
# Same trials, same row for each. A corrected trial is an EMPTY ROW: it is still
# there (grey band = "this row is a trial, and it has no spikes"), it just has
# nothing in it -- which is what your MATLAB raster shows for it.
ax = fig.add_subplot(gs[0, 1])
for r, k in enumerate(_sel):
    t = eg_corr[k]
    if t.size == 0:
        ax.axhspan(r - 0.5, r + 0.5, color="0.90", lw=0, zorder=0)
    else:
        ax.plot(t, np.full(t.size, r), "|", ms=3.4, color=col, mew=0.8)
ax.set(xlim=(-1.6, 1.6), ylim=(-1.5, len(_sel) + 0.5),
       xlabel=f"time from the event (s), {FIELD}",
       ylabel=f"trial {_sel[0]}-{_sel[-1]}  ({len(_sel)} rows)",
       title=f"A2  CORRECTED: the same {len(_sel)} rows\n"
             f"grey = a trial that is kept, and is empty")
ax.title.set_fontsize(10)
ax.text(0.015, 0.015,
        f"{len(_sel)} trials before, {len(_sel)} after\n"
        f"{int(sum(eg_corr[k].size == 0 for k in _sel))} of them now have 0 spikes",
        transform=ax.transAxes, va="bottom", ha="left", fontsize=8.5, color="0.25",
        bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
ax.spines[["top", "right"]].set_visible(False)


# ── B: spikes per trial ─────────────────────────────────────────────────────
ax = fig.add_subplot(gs[0, 2])
_n_before = np.array([t.size for _, tr in RAW for t in tr])
_n_after  = np.array([0 if _is_ph(t) else t.size for _, tr in RAW for t in tr])
_b = np.arange(-0.5, 30.5)
ax.hist(_n_before, bins=_b, color=BAD, alpha=0.75, label="as read")
ax.hist(_n_after,  bins=_b, histtype="step", lw=2, color=col,
        label="corrected (empty again)")
ax.set_yscale("log")
ax.set(xlabel="spikes in the whole trial", ylabel="trials (log)",
       title="B  as read there is NO trial with 0 spikes;\n"
             "the empty ones all sit on the bar at 2")
ax.title.set_fontsize(10)
ax.annotate("empty trials, mis-read", xy=(2, max(np.bincount(_n_before)[2], 1)),
            xytext=(9, max(np.bincount(_n_before)[2], 1) * 0.9),
            arrowprops=dict(arrowstyle="->", color=BAD), color=BAD, fontsize=9)
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

# ── C, D: the spike-time histogram at 1 ms, around each of the two values ───
_all = np.concatenate([t for _, tr in RAW for t in tr])
_keep = np.concatenate([t for _, tr in RAW for t in tr if not _is_ph(t)])
for j, (ctr, ttl) in enumerate(((0.0, "C  around t = 0"), (1.0, "D  around t = +1.000 s"))):
    ax = fig.add_subplot(gs[0, 3] if j == 0 else gs[1, 0])
    bb = np.arange(ctr - 0.01, ctr + 0.01 + 1e-9, 0.001)
    bc = bb[:-1] + 0.0005
    hb, _ = np.histogram(_all, bins=bb)
    ha, _ = np.histogram(_keep, bins=bb)
    ax.bar(bc * 1000, hb, width=0.9, color=BAD, label="as read")
    ax.bar(bc * 1000, ha, width=0.9, color=col, label="corrected")
    ax.set(xlabel=f"time from the event (ms), 1 ms bins",
           ylabel="spikes, all neurons pooled",
           title=f"{ttl}\none bin holds {hb.max()/max(np.median(hb),1):.0f}x "
                 f"the neighbouring bins")
    ax.title.set_fontsize(10)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

# ── E: what it does to a spike COUNT, which is the analysis variable ───────
ax = fig.add_subplot(gs[1, 1])
_cb = np.array([int(((t >= WIN[0]) & (t < WIN[1])).sum())
                for _, tr in RAW for t in tr])
_ca = np.array([0 if _is_ph(t) else int(((t >= WIN[0]) & (t < WIN[1])).sum())
                for _, tr in RAW for t in tr])
_bb = np.arange(-0.5, max(8, _cb.max()) + 0.5)
ax.hist(_cb, bins=_bb, color=BAD, alpha=0.75,
        label=f"as read (0 spikes in {100*np.mean(_cb == 0):.0f}% of trials)")
ax.hist(_ca, bins=_bb, histtype="step", lw=2, color=col,
        label=f"corrected (0 spikes in {100*np.mean(_ca == 0):.0f}% of trials)")
ax.set_yscale("log")
ax.set(xlabel=f"spike COUNT in [{WIN[0]*1000:+.0f}, {WIN[1]*1000:+.0f}] ms "
              f"-- the MI variable",
       ylabel="trials (log)",
       title="E  the mis-read puts a spurious count of 1 into\n"
             "every window that straddles 0")
ax.title.set_fontsize(10)
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

# ── F: the population PETH, before and after ───────────────────────────────
ax = fig.add_subplot(gs[1, 2])
_pb = np.arange(PETH_WIN[0], PETH_WIN[1] + PETH_BIN / 2, PETH_BIN)
_pc = _pb[:-1] + PETH_BIN / 2
for nm, keep_ph, c, lw in (("as read", True, BAD, 2.0),
                           ("corrected", False, col, 1.8)):
    rows = []
    for u, tr in RAW:
        use = tr if keep_ph else [t for t in tr if not _is_ph(t)]
        if not use: continue
        h, _ = np.histogram(np.concatenate(use) if use else np.array([]), bins=_pb)
        rows.append(h / (len(tr) * PETH_BIN))        # same denominator for both:
        # a silent trial still counts as a trial, it just contributes no spikes
    A = np.asarray(rows)
    ax.plot(_pc, np.nanmean(A, axis=0), color=c, lw=lw, label=nm)
ax.axvline(0, color="0.5", lw=0.8, ls=":"); ax.axvline(1.0, color="0.5", lw=0.8, ls=":")
ax.set(xlabel="time from the event (s)", ylabel="population mean rate (Hz)",
       title="F  what it does to the PETH: two fake spikes per\n"
             "empty trial, at 0 and at +1.000 s")
ax.title.set_fontsize(10)
ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

# ── G: the bookkeeping ─────────────────────────────────────────────────────
ax = fig.add_subplot(gs[1, 3]); ax.axis("off")
ax.text(0, 1, "\n".join([
    "G  what the correction does, in numbers",
    "",
    f"trials before            {n_tr:>8d}",
    f"trials after             {sum(len(tr) for _, tr in CORR):>8d}",
    f"  same count per neuron  {str(_same_n):>8s}",
    "",
    f"empty trials before      {0:>8d}",
    f"empty trials after       {_n_empty:>8d}",
    "",
    "MATLAB_empty flag vs the rule",
    f"  both                   {int((_fl & _vl).sum()):>8d}",
    f"  flag only              {int((_fl & ~_vl).sum()):>8d}",
    f"  rule only              {int((~_fl & _vl).sum()):>8d}",
    "",
    f"other trials unchanged   {str(_untouched):>8s}",
    "",
    "No trial is removed. A mis-read trial",
    "becomes an EMPTY trial in the same",
    "position, so trial k still pairs with",
    "kinematics k and label k.",
]), va="top", ha="left", fontsize=9.5, family="monospace", transform=ax.transAxes)

fig.suptitle(
    f"Empty trials mis-read as [0.0, 1.0] in {FIELD}  —  {REGION}, "
    f"{len(RAW)} neurons, {n_tr} trials"
    + ("" if SESSIONS is None else f" ({len(SESSIONS)} sessions)") + "\n"
    f"MATLAB stores a 0x1 empty array in the v7.3 file as its DIMENSIONS [0 1] "
    f"plus a MATLAB_empty flag; the python reader ignores the flag and returns "
    f"the dimensions as spike times.\nThe correction makes those trials empty "
    f"again. It does not remove them.",
    fontsize=11)
plt.tight_layout(rect=(0, 0, 1, 0.91))
if SAVE_FIGS:
    _p = config.fig_dir(OUT_SUBDIR) / f"{RUN_NAME}_{REGION}_{FIELD}.png"
    fig.savefig(_p, dpi=150, bbox_inches="tight")
    print(f"\nfigure -> {_p}")
if SHOW_FIGS: plt.show()
plt.close(fig)

# %%
