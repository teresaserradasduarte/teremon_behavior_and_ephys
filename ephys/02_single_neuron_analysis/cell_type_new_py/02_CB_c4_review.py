# %% ── 0. parameters ─────────────────────────────────────────────────────────
# STEP 02 -- what C4 said about every cerebellar-cortex neuron, INCLUDING the
# neurons it did not report.
#
# C4 (npyx/c4/predict_cell_types.py) runs an ensemble of models on every 'good'
# unit and averages their class probabilities. It then writes to the .tsv files
# ONLY the units whose confidence ratio (top probability / second probability)
# is >= 2. s04_CB_cellType.m read those .tsv files, so a neuron ended up
# unlabelled there for one of three different reasons, separated here:
#   reported         in the .tsv (ratio >= 2); s04 additionally dropped MFB with
#                    < 94% of the model votes
#   below_threshold  C4 ran, ratio < 2: it HAS a best guess and probabilities
#                    (read from cell_type_classification/ensemble_predictions...npy)
#   not_run          C4 never classified it (not 'good' when C4 ran, too few
#                    spikes, or its ensemble row could not be identified)
# Nothing is decided here: this step only gathers the numbers (cb_c4.csv) and
# shows them, so the thresholds of step 04 can be chosen.
import collections
import numpy as np
import matplotlib.pyplot as plt
import ct_config as C, ct_io as io

S04_VOTES_MFB = 0.94   # s04's criterion, drawn in the figures for reference
C4_RATIO      = 2.0    # C4's own reporting threshold

SHOW_FIGS = True
SAVE_FIGS = True
SAVE_VARS = True       # writes cb_c4.csv

CLASSES = ["PkC_ss", "PkC_cs", "MLI", "GoC", "MFB"]


# %% ── 1. gather C4's numbers per neuron ─────────────────────────────────────
FEAT = io.load_features()
CB = [r for r in FEAT["rows"] if r["region"] == "CB"]
OLD = io.old_cb_labels()
by_sess = collections.defaultdict(list)
for r in CB: by_sess[(r["mouse"], r["sess"])].append(r)

rows = []
for (mouse, sess), neus in sorted(by_sess.items()):
    folder = io.sorter_folder(mouse, sess, "CB")
    c4 = io.read_c4(folder) if folder else dict(tsv={}, ens=None)
    row_of = ({int(c): i for i, c in enumerate(c4["ens_ids"]) if c >= 0}
              if c4["ens"] is not None else {})
    n_unres = int((c4["ens_ids"] < 0).sum()) if c4["ens"] is not None else -1
    for r in neus:
        o = dict(uid=r["uid"], mouse=mouse, sess=sess, phyID=r["phyID"],
                 old_cell_type=OLD.get(r["uid"], C.UNLAB), c4_status="not_run", c4_pred="",
                 c4_prob=np.nan, c4_ratio=np.nan, c4_votes=np.nan, c4_second="",
                 **{f"p_{c}": np.nan for c in CLASSES})
        cid = r["phyID"]
        if cid in row_of:
            P = c4["ens"][row_of[cid]]                       # (n_classes, n_models)
            m = P.mean(axis=1); k = np.argsort(m)[::-1]
            o.update(c4_pred=c4["classes"][k[0]], c4_second=c4["classes"][k[1]],
                     c4_prob=float(m[k[0]]), c4_ratio=float(m[k[0]] / max(m[k[1]], 1e-12)),
                     c4_votes=float(np.mean(P.argmax(axis=0) == k[0])))
            for c, v in zip(c4["classes"], m): o[f"p_{c}"] = float(v)
            o["c4_status"] = "reported" if cid in c4["tsv"] else "below_threshold"
        elif cid in c4["tsv"]:                               # .tsv only (no ensemble file)
            t = c4["tsv"][cid]
            o.update(c4_status="reported", c4_pred=t["pred"], c4_prob=t["prob"],
                     c4_ratio=t["ratio"], c4_votes=t["votes"])
        rows.append(o)
    cnt = collections.Counter(o["c4_status"] for o in rows[-len(neus):])
    print(f"{mouse} {sess}: {len(neus):3d} CB neurons  " + "  ".join(f"{k} {v}" for k, v in cnt.items())
          + ("" if n_unres <= 0 else f"   ({n_unres} ensemble rows not identified)")
          + ("" if folder else "   NO SORTER FOLDER"))

print(f"\n{len(rows)} CB neurons")
tab = collections.Counter((o["c4_status"], o["old_cell_type"]) for o in rows)
olds = sorted({o["old_cell_type"] for o in rows})
print(f"{'':18s}" + "".join(f"{t:>11s}" for t in olds))
for st in ("reported", "below_threshold", "not_run"):
    print(f"{st:18s}" + "".join(f"{tab.get((st, t), 0):11d}" for t in olds))
if SAVE_VARS:
    io.write_csv(C.CB_C4_CSV, rows); print(f"saved -> {C.CB_C4_CSV}")


# %% ── 2. figure: where do s04's unlabelled neurons come from? ───────────────
A = lambda k, sel=rows: np.array([o[k] for o in sel], float)
fdir = C.fig_dir("02_CB_c4") if SAVE_FIGS else None
fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
ax = axes[0]; bottom = np.zeros(3); sts = ["reported", "below_threshold", "not_run"]
for t in olds:
    h = np.array([tab.get((s, t), 0) for s in sts], float)
    ax.bar(sts, h, bottom=bottom, color=C.COLORS.get(t, "0.5"), label=t, edgecolor="w")
    bottom += h
ax.set(ylabel="neurons", title="s04 label, by what C4 did with the neuron")
ax.legend(fontsize=8, frameon=False)
ax = axes[1]                                   # C4's best guess for the unlabelled ones
un = [o for o in rows if o["old_cell_type"] == C.UNLAB and o["c4_pred"]]
cc = collections.Counter((o["c4_pred"], o["c4_status"]) for o in un)
x = np.arange(len(CLASSES))
for j, st in enumerate(("reported", "below_threshold")):
    ax.bar(x + (j - 0.5) * 0.4, [cc.get((c, st), 0) for c in CLASSES], 0.4,
           color=[C.COLORS[c] for c in CLASSES], alpha=1.0 if j == 0 else 0.45,
           edgecolor="k", label=st)
ax.set_xticks(x, CLASSES)
ax.set(ylabel="neurons", title=f"C4's best guess for s04's unlabelled neurons (n={len(un)})\n"
                               f"solid: reported (dropped by s04), pale: below C4's threshold")
for a in axes: a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "unlabelled_origin.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 3. figure: confidence of every C4 prediction ─────────────────────────
# One panel per predicted class: confidence ratio vs fraction of model votes.
# Filled = kept by s04, open = unlabelled by s04. The dashed lines are C4's
# ratio threshold and s04's MFB vote threshold. Neurons just left of / below the
# lines are the candidates a relaxed threshold (step 04) would recover.
fig, axes = plt.subplots(1, len(CLASSES), figsize=(4.2 * len(CLASSES), 4.4), sharey=True)
for ax, c in zip(axes, CLASSES):
    sel = [o for o in rows if o["c4_pred"] == c]
    kept = np.array([o["old_cell_type"] != C.UNLAB for o in sel], bool)
    r, v = A("c4_ratio", sel), A("c4_votes", sel)
    ax.scatter(r[kept], v[kept], s=14, color=C.COLORS[c], alpha=0.7, lw=0, label=f"kept ({kept.sum()})")
    ax.scatter(r[~kept], v[~kept], s=16, facecolors="none", edgecolors=C.COLORS[c], lw=0.8,
               label=f"unlabelled ({(~kept).sum()})")
    ax.axvline(C4_RATIO, color="k", ls="--", lw=0.8)
    if c == "MFB": ax.axhline(S04_VOTES_MFB, color="k", ls="--", lw=0.8)
    ax.set_xscale("log"); ax.set(title=f"C4: {c} (n={len(sel)})", xlabel="confidence ratio")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("fraction of model votes")
fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "confidence_by_class.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 4. figure: how many neurons each threshold would keep ────────────────
# Survival curves: number of neurons of each predicted class with ratio >= x.
# Flat stretch = the count does not depend on the exact threshold there.
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
xs = np.geomspace(1, 50, 200)
for c in CLASSES:
    r = A("c4_ratio", [o for o in rows if o["c4_pred"] == c])
    axes[0].plot(xs, [(r >= x).sum() for x in xs], color=C.COLORS[c], lw=1.8, label=c)
    v = A("c4_votes", [o for o in rows if o["c4_pred"] == c])
    vs = np.linspace(0.3, 1, 200)
    axes[1].plot(vs, [(v >= x).sum() for x in vs], color=C.COLORS[c], lw=1.8)
axes[0].axvline(C4_RATIO, color="k", ls="--", lw=0.8); axes[0].set_xscale("log")
axes[1].axvline(S04_VOTES_MFB, color="k", ls="--", lw=0.8)
axes[0].set(xlabel="confidence ratio threshold", ylabel="neurons kept"); axes[0].legend(frameon=False, fontsize=8)
axes[1].set(xlabel="model-vote threshold")
for a in axes: a.spines[["top", "right"]].set_visible(False)
fig.suptitle("02 — neurons kept per C4 class as a function of the threshold"); fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "threshold_survival.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)
# %%
