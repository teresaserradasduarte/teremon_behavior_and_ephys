# %% ── 0. parameters ─────────────────────────────────────────────────────────
# STEP 04 -- the automatic cerebellar-cortex label, from C4 (step 02) and the
# CS-triggered pause (step 03). Every neuron gets
#   auto_cell_type   PkC_ss / PkC_cs / MLI / GoC / MFB / unlabelled
#   auto_source      which rule gave the label (below)
#   auto_conf        high / low
#   flags            things to look at in the GUI (never change the label)
#   partner          uid of its CS (for an SS) or of its SS (for a CS), if a pair passed
#
# RULES, applied in this order (a later rule can overwrite an earlier one):
#   1  c4          C4 reported it (ratio >= 2) AND it passes C4_MIN_RATIO /
#                  C4_MIN_VOTES of its class.                           -> high
#                  With the defaults this reproduces s04 exactly.
#   2  c4_relaxed  C4 ran but did not report it (ratio < 2), or s04 dropped it:
#                  accepted if ratio >= RELAX_MIN_RATIO and votes >=
#                  RELAX_MIN_VOTES of its class (MFB: 0.94, as in s04). -> low
#                  OFF while RELAX_MIN_RATIO is None -- choose it from the
#                  figures of step 02.
#   3  cs_pause    it pauses after the complex spikes of a nearby unit (a
#                  passing pair of step 03)  -> PkC_ss                  -> high
#                  This overrides C4, except when the unit is itself a CS (labelled
#                  PkC_cs, or C4's best guess is PkC_cs and it fires <= CS_FR_MAX):
#                  then the label is kept and the unit is flagged
#                  (two pieces of one CS, or a mis-sorted unit).
#   4  cs_pair     an unlabelled unit whose C4 best guess is PkC_cs, that has a
#                  pausing SS partner and fires <= CS_FR_MAX  -> PkC_cs -> high
# FLAGS
#   conflict_c4_<type>   cs_pause relabelled a unit C4 had reported as <type>
#   cs_pauses_after_cs   a PkC_cs unit pauses after another CS
#   cs_fr_high           PkC_cs firing above CS_FR_MAX (complex spikes are ~1 Hz)
#   ss_fr_low            PkC_ss firing below SS_FR_MIN
#   has_pair             its CS/SS partner was found (a positive check)
import collections
import numpy as np
import matplotlib.pyplot as plt
import ct_config as C, ct_io as io

TYPES = ["PkC_ss", "PkC_cs", "MLI", "GoC", "MFB"]
C4_MIN_RATIO = {t: 2.0 for t in TYPES}                   # rule 1
C4_MIN_VOTES = {"PkC_ss": 0.0, "PkC_cs": 0.0, "MLI": 0.0, "GoC": 0.0, "MFB": 0.94}
RELAX_MIN_RATIO = None      # rule 2: e.g. 1.5 ; None = rule off
RELAX_MIN_VOTES = {"PkC_ss": 0.8, "PkC_cs": 0.8, "MLI": 0.8, "GoC": 0.8, "MFB": 0.94}
#                             MFB stays at s04's 0.94: C4 tends to call unclear cells MFB
USE_CS_PAUSE = True         # rules 3 and 4
CS_FR_MAX = 5.0             # Hz
SS_FR_MIN = 10.0            # Hz

SHOW_FIGS = True
SAVE_FIGS = True
SAVE_VARS = True            # writes cb_auto.csv


# %% ── 1. apply the rules ────────────────────────────────────────────────────
FEAT = io.load_features()
FR = {r["uid"]: r["fr"] for r in FEAT["rows"]}
IDX = {r["uid"]: i for i, r in enumerate(FEAT["rows"])}
C4 = io.read_csv(C.CB_C4_CSV)
if not C4: raise FileNotFoundError(f"{C.CB_C4_CSV} not found -- run 02_CB_c4_review.py")
PAIRS = io.load_pkl(C.CB_PAIRS_PKL) if (USE_CS_PAUSE and C.CB_PAIRS_PKL.is_file()) else None
if USE_CS_PAUSE and PAIRS is None:
    print(f"!! {C.CB_PAIRS_PKL.name} not found: rules 3-4 skipped (run 03_CB_cs_ss_pairs.py)")
best_cs_of = PAIRS["best_cs_of"] if PAIRS else {}
best_ss_of = PAIRS["best_ss_of"] if PAIRS else {}

rows = []
for o in C4:
    u, pred = o["uid"], o["c4_pred"]
    ratio, votes = io.fnum(o["c4_ratio"]), io.fnum(o["c4_votes"])
    lab, src, conf, flags, partner = C.UNLAB, "", "", [], ""
    if (o["c4_status"] == "reported" and pred in TYPES
            and ratio >= C4_MIN_RATIO[pred] and not votes < C4_MIN_VOTES[pred]):
        lab, src, conf = pred, "c4", "high"                               # rule 1
    elif (RELAX_MIN_RATIO is not None and pred in TYPES
          and ratio >= RELAX_MIN_RATIO and votes >= RELAX_MIN_VOTES[pred]):
        lab, src, conf = pred, "c4_relaxed", "low"                        # rule 2
    if u in best_cs_of:                                                   # rule 3
        partner = best_cs_of[u]
        if lab == "PkC_cs" or (pred == "PkC_cs" and FR.get(u, np.inf) <= CS_FR_MAX):
            flags.append("cs_pauses_after_cs")
        else:
            if lab not in (C.UNLAB, "PkC_ss"): flags.append(f"conflict_c4_{lab}")
            lab, src, conf = "PkC_ss", "cs_pause", "high"
            flags.append("has_pair")
    if u in best_ss_of and "cs_pauses_after_cs" not in flags:             # rule 4
        partner = partner or best_ss_of[u]
        if lab == C.UNLAB and pred == "PkC_cs" and FR.get(u, np.inf) <= CS_FR_MAX:
            lab, src, conf = "PkC_cs", "cs_pair", "high"
        if lab == "PkC_cs": flags.append("has_pair")
    if lab == "PkC_cs" and FR.get(u, 0) > CS_FR_MAX: flags.append("cs_fr_high")
    if lab == "PkC_ss" and FR.get(u, np.inf) < SS_FR_MIN: flags.append("ss_fr_low")
    rows.append(dict(uid=u, mouse=o["mouse"], sess=o["sess"], phyID=o["phyID"],
                     old_cell_type=o["old_cell_type"], c4_status=o["c4_status"], c4_pred=pred,
                     auto_cell_type=lab, auto_source=src, auto_conf=conf,
                     flags=";".join(flags), partner=partner))

labs = TYPES + [C.UNLAB]
T = collections.Counter((r["old_cell_type"], r["auto_cell_type"]) for r in rows)
print("rows: s04 label   columns: new automatic label")
print(f"{'':12s}" + "".join(f"{t:>11s}" for t in labs) + f"{'total':>8s}")
for a in labs:
    print(f"{a:12s}" + "".join(f"{T.get((a, b), 0):11d}" for b in labs)
          + f"{sum(T.get((a, b), 0) for b in labs):8d}")
print(f"{'total':12s}" + "".join(f"{sum(T.get((a, b), 0) for a in labs):11d}" for b in labs))
print("by rule :", dict(collections.Counter(r["auto_source"] or "none" for r in rows)))
print("flags   :", dict(collections.Counter(f for r in rows for f in r["flags"].split(";") if f)))
if SAVE_VARS:
    io.write_csv(C.CB_AUTO_CSV, rows); print(f"saved -> {C.CB_AUTO_CSV}")


# %% ── 2. figure: s04 label -> new label ─────────────────────────────────────
fdir = C.fig_dir("04_CB_relabel") if SAVE_FIGS else None
M = np.array([[T.get((a, b), 0) for b in labs] for a in labs])
fig, axes = plt.subplots(1, 3, figsize=(19, 5.2), gridspec_kw=dict(width_ratios=[1.1, 1, 1.2]))
ax = axes[0]
ax.imshow(np.log1p(M), cmap="Blues")
for i in range(len(labs)):
    for j in range(len(labs)):
        if M[i, j]: ax.text(j, i, M[i, j], ha="center", va="center", fontsize=9,
                            color="w" if M[i, j] > M.max() / 3 else "k",
                            fontweight="bold" if i != j else "normal")
ax.set_xticks(range(len(labs)), labs, rotation=30, ha="right"); ax.set_yticks(range(len(labs)), labs)
ax.set(xlabel="new automatic label", ylabel="s04 label", title="off-diagonal = changed")
ax = axes[1]
x = np.arange(len(labs))
ax.bar(x - 0.2, M.sum(1), 0.4, color=[C.COLORS[l] for l in labs], alpha=0.45, edgecolor="k", label="s04 (pale)")
ax.bar(x + 0.2, M.sum(0), 0.4, color=[C.COLORS[l] for l in labs], edgecolor="k", label="new")
ax.set_xticks(x, labs, rotation=30, ha="right"); ax.set_ylabel("neurons"); ax.legend(frameon=False, fontsize=8)
ax.set_title("neurons per type")
ax = axes[2]                                              # by mouse
mice = sorted({r["mouse"] for r in rows}); bottom = np.zeros(len(mice))
for l in labs:
    h = np.array([sum(r["mouse"] == m and r["auto_cell_type"] == l for r in rows) for m in mice], float)
    ax.bar([m[2:] for m in mice], h, bottom=bottom, color=C.COLORS[l], edgecolor="w", label=l)
    bottom += h
ax.set_ylabel("neurons"); ax.set_title("new label per mouse"); ax.legend(frameon=False, fontsize=7)
for a in axes[1:]: a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "old_vs_new.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)


# %% ── 3. figure: does each type look like itself? ───────────────────────────
# Per new type: mean normalised waveform, mean ACG (each neuron scaled to its own
# 100-400 ms baseline) and the firing-rate distribution. Thin lines = neurons
# labelled by a rule other than plain C4 (the ones this pipeline added).
A = FEAT["arrays"]; t_ms = FEAT["meta"]["wf_norm_t_ms"]
lag = (np.arange(A["acg"].shape[1]) + 0.5) * FEAT["meta"]["acg_bin"] * 1e3
fig, axes = plt.subplots(3, len(labs), figsize=(3.3 * len(labs), 8.5))
for j, l in enumerate(labs):
    sel = [r for r in rows if r["auto_cell_type"] == l and r["uid"] in IDX]
    ii = np.array([IDX[r["uid"]] for r in sel], int)
    new = np.array([r["auto_source"] != "c4" and l != C.UNLAB for r in sel], bool)
    c = C.COLORS[l]
    if ii.size:
        W = A["wf_norm"][ii]; G = A["acg"][ii].astype(float)
        G = G / np.maximum(G[:, 100:400].mean(axis=1, keepdims=True), 1e-9)
        for k in np.flatnonzero(new)[:40]:
            axes[0, j].plot(t_ms, W[k], color=c, lw=0.4, alpha=0.5)
            axes[1, j].plot(lag[:100], G[k, :100], color=c, lw=0.4, alpha=0.5)
        axes[0, j].plot(t_ms, np.nanmean(W, 0), color="k" if l != "PkC_cs" else "r", lw=2)
        axes[1, j].plot(lag[:100], np.nanmean(G, 0)[:100], color="k" if l != "PkC_cs" else "r", lw=2)
        fr = np.array([FR[r["uid"]] for r in sel], float); fr = fr[np.isfinite(fr) & (fr > 0)]
        axes[2, j].hist(fr, np.geomspace(0.1, 300, 40), color=c)
    axes[2, j].set_xscale("log")
    axes[0, j].set_title(f"{l}  (n={len(sel)}, {int(new.sum())} not plain C4)", fontsize=9, color=c)
    axes[0, j].set_xlabel("ms"); axes[1, j].set_xlabel("ACG lag (ms)"); axes[2, j].set_xlabel("firing rate (Hz)")
    for a in axes[:, j]: a.spines[["top", "right"]].set_visible(False); a.tick_params(labelsize=7)
axes[0, 0].set_ylabel("normalised waveform"); axes[1, 0].set_ylabel("ACG / baseline"); axes[2, 0].set_ylabel("neurons")
fig.suptitle("04 — each new type: waveform, ACG, firing rate (bold = mean)", fontsize=11)
fig.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / "types_overview.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)
# %%
