# -*- coding: utf-8 -*-
"""
sorter_compare_weak_pause.py — Would KS4 pair more CS with a looser pause threshold, and are KS4 units
"more multi-unit" so that a real SS pause is diluted?

Uses ONLY results already computed by step0 v5 (every CS x SS pair is stored in cs_ss_pairs.csv, with its
pause ratio, asymmetry, jitter p, collision dip and distance) + label/contamination files read from the sorter
folders (Dropbox: READ ONLY). Nothing is re-sorted, step0 is not rerun, no step0 decision changes.

1. Threshold sweep.  A CS is "paired" at threshold T when some SS candidate has pause ratio < T and passes the
   other step0 criteria unchanged (asym > 0.2, jitter p < 0.05, no collision dip). T = 0.7 reproduces step0.
   Plausibility: a CS and the SS of the same Purkinje cell sit close on the probe. The share of best partners
   within NEAR_UM is compared with the share of ALL CS x SS pairs within NEAR_UM (chance).
2. Contamination.  Per SS candidate: fraction of ISIs < 1 ms (step0, same window rule in both sorters) and the
   implied contamination c: with r = viol / (1 - exp(-fr * 1 ms)) (observed / fully-Poisson violations), a unit
   made of one refractory neuron + a fraction c of Poisson spikes gives r = c (2 - c), so c = 1 - sqrt(1 - r),
   clipped to [0, 0.9]. A rough estimate (it ignores bursts and the 1 ms choice). Kilosort ContamPct is also read (same KS code in
   both sorters, but reported only for information).
   Dilution model: a unit that is a mix of the true SS (fraction 1 - c) and unrelated spikes (c) shows
       pr_obs = c + (1 - c) * pr_true      ->  pr_corr = 1 - (1 - pr_obs) / (1 - c)
   so pairing is redone with pr_corr < 0.7 (other criteria unchanged).

Output: _sorter_comparison/weak_pause_check_<BATCH_TAG>/ (tables, figure) and a copy of the figure as
FINAL_KS4_vs_IBL/figs/fig7_weak_pause_and_contamination.*
"""
# %% 1. PARAMETERS + INPUTS
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

DROPBOX_ROOT = Path(r"D:\Learning Lab Dropbox")
CODE_DIR = Path(r"C:\Users\simon\Desktop\Code\EPhys\ComplexSpike_Analysis")
LOCAL = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation")
CMP = LOCAL / "_sorter_comparison"
BATCH_TAG = os.environ.get("STEP0_BATCH_TAG", "v5batch")
RUN_DIR = f"step0_explore_bhvspan_{BATCH_TAG}"

THRESHOLDS = np.round(np.arange(0.5, 0.951, 0.05), 2)
PR_STEP0 = 0.7               # step0 candidate threshold
PR_RELAXED = 0.85            # "weak pause" upper bound used in fig3 (reported in the summary table)
ASYM_MIN, JITTER_P = 0.2, 0.05
NEAR_UM = 150.0              # CS and SS of one Purkinje cell: within ~150 um on the probe
C_MAX = 0.9
N_PERM = 10000

KS, IB = "ks4", "ibl"
COL = {KS: "#2a78d6", IB: "#eb6834"}
NAME = {KS: "Kilosort4", IB: "IBL sorter"}
INK, INK2, GRID_C = "#0b0b0b", "#52514e", "#e4e3df"


def c_from_isi(viol, fr):
    r = (viol / (1 - np.exp(-fr * 1e-3))).clip(0, 1)
    return (1 - np.sqrt(1 - r)).clip(0, C_MAX)


def safe_out_dir(p):
    p = Path(p).resolve()
    for forbidden in (DROPBOX_ROOT, CODE_DIR):
        try:
            p.relative_to(forbidden.resolve())
        except ValueError:
            continue
        raise RuntimeError(f"Refusing to write under {forbidden}: {p}")
    p.mkdir(parents=True, exist_ok=True)
    return p


OUT = safe_out_dir(CMP / f"weak_pause_check_{BATCH_TAG}")
FINAL_FIG = safe_out_dir(CMP / "FINAL_KS4_vs_IBL" / "figs")

st = pd.read_csv(CMP / "session_table.csv")
P, U = [], []
for _, r in st.iterrows():
    t = LOCAL / f"{r.mouse}_{r.R}" / RUN_DIR / "tables"
    if not (t / "cs_ss_pairs.csv").exists():
        continue
    d = Path(r.sorter_dir)                       # read only
    lab = pd.read_csv(d / "cluster_KSLabel.tsv", sep="\t").set_index("cluster_id")["KSLabel"]
    con = pd.read_csv(d / "cluster_ContamPct.tsv", sep="\t").set_index("cluster_id").iloc[:, 0]
    u = pd.read_csv(t / "unit_table.csv")
    u = u.assign(mouse=r.mouse, S=r.S, sorter=r.sorter, kslabel=u.phyID.map(lab).fillna("none"),
                 ks_contam_pct=u.phyID.map(con))
    U.append(u)
    p = pd.read_csv(t / "cs_ss_pairs.csv").assign(mouse=r.mouse, S=r.S, sorter=r.sorter)
    p = p.merge(u[["phyID", "kslabel", "ks_contam_pct", "frac_isi_lt1ms", "fr"]].rename(columns={"phyID": "ss"}),
                on="ss", how="left")
    P.append(p)
U, P = pd.concat(U, ignore_index=True), pd.concat(P, ignore_index=True)
KEY = ["sorter", "mouse", "S", "cs"]
ALL_CS = P[KEY].drop_duplicates()
N_CS = ALL_CS.groupby("sorter").size()

P["ok_other"] = (P.asym > ASYM_MIN) & (P.jitter_p < JITTER_P) & (P.collision_dip != True)  # noqa: E712
P["near"] = P.dist_um < NEAR_UM
P["c_est"] = c_from_isi(P.frac_isi_lt1ms, P.fr)
P["pr_corr"] = (1 - (1 - P.pause_ratio) / (1 - P.c_est)).clip(lower=0)
assert (((P.pause_ratio < PR_STEP0) & P.ok_other) == P.candidate).all(), "T = 0.7 must reproduce step0"
print("n CS:", N_CS.to_dict())


# %% 2. THRESHOLD SWEEP + CONTAMINATION CORRECTION
def paired_at(thr, col="pause_ratio", near_only=False):
    q = P[(P[col] < thr) & P.ok_other & (P.near if near_only else True)]
    return q.sort_values(col).groupby(KEY).head(1)          # best partner of each paired CS


sweep = []
for thr in THRESHOLDS:
    b, bn = paired_at(thr), paired_at(thr, near_only=True)
    for so in (KS, IB):
        bs = b[b.sorter == so]
        sweep.append(dict(threshold=thr, sorter=so, n_CS=N_CS[so], n_paired=len(bs), frac_paired=len(bs) / N_CS[so],
                          n_paired_near=int((bn.sorter == so).sum()), frac_paired_near=(bn.sorter == so).sum() / N_CS[so],
                          share_best_near=bs.near.mean() if len(bs) else np.nan,
                          best_dist_med=bs.dist_um.median() if len(bs) else np.nan))
sweep = pd.DataFrame(sweep)
chance_near = P.groupby("sorter").near.mean()
sweep["share_all_pairs_near_chance"] = sweep.sorter.map(chance_near)
sweep.to_csv(OUT / "threshold_sweep.csv", index=False)
print(sweep.round(3).to_string(index=False))

corr = paired_at(PR_STEP0, "pr_corr")
print("contamination-corrected pr < 0.7:", corr.groupby("sorter").size().to_dict(),
      "near share", corr.groupby("sorter").near.mean().round(2).to_dict())


# per-session fractions (session = unit of replication) + within-mouse permutation
def within_mouse_perm(df, col, seed=0):
    d = df.dropna(subset=[col])
    d = d[d.groupby("mouse")["sorter"].transform("nunique") == 2]
    groups = [(g[col].to_numpy(float), (g["sorter"] == KS).to_numpy()) for _, g in d.groupby("mouse")]

    def stat(labs):
        return np.mean([np.median(v[l]) - np.median(v[~l]) for (v, _), l in zip(groups, labs)])

    obs = stat([l for _, l in groups])
    rng = np.random.default_rng(seed)
    null = np.array([stat([rng.permutation(l) for _, l in groups]) for _ in range(N_PERM)])
    return obs, (1 + np.sum(np.abs(null) >= abs(obs) - 1e-12)) / (1 + N_PERM)


ses = ALL_CS.groupby(["sorter", "mouse", "S"]).size().rename("n_CS").reset_index()
VARIANTS = {"step0 (pr < 0.7)": paired_at(PR_STEP0), "relaxed (pr < 0.85)": paired_at(PR_RELAXED),
            "relaxed (pr < 0.95)": paired_at(0.95), "contamination-corrected pr < 0.7": corr,
            "relaxed pr < 0.85, partner < 150 um": paired_at(PR_RELAXED, near_only=True)}
rows = []
for name, b in VARIANTS.items():
    col = "f_" + str(len(rows))
    ses[col] = ses.merge(b.groupby(["sorter", "mouse", "S"]).size().rename("k").reset_index(),
                         on=["sorter", "mouse", "S"], how="left")["k"].fillna(0).to_numpy() / ses.n_CS
    k, i = ses.loc[ses.sorter == KS, col], ses.loc[ses.sorter == IB, col]
    diff, p_wm = within_mouse_perm(ses, col)
    rows.append(dict(variant=name, ks4_pooled=(b.sorter == KS).sum() / N_CS[KS], ibl_pooled=(b.sorter == IB).sum() / N_CS[IB],
                     ks4_n=(b.sorter == KS).sum(), ibl_n=(b.sorter == IB).sum(), ks4_session_median=k.median(),
                     ibl_session_median=i.median(), p_mannwhitney=mannwhitneyu(k, i).pvalue, p_within_mouse=p_wm,
                     ks4_best_near=b[b.sorter == KS].near.mean(), ibl_best_near=b[b.sorter == IB].near.mean()))
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "pairing_variants_summary.csv", index=False)
print(summary.round(3).to_string(index=False))

# contamination of the SS-candidate pool (one row per unit), by sorter x tier x curated label
pool = U[U.role.isin(["SS_cand", "SS_cand2"])].copy()
pool["tier"] = np.where(pool.role == "SS_cand", 1, 2)
pool["c_est"] = c_from_isi(pool.frac_isi_lt1ms, pool.fr)
pool["group"] = pool.group.fillna("unsorted")
cont = (pool.groupby(["sorter", "tier", "group"])
        .agg(n=("phyID", "size"), isi_viol_med=("frac_isi_lt1ms", "median"), c_est_med=("c_est", "median"),
             ks_contam_pct_med=("ks_contam_pct", "median"), fr_med=("fr", "median")).reset_index())
cont.to_csv(OUT / "ss_pool_contamination.csv", index=False)
print(cont.round(3).to_string(index=False))
c_t2 = {so: pool[(pool.sorter == so) & (pool.tier == 2)].frac_isi_lt1ms for so in (KS, IB)}
p_c_t2 = mannwhitneyu(c_t2[KS], c_t2[IB]).pvalue

# best NEAR unit of every CS (the most plausible SS of the same cell): its pause and its contamination
bn_all = P[P.near].sort_values("pause_ratio").groupby(KEY).head(1)
bn_all.to_csv(OUT / "best_near_unit_per_cs.csv", index=False)
need_c = (bn_all.pause_ratio - 0.3) / 0.7      # contamination needed to turn a true pause of 0.3 into the observed one
print("best near unit: pause median", bn_all.groupby("sorter").pause_ratio.median().round(2).to_dict(),
      "| c_est median", bn_all.groupby("sorter").c_est.median().round(2).to_dict())
n_hide = int(((bn_all.sorter == KS) & ~bn_all.candidate & (bn_all.c_est >= need_c) & (bn_all.pause_ratio >= 0.7)).sum())
n_ks_np = int(((bn_all.sorter == KS) & (bn_all.pause_ratio >= 0.7)).sum())
print("unpaired KS4 CS whose best near unit has c_est >= contamination needed to hide a 0.3 pause:", n_hide, "of", n_ks_np)

# %% 3. FIGURE
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.linewidth": 0.8, "svg.fonttype": "none", "figure.facecolor": "white"})


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, color=GRID_C, lw=0.6)
    ax.set_axisbelow(True)


def title(ax, head, sub):
    ax.set_title(head + "\n" * (1 + sub.count("\n")), loc="left", fontsize=10, fontweight="bold", color=INK)
    ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=8.2, color=INK2, va="bottom")


S = summary.set_index("variant")
fig, axs = plt.subplots(1, 3, figsize=(15, 5.4), gridspec_kw={"width_ratios": [1.1, 1, 1.15]})
fig.subplots_adjust(left=0.05, right=0.98, top=0.67, bottom=0.12, wspace=0.3)
fig.suptitle("A looser pause threshold or a contamination correction does not close the gap: "
             "the units next to KS4 CS mostly do not pause", x=0.05, ha="left", fontsize=12, fontweight="bold",
             color=INK, y=0.97)
fig.text(0.05, 0.9, "Same step0 v5 pairs, re-thresholded (asymmetry > 0.2, jitter p < 0.05, no collision dip kept). "
         f"Near = partner < {NEAR_UM:.0f} µm from the CS (same Purkinje cell plausible).", fontsize=8.5, color=INK2)

ax = axs[0]
for so in (KS, IB):
    d = sweep[sweep.sorter == so]
    ax.plot(d.threshold, d.frac_paired, color=COL[so], lw=2, marker="o", ms=4, label=f"{NAME[so]}: any distance")
    ax.plot(d.threshold, d.frac_paired_near, color=COL[so], lw=2, ls=(0, (3, 2)), label=f"{NAME[so]}: partner < {NEAR_UM:.0f} µm")
    for thr in (PR_STEP0, 0.95):
        v = d.loc[d.threshold == thr, "frac_paired"].iloc[0]
        ax.text(thr + 0.008, v + 0.008, f"{v:.0%}", fontsize=7.5, color=INK, va="bottom")
ax.axvline(PR_STEP0, color=INK2, lw=0.8, ls=(0, (2, 2)))
ax.text(PR_STEP0 + 0.005, 0.012, "step0 threshold", ha="left", fontsize=7.5, color=INK2)
ax.set_xlabel("pause-ratio threshold (pair if pause ratio < threshold)")
ax.set_ylabel("fraction of CS with an SS partner (pooled)")
ax.set_ylim(0, 0.5)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
ax.legend(frameon=False, fontsize=7.5, loc="upper left", ncol=2)
clean(ax)
_gain = int(sweep[(sweep.sorter == KS) & (sweep.threshold == 0.95)].n_paired_near.iloc[0]
            - sweep[(sweep.sorter == KS) & (sweep.threshold == PR_STEP0)].n_paired_near.iloc[0])
title(ax, "A  Relaxing the pause threshold",
      f"at < 0.95: KS4 {S.at['relaxed (pr < 0.95)', 'ks4_pooled']:.0%} vs IBL {S.at['relaxed (pr < 0.95)', 'ibl_pooled']:.0%}; "
      f"near partners only: KS4 +{_gain} CS\n"
      f"within-mouse p (sessions): < 0.7: {S.at['step0 (pr < 0.7)', 'p_within_mouse']:.2g}; "
      f"< 0.85: {S.at['relaxed (pr < 0.85)', 'p_within_mouse']:.2g};\n"
      f"< 0.85 with partner < 150 µm:{S.at['relaxed pr < 0.85, partner < 150 um', 'p_within_mouse']:.2g}")

ax = axs[1]
for so in (KS, IB):
    d = sweep[sweep.sorter == so]
    ax.plot(d.threshold, d.share_best_near, color=COL[so], lw=2, marker="o", ms=4, label=NAME[so])
    ax.axhline(chance_near[so], color=COL[so], lw=1, ls=(0, (1, 2)))
ax.text(0.505, chance_near.max() + 0.015, "chance: all CS × SS pairs", fontsize=7.5, color=INK2)
ax.set_xlabel("pause-ratio threshold")
ax.set_ylabel(f"share of best partners < {NEAR_UM:.0f} µm from the CS")
ax.set_ylim(0, 1.05)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
ax.legend(frameon=False, fontsize=8, loc="lower left")
clean(ax)
title(ax, "B  The extra KS4 pairs are mostly far away",
      "a real CS–SS pair is near; as the threshold loosens, KS4\npartners drift towards chance (likely spurious)")

ax = axs[2]
rng = np.random.default_rng(4)
bn_all["ksc"] = np.where(bn_all.candidate, "paired", "unpaired")
for so in (IB, KS):     # KS4 drawn on top
    d = bn_all[bn_all.sorter == so]
    ax.scatter(d.c_est + rng.uniform(0, 0.006, len(d)), d.pause_ratio, s=np.where(d.candidate, 28, 14),
               color=COL[so], alpha=np.where(d.candidate, 0.95, 0.5), edgecolor="white", lw=0.5, label=NAME[so])
cc = np.linspace(0, C_MAX, 50)
ax.plot(cc, cc + (1 - cc) * 0.3, color=INK, lw=1.2, ls=(0, (4, 2)))
ax.text(0.62, 0.55, "a true pause of 0.3\ndiluted by contamination c", fontsize=7.5, color=INK, va="top")
ax.axhline(PR_STEP0, color=INK2, lw=0.8, ls=(0, (2, 2)))
ax.set_xlabel("estimated contamination c of that unit (from ISI < 1 ms)")
ax.set_ylabel("pause ratio of the best unit < 150 µm")
ax.set_xlim(-0.02, 0.95)
ax.set_ylim(0, 1.3)
ax.legend(frameon=False, fontsize=8, loc="lower right")
clean(ax)
med_c = bn_all.groupby("sorter").c_est.median()
title(ax, "C  Dilution by multi-unit spikes explains little",
      f"best unit < 150 µm per CS (large = paired); median c: KS4 {med_c[KS]:.2f}, IBL {med_c[IB]:.2f}.\n"
      f"KS4 CS without a nearby pause whose unit is contaminated\nenough to hide a 0.3 pause: {n_hide} of {n_ks_np}.\n"
      f"Corrected: KS4 {S.at['contamination-corrected pr < 0.7', 'ks4_n']} (was {S.at['step0 (pr < 0.7)', 'ks4_n']}), "
      f"IBL {S.at['contamination-corrected pr < 0.7', 'ibl_n']} (was {S.at['step0 (pr < 0.7)', 'ibl_n']})")

for ext, kw in (("png", {"dpi": 200}), ("svg", {})):
    f = OUT / f"fig7_weak_pause_and_contamination.{ext}"
    fig.savefig(f, **kw)
    try:
        shutil.copy2(f, FINAL_FIG / f.name)
    except OSError as e:
        print(f"WARNING: could not copy {f.name} to FINAL ({e})")
plt.close(fig)
print("tier-2 ISI violations KS4 vs IBL: MW p =", f"{p_c_t2:.2g}")


# %% 4. CHANCE-CORRECTED PAIRING (no distance assumption needed below FAR_UM)
# Question (Simon): IBL handles drift better; KS4 may split a cell's SS so its partner sits farther away. Is IBL's
# strict 31% comparable with KS4's relaxed rate once longer distances are allowed?
# Chance level: pairs farther than FAR_UM cannot be the same Purkinje cell (CS-SS of one cell, even drift-split,
# are well within that). Their pass rate f (per session; per sorter if a session has no far pairs) is the per-pair
# false-positive rate at that threshold. For a CS with n candidates within the distance limit D, the chance of
# >= 1 spurious pass is 1 - (1 - f)^n. Excess = observed paired CS - expected spurious, per CS.
FAR_UM = 500.0
D_LIMITS = [150.0, 300.0, 500.0, np.inf]
D_NAME = {150.0: "< 150 µm", 300.0: "< 300 µm", 500.0: "< 500 µm", np.inf: "any distance"}
far_all = P[P.dist_um > FAR_UM]
ses_n = ALL_CS.groupby(["sorter", "mouse", "S"]).size().rename("n_CS")


def chance_corrected(thr, D):
    ok = (P.pause_ratio < thr) & P.ok_other
    f_ses = ok[far_all.index].groupby([far_all.sorter, far_all.mouse, far_all.S]).mean()
    f_sor = ok[far_all.index].groupby(far_all.sorter).mean()
    sub = P[P.dist_um < D]
    n = sub.groupby(KEY).size().rename("n").reset_index()
    n["f"] = [f_ses.get((a, b, c), f_sor[a]) for a, b, c in n[["sorter", "mouse", "S"]].to_numpy()]
    n["p_chance"] = 1 - (1 - n.f) ** n.n
    n["obs"] = n.set_index(KEY).index.isin(sub[ok[sub.index]].set_index(KEY).index)
    s = n.groupby(["sorter", "mouse", "S"]).agg(obs=("obs", "sum"), exp=("p_chance", "sum"))
    s = s.reindex(ses_n.index, fill_value=0).join(ses_n)
    s["excess"] = (s.obs - s.exp) / s.n_CS
    return s.reset_index(), f_sor


cc_rows, cc_ses = [], {}
for thr in THRESHOLDS:
    for D in D_LIMITS:
        s, f_sor = chance_corrected(thr, D)
        cc_ses[(thr, D)] = s
        for so in (KS, IB):
            q = s[s.sorter == so]
            cc_rows.append(dict(threshold=thr, dist_limit_um=D, sorter=so, n_CS=q.n_CS.sum(), obs=q.obs.sum(),
                                exp_chance=q.exp.sum(), frac_obs=q.obs.sum() / q.n_CS.sum(),
                                frac_excess=(q.obs.sum() - q.exp.sum()) / q.n_CS.sum(), f_far_pair=f_sor[so],
                                session_median_excess=q.excess.median()))
cc = pd.DataFrame(cc_rows)
cc.to_csv(OUT / "chance_corrected_pairing.csv", index=False)
print(cc[cc.threshold.isin([0.7, 0.85, 0.95])].round(3).to_string(index=False))


def cc_get(so, thr, D, col="frac_excess"):
    return cc[(cc.sorter == so) & np.isclose(cc.threshold, thr) & (cc.dist_limit_um == D)][col].iloc[0]


# session-level tests: same rule for both sorters, and the "fair" cross comparison asked for
# (IBL strict: < 0.7, any distance  vs  KS4 most permissive: < 0.95, < 500 um)
# The pooled numbers are dominated by one very large session per sorter (KS4 Ferrero S6: 66 of 123 CS; IBL
# Toblerone S5: 78 of 228), both with low pairing; "without largest session" drops that session in each sorter.
cmp_rows, cmp_ses = [], []
for name, (rk, ri, drop_big) in {
        "both < 0.7, any distance": ((0.7, np.inf), (0.7, np.inf), False),
        "both < 0.95, < 500 um": ((0.95, 500.0), (0.95, 500.0), False),
        "both < 0.95, any distance": ((0.95, np.inf), (0.95, np.inf), False),
        "KS4 < 0.95 & < 500 um vs IBL < 0.7 any distance": ((0.95, 500.0), (0.7, np.inf), False),
        "both < 0.7, any distance, without largest session": ((0.7, np.inf), (0.7, np.inf), True),
        "both < 0.95, < 500 um, without largest session": ((0.95, 500.0), (0.95, 500.0), True)}.items():
    sk, si = cc_ses[rk], cc_ses[ri]
    d = pd.concat([sk[sk.sorter == KS], si[si.sorter == IB]])
    if drop_big:
        d = d[d.n_CS < d.groupby("sorter").n_CS.transform("max")]
    diff, p_wm = within_mouse_perm(d, "excess")
    cmp_ses.append(d.assign(comparison=name))
    pk, pi = d[d.sorter == KS], d[d.sorter == IB]
    cmp_rows.append(dict(comparison=name, n_CS_ks4=pk.n_CS.sum(), n_CS_ibl=pi.n_CS.sum(),
                         ks4_pooled_excess=(pk.obs.sum() - pk.exp.sum()) / pk.n_CS.sum(),
                         ibl_pooled_excess=(pi.obs.sum() - pi.exp.sum()) / pi.n_CS.sum(),
                         ks4_pooled_raw=pk.obs.sum() / pk.n_CS.sum(), ibl_pooled_raw=pi.obs.sum() / pi.n_CS.sum(),
                         ks4_session_median=d[d.sorter == KS].excess.median(),
                         ibl_session_median=d[d.sorter == IB].excess.median(),
                         p_mannwhitney=mannwhitneyu(d[d.sorter == KS].excess, d[d.sorter == IB].excess).pvalue,
                         p_within_mouse=p_wm))
cmpc = pd.DataFrame(cmp_rows)
cmpc.to_csv(OUT / "chance_corrected_comparisons.csv", index=False)
cmp_ses = pd.concat(cmp_ses, ignore_index=True)
cmp_ses.to_csv(OUT / "chance_corrected_comparisons_sessions.csv", index=False)
print(cmpc.round(3).to_string(index=False))

# who are the partners gained by loosening (< 0.95, < 500 um) relative to step0 (< 0.7, any distance)?
strict_cs = set(map(tuple, P[P.candidate][KEY].drop_duplicates().to_numpy()))
gain = P[(P.pause_ratio < 0.95) & P.ok_other & (P.dist_um < 500)].sort_values("pause_ratio").groupby(KEY).head(1).copy()
gain["new"] = [tuple(r) not in strict_cs for r in gain[KEY].to_numpy()]
gain.to_csv(OUT / "partners_relaxed_095_500um.csv", index=False)
print(gain.groupby(["sorter", "new"]).agg(n=("cs", "size"), dist_med=("dist_um", "median"),
                                          pr_med=("pause_ratio", "median"), dur50_med=("pause_dur50_ms", "median"),
                                          fr_med=("ss_fr", "median"),
                                          frac_ss_like=("ss_fr", lambda x: (x >= 20).mean())).round(2).to_string())

# %% 5. FIGURE 8
fig, axs = plt.subplots(1, 3, figsize=(16, 5.6), gridspec_kw={"width_ratios": [1.05, 1.0, 1.15]})
fig.subplots_adjust(left=0.05, right=0.98, top=0.7, bottom=0.15, wspace=0.3)
fig.suptitle("Beyond chance, loosening threshold and distance lifts KS4 from 7% to 18%; IBL stays ~10–15 points ahead "
             "(also without the largest sessions)",
             x=0.05, ha="left", fontsize=12, fontweight="bold", color=INK, y=0.97)
fig.text(0.05, 0.925, f"Chance = pass rate of CS × SS pairs > {FAR_UM:.0f} µm apart (cannot be one Purkinje cell), per session, "
         "applied to every CS's candidates within the distance limit. Excess = observed − chance.\n"
         "With loose rules the session-level difference is no longer significant (C).", fontsize=8.5, color=INK2, va="top")

ax = axs[0]
for so in (KS, IB):
    for D, ls, lw in ((500.0, "-", 2.2), (np.inf, (0, (1, 1.5)), 1.6), (150.0, (0, (4, 2)), 1.6)):
        d = cc[(cc.sorter == so) & (cc.dist_limit_um == D)]
        ax.plot(d.threshold, d.frac_excess, color=COL[so], ls=ls, lw=lw, label=f"{NAME[so]}, {D_NAME[D]}")
    d = cc[(cc.sorter == so) & (cc.dist_limit_um == np.inf)]
    ax.plot(d.threshold, d.frac_obs, color=COL[so], lw=0.8, alpha=0.45)
ax.text(0.952, cc_get(KS, 0.95, np.inf, "frac_obs"), "raw", fontsize=7, color=COL[KS], va="center")
ax.text(0.952, cc_get(IB, 0.95, np.inf, "frac_obs"), "raw", fontsize=7, color=COL[IB], va="center")
ax.axvline(PR_STEP0, color=INK2, lw=0.8, ls=(0, (2, 2)))
ax.text(PR_STEP0 + 0.005, 0.012, "step0 threshold", fontsize=7.5, color=INK2)
ax.set_xlabel("pause-ratio threshold")
ax.set_ylabel("fraction of CS paired beyond chance (pooled)")
ax.set_ylim(0, 0.45)
ax.set_xlim(0.48, 0.99)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
ax.legend(frameon=False, fontsize=7, loc="upper left", ncol=2)
clean(ax)
title(ax, "A  Pairing beyond chance, by threshold and distance limit",
      f"KS4 best: {cc[cc.sorter == KS].frac_excess.max():.0%}; IBL range {cc[cc.sorter == IB].frac_excess.min():.0%}–"
      f"{cc[cc.sorter == IB].frac_excess.max():.0%}. Thin lines = raw, any distance.")

ax = axs[1]
for so, mk in ((IB, "o"), (KS, "s")):
    for new, a, s_ in ((False, 0.35, 22), (True, 0.95, 46)):
        d = gain[(gain.sorter == so) & (gain.new == new)]
        ax.scatter(d.dist_um, d.ss_fr, marker=mk, s=s_, color=COL[so] if new else "white", alpha=a if new else 1,
                   edgecolor=COL[so], lw=1, label=f"{NAME[so]}: {'gained' if new else 'step0 pair'} (n = {len(d)})",
                   zorder=3 if new else 2)
ax.axhline(20, color=INK2, lw=0.8, ls=(0, (2, 2)))
ax.text(495, 19, "20 Hz (tier-2 floor; SS fire tens of Hz)", ha="right", fontsize=7.5, color=INK2, va="top")
ax.set_yscale("log")
ax.set_ylim(0.8, 600)
ax.set_xlabel("CS–partner distance (µm)")
ax.set_ylabel("partner firing rate (Hz, log)")
ax.legend(frameon=False, fontsize=7.3, loc="upper right", ncol=2)
clean(ax)
gk = gain[(gain.sorter == KS) & gain.new]
title(ax, "B  The partners gained by loosening are not SS-like",
      f"< 0.95 and < 500 µm vs step0. KS4 gained {len(gk)}: {int((gk.ss_fr < 20).sum())} fire < 20 Hz, "
      f"{int((gk.mouse == '5_FerreroRocher').sum())} from Ferrero S6;\n"
      f"pause 50% duration median {gk.pause_dur50_ms.median():.0f} ms "
      f"(step0 pairs {gain[~gain.new].pause_dur50_ms.median():.0f} ms)")

ax = axs[2]
row = cmpc.set_index("comparison")
items = [("both < 0.7, any distance", "both step0\n(< 0.7)"),
          ("both < 0.95, < 500 um", "both loosest\n(< 0.95,\n< 500 µm)"),
          ("KS4 < 0.95 & < 500 um vs IBL < 0.7 any distance", "KS4 loosest\nvs IBL\nstep0"),
          ("both < 0.7, any distance, without largest session", "step0,\nwithout\nlargest\nsession"),
          ("both < 0.95, < 500 um, without largest session", "loosest,\nwithout\nlargest\nsession")]
for j, (key, lab) in enumerate(items):
    for k, so in enumerate((KS, IB)):
        raw = row.at[key, f"{so}_pooled_raw"]
        exc = row.at[key, f"{so}_pooled_excess"]
        x = j + (k - 0.5) * 0.36
        ax.bar(x, raw, width=0.34, color=COL[so], alpha=0.25, edgecolor="white", lw=1)
        ax.bar(x, exc, width=0.34, color=COL[so], edgecolor="white", lw=1, label=NAME[so] if j == 0 else None)
        ax.text(x, raw + 0.008, f"{exc:.0%}", ha="center", fontsize=7.8, color=INK)
    top = max(row.at[key, "ks4_pooled_raw"], row.at[key, "ibl_pooled_raw"])
    ax.text(j, top + 0.045, f"p = {row.at[key, 'p_within_mouse']:.2g}", ha="center", fontsize=7.5, color=INK2)
ax.axvline(2.5, color=GRID_C, lw=1)
ax.set_xticks(range(len(items)), [i[1] for i in items], fontsize=7.3)
ax.set_ylim(0, 0.6)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
ax.set_ylabel("fraction of CS paired (pooled)")
ax.legend(frameon=False, fontsize=8, loc="upper left")
clean(ax)
nk = int(row.at["both < 0.7, any distance, without largest session", "n_CS_ks4"])
ni = int(row.at["both < 0.7, any distance, without largest session", "n_CS_ibl"])
title(ax, "C  Fair comparisons",
      "solid = beyond chance (number), pale = raw;\np = within-mouse on session excess. Largest session:\n"
      f"KS4 Ferrero S6 (66 CS), IBL Toblerone S5 (78 CS); without: {nk} vs {ni} CS")

for ext, kw in (("png", {"dpi": 200}), ("svg", {})):
    f = OUT / f"fig8_chance_corrected_distance.{ext}"
    fig.savefig(f, **kw)
    try:
        shutil.copy2(f, FINAL_FIG / f.name)
    except OSError as e:
        print(f"WARNING: could not copy {f.name} to FINAL ({e})")
plt.close(fig)
print("outputs in", OUT)
