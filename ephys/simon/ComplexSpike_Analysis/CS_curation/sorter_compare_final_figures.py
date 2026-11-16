# -*- coding: utf-8 -*-
"""
sorter_compare_final_figures.py — final KS4 vs IBL comparison: one central folder with the key figures and tables.

Question: does one sorter (Kilosort4 vs IBL sorter, each followed by C4) do a better job at (1) detecting complex
spikes and (2) giving a CS unit whose simple-spike partner can be found?

Inputs (already computed, nothing recomputed from raw data):
  _sorter_comparison/session_summary.csv, units_all.csv, stats_session_level.csv   (pass 1: sorter_compare_CS_SS.py)
  _sorter_comparison/<STEP0_SUMMARY>/cs_units.csv, session_summary.csv            (step0 v5 batch aggregate:
                                                                                    sorter_compare_step0v3.py)
  _sorter_comparison/<ZERO_LAG>/cs_zero_lag.csv                                    (sorter_compare_zero_lag_check.py)
Population: extra-good CB units, C4 PkC_cs = CS. SS partners searched in extra-good CB units (tier 1) and any other
  CB units >= 20 Hz (tier 2). Window: win_selected, extended per unit to its stable span (step0 v3).
Statistics: session = unit of replication (each session was sorted by ONE sorter). Within-mouse permutation test
  (sorter labels shuffled within mouse; statistic = mean over mice of median KS4 - median IBL; 10 000 perms), and
  Mann-Whitney on sessions. Uncorrected, exploratory: 5 mice, 11 KS4 / 16 IBL sessions.
CAVEAT: sorter is confounded with session (no session was sorted by both), so a sorter difference could also be a
  recording difference.

Output (local only): CS_curation/_sorter_comparison/FINAL_KS4_vs_IBL/  figs/*.png|svg, tables/*.csv, README.md
GOLDEN RULE: Dropbox is READ-ONLY.
"""
# %% 1. PARAMETERS + INPUTS
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu, spearmanr

DROPBOX_ROOT = Path(r"D:\Learning Lab Dropbox")
CODE_DIR = Path(r"C:\Users\simon\Desktop\Code\EPhys\ComplexSpike_Analysis")
CMP = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation\_sorter_comparison")
STEP0_SUMMARY = "step0_v5batch"      # aggregate of step0_explore_bhvspan_v5batch (tier 2 without label gate)
STEP0_PREV = "step0v3"               # previous reference (tier 2 = cluster_group good/mua) for the v3 → v5 comparison
ZERO_LAG = "zero_lag_check_v5batch"
N_PERM = 10000
MIN_CS_SCATTER = 3                   # sessions with fewer CS give 0/1-type pairing fractions: drawn hollow, no rho

KS, IB = "ks4", "ibl"
COL = {KS: "#2a78d6", IB: "#eb6834"}  # dataviz reference palette, categorical slots 1-2
NAME = {KS: "Kilosort4", IB: "IBL sorter"}
INK, INK2, GRID_C, NEUTRAL = "#0b0b0b", "#52514e", "#e4e3df", "#d9d8d3"
MICE = ["1_CoteDor", "2_Lindt", "3_Toblerone", "4_Milka", "5_FerreroRocher"]
MARK = dict(zip(MICE, ["o", "s", "^", "D", "v"]))
MNAME = {m: m.split("_", 1)[1].replace("FerreroRocher", "Ferrero") for m in MICE}


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


OUT = safe_out_dir(CMP / "FINAL_KS4_vs_IBL")
FIG = safe_out_dir(OUT / "figs")
TAB = safe_out_dir(OUT / "tables")

s1 = pd.read_csv(CMP / "session_summary.csv")
u1 = pd.read_csv(CMP / "units_all.csv")
s3 = pd.read_csv(CMP / STEP0_SUMMARY / "session_summary.csv")
c3 = pd.read_csv(CMP / STEP0_SUMMARY / "cs_units.csv")
zl = pd.read_csv(CMP / ZERO_LAG / "cs_zero_lag.csv")

ses = s3.merge(s1[["mouse", "S", "n_eg_CB", "n_SS", "frac_SS_of_eg", "frac_CS_of_eg", "CS_fr_med",
                   "CS_footprint_span_um_med", "SS_footprint_span_um_med", "CS_amp_rcv_med", "SS_fr_med",
                   "CS_cv2_med", "win_dur_min"]], on=["mouse", "S"], how="left")
ses["n_SS_pool"] = ses["n_SS_t1"] + ses["n_SS_t2"]
# SS units that are the best partner of >= 1 paired CS (distinct units per session); 0 where the session has CS
_pss = c3[c3.paired_any].groupby(["mouse", "S"])["best_ss"].nunique().rename("n_SS_paired")
ses = ses.merge(_pss, on=["mouse", "S"], how="left")
ses.loc[ses["n_CS"] > 0, "n_SS_paired"] = ses.loc[ses["n_CS"] > 0, "n_SS_paired"].fillna(0)
zs = zl.assign(shared=zl["best_excess"] > 0.2).groupby(["mouse", "S"])["shared"].mean().rename("frac_CS_shared")
ses = ses.merge(zs, on=["mouse", "S"], how="left")
ses["sorter_name"] = ses["sorter"].map(NAME)

# Kilosort's automatic label (cluster_KSLabel.tsv: good/mua from contamination) is written by BOTH sorters with the
# same rule, unlike cluster_group (Phy-curated in KS4, unclear origin in IBL). It is used here ONLY to describe the
# tier-2 pool and the best partners (never to select). Dropbox files are read, never written.
STEP0_RUN = "step0_explore_bhvspan_" + STEP0_SUMMARY.replace("step0_", "")   # -> step0_explore_bhvspan_v5batch
stab = pd.read_csv(CMP / "session_table.csv")
_ksl, _pool = {}, []
for _, r in stab.iterrows():
    f = Path(r.sorter_dir) / "cluster_KSLabel.tsv"
    lab = pd.read_csv(f, sep="\t").set_index("cluster_id")["KSLabel"] if f.exists() else pd.Series(dtype=str)
    _ksl[(r.mouse, r.S)] = lab
    ut = CMP.parent / f"{r.mouse}_{r.R}" / STEP0_RUN / "tables" / "unit_table.csv"
    if ut.exists():
        u = pd.read_csv(ut)
        u = u[u.role == "SS_cand2"]
        _pool.append(pd.DataFrame(dict(mouse=r.mouse, S=r.S, sorter=r.sorter, phyID=u.phyID.to_numpy(),
                                       group=u.group.fillna("unsorted").to_numpy(),
                                       kslabel=u.phyID.map(lab).fillna("none").to_numpy())))
pool2 = pd.concat(_pool, ignore_index=True)   # every tier-2 SS candidate, one row per unit
c3["best_kslabel"] = [(_ksl[(m, s)].get(int(b), "none") if pd.notna(b) else np.nan)
                      for m, s, b in zip(c3.mouse, c3.S, c3.best_ss)]
_t2 = pool2.assign(ksgood=pool2.kslabel == "good").groupby(["mouse", "S"])["ksgood"].mean().rename("frac_t2_ksgood")
ses = ses.merge(_t2, on=["mouse", "S"], how="left")


# %% 2. STATISTICS
def within_mouse_perm(df, col, n_perm=N_PERM, seed=0):
    d = df.dropna(subset=[col])
    d = d[d.groupby("mouse")["sorter"].transform("nunique") == 2]
    if d.empty:
        return np.nan, np.nan, 0
    groups = [(g[col].to_numpy(float), (g["sorter"] == KS).to_numpy()) for _, g in d.groupby("mouse")]

    def stat(labs):
        return np.mean([np.median(v[l]) - np.median(v[~l]) for (v, _), l in zip(groups, labs)])

    obs = stat([l for _, l in groups])
    rng = np.random.default_rng(seed)
    null = np.array([stat([rng.permutation(l) for _, l in groups]) for _ in range(n_perm)])
    return obs, (1 + np.sum(np.abs(null) >= abs(obs) - 1e-12)) / (1 + n_perm), len(groups)


METRICS = [  # column, label, family
    ("n_CS", "extra-good CS units per session", "detection"),
    ("frac_CS_of_eg", "CS as fraction of extra-good CB units", "detection"),
    ("n_SS", "C4 PkC_ss (extra-good) per session", "detection"),
    ("frac_SS_of_eg", "PkC_ss as fraction of extra-good CB units", "detection"),
    ("n_SS_pool", "SS-partner candidates per session (tier 1 + 2)", "SS pool"),
    ("n_SS_paired", "SS units paired with a CS per session", "pairing"),
    ("frac_paired_any", "fraction of CS with an SS partner (tier 1+2)", "pairing"),
    ("frac_paired_t1", "fraction of CS with an extra-good SS partner", "pairing"),
    ("frac_strong", "fraction of CS with a strong pause (<0.3, >=5 ms)", "pairing"),
    ("frac_paired_stable", "fraction of CS paired in every 600 s block", "pairing"),
    ("best_pr_med", "pause ratio of found pairs (median)", "pairing"),
    ("n_SS_t1", "SS-partner pool: extra-good units", "SS pool"),
    ("n_SS_t2", "SS-partner pool: other CB units >= 20 Hz (any label)", "SS pool"),
    ("frac_t2_ksgood", "tier-2 pool: fraction with automatic KSLabel good", "SS pool"),
    ("frac_CS_shared", "fraction of CS with >20% spikes in a neighbour", "duplicates"),
    ("frac_CS_dup_collapsed", "fraction of CS removed as same-cell duplicates", "duplicates"),
    ("CS_fr_med", "CS firing rate (Hz, median)", "unit properties"),
    ("CS_cv2_med", "CS CV2 (median)", "unit properties"),
    ("CS_footprint_span_um_med", "CS template footprint span (um)", "unit properties"),
    ("SS_footprint_span_um_med", "SS template footprint span (um)", "unit properties"),
    ("CS_amp_rcv_med", "CS amplitude spread (IQR/median)", "unit properties"),
]
rows = []
for col, lab, fam in METRICS:
    k, i = ses.loc[ses.sorter == KS, col].dropna(), ses.loc[ses.sorter == IB, col].dropna()
    diff, p_wm, n_m = within_mouse_perm(ses, col)
    rows.append(dict(family=fam, metric=col, label=lab, ks4_median=k.median(), ibl_median=i.median(),
                     n_ks4=len(k), n_ibl=len(i), p_mannwhitney=mannwhitneyu(k, i).pvalue if len(k) and len(i) else np.nan,
                     within_mouse_diff_ks4_minus_ibl=diff, p_within_mouse=p_wm, n_mice=n_m))
stats = pd.DataFrame(rows)
stats.to_csv(TAB / "final_stats_session_level.csv", index=False)
ST = stats.set_index("metric")
print(stats[["metric", "ks4_median", "ibl_median", "n_ks4", "n_ibl", "p_mannwhitney", "p_within_mouse"]]
      .round(3).to_string(index=False))


def ptxt(col):
    return f"within-mouse p = {ST.at[col, 'p_within_mouse']:.2g}"


# %% 3. PLOT HELPERS
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.titleweight": "regular", "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.8,
                     "svg.fonttype": "none", "figure.facecolor": "white"})


def clean(ax, grid_y=True):
    ax.spines[["top", "right"]].set_visible(False)
    if grid_y:
        ax.yaxis.grid(True, color=GRID_C, lw=0.6)
        ax.set_axisbelow(True)


def title(ax, head, sub):
    ax.set_title(head + "\n" * (1 + sub.count("\n")), loc="left", fontsize=10, fontweight="bold", color=INK)
    ax.text(0, 1.02, sub, transform=ax.transAxes, fontsize=8.3, color=INK2, va="bottom")


def strip(ax, col, ylabel, fmt=None):
    rng = np.random.default_rng(1)
    for k, so in enumerate((KS, IB)):
        d = ses[(ses.sorter == so)].dropna(subset=[col])
        for m in MICE:
            dm = d[d.mouse == m]
            ax.scatter(k + rng.uniform(-0.17, 0.17, len(dm)), dm[col], marker=MARK[m], s=34, color=COL[so],
                       edgecolor="white", lw=0.8, zorder=3)
        med = d[col].median()
        ax.plot([k - 0.28, k + 0.28], [med, med], color=INK, lw=2, solid_capstyle="round", zorder=4)
        ax.text(k + 0.31, med, fmt(med) if fmt else f"{med:.2g}", va="center", fontsize=8, color=INK)
    ax.set_xticks([0, 1], [f"{NAME[KS]}\n{(ses.sorter == KS).sum()} sessions", f"{NAME[IB]}\n{(ses.sorter == IB).sum()} sessions"])
    ax.set_xlim(-0.6, 1.75)
    ax.set_ylabel(ylabel)
    clean(ax)


def mouse_legend(fig, y=0.005):
    h = [Line2D([], [], marker=MARK[m], ls="", color=INK2, markersize=6, label=MNAME[m]) for m in MICE]
    fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, fontsize=8, bbox_to_anchor=(0.5, y),
               title="mouse (marker)", title_fontsize=8)


def save(fig, name):
    for ext, kw in (("png", {"dpi": 200}), ("svg", {})):
        try:
            fig.savefig(FIG / f"{name}.{ext}", **kw)
            print("saved", FIG / f"{name}.{ext}")
        except OSError as e:     # e.g. the file is open in a viewer (Windows lock): keep going, report it
            print(f"WARNING: could not write {name}.{ext} ({e}); close it in any viewer and rerun")
    plt.close(fig)


# %% 4. FIGURE 1 — HEADLINE: detection equal, pairing better with IBL, pairs equally good
fig, axs = plt.subplots(2, 2, figsize=(11, 9.2))
fig.subplots_adjust(left=0.07, right=0.97, top=0.86, bottom=0.12, hspace=0.62, wspace=0.32)
fig.suptitle("Kilosort4 vs IBL sorter: same complex spikes found, but IBL gives their simple spikes more often",
             x=0.07, ha="left", fontsize=12.5, fontweight="bold", color=INK, y=0.975)
fig.text(0.07, 0.925, f"Extra-good CB units, C4 PkC_cs = CS; {len(MICE)} mice, one sorter per session "
         f"(sorter is confounded with session). One point = one session.", fontsize=8.5, color=INK2)

# A: detection
ax = axs[0, 0]
strip(ax, "n_CS", "extra-good CS units per session", lambda v: f"{v:g}")
nk, ni = c3[c3.sorter == KS].shape[0], c3[c3.sorter == IB].shape[0]
title(ax, "A  CS detection: no difference",
      f"median {ST.at['n_CS', 'ks4_median']:.0f} vs {ST.at['n_CS', 'ibl_median']:.1f} CS per session; "
      f"{ptxt('n_CS')}  (total {nk} vs {ni})")

# B: pairing slope by mouse
ax = axs[0, 1]
d = ses.dropna(subset=["frac_paired_any"])
rng = np.random.default_rng(2)
for k, so in enumerate((KS, IB)):
    ds = d[d.sorter == so]
    ax.scatter(k + rng.uniform(-0.1, 0.1, len(ds)), ds.frac_paired_any, s=14, color=COL[so], alpha=0.35,
               edgecolor="none", zorder=2)
n_up = n_eq = 0
lab_y = []
for m in MICE:
    dm = d[d.mouse == m]
    if dm.sorter.nunique() < 2:
        continue
    yk, yi = dm.loc[dm.sorter == KS, "frac_paired_any"].median(), dm.loc[dm.sorter == IB, "frac_paired_any"].median()
    n_up += yi > yk + 1e-9
    n_eq += abs(yi - yk) <= 1e-9
    ck, ci = int(dm.loc[dm.sorter == KS, "n_CS"].sum()), int(dm.loc[dm.sorter == IB, "n_CS"].sum())
    ax.plot([0, 1], [yk, yi], color=INK2, lw=1.2, zorder=3)
    ax.scatter([0], [yk], marker=MARK[m], s=46, color=COL[KS], edgecolor="white", lw=0.8, zorder=4)
    ax.scatter([1], [yi], marker=MARK[m], s=46, color=COL[IB], edgecolor="white", lw=0.8, zorder=4)
    lab_y.append([yi, f"{MNAME[m]}  ({ck} → {ci} CS)"])
lab_y.sort()
for j in range(1, len(lab_y)):          # keep mouse labels >= 0.055 apart
    lab_y[j][0] = max(lab_y[j][0], lab_y[j - 1][0] + 0.055)
for y, t in lab_y:
    ax.text(1.08, y, t, va="center", fontsize=7.5, color=INK2)
ax.set_xticks([0, 1], [NAME[KS], NAME[IB]])
ax.set_xlim(-0.35, 1.85)
ax.set_ylim(-0.03, 1.08)
ax.set_ylabel("fraction of CS with an SS partner")
clean(ax)
title(ax, f"B  SS partner found: IBL higher in {n_up} of 5 mice ({n_eq} equal)",
      f"session medians {ST.at['frac_paired_any', 'ks4_median']:.2f} vs {ST.at['frac_paired_any', 'ibl_median']:.2f}; "
      f"{ptxt('frac_paired_any')}\n"
      "large = mouse median of sessions (labels: CS units KS4 → IBL), faint = sessions")

# C: pooled route of pairing
ax = axs[1, 0]
cats = [("tier 1: extra-good SS", lambda c: c.paired_any & (c.best_tier == 1), 1.0),
        ("tier 2, KSLabel good", lambda c: c.paired_any & (c.best_tier == 2) & (c.best_kslabel == "good"), 0.6),
        ("tier 2, KSLabel mua", lambda c: c.paired_any & (c.best_tier == 2) & (c.best_kslabel == "mua"), 0.32),
        ("no partner", lambda c: ~c.paired_any, None)]
route = []
for yk, so in enumerate((IB, KS)):
    c = c3[c3.sorter == so]
    left = 0.0
    for name, f, a in cats:
        n = int(f(c).sum())
        w = n / len(c)
        route.append(dict(sorter=so, route=name, n=n, frac=w))
        ax.barh(yk, w, left=left, height=0.62, color=COL[so] if a else NEUTRAL, alpha=a if a else 1,
                edgecolor="white", lw=2)
        if w > 0.018:
            ax.text(left + w / 2, yk, f"{n}", ha="center", va="center", fontsize=8.5 if w > 0.045 else 7,
                    color="white" if a and a > 0.5 else INK)
        left += w
    ax.text(1.01, yk, f"{c.paired_any.mean():.0%} paired\n(n = {len(c)} CS)", va="center", fontsize=8, color=INK)
ax.set_yticks([0, 1], [NAME[IB], NAME[KS]])
ax.set_xlim(0, 1.18)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1], ["0", "25%", "50%", "75%", "100%"])
ax.set_xlabel("share of all CS units (pooled over sessions); numbers = CS units")
ax.spines[["top", "right"]].set_visible(False)
x0, y0 = 0.0, -0.69
for name, _, a in cats:   # legend: a KS4 + IBL swatch pair per category (same shading as the bars)
    w = 0.075 + 0.0125 * len(name) + 0.05
    if x0 > 0 and x0 + w > 1.15:
        x0, y0 = 0.0, y0 - 0.3
    for j, so in enumerate((KS, IB) if a else (KS,)):
        ax.add_patch(plt.Rectangle((x0 + 0.032 * j, y0 - 0.09), 0.028 if a else 0.06, 0.18, color=COL[so] if a else NEUTRAL,
                                   alpha=a if a else 1, clip_on=False, transform=ax.transData))
    ax.text(x0 + 0.075, y0, name, fontsize=7.5, va="center", color=INK2)
    x0 += w
ax.set_ylim(y0 - 0.2, 1.45)
title(ax, "C  The extra IBL pairs come from tier-2 units",
      "best SS partner of each CS; tier 2 = other CB units ≥ 20 Hz, any label, split by\n"
      "Kilosort's automatic KSLabel (same rule in both sorters; the curated label is not used)")
pd.DataFrame(route).to_csv(TAB / "pairing_route_pooled.csv", index=False)

# D: quality of the pairs that are found
ax = axs[1, 1]
for so in (KS, IB):
    v = np.sort(c3.loc[(c3.sorter == so) & c3.paired_any, "best_pr"].dropna().to_numpy())
    ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=COL[so], lw=2,
            label=f"{NAME[so]} (n = {len(v)} pairs, median {np.median(v):.2f})")
ax.axvline(0.7, color=INK2, lw=0.8, ls=(0, (2, 2)))
ax.text(0.69, 0.04, "candidate\nthreshold 0.7", ha="right", fontsize=7.5, color=INK2)
ax.set_xlim(0, 0.75)
ax.set_ylim(0, 1.02)
ax.set_xlabel("SS pause ratio after the CS (0 = full silence, 1 = no pause)")
ax.set_ylabel("cumulative fraction of paired CS")
ax.legend(frameon=False, fontsize=8, loc="upper left")
clean(ax)
title(ax, "D  When a partner is found, its pause is similar",
      f"session medians {ST.at['best_pr_med', 'ks4_median']:.2f} vs {ST.at['best_pr_med', 'ibl_median']:.2f}; "
      f"{ptxt('best_pr_med')}")
mouse_legend(fig)
save(fig, "fig1_headline_detection_vs_pairing")

# %% 5. FIGURE 2 — WHY: IBL yields more SS-like good/mua units, and pairing tracks that pool
fig, axs = plt.subplots(1, 2, figsize=(11, 4.9), gridspec_kw={"width_ratios": [1, 1.25]})
fig.subplots_adjust(left=0.07, right=0.97, top=0.76, bottom=0.2, wspace=0.28)
fig.suptitle(f"Why? IBL tends to offer more tier-2 SS candidates (median {ST.at['n_SS_t2', 'ibl_median']:.0f} vs "
             f"{ST.at['n_SS_t2', 'ks4_median']:.0f}), and pairing rises with that pool",
             x=0.07, ha="left", fontsize=12, fontweight="bold", color=INK, y=0.97)

ax = axs[0]
rng = np.random.default_rng(3)
pos = {("n_SS_t1", KS): 0, ("n_SS_t1", IB): 1, ("n_SS_t2", KS): 2.6, ("n_SS_t2", IB): 3.6}
for (col, so), x in pos.items():
    dd = ses[ses.sorter == so].dropna(subset=[col])
    for m in MICE:
        dm = dd[dd.mouse == m]
        ax.scatter(x + rng.uniform(-0.17, 0.17, len(dm)), dm[col], marker=MARK[m], s=30, color=COL[so],
                   edgecolor="white", lw=0.7, zorder=3)
    med = dd[col].median()
    ax.plot([x - 0.28, x + 0.28], [med, med], color=INK, lw=2, zorder=4, solid_capstyle="round")
    ax.text(x, med + 4, f"{med:.0f}", ha="center", va="bottom", fontsize=8, color=INK)
ax.set_xticks([0.5, 3.1], [f"extra-good units (tier 1)\n{ptxt('n_SS_t1')}",
                           f"other units ≥ 20 Hz, any label (tier 2)\n{ptxt('n_SS_t2')}"])
ax.set_ylabel("SS-partner candidates per session (CB)")
ax.set_xlim(-0.6, 4.2)
clean(ax)
for so in (KS, IB):
    ax.scatter([], [], s=34, color=COL[so], label=NAME[so])
ax.legend(frameon=False, fontsize=8, loc="upper left")
title(ax, "A  SS-partner pool per session",
      f"tier 2 higher with IBL across sessions (MW p = {ST.at['n_SS_t2', 'p_mannwhitney']:.2g}), not within mouse")

ax = axs[1]
d = ses.dropna(subset=["frac_paired_any"])
big = d[d.n_CS >= MIN_CS_SCATTER]
for _, r in d.iterrows():
    filled = r.n_CS >= MIN_CS_SCATTER
    ax.scatter(r.n_SS_t2, r.frac_paired_any, marker=MARK[r.mouse], s=18 + 2.2 * r.n_CS,
               color=COL[r.sorter] if filled else "white", edgecolor=COL[r.sorter] if not filled else "white",
               lw=1.2 if not filled else 0.8, zorder=3)
rho = spearmanr(big.n_SS_t2, big.frac_paired_any)
ax.set_xlabel("tier-2 SS-partner candidates in the session (CB, ≥ 20 Hz, any label)")
ax.set_ylabel("fraction of CS with an SS partner")
ax.set_ylim(-0.04, 1.06)
clean(ax)
for so in (KS, IB):
    ax.scatter([], [], s=40, color=COL[so], label=NAME[so])
ax.scatter([], [], s=40, color="white", edgecolor=INK2, label=f"< {MIN_CS_SCATTER} CS (not in ρ)")
ax.legend(frameon=False, fontsize=8, loc="upper right")
ax.text(0.99, 0.7, "marker size ∝ number of CS", transform=ax.transAxes, ha="right", fontsize=7.5, color=INK2)
title(ax, "B  More SS-like units → more CS paired",
      f"Spearman ρ = {rho.correlation:.2f}, p = {rho.pvalue:.2g} (sessions with ≥ {MIN_CS_SCATTER} CS, n = {len(big)})")
mouse_legend(fig, y=0.0)
save(fig, "fig2_why_ss_pool")

# %% 6. FIGURE 3 — the unpaired CS: what is missing, and is it duplication?
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [1.5, 1]})
fig.subplots_adjust(left=0.1, right=0.97, top=0.76, bottom=0.2, wspace=0.35)
fig.suptitle("Unpaired CS: mostly no unit pauses at all; CS spikes in a neighbour more often with KS4",
             x=0.07, ha="left", fontsize=12, fontweight="bold", color=INK, y=0.97)
ax = axs[0]
reason_order = [("paired", "paired"), ("near_miss", "near miss"),
                ("weak_pause(0.7-0.85)", "weak pause (0.7–0.85)"), ("no_pause_anywhere", "no pause (> 0.85)")]
alphas = [1.0, 0.65, 0.38, None]
out = []
for yk, so in enumerate((IB, KS)):
    c = c3[c3.sorter == so]
    r = c.unpaired_reason.str.replace(r"near_miss:.*", "near_miss", regex=True)
    left = 0
    for (key, lab), a in zip(reason_order, alphas):
        n = int((r == key).sum())
        w = n / len(c)
        out.append(dict(sorter=so, outcome=lab, n=n, frac=w))
        ax.barh(yk, w, left=left, height=0.62, color=COL[so] if a else "#6f6e69", alpha=a if a else 1,
                edgecolor="white", lw=2)
        if w > 0.05:
            ax.text(left + w / 2, yk, f"{w:.0%}", ha="center", va="center", fontsize=8.5,
                    color="white" if (a is None or a > 0.5) else INK)
        left += w
ax.set_yticks([0, 1], [f"{NAME[IB]}\n(n = {(c3.sorter == IB).sum()})", f"{NAME[KS]}\n(n = {(c3.sorter == KS).sum()})"])
ax.set_xlim(0, 1)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1], ["0", "25%", "50%", "75%", "100%"])
ax.spines[["top", "right"]].set_visible(False)
x0 = 0
for (key, lab), a in zip(reason_order, alphas):
    ax.add_patch(plt.Rectangle((x0, -0.8), 0.025, 0.18, color=INK2 if a else "#6f6e69", alpha=a if a else 1,
                               clip_on=False))
    ax.text(x0 + 0.032, -0.71, lab, fontsize=7.3, va="center", color=INK2)
    x0 += 0.032 + 0.0095 * len(lab) + 0.04
ax.set_ylim(-0.95, 1.45)
title(ax, "A  Outcome for each CS (lowest SS pause ratio anywhere on the probe)",
      "near miss = pause < 0.7 but one other criterion fails; categories are descriptive, not decisions")
pd.DataFrame(out).to_csv(TAB / "cs_outcome_pooled.csv", index=False)

ax = axs[1]
bars = [("frac_CS_shared", "CS with > 20% of spikes\nalso in a nearby unit"),
        ("frac_CS_dup_collapsed", "CS removed as\nsame-cell duplicates")]
pooled = {so: [(zl[zl.sorter == so].best_excess > 0.2).mean(),
               1 - ses.loc[ses.sorter == so, "n_CS_distinct_est"].sum() / ses.loc[ses.sorter == so, "n_CS"].sum()]
          for so in (KS, IB)}
for j, so in enumerate((KS, IB)):
    xs = np.arange(len(bars)) + (j - 0.5) * 0.36
    ax.bar(xs, pooled[so], width=0.34, color=COL[so], label=NAME[so], edgecolor="white", lw=1.5)
    for x, v in zip(xs, pooled[so]):
        ax.text(x, v + 0.004, f"{v:.0%}", ha="center", va="bottom", fontsize=8, color=INK)
ax.set_xticks(range(len(bars)), [b[1] for b in bars], fontsize=8)
ax.set_ylabel("share of CS units (pooled)")
ax.set_ylim(0, max(max(v) for v in pooled.values()) * 1.35 + 0.01)
ax.legend(frameon=False, fontsize=8, loc="upper right")
clean(ax)
title(ax, "B  Leakage higher with KS4; duplicates rare",
      f"session level: shared {ptxt('frac_CS_shared')}\nduplicates {ptxt('frac_CS_dup_collapsed')}")
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))
save(fig, "fig3_unpaired_and_duplicates")

# %% 7. FIGURE 4 — unit properties: firing statistics equal, template representation differs
fig, axs = plt.subplots(1, 4, figsize=(12.5, 4.4))
fig.subplots_adjust(left=0.06, right=0.98, top=0.74, bottom=0.25, wspace=0.55)
fig.suptitle("Unit properties: CS firing nearly the same; IBL templates are about twice as wide on the probe "
             "(how the sorter represents units)", x=0.06, ha="left", fontsize=11.5, fontweight="bold", color=INK, y=0.97)
for ax, (col, ylab, head) in zip(axs, [("CS_fr_med", "Hz", "CS firing rate"), ("CS_cv2_med", "CV2", "CS regularity (CV2)"),
                                       ("CS_footprint_span_um_med", "µm", "CS template footprint"),
                                       ("SS_footprint_span_um_med", "µm", "SS template footprint")]):
    strip(ax, col, ylab, (lambda v: f"{v:.0f}") if "um" in col else (lambda v: f"{v:.2f}"))
    ax.set_xticks([0, 1], [NAME[KS], NAME[IB]], fontsize=8)
    title(ax, head, ptxt(col))
mouse_legend(fig)
save(fig, "fig4_unit_properties")

# %% 7b. FIGURE 5 — simple-spike units: how many per session, how many searched, how many paired
fig, axs = plt.subplots(1, 3, figsize=(12.5, 4.9))
fig.subplots_adjust(left=0.06, right=0.98, top=0.74, bottom=0.25, wspace=0.42)
fig.suptitle("Simple-spike units per session: few C4 PkC_ss in either sorter; IBL tends higher on all three counts, "
             "none significant per session", x=0.06, ha="left", fontsize=11.5, fontweight="bold", color=INK, y=0.97)
nss = {so: int(ses.loc[ses.sorter == so, "n_SS"].sum()) for so in (KS, IB)}
for ax, (col, ylab, head, sub) in zip(axs, [
        ("n_SS", "PkC_ss units per session", "A  SS identified by C4",
         f"extra-good CB, C4 PkC_ss; total {nss[KS]} vs {nss[IB]}\n{ptxt('n_SS')}"),
        ("n_SS_pool", "units searched as SS partner", "B  SS-partner candidates searched",
         f"extra-good + other units ≥ 20 Hz (CB)\n{ptxt('n_SS_pool')}"),
        ("n_SS_paired", "SS units paired with a CS", "C  SS units paired with a CS",
         f"distinct best partners of paired CS\n{ptxt('n_SS_paired')}")]):
    strip(ax, col, ylab, lambda v: f"{v:g}")
    ax.set_xticks([0, 1], [f"{NAME[KS]}\n{ses.loc[ses.sorter == KS, col].notna().sum()} sessions",
                           f"{NAME[IB]}\n{ses.loc[ses.sorter == IB, col].notna().sum()} sessions"], fontsize=8)
    title(ax, head, sub)
mouse_legend(fig)
save(fig, "fig5_ss_counts")
print("SS counts per session (distribution):")
print(ses.groupby("sorter")[["n_SS", "n_SS_pool", "n_SS_paired"]].describe().T.round(1).to_string())

# %% 7c. FIGURE 6 — label gate: old rule (tier 2 = cluster_group good/mua, v3) vs no label gate (v5)
s_prev = pd.read_csv(CMP / STEP0_PREV / "session_summary.csv")
c_prev = pd.read_csv(CMP / STEP0_PREV / "cs_units.csv")
cmpv = s_prev[["mouse", "S", "sorter", "n_SS_t2", "frac_paired_any"]].merge(
    s3[["mouse", "S", "n_SS_t2", "frac_paired_any"]], on=["mouse", "S"], suffixes=("_gated", "_nogate"))
cmpv.to_csv(TAB / "label_gate_v3_vs_v5_sessions.csv", index=False)
fig, axs = plt.subplots(1, 3, figsize=(14, 4.9), gridspec_kw={"width_ratios": [1, 1, 1.35]})
fig.subplots_adjust(left=0.055, right=0.98, top=0.74, bottom=0.2, wspace=0.45)
fig.suptitle("Labels do not explain the sorter difference: tier 2 gated on curated good/mua (v3) vs any label (v5)",
             x=0.055, ha="left", fontsize=11.5, fontweight="bold", color=INK, y=0.97)
fig.text(0.055, 0.9, "KS4 cluster_group is the Phy-curated label (we never set mua; it mostly inherits Kilosort's "
         "automatic KSLabel); IBL cluster_group has an unclear origin. v5 ignores both.", fontsize=8.5, color=INK2)
for ax, (col, ylab, head) in zip(axs[:2], [("n_SS_t2", "tier-2 SS candidates per session", "A  Tier-2 pool"),
                                       ("frac_paired_any", "fraction of CS with an SS partner", "B  CS paired")]):
    for k, so in enumerate((KS, IB)):
        d = cmpv[cmpv.sorter == so].dropna(subset=[f"{col}_gated", f"{col}_nogate"])
        for _, r in d.iterrows():
            ax.plot([2 * k, 2 * k + 0.8], [r[f"{col}_gated"], r[f"{col}_nogate"]], color=COL[so], lw=0.8, alpha=0.45)
        for j, suf in enumerate(("gated", "nogate")):
            med = d[f"{col}_{suf}"].median()
            ax.scatter(2 * k + 0.8 * j + np.zeros(len(d)), d[f"{col}_{suf}"], marker="o", s=16, color=COL[so],
                       edgecolor="white", lw=0.6, zorder=3)
            ax.plot([2 * k + 0.8 * j - 0.22, 2 * k + 0.8 * j + 0.22], [med, med], color=INK, lw=2, zorder=4)
            ax.text(2 * k + 0.8 * j + 0.26, med, f"{med:.2g}", va="center", fontsize=8, color=INK)
    ax.set_xticks([0, 0.8, 2, 2.8], ["KS4\ngood/mua", "KS4\nany label", "IBL\ngood/mua", "IBL\nany label"], fontsize=8)
    ax.set_xlim(-0.5, 3.4)
    ax.set_ylabel(ylab)
    clean(ax)
    pk = {so: ((c_prev.sorter == so) & c_prev.paired_any).sum() / (c_prev.sorter == so).sum() for so in (KS, IB)}
    pn = {so: ((c3.sorter == so) & c3.paired_any).sum() / (c3.sorter == so).sum() for so in (KS, IB)}
    sub = (f"lines = sessions; medians shown" if col == "n_SS_t2" else
           f"pooled: KS4 {pk[KS]:.0%} → {pn[KS]:.0%}, IBL {pk[IB]:.0%} → {pn[IB]:.0%}")
    title(ax, head, sub)

# C: what the v5 tier-2 pool is made of, by curated label and by the automatic KSLabel (pooled units)
ax = axs[2]
LAB = [("good", 1.0), ("mua", 0.55), ("noise", 0.3), ("unsorted", 0.15), ("none", 0.06)]   # none = no KSLabel (unit made in Phy)
comp = []
rows_c = [(IB, "kslabel", "IBL KSLabel\n(automatic)"), (IB, "group", "IBL cluster_group"),
          (KS, "kslabel", "KS4 KSLabel\n(automatic)"), (KS, "group", "KS4 cluster_group\n(Phy-curated)")]
for yk, (so, src, lab) in enumerate(rows_c):
    p = pool2[pool2.sorter == so]
    left = 0.0
    for name, a in LAB:
        n = int((p[src] == name).sum())
        w = n / len(p)
        comp.append(dict(sorter=so, label_source=src, label=name, n=n, frac=w))
        if n:
            ax.barh(yk, w, left=left, height=0.62, color=COL[so], alpha=a, edgecolor="white", lw=2)
            if w > 0.06:
                ax.text(left + w / 2, yk, f"{n}", ha="center", va="center", fontsize=8,
                        color="white" if a > 0.5 else INK)
        left += w
    ax.text(1.01, yk, f"n = {len(p)}", va="center", fontsize=8, color=INK)
ax.set_yticks(range(4), [r[2] for r in rows_c], fontsize=8)
ax.set_xlim(0, 1.16)
ax.set_xticks([0, 0.5, 1], ["0", "50%", "100%"])
ax.set_xlabel("share of tier-2 SS candidates (pooled over sessions)")
ax.spines[["top", "right"]].set_visible(False)
for j, (name, a) in enumerate(LAB):
    ax.add_patch(plt.Rectangle(([0, 0.17, 0.34, 0.53, 0.75][j], -1.05), 0.03, 0.25, color=INK2, alpha=a, clip_on=False))
    ax.text([0, 0.17, 0.34, 0.53, 0.75][j] + 0.04, -0.93, "no KSLabel" if name == "none" else name, fontsize=7.5, va="center", color=INK2)
ax.set_ylim(-1.2, 3.5)
ks_g = ST.at["frac_t2_ksgood", "ks4_median"]
ib_g = ST.at["frac_t2_ksgood", "ibl_median"]
title(ax, "C  What the tier-2 pool is made of (v5)",
      f"KSLabel good per session: KS4 {ks_g:.0%} vs IBL {ib_g:.0%}; {ptxt('frac_t2_ksgood')}\n"
      "numbers = units; the label never selects units")
pd.DataFrame(comp).to_csv(TAB / "tier2_pool_label_composition.csv", index=False)
pool2.to_csv(TAB / "tier2_pool_units_step0v5.csv", index=False)
save(fig, "fig6_label_gate_v3_vs_v5")

# %% 8. TABLES (copies of the current per-session / per-CS results) + README
ses.to_csv(TAB / "sessions_final.csv", index=False)
c3.to_csv(TAB / "cs_units_step0v5.csv", index=False)
for src, dst in [(CMP / "session_table.csv", "session_table.csv"), (CMP / "units_all.csv", "units_pass1_all_extragood_CB.csv"),
                 (CMP / STEP0_SUMMARY / "cs_cs_evidence_all.csv", "cs_cs_evidence_step0v5.csv"),
                 (CMP / ZERO_LAG / "cs_neighbour_pairs.csv", "cs_neighbour_zero_lag_pairs.csv")]:
    shutil.copy2(src, TAB / dst)
print("tables in", TAB)
