# -*- coding: utf-8 -*-
"""
ephaptic_summary.py — across all sessions: does a CS briefly suppress the SS of NEIGHBOURING cells (0–2 ms), falling off
with distance, as in Han, Chen et al. 2020 (Nat Neurosci 23:1399, Fig. 1d–e)?

Input : CS_curation/<session>/pause_figure/ephaptic_pairs.csv and pause_figure.js (from pause_figure_data.py)
Output: CS_curation/_summary/ephaptic_all_sessions.png   (pooled + per animal: 0–2 and 3–8 ms ratios by distance group,
                                                           and the n-weighted mean rate curves per distance group)
        CS_curation/_summary/ephaptic_all_sessions.csv   (per session × group: n, median, IQR, Wilcoxon p vs 1)
Non-partner, non-flagged pairs only (partner = candidate of that CS; flagged = cs_share_flag / cs_event_flag).
Note: pairs are not independent (each SS appears with every CS of its session); p values are descriptive only.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

CUR = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation")
OUT = CUR / "_summary"
GROUPS = [(-1, 25, "≤25 µm"), (25, 50, "25–50 µm"), (50, 75, "50–75 µm"), (75, 1e9, ">75 µm")]
COLS = ["#b3261e", "#e08a2c", "#3380cc", "#8c8c8c"]

rows, curves = [], []
for f in sorted(CUR.glob("*/pause_figure/ephaptic_pairs.csv")):
    sess = f.parent.parent.name
    e = pd.read_csv(f).assign(session=sess, animal=sess.rsplit("_", 1)[0])
    rows.append(e)
    js = (f.parent / "pause_figure.js").read_text(encoding="utf-8")
    pf = json.loads(js[js.index("=") + 1:].rstrip().rstrip(";"))
    for g in pf["eph_session"]:
        if g["mean"]:
            curves.append(dict(session=sess, animal=sess.rsplit("_", 1)[0], group=g["label"], n=g["n"], mean=np.array(g["mean"], float)))
    t0, binms = pf["group_t0_ms"], pf["bin_ms"]
E = pd.concat(rows, ignore_index=True)
E["group"] = pd.cut(E.dist_um, [g[0] for g in GROUPS] + [GROUPS[-1][1]], labels=[g[2] for g in GROUPS])
N = E[~E.partner & ~E.flagged]


def stats(x):
    x = x.dropna()
    p = wilcoxon(x - 1).pvalue if len(x) >= 10 else np.nan
    return pd.Series(dict(n=len(x), median=x.median(), q25=x.quantile(.25), q75=x.quantile(.75), p_vs_1=p))


tab = []
for key, d in [("all", N)] + [(a, g) for a, g in N.groupby("animal")] + [(s, g) for s, g in N.groupby("session")]:
    for grp, gd in d.groupby("group", observed=True):
        for col in ("r_0_2", "r_3_8"):
            tab.append(dict(scope=key, group=grp, measure=col, **stats(gd[col])))
T = pd.DataFrame(tab)
OUT.mkdir(exist_ok=True)
T.round(4).to_csv(OUT / "ephaptic_all_sessions.csv", index=False)

scopes = ["all"] + sorted(N.animal.unique())
fig, ax = plt.subplots(3, len(scopes), figsize=(3.2 * len(scopes), 9), squeeze=False)
C = pd.DataFrame(curves)
for j, sc in enumerate(scopes):
    d = N if sc == "all" else N[N.animal == sc]
    for i, col in enumerate(("r_0_2", "r_3_8")):
        a = ax[i, j]
        data = [d.loc[d.group == g[2], col].dropna().values for g in GROUPS]
        bp = a.boxplot(data, showfliers=False, widths=0.6, patch_artist=True)
        for patch, c in zip(bp["boxes"], COLS):
            patch.set_facecolor(c); patch.set_alpha(0.35)
        a.axhline(1, color="k", lw=0.7, ls="--")
        a.set_xticks(range(1, 5), [g[2] for g in GROUPS], rotation=30, fontsize=8)
        for k, x in enumerate(data):
            a.text(k + 1, a.get_ylim()[1], f"n {len(x)}", ha="center", va="top", fontsize=7)
        a.set_ylabel(("0–2 ms" if col == "r_0_2" else "3–8 ms") + " ÷ baseline")
        if i == 0:
            a.set_title(f"{sc} ({d.session.nunique()} sessions)", fontsize=10)
    a = ax[2, j]
    c = C if sc == "all" else C[C.animal == sc]
    for g, col in zip(GROUPS, COLS):
        cg = c[c.group == g[2]]
        if len(cg):
            m = np.average(np.vstack(cg["mean"].values), axis=0, weights=cg.n.values)
            x = t0 + binms * (np.arange(len(m)) + 0.5)
            a.plot(x, m, color=col, lw=1.4, label=f"{g[2]} (n {cg.n.sum()})")
    a.axhline(1, color="k", lw=0.7, ls="--"); a.axvspan(0, 2, color="#f6e3e1", zorder=0); a.axvspan(3, 8, color="#e8e8e8", zorder=0)
    a.set_xlabel("time from CS (ms)"); a.set_ylabel("SS rate ÷ baseline"); a.legend(fontsize=6, loc="lower right")
fig.suptitle("Ephaptic check: SS rate of NON-partner SS after a CS, by CS–SS peak-channel distance (Han, Chen et al. 2020, Fig. 1d–e)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "ephaptic_all_sessions.png", dpi=130)
print(T[(T.scope.isin(scopes))].pivot_table(index=["scope", "group"], columns="measure", values="median", observed=True).round(3).to_string())
print("n non-partner pairs", len(N), "| sessions", N.session.nunique())
