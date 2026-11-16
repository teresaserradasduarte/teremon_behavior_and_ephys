# -*- coding: utf-8 -*-
"""
sorter_compare_step0v3.py — KS4 vs IBL comparison built on the step0 v3 (stable-span) per-session outputs, plus an inspection
of the sessions where most CS stay unpaired.

Input: CS_curation/<mouse>_<R#>/step0_explore_bhvspan_v3batch/tables/ (run_step0_batch.py, frozen step0 v3 snapshot):
per-unit stable spans (win_selected = reference, extended while FR/amplitude stay stable; a pair uses the
  overlap of its two spans; *_ref columns = win_selected only), plus everything in v2:
  - SS partner search in tier 1 (extra-good CB) + tier 2 (good/mua CB, FR >= 20 Hz), jitter test on EVERY pair,
    per-600 s block pause ratios, Lussac reference criterion;
  - CS–CS duplicate evidence lines (coinc, lag_sd, union, ss_shared, ss_cond, footprint, covar), type hints,
    evidence status and families. Nothing is merged.
Nothing here decides or merges. "Distinct-CS estimate" collapses only pairs whose evidence status is
"converges_same_cell" and is reported as a descriptive upper bound on duplication, not a curation call.

Outputs: CS_curation/_sorter_comparison/step0_<BATCH_TAG>/ (tables + figs)
Batch tag: env STEP0_BATCH_TAG (default "v5batch", shared with run_step0_batch.py and sorter_compare_zero_lag_check.py).
The first v3 results were written to _sorter_comparison/step0v3/ (input step0_explore_bhvspan_v3batch); kept as is.
"""

# %% 1. PARAMETERS
from pathlib import Path
import json
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

LOCAL = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation")
DROPBOX_ROOT = Path(r"D:\Learning Lab Dropbox")
BATCH_TAG = os.environ.get("STEP0_BATCH_TAG", "v5batch")
RUN_DIR = f"step0_explore_bhvspan_{BATCH_TAG}"   # per-session step0 outputs written by run_step0_batch.py
OUT = LOCAL / "_sorter_comparison" / f"step0_{BATCH_TAG}"
assert DROPBOX_ROOT not in OUT.resolve().parents
(OUT / "figs").mkdir(parents=True, exist_ok=True)
INSPECT = [("5_FerreroRocher", "S6"), ("3_Toblerone", "S5")]
STRONG_PR, STRONG_DUR = 0.3, 5.0      # "strong" pause: ratio < 0.3 and 50 %-recovery >= 5 ms
STABLE_BLOCK_PR = 0.7                 # pair "stable" if every 600 s block has pause ratio < this
NEAR_MISS_PR = 0.7                    # unpaired CS whose best pair still has pause ratio < this = near-miss
NO_PAUSE_PR = 0.85                    # unpaired CS whose best pair has pause ratio > this = no pause anywhere
NEAR_UM = 100.0
N_PERM = 10000
MICE = ["1_CoteDor", "2_Lindt", "3_Toblerone", "4_Milka", "5_FerreroRocher"]
SORTER_COLOR = {"ks4": "#2a78d6", "ibl": "#eb6834"}
SORTER_NAME = {"ks4": "Kilosort4", "ibl": "IBL sorter"}
MOUSE_MARKER = dict(zip(MICE, ["o", "s", "^", "D", "v"]))

sessions = pd.read_csv(LOCAL / "_sorter_comparison" / "session_table.csv")
_u1 = pd.read_csv(LOCAL / "_sorter_comparison" / "units_all.csv")
n_cs_pass1 = _u1[_u1["c4_label"] == "PkC_cs"].groupby(["mouse", "S"]).size()


# %% 2. LOAD STEP0 OUTPUTS
def uf_groups(ids, edges):
    parent = {i: i for i in ids}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for a, b in edges:
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    g = {}
    for i in ids:
        g.setdefault(find(i), []).append(i)
    return list(g.values())


cs_rows, ses_rows, ev_all, pairs_all, units_all = [], [], [], [], []
for r in sessions.itertuples():
    td = LOCAL / f"{r.mouse}_{r.R}" / RUN_DIR / "tables"
    if not (td / "cs_ss_pairs.csv").exists():
        if n_cs_pass1.get((r.mouse, r.S), 0) == 0:   # no CS: step0 not run; counts still enter the session table
            ses_rows.append(dict(mouse=r.mouse, S=r.S, sorter=r.sorter, n_CS=0))
        else:
            print(f"missing step0 output: {r.mouse} {r.S}")
        continue
    units = pd.read_csv(td / "unit_table.csv").set_index("phyID")
    pairs = pd.read_csv(td / "cs_ss_pairs.csv") if (td / "cs_ss_pairs.csv").stat().st_size > 5 else pd.DataFrame()
    try:
        ev = pd.read_csv(td / "cs_cs_evidence.csv")
    except pd.errors.EmptyDataError:
        ev = pd.DataFrame()
    prm = json.loads((td.parent / "params_used.json").read_text())
    cs_ids = units.index[units["role"] == "CS"].tolist()
    tag = dict(mouse=r.mouse, S=r.S, sorter=r.sorter)
    units_all.append(units.reset_index().assign(**tag))
    if len(pairs):
        pairs_all.append(pairs.assign(**tag))
    if len(ev):
        ev_all.append(ev.assign(**tag))

    # per CS
    same_edges = []
    if len(ev):
        same_edges = [(int(a), int(b)) for a, b, s in zip(ev.cs_a, ev.cs_b, ev.evidence_status)
                      if s == "converges_same_cell"]
    groups = uf_groups(cs_ids, same_edges)
    group_of = {c: min(g) for g in groups for c in g}
    fam_of = {}
    if len(ev) and "family" in ev:
        for a, b, f in zip(ev.cs_a, ev.cs_b, ev.family):
            if pd.notna(f):
                fam_of[int(a)] = fam_of[int(b)] = int(f)
    ss_pos = units.loc[units["role"].isin(["SS_cand", "SS_cand2"]), ["x", "y"]]
    for c in cs_ids:
        p = pairs[pairs["cs"] == c] if len(pairs) else pd.DataFrame()
        cand = p[p["candidate"]] if len(p) else p
        best = cand.sort_values("rank").iloc[0] if len(cand) else None
        top = p.sort_values("pause_ratio").iloc[0] if len(p) else None   # lowest pause ratio, candidate or not
        stable = None
        if best is not None:
            stable = bool(best["block_pr_max"] < STABLE_BLOCK_PR) if pd.notna(best["block_pr_max"]) else False
        d_near = np.hypot(ss_pos["x"] - units.at[c, "x"], ss_pos["y"] - units.at[c, "y"])
        row = dict(**tag, cs=c, fr=units.at[c, "fr"], n_spikes=units.at[c, "n_spikes"], y_um=units.at[c, "y"],
                   x_um=units.at[c, "x"], frac_isi_lt1ms=units.at[c, "frac_isi_lt1ms"],
                   paired_any=best is not None,
                   paired_ref=bool(len(p) and "candidate_ref" in p and p["candidate_ref"].fillna(False).astype(bool).any()),
                   span_min=(p["win_s"].median() / 60 if len(p) and "win_s" in p else np.nan), paired_t1=bool(len(cand) and (cand["ss_tier"] == 1).any()),
                   n_cand=len(cand), lussac_any=bool(len(p) and p["lussac_pass"].any()),
                   best_ss=best["ss"] if best is not None else np.nan,
                   best_tier=best["ss_tier"] if best is not None else np.nan,
                   best_group=best["ss_group"] if best is not None else "",
                   best_c4=best["ss_c4"] if best is not None else "",
                   best_pr=best["pause_ratio"] if best is not None else np.nan,
                   best_dur50=best["pause_dur50_ms"] if best is not None else np.nan,
                   best_jitter_z=best["jitter_z"] if best is not None else np.nan,
                   best_dist=best["dist_um"] if best is not None else np.nan,
                   best_stable=stable,
                   top_pr=top["pause_ratio"] if top is not None else np.nan,
                   top_asym=top["asym"] if top is not None else np.nan,
                   top_jitter_p=top["jitter_p"] if top is not None else np.nan,
                   top_collision=top["collision_dip"] if top is not None else np.nan,
                   top_zero_peak=top["zero_peak_ratio"] if top is not None else np.nan,
                   top_dist=top["dist_um"] if top is not None else np.nan,
                   max_zero_peak=p["zero_peak_ratio"].max() if len(p) else np.nan,
                   n_ss_within_100um=int((d_near <= NEAR_UM).sum()),
                   family=fam_of.get(c, np.nan), same_cell_group=group_of[c],
                   same_cell_group_size=len([g for g in groups if c in g][0]))
        row["strong"] = bool(row["paired_any"] and row["best_pr"] < STRONG_PR and row["best_dur50"] >= STRONG_DUR)
        cs_rows.append(row)

    # per session
    n_t1 = int((units["role"] == "SS_cand").sum())
    n_t2 = int((units["role"] == "SS_cand2").sum())
    c_s = pd.DataFrame([x for x in cs_rows if x["mouse"] == r.mouse and x["S"] == r.S])
    cand_all = pairs[pairs["candidate"]] if len(pairs) else pd.DataFrame(columns=["cs", "ss"])
    ses_rows.append(dict(**tag, win_min=prm["window_len_s"] / 60, n_CS=len(cs_ids), n_SS_t1=n_t1, n_SS_t2=n_t2,
                         n_pairs=len(pairs),
                         frac_paired_t1=c_s["paired_t1"].mean() if len(c_s) else np.nan,
                         frac_paired_any=c_s["paired_any"].mean() if len(c_s) else np.nan,
                         frac_paired_ref=c_s["paired_ref"].mean() if len(c_s) else np.nan,
                         span_med_min=c_s["span_min"].median() if len(c_s) else np.nan,
                         frac_strong=c_s["strong"].mean() if len(c_s) else np.nan,
                         frac_lussac=c_s["lussac_any"].mean() if len(c_s) else np.nan,
                         frac_paired_stable=(c_s["best_stable"] == True).mean() if len(c_s) else np.nan,  # noqa: E712
                         n_SS_claimed_multi=int((cand_all.groupby("ss")["cs"].nunique() > 1).sum()) if len(cand_all) else 0,
                         n_CS_multi_SS=int((cand_all.groupby("cs")["ss"].nunique() > 1).sum()) if len(cand_all) else 0,
                         n_dup_cand=int(ev["dup_candidate"].sum()) if len(ev) else 0,
                         n_converge_same=int((ev["evidence_status"] == "converges_same_cell").sum()) if len(ev) else 0,
                         n_CS_in_family=int(c_s["family"].notna().sum()) if len(c_s) else 0,
                         n_CS_distinct_est=len(groups),
                         best_pr_med=c_s.loc[c_s["paired_any"], "best_pr"].median() if len(c_s) else np.nan,
                         best_dur_med=c_s.loc[c_s["paired_any"], "best_dur50"].median() if len(c_s) else np.nan))

cs = pd.DataFrame(cs_rows)
ses = pd.DataFrame(ses_rows)
ev_all = pd.concat(ev_all, ignore_index=True) if ev_all else pd.DataFrame()
pairs_all = pd.concat(pairs_all, ignore_index=True)
ses["frac_CS_dup_collapsed"] = 1 - ses["n_CS_distinct_est"] / ses["n_CS"].replace(0, np.nan)
cs.to_csv(OUT / "cs_units.csv", index=False)
ses.to_csv(OUT / "session_summary.csv", index=False)
if len(ev_all):
    ev_all.to_csv(OUT / "cs_cs_evidence_all.csv", index=False)
print(f"{len(ses)} sessions, {len(cs)} CS units")


# %% 3. SORTER COMPARISON
def within_mouse_perm(df, col, n_perm=N_PERM, seed=0):
    d = df[["mouse", "sorter", col]].dropna()
    groups = [g for _, g in d.groupby("mouse") if g["sorter"].nunique() == 2]
    if not groups:
        return np.nan, np.nan, 0
    labs = [g["sorter"].values for g in groups]

    def stat(ls):
        return np.mean([np.median(g[col].values[l == "ks4"]) - np.median(g[col].values[l == "ibl"])
                        for g, l in zip(groups, ls)])
    obs = stat(labs)
    rng = np.random.default_rng(seed)
    null = np.array([stat([rng.permutation(l) for l in labs]) for _ in range(n_perm)])
    return obs, (1 + np.sum(np.abs(null) >= abs(obs) - 1e-12)) / (1 + n_perm), len(groups)


with pd.option_context("display.width", 260, "display.max_columns", 40, "display.max_rows", 100):
    print("--- per session")
    print(ses[["mouse", "S", "sorter", "win_min", "span_med_min", "n_CS", "n_SS_t1", "n_SS_t2", "frac_paired_ref",
               "frac_paired_t1", "frac_paired_any", "frac_strong",
               "frac_paired_stable", "frac_lussac", "n_SS_claimed_multi", "n_dup_cand", "n_converge_same",
               "n_CS_distinct_est"]].round(2).to_string(index=False))
    print("--- pooled CS units per sorter")
    print(cs.groupby("sorter").agg(n=("cs", "size"), paired_ref=("paired_ref", "mean"), paired_t1=("paired_t1", "mean"),
                                   paired_any=("paired_any", "mean"),
                                   strong=("strong", "mean"), stable=("best_stable", lambda v: (v == True).mean()),  # noqa: E712
                                   lussac=("lussac_any", "mean"), in_family=("family", lambda v: v.notna().mean()))
          .round(3).to_string())
    print("--- totals per sorter")
    print(ses.groupby("sorter")[["n_CS", "n_CS_distinct_est", "n_dup_cand", "n_converge_same", "n_SS_claimed_multi"]].sum())
    pa = cs[cs["paired_any"]]
    print("--- best partner tier / group (paired CS)")
    print(pd.crosstab(pa["sorter"], [pa["best_tier"], pa["best_group"]]))
    print("--- best partner pause (paired CS), medians")
    print(pa.groupby("sorter")[["best_pr", "best_dur50", "best_jitter_z", "best_dist"]].median().round(2))
    if len(ev_all):
        print("--- CS–CS evidence: status and type per sorter (pairs = dup candidates + sharp-peak-only)")
        print(pd.crosstab(ev_all["sorter"], ev_all["evidence_status"]))
        print(pd.crosstab(ev_all["sorter"], ev_all["type_hint"]).T)

st = []
for col in ["n_CS", "n_SS_t2", "frac_paired_ref", "frac_paired_t1", "frac_paired_any", "frac_strong", "frac_paired_stable", "frac_lussac",
            "frac_CS_dup_collapsed", "best_pr_med", "best_dur_med"]:
    k, i = ses.loc[ses["sorter"] == "ks4", col].dropna(), ses.loc[ses["sorter"] == "ibl", col].dropna()
    if len(k) < 2 or len(i) < 2:
        continue
    d, p, nm = within_mouse_perm(ses, col)
    st.append(dict(metric=col, ks4_median=k.median(), ibl_median=i.median(), n_ks4=len(k), n_ibl=len(i),
                   p_mannwhitney=mannwhitneyu(k, i).pvalue, within_mouse_diff=d, p_within_mouse=p, n_mice=nm))
st = pd.DataFrame(st)
st.to_csv(OUT / "stats_session_level.csv", index=False)
print("--- session-level stats (uncorrected, exploratory)")
print(st.round(4).to_string(index=False))


# %% 4. INSPECTION: why do CS stay unpaired? (every session; detailed for INSPECT)
def reason(r):
    if r["paired_any"]:
        return "paired"
    if r["same_cell_group_size"] > 1 and cs.loc[(cs["mouse"] == r["mouse"]) & (cs["S"] == r["S"])
                                                & (cs["same_cell_group"] == r["same_cell_group"]), "paired_any"].any():
        return "duplicate_of_paired_CS"
    if pd.isna(r["top_pr"]):
        return "no_SS_candidates"
    if r["top_pr"] < NEAR_MISS_PR:
        if r["top_collision"] == True:  # noqa: E712
            return "near_miss:collision_dip"
        if r["top_asym"] <= 0.2:
            return "near_miss:symmetric_dip"
        if r["top_jitter_p"] >= 0.05:
            return "near_miss:jitter_ns"
        return "near_miss:other"
    if r["top_pr"] > NO_PAUSE_PR:
        return "no_pause_anywhere"
    return "weak_pause(0.7-0.85)"


cs["unpaired_reason"] = cs.apply(reason, axis=1)
cs.to_csv(OUT / "cs_units.csv", index=False)
print("--- unpaired-CS reasons, per sorter (fraction of CS)")
print(pd.crosstab(cs["sorter"], cs["unpaired_reason"], normalize="index").round(3).T.to_string())

for m, s in INSPECT:
    c = cs[(cs["mouse"] == m) & (cs["S"] == s)]
    if not len(c):
        continue
    u = [x for x in units_all if x["mouse"].iloc[0] == m and x["S"].iloc[0] == s][0]
    print("=" * 80)
    print(f"INSPECTION {m} {s} ({c['sorter'].iloc[0]}): {len(c)} CS")
    print(c["unpaired_reason"].value_counts().to_string())
    print(f"  CS distinct estimate (converges_same_cell collapsed): {c['same_cell_group'].nunique()}; "
          f"CS in a duplicate family: {c['family'].notna().sum()}")
    print(f"  CS rate median {c['fr'].median():.2f} Hz; CS with an SS candidate within {NEAR_UM:.0f} µm: "
          f"{(c['n_ss_within_100um'] > 0).mean():.0%}; median zero-lag peak (max over SS) {c['max_zero_peak'].median():.2f}")
    print("  C4 labels of extra-good CB units:", u.loc[u["extragood"] & (u["region"] == "CB"), "c4_label"]
          .fillna("unlabelled").value_counts().to_dict())
    print("  depth (µm) CS:", np.percentile(c["y_um"], [0, 25, 50, 75, 100]).round(0).tolist(),
          "| SS tier1:", np.percentile(u.loc[u["role"] == "SS_cand", "y"], [0, 50, 100]).round(0).tolist()
          if (u["role"] == "SS_cand").any() else [],
          "| SS tier2:", np.percentile(u.loc[u["role"] == "SS_cand2", "y"], [0, 50, 100]).round(0).tolist()
          if (u["role"] == "SS_cand2").any() else [])
    with pd.option_context("display.width", 250, "display.max_rows", 120):
        print(c.sort_values("y_um")[["cs", "y_um", "fr", "paired_any", "best_ss", "best_tier", "best_pr", "best_dur50",
                                     "top_pr", "top_asym", "top_jitter_p", "top_dist", "max_zero_peak",
                                     "n_ss_within_100um", "family", "same_cell_group_size", "unpaired_reason"]]
              .round(2).to_string(index=False))

    # probe map: CS (paired / unpaired), SS candidates by tier, pairing lines
    fig, ax = plt.subplots(figsize=(6, 9))
    for role, mk, col, lab in (("SS_cand", "s", "#888888", "SS cand tier 1 (extra-good)"),
                               ("SS_cand2", ".", "#bbbbbb", "SS cand tier 2 (good/mua ≥20 Hz)")):
        q = u[u["role"] == role]
        ax.scatter(q["x"] + 4, q["y"], marker=mk, s=18, color=col, label=lab, zorder=2)
    for _, rr in c[c["paired_any"]].iterrows():
        sx, sy = u.set_index("phyID").loc[int(rr["best_ss"]), ["x", "y"]]
        ax.plot([rr["x_um"], sx + 4], [rr["y_um"], sy], color="#2a78d6", lw=0.8, zorder=1)
    ax.scatter(c.loc[c["paired_any"], "x_um"], c.loc[c["paired_any"], "y_um"], s=40, color="#2a78d6",
               edgecolor="white", label="CS paired", zorder=3)
    ax.scatter(c.loc[~c["paired_any"], "x_um"], c.loc[~c["paired_any"], "y_um"], s=40, facecolor="none",
               edgecolor="#eb6834", lw=1.5, label="CS unpaired", zorder=3)
    ax.set_xlabel("x (µm)  [SS markers drawn +4 µm right so they do not hide CS on the same site]")
    ax.set_ylabel("depth along probe (µm)")
    ax.set_title(f"{m} {s} ({SORTER_NAME[c['sorter'].iloc[0]]}): CS and SS candidates", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "figs" / f"inspect_{m}_{s}_probe_map.png", dpi=150)
    plt.close(fig)


# %% 5. FIGURES
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "font.size": 9, "savefig.dpi": 150})
rng = np.random.default_rng(0)
cols = [("frac_paired_ref", "CS paired, tier 1+2, win_selected only"), ("frac_paired_any", "CS paired, tier 1+2, stable spans"),
        ("frac_strong", f"CS strong pair (pr<{STRONG_PR}, ≥{STRONG_DUR:.0f} ms)"),
        ("frac_paired_stable", "CS pair stable in all 600 s blocks"), ("frac_CS_dup_collapsed",
                                                                       "CS removed if same-cell evidence collapsed"),
        ("n_SS_t2", "tier-2 SS candidates")]
r_st = st.set_index("metric")
fig, axs = plt.subplots(2, 3, figsize=(11, 7))
for ax, (col, t) in zip(axs.ravel(), cols):
    for xi, s in enumerate(["ks4", "ibl"]):
        sub = ses[ses["sorter"] == s].dropna(subset=[col])
        for mouse, g in sub.groupby("mouse"):
            ax.scatter(xi + rng.uniform(-0.18, 0.18, len(g)), g[col], s=28, marker=MOUSE_MARKER[mouse],
                       facecolor=SORTER_COLOR[s], edgecolor="white", linewidth=0.6, zorder=3)
        if len(sub):
            ax.plot([xi - 0.3, xi + 0.3], [sub[col].median()] * 2, color="#333333", lw=2, zorder=4)
    p = (f"MW p={r_st.at[col, 'p_mannwhitney']:.3f} | within-mouse p={r_st.at[col, 'p_within_mouse']:.3f}"
         if col in r_st.index else "")
    ax.set_title(f"{t}\n{p}", fontsize=8)
    ax.set_xticks([0, 1], [SORTER_NAME["ks4"], SORTER_NAME["ibl"]])
    ax.set_xlim(-0.6, 1.6)
h = [plt.Line2D([], [], marker=MOUSE_MARKER[m], ls="", color="#666666", label=m) for m in MICE]
fig.legend(handles=h, loc="lower center", ncol=5, frameon=False)
fig.suptitle("step0 v3 method (stable spans) on all sessions (one point per session)")
fig.tight_layout(rect=(0, 0.05, 1, 0.97))
fig.savefig(OUT / "figs" / "fig6_step0v3_pairing_by_sorter.png")
plt.close(fig)

# unpaired reasons, stacked per sorter (fractions of CS)
rs = pd.crosstab(cs["sorter"], cs["unpaired_reason"], normalize="index")
order = [c for c in ["paired", "duplicate_of_paired_CS", "near_miss:jitter_ns", "near_miss:symmetric_dip",
                     "near_miss:collision_dip", "near_miss:other", "weak_pause(0.7-0.85)", "no_pause_anywhere",
                     "no_SS_candidates"] if c in rs.columns]
fig, ax = plt.subplots(figsize=(10, 3.4))
y = np.arange(2)
left = np.zeros(2)
grays = plt.cm.Greys(np.linspace(0.25, 0.85, len(order)))
for k, c in enumerate(order):
    v = rs.reindex(["ks4", "ibl"])[c].fillna(0).values
    ax.barh(y, v, left=left, color=grays[k], edgecolor="white", linewidth=2, height=0.6)
    for yi, (l0, vv) in enumerate(zip(left, v)):
        if vv >= 0.06:
            ax.text(l0 + vv / 2, yi, f"{c}\n{vv:.0%}", ha="center", va="center", fontsize=6.5,
                    color="white" if k > len(order) / 2 else "black")
    left += v
n = cs.groupby("sorter").size()
ax.set_yticks(y, [f"{SORTER_NAME[s]} (n={n.get(s, 0)} CS)" for s in ("ks4", "ibl")])
ax.set_xlim(0, 1)
ax.set_xlabel("fraction of extra-good CS units")
ax.set_title("Pairing outcome of each CS (step0 v3; categories descriptive, not decisions)", fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "figs" / "fig7_unpaired_reasons.png")
plt.close(fig)
print(f"outputs in {OUT}")
