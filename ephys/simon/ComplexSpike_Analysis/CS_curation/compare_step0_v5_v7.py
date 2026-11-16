# -*- coding: utf-8 -*-
"""
compare_step0_v5_v7.py — what did step0 v7 change, session by session?  (read-only on the step0 outputs)

v7 = v5 + (1) pause measured after the CS event (3–8 ms instead of 0–8 ms; v6 used 4–12 ms), (2) tier-1 FR floor (toggle STEP0_SS1_FLOOR: 20 or 10 Hz),
(3) CS event decomposition table + cs_event flag (information only). (1) and (2) can change pairings; this script
separates them:
  - "floor"  : a v5 candidate SS that v7 no longer tests because it is an extra-good unit below the floor (SS_lowfr)
  - "window" : every other change of candidate status (same SS tested in both runs, different window)
Inputs : CS_curation/<mouse>_<R#>/step0_explore_bhvspan_{v5batch,v7batch}/tables/
Output : CS_curation/_sorter_comparison/step0_v7_vs_v5/  (tables + figures + README.md)
"""
# %% 1. PATHS + LOAD
from pathlib import Path
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use(os.environ.get("MPLBACKEND", "Agg"))
import matplotlib.pyplot as plt

LOCAL = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation")
# tier-1 FR floor toggle (same env as step0 / run_step0_batch): "20" (default) or "10" → reads the *_fr10 run
SS1_FLOOR = os.environ.get("STEP0_SS1_FLOOR", "20")
_FS = "" if SS1_FLOOR == "20" else f"_fr{SS1_FLOOR}"
OLD, NEW = "step0_explore_bhvspan_v5batch", os.environ.get("STEP0_V7_DIR", f"step0_explore_bhvspan_v7batch{_FS}")
OUT = LOCAL / "_sorter_comparison" / f"step0_v7{_FS}_vs_v5"
OUT.mkdir(parents=True, exist_ok=True)
tab = pd.read_csv(LOCAL / "_sorter_comparison" / "session_table.csv")

rows, pair_rows, dec_all, low_all = [], [], [], []
for _, s in tab.iterrows():
    sess = LOCAL / f"{s['mouse']}_{s['R']}"
    a, b = sess / OLD / "tables", sess / NEW / "tables"
    if not (a / "cs_ss_pairs.csv").exists() or not (b / "cs_ss_pairs.csv").exists():
        continue
    p5, p6 = pd.read_csv(a / "cs_ss_pairs.csv"), pd.read_csv(b / "cs_ss_pairs.csv")
    u6 = pd.read_csv(b / "unit_table.csv").set_index("phyID")
    lowfr = set(u6.index[u6["role"] == "SS_lowfr"])
    cs_ids = sorted(u6.index[u6["role"] == "CS"])
    key = dict(mouse=s["mouse"], R=s["R"], S=s["S"], sorter=s.get("sorter_type", s.get("sorter", "")))
    m = p5[["cs", "ss", "candidate", "pause_ratio", "rank"]].merge(
        p6[["cs", "ss", "candidate", "pause_ratio", "pause_ratio_0_8", "rank", "cs_event_excess", "pause_dur50_ms",
            "pause_min_lag_ms"]], on=["cs", "ss"], how="outer", suffixes=("_v5", "_v7"))
    m["candidate_v5"] = m["candidate_v5"].fillna(False).astype(bool)
    m["candidate_v7"] = m["candidate_v7"].fillna(False).astype(bool)
    m["change"] = np.select([m["candidate_v5"] & ~m["candidate_v7"] & m["ss"].isin(lowfr),
                             m["candidate_v5"] & ~m["candidate_v7"], ~m["candidate_v5"] & m["candidate_v7"]],
                            ["lost (floor)", "lost (window)", "gained (window)"], "")
    pair_rows.append(m[m["change"] != ""].assign(**key))
    m5 = set(p5.loc[p5["candidate"], "cs"]); m6 = set(p6.loc[p6["candidate"], "cs"])
    r5 = p5[p5["rank"] == 1].set_index("cs")["ss"]; r6 = p6[p6["rank"] == 1].set_index("cs")["ss"]
    both = sorted(m5 & m6)
    rows.append(dict(key, n_cs=len(cs_ids), matched_v5=len(m5), matched_v7=len(m6),
                     lost_cs=len(m5 - m6), gained_cs=len(m6 - m5),
                     rank1_changed=int(sum(r5.get(c) != r6.get(c) for c in both)),
                     lost_pairs_floor=int((m["change"] == "lost (floor)").sum()),
                     lost_pairs_window=int((m["change"] == "lost (window)").sum()),
                     gained_pairs_window=int((m["change"] == "gained (window)").sum()),
                     n_lowfr_units=len(lowfr)))
    if (b / "cs_event_decomposition.csv").exists():
        d = pd.read_csv(b / "cs_event_decomposition.csv").assign(**key)
        if "strong" not in d:                       # early v7 test output: same rule as step0 v7
            d["strong"] = (d["excess"] > 0.2) | d["mostly_cs_pieces"]
        dec_all.append(d)
    if (b / "cs_ss_pairs_lowfr_excluded.csv").exists() and (b / "cs_ss_pairs_lowfr_excluded.csv").stat().st_size > 5:
        low_all.append(pd.read_csv(b / "cs_ss_pairs_lowfr_excluded.csv").assign(**key))

S = pd.DataFrame(rows)
P = pd.concat(pair_rows, ignore_index=True) if pair_rows else pd.DataFrame()
DEC = pd.concat(dec_all, ignore_index=True) if dec_all else pd.DataFrame()
LOW = pd.concat(low_all, ignore_index=True) if low_all else pd.DataFrame()
print(f"{len(S)} sessions with both runs")

# %% 2. TABLES
S.to_csv(OUT / "per_session.csv", index=False)
P.to_csv(OUT / "changed_pairs.csv", index=False)
if len(LOW):
    LOW.to_csv(OUT / "lowfr_excluded_pairs_all.csv", index=False)
if len(DEC):
    DEC.to_csv(OUT / "cs_event_decomposition_all.csv", index=False)
    dsum = DEC.groupby(["mouse", "R", "sorter"]).agg(n_pieces=("unit", "size"), n_strong=("strong", "sum"),
                                                     n_zero_lag=("piece", lambda x: (x == "zero-lag copy").sum()),
                                                     n_ss_cand_in_event=("is_ss_candidate", "sum")).reset_index()
    dsum.to_csv(OUT / "decomposition_per_session.csv", index=False)

# %% 3. FIGURES
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
x = np.arange(len(S)); lab = [f"{m.split('_')[1][:6]} {r}" for m, r in zip(S["mouse"], S["R"])]
ax[0].bar(x - 0.2, S["matched_v5"] / S["n_cs"], 0.4, label="v5", color="#bbbbbb")
ax[0].bar(x + 0.2, S["matched_v7"] / S["n_cs"], 0.4, label="v7", color="#9a7bc1")
ax[0].set_xticks(x, lab, rotation=70, fontsize=7); ax[0].set_ylabel("fraction of CS with an SS candidate")
ax[0].legend(frameon=False); ax[0].set_title("Matched CS per session")
for k, (c, col) in enumerate([("lost_pairs_floor", "#d1495b"), ("lost_pairs_window", "#e8a33d"),
                              ("gained_pairs_window", "#3a8a5c")]):
    ax[1].bar(x + (k - 1) * 0.27, S[c] * (1 if "gained" in c else -1), 0.27, color=col,
              label=c.replace("_pairs_", " pairs, ").replace("_", " "))
ax[1].axhline(0, color="#333", lw=0.8); ax[1].set_xticks(x, lab, rotation=70, fontsize=7)
ax[1].set_ylabel("CS×SS candidate pairs (lost < 0 < gained)"); ax[1].legend(frameon=False, fontsize=8)
ax[1].set_title("Which v7 change moved the candidates")
if len(P):
    q = P[P["change"] != "lost (floor)"]
    col = q["change"].map({"lost (window)": "#e8a33d", "gained (window)": "#3a8a5c"})
    ax[2].scatter(q["pause_ratio_v5"], q["pause_ratio_v7"], c=col, s=22, edgecolor="none")
    ax[2].plot([0, 1.5], [0, 1.5], color="#999", lw=0.8); ax[2].axhline(0.7, color="#bbb", ls=":"); ax[2].axvline(0.7, color="#bbb", ls=":")
    ax[2].set_xlabel("pause ratio v5 (0–8 ms)"); ax[2].set_ylabel("pause ratio v7 (3–8 ms)")
    ax[2].set_title("Pairs whose candidate status changed by the window")
fig.tight_layout(); fig.savefig(OUT / "fig1_v5_vs_v7_per_session.png", dpi=150); plt.close(fig)

if len(P):
    q = P[P["change"] == "lost (window)"]
    if len(q):
        fig, ax = plt.subplots(figsize=(5.5, 4.2))
        ax.hist(q["pause_dur50_ms"].dropna(), bins=np.arange(0, 31, 1), color="#e8a33d")
        ax.set_xlabel("pause 50% recovery (ms, v7 CCG)"); ax.set_ylabel("pairs lost by the window")
        ax.set_title("Pairs lost by the 3–8 ms window: short pauses?")
        fig.tight_layout(); fig.savefig(OUT / "fig2_lost_by_window_pause_duration.png", dpi=150); plt.close(fig)

if len(DEC):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    st = DEC[DEC["strong"]]
    ax[0].hist(st["lag_med_ms"], bins=np.arange(-0.5, 4.01, 0.1), color="#9a7bc1")
    ax[0].set_xlabel("median lag after the CS spike (ms)"); ax[0].set_ylabel("strong pieces (all sessions)")
    ax[0].set_title("When do the strong CS pieces occur?")
    ax[1].scatter(DEC["excess"], DEC["unit_locked_frac"], s=8, c=np.where(DEC["strong"], "#9a7bc1", "#bbbbbb"))
    ax[1].set_xlabel("excess (share of CS events with a spike, minus chance)")
    ax[1].set_ylabel("share of the cluster's own spikes that are CS-locked"); ax[1].set_title("All CS × cluster pieces")
    fig.tight_layout(); fig.savefig(OUT / "fig3_event_decomposition_all_sessions.png", dpi=150); plt.close(fig)

# %% 4. README
tot = S[["n_cs", "matched_v5", "matched_v7", "lost_cs", "gained_cs", "rank1_changed", "lost_pairs_floor",
         "lost_pairs_window", "gained_pairs_window"]].sum()
lines = [f"# step0 v7 vs v5 ({len(S)} sessions; tier-1 FR floor {SS1_FLOOR} Hz)", "",
         f"- Matched CS (>= 1 SS candidate): v5 {tot.matched_v5} → v7 {tot.matched_v7} of {tot.n_cs} CS "
         f"(lost {tot.lost_cs}, gained {tot.gained_cs}); rank-1 SS changed for {tot.rank1_changed} CS matched in both.",
         f"- Candidate pairs lost by the tier-1 FR floor: {tot.lost_pairs_floor}; lost by the 3–8 ms window: "
         f"{tot.lost_pairs_window}; gained by the window: {tot.gained_pairs_window}.",
         ]
if len(LOW):
    w = LOW[LOW["would_be_candidate"]]
    lines.append(f"- Pairs with a sub-{SS1_FLOOR} Hz extra-good unit that WOULD be candidates in v7 without the floor: {len(w)} "
                 f"({w[['mouse', 'R']].drop_duplicates().shape[0]} sessions; table lowfr_excluded_pairs_all.csv).")
if len(DEC):
    lines.append(f"- Event decomposition: {len(DEC)} CS × cluster pieces, {int(DEC['strong'].sum())} strong; "
                 f"SS candidates carrying pieces of their CS's event: {int(DEC['is_ss_candidate'].sum())}.")
if len(P):
    q = P[P["change"] == "lost (window)"]
    if len(q):
        lines.append(f"- Pairs lost by the window: median 50% recovery {q['pause_dur50_ms'].median():.0f} ms "
                     f"median v7 pause ratio "
                     f"{q['pause_ratio_v7'].median():.2f}.")
lines += ["", "Files: per_session.csv, changed_pairs.csv, lowfr_excluded_pairs_all.csv, cs_event_decomposition_all.csv,",
          "decomposition_per_session.csv, fig1–fig3. Script: Code/EPhys/ComplexSpike_Analysis/CS_curation/compare_step0_v5_v7.py"]
(OUT / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(lines))
with pd.option_context("display.width", 200):
    print(S.to_string(index=False))
