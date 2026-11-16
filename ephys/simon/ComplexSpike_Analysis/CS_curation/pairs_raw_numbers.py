# -*- coding: utf-8 -*-
"""
pairs_raw_numbers.py — raw counts of CS, paired CS, SS and CS→SS pairs per session (step0 v7, both FR floors),
BEFORE and AFTER collapsing potential duplicates. Counting only: nothing is merged in any data file.

Duplicates used for the "after" columns:
  * CS: step0's duplicate families (tables/cs_dup_families.json = connected components of the CS–CS pairs with
    dup_candidate = True in cs_cs_evidence.csv). A family counts as ONE cell; it is paired if any member is paired.
    `cs_dup_pairs_conflicting` = dup candidates whose evidence lines disagree (evidence_status "conflicting").
  * SS: ss_same_cell_check.py verdicts (raw layer only): SS pairs called "possibly the same cell" / "same cell" are
    collapsed (connected components). Sessions without the raw layer have no SS collapse (`raw_layer` = False).
  Pairs after = distinct (CS cell, SS cell) combinations among the candidate pairs.

Output: CS_curation/_summary/pairs_raw_numbers.csv  (one row per session × floor, plus per-animal and total rows)
"""
from pathlib import Path
import json

import pandas as pd

CUR = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation")
RUNS = {"20 Hz": "step0_explore_bhvspan_v7batch", "10 Hz": "step0_explore_bhvspan_v7batch_fr10"}


def components(edges, nodes):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in edges:
        parent.setdefault(a, a); parent.setdefault(b, b)
        parent[find(a)] = find(b)
    return {n: find(n) for n in parent}


def ss_same_cell_edges(sess):
    out = []
    for f in (CUR / sess / "ss_same_cell").glob("*/summary.json"):
        j = json.loads(f.read_text(encoding="utf-8"))
        v = str(j.get("verdict", ""))
        if v.startswith("possibly the same cell") or v.startswith("same cell"):
            out.append((int(j["ss_a"]), int(j["ss_b"])))
    return out


rows = []
for d in sorted(CUR.glob(f"*/{RUNS['20 Hz']}")):
    sess = d.parent.name
    raw = (CUR / sess / "raw_waveforms" / "_raw_batch_done.json").exists() or (CUR / sess / "ss_same_cell" / "ss_same_cell.js").exists()
    ss_edges = ss_same_cell_edges(sess) if raw else []
    for floor, run in RUNS.items():
        t = CUR / sess / run / "tables"
        if not t.exists():
            continue
        units = pd.read_csv(t / "unit_table.csv")
        cs_ids = units.loc[units.role == "CS", "phyID"].astype(int).tolist()
        p = pd.read_csv(t / "cs_ss_pairs.csv") if (t / "cs_ss_pairs.csv").stat().st_size > 5 else pd.DataFrame(columns=["cs", "ss", "candidate"])
        c = p[p.candidate.astype(bool)][["cs", "ss"]].astype(int)
        fam_f = t / "cs_dup_families.json"
        fams = json.loads(fam_f.read_text()) if fam_f.exists() and fam_f.stat().st_size > 2 else {}
        cs_edges = [(int(m[0]), int(x)) for m in fams.values() for x in m[1:]]
        ev_f = t / "cs_cs_evidence.csv"
        ev = pd.read_csv(ev_f) if ev_f.exists() and ev_f.stat().st_size > 5 else pd.DataFrame(columns=["dup_candidate", "evidence_status"])
        dup = ev[ev.dup_candidate.astype(bool)] if len(ev) else ev
        cs_cell = components(cs_edges, cs_ids)
        ss_cell = components(ss_edges, set(c.ss))
        c = c.assign(cs_cell=c.cs.map(cs_cell), ss_cell=c.ss.map(lambda s: ss_cell.get(s, s)))
        rows.append(dict(
            session=sess, animal=sess.rsplit("_", 1)[0], floor=floor, raw_layer=raw,
            n_cs=len(cs_ids), n_cs_paired=c.cs.nunique(), n_ss_paired=c.ss.nunique(), n_pairs=len(c),
            n_cs_multi_ss=int((c.groupby("cs").size() > 1).sum()),
            cs_dup_pairs=len(dup), cs_dup_pairs_converging=int(dup.evidence_status.isin(["converges_same_cell", "leans_same_cell"]).sum()) if len(dup) else 0,
            cs_dup_pairs_conflicting=int((dup.evidence_status == "conflicting").sum()) if len(dup) else 0,
            cs_dup_families=len(fams), cs_in_families=sum(len(m) for m in fams.values()),
            n_cs_after=len(set(cs_cell.values())), n_cs_paired_after=c.cs_cell.nunique(),
            ss_same_cell_pairs=len(ss_edges), n_ss_paired_after=c.ss_cell.nunique(),
            n_pairs_after=len(c[["cs_cell", "ss_cell"]].drop_duplicates()),
        ))
R = pd.DataFrame(rows)
num = [k for k in R.columns if k.startswith(("n_", "cs_", "ss_"))]
agg = []
for floor, g in R.groupby("floor", sort=False):
    for a, ga in g.groupby("animal"):
        agg.append(dict(session=f"{a} (all)", animal=a, floor=floor, raw_layer=f"{int(ga.raw_layer.sum())}/{len(ga)}", **ga[num].sum()))
    agg.append(dict(session="TOTAL", animal="all", floor=floor, raw_layer=f"{int(g.raw_layer.sum())}/{len(g)}", **g[num].sum()))
R = pd.concat([R, pd.DataFrame(agg)], ignore_index=True)
R["pct_cs_paired"] = (100 * R.n_cs_paired / R.n_cs).round(1)
R["pct_cs_paired_after"] = (100 * R.n_cs_paired_after / R.n_cs_after).round(1)
out = CUR / "_summary"; out.mkdir(exist_ok=True)
R.to_csv(out / "pairs_raw_numbers.csv", index=False)
show = ["session", "floor", "raw_layer", "n_cs", "n_cs_after", "n_cs_paired", "n_cs_paired_after", "n_ss_paired",
        "n_ss_paired_after", "n_pairs", "n_pairs_after", "cs_dup_pairs", "cs_dup_pairs_conflicting", "ss_same_cell_pairs"]
print(R[R.floor == "20 Hz"][show].to_string(index=False))
print(R[(R.floor == "10 Hz") & R.session.str.contains("all|TOTAL")][show].to_string(index=False))
