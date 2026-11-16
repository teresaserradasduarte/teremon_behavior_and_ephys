# -*- coding: utf-8 -*-
"""
check_v8_reproduces_v7.py — does step0 v8 with both extra CS sources OFF (env STEP0_V8_EXTRA_CS=0) reproduce v7?

Compares, table by table, only the columns v7 has (v8 adds columns), against the frozen v7batch output of the same
session. With the sources ON, use --c4 to compare only the rows of C4 CS (pairs / CS–CS / decomposition), which v8
must also leave unchanged. Read-only; prints a report.
Usage: python check_v8_reproduces_v7.py <session> <v8 folder name> [--c4]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CUR = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation")
sess, v8name = sys.argv[1], sys.argv[2]
C4_ONLY = "--c4" in sys.argv
A = CUR / sess / "step0_explore_bhvspan_v7batch" / "tables"
B = CUR / sess / v8name / "tables"
KEYS = {"cs_ss_pairs.csv": ["cs", "ss"], "cs_ss_pairs_lowfr_excluded.csv": ["cs", "ss"], "cs_cs_pairs.csv": ["cs_a", "cs_b"],
        "cs_cs_evidence.csv": ["cs_a", "cs_b"], "cs_event_decomposition.csv": ["cs", "unit"], "conflicts.csv": ["type", "unit"],
        "unit_table.csv": ["phyID"]}
ok_all = True
c4 = None
if C4_ONLY:
    u = pd.read_csv(B / "unit_table.csv")
    c4 = set(u.loc[u["cs_c4"].astype(bool), "phyID"])
for name, key in KEYS.items():
    fa, fb = A / name, B / name
    if fa.stat().st_size < 6 and fb.stat().st_size < 6:
        print(f"{name:34s} both empty"); continue
    a = pd.read_csv(fa) if fa.stat().st_size > 5 else pd.DataFrame()
    b = pd.read_csv(fb) if fb.stat().st_size > 5 else pd.DataFrame()
    if C4_ONLY and name not in ("unit_table.csv", "conflicts.csv", "cs_ss_pairs_lowfr_excluded.csv"):
        cols = [k for k in key if k in ("cs", "cs_a", "cs_b")]
        b = b[np.logical_and.reduce([b[k].isin(c4) for k in cols])] if len(b) else b
        a = a[np.logical_and.reduce([a[k].isin(c4) for k in cols])] if len(a) else a
    if name == "unit_table.csv" and C4_ONLY:
        cols_skip = {"role"}          # role legitimately changes for units that became CS
    else:
        cols_skip = set()
    common = [c for c in a.columns if c in b.columns and c not in cols_skip]
    missing = [c for c in a.columns if c not in b.columns]
    a2, b2 = a[common].set_index(key).sort_index(), b[common].set_index(key).sort_index()
    rows_only_a, rows_only_b = a2.index.difference(b2.index), b2.index.difference(a2.index)
    both = a2.index.intersection(b2.index)
    diffs = []
    for c in a2.columns:
        x, y = a2.loc[both, c], b2.loc[both, c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            bad = ~(np.isclose(x.astype(float), y.astype(float), rtol=1e-9, atol=1e-9, equal_nan=True))
        else:
            bad = ~((x.astype(str) == y.astype(str)) | (x.isna() & y.isna()))
        if bad.any():
            diffs.append(f"{c} ({int(bad.sum())})")
    ok = not (len(rows_only_a) or len(rows_only_b) or diffs or missing)
    ok_all &= ok
    print(f"{name:34s} {'IDENTICAL' if ok else 'DIFF'}: rows v7 {len(a2)} / v8 {len(b2)}, only-v7 {len(rows_only_a)}, "
          f"only-v8 {len(rows_only_b)}, columns compared {len(common)}" + (f", differing: {', '.join(diffs)}" if diffs else "")
          + (f", missing in v8: {missing}" if missing else ""))
fa, fb = json.loads((A / "cs_dup_families.json").read_text()), json.loads((B / "cs_dup_families.json").read_text())
print(f"{'cs_dup_families.json':34s} {'IDENTICAL' if fa == fb else 'DIFF'}: v7 {fa} | v8 {fb}")
ok_all &= fa == fb
print("RESULT:", "v8 reproduces v7 on every v7 column" if ok_all else "differences found (see above)")
