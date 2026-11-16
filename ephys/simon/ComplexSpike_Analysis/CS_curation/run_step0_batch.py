# -*- coding: utf-8 -*-
"""
run_step0_batch.py — run step0_explore_session.py on every session of the sorter-comparison table.

Runs a FROZEN SNAPSHOT of step0 (STEP0 below), so edits to the live script cannot change the method mid-batch.
Each session runs in its own process (STEP0_SESSION="<mouse>:<S#>", STEP0_RUN_TAG=RUN_TAG, Agg backend) and writes to
CS_curation/<mouse>_<R#>/<OUT_SUBDIR>/. Sessions whose tables/cs_ss_pairs.csv already exists are skipped (FORCE to
redo); sessions with no extra-good CB PkC_cs unit (from pass-1 units_all.csv) are skipped (nothing to pair).
Failed sessions are retried once (Dropbox hydration of online-only files can fail under parallel load).
Logs: CS_curation/_sorter_comparison/step0_batch_logs_<RUN_TAG>/.

Snapshot rule: after any change to step0_explore_session.py, copy it UNCHANGED to a new folder
_sorter_comparison/step0_<version>_snapshot/ and point STEP0 / RUN_TAG here at it (no hand patches).
RUN_TAG = env STEP0_BATCH_TAG (default "v7batch"; v6batch = step0_v6_snapshot, v5batch = step0_v5_snapshot); the downstream scripts sorter_compare_step0v3.py and
sorter_compare_zero_lag_check.py read the same variable. History: v3batch = step0_v3_snapshot (hand-patched). Old snapshots (v3, v4) are in _sorter_comparison/_archive/.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import os
import subprocess
import sys
import time

import pandas as pd

LOCAL = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation")
STEP0 = LOCAL / "_sorter_comparison" / "step0_v7_snapshot" / "step0_explore_session.py"
RUN_TAG = os.environ.get("STEP0_BATCH_TAG", "v7batch")
# tier-1 FR floor toggle: env STEP0_SS1_FLOOR = "20" (default) or "10"; passed to every session process (inherited
# environment). step0 appends "_fr10" to its output folder for the 10 Hz floor, so the batch does the same here.
SS1_FLOOR = os.environ.get("STEP0_SS1_FLOOR", "20")
assert SS1_FLOOR in ("20", "10"), SS1_FLOOR
_FLOOR_SUFFIX = "" if SS1_FLOOR == "20" else f"_fr{SS1_FLOOR}"
OUT_SUBDIR = f"step0_explore_bhvspan_{RUN_TAG}{_FLOOR_SUFFIX}"      # step0 v7 (pause 3–8 ms, tier-1 FR floor, event decomposition); same rule as step0 itself
assert STEP0.exists(), f"snapshot missing: {STEP0}"
N_WORKERS = 3
FORCE = False
FIRST = [("5_FerreroRocher", "S6"), ("3_Toblerone", "S5")]   # inspection sessions first

tab = pd.read_csv(LOCAL / "_sorter_comparison" / "session_table.csv")
u1 = pd.read_csv(LOCAL / "_sorter_comparison" / "units_all.csv")
n_cs = u1[u1["c4_label"] == "PkC_cs"].groupby(["mouse", "S"]).size()
LOG_DIR = LOCAL / "_sorter_comparison" / f"step0_batch_logs_{RUN_TAG}{_FLOOR_SUFFIX}"
LOG_DIR.mkdir(parents=True, exist_ok=True)
order = FIRST + [(m, s) for m, s in zip(tab["mouse"], tab["S"]) if (m, s) not in FIRST]


def run(ms):
    m, s = ms
    if n_cs.get((m, s), 0) == 0:
        return f"skip   {m} {s} (0 CS)"
    r = tab.loc[(tab["mouse"] == m) & (tab["S"] == s), "R"].iloc[0]
    done = LOCAL / f"{m}_{r}" / OUT_SUBDIR / "tables" / "cs_ss_pairs.csv"
    if done.exists() and not FORCE:
        return f"skip   {m} {s} (done)"
    env = dict(os.environ, STEP0_SESSION=f"{m}:{s}", STEP0_RUN_TAG=RUN_TAG, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
    t0 = time.time()
    for attempt in (1, 2):
        with open(LOG_DIR / f"{m}_{s}.log", "w", encoding="utf-8") as f:
            rc = subprocess.run([sys.executable, "-u", str(STEP0)], env=env, stdout=f,
                                stderr=subprocess.STDOUT).returncode
        if rc == 0 and done.exists():
            break
    ok = rc == 0 and done.exists()
    return f"{'ok   ' if ok else 'FAIL ' + str(rc)}  {m} {s}  {(time.time() - t0) / 60:.1f} min (attempts {attempt})"


with ThreadPoolExecutor(N_WORKERS) as ex:
    for msg in ex.map(run, order):
        print(msg, flush=True)
