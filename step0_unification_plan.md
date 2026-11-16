# step0 unification plan (handoff for the chat that will write the final script)

Written 2026-11-05 by the "CS curation / browser" chat (session 7f49ba1e). This is a read-only review: no code was changed while preparing it.
The other chat is the "sorter comparison" chat (session f61de47e).

Code folder: `C:\Users\simon\Desktop\Code\EPhys\ComplexSpike_Analysis\CS_curation\`
Data folder: `C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation\`
The golden rule still applies: nothing under `D:\Learning Lab Dropbox\` is ever written.

## 1. Where things stand

**There is only one live `step0_explore_session.py`.** Both chats edited it in place, one after the other, and their edits don't collide textually. There is no git repo in the code folder, so this history comes from the two transcripts.

| time (11-05) | chat | change to `step0_explore_session.py` |
|---|---|---|
| 12:14 | sorter cmp | `STEP0_SESSION="<mouse>:<S#>"` env override: reads `_sorter_comparison/session_table.csv` and sets `RUN_TAG = "v2batch"` (section 1, lines ~53-65) |
| 12:49-13:13 | curation | **v3**: per-unit stable spans (`STABLE_SPAN`, `STAB_*`, section 3b), pair windows, `*_ref` reference-window columns, template overlay (8b), stability charts. `RUN_TAG = "v3"` |
| 13:14 | curation | `browser_template.html`: overlay, stability cards, pair-window columns, shading-order fix |
| 13:19 | sorter cmp | **frozen snapshot**: `_sorter_comparison/step0_v3_snapshot/step0_explore_session.py` = live v3 + 2 patches (`STEP0_RUN_TAG` env, empty-`dups` guard) |
| 13:21-14:46 | sorter cmp | `run_step0_batch.py` runs the snapshot on all sessions → `step0_explore_bhvspan_v3batch` (all 20 sessions with CS succeeded) |
| 14:47 | sorter cmp | overview x-label note: SS markers drawn +4 µm right |
| 16:11 | sorter cmp | **CS-spike sharing flag**: `CS_SHARE_FLAG = 0.2`, `cs_spike_sharing()`, self-test (d), columns `cs_share_frac/null/excess`, `ss_from_cs_frac`, `cs_share_flag` in `cs_ss_pairs`, printout. Live script now also reads `STEP0_RUN_TAG` |

**Live vs snapshot (diff checked):**
- The live file has the sharing flag, its self-test and the x-label note. The snapshot doesn't.
- The snapshot has the empty-`dups` guard. The live file doesn't.
- Otherwise they are identical.

**Method check:** R5 `step0_explore_bhvspan_v3` (my run) and R5 `step0_explore_bhvspan_v3batch` (snapshot run) are numerically identical:
- `cs_ss_pairs.csv` and `unit_table.csv` have the same shape and columns.
- The maximum absolute difference is 0.

So **the v3batch outputs are the current reference results for every session.**

**When this was written, two python processes were running:** the sorter chat's `STEP0_RUN_TAG=shareflag` test of the live script on Lindt S6 and Toblerone S3. They write to `…/step0_explore_bhvspan_shareflag`. Let them finish before touching the live file.

## 2. What the unified script must fix or decide

1. **Session selection and output naming are split, and one label is wrong.**
   - The default session is hard-coded (R5) with `RUN_TAG = "v3"`.
   - The env override sets `RUN_TAG = os.environ.get("STEP0_RUN_TAG", "v2batch")`. Without `STEP0_RUN_TAG`, v3 code writes into folders named `…bhvspan_v2batch`.
     - That already happened for Toblerone R4/R6, Milka R1/R5/R6/R7 and Ferrero R1/R4/R5.
     - Those folders are mislabelled, early v3 runs. Treat them as obsolete. Don't delete them; just ignore them, or ask Simon.
   - **Proposal:** one `RUN_TAG` constant (e.g. `"v4"`) that is overridable only by `STEP0_RUN_TAG`, and the same naming rule for single and batch runs: `<mouse>_<R#>/step0_explore_bhvspan_<RUN_TAG>`.
     - Keep `STEP0_SESSION` as the only way to pick a session; with no variable, the default is R5.
2. **Port the empty-`dups` guard from the snapshot** (section 6, `dups = pd.DataFrame(rows)`). It is needed for sessions with fewer than 2 CS.
3. **Fix `NameError: lc` in section 5.**
   - `LAGS_COARSE, LAGS_FINE = lc, lf` (live line ~1291) fails when `SS_IDS` is empty, because the loop never runs.
   - Fix: take the lags from the CCG constants, or guard with `if rows:`. Also make section 5 onwards tolerate an empty `pairs` table.
   - The batch never hit this only because 0-CS sessions are skipped. A session with CS but no SS candidates would crash.
4. **Show the CS-sharing flag in the HTML browser.** It currently exists only in the CSV and the printout.
   - Add `cs_share_excess`, `ss_from_cs_frac` and `cs_share_flag` to `pcols` (section 8b, line ~2148).
   - In `browser_template.html`, add:
     - a column or badge in the Pairing table;
     - a line in the pair detail, e.g. "x% of CS spikes also in this SS at zero lag (flag)".
   - It stays **flag only**: `candidate` is unchanged, as both chats agreed.
5. **Optional, ask Simon first:** the same zero-lag sharing measure on CS–CS pairs in 6b.
   - The 585/874 overlay already shows a fixed-lag spikelet-like component.
   - Sharing is close to what `coinc_frac_peak` already measures, so it may be redundant.
6. **Batch runner (`run_step0_batch.py`):**
   - Keep the frozen-snapshot practice (good: the method can't change mid-batch).
   - After unification: copy the unified script to a new snapshot folder (e.g. `_sorter_comparison/step0_v4_snapshot/`) and set `STEP0`, `RUN_TAG` and `OUT_SUBDIR` to the new tag.
   - The snapshot then needs **no hand patches**.
7. **Downstream scripts hard-code the input folder:**
   - `sorter_compare_step0v3.py` (`RUN_DIR`, line 33);
   - `sorter_compare_zero_lag_check.py` (`RUN_DIR`, line 36).
   - Point them at the new tag, preferably one shared constant or an env variable. Write their outputs to new folders (`step0v4`, …); don't overwrite `step0v3`.
8. **Update `CLAUDE.md`** (in `ComplexSpike_Analysis\`):
   - stable spans and pair windows, `*_ref` columns, template overlay;
   - the CS-sharing flag;
   - the batch runner and snapshot rule;
   - the output-folder naming, with a short table of which folder is current and which is obsolete.
9. **Method questions left open for Simon.** These are not part of the merge; keep current values:
   - `STAB_MAX_GAP = 1` (two bad bins stop a span) may be conservative. CS 622 stops at 2205 s, although later bins look stable again.
   - `CS_SHARE_FLAG = 0.2`.

## 3. Suggested order of work

1. Wait for the running `shareflag` tests. Then take the **live** file as the base: it already contains everything except the empty-`dups` guard.
2. Apply items 1-4 (+ 8). Bump `RUN_TAG` to `"v4"` so no earlier folder is touched.
3. Verify:
   - Self-tests (a)-(d) and the evidence self-test print PASSED.
   - R5 single run (no env var):
     - all v3 columns of `cs_ss_pairs.csv`, `unit_table.csv` and `cs_cs_evidence.csv` equal `step0_explore_bhvspan_v3` exactly, with only the new `cs_share_*` columns added;
     - candidates 24, families `{585: [585, 734, 874], 622: [622, 876]}`.
   - One session with fewer than 2 CS, and one with CS but no SS candidate (if one exists; otherwise a synthetic forced `SS_IDS = []` check), both run to the end.
   - `browser.html` opens with no console errors and shows the sharing flag.
4. Snapshot → batch (item 6) → downstream scripts (item 7).

## 4. Output folder inventory (data side, nothing deleted)

| folder | produced by | status |
|---|---|---|
| `1_CoteDor_R5/step0_explore`, `step0_explore_bhvwin` | v1 | old, keep |
| `1_CoteDor_R5/step0_explore_bhvwin_v2` | v2 single run | old, keep |
| `*/step0_explore_bhvwin_v2batch` | v2 batch | old, keep. Some sessions failed (`coinc_frac_peak` KeyError, `lc` NameError, Dropbox OSError) |
| `*/step0_explore_bhvspan_v2batch` | v3 code, wrong label | **obsolete / mislabelled**, ignore |
| `1_CoteDor_R5/step0_explore_bhvspan_v3` | v3 single run | current R5 reference (= v3batch) |
| `*/step0_explore_bhvspan_v3batch` | v3 snapshot batch | **current results, all sessions** |
| `2_Lindt_R6`, `3_Toblerone_R3` `/step0_explore_bhvspan_shareflag` | live script + sharing flag test | test only |
| `_sorter_comparison/step0v2`, `step0v3`, `zero_lag_check` | sorter-comparison scripts | step0v3 + zero_lag_check are current |
