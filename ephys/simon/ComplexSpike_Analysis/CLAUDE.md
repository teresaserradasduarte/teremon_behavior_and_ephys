# CLAUDE.md — Cerebellar complex-spike analysis (push/pull lever + reach task)

Context file for Claude sessions working on this project. Read it fully before touching code or data.
It records what is known, how the data are produced, the constraints on working with them, and the
questions that are still open for Simon. **Do not start an analysis plan from this file alone.**
The open questions in §9 come first.

## 0. GOLDEN RULE — DROPBOX IS READ-ONLY. NEVER, EVER MODIFY IT.

**Never write, modify, move, rename, overwrite or delete anything under
`D:\Learning Lab Dropbox\` — not a data file, not a label `.tsv`, not a folder, nothing. Ever.**
- This applies to every script, every shell command and every Claude action, with no exceptions.
- Dropbox is a **read-only source**. All outputs go to local folders (e.g. `C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation\…`).
- Scripts must guard against it: refuse any output path under the Dropbox root, and open data files read-only (`np.load(..., mmap_mode='r')`, `pd.read_csv`).
- Reading a cloud-only file makes Dropbox download it; that is allowed (it is a read).

> Sister project: the temperature / Fixed Interval reach pipeline (`ReachAnalysis`, with its own CLAUDE.md).
> That is a **different task and dataset**. Do not reuse its assumptions (block structure, trial
> definition, reach classes, px→mm factors) without checking.

> **Source of §3–§7.** These sections come from reading the lab repo
> `teresaserradasduarte/teremon_behavior_and_ephys` (GitHub, public; commit `8d400de`, 2026-11-03;
> local copy `C:\Users\simon\Documents\GitHub\teremon_behavior_and_ephys`). The dataset
> ("Teremon", ChocolateGroup) was collected and preprocessed by Teresa Serradas Duarte, mostly
> in MATLAB. Everything here was **inferred from code, not checked on the data**. Treat it as
> "expected format" until a script confirms it on one session.

---

## 1. Scientific goal

Understand **when and where (in kinematic space)** cerebellar **complex spikes (CS)** occur during a
reach, and **what they encode**. The task mixes two problems:

- a **discrete** problem: choosing push vs. pull under an uncued, block-switching rule (trial and error);
- a **continuous** problem: reaching to a droplet whose location (L / C / R) also changes in blocks.

**Step 0, before any behavioural analysis, is dataset curation:**
1. Count the CS units per session and per animal.
2. Estimate what proportion of CS clusters are **duplicates** of one another (Kilosort can duplicate clusters).
3. Count how many CS units have a **matching simple-spike (SS) unit**, identified by an SS pause in the
   CS→SS cross-correlogram.

Only the curated set goes downstream: one row per identified Purkinje cell, with or without a paired SS.

## 2. Task design

- **Uncued two-forced-choice, self-initiated.** When a trial is available, the push or pull is both the trial initiation and the choice.
- **Correct choice:** a water droplet appears at L / C / R and the mouse reaches for it.
- **Wrong choice:** buzz + light, then a timeout.
- **Rule blocks (push vs. pull)** switch without a cue; the MI project's notes say blocks are 20–30 trials. Droplet-location blocks switch **independently** of the rule blocks.
- **Effectors:** both paws push or pull the lever; **one dominant paw reaches**.
  - Dominant paw is **R** for 1_CoteDor, 4_Milka and 5_FerreroRocher.
  - Dominant paw is **L** for 2_Lindt and 3_Toblerone.
  - Source: `s01_NeuronStructAndPlot_v5.m`; `bhv.paw_pref` holds it per session.
- **Lever = load cell.** The raw stream is `behavior.inputs.paws_all`: (n, 2) = Harp timestamp + force in mV/V, at about 1 kHz. This load cell defines push and pull.
- Setup name in paths: `headfixed_dynamicTarget`.

## 3. Dataset and folder layout

### Animals and sessions
- Group `20230801_ChocolateGroup`: 5 mice, `1_CoteDor`, `2_Lindt`, `3_Toblerone`, `4_Milka`, `5_FerreroRocher`.
- Up to 7 sessions per mouse.
  - Ephys folders use **`S#`**, e.g. `18082023_Milka_StrCer_S4_g0`.
  - Behaviour and derived `.mat` folders use **`R#`** for the same session (S4 ↔ R4).
- 27 cerebellar sessions are listed in `paths_ephysTeremon_windows(_mat).m/.py`.
- Milka S4 is the "reference session" used for development.

### Roots (Windows)
- Lab root: `D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\`. The alternatives in the code are `Y:\shared-paton\teremon\ephys_curated\` and `F:\curated_ephys_data\`.
- **Working copy (decided): `C:\Users\simon\Desktop\Data\ephys_and_behavior\`.** Simon places the needed
  files there. The Dropbox tree (`D:\Learning Lab Dropbox\…`) is the original and is always available as a fallback.
  Scripts should take both roots from one config and prefer the local copy.
- `C:\Users\SimonZ\Documents\Data\EPhys\DTATA_2CURATE` appears in old code and is not used.

### Folder structure under the lab root

```
ephys_data\20230801_ChocolateGroup\<n_Mouse>\<ddmmyyyy_Mouse_StrCer_S#_g0>\<..._g0_imecN>\
    *.ap.bin / *.ap.meta / *.lf.bin / *.lf.meta        (SpikeGLX raw — never load whole)
    catGT\kilosort4\            ← sorter folder, OR
    ibl_sorter_results[_driftAdaptNN]\
behavior_data\raw_data\20230511_ChocolateGroup\headfixed_dynamicTarget\<mouse>\R#\behavior_session.mat
behavior_data\analyzed_data\mat_files\20230511_ChocolateGroup\headfixed_dynamicTarget\<mouse>\R#\session_reaching_data_paw.mat
ephys_and_behavior\mat_files\20230801_ChocolateGroup\<mouse>\R#\  behavior_fundamentals.mat, eg_neurons_CB.mat, eg_neurons_BG.mat, cell_type_CB.mat
ephys_and_behavior\mat_files\20230801_ChocolateGroup\  st_CB_DCN.mat, st_eg_neurons.mat, eg_neurons.mat, eg_neurons_CB_cellType.mat, teremon_paths.mat
ephys_and_behavior\out_files\20230801_ChocolateGroup\...   (figures)
histology\20230801_ChocolateGroup\<mouse>\LineFit\CB\outputForIBL\<session_imecN>\  (channel_locations)
```

The behaviour raw tree uses a **different group name** (`20230511_…`) from the ephys tree (`20230801_…`). The mouse and session names underneath are the same.

## 4. Neural data: acquisition, sorting, labels

### Acquisition
- **Neuropixels 1.0.** SpikeGLX writes 385 channels per file, and channel 385 is the sync channel.
- **Two probes per session**, `imec0` and `imec1`, one in CB and one in BG, each **sorted separately**.
- Which imec index is CB depends on the mouse:
  - R-pawed mice (1, 4, 5): **imec0 = CB**.
  - L-pawed mice (2, 3): **imec1 = CB**.
  - Toblerone S2 follows the rule: **CB = imec1** (confirmed by Simon; the `.py` list is wrong, see §8.3).
- Some early sessions are single-region, e.g. `CoteDorS1CB` and `CoteDor_Cer_S5`.
- Sampling rate is about 30 kHz. Some overview scripts hard-code `30000`. **Always read `sample_rate` from `params.py`**, or `imSampRate` from `.ap.meta`.

### Sorting: two pipelines are mixed across sessions
- **CatGT → Kilosort 4** (`catGT\kilosort4`). Extra files: `ops.npy` (contains `dshift` drift), `spike_positions.npy`, `cluster_ContamPct.tsv`.
- **IBL sorter** (`ibl_sorter_results`, some with `_driftAdapt2/20/40`). Extra files: `drift.times.npy`, `drift.um.npy`, `drift_depths.um.npy`, `pc_features.npy`, `pc_feature_ind.npy`.
- Both are phy-compatible. Shared files: `spike_times.npy` (samples), `spike_clusters.npy`, `spike_templates.npy`, `templates.npy` (whitened), `whitening_mat_inv.npy`, `amplitudes.npy`, `channel_map.npy`, `channel_positions.npy`, `params.py`, `cluster_info.tsv`.
- Curation code must not assume one sorter. Record which sorter produced each session in the curation output: sorter identity affects duplicate rates and collision dips.

### Phy manual-curation label files (in the sorter folder)
- `cluster_group.tsv`: `noise` / `mua` / `good`. The MATLAB loader codes these as 0 = noise, 1 = mua, 2 = good.
- `cluster_Extragood.tsv` (sometimes lowercase): column `Extragood`, value `y` marks **extra-good** units. Everything downstream in MATLAB uses **extra-good units only**.
- `cluster_CS.tsv` **or** `cluster_potential_CS.tsv`: value `y` marks a **manually flagged potential CS** cluster. This answers the old question "were CS clusters labelled manually?". It does exist, at least for some sessions. **Not yet known:** whether every session has it, and how it relates to the C4 `PkC_cs` labels.

### C4 cell-type output (in the sorter folder)
- `cluster_predicted_cell_type.tsv`, `cluster_confidence_ratio.tsv`, `cluster_pred_probability.tsv`, and `cluster_model_votes.tsv` (votes stored as the string `"k/n"`). All are keyed by `cluster_id`.
- Labels seen in the code: `PkC_cs`, `PkC_ss`, `MLI`, `MFB` (C4 also has `GoC`).
- Lab post-filter (`s04_CB_cellType.m`):
  - Only `MFB` with vote ratio < 0.94 is set to NaN.
  - **No extra threshold is applied to `PkC_cs` / `PkC_ss`.**
  - In DCN, valid units need vote ratio ≥ 0.94 and pred_prob > 0.6. They are relabelled `DCN_MF` / `DCN_PCproj`. A `PkC_cs` unit found in DCN triggers a "revisit tracking" warning.
- C4 does not pair CS with SS and does not detect duplicates (see §6).

### Brain region per cluster
- Histology was aligned with the IBL ephys-alignment workflow and gives `channel_locations` (json/npy). `histology/simplify_region_map.py` reduces it to a TSV with columns `channel, brain_region, …`.
- `histology/mapping_channel_to_regions.py` writes **`cluster_channel_region_map.txt`** in the sorter folder, with columns `cluster_id, channel, brain_region`. Each cluster gets the channel closest to its **median spike depth**.
- CB-probe regions include `CB` (cortex) and `DCN`; the code also mentions `BS`, `MB`, `P`. The **current scope is `brain_region == "CB"`.**

## 5. MATLAB pipeline that produced the derived files

1. **`loadKSdir_TD.m`** loads a sorter folder.
   - `st = spike_times / sample_rate`, in seconds on the **ephys (imec) clock**.
   - Removes noise clusters and un-whitens templates.
   - Reads the Extragood and potential-CS files.
   - Peak channel = largest peak-to-peak of the un-whitened template.
2. **`s01_NeuronStructAndPlot_v5.m`** writes `<sorter folder>\neurons_session.mat` (`neurons`, `neurons_params`). It runs once per probe-session.
   - **Sync:**
     - Reads the SpikeGLX sync channel from the **LF** `.bin` (`extractSyncChannel`, ch 385) and parses digital lines (`spikeGLXdigitalParse`); bit 10 is the per-trial TTL (`eventTimes_water_giv`).
     - Fits a line from these times to `behavior.sync.time_newTrial_log` (the Harp log).
     - The resulting function is misleadingly named `tm_bhv2ephys`, but it maps **ephys → Harp/behaviour time**.
     - Manual flags handle crashed sessions (`cut_sess_sync_flag`) and a dropped first behaviour trial (`flag_remove_1st_behavior_frame`, saved back into `behavior_session.mat`).
     - The fit coefficients are **not saved explicitly**. They can be recovered from the paired fields `st_raw_all` (ephys clock) and `st_all` (Harp clock).
   - **Per unit:**
     - `st_raw_all` / `st_all`: all spikes on the ephys / Harp clock.
     - `st_raw` / `st`: only spikes inside the behaviour epoch.
     - `phyID`, `quality` (1 mua / 2 good), `extraGood`, `potentialCS`.
     - Templates (IDs, weights, waveforms, peak channel), computed **from behaviour-epoch spikes**.
     - `amplitudes`, `depth` (median spike depth).
     - `CCG`: autocorrelogram, 1 ms bins, ±150 ms, `CCGBz`, rate-normalised.
     - ISI histograms, using 0–5 s bins when `potentialCS == 1`.
     - `st_init{trial}` and `st_reach{trial}`.
3. **`s02_EGneurons_FR_all.m`** writes `…\mat_files\…\<mouse>\R#\eg_neurons_CB.mat` (`eg_neurons`, `figProp`, `eg_neu_FR_params`).
   - **Extra-good units only.**
   - `st` = **`st_all`** (the whole recording on the Harp clock, *not* restricted to the behaviour epoch).
   - `st_reach{k}` = spikes relative to `bhv.reach_times_inVec(k)`, in a [-3, 3] s window plus kernel padding; `st_init{k}` is the same relative to `bhv.pp_times(k)`.
   - `FR_reach` / `FR_init`: instantaneous 1/ISI rate at 2 ms bins, then a Gaussian with σ = 20 ms.
   - **`potentialCS` is not copied over.**
4. **`s03_allEGneurons_perRegTracked.m`** adds `brain_region` from `cluster_channel_region_map.txt`.
   **"Tracked" means the probe track reconstructed from histology** (IBL alignment).
   - It replaces the older `s03_allEGneurons_perReg.m`, which assigned regions from manual depth cut-offs (`regions separation - <mouse>.csv`).
   - It has **nothing to do with tracking cells across sessions.**
   - It also records `is_nreg_matching`, the number of units with no region match per mouse/session/probe. It also writes the group-level `eg_neurons.mat`, which includes the count arrays `n_eg_neu`, `n_good_neu`, `n_mua_neu`, sized [mouse, session, BG/CB].
5. **`s04_CB_cellType.m`** joins the C4 tables onto `eg_neurons` by `phyID` and adds `cell_type`. It saves `cell_type_CB.mat` per session and `eg_neurons_CB_cellType.mat` for the group.
6. **`generateRegStruct_selectedFields_CB.m`** writes `st_CB_DCN.mat` (`all_CB_sts`, `all_DCN_sts`). **`s05_allEGneurons_sts.m`** writes `st_eg_neurons.mat` for all regions.
   - One struct entry per neuron, with 16 fields: `phyID, st_init, st_reach, brain_region, cell_type, idx_init_vPP_invPP, DCnD_init, idx_reach_cat, idx_reach_hit, idx_reach_succ, idx_reach_hDCnD, idx_reach_LCR, idx_reach_PP, reach_px, mouse, sess`.
   - Behaviour fields are **copied identically into every neuron of a session**.
   - uid convention in the Python projects: `f"{REG}|{mouse}|{sess}|{phyID}"`.
   - Files are MATLAB **v7.3** (HDF5). Simon's earlier scripts read them with `h5py` plus a recursive reference resolver. Strings come back as uint16 char codes and must be decoded.

**Existing CS/SS code in the lab repo**
- `ephys/02_single_neuron_analysis/cerebellum_neurons/SS_and_CCs_pair_histogram.m` plots the CS×SS CCG for **one hand-picked pair** (Milka R4, `neurons` indices 86 and 89; 1 ms bins, ±150 ms). There are no criteria and no batch processing.
- `aux_fig_CS_v2_zoom_Pcs_v2.m` is a figure script.
- No duplicate detection exists anywhere in the repo.

**Simon's earlier CS code** (`ephys_and_behavior/Simon_CS_analysis/`)
- `CS_analysis_workspace.py` loads `st_CB_DCN.mat` and selects `cell_type == "PkC_cs"`. It caches to `filtered_data.pkl`, which is pickled (the parquet attempt failed). It plots rasters/PSTHs by reach category, success and L/C/R, and detects onsets.
- `Whole_Session_Analysis_5_R5.py` works on Côte d'Or R5.
  - Inputs: `session_reaching_data_paw.mat` (`reaches`, `session`, `trials`), session-long paw traces `ppDom_paw_sess` / `ppNonDom_paw_sess` (3 × n_frames, **px**), and `eg_neurons_CB.mat` (`st`, `cell_type`).
  - It reconstructs a continuous frame-timestamp vector from `trials.sync.frames_timestamps` and `trials.size_trial`.
  - It segments reaches itself (peak/onset, speed and PC1 arbitration) and does grooming / small-movement detection.
- These are exploratory scripts, not a pipeline. Reuse the loaders, not the assumptions.

## 6. Behaviour and kinematics

- **Camera: 400 fps** (2.5 ms per frame), DLC-tracked 3D paw, with one video per trial.
  - Frame timestamps come from a per-trial **linear interpolation** between consecutive Harp sync times (`step1_find_reaches…m`). This assumes no dropped frames within a trial.
    Exactly: frames 1…siz(tt) of video trial tt are mapped onto [sync(tt), sync(tt+1)], with
    `sync = [behavior.sync.time_newTrial_log; real_log_all(end,1)]` → stored as `trials.sync.frames_timestamps` (n_trials × max siz).
    Consequences: implied frame period = Δsync/(siz−1) (≈399.9 fps, not 400); the last frame of trial tt and frame 1 of
    tt+1 get the **same** timestamp (drop one when stitching); any video that does not span TTL→TTL is stretched/shifted.
  - Frame times are on the Harp clock (see §14.2 for the full clock chain).
- **PX2MM = (0.1, 0.07, 0.1)** for (x, y, z). **y = medio-lateral.** Source: `config.py`, verified there against `bhv.reaches_inVec_mm`. Whether these are the same values as the cooling rig is still open.
- **`behavior_fundamentals.mat` (`bhv`)** is written by `behavior/paw_and_body_trajectories/step3_bhv_struct.m`. Key fields:
  - `pp_times`, the push/pull event times. It is **concatenated as [valid initiations; invalid pushes; invalid pulls]. It is NOT chronological.**
  - `init_invPush_invPull_idx` (1 valid / 2 invalid push / 3 invalid pull) and `valPushPull_invalPushPull_idx`.
  - `push_idx`, `pull_idx`, `left/center/right_idx` (trial indices).
  - `reach_times_inVec`, `reach_trials_inVec`, `cat_reach_inVec`, `hit_inVec`, `suc_inVec`, and `LCR_inVec_idx` (1 L / 2 C / 3 R).
  - `pp_inVec_idx` (0 push / 1 pull), `DCnD_*` (dominant vs. non-dominant), and `paw_pref`.
  - `reaches_inVec_px` / `_mm`: (201 × 3 × n_reaches).
  - **`win_selected`: the lab's stable good-behaviour + ephys window** (`trials_vec`, `select_start_time`,
    `time_sess_end`, Harp s). Set in `step3_bhv_struct.m`: start = trial-init time of the first selected
    trial (manual `cut_trials_start`), end = `behavior.behavior_duration` end (manual `shorten_display` in
    `daily_figure_PushPull_MovTarget_Free_v2.m`). `all_trials` holds the unrestricted events.
    Curation metrics use this window (step0 `ANALYSIS_WINDOW = "bhv_win_selected"`). To map it to imec:
    harp = a·imec + b, fitted from `st_raw_all` / `st_all` in the sorter folder's `neurons_session.mat`
    (v7.3, read with h5py). Example, CoteDor R5: trials 10–110 → imec 406.4–2025.0 s of 3481 s.
- **Reach records are 201 frames.** That is 150 before the **max-reach (endpoint) frame** plus 50 after, so the endpoint sits at 0-based index 150 and the window runs from −375 to +125 ms.
  - **`st_reach` t = 0 is the max-reach frame, not reach onset.** This is a property of the data, not an analysis choice.
- **Reach categories:** 1 = ballistic from the resting bar, 2 = lifted paw, 3 = grooming-like, 4 = lifted paw drink/reach.
  - "Purposeful" (`purpose_reach`): max-reach frame before the droplet disappears from DLC tracking (+10 frames; +40 for Toblerone).
  - "Hit" (`hit_reach`): among purposeful reaches, the **last one inside the water window** (the reach that took the water). Non-hit purposeful = earlier attempts.
  - "Success" (`success_reach`, step1 ~l.1250–1328, defined only for hits): the tracked droplet did **not** leave an ROI / pixel limits
    in a window [−100, +300] frames around the reach, i.e. the water was not knocked away. NaN for non-hits.
- **Reach segmentation exists** (step1/step2/step3 in MATLAB). Simon also had a separate onset detector in Python.

## 7. Size and access constraints (IMPORTANT)

- Raw electrophysiology (`.ap.bin` / `.lf.bin`, hundreds of GB) is **far too heavy for Claude to load or browse**.
  - Never write scripts that read a whole raw binary into RAM.
  - If raw data are needed (waveforms around CS times, or the sync channel), use `np.memmap` (385 int16 channels, interleaved), read short snippets only, process **one session at a time**, and cache results to small files.
  - The lab's `extractSyncChannel` reads the full LF file. Do not port that pattern.
- Claude normally works only from **lightweight files**: the sorter-folder `.npy`/`.tsv` files (§4), `params.py`, C4 tables, `cluster_channel_region_map.txt`, and small exports.
- The derived `.mat` files (`st_CB_DCN.mat` and similar) can be large: one comment mentions 67 GB of strings. Read them selectively with `h5py`, one session at a time.
- Data live on **Dropbox**, and cloud-only placeholders may need hydrating. `mapping_channel_to_regions.py` has an `ensure_file_available` helper for this.
- Scripts run **on Simon's workstation**. Claude usually sees only the code, printed outputs, small tables and figures. Scripts should print compact summaries and save diagnostic figures.
- Test new code on **one session first** (Milka S4/R4 is the natural choice), or on synthetic data.

## 8. Known issues found in the lab code (check before relying on it)

1. **Missing commas** in `paths_ephysTeremon_windows.py` (root and `ephys/` copies) after **Milka S1** and **Ferrero S5**, in both the CB and BG lists. Python silently joins each with the next string, so each list holds 25 entries instead of 27, and two of them are invalid paths. The MATLAB `.m` version is correct.
2. **Leading backslash** on some entries (`\2_Lindt\…S3`, `\4_Milka\…S1`). With `pathlib`, `root / r"\2_Lindt…"` resets to the drive root and discards `root_path`. Strip leading separators.
3. **Toblerone S2 imec mismatch.** The `.py` lists give `imec0` for *both* CB and BG. **Resolved: CB = imec1**, as in the `.m` list.
   Never take the `.py` path list as-is. Build our own session table (§9).
4. **Different spike sets across files.**
   - `neurons(i).st` is restricted to the behaviour epoch, but `eg_neurons(i).st` is the whole recording.
   - Templates were computed from behaviour-epoch spikes only.
   - For curation, compute every CCG, firing rate and shared-spike fraction on **one explicitly defined spike set and time window**, the same one for both units of a pair.
5. **Two clocks.** `spike_times.npy` is on the imec clock; every `.mat` spike time is on the Harp clock.
   - For pairwise CCGs within one probe the clock does not matter.
   - Never mix the clocks when aligning to behaviour, and never mix them across the two probes.
6. **The curation population is wider than the `.mat` files.** The `.mat` files keep only extra-good units and drop `potentialCS`. SS partners, and duplicate copies of a CS, may be good/mua units that are not extra-good. **Curation should therefore start from the sorter folders, not from `eg_neurons_CB.mat` / `st_CB_DCN.mat`.**
7. `s01` PSTH SEM divides by n instead of √n. It does not affect curation, but do not reuse those SEMs.
8. `step1_find_reaches…m` is committed in a debug state (`tt=1`, trial loop commented out).
9. `s04` warning lines reference the undefined `animals{m}`. They are harmless unless a path is missing.
10. **`step1` reach-timestamp bug (~l.956).** When a reach's +50 post-max frames run into the next video trial,
    `reach_timestamps = cat(2, frames_timestamps(tt,reach_curr), squeeze(ppDom_paw(tt+1,reach_nex)))` — the tail columns of
    `reaches.reach_timestamps_mat` are **paw values, not times**. Column 151 (max reach) is fine. Rebuild the tail as
    t_max + k/fps (k = frames after max) and count affected reaches. (The pre-max branch, reach in the first 150 frames, is correct.)
11. **Video-trial compaction risk (step1 ~l.126–137).** Empty videos (siz==1) and the following one are removed (`trials_w_vid`) and
    `siz` re-indexed, but `sync_from_log_all(tt)` is still indexed by the compacted tt → any removed video would shift all later
    frame↔Harp mappings. Not known whether it ever happened. Check: per-trial implied fps (≈400) and Harp `timestamp_reached` vs
    kinematic max reach (constant lag, no jumps).
12. **Simon's old `Whole_Session_Analysis_5_R5.py` timestamps.** It concatenates `frames_timestamps` over *all* `size_trial` entries,
    but `session.ppDom_paw_sess` starts at `start_trials` (trial 1 skipped when its video is > 7 min). Positions and times can be offset
    by one trial. Pair `trials.ppDom_paw(1:siz(tt),:,tt)` with `frames_timestamps(tt,1:siz(tt))` trial by trial instead.
13. **`step1` speeds are in px/frame** (`speed_DomPaw_sess`, SG derivative N=1, K=3, F=11, computed on the concatenated session across
    trial borders). px is anisotropic (y 0.07 vs x/z 0.1 mm/px), so px speed magnitudes are distorted → recompute speed from mm positions.
14. `step1` also drops the last trial when `n_trials == length(sync_from_log_all)` (l.711), after the position arrays were built:
    the last trial has positions but NaN `frames_timestamps`.

## 9. Decisions (Simon, 2026-11-03) and remaining open points

### Decided
- **Session table.** Only sessions in the CB curation list are included, using the sorter folder that carries the Phy + C4 labels. Build our own explicit session table (mouse, R#, S#, imec, sorter folder, sorter type) instead of reusing the buggy `.py` path list (§8.1–8.3).
- **Cell-type truth = C4.**
  - C4 was run on **all extra-good units** (and only those).
  - The manual Phy flag `cluster_CS.tsv` / `cluster_potential_CS.tsv` is *not* the classification. Keep it only as a diagnostic: report "manual CS but C4 ≠ PkC_cs" (e.g. unlabelled) and the reverse, but do not use it for selection.
- **First-pass population = extra-good CB units.**
  - CS candidates are C4 `PkC_cs`.
  - SS partner candidates are **extra-good units only** (`PkC_ss`, plus unlabelled extra-good units as candidates; whether to admit other C4 classes is still to decide).
  - Goal: measure how good the matching is. If the matching rate is poor, a second pass searches **good + mua** for partners that did not survive the extra-good selection. C4 labels do not exist for those units.
- **CS/SS mixing inside clusters** is accepted as a known, partly unfixable sorting limitation. Curation should *detect and report* it (e.g. zero-lag CS×SS peaks, CS-like events in SS clusters), not try to fix it.
- **Aggregation order:** curate per session, then pool sessions within an animal, then pool animals. No cross-session cell tracking.
  **Each session is a fresh acute insertion**, so the same PC cannot appear in two sessions and pooling sessions needs no de-duplication.
- **One-to-one CS↔SS is enforced.**
  - Two CS units pausing the *same* SS unit are treated as duplicates.
  - CS units pausing *different* SS units are treated as distinct, possibly synchronous, cells.
- **Method: candidates first, adjudication second.**
  - No hard thresholds or early decisions.
  - First, enumerate *all* candidate relations per session: CS–SS pause candidates and CS–CS duplicate candidates, each with its evidence scores.
  - Then adjudicate: resolve conflicts to the best-fitting one-to-one assignment, and record the decision and the competing candidates for every unit.
  - Keep rejected candidates and the reason for rejection in the outputs so that every call can be audited.
- **No automatic merging of CS duplicates (Simon, 2026-11-04).** First identify the *source* of each duplication
  (re-detection of the same event, fixed-lag spikelet/component, partial split, or true synchrony) from several
  independent evidence lines that agree. SS co-silencing is a key line: a shared SS partner, and spikes of one unit
  that the other missed still silencing that partner. step0 v2 writes `cs_cs_evidence.csv` (per-line votes, type
  hint, families) for this; it never merges.
- **SS partner search includes a tier 2** of non-extra-good CB units with FR ≥ 20 Hz (step0 v2). On R5 it doubled
  the matched CS (11 → 22 of 34).
- **Tier 2 has NO label gate (Simon, 2026-11-05; step0 v5).** `cluster_group.tsv` is not comparable across sorters:
  KS4 = Phy label (mostly KSLabel + manual edits; CoteDor S1 never labelled, 392 `unsorted`); IBL = a different,
  unclear origin (matches its own classifier argmax for only 18–70% of clusters) and ~2× more clusters. v2–v4 gated
  tier 2 on good/mua; v5 admits any label. Removing the gate added 91 KS4 tier-2 units but changed no pairing.
- **Output:** a curation table (CSV/parquet) keyed by `CB|mouse|sess|phyID`, plus sorter folder and imec. Saved **locally only** for now.
- **Environment:** a new, dedicated conda env (proposed name `teremon-cs`; Simon to confirm).
- **Thresholds** for duplicates and SS pause are not fixed in advance. They will be worked out in a dedicated exploratory script, informed by the literature (§12).
- **Blocks:** rule and droplet-location block information is in the st/eg `.mat` files and in `behavior_session.mat`. Field names are still to be identified when needed.
- **Kinematics:** same rig as the cooling experiment, but **keep this dataset's PX2MM = (0.1, 0.07, 0.1)** as is.
- **Reach alignment, first pass:** use the existing lab labels and detected reaches (max-reach-aligned `st_reach`, existing categories). See §13 for the planned extension.

### Still open
- Exact block field names (needed only for behavioural steps).
- Which C4 classes (if any) besides `PkC_ss` and unlabelled may serve as SS partners. For example, a `PkC_cs`→`MLI`-labelled unit with a clear pause could be an MLI-labelled SS.
- `STAB_MAX_GAP = 1` (two unstable bins end a span) may be conservative (R5 CS 622 stops at 2205 s, later bins look stable).
- `CS_SHARE_FLAG = 0.2` (provisional).
- **Tier 1 FR floor — DECIDED (Simon, 2026-11-06): tier 1 also needs ≥ `SS2_FR_MIN` (20 Hz)**, as tier 2 does.
  Reason: a loose pause threshold let 1–8 Hz extra-good units (often C4 MFB) become "SS partners".
  **Implemented in step0 v6 (2026-11-07)**: role `SS_lowfr`, audit table `cs_ss_pairs_lowfr_excluded.csv`
  (`would_be_candidate`). Consequence: CoteDor R5 SS 938 (14 Hz) is excluded (585/734/874→938 would be candidates);
  In v7 (3–8 ms) 585/734→938 pause at 0.07/0.06, i.e. among the best pairs of R5, and Ferrero R6 423→344 (13 Hz)
  is also cut; the other would-be pairs are 0.3–5 Hz units. A 10 Hz floor would keep 938/344 — OPEN for Simon.
  **Toggle (2026-11-08):** `SS1_FLOOR = "20" | "10"` at the top of step0 (or env `STEP0_SS1_FLOOR`, also read by
  `run_step0_batch.py` and `compare_step0_v5_v7.py`). "20" = default, folders unchanged; "10" writes to
  `step0_explore_bhvspan_<tag>_fr10` (batch logs / comparison likewise `_fr10`), so both settings coexist.
  R5 at 10 Hz (`1_CoteDor_R5/step0_explore_bhvspan_v7batch_fr10`): 938 becomes rank 2 for CS 585/734/874 (pr 0.07/0.06/0.40),
  behind 925; matched CS unchanged (22), candidate pairs 25 → 28.
  Full 10 Hz batch (2026-11-09, all 21 sessions, `v7batch_fr10`; comparison `_sorter_comparison/step0_v7_fr10_vs_v5/`):
  vs 20 Hz only 4 new candidate pairs (R5 585/734/874→938 at rank 2; Ferrero R6 423→344, pr 0.61, rank 1, newly
  matched CS) → matched CS 91 → 92. 938 carries a weak late CS-event piece (excess 0.07–0.12, lag ~1.8 ms; below
  the 0.2 flag). Other rank-1 changes are non-candidate units only. `fr10_vs_fr20_new_pairs.csv` lists the 4.
- **Pause window — DECIDED (Simon, 2026-11-07): 3–8 ms (step0 v7, current).** v6's 4–12 ms dropped the clearest v5 pairs
  (Toblerone R1 379→349 pr 0.07→0.71, R5 871→841 0.12→0.76): silence 0–7 ms, then a 1.3–1.8× rebound at 8–10 ms that
  the 4–12 window averages in. Read-only what-if (`_sorter_comparison/step0_v6_vs_v5/window_whatif/`): 3–8 or 4–8 ms
  keeps them (≤ 0.23) and gives 103 passing pairs vs 93 (4–12) / 94 (0–8). v7batch vs v5batch: matched CS 88 → 91 of 351,
  window lost 2 borderline pairs (Lindt R6 1190→921, Toblerone R1 412→401) and gained 13; floor lost 7
  (`_sorter_comparison/step0_v7_vs_v5/`).
- **KS4 vs IBL: parked (Simon, 2026-11-06).** Treated as doing equally well for now; no more time on the comparison.
  Background: it is not settled (sorter is confounded with session: each session was sorted by one sorter only). The
  proposed decisive test is to re-sort 2–3 sessions with the other sorter and run step0 on both outputs:
  e.g. KS4 Ferrero S6 → IBL; IBL CoteDor S5 or Toblerone S5 → KS4. Not started. Outputs must go outside Dropbox.

### KS4 vs IBL comparison: conclusions so far (2026-11-06; details in `FINAL_KS4_vs_IBL/README.md`)
- **CS detection is equal** (median 4 vs 4.5 per session).
- **IBL pairs more CS with an SS:** 31% vs 14% pooled; within-mouse p = 0.045 at the strict threshold.
  - The label gate is not the cause (v5, fig6).
  - Multi-unit dilution is not the cause (fig7): KS4 tier-2 units do have ~3× more ISI violations, but the units next
    to unpaired KS4 CS are clean and do not pause.
  - Beyond chance, loosening threshold and distance lifts KS4 from 7% to 18%; IBL stays ~30% (fig8). At the session
    level the loose comparison is not significant (p = 0.2–0.25).
- **More SS units with IBL, not significant** (C4 PkC_ss median 2 vs 4.5 per session; totals 35 vs 97).
- **CS spikes leak into neighbouring units more with KS4** (13% vs 5% of CS).
- Labels: compare sorters with the automatic `cluster_KSLabel.tsv` (same rule in both sorters), never with
  `cluster_group`.

### step0 script state (v5, 2026-11-05; v4 = unified, v5 = v4 + tier 2 without label gate, `SS2_GROUPS = None`)
`CS_curation/step0_explore_session.py` is the ONE live script (two chats had edited it; merged per
`step0_unification_plan.md`). What it does beyond v2 (tier 2 SS, duplicate evidence):
- **Per-unit stable spans** (`STABLE_SPAN = True`, section 3b). `win_selected` is the reference and is always kept; each CS/SS
  window is extended in 60 s bins inside the behaviour period while FR stays within ×/÷1.5 and amplitude ×/÷1.3 of the
  unit's reference values. A pair (CS×SS or CS–CS) uses the **overlap** of its two spans. `*_ref` columns
  (`candidate_ref`, `pause_ratio_ref`, `jitter_p_ref`, …) give the same metrics on `win_selected` only. In the v3 batch
  (21 sessions) the fraction of paired CS per session was identical with spans and with `win_selected` only.
- **Template overlay** (browser): CS and SS templates drawn on shared channels.
- **CS-spike sharing flag (flag only, never changes `candidate`).** Per CS×SS pair: `cs_share_frac` (CS spikes with an SS
  spike within ±0.5 ms), `cs_share_null` (Poisson chance for that SS rate), `cs_share_excess`, `ss_from_cs_frac`,
  `cs_share_flag = cs_share_excess > CS_SHARE_FLAG`. A flagged "SS" may carry CS spikes/spikelets or be a split piece
  of the CS (seen e.g. Lindt S6 CS 950 / unit 949, Toblerone S3 476 / 754). Shown in the browser (Pairing grid
  "CS shared", pair detail line, Tables).
- Robust to < 2 CS (empty `cs_cs_pairs`) and to no SS candidate (empty `cs_ss_pairs`); test hook `STEP0_TEST_NO_SS=1`.
- **Session / output naming (same rule for single and batch runs):** `RUN_TAG` constant (now `"v5"`), overridable only
  by env `STEP0_RUN_TAG`; env `STEP0_SESSION="<mouse>:<S#>"` picks the session from
  `_sorter_comparison/session_table.csv` (default without it: CoteDor R5). Output:
  `<mouse>_<R#>/step0_explore_bhvspan_<RUN_TAG>/`.
- **Batch runner** `run_step0_batch.py` runs a **frozen snapshot** (`_sorter_comparison/step0_<version>_snapshot/`),
  copied UNCHANGED from the live script (no hand patches; only the current one stays there, older ones go to
  `_sorter_comparison/_archive/`), one process per session, retry once (Dropbox
  hydration). Batch tag = env `STEP0_BATCH_TAG` (default `v5batch`), shared with the downstream scripts
  `sorter_compare_step0v3.py` (→ `_sorter_comparison/step0_<tag>/`) and `sorter_compare_zero_lag_check.py`
  (→ `zero_lag_check_<tag>/`). `sorter_compare_weak_pause.py` (same env) re-thresholds the stored step0 pairs:
  pause-threshold sweep, contamination, chance correction (→ `weak_pause_check_<tag>/`, fig7–8 copied into FINAL).
  `sorter_compare_final_figures.py` builds the central summary folder
  `_sorter_comparison/FINAL_KS4_vs_IBL/` (fig1–6, tables, README with the KS4 vs IBL conclusions).

**Output folders (data side; nothing is deleted).** Rule (Simon, 2026-11-06): keep only the CURRENT version in place;
every refinement step goes to an `_archive/` folder as soon as it is superseded, to avoid confusion.

| folder | produced by | status |
|---|---|---|
| `<mouse>_<R#>/step0_explore_bhvspan_v5batch` (21 sessions with CS) | v5 snapshot batch | **current reference** — still read by `raw_waveforms_*`, `FINAL_KS4_vs_IBL` and the downstream sorter-comparison scripts; archive only on Simon's OK |
| `<mouse>_<R#>/step0_explore_bhvspan_v7batch` (21 sessions with CS) | v7 snapshot batch (`_sorter_comparison/step0_v7_snapshot`, 2026-11-07): pause 3–8 ms + event decomposition + 20 Hz tier-1 floor | **current** |
| `_sorter_comparison/step0_v7_vs_v5`, `step0_v7_snapshot`, `step0_batch_logs_v7batch` | `compare_step0_v5_v7.py` | **current** |
| `_sorter_comparison/step0_v6_vs_v5` | `compare_step0_v5_v6.py`; README review notes + `window_whatif/` = the evidence for the v7 window | kept as reference |
| `CS_curation/_archive/step0_runs/<mouse>_<R#>/step0_explore_bhvspan_{v6batch,v6test,v7test}`, `_sorter_comparison/_archive/{step0_v6_snapshot,step0_batch_logs_v6batch}` | v6 (4–12 ms) and test runs | archived 2026-11-07 |
| `CS_curation/_archive/step0_runs/<mouse>_<R#>/` | all earlier step0 runs (v1, v2/`bhvwin`, mislabelled `bhvspan_v2batch`, v3, v3batch, v4, partial v4batch, `shareflag`, `*_v4test_*`, `*_v5test`); the 6 sessions without CS are moved there whole | archived 2026-11-06, not read by any script |
| `_sorter_comparison/step0v3` | sorter comparison on v3batch | previous reference (fig6 only) |
| `_sorter_comparison/step0_v5batch`, `zero_lag_check_v5batch`, `weak_pause_check_v5batch`, `FINAL_KS4_vs_IBL`, `step0_v5_snapshot`, `step0_batch_logs_v5batch` | sorter comparison on v5batch | **current** |
| `_sorter_comparison/*.csv` (`session_table`, `units_all`, `session_summary`, `stats_session_level`) | `sorter_compare_CS_SS.py` (step 1) | **current**, read by every downstream script |
| `_sorter_comparison/_archive/` | old snapshots (v3, v4), early-pass `figs`/`sessions`, `pass2_*`, `step0v2`, v3 `zero_lag_check`, old batch logs | archived 2026-11-06, not read by any script |

### Raw-data layer (CoteDor R5 only, 2026-11-06): CS events seen in the raw `.ap.bin`
Purpose: judge spikelets and CS pieces absorbed by other clusters, which the ±1.3 ms sorter templates cannot show.
- **2026-11-09: rebuilt on step0 v7batch** (`STEP0_DIR` in `raw_waveforms_session.py` / `raw_waveforms_figures.py`);
  pairing re-check windows now 0–8 (v5) / 3–8 (= v7) / 8–16 ms, `post_rank` by 3–8 ms. Previous (v5-based) output copied to
  `CS_curation/_archive/raw_waveforms_runs/1_CoteDor_R5/raw_waveforms_v5basis_2026-11-06/`.
- **`raw_dup_events.py`** (add-on, ~2.5 min): single raw events of every step0 CS–CS pair, ordered first → second unit,
  −3…+10 ms around the first unit's spike, 2 top channels per unit, every spike of both units marked; 40 'both' + 15
  first-only + 15 second-only events, lag histogram. Writes `raw_waveforms/raw_dup_events.js` (window.CS_RAW_DUPEV); the
  browser shows it as "Single raw events of CS A → B" in Duplicates and Raw data.
- **`ss_same_cell_check.py`** (2026-11-09, ~1 min, R5): are two SS units the same cell (soma + axon)? SS–SS CCG vs both
  ACGs, coincidences vs chance, cross-triggered raw means on a channel strip covering both, per-event "signature" fit
  (B's waveform at A + lag and vice versa) vs random times. Runs every pair of SS that are both candidates of one CS
  (v7batch or v7batch_fr10; ~3 min for 5 pairs); env `SSPAIR_A/SSPAIR_B` = one pair only. Output
  `<session>/ss_same_cell/<A>_<B>/` (png + summary.json) and `ss_same_cell/ss_same_cell.js` (window.CS_SSPAIR) →
  browser card "Same cell? SS A vs SS B" in the Pairing tab under the CCG grid (selector if several pairs).
  R5 results: 925–938 two cells (flat CCG, chance coincidences, no 938 shape on 925 spikes: 0.20 = random 0.20);
  510–628 two cells; 615/616/617 = CCG copies the ACGs (≈0.18× at ±1 ms), same waveform on ch 245 → likely ONE cell
  split in three by the sorter (not merged; Simon to decide).
- **Browser (2026-11-09):** every card has a hover ⓘ (what / how / read it; table `INFO` in browser_template.html,
  matched on the card title — a new card needs a new INFO entry); per-channel row plots (raw means, with/without,
  locked/other) use two columns. Python http.server pages are cached: after rebuilding data, hard-reload (Ctrl+F5).
- **Raw copy (local, read-only use):** `Data/ephys_and_behavior/raw_copies/1_CoteDor/04082023_CoteDor_Cer_S5_g0_imec0/`
  (80.4 GB, copied with `CS_curation/copy_raw_local.py`). Online-only Dropbox files cannot be copied (Windows times out
  the hydration): Simon first sets the file to "Make available offline" himself, then the copy runs at ~80 MB/s.
- `CS_curation/raw_waveforms_session.py` (~35 min) → `<mouse>_<R#>/raw_waveforms/`; then
  `raw_waveforms_figures.py` (~1 min) → `raw_waveforms/figures/` (session summary, per-CS events heat map, summary page,
  with/without grid, locked-vs-other grid, CS–CS pairs). The browser loads `../raw_waveforms/raw_data.js` if present
  ("Raw data" tab, sidebar badges, `pr 4–12` column in Pairing, raw card in Duplicates).
- Preprocessing: per-channel mean removed, common median reference, 300 Hz high-pass (order 3, zero-phase), µV.
  Spike times = template trough (checked: CS 585 raw trough at −0.03 ms).
- **Findings R5:** (1) 585's event is split: 734 = remote re-detection of the main spike (raw cosine 0.99), 874 = the
  spikelet at ~1.8 ms (present only in events where 874 fires), smaller late pieces in 601/938/589. (2) Pause measured
  after the event (4–12 ms) instead of 0–8 ms: 585/734/874 → SS 938 (0.30) instead of 925 (0.69); 610 → 617 vs 615
  (weak). (3) Single-channel per-event spikelet counts are weak in this session (background spikes): excess over a
  random-time null is ≤ 0.2/event beyond 1.2 ms except CS 564; use the with/without (conditional) means instead.
  (4) Raw CS–CS similarities agree with step0's duplicate verdicts on all 18 pairs.
- **All sessions (`raw_all_sessions.py`, 2026-11-09…11):** per session: pin the raw .ap.bin in Dropbox (`attrib +P`,
  narrow exception, see memory) + one read to start the download → copy to `raw_copies/<session>/` → `attrib -P +U` →
  the 4 scripts above (env `RAW_SESSION` / `RAW_BIN` / `RAW_META`) → delete the C: copy → `_raw_batch_done.json`.
  Logs `CS_curation/_raw_batch_logs/`. Done: CoteDor R5 (manual), R1, FerreroRocher R4. The other 18 failed on the
  Dropbox download (stalled from the 3rd session on, D: free space frozen); rerun when Dropbox downloads again (done
  sessions are skipped). Fixed on the way: empty `cs_cs_evidence.csv` (9 sessions), no SS pairs → no output folder,
  a failed download left the file pinned. Kilosort4 sessions have 383 sorter channels (`len(chan_map)`).
- **2026-11-16 rerun of the 18 (Simon: "go ahead")**, detached process, same log. Added `raw_cs_overlay.py` as a 5th
  step and a D: guard (`MIN_FREE_D` 40 GB: never pin a download that would leave D: below it — files set back to
  online-only have NOT freed D: so far: CoteDor R1, FerreroRocher R4 still local at 0x100020; CoteDor R4 had finished
  downloading and is left as it is). D: had 378 GB free for ~1.4 TB of downloads, so the guard will stop the batch
  after ~4 downloads unless Dropbox frees space. Old driver in `_archive/raw_batch_before_overlay_dguard_2026-11-16/`.
- **Dropbox sync state, found 2026-11-16:** this Dropbox uses its own Smart Sync driver (reparse tag 0x44), so the
  `attrib +P` / `-P +U` the batch used were ignored (never freed D:). Now: download = one read; afterwards Dropbox's
  own right-click "Make online-only" via `CS_curation/dropbox_free.py` (shell verb by exact name, only when Dropbox is
  not downloading — it was ignored during a download). The bytes then sit in `.dropbox.cache/old_files` until Dropbox
  clears it (~3 days). **Nothing under Dropbox is ever deleted, cache included (Simon).** The batch therefore WAITS for
  D: space (up to 5 days) instead of failing. Restarted 15:08 with this code (17 sessions; pre-change driver in
  `_archive/raw_batch_before_dropbox_free_2026-11-16/`). Expect ~4 sessions per ~3 days unless Dropbox clears sooner.
- **`pairs_raw_numbers.py`** → `_summary/pairs_raw_numbers.csv`: per session × floor (+ per animal, total) CS, paired
  CS, paired SS, pairs, before and after collapsing duplicates (counting only, nothing merged): CS = step0's
  `cs_dup_families.json`; SS = `ss_same_cell` verdicts "possibly the same cell" (raw layer only). Rerun after the raw batch.

### Pause figure + ephaptic check (2026-11-16; spike times only, all sessions)
After Han, Chen et al. 2020 (Nat Neurosci 23:1399): a CS pauses its own PC's SS for ~10+ ms AND briefly (~0–2 ms)
suppresses the SS of NEIGHBOURING PCs via its extracellular field (strong ≤ 25 µm, gone > 75 µm). So a 0–2 ms dip is
not evidence of a partner; step0's 3–8 ms pause window already avoids it (the old 0–8 ms window does not).
- `pause_figure_data.py` (~1–3 min/session) → `<session>/pause_figure/pause_figure.js` (window.CS_PAUSEFIG) +
  `ephaptic_pairs.csv`. Exported pairs = the browser's shown set (candidate | step0 rank ≤ 8 | combined top 8 | #1/#2 on
  pause, asym, jitter z, distance) + each CS's control SS (nearest non-partner, non-flagged). Per pair: raster (300
  events), PSTH 0.5 ms −10…+25 ms ÷ baseline (−50…−10 ms) with 95% bootstrap band, zoom −2…+5 ms at 0.2 ms, 4 quarters
  of the pair window. Every CS×SS pair: r_0_2 / r_3_8 = rate ÷ baseline, distance, partner, flagged.
- Browser: "Pause figure" + "Pause onset" in the pair detail (control SS menu, hash `ctrl`), "Ephaptic check" per CS in
  Pairing, "Ephaptic check (session)" in Tables. `ephaptic_summary.py` → `CS_curation/_summary/ephaptic_all_sessions.png/.csv`.
- Strict ephaptic measure (2026-11-16): `r_05_2` = 0.5–2 ms ÷ baseline (skips the sorter's overlap blind spot at 0 ms);
  `strict` = non-partner, non-flagged AND `cs_event_excess` and `cs_share_excess` ≤ 0.05. Caveat: `cs_event_excess`
  counts spikes −0.5…4 ms after the CS, so the filter can itself lower the early ratio; the >75 µm group is the check.
  Old scripts in `CS_curation/_archive/pause_figure_before_strict_2026-11-16/`. Browser cards still show the 0–2 ms set.
- Raw CS overlay (2026-11-16): `raw_cs_overlay.py` → `<session>/raw_overlay/raw_overlay.js` + `figures/cs<CS>_ss<SS>.png`.
  40 random CS events in the best partner's pair window, −10…+30 ms, CAR + 300–6000 Hz, on the CS and SS peak channels,
  with sorter spike ticks below. Per CS: partners (≤ 4; else the step0 #1 SS) + the control SS. Only sessions with a raw
  file on this computer (local copy, or a Dropbox file already on the device — never triggers a download): now CoteDor
  R5, R1, FerreroRocher R4. Browser: "Raw CS overlay" card in the pair detail, after the Pause figure.
- Pairing table (2026-11-16): "rank by: combined" = mean of the ranks on pause, asym, jitter z, distance; the top N is
  followed by any SS that is #1/#2 on one criterion (tagged). Display only: step0 `rank` unchanged.

## 10. Working conventions

- **Python** on a **Windows** workstation, in **VS Code with `# %%` cell-style scripts**. Use `pathlib` for paths, strip leading separators from relative path strings (§8.2), and never hard-code POSIX paths.
- Pipeline style: numbered steps (`step1_…`, `step2_…`). Each step reads the previous step's saved outputs and can be re-run on its own.
- **Never delete existing analysis steps or data files.** Write new outputs to new, clearly named folders, and keep data out of the code folder (lab convention).
- **Think through the methodology before writing code.** Simon often wants a design discussion or a well-specified prompt first.
- Be alert to **estimator / reference mismatches**, e.g. a baseline computed on different trials, time windows or spike sets than the quantity it normalises (see §8.4).
- Units:
  - Times are in seconds internally. State the clock (imec vs. Harp) and the sampling-rate source (`params.py`).
  - Kinematics are in mm (PX2MM above), not pixels.
- Lab colours (MI project): CB purple `[154,123,193]/255`, CP green `[106,176,134]/255`.

## 11. Curation pitfalls to keep in mind (not a plan)

- **CS duplicates vs. real climbing-fibre synchrony.** Neighbouring PCs in the same zone fire synchronous CS (±few ms), so a zero-lag peak alone does **not** prove a duplicate. Duplicates show a very sharp zero-lag peak, a large shared-spike fraction, nearby peak channels and similar templates (`similar_templates.npy` if present).
- **Within-cluster double detection.** One CS, with its spikelets, can be counted twice within a few ms. This shows up in the ACG as a peak below 5 ms.
- **CS hidden inside SS clusters**, or the reverse. This shows up as a zero-lag peak in the CS×SS CCG.
- **Sorter collision dips.** Overlapping spikes are removed, which leaves a narrow (~1 ms) zero-lag dip between nearby units. A true CS-induced SS pause lasts about 10–20 ms. The size of the dip may differ between KS4 and the IBL sorter.
- **One-to-one constraint.** A CS unit should match at most one SS unit, and vice versa. Violations point to duplicates.
- **C4-unlabelled units** may still be PC SS units. Candidate SS partners must include unlabelled, and possibly non-extra-good, units (§8.6).
- **Drift** affects every CCG-based criterion. Drift traces are available: IBL `drift.um.npy` and KS4 `ops['dshift']`.

## 12. Literature / prior art for curation (starting point, not decided)

Citations below are from memory, except the two code bases that were read directly. Verify the references before relying on them.

**CS–SS pairing (same Purkinje cell).**
- The standard criterion is a CS-triggered SS pause in the CCG, at positive lags only. The pause is typically about 10–20 ms and sometimes longer.
- **Lussac** (BarbourLab/lussac, `modules/find_purkinje_cells.py`, code read):
  - Putative SS units: FR ≥ 40 Hz. Putative CS units: 0.5–3 Hz.
  - CCG with 1 ms bins over ±25 ms.
  - Pause test: median CCG in [0, 8) ms < 0.4 × the median of the negative-lag baseline.
  - It then rejects symmetric dips, which removes collision/refractory artefacts.
  - It does not enforce one-to-one pairing and does not compare templates (there is a TODO in the code).
- **NeuroPyxels / npyx** (m-beau, the same group as C4; code read) has CCG significance tests, Stark & Abeles 2009 and Kopelowitz & Cohen 2014, as well as `ccg_vs_fr` and synchrony indices (`synchrony_zscore`, `synchrony_regehr`). The C4 paper (Beau et al., Cell 2025) used CS→SS pauses as part of its ground-truth PC identification.
- **P-sort** (Sedaghat-Nejad et al., J Neurophysiol 2022, Shadmehr lab) is semi-manual CS/SS sorting and confirms pairs by CS-triggered SS suppression.

**Duplicates.**
- Generic tools:
  - Shared-spike ("agreement") fraction within a tight window. SpikeInterface `remove_redundant_units`, as used by Lussac, defaults to delta = 0.3 ms with a duplicate threshold of 0.7.
  - Template similarity (`similar_templates.npy`).
  - KS4's own CCG-plus-similarity merging.
- **CS-specific problem:** CS waveforms are multiphasic and last several ms (spikelets) with a large spatial footprint. Two clusters can therefore lock onto different phases or channels of the *same* event and coincide at a **small, constant, non-zero lag**. A 0.3 ms window centred at 0 would miss them, so coincidence should be measured at the CCG peak lag.

**Duplicate vs. true climbing-fibre synchrony.** Candidate discriminants, to be combined:
- **Coincident fraction.** Duplicates: most spikes of the smaller unit are coincident. Synchronous neighbours: a minority.
- **Peak shape.** Duplicates show a sub-ms, jitter-free peak, possibly at a fixed offset. Synchrony shows a broader peak of ± a few ms.
- **Union-train refractoriness.** Merge the two trains after removing coincident pairs. Duplicates leave no new short ISIs; two distinct cells produce CS ISIs that are physiologically too short.
- **Amplitude co-variation** of coincident spikes, plus peak-channel distance and template similarity.
- **SS-partner logic, probably the most decisive and linking both problems.** Two CS units that pause the *same* SS unit point to a duplicate (one-to-one). Two CS units pausing *different* SS units are two synchronous cells. This fails when the partner is missing, which is likely for some units given the extra-good-only first pass.

**Pause vs. collision dip.**
- A collision dip is about 1 ms wide and symmetric around 0.
- An SS pause is asymmetric (positive lags), about 10–20 ms, and its depth is relative to the SS baseline rate.
- The baseline must come from the **same spike set and time window** as the pause estimate. This is an estimator/reference concern (§8.4).

**Statistics and drift.**
- Use a jitter or partially-hollowed-convolution null (Stark & Abeles 2009) instead of fixed ratios alone. Statistical power scales with the number of CS events (≈ 1 Hz × session length).
- Compute the criteria in time blocks and require stability, because SS rate and CS detectability drift. Drift traces: IBL `drift.um.npy`, KS4 `ops['dshift']`.

**Later (cross-session):** UnitMatch (van Beest et al., Nat Methods 2024), only if cross-session tracking is ever needed.

## 13. Deferred / planned extensions (keep in mind)

- **Reach segmentation consistent with the temperature project.**
  - Later, re-run that project's peak detection and reach classification on this dataset's kinematics.
  - Use the existing lab labels and detected reaches as seeds for that pipeline's **DTW-based onset refinement**.
  - Then **recompute spike times relative to the new windows** (onset-aligned instead of max-reach-aligned).
  - The first pass uses the existing lab reaches and labels only.
- **Second-pass SS partner search** in good + mua units if the extra-good matching rate is poor (§9).
- Pooling across sessions, then animals, after per-session curation.
- CS × reach-kinematics chapter: see §14.

## 14. CS × reach kinematics chapter (started 2026-11-16)

**Working mode (Simon):** the claude.ai chat is for design, literature and writing prompts; **code is written by the
Claude Code agent** from those prompts. Everything below comes from reading the lab MATLAB code (repo HEAD `4d6cf0e`,
2026-11-16) and was **not yet checked on the data**.

### 14.1 Goal and decisions
- Question: **when** (reach phase) and **where** (paw position) CS occur during the reach, and what they encode.
  Foundation paper: Hage et al. 2025, PLoS Biol 23:e3003110 (marmoset tongue, lingual vermis; see 14.5).
- **First pass: one mouse, one session — 1_CoteDor R5** (many CS units). Curation runs in parallel (§9); use one unit per
  duplicate family (step0 `cs_dup_families.json`) and keep the SS-partner information.
- **Space = absolute paw coordinates in mm** (dominant paw, x forward, y medio-lateral, z vertical; PX2MM §6).
- **Reach snippets are sufficient** (201 frames, −375…+125 ms around max reach, existing lab detection/labels). (Simon)
- **Droplet:** L/C/R from trial info; droplet position also tracked by DLC (see 14.3). Use both.
- **Licks: out of scope for now** (not tracked).
- **Deliverable for exploration:** an interactive HTML page to inspect kinematics and CS jointly (same pattern as the step0
  browser: Python writes a `.js` data file, one HTML template reads it).
- Scope of "CS units": C4 `PkC_cs`, extra-good → all present in `eg_neurons_CB.mat` (fields `st`, `phyID`, `cell_type`; v7.3).

### 14.2 Clock chain — everything is on the Harp clock
- **Spikes** (`eg_neurons.st` = `neurons.st_all`): `s01` fits `polyfit(eventTimes_water_giv(2:end) [imec, SpikeGLX bit 10],
  time_newTrial_log(2:end) [Harp], 1)` (or `(3:end)` with `flag_remove_1st_behavior_frame`) and applies it to all imec spike times.
- **Video frames:** `trials.sync.frames_timestamps`, interpolated between the same `time_newTrial_log` entries (§6).
- **Events:** `bhv.pp_times`, `trials.sync.timestamp_delivered/reached/iti_*` are `behavior.inputs.read_log` times (Harp).
- `bhv.reach_times_inVec` = `reaches.reach_timestamps_mat(:, 151)` restricted to `win_selected` trials.
- So spike ↔ kinematics alignment needs **no conversion** when spikes come from the `.mat` files. Spikes taken from the
  sorter folder (imec) need harp = a·imec + b (§6, `win_selected`).
- Checks to run on the first session (proposed): per-trial implied fps; border duplicates/gaps after stitching; `reach_mat`
  vs the continuous trace sampled at `reach_timestamps_mat`; Harp `timestamp_reached` − kinematic max reach of the hit reach
  (constant? drift over the session?); spike range inside video coverage; optional sorter imec→Harp refit per CS unit.

### 14.3 Where the needed fields are (`session_reaching_data_paw.mat`, v7.3; h5py arrays come transposed vs MATLAB)
- `trials.ppDom_paw` / `ppNonDom_paw`: max_siz × 3 × n_trials, **px**, pchip-filled after outlier/jump removal. `trials.size_trial` = siz.
- `trials.water`: max_siz × 3 × n_trials, DLC droplet ("water", bodypart 11) position in px, NaN when not tracked.
  `reaches.water_reach_mat`: same for each 201-frame reach record.
- `trials.water_location_trial` (1 L / 2 C / 3 R), `trials.is_push_trial`, indexed by **video trial tt**; behaviour trial = tt−1.
- `trials.sync.timestamp_delivered` = `read_log(trial_init_ind(tt−1)+1)` (droplet delivered); `timestamp_reached` =
  `read_log(reached_ind(tt−1))`; `timestamp_iti_start/end`; invalid push/pull times per trial.
- `reaches.reach_trial` = video trial tt; `reach_mat` 201×3×n (px); `reach_timestamps_mat` n×201 (Harp, see §8.10).
- Reach onsets (1-based indices inside the 201-frame record):
  - `start_forw.start_forw`: start of the last stable forward x-movement (x-speed > 1 px/frame for > 6 frames; at least 28 frames before max).
  - `start_stopX.start_stopX` (+ flag `onset_stop_mat`): movement start from rest in x (|vx| < 2 for ≥ 8 frames).
  - `start_realStop.start_realStop` (+ flag `onset_REALstop_mat`): start from rest in 3-D speed (< 200).
- Detection: `findpeaks` on x, height 0.84 × median droplet x, min distance 55 frames, prominence 40 (R-pawed) / 20 (L-pawed).
- `cat_reach`, `purpose_reach`, `hit_reach`, `success_reach` (§6), `nonDomPaw_resting`.

### 14.4 Analysis principles agreed (first pass, R5)
- Several alignments (droplet delivered, onset, peak speed, max reach) and a phase split in the spirit of Hage
  (pre-onset, acceleration, deceleration, post-max); compare with groups: L/C/R, hit vs non-hit, success vs fail,
  cat 1 vs grooming-like (cat 3).
- **Occupancy normalisation is mandatory** for "where": Pr(CS|x) ∝ p(x|CS)/p(x) (CS count ÷ time at x), as in Hage Fig 8.
- **Lag sweep** (e.g. −300…+300 ms) between CS time and paw position, instead of Hage's single −25 ms, to tell
  predictive from feedback CS (Streng et al. 2017).
- Null: circular shift of each spike train within the window (keeps ISIs); pooled null shifts units independently.
- CS are rare (~1 Hz, ~0.2–0.4 CS per reach per cell): show counts everywhere, coarse bins/smoothing.
- Confounds to keep visible: push/pull, droplet delivery, error feedback (buzz), reward; climbing-fibre reward signals.
- Later: Hage-style CS→SS-pause→kinematics test on paired cells; CS-on direction per cell; synchrony; BG probe.

### 14.5 Literature (verified 2026-11-16 unless marked)
- Hage et al. 2025 PLoS Biol — CS-induced SS suppression during deceleration → hypermetria; Pr(CS|x) maps, error CS after a missed tube.
- Kitazawa, Kimura & Yin 1998 Nature — monkey reach: early CS encode destination, late CS encode endpoint error.
- Streng, Popa & Ebner 2017 J Neurophysiol — CS encode kinematics and position error, mostly predictively.
- Wagner et al. 2021 Cell — mouse forelimb reach: synchronous CS bursts at movement onset emerge with learning.
- Calame, Becker & Person 2023 Nat Neurosci — mouse reach, PC population response scales inversely with reach velocity.
- Kostadinov et al. 2019 Nat Neurosci; Heffley et al. 2018 Nat Neurosci — reward-related CF signals in forelimb tasks.
- Ikezoe et al. 2023 Commun Biol — CF in a self-initiated lever-pull task multiplex lever movement, licking, reward.
- Tsutsumi et al. 2019 eLife — modular (zebrin) organisation of CF signals in a go/no-go task.
- Unverified (from memory): Herzfeld et al. 2018 Nat Neurosci (CS-on direction); Becker & Person 2019 Neuron (interpositus, reach deceleration).

### 14.6 Open
- Exact step0 file formats to read (`cs_dup_families.json` structure, `cs_ss_pairs.csv` columns) — check before writing the loader.
- Which onset definition to use as default (`start_stopX` with `start_forw` fallback was proposed).
- Phase boundaries and lag range: to fix after seeing R5.
