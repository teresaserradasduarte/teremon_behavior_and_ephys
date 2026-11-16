# MI project — context for Claude Code

Persistent brief so any Claude Code session starts proficient. `CLAUDE.md` imports this
via `@MI_PROJECT_CONTEXT.md`. Keep this high-signal; put the live checklist in `roadmap.md`.

## Layout — code and data are SEPARATE

- This repo holds CODE only. Raw data and all derived outputs live under a DATA ROOT set in
  `config.py`, mirroring the MATLAB `rootdir / group / out_files` layout. Never write data into
  the code folder.
- `config.py` defines `ROOT_DIR`, `GROUP`, `SETUP`, `ANIMALS`, `REGION_STRUCT`, `REGIONS`, and
  derives `MAT_DIR`, `OUT_DIR`, `DATA_MAT` (default `MAT_DIR/st_eg_neurons.mat`), `VAR_DIR`,
  `CACHE_DIR`. Override per-machine without editing via env vars `MI_ROOT_DIR` and `MI_DATA_MAT`.
  - intermediate variables / pickles: `VAR_DIR = MAT_DIR/MutualInformation_py` (`config.var_path(name)`),
    parsed cache in `VAR_DIR/python_cache` (`config.cache_path()`);
  - figures: `OUT_DIR = .../out_files/<GROUP>/group/MutualInformation_NEW` (`config.fig_dir(*parts)`).
- Behaviour lives in OTHER trees: `BEHAV_MAT_DIR` (group `BEHAV_GROUP` = `20230511_ChocolateGroup`,
  water position), `BEHAV_RAW_DIR` (raw load cell + initiation times), and per-session
  `behavior_fundamentals.mat` in `MAT_DIR`. Same mouse/session names underneath.

## Scientific goal

Compare how, and *when*, basal ganglia (striatum / CP) vs cerebellum (CB) neurons carry
information about forelimb reach kinematics, in a head-fixed "select-and-collect" task
(push/pull to trigger water = select; reach to a spout = collect; rules alternate in uncued
20–30 trial blocks). Prior result to re-test with a better estimator: on average CB *leads* and
BG *lags* reach-endpoint kinematics. Claude owns the plumbing; the human owns the inference.

## Data

- `config.DATA_MAT` = a .mat with variables `all_<REG>_sts`, one struct array per region, one
  entry per neuron. **These `_sts` structs carry `st_reach` spike times** — NOT the MATLAB
  `all_<REG>_neurons` variant (spikes removed there). Regions of interest: CP (striatum, low-FR
  MSNs), CB (cerebellum, high-FR). Others (P, DCN, BS, MB, CTX) secondary.
- Load through `mi_io.parse_mat`, which handles both formats: classic .mat via
  `scipy.io.loadmat(path, squeeze_me=True, struct_as_record=False)` (squeeze collapses singletons →
  always normalise) and v7.3 via h5py. **v7.3 trap:** an empty MATLAB array is stored as its
  dimensions `[0, 1]` with attribute `MATLAB_empty = 1` — always check the flag, or an empty trial
  reads as spikes at 0 s and 1 s (bug fixed 2026-10-12, see `19_spike_placeholder.py`).
- `reach_px` frame 150 (0-based) = the REACH ENDPOINT = st_reach t = 0 (`config.REACH_ENDPOINT_FRAME`,
  a property of how MATLAB cut the reaches). px → mm = `config.PX2MM` = (0.1, 0.07, 0.1) for x, y, z.
- `st_init{k}` = spike times relative to push/pull `k` (t = 0 = rig detection = cue; the paw
  movement is BEFORE 0, water follows at a fixed +222 ms). Codes `idx_init_vPP_invPP`:
  0 valid push, 1 valid pull, 2 invalid push, 3 invalid pull.
- `cell_type` exists only for CB (PkC_cs, PkC_ss, MFB, GoC, MLI, NaN = unlabelled); `mi_io` does not
  keep it (24 reads it directly).
- Per-neuron fields (16): `phyID`, `st_init`, `st_reach`, `brain_region`, `cell_type`,
  `idx_init_vPP_invPP`, `DCnD_init`, `idx_reach_cat`, `idx_reach_hit`, `idx_reach_succ`,
  `idx_reach_hDCnD`, `idx_reach_LCR`, `idx_reach_PP`, `reach_px`, `mouse`, `sess`.
- `st_reach{k}` = spike times (s) relative to reach `k`, spanning ~[-3, 3].
- `idx_reach_*` = 1×n_trials per-trial labels (hit, success, L/C/R, push/pull, category).
- `reach_px` = `201 × 3 × n_trials` (frames × x,y,z × trials), **identical for every neuron of a
  session**, so `reach_px[:, c, k]` aligns with `st_reach{k}`. **Med-lateral = y = column index 1**
  (0-indexed). **Endpoint window = frames 147–167** (MATLAB 1-indexed → Python `146:167`).
- Struct concatenates sessions/mice; neurons of one `(mouse, sess)` share trials.
  uid = `f"{REG}|{mouse}|{sess}|{phyID}"`.

## Analysis design — fixed decisions (do not silently revert)

- Neural variable = **spike COUNT in the integration window `[start, stop]`** (the old MI
  landscape axes). NOT instantaneous 1/ISI rate. ISI-artefact-free, degrades gracefully.
- Neural window width = **70 ms** (step 25 ms), chosen from the width sweep in `20`.
- Bin at **2.5 ms** then sum over the window (matches 400 fps camera grid). **No smoothing.**
  Binarising bins to 0/1 is ~irrelevant — keep raw counts.
- **Empty trials** (0 spikes in window) are a (neuron × window) property, expected from Poisson:
  `P(empty)=exp(-r·T)`. Diagnose as observed-vs-Poisson. **IMPORTANT: do NOT drop empty trials**
  (it conditions on firing and biases MI; silence is informative). Prefer widening the window.
- **Estimator (load-bearing):** primary = **KSG both-continuous**
  (`mutual_info_regression`, count→float + tiny jitter) — uniform across CP low and CB high
  counts, so the cross-region comparison stays fair. Cross-check only = **Ross**
  (`mutual_info_classif`, count discrete) — destabilises for high-FR CB (kNN within sparse
  count-levels starves). **Avoid the binning plug-in for headline numbers** (positive bias scales
  with count-levels → worse for CB, confounds comparison). Report MI in **bits** (nats ÷ ln 2).
- `mutual_info_classif/regression` give **per-feature** MI, not **joint**. For multi-D kinematics
  use a true multivariate KSG on the full continuous matrix (`09`'s `_ksg_md`, Chebyshev norm;
  batched GPU version in `gpu_mi.py`, validated in `validate_gpu_mi.py`).
- Multi-D feature set = **x, y, z, speed, dist_water** (mm; speed by Savitzky-Golay derivative),
  chosen by excess MI in `11` (originally planned [PC1, PC2, speed, accel]).
- Compare and test on **EXCESS MI = observed − mean of trial-shuffle nulls** (raw KSG MI rises with
  dimension through bias alone). Statistics must be measured against the shuffle null, not within
  a map (`_validate_contrast.py`).

## MI formulas (reference)

- Old plug-in (avoid): `I = Σ p(x,y) log2[p(x,y)/(p(x)p(y))] = H(X)+H(Y)-H(X,Y)`; the three
  histograms are *only* what binning consumes — KNN estimators build none.
- KSG: `I ≈ ψ(k) − ⟨ψ(n_x+1)+ψ(n_y+1)⟩ + ψ(N)`.
- Ross: `I ≈ ψ(N) − ⟨ψ(N_d)⟩ + ψ(k) − ⟨ψ(m)⟩`, kNN restricted to same-count trials.

## Estimator wrapper convention (reuse)

```python
NAT2BIT = 1.0 / np.log(2)
def mi_bits(kin, count, rng):                      # kin: (N,) or (N,d); count: (N,) int
    X = np.atleast_2d(kin).T if kin.ndim == 1 else kin
    ksg  = mutual_info_regression(X, count + rng.normal(0,1e-9,count.size),
                                  discrete_features=False, n_neighbors=3, random_state=0) * NAT2BIT
    ross = mutual_info_classif(X, count.astype(int),
                               discrete_features=False, n_neighbors=3, random_state=0) * NAT2BIT
    return ksg, ross                               # per-feature; sum is NOT joint MI
```

## Repo files

Shared modules
- `config.py` — paths/constants; data + outputs live under the data root (NOT here).
- `mi_io.py` — loader (scipy + h5py): `parse_mat(path, struct_names)`, `counts(neuron, win)`,
  `save/load`. Keeps spike times so counts recompute for any window.
- `mi_rates.py` — neural variable modes: `count` (default), `isi` (MATLAB 1/ISI port), `conv`.
- `gpu_mi.py` — batched torch KSG / Ross + shuffle-null excess MI (CPU fallback).
- `mi2d_figures.py` — analysis + figures of the TIME × TIME maps (12's cells 4b–7 as functions),
  used by 18. 12 still has its own copy: **a change in one must be made in the other.**

**Core scripts: `04`, `09`, `12`, `15`, `18`, `28`. The most important is `12`** (the TIME × TIME
map); 15 and 18 are its onset and push/pull versions, so changes to 12 must stay comparable with them.
The others are diagnostics, checks, or steps toward these.

Scripts (numbered = run order; later ones load earlier pickles)
- `01` .mat → cache. `02` empty trials by category. `03` MI three ways on example cells.
- `04` 1-D MI curve per neuron (endpoint y; binning / KSG / Ross), lead/lag. `05` neural start × stop map (GPU).
- `06`/`07` count vs isi vs conv. `08` 04 per mouse.
- `09` joint multi-D KSG curve. `10` water-referenced features. `11` compare feature sets.
- `12` **TIME × TIME map**: MI(kinematics in kin window ; count in neural window), t = 0 = reach endpoint.
  `15` same, aligned to reach onset. `16` push/pull kinematics. `21` clean push/pulls + force-onset
  alignment. `18` TIME × TIME around the push/pull (`st_init`).
- `13` standalone lead/lag figure. `14` why binning is biased. `17` MI curves next to raster + scatter.
- `19` the v7.3 empty-trial bug, before/after. `20` neural window width sweep (knee).
- Map shapes (what makes up the mean bands): `22` split-half profiles, `23` hot cluster, `24` CB cell
  types, `25` template bank, `26` pick maps to label, `27` browser labelling GUI (run from a terminal),
  `28` classifier trained on the hand labels.
- Checks: `validate_gpu_mi.py`, `_validate_contrast.py`.

If a file is missing, regenerate from the decisions above, not from memory of older versions.

## Roadmap (live checklist in `roadmap.md`)

1–5 done: empty trials, deep-dive, 1-D bridge, multi-D joint KSG, temporal landscapes (05, 12/15/18).
6 in progress: CP-vs-CB lead/lag — currently classifying single-neuron map shapes (22–28).
Open: 09/12 use frame 147 as the kinematic zero, 10 ms before the endpoint (frame 150);
09 still uses a 100 ms neural window (chosen width is 70 ms).

## Conventions

- Colours: CP green `[106,176,134]/255`, CB purple `[154,123,193]/255`.
- Cell-based scripts (`# %%`) for the VS Code interactive window; small composable scripts that
  import `config` and `mi_io`. Run `01` once per data update; later scripts reuse the in-kernel
  `NEURONS` or the cache. Conda env: `neuro_mi` (`neuro_mi.yml`).
- Validate any new estimator/loader on tiny synthetic data with a KNOWN answer first.
- Show evidence (numbers, figures), don't assert. Ask before large refactors, before changing a
  fixed decision, before deleting files.
