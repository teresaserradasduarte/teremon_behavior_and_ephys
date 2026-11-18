# Granger analysis: directional neural–kinematic coupling

**Question:** for each neuron (CP, CB, cortex, ...), does its past improve prediction of the kinematics, and/or the reverse? At what lag?
**Interpretation rule (Stokes & Purdon 2017):** results are *predictive precedence*, not physiological "driving". No spectral (Geweke) GC.

Sibling project of the MI pipeline, in the same parent repo, separate folder.

---

## 1. Folder layout

```
ephys/04_population_linear_dependence_py/   # (originally planned as granger_analysis/)
├── CLAUDE.md
├── PLAN.md                   # design (this file)
├── roadmap.md                # live checklist + open issues
├── gctools/                  # shared package; scripts only set parameters and call it
│   ├── config.py             # paths (roots from the MI config), constants, lag_tag, var_path / fig_dir
│   ├── io.py                 # load data through the MI loader (mi_io), trial table, neuron_table, REACH_LONG
│   ├── cache.py              # parameter-hashed disk cache for intermediate results
│   ├── preprocess.py         # binning, causal gamma smoothing, ZOH resampling, ensemble-mean removal
│   ├── design.py             # per-trial lagged matrices, discrete lag windows, raised-cosine bases, event regressors
│   ├── kernels.py            # 03's model: spec, causal bases (cross / own history / event), design,
│   │                         #   per-neuron fit both directions, CV over H, lag bootstrap
│   ├── models.py             # OLS VAR + BIC; Poisson GLM (IRLS); Gaussian AR; LR tests; trial-wise CV
│   ├── nulls.py              # trial-shuffle null (within condition)
│   ├── stats.py              # BH-FDR, trial bootstrap, bits/spike, lag readout from kernels
│   ├── runs.py               # shared run plumbing of 01-03: session list, run names, run_info, adaptive
│   │                         #   shuffles, FDR + direction classes, summaries, tables, common figures
│   ├── sim.py                # simulator (real kinematics + known ground truth)
│   └── plotting.py
├── scripts/                  # numbered = run order
│   ├── 00_core_checks.py         # milestone-1 evidence: causal kernels, raw vs matched delay, rasters, sims
│   ├── 00b_extract_reach_long.py # REACH_LONG: kinematics -3..+3 s around max reach from the session tracking
│   ├── 01_var_continuous.py      # Gaussian VAR baseline
│   ├── 02_ppglm_discrete_lags.py # Poisson GLM, discrete lag windows
│   ├── 02a_dispersion_check.py   # Poisson dispersion of every real neuron (02's reduced model)
│   ├── 03_ppglm_kernels.py       # MAIN: Poisson GLM, smooth causal kernels
│   └── 03b_explain_kernels.py    # explainer figures for 03 (how the model works; not results)
│                                 # 04_simulation_validation.py: NOT BUILT (sims run inside each script, §4)
└── tests/                    # plain-assert scripts; tests/run_tests.py runs all
                              #   (cache, causality, glm, io_real, kernels, models, nulls, sim, stats, trial_boundaries)
```
Data and outputs are **outside the repo** (see `gctools/config.py` and `CLAUDE.md` › Layout):
- `<MAT_DIR>/LinearDependence_py/<lag_tag>/<OUT_SUBDIR>/` pkls + `_neurons.csv` + `_run_info.txt`;
  `<MAT_DIR>/LinearDependence_py/python_cache/` parameter-hashed cache; `.../reach_long/` REACH_LONG (00b)
- `<ROOT>/ephys_and_behavior/out_files/<GROUP>/group/LinearDependence/<lag_tag>/<OUT_SUBDIR>/<RUN_NAME>/` figures and tables
- `OUT_SUBDIR`: `00_core_checks`, `01_VAR`, `02_GLM` (02 and 02a), `03_KERNELS` (03b in `03_KERNELS/explainers/`)

Scripts are cell-based (`# %%`). Cell 1 holds parameters and initialization; the remaining cells follow the MI project's style.

---

## 2. Conventions shared by all scripts

- **Causality:** everything that enters a predictor must use only the past.
  - Causal filters only.
  - Zero-order hold (not linear interpolation) when upsampling kinematics.
  - Backward differences for any derivative.
  - Event regressors only after the event.
- **Trial boundaries:** lagged matrices are built per trial and never cross trial boundaries. The first `max_lag` samples of each trial are dropped from the fit, but still used as history.
- **Analysis window:** `TRIAL_WINDOW_MS` is aligned to `ALIGN_EVENT` (e.g., reach onset). Set it once and share it across scripts.
- **Kinematic variable:** med-lat position (`KIN_VAR = 'y'`, coordinate 1 of reach_px, in mm). Others can be added later.
- **Trial-shuffle null:** pair neuron trial *i* with kinematics from trial *j ≠ i*, drawn **from the same condition** (`SHUFFLE_WITHIN`, e.g., spout position / action type).
  - This keeps each signal's autocorrelation and the task-locked structure.
  - What survives the shuffle is trial-by-trial coupling.
  - `N_SHUFFLE` is a parameter.
- **Both directions** are tested wherever possible:
  - **x→y (kinematics → neuron):** the neuron equation.
  - **y→x (neuron → kinematics):** the kinematics equation.
- **Comparing directions:** use shuffle p-values or normalized cross-validated gains, never raw deviances. The two equations use different likelihoods (Poisson vs Gaussian), so raw deviances aren't comparable.
- **FDR:** BH-FDR across neurons, applied per direction and per region.
- **Data source:** `DATA_SOURCE = 'real' | 'sim'`, plus `SIM_SCENARIO`, in every script's parameter cell (see §4).
- **Caching:** preprocessed inputs and fitted results are cached under `CACHE_DIR`, keyed by a hash of the relevant parameters. Load the preprocessed data from the MI pipeline wherever possible; don't recompute it from raw.

---

## 3. Scripts

### 01: Bivariate Gaussian VAR on continuous signals (didactic baseline)
This is not for final inference, because its known issues (Gaussian model on spike counts, stationarity) remain. Its purpose is to build and debug the machinery.

- **Neural signal:** spike counts smoothed with a **causal gamma kernel** (`FR_KERNEL_SHAPE`, `FR_KERNEL_SCALE_MS`). Report the kernel's mean delay.
- **Kinematics:** `KIN_SMOOTH = 'raw' | 'matched'`. `'matched'` applies the same causal gamma kernel, so both signals get the same delay.
- **Binning:** `BIN_MS` (default 10) after smoothing.
- **Optional:** `REMOVE_ENSEMBLE_MEAN` subtracts the across-trial mean trajectory and PSTH, per condition.
- **Fitting:** stacked per-trial OLS. Choose p by BIC up to `P_MAX`.
- **Output per neuron, both directions:**
  - GC as a log variance ratio
  - F-test p
  - Trial-shuffle p
  - Selected p
- **Figures:**
  - GC x→y vs y→x per neuron
  - `raw` vs `matched` comparison, which demonstrates the filtering artifact

### 02: Point-process GLM with discrete lag windows (sanity check)
- **Spikes:** binned at `BIN_MS = 2.5` (the same as the MI pipeline).
- **Kinematics:** med-lat, **raw**, upsampled with zero-order hold.
- **x→y:** Poisson GLM (`GLM_FAMILY = 'poisson' | 'negbin'`).
  - Fixed own-history windows `HIST_WINDOWS_MS`, e.g. `[(1,2),(3,5),(6,10),(11,20),(21,40),(41,80)]`.
  - Cross lag windows `CROSS_WINDOWS_MS`, e.g. `[(0,10),(10,20),(20,40),(40,80),(80,160),(160,300)]`, with the kinematic mean in each window.
- **y→x:** Gaussian AR on kinematics, with the same windows applied to the spike counts.
- **Ridge:** optional small penalty `RIDGE_ALPHA`, identical in the full and reduced models.
- **Tests:**
  - Likelihood-ratio test (full vs reduced) and trial-shuffle null, in both directions
  - BH-FDR across neurons
- **Output:**
  - Coarse lag profile (β per window)
  - Significance per direction

### 03: Point-process GLM with smooth causal kernels (main analysis)
- **Cross-kernels:** raised-cosine basis, **strictly causal**, fixed and identical for all neurons.
  - `BASIS_SPACING = 'log' | 'uniform'`
  - `CROSS_LAG_RANGE_MS = (0, 300)`
  - `N_CROSS_BASES`
- **Own-history kernel:** selected **per neuron on the reduced model**, by trial-wise cross-validation over `HIST_LENGTH_CANDIDATES_MS` (e.g., 20, 50, 100, 200). This guards against the kinematic lags absorbing unmodelled history.
- **Event kernels:** `EVENT_KERNELS = []` or `['endpoint']`.
  - All causal: they are active only after the event, over `EVENT_KERNEL_RANGE_MS`.
  - The list is extensible, so other events can be added later.
- **Cross-validation:** trial-wise K-fold (`N_FOLDS`), in both directions.
  - x→y gain is reported in bits/spike.
  - y→x gain is reported as the fraction of held-out variance explained, beyond the reduced model.
- **Tests:**
  - Likelihood-ratio test
  - Trial-shuffle null (`N_SHUFFLE`)
  - BH-FDR
- **Lag readout:** reconstruct β(τ) = B·w.
  - Report the lag of peak |β| and the kernel centroid.
  - 95% confidence intervals come from a **trial bootstrap** (`N_BOOT`): resample trials and refit.
- **Robustness checks (flags):**
  - log vs uniform bases
  - with/without the endpoint kernel
  - Poisson vs NB
- **Summaries:**
  - Per region (CB vs CP vs cortex): fraction of significant neurons per direction
  - Distributions of lags

### 04: Simulation-based validation
Script 04 is the systematic validation harness. It calls the **same `gctools` functions** as 03 (and optionally 01/02, to document their artifacts). See §4.

---

## 4. Simulation strategy (recommended)

**One simulator (`gctools/sim.py`), used two ways:**

1. **Quick check in any script:** set `DATA_SOURCE = 'sim'` and choose `SIM_SCENARIO` in the parameter cell. The script then runs unchanged on simulated data with known ground truth.
   - This lets you sanity-check each script with one flag change.
   - It avoids a separate simulation cell in every script.
2. **Script 04, the systematic sweep:** scenarios × rate regimes × trial counts × repeats. It produces:
   - False-positive rates
   - Power (sensitivity)
   - Lag-recovery error

**Real kinematics throughout:** simulations use real med-lat trials from a chosen session, so the autocorrelation and task structure are realistic.

**Scenarios:**

| Scenario | Construction | Expected result |
|---|---|---|
| `null_independent` | Spikes from a history-only GLM | No GC in either direction |
| `null_shared_drive` | Rate driven by the **condition-mean** trajectory, not the trial's own trajectory | Classic LR test may be significant; **the shuffle null must not be** |
| `x_to_y` | Rate driven by the trial's own kinematics through a known causal kernel peaking at lag `L` | x→y only; lag ≈ L |
| `y_to_x` | Real spikes (or simulated ones) with a known spike-driven term added to the real kinematics | y→x only; lag ≈ L |
| `bidirectional` | Both of the above | Both directions; both lags recovered |
| `common_latent` (limitation demo) | Trial-by-trial latent (e.g., vigor) scales both the rate and the kinematics | Documents what bivariate GC cannot rule out |

**Parameters swept:**
- Rate regime: MSN-like (sparse) vs Purkinje-like (high rate)
- Number of trials
- Coupling strength
- Lag L

**Pass criteria (to fix before running 03 on real data):**
- False-positive rate ≤ nominal α under both null scenarios
- Correct direction in the directional scenarios
- Median lag error less than about one basis width

---

## 5. Build order (plan review at each milestone)

Status as of 2026-11-16 (code exists and has been run; the live checklist is `roadmap.md`).

1. **`gctools` core:** io + cache + preprocess + design + sim (scenarios `null_independent` and `x_to_y` only). Check: lagged matrices never cross trials; filters are causal (unit test with an impulse). — **done** (`tests/`, evidence figures in `00`).
2. **01** on sim, then on real data. Demonstrate the `raw` vs `matched` difference. — **done** (sim; real 22 sessions).
3. **02** on sim, then on real data. — **done** (sim incl. 100–150 Hz null; real 22 sessions; dispersion in `02a`).
4. **03** together with **04**: add the remaining scenarios, then sweep. Fix parameter defaults from the 04 results before the final real-data run. — **partly**: 03 run on sim (no coupling, kin→neuron, neuron→kin, shared drive) and on the real 22 sessions; the remaining scenarios (`bidirectional`, `common_latent`) and the systematic 04 sweep (FP rate, power, lag error) are **not built**, so the pass criteria of §4 are not formally checked yet.

## 6. Decisions (2026-09-29) and open items
- `ALIGN_EVENT = 'max_reach'`: t=0 of `st_reach` is the max-reach frame (frame 151, 1-indexed; `step1_find_reaches…m`: `reach_goBack_frames=150`, `max_reach=151`). It is **not** reach onset.
- `TRIAL_WINDOW_MS = (-375, 125)`: the span covered by the tracked kinematics (201 frames at 400 fps). Spike history may reach earlier, within the same trial record.
- `CROSS_LAG_RANGE_MS = (0, 150)` for the first pass. Lags start at ≥ 1 bin (no lag 0). **To do later:** longer lags and windows (motivated by the MI results). This needs kinematics from the raw tracking. Keep both result sets, separated by the lag tag.
- `SHUFFLE_WITHIN = ('LCR',)`: spout side. `None` gives the plain permutation, kept as a secondary readout (it mixes side selectivity with within-side coupling).
- `REACH_CAT = (1,)` by default (parameter).
- First-pass session: `5_FerreroRocher R1` (339 cat-1 reaches, L/C/R 133/123/83, 37 CP / 49 CB). Second session: `1_CoteDor R6`.
- Regions: CP and CB only. CTX later.
- No per-region firing-rate assumptions: every region has mixed cell types, so rate-dependent decisions are made per neuron (`io.neuron_table`). Simulated rates are given in Hz, never by region.
- CB `cell_type` (in the .mat, not parsed by `mi_io`): skipped for now.
- Count model (01/10): dispersion was measured on all 2,009 real neurons (`02a_dispersion_check.py`). 92% of CP and 79% of CB are Poisson-like; neurons above 100 Hz (CB only) are UNDER-dispersed (median 0.69); over-dispersion is rare (~1%). Decision: keep 2.5 ms bins (= camera frames); **no negative-binomial** (it cannot model under-dispersion); the shuffle test is validated at 100–150 Hz. 03 reports effect sizes with and without neurons > 100 Hz; a 1 ms bin option is built only if that changes a conclusion.
- Own history in 03 (01/10): kinematics→neuron uses per-neuron CV over H ∈ {150, 200, 300} ms (never shorter than the cross range, after the 01 artifact); neuron→kinematics uses a fixed 150 ms (longer would cut rows); `HIST_MODE = "fixed"` uses 150 ms for both, as a check; an adequacy test flags neurons whose residuals still depend on older own history. *(The neuron→kinematics part was superseded by REACH_LONG, below.)*

### Decisions added since (2026-10-01; current defaults in the scripts' cell 0)
- Sessions and neurons: `SESSIONS = "all"` = every session with ≥ `MIN_STRATUM` = 20 trials per spout side (22 of 29 with reach category 1); neurons ≥ `MIN_RATE_HZ` = 1 Hz. BH-FDR q < 0.05 per region × direction, pooled over sessions. `N_SHUFFLE = 1000`, with **adaptive shuffles** (neurons that beat all 1000 shuffles but are not yet BH-significant get more, sequential p of Besag & Clifford, capped where more cannot change the BH decision; `runs.adaptive_extend`).
- **REACH_LONG** (`00b`): the same reaches cut from the continuous session tracking over −3..+3 s around max reach; cut back to −375..+125 ms it reproduces `reaches.reach_mat` and the MI `reach_px` exactly (checked per session in 00b). Timing at video-trial boundaries is uncertain by ~1–2 frames and flagged per reach. `KIN_SOURCE = "long"` gives **real kinematic history before −375 ms**; the analysis window is unchanged.
- 03 defaults (real 22-session run `real_22sessions_kernels_kinLong_kinHmax900_from-225`):
  - neuron own history: per-neuron CV over H ∈ {150, 200, 300} ms (`HIST_MODES = ("cv",)`).
  - kinematic own history: `KIN_HIST_MODES = ("max",)` over `H_KIN_LIST = (450, 600, 900)` → 900 ms (in sims CV always picked the longest; 900 ms clears the adequacy flag). Needs `KIN_SOURCE = "long"`.
  - `FIT_FROM_MS = -225`: fitted rows start at −225 ms, the same rows as the short-kinematics runs (`None` would start at −375 with REACH_LONG).
  - own-history log bumps scale with H (`SCALE_BASES`); 6 log cross bases over (0, 150] ms; ridge 0.1; 5 trial folds stratified by spout; adequacy over (H, H + 150] ms, flag p < 0.01; lag bootstrap 200 (FDR-significant neurons only).
  - `EVENT_KERNEL = False`, `TASK_MEAN_NUISANCE = False` (both built; the nuisance term was tried on sims).
- 03 effect size: **primary = EXCESS over the shuffle** (observed − mean of the neuron's own shuffle null; bits/spike for kin→neuron, % partial R² for neuron→kin) = the trial-specific part, ~0 under the shared-drive null. Secondary = cross-validated gain of full over reduced (includes the shared reach shape of the spout side).
- 01 lag readout is a PREVIEW (smooth causal cross-kernel on the VAR; use the centroid, the peak is unreliable — `models.LAG_NOTE`).

### Open issues
- **Long lags not run yet**: everything on disk is `lag150_win-375to125`. Longer lags / windows (motivated by the MI results) are the next step; keep both result sets under their lag tags.
- `<LAG_TAG>/01_var_continuous/` (VAR_DIR and OUT_DIR, 29–30 Sep) holds early 01 runs from before the output folder was renamed `01_VAR`; the current 01 runs are in `01_VAR/`.
- Script 04 (systematic simulation sweep) and the `bidirectional` / `common_latent` scenarios are not built; the §4 pass criteria are not formally checked.
- CTX not analysed; CB `cell_type` not used.
- 03 prints its summaries with and without neurons > 100 Hz (decision above); a 1 ms bin option is built only if that changes a conclusion — not yet assessed.
