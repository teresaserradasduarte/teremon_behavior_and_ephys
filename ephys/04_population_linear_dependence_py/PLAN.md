# Granger analysis: directional neural–kinematic coupling

**Question:** for each neuron (CP, CB, cortex, ...), does its past improve prediction of the kinematics, and/or the reverse? At what lag?
**Interpretation rule (Stokes & Purdon 2017):** results are *predictive precedence*, not physiological "driving". No spectral (Geweke) GC.

Sibling project of the MI pipeline, in the same parent repo, separate folder.

---

## 1. Folder layout

```
granger_analysis/
├── CLAUDE.md
├── PLAN.md
├── gctools/                  # shared package; scripts only set parameters and call it
│   ├── io.py                 # load preprocessed data (reuse the MI project's loaders), trial table
│   ├── cache.py              # parameter-hashed disk cache for intermediate results
│   ├── preprocess.py         # binning, causal gamma smoothing, ZOH resampling, ensemble-mean removal
│   ├── design.py             # per-trial lagged matrices, discrete lag windows, raised-cosine bases, event regressors
│   ├── models.py             # OLS VAR + BIC; Poisson/NB GLM; Gaussian AR; LR tests; trial-wise CV
│   ├── nulls.py              # trial-shuffle null (within condition)
│   ├── stats.py              # BH-FDR, trial bootstrap, bits/spike, lag readout from kernels
│   ├── sim.py                # simulator (real kinematics + known ground truth)
│   └── plotting.py
├── scripts/
│   ├── 01_var_continuous.py
│   ├── 02_ppglm_discrete_lags.py
│   ├── 03_ppglm_kernels.py
│   └── 04_simulation_validation.py
└── tests/                    # plain-assert scripts; tests/run_tests.py runs all
```
Data and outputs are **outside the repo** (see `gctools/config.py`):
- `<MAT_DIR>/LinearDependence_py/<lag_tag>/` pkls, `<MAT_DIR>/LinearDependence_py/python_cache/` parameter-hashed cache
- `<ROOT>/ephys_and_behavior/out_files/<GROUP>/group/LinearDependence/<lag_tag>/<OUT_SUBDIR>/<RUN_NAME>/` figures and tables

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

1. **`gctools` core:** io + cache + preprocess + design + sim (scenarios `null_independent` and `x_to_y` only). Check: lagged matrices never cross trials; filters are causal (unit test with an impulse).
2. **01** on sim, then on real data. Demonstrate the `raw` vs `matched` difference.
3. **02** on sim, then on real data.
4. **03** together with **04**: add the remaining scenarios, then sweep. Fix parameter defaults from the 04 results before the final real-data run.

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
- Own history in 03 (01/10): kinematics→neuron uses per-neuron CV over H ∈ {150, 200, 300} ms (never shorter than the cross range, after the 01 artifact); neuron→kinematics uses a fixed 150 ms (longer would cut rows); `HIST_MODE = "fixed"` uses 150 ms for both, as a check; an adequacy test flags neurons whose residuals still depend on older own history.
