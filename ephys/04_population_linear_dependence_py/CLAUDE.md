# CLAUDE.md: granger_analysis

This is a directional (Granger-style) analysis of neural–kinematic coupling, using Neuropixels recordings (CP, CB, cortex) in head-fixed mice performing a select-and-collect task. The full design is in `PLAN.md`; read it first. The live checklist and open issues are in `roadmap.md`.

## Working rules
- **Plan review before any non-trivial change.** Use plan mode, and wait for approval.
- **No large refactors without discussion.** This applies to `gctools/` and to the sibling MI project.
- **Iterative, evidence-driven development.** Validate every new component on simulated data (`DATA_SOURCE = 'sim'`) before running it on real data.
- **Match the MI project's coding style.** Before writing code, read its scripts and shared package. Mirror its:
  - cell layout (`# %%`, with the first cell holding parameters and initialization)
  - naming
  - docstrings
  - config style
  - plotting conventions
- **Reuse the MI project's data loaders and preprocessed outputs.** Don't reimplement loading. Import from the MI project. If that isn't possible, raise it instead of copying code.
- **Cache expensive steps** under `gctools.config.CACHE_DIR`, keyed by a hash of the parameters. Never recompute preprocessed data that already exists.

## Non-negotiable analysis rules
1. **Causality of predictors.** Any predictor may use only the past:
   - Causal filters only; no `filtfilt`, no Savitzky–Golay, no centred Gaussian.
   - Zero-order-hold resampling.
   - Backward differences.
   - Event kernels active only after the event.
   - Every new filter or resampler gets an impulse test that proves it is causal.
2. **No lag crosses a trial boundary.** Lagged matrices are built per trial.
3. **Full and reduced models are identical except for the tested terms.** Same data rows, same own-history, same penalty.
4. **Significance comes from the trial-shuffle null** (shuffled within condition), with BH-FDR across neurons. The χ² likelihood-ratio p is reported, but it is secondary.
5. **Don't compare raw deviances across directions.** Use shuffle p-values or normalized cross-validated gains.
6. **Wording:** "predictive precedence", not "drives" or "causes". No spectral GC.

## Layout
- `gctools/`: shared code. `tests/`: plain-assert scripts (`tests/run_tests.py` runs all).
- `scripts/0X_*.py`: parameters and calls only. The full file list is in `PLAN.md` §1.
- **Data and outputs live outside the repo** (mirrors the MI project; paths in `gctools/config.py`).
  The roots (`ROOT_DIR`, `GROUP`, `MAT_DIR`) come from the MI `config.py`, loaded by file path as
  `mi_config` (and `mi_io`), so the MI env overrides `MI_ROOT_DIR` / `MI_DATA_MAT` apply here too:
  - pkls: `VAR_DIR = <MAT_DIR>/LinearDependence_py/`, next to `MutualInformation_py/` (`config.var_path(*parts)`).
  - cache: `CACHE_DIR = <VAR_DIR>/python_cache/`, keyed by a parameter hash (`gctools/cache.py`).
  - long-window kinematics: `<VAR_DIR>/reach_long/reach_long_win-3000to3000.pkl`, made by `00b`
    (`io.load_reach_long()`).
  - figures / tables: `OUT_DIR = <ROOT>/ephys_and_behavior/out_files/<GROUP>/group/LinearDependence/` (`config.fig_dir(*parts)`).
  - Folder levels, under both `VAR_DIR` and `OUT_DIR`: `<LAG_TAG>/<OUT_SUBDIR>/` (figures one more level: `<RUN_NAME>/`).
    - `LAG_TAG` (e.g. `lag150_win-375to125`, from `config.lag_tag()`), so short- and long-lag runs coexist.
      Keep the tag in filenames where relevant.
    - `OUT_SUBDIR` = `00_core_checks`, `01_VAR`, `02_GLM` (02 and 02a), `03_KERNELS`.
    - `RUN_NAME` is built by `runs.data_name` from the data source, sessions and settings
      (e.g. `real_22sessions_poissonWindows`); each run also writes `<RUN_NAME>_run_info.txt` and `_neurons.csv`.
  - Windows paths are long: `config.long_path` / `warn_long` handle names over 259 characters.
- **Environment:** conda env `neuro_mi` (as the MI project), with numpy / scipy / statsmodels / pandas /
  matplotlib / joblib. Ask before adding new dependencies.
