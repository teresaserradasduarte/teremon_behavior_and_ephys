# CLAUDE.md: granger_analysis

This is a directional (Granger-style) analysis of neural–kinematic coupling, using Neuropixels recordings (CP, CB, cortex) in head-fixed mice performing a select-and-collect task. The full design is in `PLAN.md`; read it first.

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
- `gctools/`: shared code.
- `scripts/0X_*.py`: parameters and calls only.
- **Data and outputs live outside the repo** (mirrors the MI project; paths in `gctools/config.py`):
  - pkls / cache: `<MAT_DIR>/LinearDependence_py/` (cache in `python_cache/`), next to `MutualInformation_py/`.
  - figures / tables: `<ROOT>/ephys_and_behavior/out_files/<GROUP>/group/LinearDependence/`.
  - First folder level = lag tag (e.g. `lag150_win-375to125`), so short- and long-lag runs coexist. Keep the tag in filenames where relevant.
- **Environment:** Python, with numpy / scipy / statsmodels. Ask before adding new dependencies.
