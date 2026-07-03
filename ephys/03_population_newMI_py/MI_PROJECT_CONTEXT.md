# MI project — context for Claude Code

Persistent brief so any Claude Code session starts proficient. `CLAUDE.md` imports this
via `@MI_PROJECT_CONTEXT.md`. Keep this high-signal; put the live checklist in `docs/roadmap.md`.

## Layout — code and data are SEPARATE

- This repo holds CODE only. Raw data and all derived outputs live under a DATA ROOT set in
  `config.py`, mirroring the MATLAB `rootdir / group / out_files` layout. Never write data into
  the code folder.
- `config.py` defines `ROOT_DIR`, `GROUP`, `SETUP`, `ANIMALS`, `REGION_STRUCT`, `REGIONS`, and
  derives `MAT_DIR`, `OUT_DIR`, `DATA_MAT`, `CACHE_DIR`. Override per-machine without editing via
  env vars `MI_ROOT_DIR` and `MI_DATA_MAT`. Parsed cache + figures go under
  `OUT_DIR/python_cache` and `OUT_DIR/<name>`.

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
- Load with `scipy.io.loadmat(path, squeeze_me=True, struct_as_record=False)`; fields are
  attributes; squeeze collapses singletons → always normalise (helpers in `mi_io.py`).
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
- `mutual_info_classif/regression` give **per-feature** MI, not **joint**. For multi-D
  [PC1, PC2, speed, accel] use a true multivariate KSG/Ross on the full continuous matrix.

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

- `config.py` — paths/constants; data + outputs live under the data root (NOT here).
- `mi_io.py` — squeeze-safe loader: `parse_mat(path, struct_names)`, `counts(neuron, win)`,
  `save/load`. Keeps spike times so counts recompute for any window.
- `01_load_data.py` — read .mat once → cache `parsed_<REGIONS>.pkl`; print neurons/sessions/trials.
- `02_empty_trials_by_category.py` — % empty by hit / push-pull / L-C-R / reach cat, pooled and within-neuron.
- `03_mi_demo.py` — marginals + joint + MI three ways (bin/KSG/Ross) on representative cells.

If a file is missing, regenerate from the decisions above, not from memory of older versions.

## Roadmap (live checklist in `docs/roadmap.md`)

1. (done) Empty-trial diagnostic. 2. (done) Representative-neuron deep-dive.
3. (next) 1-D bridge: reproduce old endpoint-y result with count + KSG.
4. Multi-D kinematics [PC1, PC2, speed, accel], joint KSG.
5. Full temporal landscape without redundant sliding.
6. CP-vs-CB lead/lag + entropy with new estimator + proper stats.

## Conventions

- Colours: CP green `[106,176,134]/255`, CB purple `[154,123,193]/255`.
- Cell-based scripts (`# %%`) for the VS Code interactive window; small composable scripts that
  import `config` and `mi_io`. Run `01` once per data update; `02`/`03` reuse the in-kernel
  `NEURONS` or the cache.
- Validate any new estimator/loader on tiny synthetic data with a KNOWN answer first.
- Show evidence (numbers, figures), don't assert. Ask before large refactors, before changing a
  fixed decision, before deleting files.
