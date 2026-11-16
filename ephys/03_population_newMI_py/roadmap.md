# Roadmap (live checklist)

Status = whether the code exists and has been run; the scientific conclusions are the human's.

Core scripts: `04`, `09`, `12`, `15`, `18`, `28` — most important: **`12`**.

- [x] 1. Empty-trial diagnostic — how many, and is emptiness structured? (`02`)
- [x] 2. Representative-neuron deep-dive — raw → discrete variable → what MI sees (`03`, `14` for the old binning method)
- [x] 3. 1-D bridge: MI(endpoint y ; spike count) per sliding window, binning / KSG / Ross (`04`);
      neural-variable modes count / isi / conv compared (`06`, `07`); per mouse (`08`); curves with raw data (`17`)
- [x] 4. Multi-dimensional kinematics, joint multivariate KSG (`09`). Feature set chosen by excess MI in `11`
      (now x, y, z, speed, dist_water; water features from `10`), replacing the planned [PC1, PC2, speed, accel]
- [x] 5. Temporal landscapes
  - [x] neural window start × stop, kinematics fixed (`05`, GPU)
  - [x] TIME × TIME maps: kinematic window × neural window, t = 0 = reach endpoint (`12`), reach onset (`15`),
        push/pull (`18`, kinematics from `16`, cleaned and force-onset alignment from `21`)
  - [x] neural window width measured instead of chosen by eye (`20`)
- [ ] 6. CP-vs-CB lead/lag, with the new estimator + proper statistics
  - [x] h(lag) per mouse, diagonal (lag) vs neural-time R² per neuron, per-mouse maps (`12`, `mi2d_figures.py`)
  - [x] standalone lead/lag figure (`13`)
  - [x] CB maps split by cell type (`24`)
  - [ ] **in progress — what single-neuron map shapes make up the mean bands** (CP vertical band, CB diagonal band):
        split-half profiles (`22`, too noisy) → hot cluster (`23`) → template bank (`25`) → hand labels
        (`26` selection, `27` GUI) → classifier trained on the labels (`28`)
  - [ ] entropy comparison with the new estimator

## Open issues
- `mi2d_figures.py` holds 12's cells 4b–7 as functions and is used by 18; 12 still runs its own copy.
  Until 12 is switched over, a change in one must be made in the other.
- 09 / 12 use frame 147 as the kinematic zero; the true max-reach frame (st_reach t = 0) is index 150,
  10 ms later (found in `15`). 15 avoids it by aligning both axes per trial to onset.
- 09 still uses a 100 ms neural window; the chosen width is 70 ms (`20`), as in 12 / 18 / 22.

## Decisions log (append as they solidify)
- Endpoint kinematic for the 1-D bridge: mean med-lateral y over frames 147–167.
- MI reported in bits; primary estimator KSG (both continuous); Ross as sparse-cell cross-check.
- Analysis window: spike count in a sliding **70 ms** window, step 25 ms. Width chosen from `20`
  (knee = smallest width reaching KNEE_FRAC of the plateau of excess MI).
- Significance / comparisons use EXCESS MI = observed − mean of trial-shuffle nulls (needed to compare
  feature sets of different dimension, `11`). `_validate_contrast.py`: statistics measured against the
  shuffle null separate real from shuffled maps; within-map contrasts do not.
- Multi-D kinematic features: x, y, z, speed, dist_water (`11`); speed by Savitzky-Golay derivative
  (11 frames, order 3) on mm traces (`config.PX2MM`).
- TIME × TIME maps (12): kinematic windows 21 frames (52.5 ms), step 4 frames (10 ms), range −0.34..0.10 s
  re endpoint; neural windows 70 ms, step 25 ms, centres −0.5..0.5 s; 20 shuffles per cell.
- 2026-10-12: empty trials in v7.3 .mat were read as spikes at 0 s and 1 s (MATLAB_empty flag ignored);
  fixed in `mi_io`, documented in `19`. Caches / results made before the fix may need recomputing.
- Map shapes: automatic rules (23 hot cluster, 25 templates) mislabel too many maps → hand labels are the
  ground truth (classes: column, diagonal, tiling, 2 tiles, broad, nothing, mixed/unsure), classifier in `28`.
