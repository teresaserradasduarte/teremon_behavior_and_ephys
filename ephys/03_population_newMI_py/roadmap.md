# Roadmap (live checklist)

- [x] 1. Empty-trial diagnostic — how many, and is emptiness structured? (`02_empty_trials_by_category.py`)
- [x] 2. Representative-neuron deep-dive — raw → discrete variable → what MI sees (`03_mi_demo.py`)
- [ ] 3. 1-D bridge: reproduce the old endpoint-y result with the new count variable + KSG
- [ ] 4. Multi-dimensional kinematics [PC1, PC2, speed, accel], joint KSG
- [ ] 5. Full temporal landscape (sweep neuron window start/stop) without redundant sliding
- [ ] 6. CP-vs-CB lead/lag and entropy comparisons with the new estimator + proper statistics

## Decisions log (append as they solidify)
- Endpoint kinematic for the 1-D bridge: mean med-lateral y over frames 147–167.
- MI reported in bits; primary estimator KSG (both continuous); Ross as sparse-cell cross-check.
- Analysis window [start, stop]: TBD from real empty-trial numbers (run 02 first).
