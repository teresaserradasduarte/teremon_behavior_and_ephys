# Roadmap (live checklist)

Status = whether the code exists and has been run; the scientific conclusions are the human's.
Design and decisions: `PLAN.md` (§5 build order, §6 decisions and open issues).

Main script: **`03`** (point-process GLM with smooth causal kernels). `01` and `02` are the baseline and
sanity check; `00`, `00b`, `02a`, `03b` are evidence, data preparation, diagnostics and explainers.

All runs so far: lag tag `lag150_win-375to125` (lags (0, 150] ms, window −375..+125 ms re max reach),
kinematics = med-lat y, reach category 1, CP and CB.

- [x] 1. `gctools` core: io, cache, preprocess, design, sim; causality and trial-boundary tests (`tests/`), evidence figures (`00`)
- [x] 2. `01` Gaussian VAR baseline — sim (no coupling, kin→neuron, neuron→kin), raw vs matched smoothing, real 22 sessions
- [x] 3. `02` Poisson GLM with discrete lag windows — sim (incl. no-coupling null at 100–150 Hz), real 22 sessions
  - [x] Poisson dispersion of all real neurons (`02a`) → keep 2.5 ms bins, no negative binomial
- [ ] 4. `03` Poisson GLM with smooth causal kernels (main)
  - [x] REACH_LONG: kinematics −3..+3 s from the session tracking, validated against reach_px (`00b`)
  - [x] sim: no coupling, kin→neuron, neuron→kin, shared drive (incl. task-mean nuisance and kinematic-history tests)
  - [x] real: FerreroRocher R1, then all 22 sessions (`real_22sessions_kernels_kinLong_kinHmax900_from-225`)
  - [x] explainer figures (`03b`)
  - [ ] remaining sim scenarios: `bidirectional`, `common_latent`
  - [ ] `04` systematic simulation sweep (FP rate, power, lag error) → check the §4 pass criteria, fix defaults
  - [ ] final real-data run with the defaults fixed by 04
- [ ] 5. **Longer lags / windows** (motivated by the MI results) — new lag tag, keep the first pass
- [ ] 6. CP vs CB summaries: fraction significant per direction, lag distributions, with / without neurons > 100 Hz
- [ ] 7. Later: CTX; CB cell types

## Open issues
- `<LAG_TAG>/01_var_continuous/` holds early 01 runs (29–30 Sep, before the folder was renamed `01_VAR`).
- 1 ms bins only if excluding neurons > 100 Hz changes a conclusion — not yet assessed.

## Literature
Organised by the step it motivates. **KEY** = the paper the step is built on; the rest support it.
References written from memory (2026-11-16): check volume / pages before citing them in a manuscript.

### Step 1 — Granger causality (`01`, Gaussian VAR baseline)
- **KEY** Granger CWJ (1969). Investigating causal relations by econometric models and cross-spectral methods. *Econometrica* 37:424–438. — the definition: the past of X improves the prediction of Y beyond Y's own past (full vs reduced model).
- **KEY** Stokes PA, Purdon PL (2017). A study of problems encountered in Granger causality analysis from a neuroscience perspective. *PNAS* 114:E7063–E7072. — the interpretation rule ("predictive precedence", not "drives"); no spectral GC; filtering distorts GC (01 raw vs matched smoothing); the own-history model must be adequate (step 3).
- Supporting: Seth AK, Barrett AB, Barnett L (2015). Granger causality analysis in neuroscience and neuroimaging. *J Neurosci* 35:3293–3297 (review, latent drivers → `common_latent`); Barnett L, Seth AK (2014). The MVGC multivariate Granger causality toolbox. *J Neurosci Methods* 223:50–68 (VAR GC, order by BIC); Geweke J (1982). *JASA* 77:304–313 (the spectral GC we do NOT use); Barnett L, Barrett AB, Seth AK (2009). Granger causality and transfer entropy are equivalent for Gaussian variables. *Phys Rev Lett* 103:238701 (bridge to the MI project).

### Step 2 — moving to point processes (`02`, Poisson GLM on raw 2.5 ms counts)
- **KEY** Truccolo W, Eden UT, Fellows MR, Donoghue JP, Brown EN (2005). A point process framework for relating neural spiking activity to spiking history, neural ensemble, and extrinsic covariate effects. *J Neurophysiol* 93:1074–1089. — THE framework: conditional intensity with the neuron's own spiking history plus extrinsic covariates (here kinematics), fitted as a Poisson GLM on fine bins; replaces the Gaussian model on smoothed rates of step 1.
- Supporting: Kim S, Putrino D, Ghosh S, Brown EN (2011). A Granger causality measure for point process models of ensemble neural spiking activity. *PLoS Comput Biol* 7:e1001110 (GC as a likelihood ratio between full and reduced point-process GLMs; rule 3); Paninski L (2004). Maximum likelihood estimation of cascade point-process neural encoding models. *Network* 15:243–262 (concave likelihood → IRLS fit); Kass RE, Eden UT, Brown EN (2014). *Analysis of Neural Data*. Springer (GLMs, dispersion → `02a`).

### Step 3 — point-process GLM with smooth temporal kernels (`03`, main analysis)
- **KEY** Park IM, Meister MLR, Huk AC, Pillow JW (2014). Encoding and decoding in parietal cortex during sensorimotor decision-making. *Nat Neurosci* 17:1395–1403. — smooth temporal kernels for every covariate, each a weighted sum of raised-cosine basis functions, in a Poisson GLM; the kernel shape β(τ) is the readout (our lag readout).
- How 03 applies it:
  - **strictly causal kernels** — every basis starts at lag ≥ 1 bin; event kernels active only after the event (CLAUDE.md rule 1).
  - **cross-kernels on log-spaced raised-cosine bases** (fine at short lags, broad at long lags) — Pillow JW, Shlens J, Paninski L, Sher A, Litke AM, Chichilnisky EJ, Simoncelli EP (2008). Spatio-temporal correlations and visual signalling in a complete neuronal population. *Nature* 454:995–999; Park et al. 2014.
  - **cross-validation in both directions** — trial-wise K-fold held-out gain of full over reduced: bits/spike for kin→neuron (Poisson), held-out variance for neuron→kin (Gaussian).
  - **own-history kernel chosen per neuron** on the reduced model, by CV, so the cross-kernel cannot absorb unmodelled history — as recommended by **Stokes & Purdon 2017** (step 1).

### Statistics used throughout (01–03)
- Besag J, Clifford P (1991). Sequential Monte Carlo p-values. *Biometrika* 78:301–304. — adaptive shuffles (`stats.sequential_p`, `runs.adaptive_extend`).
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS B* 57:289–300. — BH-FDR across neurons.
- Storey JD (2002). A direct approach to false discovery rates. *JRSS B* 64:479–498. — π₁, the fraction of non-null neurons (`stats.py`).

## Decisions log (append as they solidify)
- 2026-09-29: t = 0 = max-reach frame (python index 150); window −375..+125 ms; lags (0, 150] ms first pass;
  shuffle within spout side (LCR); CP and CB only. Details: `PLAN.md` §6.
- 2026-10-01: Poisson only (no NB), 2.5 ms bins (`02a`). Own history never shorter than the cross range (01 artifact).
- 2026-10-01: 03 uses REACH_LONG, kinematic own history 900 ms, neuron own history by CV (150/200/300 ms),
  fitted rows from −225 ms; primary effect size = excess over the shuffle null.

## Reference list (as provided, 2026-11-16)
Four entries corrected the same day, checked online: Hatsopoulos & Suminski, Garenne et al., Muzzu et al., and Gretton et al.

### Mutual Information estimation (used estimator)
- Kraskov, Stögbauer & Grassberger (2004). Phys Rev E 69:066138.

### MI in Neuroscience
- Paninski, Fellows, Hatsopoulos & Donoghue (2004). J Neurophysiol 91:515–532. Spatiotemporal tuning of motor cortical neurons for hand position and velocity.
- Paninski, Shoham, Fellows, Hatsopoulos & Donoghue (2004). J Neurosci 24:8551–8561. Superlinear population encoding of dynamic hand trajectory in M1.
- Hatsopoulos & Suminski (2011). Neuron 72(3):477–487. Sensing with the motor cortex (review). https://doi.org/10.1016/j.neuron.2011.10.020
- Garenne, Pasquereau, Guthrie, Bioulac & Boraud (2011). Frontiers in Systems Neuroscience 5:23. Basal ganglia preferentially encode context dependent choice in a two-armed bandit task. https://doi.org/10.3389/fnsys.2011.00023
- Muzzu, Mitolo, Gava & Schultz (2018). PLOS ONE 13(9):e0203900. Encoding of locomotion kinematics in the mouse cerebellum. https://doi.org/10.1371/journal.pone.0203900
- Yamamoto, Kawato et al. (2007). J Neurophysiol / PMID 17079350. Encoding of movement dynamics by Purkinje cell simple spikes under resistive and assistive force fields.
- Popa, Streng & Ebner (2017). eNeuro 4(2). Long-term predictive and feedback encoding in Purkinje cell simple spikes.

### Dependence measures
- Rényi, A. (1959). On measures of dependence. Acta Mathematica Hungarica, 10(3–4), 441–451. https://doi.org/10.1007/bf02024507
- Gretton, Smola, Bousquet, Herbrich, Belitski, Augath, Murayama, Pauls, Schölkopf & Logothetis (2005). Kernel constrained covariance for dependence measurement. Proceedings of the Tenth International Workshop on Artificial Intelligence and Statistics (AISTATS), PMLR R5:112–119. http://proceedings.mlr.press/r5/gretton05a/gretton05a.pdf
- Yamada, M., Jitkrittum, W., Sigal, L., Xing, E. P., & Sugiyama, M. (2014). High-dimensional feature selection by feature-wise kernelized Lasso. Neural Computation, 26(1), 185–207. https://doi.org/10.1162/NECO_a_00537
- Stokes, P. A., & Purdon, P. L. (2017). A study of problems encountered in Granger causality analysis from a neuroscience perspective. Proceedings of the National Academy of Sciences of the United States of America, 114(34), E7063–E7072. https://doi.org/10.1073/pnas.1704663114

### Memming's old work
- Poisson extension of autoregressive modeling: Park, I. M., Meister, M. L. R., Huk, A. C. and Pillow, J. W. (2014). Encoding and decoding in parietal cortex during sensorimotor decision-making. Nature Neuroscience, 17(10), 1395-1403. https://doi.org/10.1038/nn.3800
- Entropy estimator: Archer, E., Park, I. M. and Pillow, J. (2014). Bayesian Entropy Estimation for Countable Discrete Distributions. Journal of Machine Learning Research, 15, 2833-2868. http://jmlr.org/papers/v15/archer14a.html
- Mutual information estimator: Archer, E., Park, I. M. and Pillow, J. (2013). Bayesian and Quasi-Bayesian Estimators for Mutual Information from Discrete Data. Entropy, 15(5), 1738-1755. https://doi.org/10.3390/e15051738
- Spike train kernels for low-firing rate neurons: Park, I. M., Seth, S., Rao, M. and Príncipe, J. C. (2012). Strictly positive definite spike train kernels for point process divergences. Neural Computation, 24(8), 2223-2250. https://doi.org/10.1162/NECO_a_00309

### On bivariate VAR
- Granger, C.W.J. (1969). Investigating Causal Relations by Econometric Models and Cross-spectral Methods. Econometrica — the original bivariate VAR/Granger framework.
- Geweke, J. (1982). Measurement of Linear Dependence and Feedback Between Multiple Time Series. JASA — the frequency-domain decomposition.
- Bressler, S.L. & Seth, A.K. (2011). Wiener–Granger Causality: A well established methodology. NeuroImage — a good practical review for neuroscience applications, including the stationarity caveats.
- Truccolo, W., Eden, U.T., Fellows, M.R., Donoghue, J.P. & Brown, E.N. (2005). A point process framework for relating neural spiking activity to spiking history, neural ensemble, and extrinsic covariate effects. J. Neurophysiology — the foundational point-process GLM paper, exactly built for spike-count-vs-covariate autoregression.
- Kim, S., Putrino, D., Ghosh, S. & Brown, E.N. (2011). A Granger causality measure for point process models of ensemble neural spiking activity. PLoS Computational Biology — point-process Granger causality.
