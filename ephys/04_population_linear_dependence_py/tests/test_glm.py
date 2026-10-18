"""Poisson GLM machinery (script 02): agreement with statsmodels, and calibration of
the shuffle / chi2 tests under the null."""
from _harness import run
import numpy as np
import statsmodels.api as sm
from gctools import models, nulls, design as ds


def _poisson_data(rng, n=5000, k=4, kc=3, bc=0.0):
    Xr = np.column_stack([np.ones(n), rng.normal(size=(n, k - 1)) * 0.3])
    Xc = rng.normal(size=(n, kc)) * 0.3
    b = np.r_[np.log(0.2), rng.normal(size=k - 1) * 0.5]
    y = rng.poisson(np.exp(Xr @ b + Xc @ (np.full(kc, bc))))
    return Xr, Xc, y


def test_irls_matches_statsmodels():
    rng = np.random.default_rng(0)
    Xr, Xc, y = _poisson_data(rng, bc=0.3)
    X = np.hstack([Xr, Xc])
    ours = models.poisson_irls(X, y)
    ref = sm.GLM(y, X, family=sm.families.Poisson()).fit()
    print(f"  max |b - b_sm| = {np.max(np.abs(ours['b'] - ref.params)):.2e}; "
          f"deviance {ours['deviance']:.6f} vs {ref.deviance:.6f}")
    assert np.allclose(ours["b"], ref.params, atol=1e-6)
    assert np.isclose(ours["deviance"], ref.deviance, rtol=1e-9)
    assert np.allclose(np.sqrt(np.diag(ours["cov"])), ref.bse, rtol=1e-5)


def test_score_matches_statsmodels():
    rng = np.random.default_rng(1)
    Xr, Xc, y = _poisson_data(rng, bc=0.1)
    sc = models.PoissonScore(Xr, y)
    ref = sm.GLM(y, Xr, family=sm.families.Poisson()).fit()
    chi2_sm = np.asarray(ref.score_test(exog_extra=Xc)[0]).ravel()[0]
    print(f"  score statistic ours {sc.stat(Xc):.6f} vs statsmodels {float(chi2_sm):.6f}")
    assert np.isclose(sc.stat(Xc), float(chi2_sm), rtol=1e-6)


def test_ridge_objective_decreases_and_shrinks():
    rng = np.random.default_rng(2)
    Xr, Xc, y = _poisson_data(rng, bc=0.3)
    X = np.hstack([Xr, Xc])
    b0 = models.poisson_irls(X, y)["b"]; b1 = models.poisson_irls(X, y, ridge=500.0)["b"]
    assert np.all(np.abs(b1[1:]) <= np.abs(b0[1:]) + 1e-9), "ridge did not shrink"


def test_glm_neuron_null_calibrated():
    """Independent spikes and 'kinematics' on trials: shuffle p and chi2 p calibrated for
    x->y; shuffle p calibrated for y->x (smooth AR kinematics)."""
    rng = np.random.default_rng(3)
    hist = ds.ms_windows_to_lags([(0, 2.5), (2.5, 5), (5, 10), (10, 20), (20, 40)], 2.5)
    cross = ds.ms_windows_to_lags([(0, 10), (10, 20), (20, 40)], 2.5)
    centers = np.array([5, 15, 30.])
    inw = np.arange(80) >= 20
    cond = np.repeat([1, 2], 20)
    pS, pC, pSx = [], [], []
    for i in range(150):
        yc = rng.poisson(0.05, size=(40, 80)).astype(float)
        x = np.cumsum(rng.normal(size=(40, 80)), axis=1) * 0.1
        perms = nulls.make_perms(cond, 99, seed=i)
        r = models.glm_neuron(yc, x, inw, hist, cross, centers, perms)
        pS.append(r["xy"]["p_shuf"]); pC.append(r["xy"]["p_chi2"]); pSx.append(r["yx"]["p_shuf"])
    fS, fC, fX = (np.mean(np.array(v) < 0.05) for v in (pS, pC, pSx))
    print(f"  null FPR at 0.05 (n=150): x->y shuffle {fS:.3f}, x->y chi2 {fC:.3f}, y->x shuffle {fX:.3f}")
    assert fS <= 0.11 and fX <= 0.11 and fC <= 0.12


def test_glm_neuron_detects_coupling_window():
    """Spikes driven by x at lag (10, 20] ms: x->y detected, best window = (10, 20]."""
    rng = np.random.default_rng(4)
    hist = ds.ms_windows_to_lags([(0, 2.5), (2.5, 5), (5, 10), (10, 20), (20, 40)], 2.5)
    cross = ds.ms_windows_to_lags([(0, 10), (10, 20), (20, 40)], 2.5)
    centers = np.array([5, 15, 30.])
    x = rng.normal(size=(60, 100))
    drive = np.zeros_like(x); drive[:, 6:] = x[:, :-6]           # lag 6 bins = 15 ms
    yc = rng.poisson(np.exp(np.log(0.1) + 0.8 * drive)).astype(float)
    perms = nulls.make_perms(np.repeat([1, 2, 3], 20), 99, seed=0)
    r = models.glm_neuron(yc, x, np.arange(100) >= 20, hist, cross, centers, perms)
    print(f"  x->y p_shuf {r['xy']['p_shuf']:.3f}, best window {r['xy']['lag']['best_ms']} ms, "
          f"bits/spike {r['xy']['bits_per_spike']:.3f}; y->x p_shuf {r['yx']['p_shuf']:.2f}")
    assert r["xy"]["p_shuf"] <= 0.01 and r["xy"]["lag"]["best_ms"] == 15.0


def test_glmnull_reproduces_stored_null():
    rng = np.random.default_rng(6)
    hist = ds.ms_windows_to_lags([(0, 2.5), (2.5, 5), (5, 10), (10, 20)], 2.5)
    cross = ds.ms_windows_to_lags([(0, 10), (10, 20)], 2.5)
    x = rng.normal(size=(30, 60))
    yc = rng.poisson(0.2, size=(30, 60)).astype(float)
    inw = np.arange(60) >= 10
    perms = nulls.make_perms(np.repeat([1, 2], 15), 20, seed=3)
    r = models.glm_neuron(yc, x, inw, hist, cross, np.array([5, 15.]), perms, 0.1)
    for d in ("xy", "yx"):
        g = models.GlmNull(yc, x, inw, hist, cross, 0.1, d)
        assert np.array_equal(g.stats(perms), r[d]["null"]), d


if __name__ == "__main__":
    run(globals())
