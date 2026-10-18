"""VAR machinery (script 01): agreement with statsmodels, recovery of a known
Gaussian VAR, and F-test calibration where its assumptions hold."""
from _harness import run
import numpy as np
import statsmodels.api as sm
from gctools import models, nulls


def _var2(n_tr, n_b, b_xy, rng, b_yx=0.0):
    """Trials of a stationary Gaussian VAR(2): x AR(2); y AR(2) + b_xy * x[t-2]
    (+ b_yx * y[t-1] into x). Returns y, x (n_tr, n_b) after a 100-sample burn-in."""
    T = n_b + 100
    y = np.zeros((n_tr, T)); x = np.zeros((n_tr, T))
    ex, ey = rng.normal(size=(2, n_tr, T))
    for t in range(2, T):
        x[:, t] = 0.5 * x[:, t-1] - 0.2 * x[:, t-2] + b_yx * y[:, t-1] + ex[:, t]
        y[:, t] = 0.4 * y[:, t-1] - 0.3 * y[:, t-2] + b_xy * x[:, t-2] + ey[:, t]
    return y[:, 100:], x[:, 100:]


def test_matches_statsmodels():
    rng = np.random.default_rng(0)
    y, x = _var2(40, 80, 0.3, rng)
    inw = np.ones(80, bool)
    res = models.var_neuron(y, x, inw, p_max=6, p_fixed=3)
    Ly, Lx, rows = models.var_lags(y, x, inw, 6)
    yt = y[:, rows].ravel(); p = 3
    Yl = Ly[:, rows, :p].reshape(-1, p); Xl = Lx[:, rows, :p].reshape(-1, p)
    full = sm.OLS(yt, sm.add_constant(np.column_stack([Yl, Xl]))).fit()
    red = sm.OLS(yt, sm.add_constant(Yl)).fit()
    F, pF, _ = full.compare_f_test(red)
    gc = np.log(red.ssr / full.ssr)
    print(f"  GC ours {res['xy']['gc']:.10f} vs statsmodels {gc:.10f}; "
          f"F {res['xy']['F']:.6f} vs {F:.6f}")
    assert np.isclose(res["xy"]["gc"], gc, rtol=0, atol=1e-10)
    assert np.isclose(res["xy"]["F"], F, rtol=1e-9)
    assert np.isclose(res["xy"]["p_F"], pF, rtol=1e-8, atol=1e-300)
    assert rows[:6].sum() == 0 and rows[6:].all(), "rows must start after p_max lags"


def test_recovers_direction_and_order():
    rng = np.random.default_rng(1)
    ps, sig_xy, sig_yx = [], 0, 0
    for _ in range(20):
        y, x = _var2(30, 100, 0.4, rng)
        r = models.var_neuron(y, x, np.ones(100, bool), p_max=8)
        ps.append(r["p"]); sig_xy += r["xy"]["p_F"] < 0.05; sig_yx += r["yx"]["p_F"] < 0.05
    print(f"  BIC p: {np.bincount(ps).tolist()}  | x->y detected {sig_xy}/20, "
          f"y->x (null) {sig_yx}/20")
    assert np.mean(np.array(ps) == 2) >= 0.7
    assert sig_xy == 20 and sig_yx <= 3


def test_ftest_and_shuffle_calibrated_under_null():
    rng = np.random.default_rng(2)
    reps, pF, pS = 300, [], []
    for i in range(reps):
        y, x = _var2(20, 60, 0.0, rng)
        cond = np.repeat([1, 2], 10)
        perms = nulls.make_perms(cond, 99 if i < 200 else 0, seed=i)
        r = models.var_neuron(y, x, np.ones(60, bool), p_max=4, p_fixed=2,
                              perms=perms if len(perms) else None)
        pF.append(r["xy"]["p_F"]); pS.append(r["xy"]["p_shuf"])
    fF = np.mean(np.array(pF) < 0.05)
    pS = np.array(pS)[:200]; fS = np.mean(pS < 0.05)
    print(f"  null FPR at 0.05: F-test {fF:.3f} (n={reps}), shuffle {fS:.3f} (n=200)")
    assert 0.02 <= fF <= 0.09          # binomial 99% CI around 0.05 for n=300
    assert fS <= 0.095                 # n=200: 99% CI upper bound for 0.05


def test_varnull_reproduces_stored_null():
    """Adaptive shuffles use VarNull on the continued sequence; on the first shuffles it
    must reproduce var_neuron's stored null bit-for-bit (both directions, BIC and fixed p)."""
    rng = np.random.default_rng(5)
    y, x = _var2(30, 70, 0.3, rng)
    inw = np.ones(70, bool)
    cond = np.repeat([1, 2, 3], 10)
    perms = nulls.make_perms(cond, 40, seed=2)
    for p_fixed in (None, 6):
        r = models.var_neuron(y, x, inw, 6, perms[:25], p_fixed)
        for d in ("xy", "yx"):
            vn = models.VarNull(y, x, inw, 6, r["p"], d)
            assert np.array_equal(vn.stats(perms[:25]), r[d]["null"]), (p_fixed, d)
            assert np.array_equal(vn.stats(perms[25:40]),
                                  models.var_neuron(y, x, inw, 6, perms, p_fixed)[d]["null"][25:40])


def test_lag_readout_recovers_known_lag():
    """VAR(2) with x->y only at lag 2 (x[t-2]): basis kernel peak at 2 bins +-1, and the
    bootstrap CI covers it in most repetitions."""
    from gctools import design as ds
    rng = np.random.default_rng(3)
    B, lags, _ = ds.raised_cosine_basis(6, (0, 80), 10.0, "log")
    peaks, cover = [], 0
    for i in range(10):
        y, x = _var2(40, 80, 0.4, rng)
        lk = models.var_lag_kernel(y, x, np.ones(80, bool), 8, B, lags, 10.0, n_boot=50, seed=i)
        pk, ci = lk["xy"]["peak_ms"], lk["xy"]["peak_ci"]
        peaks.append(pk); cover += ci[0] <= 20.0 <= ci[1]
    print(f"  peak lags (true 20 ms): {peaks}; bootstrap CI covers truth {cover}/10")
    assert np.all(np.abs(np.array(peaks) - 20.0) <= 10.0)
    assert cover >= 7


if __name__ == "__main__":
    run(globals())
