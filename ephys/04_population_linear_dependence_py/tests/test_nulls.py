"""Trial-shuffle permutations: within stratum, never self, reproducible."""
from _harness import run
import numpy as np
from gctools import nulls


def test_derangement_within_strata():
    cond = np.array(["L"] * 7 + ["C"] * 2 + ["R"] * 11 + ["X"])      # "X": singleton
    P = nulls.make_perms(cond, 500, seed=3)
    idx = np.arange(cond.size)
    single = cond == "X"
    assert np.all(cond[P] == cond[None, :]), "a partner came from another stratum"
    assert not np.any(P[:, ~single] == idx[~single]), "a trial was paired with itself"
    assert np.all(P[:, single] == idx[single])
    assert all(np.array_equal(np.sort(p), idx) for p in P), "not a permutation"
    # every trial of a >=3 stratum gets several different partners
    n_part = [np.unique(P[:, i]).size for i in np.flatnonzero(cond == "R")]
    print(f"  distinct partners per R trial over 500 shuffles: min {min(n_part)}")
    assert min(n_part) >= 5


def test_reproducible():
    cond = np.repeat([1, 2, 3], 20)
    assert np.array_equal(nulls.make_perms(cond, 50, 7), nulls.make_perms(cond, 50, 7))
    assert not np.array_equal(nulls.make_perms(cond, 50, 7), nulls.make_perms(cond, 50, 8))


def test_longer_sequence_extends_shorter():
    cond = np.repeat([1, 2, 3], 15)
    a = nulls.make_perms(cond, 300, 5)
    assert np.array_equal(a[:100], nulls.make_perms(cond, 100, 5))
    assert np.array_equal(nulls.perm_batch(cond, 100, 300, 5), a[100:300])


def test_sequential_p_calibrated():
    """20,000 null 'neurons' whose observed value is one more draw from the same
    distribution: P(p <= t) <= t within Monte Carlo error, down to small t."""
    from gctools import stats
    rng = np.random.default_rng(0)
    n0, cap, n = 200, 5000, 20000
    ps = np.empty(n)
    for i in range(n):
        obs = rng.normal()
        first = rng.normal(size=n0)
        g0 = int(np.sum(first >= obs))
        ps[i] = stats.sequential_p(obs, g0, n0, lambda a, b: rng.normal(size=b - a), cap, 10, 1000)[0]
    for t in (0.05, 0.01, 0.001, 0.0002):
        f = np.mean(ps <= t)
        se = np.sqrt(t * (1 - t) / n)
        print(f"  P(p <= {t:g}) = {f:.5f}  (nominal {t:g}, +3 SE {t + 3*se:.5f})")
        assert f <= t + 3 * se


def test_sequential_p_resolution():
    from gctools import stats
    rng = np.random.default_rng(1)
    more = lambda a, b: rng.normal(size=b - a)
    p, n, m = stats.sequential_p(50.0, 0, 1000, more, 24000, 10, 1000)   # never reached
    assert m == "cap" and n == 24000 and np.isclose(p, 1 / 24001)
    p, n, m = stats.sequential_p(3.3, 0, 1000, more, 24000, 10, 1000)    # ~1 in 2000 above
    print(f"  obs 3.3 sd: p = {p:.2e} after {n} shuffles ({m}); true tail 4.8e-04")
    assert m == "sequential" and p == 10 / n and 1e-4 < p < 2e-3
    assert stats.sequential_p(0.0, 3, 1000, more, 24000)[2] == "fixed"   # not at floor


if __name__ == "__main__":
    run(globals())
