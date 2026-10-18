"""Population effect measures: Storey pi1 (fraction of neurons with an effect) and the
excess (shuffle-corrected) effect sizes."""
from _harness import run
import numpy as np
from gctools import stats, models


def test_storey_pi1_uniform_and_mixture():
    rng = np.random.default_rng(0)
    p0 = rng.uniform(size=2000)
    pi, lo, hi = stats.storey_pi1(p0)
    print(f"  uniform p: pi1 = {pi:.3f} [{lo:.3f}, {hi:.3f}] (true 0)")
    assert pi < 0.05 and lo < 0.01                       # pi1 >= 0 by construction
    p = np.r_[rng.uniform(size=1600), rng.beta(0.2, 20, size=400)]      # 20 % with an effect
    pi, lo, hi = stats.storey_pi1(p)
    print(f"  20% mixture: pi1 = {pi:.3f} [{lo:.3f}, {hi:.3f}] (true 0.20)")
    assert lo <= 0.20 <= hi
    # discrete shuffle p-values (floor 1/201) behave the same
    pd_ = np.maximum(np.round(p * 201) / 201, 1 / 201)
    assert abs(stats.storey_pi1(pd_)[0] - pi) < 0.03


def test_storey_session_bootstrap_honest():
    """p-values correlated within sessions (shared shift, all null in truth): the
    session-level CI is wider than the neuron-level one and covers pi1 = 0 far more often."""
    from scipy.stats import norm
    rng = np.random.default_rng(3)
    cover_n, cover_s, w_n, w_s = 0, 0, [], []
    for rep in range(60):
        g = np.repeat(np.arange(20), 40)                       # 20 sessions x 40 neurons
        shift = rng.normal(0, 0.5, 20)[g]                      # session-shared nuisance
        z = (rng.normal(size=g.size) + shift) / np.sqrt(1.25)  # marginally N(0,1): null
        p = norm.sf(z)
        _, lo1, hi1 = stats.storey_pi1(p, n_boot=300, seed=rep)
        _, lo2, hi2 = stats.storey_pi1(p, n_boot=300, seed=rep, groups=g)
        cover_n += lo1 <= 0 <= hi1 + 1e-12; cover_s += lo2 <= 0 <= hi2 + 1e-12
        w_n.append(hi1 - lo1); w_s.append(hi2 - lo2)
    print(f"  CI covers true pi1=0: neuron-bootstrap {cover_n}/60, session-bootstrap {cover_s}/60; "
          f"mean width {np.mean(w_n):.3f} vs {np.mean(w_s):.3f}")
    assert np.mean(w_s) > 1.5 * np.mean(w_n) and cover_s >= 54
    assert np.isnan(stats.storey_pi1(p, groups=np.zeros(p.size))[1])   # < 3 sessions -> no CI


def test_excess_effects_centred_under_null():
    """If the observed GC is just another draw from its null, the excess partial R^2
    averages ~0 across neurons (and is clearly > 0 when it is shifted up)."""
    rng = np.random.default_rng(1)
    ex0, ex1 = [], []
    for _ in range(500):
        null = rng.chisquare(15, size=200) / 12000          # GC-like null (15 params, 12k rows)
        ex0.append(models.excess_partial_r2(rng.chisquare(15) / 12000, null))
        ex1.append(models.excess_partial_r2(rng.chisquare(15) / 12000 + 0.01, null))
    print(f"  excess partial R2 (%): null mean {np.mean(ex0):.4f}, shifted mean {np.mean(ex1):.3f}")
    assert abs(np.mean(ex0)) < 0.005 and np.mean(ex1) > 0.9
    b = [models.excess_bits_per_spike(rng.chisquare(5), rng.chisquare(5, 200), 800) for _ in range(500)]
    assert abs(np.mean(b)) < 0.002


if __name__ == "__main__":
    run(globals())
