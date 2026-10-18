"""stats.py - shuffle p-values, null-standardised effect sizes, BH-FDR.

Rule 5: effect sizes are compared across directions only after standardising each
direction against ITS OWN shuffle null (null_z) or via shuffle p-values -- never
as raw deviances / variance ratios.
"""
import numpy as np
from statsmodels.stats.multitest import multipletests


def shuffle_p(obs, null):
    """One-sided permutation p = (1 + #{null >= obs}) / (1 + n). NaN if obs is NaN."""
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size == 0:
        return np.nan
    return (1.0 + np.sum(null >= obs)) / (1.0 + null.size)


def null_z(obs, null):
    """(obs - mean(null)) / sd(null): effect in units of the direction's own null."""
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size < 2 or np.std(null) == 0:
        return np.nan
    return (obs - null.mean()) / null.std(ddof=1)


def bh_fdr(p, q=0.05):
    """Benjamini-Hochberg across the finite entries of p -> (reject bool, p_adj).
    NaN p -> not rejected, p_adj NaN."""
    p = np.asarray(p, float)
    rej = np.zeros(p.shape, bool)
    padj = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        r, pa, _, _ = multipletests(p[ok], alpha=q, method="fdr_bh")
        rej[ok], padj[ok] = r, pa
    return rej, padj


def sequential_p(obs, g0, n0, more_fn, cap, h=10, batch=1000):
    """Adaptive (sequential) permutation p-value for ONE neuron (Besag & Clifford 1991).
      g0, n0   exceedances (null >= obs) among the first n0 shuffles
      more_fn  more_fn(a, b) -> null statistics of shuffles a..b-1 (same sequence)
    Only a neuron that beat every one of the first n0 shuffles (g0 == 0) is extended, in
    batches, until h shuffles reach it (p = h / L at the L-th shuffle) or the cap is hit
    (p = (g + 1) / (cap + 1)). Otherwise p = (g0 + 1) / (n0 + 1) as usual.
    -> (p, n_used, method in {"fixed", "sequential", "cap"})"""
    if g0 > 0 or n0 >= cap:
        return (g0 + 1.0) / (n0 + 1.0), n0, "fixed"
    g, n = 0, n0
    while n < cap:
        k = min(batch, cap - n)
        exceed = np.asarray(more_fn(n, n + k)) >= obs
        cum = g + np.cumsum(exceed)
        if cum.size and cum[-1] >= h:
            L = n + int(np.argmax(cum >= h)) + 1
            return h / L, L, "sequential"
        g, n = int(cum[-1]) if cum.size else g, n + k
    return (g + 1.0) / (cap + 1.0), cap, "cap"


def storey_pi1(p, lam=0.5, n_boot=1000, seed=0, groups=None):
    """Estimated FRACTION OF NEURONS WITH AN EFFECT from the whole p-value distribution
    (Storey 2002): null p-values are uniform, so the density above lam is almost all null:
        pi0 = min(1, #{p > lam} / ((1 - lam) n)),   pi1 = 1 - pi0.
    Counts effects too weak to pass FDR one by one. 95% CI by bootstrap:
      groups=None  -> resample neurons (assumes independent p-values)
      groups=array -> resample whole GROUPS (sessions): neurons of a session share the
                      kinematics and the shuffles, so their p-values are correlated and a
                      neuron-level bootstrap is too narrow. With < 3 groups the CI is NaN.
    -> (pi1, lo, hi); NaN if no finite p."""
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    p = p[ok]
    if p.size == 0:
        return np.nan, np.nan, np.nan
    f = lambda q: 1.0 - min(1.0, np.mean(q > lam) / (1.0 - lam)) if q.size else np.nan
    rng = np.random.default_rng(seed)
    if not n_boot:
        return f(p), np.nan, np.nan
    if groups is None:
        boot = [f(p[rng.integers(0, p.size, p.size)]) for _ in range(n_boot)]
    else:
        g = np.asarray(groups)[ok]
        ug = np.unique(g)
        if ug.size < 3:
            return f(p), np.nan, np.nan
        members = [p[g == u] for u in ug]
        boot = [f(np.concatenate([members[k] for k in rng.integers(0, ug.size, ug.size)]))
                for _ in range(n_boot)]
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    return f(p), float(lo), float(hi)


def binom_ci(k, n, level=0.95):
    """Clopper-Pearson interval for a rate k/n."""
    from scipy.stats import beta
    a = (1 - level) / 2
    lo = beta.ppf(a, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(1 - a, k + 1, n - k) if k < n else 1.0
    return lo, hi
