"""models.py - model fitting.

M2: bivariate Gaussian VAR on continuous signals (script 01, didactic baseline).
    y equation:  y_t = c + sum_{j=1..p} a_j y_{t-j} + sum_{j=1..p} b_j x_{t-j} + e
    x equation:  the same with x and y swapped
  GC (x->y) = ln(RSS_reduced / RSS_full) of the y equation (log variance ratio);
  F-test p reported as SECONDARY; significance comes from the trial-shuffle null.

ROWS (rule 3). One set of rows per neuron: bins where every trial has y, x and all
P_MAX lags of both finite, inside the analysis window. Every order p candidate, both
equations, full and reduced models, and every shuffle use exactly these rows. The
row set is the SAME bins in every trial, so a trial-shuffle only swaps which trial
the cross-lags come from.

Fitting uses the reduced model's QR: the full-model RSS is the reduced residual
minus its projection on the cross-lag columns residualised on the reduced design
(Frisch-Waugh-Lovell) -- identical to a full lstsq fit, but the reduced part is
computed once and reused by every shuffle. Checked against statsmodels in tests.
"""
import numpy as np
from scipy.stats import f as _f
from .design import lagged_block, lag_matrix
from . import stats


# ---------------------------------------------------------------------------
# OLS pieces
# ---------------------------------------------------------------------------
def ols_rss(Z, t):
    """RSS of the least-squares fit of t on Z (Z should contain the constant)."""
    b, *_ = np.linalg.lstsq(Z, t, rcond=None)
    r = t - Z @ b
    return float(r @ r)


def _reduced(t, O):
    """QR of the reduced design [1, O] and the reduced residual of t."""
    Z = np.column_stack([np.ones(len(t)), O])
    Q, _ = np.linalg.qr(Z)
    r = t - Q @ (Q.T @ t)
    return Q, r


def _full_rss(Q, r, C):
    """RSS after adding columns C to the reduced design (FWL)."""
    Ct = C - Q @ (Q.T @ C)
    b, *_ = np.linalg.lstsq(Ct, r, rcond=None)
    e = r - Ct @ b
    return float(e @ e)


def excess_partial_r2(gc, null_gc):
    """Effect STRENGTH, in % of the variance left after own history: partial R^2 =
    1 - RSS_full/RSS_red = 1 - exp(-GC), minus its mean over the neuron's own shuffles
    (removes the positive bias of the extra cross parameters; ~0 without an effect)."""
    null_gc = np.asarray(null_gc, float)
    null_gc = null_gc[np.isfinite(null_gc)]
    if not np.isfinite(gc) or null_gc.size == 0:
        return np.nan
    return 100.0 * ((1 - np.exp(-gc)) - np.mean(1 - np.exp(-null_gc)))


def excess_bits_per_spike(lr, null_score, n_spikes):
    """Poisson effect STRENGTH: (LR - mean shuffle score statistic) / (2 n_spikes ln 2).
    The score statistic and the LR agree for large samples, so the shuffle mean
    estimates the chance-level LR. ~0 without an effect."""
    null_score = np.asarray(null_score, float)
    null_score = null_score[np.isfinite(null_score)]
    if not np.isfinite(lr) or null_score.size == 0 or n_spikes <= 0:
        return np.nan
    return (lr - null_score.mean()) / (2 * n_spikes * np.log(2))


def gc_ftest(rss_r, rss_f, n, k_full, q):
    """GC = ln(rss_r / rss_f) and the nested F-test for q added columns."""
    gc = np.log(rss_r / rss_f)
    df2 = n - k_full
    F = ((rss_r - rss_f) / q) / (rss_f / df2)
    return gc, F, float(_f.sf(F, q, df2)), df2


# ---------------------------------------------------------------------------
# VAR
# ---------------------------------------------------------------------------
def var_lags(y, x, in_window, p_max):
    """Per-trial lag blocks (lags 1..p_max) and the common row mask.
    y, x (n_trials, n_bins) -> Ly, Lx (n_trials, n_bins, p_max), rows (n_bins,) bool."""
    lags = np.arange(1, p_max + 1)
    Ly = lagged_block(y, lag_matrix, lags)
    Lx = lagged_block(x, lag_matrix, lags)
    fin = (np.isfinite(y) & np.isfinite(x)
           & np.isfinite(Ly).all(axis=2) & np.isfinite(Lx).all(axis=2))
    rows = np.asarray(in_window, bool) & fin.all(axis=0)
    return Ly, Lx, rows


def select_order_bic(yt, xt, Yl, Xl, p_max):
    """BIC of the bivariate system over p = 1..p_max on fixed rows:
    ln det(Sigma_p) + (4p+2) ln N / N.  -> (p_best, bic (p_max,))."""
    N = len(yt)
    bic = np.full(p_max, np.nan)
    for p in range(1, p_max + 1):
        Z = np.column_stack([np.ones(N), Yl[:, :p], Xl[:, :p]])
        B, *_ = np.linalg.lstsq(Z, np.column_stack([yt, xt]), rcond=None)
        E = np.column_stack([yt, xt]) - Z @ B
        sign, logdet = np.linalg.slogdet(E.T @ E / N)
        bic[p - 1] = logdet + (4 * p + 2) * np.log(N) / N if sign > 0 else np.inf
    return int(np.nanargmin(bic)) + 1, bic


LAG_NOTE = ("PREVIEW LAG (biased; grain of salt): Gaussian VAR on smoothed rates, 10 ms bins, "
            "lags 10-150 ms only. The cross-kernel is CONDITIONAL on the target's own past, "
            "not the raw coupling kernel. Use the CENTROID: in simulations with a true lag of "
            "50 ms it read 52 ms (kin->neuron) and 43 ms (neuron->kin) with 'matched' kinematics "
            "(42 ms with 10 leaders in one session); 'raw' adds ~+15 ms to kin->neuron. The PEAK "
            "is unreliable (10-30 ms regardless of the true lag). Final lags come from script 03.")


def var_lag_kernel(y, x, in_window, p_max, B, lags, bin_ms, n_boot=100, seed=0):
    """LAG READOUT (preview; see LAG_NOTE) for one neuron, both directions.
    Same rows and target as var_neuron; own history = p_max raw lags; the cross term is
    the partner's lags projected on the causal basis B (n_lags, n_bases), so the
    cross-kernel is beta(tau) = B @ w -- smooth and far less collinear than 15 raw lags.
    Significance does NOT come from here (it comes from var_neuron + shuffle).
    -> {"xy": dict(tau_ms, beta, peak_ms, peak_sign, centroid_ms, peak_ci, centroid_ci,
                   beta_lo, beta_hi), "yx": ...}   (CIs: trial bootstrap, 2.5/97.5 %)"""
    from .design import basis_features
    y = np.asarray(y, float); x = np.asarray(x, float)
    Ly, Lx, rows = var_lags(y, x, in_window, p_max)
    Fy = lagged_block(y, basis_features, B, lags)
    Fx = lagged_block(x, basis_features, B, lags)
    if not (np.isfinite(Fy[:, rows]).all() and np.isfinite(Fx[:, rows]).all()):
        raise ValueError("basis lags exceed p_max rows")
    tau = np.asarray(lags) * bin_ms
    n_tr = y.shape[0]
    rng = np.random.default_rng(seed)
    boots = [np.arange(n_tr)] + [rng.integers(0, n_tr, n_tr) for _ in range(n_boot)]

    def _summ(beta):
        a = np.abs(beta)
        i = int(np.argmax(a))
        return tau[i], np.sign(beta[i]), float(np.sum(tau * a) / np.sum(a)) if a.sum() else np.nan

    out = {}
    for key, tgt, own, cross in (("xy", y, Ly, Fx), ("yx", x, Lx, Fy)):
        T = tgt[:, rows]; O = own[:, rows, :]; C = cross[:, rows, :]
        betas = []
        for bi in boots:
            t_ = T[bi].ravel()
            Z = np.column_stack([np.ones(t_.size), O[bi].reshape(-1, O.shape[2]),
                                 C[bi].reshape(-1, C.shape[2])])
            w, *_ = np.linalg.lstsq(Z, t_, rcond=None)
            betas.append(B @ w[-B.shape[1]:])
        betas = np.array(betas)
        pk, sg, ce = _summ(betas[0])
        if n_boot:
            bs = np.array([_summ(b) for b in betas[1:]])
            pci = tuple(np.nanpercentile(bs[:, 0], [2.5, 97.5]))
            cci = tuple(np.nanpercentile(bs[:, 2], [2.5, 97.5]))
        else:
            pci = cci = (np.nan, np.nan)
        out[key] = dict(tau_ms=tau, beta=betas[0], peak_ms=float(pk), peak_sign=float(sg),
                        centroid_ms=ce, peak_ci=pci, centroid_ci=cci,
                        beta_lo=np.percentile(betas[1:], 2.5, axis=0) if n_boot else None,
                        beta_hi=np.percentile(betas[1:], 97.5, axis=0) if n_boot else None)
    return out


def var_neuron(y, x, in_window, p_max, perms=None, p_fixed=None, lag_readout=None):
    """Bivariate VAR GC for one neuron, both directions, with the trial-shuffle null.
      y, x      (n_trials, n_bins) neural signal and kinematics (same grid)
      perms     (n_shuffle, n_trials) partner indices (nulls.make_perms) or None
      p_fixed   use this order instead of BIC
      lag_readout  None, or dict(B, lags, bin_ms, n_boot, seed) -> adds out[d]["lag"]
                   from var_lag_kernel (preview lag estimate, see LAG_NOTE)
    -> dict(p, bic, n_rows, rows_per_trial,
            xy=dict(gc, F, p_F, df2, null, p_shuf, z),   # x -> y (y equation)
            yx=dict(...))                                # y -> x (x equation)"""
    y = np.asarray(y, float); x = np.asarray(x, float)
    Ly, Lx, rows = var_lags(y, x, in_window, p_max)
    n_tr, n_r = y.shape[0], int(rows.sum())
    if n_r == 0:
        raise ValueError("no usable rows (window too short for p_max?)")
    yt, xt = y[:, rows].ravel(), x[:, rows].ravel()
    Yl = Ly[:, rows, :].reshape(-1, p_max)
    Xl = Lx[:, rows, :].reshape(-1, p_max)
    N = yt.size

    p, bic = (p_fixed, None) if p_fixed else select_order_bic(yt, xt, Yl, Xl, p_max)
    out = dict(p=p, bic=bic, n_rows=N, rows_per_trial=n_r)
    for key, t, own, cross_blk in (("xy", yt, Yl, Lx), ("yx", xt, Xl, Ly)):
        Q, r = _reduced(t, own[:, :p])
        rss_r = float(r @ r)
        C = cross_blk[:, rows, :p].reshape(-1, p)
        rss_f = _full_rss(Q, r, C)
        gc, F, pF, df2 = gc_ftest(rss_r, rss_f, N, 2 * p + 1, p)
        null = np.full(0 if perms is None else len(perms), np.nan)
        if perms is not None:
            for s, pi in enumerate(perms):
                Cs = cross_blk[pi][:, rows, :p].reshape(-1, p)
                null[s] = np.log(rss_r / _full_rss(Q, r, Cs))
        out[key] = dict(gc=gc, F=F, p_F=pF, df2=df2, null=null,
                        p_shuf=stats.shuffle_p(gc, null) if null.size else np.nan,
                        z=stats.null_z(gc, null) if null.size else np.nan)
    if lag_readout is not None:
        lk = var_lag_kernel(y, x, in_window, p_max, **lag_readout)
        for key in ("xy", "yx"):
            out[key]["lag"] = lk[key]
    return out


class VarNull:
    """Shuffle statistics of ONE neuron and direction for script 01, on any set of shuffles
    (adaptive extension). Same operations as var_neuron, so perms[a:b] reproduces
    var_neuron(...)[d]["null"][a:b] exactly. Setup is lazy (done in the worker)."""

    def __init__(self, y, x, in_window, p_max, p, direction):
        self.args = (np.asarray(y, float), np.asarray(x, float), np.asarray(in_window, bool),
                     int(p_max), int(p), direction)
        self._ready = False

    def _setup(self):
        y, x, inw, p_max, p, d = self.args
        Ly, Lx, rows = var_lags(y, x, inw, p_max)
        yt, xt = y[:, rows].ravel(), x[:, rows].ravel()
        own = (Ly if d == "xy" else Lx)[:, rows, :].reshape(-1, p_max)
        self.cross_blk, self.rows, self.p = (Lx if d == "xy" else Ly), rows, p
        self.Q, self.r = _reduced(yt if d == "xy" else xt, own[:, :p])
        self.rss_r = float(self.r @ self.r)
        self._ready = True

    def stats(self, perms):
        if not self._ready:
            self._setup()
        p, rows = self.p, self.rows
        return np.array([np.log(self.rss_r / _full_rss(self.Q, self.r,
                                                        self.cross_blk[pi][:, rows, :p].reshape(-1, p)))
                         for pi in perms])


# ---------------------------------------------------------------------------
# Point-process GLM with discrete lag windows (script 02)
# ---------------------------------------------------------------------------
GLM_LAG_NOTE = ("LAG WINDOWS (coarse; grain of salt): cross coefficients per lag window, "
                "CONDITIONAL on the target's own history (same caveat as 01). Read the window "
                "with the largest |beta/SE| and the t-weighted centroid over windows. In "
                "simulations with a true lag of 50 ms: kin->neuron strongest window 40-80 ms "
                "(centroid 65 ms); neuron->kin strongest window 20-40 ms (centroid 37 ms), i.e. "
                "neuron->kin lags read SHORT -- and SHORTER when several neurons lead the same "
                "kinematics (10 leaders in one session: strongest window 10-20 ms, centroid 27 ms), "
                "because the kinematics' own history then holds the older contributions of all "
                "leaders. Treat neuron->kin windows as an UNDERESTIMATE. Final lags come from 03.")


def _pen(k, ridge):
    """Ridge penalty matrix: ridge on every column except the intercept (column 0)."""
    P = np.eye(k) * ridge
    P[0, 0] = 0.0
    return P


def poisson_irls(X, y, ridge=0.0, maxit=100, tol=1e-9):
    """Poisson GLM (log link) by damped Newton / IRLS, optional ridge (not on column 0,
    which must be the constant).
    -> dict(b, mu, deviance, fisher (X'WX + P), cov (fisher^-1), converged, n_iter)."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    n, k = X.shape
    P = _pen(k, ridge)
    b = np.zeros(k); b[0] = np.log(max(y.mean(), 1e-8))

    def _obj(b):                                   # negative penalised log-likelihood
        eta = X @ b
        return float(np.sum(np.exp(eta) - y * eta) + 0.5 * b @ P @ b)

    f_old, conv = _obj(b), False
    for it in range(maxit):
        mu = np.exp(X @ b)
        g = X.T @ (y - mu) - P @ b
        H = (X * mu[:, None]).T @ X + P
        step = np.linalg.solve(H + 1e-10 * np.eye(k), g)
        t, ok = 1.0, False
        with np.errstate(over="ignore", invalid="ignore"):
            while t >= 1e-10:                     # step halving: finite AND objective falling
                bn = b + t * step
                fn = _obj(bn)
                if np.isfinite(fn) and fn <= f_old + 1e-12:
                    ok = True
                    break
                t *= 0.5
        if not ok:                                # no improving step left: at the optimum
            conv = np.max(np.abs(g)) < 1e-6 * max(1.0, np.abs(y).sum())
            break
        b, f_old = bn, fn
        if np.max(np.abs(t * step)) < tol:
            conv = True
            break
    mu = np.exp(X @ b)
    H = (X * mu[:, None]).T @ X + P
    return dict(b=b, mu=mu, deviance=poisson_deviance(y, mu), fisher=H,
                cov=np.linalg.inv(H), converged=conv, n_iter=it + 1)


def poisson_deviance(y, mu):
    """2 * sum[y log(y/mu) - (y - mu)], with y log y = 0 at y = 0."""
    y = np.asarray(y, float)
    t = np.where(y > 0, y * np.log(np.where(y > 0, y, 1.0) / mu), 0.0)
    return float(2 * np.sum(t - (y - mu)))


class PoissonScore:
    """Score (Rao) statistic for ADDING columns Xc to a fitted reduced Poisson model.
    The reduced fit is done once; each call only needs the new columns, so the trial
    shuffle (which only changes Xc) costs a few matrix products per shuffle.
      S = U' I^-1 U,  U = Xc'(y - mu_r),
      I = Xc'W Xc - (Xr'W Xc)' (Xr'W Xr + P_r)^-1 (Xr'W Xc) + ridge * Id,  W = diag(mu_r)."""

    def __init__(self, Xr, y, ridge=0.0):
        self.Xr, self.y, self.ridge = np.asarray(Xr, float), np.asarray(y, float), ridge
        self.fit = poisson_irls(self.Xr, self.y, ridge)
        self.mu = self.fit["mu"]
        self.resid = self.y - self.mu
        self.WXr = self.Xr * self.mu[:, None]
        self.Ainv = np.linalg.inv(self.WXr.T @ self.Xr + _pen(self.Xr.shape[1], ridge))

    def stat(self, Xc):
        Xc = np.asarray(Xc, float)
        U = Xc.T @ self.resid
        G = self.WXr.T @ Xc
        I = (Xc * self.mu[:, None]).T @ Xc - G.T @ self.Ainv @ G + self.ridge * np.eye(Xc.shape[1])
        return float(U @ np.linalg.solve(I, U))


def _window_lag_summary(beta, se, centers_ms):
    """Coarse lag readout: window with the largest |beta/SE| and the |t|-weighted centroid."""
    t = np.abs(beta / se)
    i = int(np.nanargmax(t))
    cen = float(np.sum(centers_ms * t) / np.sum(t)) if np.sum(t) > 0 else np.nan
    return dict(beta=beta, se=se, t=beta / se, best_ms=float(centers_ms[i]),
                best_sign=float(np.sign(beta[i])), centroid_ms=cen)


def glm_neuron(y, x, in_window, hist_sets, cross_sets, cross_centers_ms, perms=None,
               ridge=0.0):
    """Script 02, one neuron, both directions, on raw 2.5 ms spike counts y and raw
    kinematics x (n_trials, n_bins). Every predictor is a WINDOW MEAN of past lags
    (design.window_means); own history and cross terms use the lag sets given.
      x -> y (neuron eq.):     Poisson GLM   y ~ 1 + own(y) [+ cross(x)]
      y -> x (kinematics eq.): Gaussian LM   x ~ 1 + own(x) [+ cross(y)]
    Rows: in_window bins where every window of both signals exists in every trial
    (same rows for full, reduced and every shuffle; rule 3).
    Shuffle statistic: Poisson score statistic (x->y), ln RSS ratio (y->x). Observed
    LR deviance / F-test p are SECONDARY.
    ridge: L2 penalty for the Poisson equation (identical in reduced, full and every
    shuffle; not on the intercept). The Gaussian equation is unpenalised.
    -> dict(n_rows, rows_per_trial, n_spikes,
            xy=dict(lr, p_chi2, score, null, p_shuf, z, dev_expl_red, dev_expl_full,
                    bits_per_spike, dispersion, lag=dict(...)),
            yx=dict(gc, F, p_F, null, p_shuf, z, r2_red, r2_full, partial_r2, lag=dict(...)))"""
    from scipy.stats import chi2
    from .design import window_means
    y = np.asarray(y, float); x = np.asarray(x, float)
    Hy = lagged_block(y, window_means, hist_sets)
    Hx = lagged_block(x, window_means, hist_sets)
    Cy = lagged_block(y, window_means, cross_sets)
    Cx = lagged_block(x, window_means, cross_sets)
    fin = (np.isfinite(y) & np.isfinite(x) & np.isfinite(Hy).all(2) & np.isfinite(Hx).all(2)
           & np.isfinite(Cy).all(2) & np.isfinite(Cx).all(2))
    rows = np.asarray(in_window, bool) & fin.all(axis=0)
    n_r = int(rows.sum())
    if n_r == 0:
        raise ValueError("no usable rows")
    yt, xt = y[:, rows].ravel(), x[:, rows].ravel()
    q = len(cross_sets)
    flat = lambda A: A[:, rows, :].reshape(-1, A.shape[2])
    out = dict(n_rows=yt.size, rows_per_trial=n_r, n_spikes=float(yt.sum()))

    # ---- x -> y: Poisson ------------------------------------------------------
    # Columns are standardised (numerics: counts vs mm differ by orders of magnitude).
    # LR, score statistic and t = beta/SE are invariant to this; a trial shuffle only
    # reassigns rows between trials, so shuffled columns have the SAME mean/sd and get
    # the identical transform. Ridge (if any) acts on the standardised coefficients.
    def _z(A, m=None, s=None):
        m = A.mean(0) if m is None else m
        s = A.std(0) if s is None else s
        s = np.where(s > 0, s, 1.0)
        return (A - m) / s, m, s
    one = np.ones((yt.size, 1))
    Hz, _, _ = _z(flat(Hy))
    Xr = np.hstack([one, Hz])
    sc = PoissonScore(Xr, yt, ridge)
    Xc, mc, sdc = _z(flat(Cx))
    full = poisson_irls(np.hstack([Xr, Xc]), yt, ridge)
    dev0 = poisson_deviance(yt, np.full(yt.size, yt.mean()))
    lr = sc.fit["deviance"] - full["deviance"]
    s_obs = sc.stat(Xc)
    null = (np.array([sc.stat(_z(flat(Cx[pi]), mc, sdc)[0]) for pi in perms])
            if perms is not None else np.zeros(0))
    k_r = Xr.shape[1]
    b_c = full["b"][k_r:]; se_c = np.sqrt(np.diag(full["cov"])[k_r:])
    pear = float(np.sum(sc.resid ** 2 / sc.mu) / (yt.size - k_r))
    out["xy"] = dict(lr=lr, p_chi2=float(chi2.sf(lr, q)), score=s_obs, null=null,
                     p_shuf=stats.shuffle_p(s_obs, null) if null.size else np.nan,
                     z=stats.null_z(s_obs, null) if null.size else np.nan,
                     dev_expl_red=1 - sc.fit["deviance"] / dev0,
                     dev_expl_full=1 - full["deviance"] / dev0,
                     bits_per_spike=(lr / 2) / (max(yt.sum(), 1) * np.log(2)),
                     dispersion=pear, converged=bool(full["converged"] and sc.fit["converged"]),
                     lag=_window_lag_summary(b_c, se_c, np.asarray(cross_centers_ms, float)))

    # ---- y -> x: Gaussian (same machinery as 01) --------------------------------
    O = flat(Hx); C = flat(Cy)
    Q, r = _reduced(xt, O)
    rss_r = float(r @ r); rss_f = _full_rss(Q, r, C)
    gc, F, pF, df2 = gc_ftest(rss_r, rss_f, xt.size, 1 + O.shape[1] + q, q)
    nullx = (np.array([np.log(rss_r / _full_rss(Q, r, flat(Cy[pi]))) for pi in perms])
             if perms is not None else np.zeros(0))
    Z = np.column_stack([np.ones(xt.size), O, C])
    w, *_ = np.linalg.lstsq(Z, xt, rcond=None)
    sig2 = rss_f / df2
    cov = sig2 * np.linalg.pinv(Z.T @ Z)
    tss = float(np.sum((xt - xt.mean()) ** 2))
    out["yx"] = dict(gc=gc, F=F, p_F=pF, null=nullx,
                     p_shuf=stats.shuffle_p(gc, nullx) if nullx.size else np.nan,
                     z=stats.null_z(gc, nullx) if nullx.size else np.nan,
                     r2_red=1 - rss_r / tss, r2_full=1 - rss_f / tss, partial_r2=1 - rss_f / rss_r,
                     lag=_window_lag_summary(w[-q:], np.sqrt(np.diag(cov)[-q:]),
                                             np.asarray(cross_centers_ms, float)))
    return out


class GlmNull:
    """Shuffle statistics of ONE neuron and direction for script 02 (x->y: Poisson score
    statistic; y->x: ln RSS ratio) on any set of shuffles. Same operations as glm_neuron,
    so perms[a:b] reproduces glm_neuron(...)[d]["null"][a:b] exactly. Lazy setup."""

    def __init__(self, y, x, in_window, hist_sets, cross_sets, ridge, direction):
        self.args = (np.asarray(y, float), np.asarray(x, float), np.asarray(in_window, bool),
                     hist_sets, cross_sets, ridge, direction)
        self._ready = False

    def _setup(self):
        from .design import window_means
        y, x, inw, hist_sets, cross_sets, ridge, d = self.args
        Hy = lagged_block(y, window_means, hist_sets)
        Hx = lagged_block(x, window_means, hist_sets)
        Cy = lagged_block(y, window_means, cross_sets)
        Cx = lagged_block(x, window_means, cross_sets)
        fin = (np.isfinite(y) & np.isfinite(x) & np.isfinite(Hy).all(2) & np.isfinite(Hx).all(2)
               & np.isfinite(Cy).all(2) & np.isfinite(Cx).all(2))
        self.rows = rows = inw & fin.all(axis=0)
        yt, xt = y[:, rows].ravel(), x[:, rows].ravel()
        flat = self._flat
        if d == "xy":
            def _z(A, m=None, s=None):
                m = A.mean(0) if m is None else m
                s = A.std(0) if s is None else s
                s = np.where(s > 0, s, 1.0)
                return (A - m) / s, m, s
            self._zf = _z
            Hz, _, _ = _z(flat(Hy))
            self.sc = PoissonScore(np.hstack([np.ones((yt.size, 1)), Hz]), yt, ridge)
            _, self.mc, self.sdc = _z(flat(Cx))
            self.cross = Cx
        else:
            self.Q, self.r = _reduced(xt, flat(Hx))
            self.rss_r = float(self.r @ self.r)
            self.cross = Cy
        self._ready = True

    def _flat(self, A):
        return A[:, self.rows, :].reshape(-1, A.shape[2])

    def stats(self, perms):
        if not self._ready:
            self._setup()
        if self.args[-1] == "xy":
            return np.array([self.sc.stat(self._zf(self._flat(self.cross[pi]), self.mc, self.sdc)[0])
                             for pi in perms])
        return np.array([np.log(self.rss_r / _full_rss(self.Q, self.r, self._flat(self.cross[pi])))
                         for pi in perms])
