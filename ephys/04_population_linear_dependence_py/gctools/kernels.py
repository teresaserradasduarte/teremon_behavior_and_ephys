"""kernels.py - script 03: point-process GLM with smooth causal kernels (main analysis).

Per neuron, both directions, on raw 2.5 ms spike counts y and raw med-lat x
(n_trials, n_bins), with predictors built from causal bases (design.py):

  kin -> neuron (x->y), Poisson GLM, log link:
      log lam_t = b0 + own_y(H) [+ event] [+ nuisance] + cross(x)
      own history H: hist_mode "cv" = chosen PER NEURON by trial-wise CV on the REDUCED
      model among h_list (never shorter than the cross range -- the 01 lesson), "fixed" =
      h_fixed (150 ms), "max" = the longest candidate.
  neuron -> kin (y->x), Gaussian LM:
      x_t = c + own_x(H_kin) [+ event] [+ nuisance] + cross(y)
      H_kin: kin_hist_mode "cv" (CV on the reduced Gaussian model, held-out SSE, among
      h_kin_list) | "fixed" (h_kin, 150 ms) | "max". Kinematic history before -375 ms only
      exists with io kin_source="long" (REACH_LONG); with "short" a long H_kin cuts rows.
  Own-history bases: 2 single lags + design.n_log_for(H) log bumps (~2 more per doubling of
  H; scale_bases=False keeps 6 for every H, as before).

ROWS (rule 3): bins in the window where EVERY candidate term exists in every trial --
identical for all H candidates, reduced/full, every CV fold and every shuffle.

Significance: trial-shuffle (Poisson score statistic / ln RSS ratio), as in 02; H is
chosen on the reduced model (no cross terms), so shuffles never re-select it.
Effect size: trial-wise K-fold CROSS-VALIDATED gain of full over reduced (bits/spike;
fraction of held-out med-lat variance) -- ~<= 0 by chance, comparable across directions
in the normalised sense of rule 5.
Lag: kernel beta(tau) = B_cross @ w (raw units), peak |beta| and centroid; CIs by trial
bootstrap (kernel_lag_bootstrap, run for FDR-significant neurons only).
Adequacy: does the reduced model's residual still depend on the target's own past
beyond H? (x->y: LR test adding own-history bases over (H, H+150] ms; y->x: F-test of
x bases over (H_kin, H_kin+150] ms on the rows where that history exists). Flagged, not removed.
"""
import numpy as np
from scipy.stats import chi2, f as _f
from . import design as ds
from . import stats
from .models import (poisson_irls, poisson_deviance, PoissonScore, _reduced, _full_rss,
                     gc_ftest)

KERNEL_LAG_NOTE = ("KERNEL LAGS (script 03): beta(tau) = B @ w is CONDITIONAL on the target's own "
                   "history, as in 01/02. Use the CENTROID (sims, true 50 ms: kin->neuron 57 ms; "
                   "neuron->kin 51 ms with 1 leader, 38 ms with 5 -- reads short with several "
                   "leaders); the PEAK is unreliable (often lands on the first narrow bump). "
                   "Full bias table: script 04.")


# ---------------------------------------------------------------------------
# design
# ---------------------------------------------------------------------------
def _candidates(mode, cand, fixed):
    if mode == "cv":
        return tuple(cand)
    if mode == "fixed":
        return (fixed,)
    if mode == "max":
        return (max(cand),)
    raise ValueError(f"history mode must be 'cv' | 'fixed' | 'max', got {mode!r}")


def make_spec(bin_ms=2.5, h_list=(150, 200, 300), hist_mode="cv", h_kin=150,
              n_cross=6, cross_range=(0, 150), spacing="log", event=False,
              event_range=(0, 125), n_event=4, nuisance=False, ridge=0.1, n_folds=5,
              adequacy_ms=150, seed=0, h_fixed=150, h_kin_list=(150,), kin_hist_mode="fixed",
              scale_bases=True, fit_from_ms=None):
    """All modelling choices of script 03 in one dict (also the cache key).
      h_list / hist_mode           neuron own history: candidates, "cv" | "fixed" (h_fixed) | "max"
      h_kin_list / kin_hist_mode   kinematics own history: candidates, "cv" | "fixed" (h_kin) | "max"
      scale_bases                  log bumps per history = design.n_log_for(H) (else 6)
      fit_from_ms                  fitted rows start no earlier than this (ms re max reach;
                                   None = wherever every term exists)"""
    Bc, Lc, Cc = ds.raised_cosine_basis(n_cross, cross_range, bin_ms, spacing)
    spec = dict(bin_ms=bin_ms, h_list=_candidates(hist_mode, h_list, h_fixed),
                hist_mode=hist_mode, h_kin=h_kin, h_fixed=h_fixed,
                h_kin_list=_candidates(kin_hist_mode, h_kin_list, h_kin), kin_hist_mode=kin_hist_mode,
                scale_bases=bool(scale_bases), fit_from_ms=fit_from_ms,
                n_cross=n_cross, cross_range=tuple(cross_range),
                spacing=spacing, B_cross=Bc, L_cross=Lc, cross_centers=Cc, event=event,
                nuisance=nuisance, ridge=ridge, n_folds=n_folds, adequacy_ms=adequacy_ms, seed=seed)
    if event:
        spec["B_ev"], spec["L_ev"], _ = ds.raised_cosine_basis(n_event, event_range, bin_ms, "log")
    return spec


def pad_needed_ms(spec):
    """History the binning must provide before the window (own history + adequacy)."""
    return max(max(spec["h_list"]), max(spec["h_kin_list"])) + spec["adequacy_ms"]


def own_bases(H, spec):
    """Own-history basis of length H ms (2 single lags + log bumps)."""
    n_log = ds.n_log_for(H) if spec.get("scale_bases", False) else 6
    return ds.history_bases(H, spec["bin_ms"], n_log=n_log)


def _feat(Z, B, L):
    return ds.lagged_block(Z, ds.basis_features, B, L)


def _extra_bases(H, extra, bin_ms):
    """Own-history bases over (H, H+extra] only (adequacy test): 3 log bumps."""
    B, L, _ = ds.raised_cosine_basis(3, (H, H + extra), bin_ms, "log")
    return B, L


def build_design(y, x, in_window, cond, t, spec):
    """Per-trial features for every candidate term, then the common rows.
    -> dict of flattened (n_rows, k) blocks + targets + row trial ids + per-trial blocks
       (the per-trial cross blocks are kept un-flattened for shuffles)."""
    y = np.asarray(y, float); x = np.asarray(x, float)
    b = spec["bin_ms"]
    blk = {}
    for H in spec["h_list"]:
        blk[f"own_y{H}"] = _feat(y, *own_bases(H, spec))
        blk[f"adq_y{H}"] = _feat(y, *_extra_bases(H, spec["adequacy_ms"], b))
    for H in spec["h_kin_list"]:
        blk[f"own_x{H}"] = _feat(x, *own_bases(H, spec))
    blk["cross_x"] = _feat(x, spec["B_cross"], spec["L_cross"])
    blk["cross_y"] = _feat(y, spec["B_cross"], spec["L_cross"])
    if spec["event"]:
        ev = ds.event_regressors(t, 0.0, spec["B_ev"], spec["L_ev"])
        blk["event"] = np.broadcast_to(ev, (y.shape[0],) + ev.shape).copy()
    if spec["nuisance"]:
        blk["nui_y"] = np.log(ds.nuisance_mean(y, cond) + 1e-3)[:, :, None]
        blk["nui_x"] = ds.nuisance_mean(x, cond)[:, :, None]
    # rows: every term except the kinematic adequacy blocks must exist (those are tested
    # on their own subset of rows, where H_kin + adequacy_ms of kinematic history exists)
    fin = np.isfinite(y) & np.isfinite(x)
    for A in blk.values():
        fin &= np.isfinite(A).all(axis=2)
    rows = np.asarray(in_window, bool) & fin.all(axis=0)
    if spec.get("fit_from_ms") is not None:
        rows &= np.asarray(t) >= spec["fit_from_ms"] / 1000.0 - 1e-9
    n_tr = y.shape[0]
    flat = lambda A: A[:, rows, :].reshape(-1, A.shape[2])
    D = {k: flat(v) for k, v in blk.items()}
    D.update(yt=y[:, rows].ravel(), xt=x[:, rows].ravel(), rows=rows,
             trial=np.repeat(np.arange(n_tr), int(rows.sum())),
             cross_x_trials=blk["cross_x"][:, rows, :], cross_y_trials=blk["cross_y"][:, rows, :])
    for H in spec["h_kin_list"]:
        A = flat(_feat(x, *_extra_bases(H, spec["adequacy_ms"], b)))
        D[f"adq_x{H}"], D[f"adq_x{H}_ok"] = A, np.isfinite(A).all(axis=1)
    return D


def _z(A, m=None, s=None):
    if A.shape[1] == 0:
        return A, np.zeros(0), np.ones(0)
    m = A.mean(0) if m is None else m
    s = A.std(0) if s is None else s
    s = np.where(s > 0, s, 1.0)
    return (A - m) / s, m, s


def _nuis(D, which):
    keys = (["event"] if "event" in D else []) + ([f"nui_{which}"] if f"nui_{which}" in D else [])
    return np.hstack([D[k] for k in keys]) if keys else np.zeros((D["yt"].size, 0))


def _poisson_ll(y, mu):
    return float(np.sum(y * np.log(np.maximum(mu, 1e-300)) - mu))


def folds_by_trial(cond, n_folds, seed):
    """Trial-wise K-fold assignment, stratified by condition (spout side)."""
    cond = np.asarray(cond)
    rng = np.random.default_rng(seed)
    f = np.empty(cond.size, int)
    for c in np.unique(cond):
        idx = rng.permutation(np.flatnonzero(cond == c))
        f[idx] = np.arange(idx.size) % n_folds
    return f


# ---------------------------------------------------------------------------
# per neuron
# ---------------------------------------------------------------------------
def kernel_neuron(y, x, in_window, cond, t, spec, perms=None):
    """Script 03 for one neuron, both directions (see module docstring).
    -> dict(n_rows, rows_per_trial, n_spikes, H, cv_ll_H, H_kin, cv_sse_Hkin,
            xy=dict(lr, p_chi2, score, null, p_shuf, z, cv_gain_bits, beta, tau_ms, peak_ms,
                    peak_sign, centroid_ms, adequacy_p, dispersion, converged),
            yx=dict(gc, F, p_F, null, p_shuf, z, cv_gain_r2, beta, tau_ms, peak_ms, peak_sign,
                    centroid_ms, adequacy_p))"""
    D = build_design(y, x, in_window, cond, t, spec)
    yt, xt, tr = D["yt"], D["xt"], D["trial"]
    if yt.size == 0:
        raise ValueError("no usable rows")
    fold = folds_by_trial(cond, spec["n_folds"], spec["seed"])[tr]
    ridge, K = spec["ridge"], spec["n_folds"]
    tau = spec["L_cross"] * spec["bin_ms"]
    out = dict(n_rows=int(yt.size), rows_per_trial=int(D["rows"].sum()), n_spikes=float(yt.sum()))
    one = np.ones((yt.size, 1))

    # ---- kin -> neuron: choose H by CV on the REDUCED model ---------------------
    nuis_y = _z(_nuis(D, "y"))[0]
    cv_ll = {}
    for H in spec["h_list"]:
        Xr = np.hstack([one, _z(D[f"own_y{H}"])[0], nuis_y])
        if len(spec["h_list"]) == 1:
            cv_ll[H] = np.nan
            break
        ll = 0.0
        for k in range(K):
            te = fold == k
            fit = poisson_irls(Xr[~te], yt[~te], ridge)
            ll += _poisson_ll(yt[te], np.exp(Xr[te] @ fit["b"]))
        cv_ll[H] = ll
    H = (max(cv_ll, key=lambda h: cv_ll[h]) if len(spec["h_list"]) > 1 else spec["h_list"][0])
    out["H"], out["cv_ll_H"] = int(H), cv_ll
    Xr = np.hstack([one, _z(D[f"own_y{H}"])[0], nuis_y])
    Xc, mc, sc = _z(D["cross_x"])
    sco = PoissonScore(Xr, yt, ridge)
    full = poisson_irls(np.hstack([Xr, Xc]), yt, ridge)
    lr = sco.fit["deviance"] - full["deviance"]
    s_obs = sco.stat(Xc)
    null = (np.array([sco.stat(_z(D["cross_x_trials"][pi].reshape(-1, Xc.shape[1]), mc, sc)[0])
                      for pi in perms]) if perms is not None else np.zeros(0))
    # CV gain (same folds, same rows): held-out bits/spike of full over reduced
    dll, nsp = 0.0, 0.0
    for k in range(K):
        te = fold == k
        fr = poisson_irls(Xr[~te], yt[~te], ridge)
        ff = poisson_irls(np.hstack([Xr, Xc])[~te], yt[~te], ridge)
        dll += (_poisson_ll(yt[te], np.exp(np.hstack([Xr, Xc])[te] @ ff["b"]))
                - _poisson_ll(yt[te], np.exp(Xr[te] @ fr["b"])))
        nsp += yt[te].sum()
    k_r = Xr.shape[1]
    beta = spec["B_cross"] @ (full["b"][k_r:] / sc)               # raw units (per mm)
    # adequacy: does own history beyond H still predict the neuron?
    Xa = np.hstack([Xr, _z(D[f"adq_y{H}"])[0]])
    adq = poisson_irls(Xa, yt, ridge)
    lr_a = sco.fit["deviance"] - adq["deviance"]
    pear = float(np.sum(sco.resid ** 2 / sco.mu) / (yt.size - k_r))
    out["xy"] = dict(lr=lr, p_chi2=float(chi2.sf(lr, Xc.shape[1])), score=s_obs, null=null,
                     p_shuf=stats.shuffle_p(s_obs, null) if null.size else np.nan,
                     z=stats.null_z(s_obs, null) if null.size else np.nan,
                     cv_gain_bits=dll / (max(nsp, 1.0) * np.log(2)),
                     **_lag_summary(beta, tau),
                     adequacy_p=float(chi2.sf(lr_a, Xa.shape[1] - k_r)),
                     dispersion=pear, converged=bool(full["converged"] and sco.fit["converged"]))

    # ---- neuron -> kin: Gaussian; choose H_kin by CV on the REDUCED model (held-out SSE)
    cv_sse = {}
    if len(spec["h_kin_list"]) > 1:
        for Hk in spec["h_kin_list"]:
            Zr = np.column_stack([np.ones(xt.size), D[f"own_x{Hk}"], _nuis(D, "x")])
            sse = 0.0
            for k in range(K):
                te = fold == k
                w, *_ = np.linalg.lstsq(Zr[~te], xt[~te], rcond=None)
                sse += float(np.sum((xt[te] - Zr[te] @ w) ** 2))
            cv_sse[Hk] = sse
        Hk = min(cv_sse, key=lambda h: cv_sse[h])
    else:
        Hk = spec["h_kin_list"][0]
        cv_sse[Hk] = np.nan
    out["H_kin"], out["cv_sse_Hkin"] = int(Hk), cv_sse
    O = np.hstack([D[f"own_x{Hk}"], _nuis(D, "x")])
    Q, r = _reduced(xt, O)
    rss_r = float(r @ r)
    C = D["cross_y"]
    rss_f = _full_rss(Q, r, C)
    gc, F, pF, df2 = gc_ftest(rss_r, rss_f, xt.size, 1 + O.shape[1] + C.shape[1], C.shape[1])
    nullx = (np.array([np.log(rss_r / _full_rss(Q, r, D["cross_y_trials"][pi].reshape(-1, C.shape[1])))
                       for pi in perms]) if perms is not None else np.zeros(0))
    Zr = np.column_stack([np.ones(xt.size), O]); Zf = np.column_stack([Zr, C])
    sse_r = sse_f = 0.0
    for k in range(K):
        te = fold == k
        wr, *_ = np.linalg.lstsq(Zr[~te], xt[~te], rcond=None)
        wf, *_ = np.linalg.lstsq(Zf[~te], xt[~te], rcond=None)
        sse_r += float(np.sum((xt[te] - Zr[te] @ wr) ** 2))
        sse_f += float(np.sum((xt[te] - Zf[te] @ wf) ** 2))
    wf, *_ = np.linalg.lstsq(Zf, xt, rcond=None)
    beta_x = spec["B_cross"] @ wf[-C.shape[1]:]
    # adequacy: residual of the reduced model vs kinematics (H_kin, H_kin+150] back (where available)
    ok, Ax = D[f"adq_x{Hk}_ok"], D[f"adq_x{Hk}"]
    if ok.sum() > 50:
        _, r0 = _reduced(xt[ok], O[ok])                               # own history only
        _, r1 = _reduced(xt[ok], np.hstack([O[ok], Ax[ok]]))          # + older history
        rss0, rss1 = float(r0 @ r0), float(r1 @ r1)
        q_a = Ax.shape[1]
        Fa = ((rss0 - rss1) / q_a) / (rss1 / (ok.sum() - 1 - O.shape[1] - q_a))
        p_adq_x = float(_f.sf(Fa, q_a, ok.sum() - 1 - O.shape[1] - q_a))
    else:
        p_adq_x = np.nan
    out["yx"] = dict(gc=gc, F=F, p_F=pF, null=nullx,
                     p_shuf=stats.shuffle_p(gc, nullx) if nullx.size else np.nan,
                     z=stats.null_z(gc, nullx) if nullx.size else np.nan,
                     cv_gain_r2=1.0 - sse_f / sse_r, **_lag_summary(beta_x, tau),
                     adequacy_p=p_adq_x)
    return out


def _lag_summary(beta, tau):
    a = np.abs(beta)
    i = int(np.argmax(a))
    return dict(beta=beta, tau_ms=tau, peak_ms=float(tau[i]), peak_sign=float(np.sign(beta[i])),
                centroid_ms=float(np.sum(tau * a) / a.sum()) if a.sum() > 0 else np.nan)


class KernelNull:
    """Shuffle statistics of one neuron and direction for script 03 on any shuffles
    (adaptive extension); same operations as kernel_neuron, so perms[a:b] reproduces its
    stored null exactly. H = own history of the TARGET of this direction: the neuron's H
    for "xy", the session's H_kin for "yx". Lazy setup."""

    def __init__(self, y, x, in_window, cond, t, spec, H, direction):
        self.args = (y, x, in_window, cond, t, spec, int(H), direction)
        self._ready = False

    def _setup(self):
        y, x, inw, cond, t, spec, H, d = self.args
        D = build_design(y, x, inw, cond, t, spec)
        self.d = d
        if d == "xy":
            one = np.ones((D["yt"].size, 1))
            Xr = np.hstack([one, _z(D[f"own_y{H}"])[0], _z(_nuis(D, "y"))[0]])
            self.sco = PoissonScore(Xr, D["yt"], spec["ridge"])
            _, self.mc, self.sc = _z(D["cross_x"])
            self.cross = D["cross_x_trials"]
        else:
            O = np.hstack([D[f"own_x{H}"], _nuis(D, "x")])
            self.Q, self.r = _reduced(D["xt"], O)
            self.rss_r = float(self.r @ self.r)
            self.cross = D["cross_y_trials"]
        self._ready = True

    def stats(self, perms):
        if not self._ready:
            self._setup()
        k = self.cross.shape[2]
        if self.d == "xy":
            return np.array([self.sco.stat(_z(self.cross[pi].reshape(-1, k), self.mc, self.sc)[0])
                             for pi in perms])
        return np.array([np.log(self.rss_r / _full_rss(self.Q, self.r, self.cross[pi].reshape(-1, k)))
                         for pi in perms])


def kernel_lag_bootstrap(y, x, in_window, cond, t, spec, H, direction, n_boot=200, seed=0):
    """Trial bootstrap of the lag kernel (refit the FULL model on resampled trials).
    H = own history of the TARGET (the neuron's H for "xy", H_kin for "yx").
    -> dict(beta_lo, beta_hi (2.5/97.5 % per lag), peak_ci, centroid_ci, n_boot)."""
    D = build_design(y, x, in_window, cond, t, spec)
    tr, n_tr = D["trial"], np.asarray(y).shape[0]
    rows_of = [np.flatnonzero(tr == i) for i in range(n_tr)]
    tau = spec["L_cross"] * spec["bin_ms"]
    rng = np.random.default_rng(seed)
    if direction == "xy":
        one = np.ones((D["yt"].size, 1))
        Xr = np.hstack([one, _z(D[f"own_y{H}"])[0], _z(_nuis(D, "y"))[0]])
        Xc, _, sc = _z(D["cross_x"])
        X, tgt, k_r = np.hstack([Xr, Xc]), D["yt"], Xr.shape[1]
    else:
        Zr = np.column_stack([np.ones(D["xt"].size), D[f"own_x{H}"], _nuis(D, "x")])
        X, tgt, k_r = np.column_stack([Zr, D["cross_y"]]), D["xt"], Zr.shape[1]
    betas, peaks, cents = [], [], []
    for _ in range(n_boot):
        idx = np.concatenate([rows_of[i] for i in rng.integers(0, n_tr, n_tr)])
        if direction == "xy":
            w = poisson_irls(X[idx], tgt[idx], spec["ridge"])["b"][k_r:] / sc
        else:
            w = np.linalg.lstsq(X[idx], tgt[idx], rcond=None)[0][k_r:]
        b = spec["B_cross"] @ w
        s = _lag_summary(b, tau)
        betas.append(b); peaks.append(s["peak_ms"]); cents.append(s["centroid_ms"])
    betas = np.array(betas)
    return dict(beta_lo=np.percentile(betas, 2.5, 0), beta_hi=np.percentile(betas, 97.5, 0),
                peak_ci=tuple(np.percentile(peaks, [2.5, 97.5])),
                centroid_ci=tuple(np.nanpercentile(cents, [2.5, 97.5])), n_boot=n_boot)
