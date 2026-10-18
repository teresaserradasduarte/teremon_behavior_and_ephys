"""design.py - per-trial lagged design matrices.

LAG CONVENTION. Lag j (in bins) at row t reads bin t-j. Every predictor has j >= 1:
a lag-0 term would be instantaneous coupling, not predictive precedence. On a
bin_ms grid, lag j covers delays in ((j-1)*bin, j*bin]; a window (lo, hi] ms maps to
lags lo/bin+1 .. hi/bin.

TRIAL BOUNDARIES. Lagged features are built from ONE trial's own series at a time
(lagged_block loops over trials); a lag that would reach before the trial's first
bin is NaN, and any row with a NaN feature is dropped by build_design. Stacking
across trials happens only after that, so no lag ever crosses a trial boundary.

FULL vs REDUCED. build_design computes the row mask ONCE, for the full model;
reduced() only removes column groups, so both models use exactly the same rows
and the same y (rule 3).
"""
import numpy as np


def lag_matrix(z, lags, fill=np.nan):
    """ONE trial. z (n_bins,), lags (n_lags,) ints >= 1 -> (n_bins, n_lags) with
    out[t, i] = z[t - lags[i]]; `fill` where t - lag < 0."""
    z = np.asarray(z, float)
    lags = np.asarray(lags, int)
    if lags.ndim != 1 or np.any(lags < 1):
        raise ValueError("all lags must be >= 1 bin (no contemporaneous terms)")
    n = z.shape[0]
    out = np.full((n, lags.size), fill, dtype=float)
    for i, L in enumerate(lags):
        if L < n:
            out[L:, i] = z[:n - L]
    return out


def ms_windows_to_lags(windows_ms, bin_ms):
    """[(lo, hi], ...] ms -> list of lag-index arrays (bins). Raises if a window holds
    no whole lag at this bin size."""
    out = []
    for lo, hi in windows_ms:
        a = int(np.floor(lo / bin_ms + 1e-9)) + 1
        b = int(np.floor(hi / bin_ms + 1e-9))
        if b < a:
            raise ValueError(f"window ({lo}, {hi}] ms holds no lag at {bin_ms} ms bins")
        out.append(np.arange(a, b + 1))
    return out


def window_means(z, lag_sets):
    """ONE trial. Mean of z over each set of lags (discrete lag windows, script 02)
    -> (n_bins, n_windows). NaN if any lag in the window is unavailable."""
    return np.column_stack([lag_matrix(z, L).mean(axis=1) for L in lag_sets])


def raised_cosine_basis(n_bases, lag_range_ms, bin_ms, spacing="log", log_offset_ms=5.0):
    """Raised-cosine bumps over lags (lo, hi] ms (Pillow et al. 2008 style).
    Lags start at >= 1 bin. In a warped axis u (u = tau for "uniform",
    u = log(tau + offset) for "log") the centres are evenly spaced from u(first lag)
    to u(last lag) with spacing dc, and each bump is
        0.5 * (1 + cos(clip((u - c) * pi / (2 dc), -pi, pi)))
    -> (B (n_lags, n_bases), lags (n_lags,) bins, centres_ms (n_bases,))."""
    lo, hi = lag_range_ms
    lags = np.arange(int(np.floor(lo / bin_ms + 1e-9)) + 1,
                     int(np.floor(hi / bin_ms + 1e-9)) + 1)
    if lags.size == 0:
        raise ValueError("empty lag range")
    tau = lags * bin_ms
    if spacing == "log":
        warp, unwarp = (lambda x: np.log(x + log_offset_ms)), (lambda u: np.exp(u) - log_offset_ms)
    elif spacing == "uniform":
        warp, unwarp = (lambda x: x), (lambda u: u)
    else:
        raise ValueError(f"spacing must be 'log' or 'uniform', got {spacing!r}")
    u = warp(tau)
    if n_bases == 1:
        c, dc = np.array([u.mean()]), (u[-1] - u[0]) / 2 or 1.0
    else:
        c = np.linspace(u[0], u[-1], n_bases)
        dc = c[1] - c[0]
    arg = np.clip((u[:, None] - c[None, :]) * np.pi / (2 * dc), -np.pi, np.pi)
    B = 0.5 * (1 + np.cos(arg))
    return B, lags, unwarp(c)


def n_log_for(H_ms, n_at_150=6, per_doubling=2):
    """Log bumps for an own history of H ms: n_at_150 at 150 ms, + per_doubling for every
    doubling of H (150 -> 6, 200 -> 7, 300 -> 8, 450 -> 9, 600 -> 10). Log bumps widen with
    lag, so a longer history needs only a few more of them."""
    return max(1, int(round(n_at_150 + per_doubling * np.log2(H_ms / 150.0))))


def history_bases(H_ms, bin_ms, n_fine=2, n_log=6, log_offset_ms=5.0):
    """OWN-HISTORY basis up to H ms (script 03), strictly causal (lags >= 1 bin):
      - the first n_fine lags each get their OWN column (identity): refractoriness and
        bursts live in the first few ms and need full resolution;
      - the rest, (n_fine * bin, H], is covered by n_log log raised-cosine bumps
        (fine at short lags, coarse at long ones), so a long H costs few parameters.
    -> (B (n_lags, n_fine + n_log), lags (n_lags,) bins 1..H/bin)."""
    lags = np.arange(1, int(np.floor(H_ms / bin_ms + 1e-9)) + 1)
    B = np.zeros((lags.size, n_fine + n_log))
    B[np.arange(n_fine), np.arange(n_fine)] = 1.0
    Bl, ll, _ = raised_cosine_basis(n_log, (n_fine * bin_ms, H_ms), bin_ms, "log", log_offset_ms)
    B[ll - 1, n_fine:] = Bl
    return B, lags


def nuisance_mean(Z, cond):
    """TASK-LOCKED NUISANCE regressor (script 03, optional): for every trial and time t,
    the mean of Z at the SAME t over the OTHER trials of the same condition (leave-one-
    out). A time-varying intercept: it never uses the trial's own values, and only the
    other trials' values at time t (no temporal mixing). Z (n_trials, n_bins) -> same."""
    Z = np.asarray(Z, float)
    cond = np.asarray(cond)
    out = np.full_like(Z, np.nan)
    for c in np.unique(cond):
        idx = np.flatnonzero(cond == c)
        if idx.size < 2:
            continue
        Zi = Z[idx]
        fin = np.isfinite(Zi)
        s = np.where(fin, Zi, 0.0).sum(0, keepdims=True)
        n = fin.sum(0, keepdims=True)
        own = np.where(fin, Zi, 0.0)
        m = (s - own) / np.maximum(n - fin, 1)
        out[idx] = np.where((n - fin) > 0, m, np.nan)
    return out


def basis_features(z, B, lags):
    """ONE trial. lag_matrix(z, lags) @ B -> (n_bins, n_bases). NaN if any lag in the
    basis support is unavailable (so the row is dropped, never partially filled)."""
    return lag_matrix(z, lags) @ B


def event_regressors(t, event_t, B, lags):
    """ONE trial. Impulse at the bin containing event_t, seen through the causal basis:
    non-zero only in bins STRICTLY after the event bin. Before the event (or the trial
    start) the value is 0, not NaN -- there is no event there to be missing.
    t: (n_bins,) bin left edges; event_t: scalar (s) or NaN (no event)."""
    t = np.asarray(t, float)
    e = np.zeros(t.size)
    if np.isfinite(event_t):
        dt = t[1] - t[0]
        b = int(np.floor((event_t - t[0]) / dt + 1e-9))
        if 0 <= b < t.size:
            e[b] = 1.0
    return lag_matrix(e, lags, fill=0.0) @ B


def lagged_block(Z, fn, *args, **kw):
    """Apply a ONE-trial feature builder to every trial separately.
    Z (n_trials, n_bins) -> (n_trials, n_bins, p). This is the only way features are
    built, which is what guarantees that no lag crosses a trial boundary."""
    Z = np.asarray(Z, float)
    return np.stack([np.atleast_2d(fn(Z[i], *args, **kw).T).T for i in range(Z.shape[0])])


def build_design(blocks, y, in_window):
    """Stack per-trial features into a design.
      blocks    dict name -> (n_trials, n_bins, p) (from lagged_block / event_regressors)
      y         (n_trials, n_bins) target
      in_window (n_bins,) bool, rows eligible for the fit
    Row mask (computed ONCE, for this full model): in_window & y finite & every
    column finite.
    -> D = dict(X (n_rows, p), y (n_rows,), trial (n_rows,), tbin (n_rows,),
                groups {name: slice}, names, row_mask (n_trials, n_bins))."""
    y = np.asarray(y, float)
    n_tr, n_b = y.shape
    names, groups, cols, j = list(blocks), {}, [], 0
    for nm in names:
        A = np.asarray(blocks[nm], float)
        if A.ndim == 2:
            A = A[:, :, None]
        if A.shape[:2] != (n_tr, n_b):
            raise ValueError(f"block {nm!r} has shape {A.shape}, expected ({n_tr}, {n_b}, p)")
        groups[nm] = slice(j, j + A.shape[2]); j += A.shape[2]
        cols.append(A)
    F = np.concatenate(cols, axis=2) if cols else np.zeros((n_tr, n_b, 0))
    row_mask = (np.asarray(in_window, bool)[None, :] & np.isfinite(y)
                & np.isfinite(F).all(axis=2))
    tr, tb = np.nonzero(row_mask)
    return dict(X=F[tr, tb], y=y[tr, tb], trial=tr, tbin=tb, groups=groups,
                names=names, row_mask=row_mask)


def reduced(D, drop):
    """Same rows and y as D, without the column groups in `drop`."""
    drop = [drop] if isinstance(drop, str) else list(drop)
    missing = [d for d in drop if d not in D["groups"]]
    if missing:
        raise KeyError(f"unknown groups {missing}; have {list(D['groups'])}")
    keep_names = [n for n in D["names"] if n not in drop]
    idx, groups, j = [], {}, 0
    for n in keep_names:
        s = D["groups"][n]
        idx.extend(range(s.start, s.stop))
        groups[n] = slice(j, j + s.stop - s.start); j += s.stop - s.start
    return dict(D, X=D["X"][:, idx], groups=groups, names=keep_names)
