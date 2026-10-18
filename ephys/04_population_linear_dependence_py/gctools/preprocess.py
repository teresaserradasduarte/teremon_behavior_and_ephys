"""preprocess.py - binning, causal smoothing, zero-order-hold resampling and
ensemble-mean removal. EVERYTHING here is causal and per trial (time = last axis):
an output sample depends only on inputs at the same or earlier time.

TIME GRID. Anchored at t=0 (the align event): bin edges are integer multiples of
bin_ms, so 2.5 ms bins coincide with camera frames. Bin b covers [t_b, t_b + bin).
The analysis window is snapped INWARD to the grid; `pad_ms` of extra bins before
it hold the history that lagged regressors need (spike history is real data from
the same trial record; kinematics are NaN where they were not tracked).

KINEMATICS -> BINS: zero-order hold, queried at the bin's right edge minus epsilon,
i.e. a bin carries the LAST frame at or before its end. Never interpolated.
At 2.5 ms bins this is the identity; at coarser bins it DECIMATES (the other frames
of the bin are dropped -- mild aliasing for smooth med-lat position, no leakage).

SPIKES: counted on a fine 2.5 ms grid, optionally smoothed there with a CAUSAL gamma
kernel, then aggregated into bin_ms (counts summed, rates averaged over the sub-bins
of the bin -- no sample outside the bin is used).
Every function has an impulse test in tests/test_causality.py.
"""
import numpy as np
from scipy.signal import lfilter, lfilter_zi
from scipy.stats import gamma as _gamma
from . import config

FINE_MS = config.FRAME_MS      # 2.5 ms: spike grid before aggregation == frame period
_EPS_S  = 1e-7                 # s; query just inside a bin's right edge


# ---------------------------------------------------------------------------
# grid + binning
# ---------------------------------------------------------------------------
def bin_edges(window_ms, bin_ms, pad_ms=0.0):
    """-> (edges_s (n_bins+1,), t_s (n_bins,) bin LEFT edges, in_window (n_bins,) bool).
    Window snapped inward to multiples of bin_ms; pad_ms (rounded up to whole bins)
    of history bins precede it."""
    k0 = int(np.ceil(window_ms[0] / bin_ms - 1e-9))
    k1 = int(np.floor(window_ms[1] / bin_ms + 1e-9))
    kp = int(np.ceil(pad_ms / bin_ms - 1e-9))
    k  = np.arange(k0 - kp, k1 + 1)
    edges = k * bin_ms / 1000.0
    t     = edges[:-1]
    in_window = k[:-1] >= k0
    return edges, t, in_window


def bin_spikes(spikes, edges):
    """spikes: list (n_trials) of spike-time arrays (s) -> (n_trials, n_bins) int counts,
    bins [edges[i], edges[i+1])."""
    n_b = len(edges) - 1
    C = np.zeros((len(spikes), n_b), dtype=np.int32)
    for i, s in enumerate(spikes):
        s = np.asarray(s, float)
        idx = np.searchsorted(edges, s, side="right") - 1
        idx = idx[(idx >= 0) & (idx < n_b)]
        if idx.size:
            C[i] = np.bincount(idx, minlength=n_b)
    return C


def zoh_resample(t_src, X_src, t_query, max_hold_s=None):
    """Zero-order hold: value of the LATEST source sample with t_src <= t_query.
    NaN before the first sample, and (if max_hold_s) NaN when the latest sample is
    older than max_hold_s (i.e. beyond the end of the source coverage).
    X_src (..., n_src) -> (..., n_query)."""
    t_src, t_query = np.asarray(t_src, float), np.asarray(t_query, float)
    idx   = np.searchsorted(t_src, t_query, side="right") - 1
    valid = idx >= 0
    if max_hold_s is not None:
        valid &= (t_query - t_src[np.clip(idx, 0, None)]) <= max_hold_s + 1e-12
    out = np.asarray(X_src, float)[..., np.clip(idx, 0, None)]
    out[..., ~valid] = np.nan
    return out


def backward_diff(X, dt):
    """(x[t] - x[t-1]) / dt along the last axis; first sample NaN."""
    X = np.asarray(X, float)
    out = np.full_like(X, np.nan)
    out[..., 1:] = np.diff(X, axis=-1) / dt
    return out


# ---------------------------------------------------------------------------
# causal smoothing
# ---------------------------------------------------------------------------
def gamma_kernel(shape, scale_ms, bin_ms, q=0.999):
    """Causal gamma kernel on a grid of bin_ms. k[j] = gamma probability mass in
    [j*bin, (j+1)*bin) -- weight given to the input j bins in the PAST (j=0 is the
    current bin; nothing from the future). Normalised to sum 1.
    -> (k, mean_delay_ms) with mean_delay_ms = shape * scale_ms (continuous mean)."""
    t_max = _gamma.ppf(q, shape, scale=scale_ms)
    n = int(np.ceil(t_max / bin_ms)) + 1
    k = np.diff(_gamma.cdf(np.arange(n + 1) * bin_ms, shape, scale=scale_ms))
    return k / k.sum(), shape * scale_ms


def causal_filter(X, k, pre="hold"):
    """FIR filter along the last axis: y[t] = sum_j k[j] x[t-j], j >= 0.
    pre = "hold" -> the unseen past is the FIRST sample held constant (causal: x[0]
                    is never later than any t it feeds)
          "zero" -> unseen past = 0
          "nan"  -> the first len(k)-1 outputs are NaN (incomplete kernel)."""
    X = np.asarray(X, float)
    k = np.asarray(k, float)
    if pre == "zero" or len(k) == 1:
        return lfilter(k, [1.0], X, axis=-1)
    if pre == "hold":
        zi = lfilter_zi(k, [1.0]) * X[..., :1]
        y, _ = lfilter(k, [1.0], X, axis=-1, zi=zi)
        return y
    if pre == "nan":
        y = lfilter(k, [1.0], X, axis=-1)
        y[..., :len(k) - 1] = np.nan
        return y
    raise ValueError(f"pre must be 'hold', 'zero' or 'nan', got {pre!r}")


def remove_ensemble_mean(X, cond, loo=True):
    """Subtract the across-trial mean trajectory per condition. X (..., n_trials, n_bins),
    cond (n_trials,). loo=True uses the leave-one-out mean (a trial's own values
    never enter its own template). Contemporaneous only: no mixing across time."""
    X = np.asarray(X, float)
    out = np.array(X, copy=True)
    cond = np.asarray(cond)
    for c in np.unique(cond):
        idx = np.flatnonzero(cond == c)
        Xi  = X[..., idx, :]
        fin = np.isfinite(Xi)
        s   = np.nansum(Xi, axis=-2, keepdims=True)
        n   = fin.sum(axis=-2, keepdims=True)
        if loo:
            m = (s - np.where(fin, Xi, 0.0)) / np.maximum(n - fin, 1)
            m = np.where((n - fin) > 0, m, np.nan)
        else:
            m = s / np.maximum(n, 1)
        out[..., idx, :] = Xi - m
    return out


# ---------------------------------------------------------------------------
# session -> binned arrays
# ---------------------------------------------------------------------------
def bin_session(sess, bin_ms, window_ms=None, hist_pad_ms=0.0,
                kin_smooth="raw", fr_kernel=None, remove_mean=False):
    """Bin a session dict (io.load_session / sim.simulate_session) on a common grid.
      bin_ms       multiple of 2.5 ms
      hist_pad_ms  history bins before the window (must stay inside sess["spike_span"])
      kin_smooth   "raw" | "matched" (same causal gamma kernel as the spikes)
      fr_kernel    None (Y = counts) or (shape, scale_ms) -> Y = causal gamma rate (Hz)
      remove_mean  subtract the leave-one-out mean per trial condition from Y and X
    -> dict(t, edges, in_window, kin_avail, bin_ms, Y (n_neu, n_trials, n_bins),
            X (n_trials, n_bins; NaN where untracked), uids, regions, cond, trials,
            fr_kernel_mean_delay_ms, params)"""
    window_ms = window_ms or config.TRIAL_WINDOW_MS
    r = bin_ms / FINE_MS
    if abs(r - round(r)) > 1e-9 or round(r) < 1:
        raise ValueError(f"bin_ms must be a multiple of {FINE_MS} ms, got {bin_ms}")
    r = int(round(r))
    if kin_smooth == "matched" and fr_kernel is None:
        raise ValueError("kin_smooth='matched' needs fr_kernel=(shape, scale_ms)")

    edges, t, in_window = bin_edges(window_ms, bin_ms, hist_pad_ms)
    k, mean_delay = (gamma_kernel(*fr_kernel, FINE_MS) if fr_kernel is not None
                     else (None, 0.0))

    # fine spike grid: extra whole bins in front so the smoothing kernel is fully
    # formed from the first output bin on (real spikes, no edge assumption)
    n_extra = 0 if k is None else int(np.ceil((len(k) - 1) / r))
    k_first = int(round(edges[0] * 1000 / FINE_MS)) - n_extra * r
    k_last  = int(round(edges[-1] * 1000 / FINE_MS))
    fine_edges = np.arange(k_first, k_last + 1) * FINE_MS / 1000.0
    if fine_edges[0] < sess["spike_span"][0] - 1e-9:
        raise ValueError(f"history reaches {fine_edges[0]:.3f} s, before the spike span "
                         f"{sess['spike_span']} -- reduce hist_pad_ms / kernel")

    Y = []
    for n in sess["neurons"]:
        c = bin_spikes(n["spikes"], fine_edges).astype(float)
        if k is not None:
            c = causal_filter(c / (FINE_MS / 1000.0), k, pre="zero")
        c = c[:, n_extra * r:]                                  # drop kernel burn-in
        c = c.reshape(c.shape[0], -1, r)
        Y.append(c.mean(axis=2) if k is not None else c.sum(axis=2))
    Y = np.stack(Y) if Y else np.zeros((0, len(sess["trials"]), len(t)))
    if k is None:
        Y = Y.astype(np.int32)

    kin_t, X = sess["kin"]["t"], np.asarray(sess["kin"]["X"], float)
    if kin_smooth == "matched":
        X = causal_filter(X, k, pre="hold")
    elif kin_smooth != "raw":
        raise ValueError(f"kin_smooth must be 'raw' or 'matched', got {kin_smooth!r}")
    Xb = zoh_resample(kin_t, X, edges[1:] - _EPS_S, max_hold_s=config.FRAME_MS / 1000.0)

    cond = sess["trials"]["cond"].to_numpy()
    if remove_mean:
        Y  = remove_ensemble_mean(Y, cond)
        Xb = remove_ensemble_mean(Xb, cond)

    return dict(t=t, edges=edges, in_window=in_window,
                kin_avail=np.isfinite(Xb).all(axis=0), bin_ms=bin_ms,
                Y=Y, X=Xb, uids=[n["uid"] for n in sess["neurons"]],
                regions=[n["region"] for n in sess["neurons"]], cond=cond,
                trials=sess["trials"], fr_kernel_mean_delay_ms=mean_delay,
                params=dict(bin_ms=bin_ms, window_ms=tuple(window_ms),
                            hist_pad_ms=hist_pad_ms, kin_smooth=kin_smooth,
                            fr_kernel=fr_kernel, remove_mean=remove_mean))
