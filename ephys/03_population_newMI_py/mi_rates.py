"""mi_rates.py - three ways to turn spike times into a per-trial neural variable.

MODES
  "count" : spike COUNT in the integration window (the project default; exact,
            ISI-artefact-free, integer-valued -> Ross/plug-in friendly).
  "isi"   : instantaneous FR = 1/ISI sampled on a fine grid, then lightly smoothed
            per trial with a Gaussian. Ported to match the MATLAB pipeline exactly
            (auxiliar_functions/calculate_instant_FR_reach.m + gausskernel.m,
            called from s02_EGneurons_FR_all.m with bin_width=2 ms, sig_pk=20 ms).
  "conv"  : FR = spike train binned on the same fine grid and convolved with a
            wider Gaussian (default sigma 70 ms). No ISI step artefacts.

"isi" and "conv" return a RATE in Hz averaged over the integration window;
"count" returns an integer spike count. Downstream estimators must treat the
first two as CONTINUOUS (quantile-bin them for the plug-in estimator; Ross is
not applicable to them at all).

MATLAB parity notes for "isi":
  - 1/ISI is a step function: for a sample t falling in (s[j-1], s[j]] the rate is
    1/(s[j]-s[j-1]).
  - MATLAB reaches OUTSIDE the trial window (into the full spike train) to define
    the ISI before the first and after the last spike. Here st_reach is already
    trimmed to ~[-3,3] per trial, so those two edge segments are set to 0 (the same
    branch MATLAB takes when the spike is the first/last of the whole train).
    Irrelevant for the usual [-1,1] analysis range, which sits well inside the data.
  - Smoothing is a normalised Gaussian applied per trial. The raw trace is built on
    a PADDED grid and the convolution result is cropped back, reproducing MATLAB's
    win_padded + conv2(...,"valid") scheme (no edge roll-off inside the window).
"""
import numpy as np
from scipy.signal import fftconvolve

MODES = ("count", "isi", "conv")

# ---- defaults (match the MATLAB pipeline) ----------------------------------
RATE_BIN   = 0.002    # s, fine sampling grid            (MATLAB bin_width)
ISI_SIGMA  = 0.020    # s, light smoothing of 1/ISI      (MATLAB sig_pk)
CONV_SIGMA = 0.070    # s, wider kernel for "conv" mode  (70 ms, as requested)

MODE_LABEL = {
    "count": "spike count",
    "isi":   f"1/ISI  (Gauss {ISI_SIGMA*1000:.0f} ms)",
    "conv":  f"conv FR (Gauss {CONV_SIGMA*1000:.0f} ms)",
}
MODE_UNIT = {"count": "spikes", "isi": "Hz", "conv": "Hz"}


def gauss_kernel(sigma, bin_width=RATE_BIN):
    """Normalised Gaussian kernel, matching MATLAB gausskernel.m exactly:
    nbins = 2*round(sigma/bin_width * 9/2)  ->  +-4.5 sigma support."""
    nbins = 2 * int(round(sigma / bin_width * 9 / 2))
    bins  = np.arange(-nbins // 2 + 1, nbins // 2 + 1)
    x     = bins * bin_width
    pdf   = np.exp(-0.5 * (x / sigma) ** 2)
    return pdf / pdf.sum()


def _grid(t_lo, t_hi, bin_width, pad):
    lo = t_lo - pad
    n  = int(np.ceil((t_hi + pad - lo) / bin_width)) + 1
    return lo + np.arange(n) * bin_width


def isi_trace(spikes, t):
    """1/ISI step function sampled at t. Returns (n_trials, len(t)) in Hz.
    Sample t in (s[j-1], s[j]]  ->  1/(s[j]-s[j-1]);  before first / after last
    spike of the trial -> 0 (see module docstring)."""
    out = np.zeros((len(spikes), t.size))
    for i, s in enumerate(spikes):
        s = np.sort(np.asarray(s, float))
        if s.size < 2:
            continue
        idx   = np.searchsorted(s, t, side="left")   # s[idx-1] < t <= s[idx]
        inner = (idx > 0) & (idx < s.size)
        j     = idx[inner]
        isi   = s[j] - s[j - 1]
        out[i, inner] = np.where(isi > 0, 1.0 / np.maximum(isi, 1e-12), 0.0)
    return out


def binned_rate(spikes, t, bin_width):
    """Spike train binned on t, expressed in Hz. Returns (n_trials, len(t)).
    One bincount over all trials at once -- ~2x faster than per-trial histogram
    and numerically identical."""
    edges = np.append(t - bin_width / 2, t[-1] + bin_width / 2)
    n_tr, n_b = len(spikes), t.size
    rows, cols = [], []
    for i, s in enumerate(spikes):
        s = np.asarray(s, float)
        s = s[(s >= edges[0]) & (s < edges[-1])]
        if s.size:
            cols.append(np.searchsorted(edges, s, side="right") - 1)
            rows.append(np.full(s.size, i))
    if not rows:
        return np.zeros((n_tr, n_b))
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    flat = np.bincount(rows * n_b + cols, minlength=n_tr * n_b)
    return flat.reshape(n_tr, n_b) / bin_width


def _smooth_rows(R, k):
    """Gaussian smoothing of every trial at once (FFT convolution along time).
    Matches a per-trial np.convolve(..., "same") to ~1e-14."""
    return fftconvolve(R, k[None, :], mode="same", axes=1)


def rate_traces(spikes, mode, t_lo, t_hi, bin_width=RATE_BIN,
                isi_sigma=ISI_SIGMA, conv_sigma=CONV_SIGMA):
    """Continuous per-trial trace on [t_lo, t_hi].  Returns (t, R) with R
    (n_trials, len(t)).  "count" returns the raw binned rate (unsmoothed) -- it
    exists for VISUALISATION only; the count mode's analysis variable is the exact
    integer count from windowed()."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")

    sigma = {"count": 0.0, "isi": isi_sigma, "conv": conv_sigma}[mode]
    k     = gauss_kernel(sigma, bin_width) if sigma > 0 else None
    pad   = (len(k) // 2 + 1) * bin_width if k is not None else 0.0

    t_pad = _grid(t_lo, t_hi, bin_width, pad)
    R = isi_trace(spikes, t_pad) if mode == "isi" else binned_rate(spikes, t_pad, bin_width)
    if k is not None:
        R = _smooth_rows(R, k)

    keep = (t_pad >= t_lo - bin_width / 2) & (t_pad <= t_hi + bin_width / 2)
    return t_pad[keep], R[:, keep]


def windowed(spikes, win_starts, win_width, mode, bin_width=RATE_BIN,
             isi_sigma=ISI_SIGMA, conv_sigma=CONV_SIGMA):
    """Per-trial neural variable for every sliding window.
    Returns (n_trials, n_wins):
      "count" -> int spike count in [start, start+width)   (exact, no binning)
      "isi"/"conv" -> mean rate (Hz) over the window."""
    win_starts = np.asarray(win_starts, float)

    if mode == "count":
        hi = win_starts + win_width
        C  = np.zeros((len(spikes), win_starts.size), dtype=np.int32)
        for ti, spk in enumerate(spikes):
            if len(spk) == 0:
                continue
            s = np.sort(np.asarray(spk, float))
            C[ti] = np.searchsorted(s, hi) - np.searchsorted(s, win_starts)
        return C

    t, R = rate_traces(spikes, mode, win_starts[0], win_starts[-1] + win_width,
                       bin_width=bin_width, isi_sigma=isi_sigma, conv_sigma=conv_sigma)
    csum = np.concatenate([np.zeros((R.shape[0], 1)), np.cumsum(R, axis=1)], axis=1)
    # Half-bin tolerance: a sample belongs to [start, start+width) when its CENTRE
    # does. Without it, float noise in win_starts (np.arange never lands exactly on
    # 0.0) shifts searchsorted by a whole sample, so window sample counts vary
    # (51 vs 50) and a single-window call disagrees with the full-landscape call.
    tol = bin_width / 2
    i0 = np.searchsorted(t, win_starts - tol, side="left")
    i1 = np.maximum(np.searchsorted(t, win_starts + win_width - tol, side="left"),
                    i0 + 1)
    return (csum[:, i1] - csum[:, i0]) / (i1 - i0)[None, :]


def is_discrete(mode):
    """True only for "count" -- the one mode Ross / unique-value binning applies to."""
    return mode == "count"


# ---------------------------------------------------------------------------
# Cache helpers for the WINDOWED variable
# ---------------------------------------------------------------------------
# The windowed neural variable depends on BOTH the smoothing params and the
# sliding window, so the window params belong in the cache key -- they are put in
# the FILENAME so different window configs coexist instead of overwriting one
# another. Only "isi"/"conv" are worth caching; "count" is a searchsorted away.
def mode_sigma(mode, isi_sigma=ISI_SIGMA, conv_sigma=CONV_SIGMA):
    """The smoothing sigma that actually affects `mode` (0 for "count")."""
    return {"isi": isi_sigma, "conv": conv_sigma}.get(mode, 0.0)


def cache_name(mode, win_width, win_step, win_range,
               rate_bin=RATE_BIN, isi_sigma=ISI_SIGMA, conv_sigma=CONV_SIGMA):
    """Filename for a cached windowed variable.

    EVERY parameter the values depend on is in the NAME -- the window config AND
    the smoothing that applies to this mode -- so caches for different settings
    coexist instead of silently overwriting one another. Only the sigma relevant
    to `mode` is included: "isi" does not depend on CONV_SIGMA and vice versa."""
    sig = mode_sigma(mode, isi_sigma, conv_sigma)
    return (f"neuralvar_{mode}_w{win_width*1000:.0f}ms_s{win_step*1000:.0f}ms"
            f"_r{win_range[0]:g}to{win_range[1]:g}"
            f"_sig{sig*1000:.0f}ms_bin{rate_bin*1000:.0f}ms.pkl")


def cache_matches(res, uids, win_starts, win_width,
                  rate_bin=RATE_BIN, isi_sigma=ISI_SIGMA, conv_sigma=CONV_SIGMA):
    """True only if the cache was built from the same neurons, the same window and
    the same smoothing FOR THIS MODE. Deliberately ignores the other mode's sigma:
    an "isi" cache must not be thrown away because CONV_SIGMA changed."""
    try:
        want = mode_sigma(res["mode"], isi_sigma, conv_sigma)
        have = mode_sigma(res["mode"], res["ISI_SIGMA"], res["CONV_SIGMA"])
        return (list(res["uids"]) == list(uids)
                and len(res["win_starts"]) == len(win_starts)
                and np.allclose(res["win_starts"], win_starts)
                and np.isclose(res["WIN_WIDTH"], win_width)
                and np.isclose(res["RATE_BIN"], rate_bin)
                and np.isclose(have, want))
    except KeyError:
        return False
