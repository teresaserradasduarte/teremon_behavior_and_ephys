"""ct_features.py - waveform and spike-train features. No file access here."""
import numpy as np
from scipy.interpolate import CubicSpline

UPS = 10            # waveform upsampling factor (30 kHz -> 300 kHz, 3.3 us steps)


# ── correlograms ─────────────────────────────────────────────────────────────
def ccg(ref, tgt, bin_s=0.001, win_s=0.05, chunk=200_000):
    """Spikes of `tgt` around each spike of `ref`: counts per lag bin.
    -> (counts, lag centres in s). lag > 0 = tgt AFTER ref. Both sorted, seconds."""
    nb = int(round(win_s / bin_s))
    edges = np.arange(-nb, nb + 1) * bin_s
    counts = np.zeros(2 * nb, np.int64)
    ref = np.asarray(ref, float); tgt = np.asarray(tgt, float)
    for i in range(0, ref.size, chunk):
        r = ref[i:i + chunk]
        lo = np.searchsorted(tgt, r - win_s, "left")
        hi = np.searchsorted(tgt, r + win_s, "right")
        n = hi - lo
        if n.sum() == 0: continue
        idx = np.repeat(lo, n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
        counts += np.histogram(tgt[idx] - np.repeat(r, n), edges)[0]
    return counts, (edges[:-1] + edges[1:]) / 2


def acg(st, bin_s=0.001, win_s=0.5):
    """Autocorrelogram in spikes/s, one side only (lags bin_s/2 ... win_s).
    -> (rate per lag, lag centres in s)."""
    st = np.asarray(st, float)
    c, lags = ccg(st, st, bin_s, win_s)
    nb = c.size // 2
    pos = c[nb:].astype(float)
    # each spike paired with itself has lag exactly 0 -> lands in the first positive bin
    pos[0] = max(pos[0] - st.size, 0.0)
    return pos / max(st.size, 1) / bin_s, lags[nb:]


# ── spike-train features ─────────────────────────────────────────────────────
def train_features(st, acg_rate, acg_lags, long_isi=2.0):
    """Features of one spike train (st sorted, s). acg_rate: acg() with 1 ms bins, 0.5 s."""
    st = np.asarray(st, float)
    out = dict(n_spikes=int(st.size), fr=np.nan, cv2=np.nan, isi_median_ms=np.nan,
               prop_long_isi=np.nan, post_spike_supp_ms=np.nan, refractory_ms=np.nan,
               burst_index=np.nan, acg_peak_ms=np.nan)
    if st.size < 20: return out
    dur = st[-1] - st[0]
    isi = np.diff(st)
    out["fr"] = st.size / dur
    out["isi_median_ms"] = 1e3 * np.median(isi)
    out["cv2"] = float(np.mean(2 * np.abs(isi[1:] - isi[:-1]) / (isi[1:] + isi[:-1] + 1e-12)))
    # bombcell propLongISI: fraction of the recording spent in ISIs longer than 2 s
    out["prop_long_isi"] = float(isi[isi > long_isi].sum() / dur)
    # bombcell postSpikeSuppression: first lag at which the ACG is back at its
    # baseline (mean of the 100-400 ms lags). bombcell's ACG is 1 ms / +-500 ms and it
    # searches indices 500:1000, i.e. from lag 0 -- same thing on a one-sided ACG.
    base = np.mean(acg_rate[(acg_lags >= 0.1) & (acg_lags < 0.4)])
    if base > 0:
        k = np.flatnonzero(acg_rate >= base)
        out["post_spike_supp_ms"] = 1e3 * acg_lags[k[0]] if k.size else np.nan
        # refractory period: first lag where the ACG passes 10% of its baseline
        k = np.flatnonzero(acg_rate >= 0.1 * base)
        out["refractory_ms"] = 1e3 * acg_lags[k[0]] if k.size else np.nan
        # burst index: ACG at 1-10 ms relative to the baseline (>1 = bursty)
        out["burst_index"] = float(np.mean(acg_rate[acg_lags < 0.010]) / base)
    out["acg_peak_ms"] = 1e3 * acg_lags[int(np.argmax(acg_rate))]
    return out


# ── waveform features ────────────────────────────────────────────────────────
def peak_channel(W):
    """W: (n_channels, n_samples) template -> index of the largest-amplitude channel."""
    return int(np.argmax(W.max(axis=1) - W.min(axis=1)))


def upsample(w, ups=UPS):
    x = np.arange(w.size)
    xu = np.arange(0, w.size - 1 + 1e-9, 1 / ups)
    return CubicSpline(x, w)(xu)


def waveform_features(w, fs=30000.0):
    """Features of ONE channel's mean waveform (1-D, any units).
    duration: bombcell's waveformDuration_peakTrough -- time between the main
    deflection and the largest opposite deflection that FOLLOWS it."""
    out = dict(wf_duration_us=np.nan, wf_halfwidth_us=np.nan, wf_pt_ratio=np.nan,
               wf_polarity=0, wf_repol_slope=np.nan, wf_amp=np.nan)
    w = np.asarray(w, float)
    if not np.any(np.isfinite(w)) or np.ptp(w) == 0: return out
    w = w - np.median(w[: max(5, w.size // 8)])            # baseline = the first samples
    u = upsample(w); dt_us = 1e6 / fs / UPS
    i0 = int(np.argmax(np.abs(u)))
    sgn = -1.0 if u[i0] < 0 else 1.0
    v = -sgn * u                                            # main deflection is now a TROUGH
    i1 = i0 + int(np.argmax(v[i0:]))                        # the following peak
    out["wf_polarity"] = int(sgn)
    out["wf_amp"] = float(np.ptp(w))
    out["wf_duration_us"] = (i1 - i0) * dt_us
    half = v[i0] / 2
    a = i0
    while a > 0 and v[a] < half: a -= 1
    b = i0
    while b < v.size - 1 and v[b] < half: b += 1
    out["wf_halfwidth_us"] = (b - a) * dt_us
    out["wf_pt_ratio"] = float(v[i1] / -v[i0]) if v[i0] != 0 else np.nan
    if i1 > i0 + 1:                                         # steepest return from the trough
        out["wf_repol_slope"] = float(np.max(np.diff(v[i0:i1 + 1])) / -v[i0] / (dt_us * 1e-3))
    return out


def norm_waveform(w, n=82):
    """Peak-channel waveform for embedding: baseline-subtracted, trough-aligned,
    flipped so the main deflection is negative, scaled to unit amplitude,
    resampled to `n` points (so 61-sample KS4 and 82-sample templates are comparable).
    The window is -1.0 .. +1.7 ms around the main deflection."""
    w = np.asarray(w, float)
    if not np.any(np.isfinite(w)) or np.ptp(w) == 0: return np.full(n, np.nan)
    w = w - np.median(w[: max(5, w.size // 8)])
    u = upsample(w)
    i0 = int(np.argmax(np.abs(u)))
    u = u / -u[i0]                                          # trough = -1
    t = (np.arange(u.size) - i0) / (30.0 * UPS)             # ms
    return np.interp(np.linspace(-1.0, 1.7, n), t, u, left=0.0, right=0.0)


# ── CS-triggered pause ───────────────────────────────────────────────────────
def pause_stats(counts, lags, n_ref, bin_s, base_win=(-0.05, -0.01), pause_win=(0.003, 0.008)):
    """CS-triggered histogram of another unit -> how strongly it pauses after the CS.
    ratio      rate in pause_win / rate in base_win   (0 = complete pause)
    mirror     the same in the mirrored window BEFORE the CS (-8..-3 ms)
    asym       mirror - ratio: a real pause is one-sided (~1); a sorter collision
               dip or a shared refractory period is symmetric (~0)
    z          Poisson z-score of the drop
    pause_ms   how long after the CS the rate stays below half of the baseline"""
    rate = counts / max(n_ref, 1) / bin_s
    b = (lags >= base_win[0]) & (lags < base_win[1])
    p = (lags >= pause_win[0]) & (lags < pause_win[1])
    q = (lags >= -pause_win[1]) & (lags < -pause_win[0])
    base = rate[b].mean()
    out = dict(base_hz=float(base), ratio=np.nan, mirror=np.nan, asym=np.nan, z=np.nan,
               pause_ms=0.0, win_hz=float(rate[p].mean()))
    if base <= 0: return out
    exp = base * bin_s * n_ref * p.sum()                    # expected count in the window
    out["ratio"] = float(rate[p].mean() / base)
    out["mirror"] = float(rate[q].mean() / base)
    out["asym"] = out["mirror"] - out["ratio"]
    out["z"] = float((counts[p].sum() - exp) / np.sqrt(max(exp, 1e-9)))
    post = np.flatnonzero(lags >= pause_win[0])
    sm = np.convolve(rate, np.ones(3) / 3, "same")          # 3-bin smoothing
    k = next((j for j in post if sm[j] >= 0.5 * base), post[-1])
    out["pause_ms"] = float(1e3 * (lags[k] - pause_win[0]))
    return out
