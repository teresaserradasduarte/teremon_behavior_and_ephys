"""Causality of every predictor builder (CLAUDE.md rule 1).

Two complementary checks:
  IMPULSE        a unit impulse at sample n must leave every output before n exactly
                 zero (or NaN) -- and, for strictly lagged builders, also AT n.
  FUTURE-PERTURB randomise all input after t*; outputs at <= t* must be bit-identical
                 (for strictly lagged builders, outputs at t*+1 too, since they read
                 only up to t*).
NEGATIVE CONTROL: a centred Gaussian (filtfilt-like) must FAIL the future-perturb
check, which proves the check can detect leakage.
"""
from _harness import run
import numpy as np
from scipy.ndimage import gaussian_filter1d
from gctools import preprocess as pp, design as ds, sim, config

RNG = np.random.default_rng(0)
N = 200


def _impulse(n=N, at=80):
    z = np.zeros(n); z[at] = 1.0
    return z, at


def _zero_before(y, at, strict):
    """y: (n,) or (n, p) output over time (axis 0). strict -> also zero at `at`."""
    y = np.asarray(y, float)
    y = y if y.ndim > 1 else y[:, None]
    stop = at + 1 if strict else at
    pre = y[:stop]
    return np.all((pre == 0) | np.isnan(pre))


def _future_invariant(builder, strict, n=N, t_star=90, reps=5):
    """builder: z (n,) -> (n,) or (n, p). Randomise z[t*+1:], compare outputs at <= t*
    (strict: <= t*+1)."""
    z = RNG.normal(size=n)
    base = np.asarray(builder(z), float)
    upto = t_star + 2 if strict else t_star + 1
    for _ in range(reps):
        z2 = z.copy(); z2[t_star + 1:] = RNG.normal(size=n - t_star - 1)
        out = np.asarray(builder(z2), float)
        if not np.array_equal(base[:upto], out[:upto], equal_nan=True):
            return False
    return True


# ---------------------------------------------------------------------------
def test_gamma_kernel_causal_and_delay():
    for shape, scale, bin_ms in [(2.0, 10.0, 2.5), (4.0, 5.0, 2.5), (3.0, 20.0, 10.0)]:
        k, mean_delay = pp.gamma_kernel(shape, scale, bin_ms)
        z, at = _impulse(400, 50)
        y = pp.causal_filter(z, k, pre="zero")
        assert _zero_before(y, at, strict=False), "gamma filter responds before the impulse"
        # discrete mean lag: mass of [j, j+1) bins sits at lag j -> E = mean - bin/2
        emp = np.sum(np.arange(k.size) * k) * bin_ms
        print(f"  gamma a={shape} th={scale} bin={bin_ms}: mean delay {mean_delay:.2f} ms, "
              f"discrete {emp + bin_ms/2:.2f} ms (+bin/2)")
        assert abs(emp + bin_ms / 2 - mean_delay) <= bin_ms / 2, "kernel mean delay off"
        assert abs(y.sum() - 1) < 1e-9


def test_causal_filter_future_invariance():
    k, _ = pp.gamma_kernel(3.0, 8.0, 2.5)
    for pre in ("hold", "zero", "nan"):
        assert _future_invariant(lambda z: pp.causal_filter(z, k, pre=pre), strict=False), pre


def test_negative_control_centred_gaussian_fails():
    leaky = lambda z: gaussian_filter1d(z, sigma=5)            # centred -> uses the future
    assert not _future_invariant(leaky, strict=False), "the check failed to detect leakage"
    z, at = _impulse()
    assert not _zero_before(leaky(z), at, strict=False)
    print("  centred Gaussian correctly flagged as non-causal")


def test_zoh_resample():
    # Downsampling by ZOH keeps only the LAST frame of each bin, so a one-frame
    # impulse can vanish (aliasing, not leakage). The causal probe is therefore a
    # STEP: it must appear first in the bin that contains it, never earlier.
    t_src = np.arange(N) * 0.0025
    edges = np.arange(0, N * 0.0025 + 1e-12, 0.010)
    for at in (80, 81, 82, 83):                            # every phase within a 10 ms bin
        z = (np.arange(N) >= at).astype(float)
        y = pp.zoh_resample(t_src, z, edges[1:] - 1e-7, max_hold_s=0.0025)
        first = np.flatnonzero(y > 0)[0]
        b_true = int(np.floor(t_src[at] / 0.010 + 1e-9))
        print(f"  step at t={t_src[at]*1000:.1f} ms appears first in bin "
              f"[{edges[first]*1000:.0f}, {edges[first+1]*1000:.0f}) ms")
        assert first == b_true, "step appears in the wrong bin"
        assert np.all(y[:first] == 0)
    # future invariance on an arbitrary query grid
    q = np.linspace(0, t_src[-1], 137)
    t_star = 90
    base = pp.zoh_resample(t_src, RNG.normal(size=N), q)
    z = RNG.normal(size=N); b = pp.zoh_resample(t_src, z, q)
    z2 = z.copy(); z2[t_star + 1:] = 0
    b2 = pp.zoh_resample(t_src, z2, q)
    m = q <= t_src[t_star] + 1e-12
    assert np.array_equal(b[m], b2[m])
    # identity at 2.5 ms (frames == bins)
    e25 = np.arange(N + 1) * 0.0025
    assert np.array_equal(pp.zoh_resample(t_src, z, e25[1:] - 1e-7), z)
    # NaN before the first sample and beyond the coverage
    y = pp.zoh_resample(t_src, z, np.array([-0.001, t_src[-1] + 0.01]), max_hold_s=0.0025)
    assert np.all(np.isnan(y))


def test_backward_diff():
    z, at = _impulse()
    y = pp.backward_diff(z, 0.0025)
    assert _zero_before(y, at, strict=False)
    assert set(np.flatnonzero((y != 0) & ~np.isnan(y))) == {at, at + 1}
    assert _future_invariant(lambda s: pp.backward_diff(s, 0.0025), strict=False)


def test_lag_matrix_strict():
    lags = np.arange(1, 13)
    z, at = _impulse()
    y = ds.lag_matrix(z, lags)
    assert _zero_before(y, at, strict=True), "lag matrix sees the current/future sample"
    assert _future_invariant(lambda s: ds.lag_matrix(s, lags), strict=True)
    try:
        ds.lag_matrix(z, [0, 1]); raise AssertionError("lag 0 was accepted")
    except ValueError:
        pass


def test_raised_cosine_basis_strict():
    for spacing in ("log", "uniform"):
        B, lags, c = ds.raised_cosine_basis(6, (0, 150), 2.5, spacing=spacing)
        assert lags[0] >= 1
        z, at = _impulse(400, 100)
        y = ds.basis_features(z, B, lags)
        assert _zero_before(y, at, strict=True), spacing
        assert _future_invariant(lambda s: ds.basis_features(s, B, lags), strict=True, n=400,
                                 t_star=200)
        # impulse response of the features = the basis itself, shifted by the lag
        resp = y[at + lags]
        assert np.allclose(resp, B)
        print(f"  {spacing}: {len(lags)} lags ({lags[0]*2.5:.1f}-{lags[-1]*2.5:.1f} ms), "
              f"centres {np.round(c, 1).tolist()} ms")


def test_history_bases_strict():
    for H in (150, 200, 300):
        B, lags = ds.history_bases(H, 2.5)
        assert lags[0] == 1 and lags[-1] * 2.5 == H
        assert np.allclose(B[0, :2], [1, 0]) and np.allclose(B[1, :2], [0, 1])
        assert np.all(B[:2, 2:] == 0) and np.all(B[2:, :2] == 0), "fine and log parts overlap"
        assert np.all(B.sum(axis=1) > 0), f"H={H}: some lag is not covered"
        z, at = _impulse(400, 120)
        y = ds.basis_features(z, B, lags)
        assert _zero_before(y, at, strict=True)
        assert _future_invariant(lambda s: ds.basis_features(s, B, lags), strict=True, n=400, t_star=200)
    print(f"  H=300: {B.shape[1]} columns for {lags.size} lags (2 single-lag + 6 log bumps)")


def test_history_bases_scaled_strict():
    """Own-history bases with n_log_for(H) bumps (script 03, scale_bases): still strictly
    causal and covering every lag up to H."""
    for H in (150, 200, 300, 450, 600):
        B, lags = ds.history_bases(H, 2.5, n_log=ds.n_log_for(H))
        assert B.shape[1] == 2 + ds.n_log_for(H)
        assert lags[0] == 1 and lags[-1] * 2.5 == H and np.all(B.sum(axis=1) > 0)
        z, at = _impulse(600, 300)
        assert _zero_before(ds.basis_features(z, B, lags), at, strict=True)
        assert _future_invariant(lambda s: ds.basis_features(s, B, lags), strict=True, n=600, t_star=300)
    print("  log bumps per H:", {H: ds.n_log_for(H) for H in (150, 200, 300, 450, 600)})


def test_nuisance_mean_never_uses_own_trial_or_future():
    Z = RNG.normal(size=(30, N)); cond = np.repeat([1, 2, 3], 10)
    base = ds.nuisance_mean(Z, cond)
    Z2 = Z.copy(); Z2[4] += 100.0                          # change trial 4 only
    d = ds.nuisance_mean(Z2, cond) - base
    assert np.allclose(d[4], 0), "trial's own values entered its nuisance regressor"
    assert np.allclose(d[[i for i in range(10) if i != 4]], 100.0 / 9)   # others in its condition
    assert np.allclose(d[10:], 0)                                        # other conditions untouched
    Z3 = Z.copy(); Z3[:, 101:] = RNG.normal(size=(30, N - 101))          # change the future of ALL
    assert np.array_equal(ds.nuisance_mean(Z3, cond)[:, :101], base[:, :101])


def test_window_means_strict():
    lag_sets = ds.ms_windows_to_lags([(0, 10), (10, 20), (20, 40)], 2.5)
    assert lag_sets[0].tolist() == [1, 2, 3, 4]
    z, at = _impulse()
    assert _zero_before(ds.window_means(z, lag_sets), at, strict=True)
    assert _future_invariant(lambda s: ds.window_means(s, lag_sets), strict=True)


def test_event_regressors_only_after_event():
    B, lags, _ = ds.raised_cosine_basis(4, (0, 100), 2.5)
    t = (np.arange(200) - 100) * 0.0025
    y = ds.event_regressors(t, 0.0, B, lags)
    ev_bin = 100
    assert np.all(y[:ev_bin + 1] == 0), "event kernel active at/before the event"
    assert np.any(y[ev_bin + 1:] != 0)
    assert np.all(ds.event_regressors(t, np.nan, B, lags) == 0)


def test_remove_ensemble_mean_no_temporal_mixing():
    X = RNG.normal(size=(30, N)); cond = np.repeat([1, 2, 3], 10)
    base = pp.remove_ensemble_mean(X, cond)
    X2 = X.copy(); X2[:, 101:] = RNG.normal(size=(30, N - 101))
    out = pp.remove_ensemble_mean(X2, cond)
    assert np.array_equal(base[:, :101], out[:, :101])
    # leave-one-out: a trial's own value never enters its own template
    X3 = X.copy(); X3[0, 50] += 100.0
    d = pp.remove_ensemble_mean(X3, cond) - base
    assert np.isclose(d[0, 50], 100.0) and np.allclose(d[1:10, 50], -100.0 / 9)
    assert np.allclose(d[10:, 50], 0) and np.count_nonzero(d[:, np.r_[0:50, 51:N]]) == 0


def test_sim_drive_causal():
    k, k_t = sim.coupling_kernel(50.0, 5.0, 0.5)
    assert k[0] == 0.0, "coupling kernel has a contemporaneous term"
    z, at = _impulse(1000, 300)
    y = pp.causal_filter(z, k, pre="hold")
    assert _zero_before(y, at, strict=True)
    mode = k_t[np.argmax(k)]
    print(f"  coupling kernel mode {mode:.1f} ms (target 50)")
    assert abs(mode - 50.0) <= 1.0


def test_sim_spike_to_kin_drive_causal():
    """y_to_x: a spike at time s may change the kinematics only at frames after s."""
    k, _ = sim.coupling_kernel(50.0, 5.0, config.FRAME_MS)
    t_f = (np.arange(config.N_FRAMES) - config.ALIGN_FRAME) / config.FPS
    for s in (-0.2013, -0.2000, -0.19999):            # mid-frame, on a frame, just after
        D = sim.spike_drive([np.array([s])], t_f, k, t_start=-0.675)[0]
        first = np.flatnonzero(D != 0)[0]
        print(f"  spike at {s*1000:.2f} ms -> first kinematic change at frame "
              f"{t_f[first]*1000:.1f} ms (peak {t_f[np.argmax(D)]*1000 - s*1000:.1f} ms later)")
        assert t_f[first] > s and np.all(D[t_f <= s] == 0)
    # future perturbation: spikes after t* never change frames <= t*
    rng = np.random.default_rng(0)
    base = np.sort(rng.uniform(-0.6, 0.12, 40))
    t_star = -0.1
    a = sim.spike_drive([base], t_f, k, -0.675)[0]
    extra = np.sort(np.r_[base[base <= t_star], rng.uniform(t_star + 1e-6, 0.12, 30)])
    b = sim.spike_drive([extra], t_f, k, -0.675)[0]
    assert np.array_equal(a[t_f <= t_star], b[t_f <= t_star])


def test_bin_session_spikes_and_kin_causal():
    """End-to-end: one spike and one kinematic STEP (see test_zoh_resample for why a
    step); the binned outputs must not see them before the bin that contains them
    (and 'matched' smoothing must not either)."""
    S = sim.synthetic_template(n_trials=2, seed=1)
    f_imp = config.ALIGN_FRAME - 39                      # frame at -97.5 ms (mid-bin at 10 ms)
    S["kin"]["X"][:] = (np.arange(config.N_FRAMES) >= f_imp).astype(float)
    t_spk = -0.2013
    S["neurons"] = [dict(uid="T|x|y|0", region="SIM", phyID=0, spikes=[np.array([t_spk])] * 2)]
    for bin_ms in (2.5, 10.0):
        for kin_smooth, frk in (("raw", None), ("matched", (3.0, 5.0))):
            B = pp.bin_session(S, bin_ms, hist_pad_ms=100, kin_smooth=kin_smooth, fr_kernel=frk)
            t = B["t"]
            b_spk = np.flatnonzero((t <= t_spk) & (t + bin_ms / 1000 > t_spk))[0]
            y = B["Y"][0, 0]
            assert np.all(y[:b_spk] == 0) and y[b_spk] > 0, (bin_ms, kin_smooth)
            t_imp = S["kin"]["t"][f_imp]
            b_kin = np.flatnonzero((t <= t_imp + 1e-9) & (t + bin_ms / 1000 > t_imp + 1e-9))[0]
            x = np.nan_to_num(B["X"][0])
            assert np.all(x[:b_kin] == 0), (bin_ms, kin_smooth)
            print(f"  bin {bin_ms} ms, {kin_smooth}: spike bin {t[b_spk]*1000:.1f} ms, "
                  f"kin first non-zero {t[np.flatnonzero(x)[0]]*1000:.1f} ms "
                  f"(step at {t_imp*1000:.1f} ms)")


if __name__ == "__main__":
    run(globals())
