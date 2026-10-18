"""Script 03 model (gctools.kernels): shuffle reproducibility, CV gain under the null,
direction + lag recovery, adequacy-test calibration, lag bootstrap."""
from _harness import run
import numpy as np
from gctools import sim, kernels, nulls, preprocess as pp

LAG = 50.0


def _binned(scenario, n_tr=200, seed=0, kin="white", **kw):
    T = sim.synthetic_template(n_trials=n_tr, seed=seed, kin=kin)
    spec = kernels.make_spec(h_list=(150, 200), n_folds=4)
    pad = kernels.pad_needed_ms(spec)
    S = sim.simulate_session(T, scenario, seed=seed + 1, lag_ms=LAG, hist_pad_ms=pad + 100, **kw)
    B = pp.bin_session(S, 2.5, hist_pad_ms=pad)
    return S, B, spec


def test_kernelnull_reproduces_stored_null():
    S, B, _ = _binned("x_to_y", n_tr=60, rate_hz=(20.0,), n_per_rate=1, coupling=1.0)
    perms = nulls.make_perms(B["cond"], 15, seed=2)
    for spec in (kernels.make_spec(h_list=(150, 200), n_folds=3),
                 kernels.make_spec(hist_mode="fixed", event=True, nuisance=True, n_folds=3)):
        r = kernels.kernel_neuron(B["Y"][0], B["X"], B["in_window"], B["cond"], B["t"], spec, perms)
        for d in ("xy", "yx"):
            kn = kernels.KernelNull(B["Y"][0], B["X"], B["in_window"], B["cond"], B["t"], spec, r["H"], d)
            assert np.array_equal(kn.stats(perms), r[d]["null"]), (spec["hist_mode"], d)
    print(f"  chosen H = {r['H']} ms; nulls reproduced bit-for-bit (cv and fixed+event+nuisance)")


def test_null_cv_gain_and_adequacy_calibrated():
    S, B, spec = _binned("null_independent", n_tr=150, rate_hz=(10.0, 40.0), n_per_rate=6)
    g_xy, g_yx, padq = [], [], []
    for i in range(len(S["neurons"])):
        r = kernels.kernel_neuron(B["Y"][i], B["X"], B["in_window"], B["cond"], B["t"], spec)
        g_xy.append(r["xy"]["cv_gain_bits"]); g_yx.append(r["yx"]["cv_gain_r2"])
        padq.append(r["xy"]["adequacy_p"])
    print(f"  null CV gain: kin->neuron median {np.median(g_xy):+.5f} bits/spike, "
          f"neuron->kin median {np.median(g_yx):+.5f}; adequacy p<0.05 in {np.mean(np.array(padq)<0.05):.0%}")
    assert np.median(g_xy) <= 0.001 and np.median(g_yx) <= 0.001
    assert np.mean(np.array(padq) < 0.05) <= 0.25        # n=12: loose bound around 5%


def test_x_to_y_detected_with_right_lag():
    S, B, spec = _binned("x_to_y", n_tr=200, rate_hz=(25.0,), n_per_rate=3, coupling=1.0, kernel_shape=8.0)
    perms = nulls.make_perms(B["cond"], 49, seed=1)
    for i in range(3):
        r = kernels.kernel_neuron(B["Y"][i], B["X"], B["in_window"], B["cond"], B["t"], spec, perms)
        print(f"  neuron {i}: p_shuf {r['xy']['p_shuf']:.3f}, CV gain {r['xy']['cv_gain_bits']:.3f} bits/spike, "
              f"peak {r['xy']['peak_ms']:.1f} ms (sign {r['xy']['peak_sign']:+.0f}), centroid "
              f"{r['xy']['centroid_ms']:.1f} ms; neuron->kin p {r['yx']['p_shuf']:.2f}")
        assert r["xy"]["p_shuf"] <= 0.02 and r["xy"]["cv_gain_bits"] > 0
        assert abs(r["xy"]["peak_ms"] - LAG) <= 15 and r["xy"]["peak_sign"] > 0


def test_lag_bootstrap_covers_truth():
    S, B, spec = _binned("x_to_y", n_tr=200, rate_hz=(25.0,), n_per_rate=1, coupling=1.0, kernel_shape=8.0)
    r = kernels.kernel_neuron(B["Y"][0], B["X"], B["in_window"], B["cond"], B["t"], spec)
    bt = kernels.kernel_lag_bootstrap(B["Y"][0], B["X"], B["in_window"], B["cond"], B["t"], spec,
                                      r["H"], "xy", n_boot=40, seed=0)
    print(f"  peak {r['xy']['peak_ms']:.1f} ms, bootstrap 95% CI {bt['peak_ci']}; centroid CI "
          f"({bt['centroid_ci'][0]:.1f}, {bt['centroid_ci'][1]:.1f})")
    assert bt["peak_ci"][0] - 10 <= LAG <= bt["peak_ci"][1] + 10


def test_kinematic_history_choice():
    """kin_hist_mode: CV picks the H_kin with the smallest held-out SSE (all candidates on the
    same folds and rows); "max" uses the longest; KernelNull reproduces the stored null for
    both directions with the chosen histories."""
    T = sim.synthetic_template(n_trials=60, seed=3, kin="ar1", pre_ms=700)
    S = sim.simulate_session(T, "y_to_x", seed=4, rate_hz=(20.0,), n_per_rate=2, driver_rate_hz=30.0,
                             coupling=1.0, hist_pad_ms=800)
    perms = nulls.make_perms(S["trials"]["cond"].to_numpy(), 9, seed=1)
    for mode in ("cv", "max"):
        spec = kernels.make_spec(h_list=(150, 200), n_folds=3, h_kin_list=(150, 300, 450), kin_hist_mode=mode)
        B = pp.bin_session(S, 2.5, hist_pad_ms=kernels.pad_needed_ms(spec))
        r = kernels.kernel_neuron(B["Y"][0], B["X"], B["in_window"], B["cond"], B["t"], spec, perms)
        D = kernels.build_design(B["Y"][0], B["X"], B["in_window"], B["cond"], B["t"], spec)
        assert D["rows"].sum() == B["in_window"].sum(), "long kinematics: rows should start at the window start"
        if mode == "cv":
            sse = r["cv_sse_Hkin"]
            assert set(sse) == {150, 300, 450} and r["H_kin"] == min(sse, key=sse.get)
        else:
            assert r["H_kin"] == 450
        for d, H in (("xy", r["H"]), ("yx", r["H_kin"])):
            kn = kernels.KernelNull(B["Y"][0], B["X"], B["in_window"], B["cond"], B["t"], spec, H, d)
            assert np.array_equal(kn.stats(perms), r[d]["null"]), (mode, d)
        print(f"  {mode}: H_kin = {r['H_kin']} ms, cv SSE {({k: round(v, 1) for k, v in r['cv_sse_Hkin'].items()})}, "
              f"rows/trial {int(D['rows'].sum())}, adequacy p (kin) {r['yx']['adequacy_p']:.2g}")


def test_y_to_x_detected():
    S, B, spec = _binned("y_to_x", n_tr=200, rate_hz=(20.0,), n_per_rate=3, driver_rate_hz=30.0,
                         coupling=1.0, n_drivers=1)
    perms = nulls.make_perms(B["cond"], 49, seed=1)
    ps = [kernels.kernel_neuron(B["Y"][i], B["X"], B["in_window"], B["cond"], B["t"], spec, perms)["yx"]["p_shuf"]
          for i in range(3)]
    print(f"  neuron->kin p_shuf: leader {ps[0]:.3f}, null neurons {ps[1]:.2f}, {ps[2]:.2f}")
    assert ps[0] <= 0.02


if __name__ == "__main__":
    run(globals())
