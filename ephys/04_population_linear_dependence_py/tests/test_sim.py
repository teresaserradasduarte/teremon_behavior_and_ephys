"""Simulator sanity (NOT the GC test itself -- that comes with the GLMs in M2+).
Synthetic WHITE-noise kinematics are used for the kernel check: with white input the
spike-triggered average of past kinematics is proportional to the coupling kernel,
so its peak must sit at the true lag. On real (smooth, task-locked) kinematics the
STA is smeared -- that is shown, not tested, in scripts/00_core_checks.py.
"""
from _harness import run
import numpy as np
import pandas as pd
from gctools import sim, preprocess as pp

LAG_MS, BIN_MS = 50.0, 2.5
LAGS = np.arange(1, 61)                                    # 2.5 .. 150 ms


def _pooled_sta(S, pair=None):
    B = pp.bin_session(S, BIN_MS, hist_pad_ms=160)
    return np.mean([sim.sta_kin(B["Y"][i], B["X"], B["in_window"], LAGS, pair=pair)
                    for i in range(B["Y"].shape[0])], axis=0), B


def test_rates_hit_targets():
    T = sim.synthetic_template(n_trials=300, seed=2)
    S = sim.simulate_session(T, "null_independent", rate_hz=(2.0, 20.0, 60.0), n_per_rate=3,
                             seed=3)
    w0, w1 = -0.375, 0.125
    for n in S["neurons"]:
        emp = np.mean([np.sum((s >= w0) & (s < w1)) for s in n["spikes"]]) / (w1 - w0)
        rel = emp / n["rate_target_hz"] - 1
        print(f"  target {n['rate_target_hz']:5.1f} Hz  empirical {emp:6.2f} Hz  ({rel:+.1%})")
        assert abs(rel) < 0.15


def test_deterministic_seed():
    T = sim.synthetic_template(n_trials=20, seed=2)
    a = sim.simulate_session(T, "x_to_y", rate_hz=(20.0,), n_per_rate=2, seed=7)
    b = sim.simulate_session(T, "x_to_y", rate_hz=(20.0,), n_per_rate=2, seed=7)
    for na, nb in zip(a["neurons"], b["neurons"]):
        assert all(np.array_equal(x, y) for x, y in zip(na["spikes"], nb["spikes"]))


def test_x_to_y_kernel_recovered_null_flat():
    T = sim.synthetic_template(n_trials=400, seed=4)
    kw = dict(rate_hz=(30.0,), n_per_rate=4, coupling=1.5, lag_ms=LAG_MS, kernel_shape=20.0)
    S1 = sim.simulate_session(T, "x_to_y", seed=5, **kw)
    S0 = sim.simulate_session(T, "null_independent", seed=5, **kw)
    sta1, B = _pooled_sta(S1)
    sta0, _ = _pooled_sta(S0)
    rng = np.random.default_rng(0)
    shuf = np.array([_pooled_sta(S0, pair=sim.perm_within(B["cond"], rng))[0]
                     for _ in range(30)])
    band = np.percentile(np.abs(shuf), 99, axis=0)
    pk = LAGS[np.argmax(sta1)] * BIN_MS
    frac_out0 = np.mean(np.abs(sta0) > band)
    print(f"  x_to_y: STA peak at {pk:.1f} ms (true mode {LAG_MS} ms), "
          f"peak/band = {sta1.max() / band[np.argmax(sta1)]:.1f}")
    print(f"  null  : {frac_out0:.0%} of lags outside the 99% within-LCR shuffle band, "
          f"max|STA|/band = {np.max(np.abs(sta0) / band):.2f}")
    assert abs(pk - LAG_MS) <= 2 * BIN_MS
    assert sta1.max() > 3 * band[np.argmax(sta1)]
    assert frac_out0 <= 0.10


def test_y_to_x_kinematics_follow_driver():
    """y_to_x: kinematics AFTER a driver spike deviate, peaking ~lag_ms later; a
    non-driver neuron shows nothing beyond its shuffle band."""
    from gctools import design as ds
    T = sim.synthetic_template(n_trials=300, seed=6)
    S = sim.simulate_session(T, "y_to_x", rate_hz=(20.0,), n_per_rate=3, driver_rate_hz=30.0,
                             coupling=1.0, lag_ms=LAG_MS, seed=7)
    assert S["truth"]["expect"]["y_to_x"].tolist() == [True, False, False]
    B = pp.bin_session(S, BIN_MS, hist_pad_ms=160)

    def kin_after_spike(Y, X):                    # sum y[t-j] x[t] / sum y[t-j] - mean x
        L = ds.lagged_block(Y.astype(float), ds.lag_matrix, LAGS)
        ok = B["in_window"][None, :, None] & np.isfinite(L) & np.isfinite(X)[:, :, None]
        Xb = np.where(ok, X[:, :, None], 0.0); Lz = np.where(ok, L, 0.0)
        return ((Lz * Xb).sum((0, 1)) / np.maximum(Lz.sum((0, 1)), 1e-12)
                - Xb.sum((0, 1)) / np.maximum(ok.sum((0, 1)), 1))

    drv = kin_after_spike(B["Y"][0], B["X"])
    oth = kin_after_spike(B["Y"][1], B["X"])
    pk = LAGS[np.argmax(drv)] * BIN_MS
    print(f"  driver: kinematics peak {pk:.1f} ms after its spikes (true mode {LAG_MS}); "
          f"peak {drv.max():.3f} vs non-driver max |.| {np.abs(oth).max():.3f}")
    assert abs(pk - LAG_MS) <= 2 * BIN_MS
    assert drv.max() > 5 * np.abs(oth).max()


def test_y_to_x_several_leaders():
    """n_drivers=3 in ONE session: the kinematics follow each leader ~lag_ms later, each
    by about as much as a single leader would; non-leaders show nothing. n_drivers=1 is
    bit-identical to the default."""
    from gctools import design as ds
    T = sim.synthetic_template(n_trials=300, seed=6)
    kw = dict(rate_hz=(20.0,), n_per_rate=6, driver_rate_hz=30.0, coupling=1.0, lag_ms=LAG_MS, seed=7)
    S1 = sim.simulate_session(T, "y_to_x", **kw)
    S1b = sim.simulate_session(T, "y_to_x", n_drivers=1, **kw)
    assert np.array_equal(S1["kin"]["X"], S1b["kin"]["X"])
    S3 = sim.simulate_session(T, "y_to_x", n_drivers=3, **kw)
    assert S3["truth"]["drivers"] == [0, 1, 2]
    assert S3["truth"]["expect"]["y_to_x"].tolist() == [True] * 3 + [False] * 3

    def kin_after_spike(S, i):
        B = pp.bin_session(S, BIN_MS, hist_pad_ms=160)
        L = ds.lagged_block(B["Y"][i].astype(float), ds.lag_matrix, LAGS)
        X = B["X"]
        ok = B["in_window"][None, :, None] & np.isfinite(L) & np.isfinite(X)[:, :, None]
        Xb = np.where(ok, X[:, :, None], 0.0); Lz = np.where(ok, L, 0.0)
        return ((Lz * Xb).sum((0, 1)) / np.maximum(Lz.sum((0, 1)), 1e-12)
                - Xb.sum((0, 1)) / np.maximum(ok.sum((0, 1)), 1))

    single = kin_after_spike(S1, 0).max()
    for i in range(6):
        k = kin_after_spike(S3, i)
        pk = LAGS[np.argmax(k)] * BIN_MS
        print(f"  neuron {i} ({'leader' if i < 3 else 'null  '}): peak {k.max():.3f} at {pk:.1f} ms "
              f"(single leader alone: {single:.3f})")
        if i < 3:
            assert abs(pk - LAG_MS) <= 2 * BIN_MS and 0.6 < k.max() / single < 1.4
        else:
            assert np.abs(k).max() < 0.2 * single
    print(f"  added signal: {S3['truth']['added_sd_frac']:.2f} x sd (1 leader: "
          f"{S1['truth']['added_sd_frac']:.2f})")


def test_null_shared_drive_task_locked_not_trial_specific():
    """null_shared_drive: rate follows the spout-side MEAN kinematics -> PSTHs differ by side
    (task-locked), but within a side a trial's spike count does not track its OWN kinematics
    (x_to_y does)."""
    T = sim.synthetic_template(n_trials=240, seed=7, kin="ar1", n_cond=3)
    side = T["trials"]["LCR"].to_numpy()
    T["kin"]["X"] = T["kin"]["X"] + 2.0 * (side[:, None] - 2) * np.linspace(0, 1, T["kin"]["X"].shape[1])  # side-specific reach
    out = {}
    for scen in ("null_shared_drive", "x_to_y"):
        S = sim.simulate_session(T, scen, rate_hz=(30.0,), n_per_rate=4, coupling=1.0, seed=3)
        assert not S["truth"]["expect"]["x_to_y"].any() if scen == "null_shared_drive" else S["truth"]["expect"]["x_to_y"].all()
        B = pp.bin_session(S, 10.0, hist_pad_ms=0)
        w = B["in_window"]
        late = w & (B["t"] > 0)
        cnt = B["Y"][:, :, late].sum(axis=2)                                  # (n_neu, n_trials)
        xm = np.nanmean(B["X"][:, late], axis=1)
        side_effect = np.mean([np.ptp([cnt[i, side == s].mean() for s in (1, 2, 3)]) / cnt[i].mean() for i in range(4)])
        r_within = np.mean([np.corrcoef(cnt[i] - pd.Series(cnt[i]).groupby(side).transform("mean"),
                                        xm - pd.Series(xm).groupby(side).transform("mean"))[0, 1] for i in range(4)])
        out[scen] = (side_effect, r_within)
        print(f"  {scen:18s}: rate difference between sides {side_effect:.2f} (rel.), "
              f"within-side corr(count, own kinematics) {r_within:+.3f}")
    assert out["null_shared_drive"][0] > 0.2, "shared drive should make the rate task-locked"
    assert abs(out["null_shared_drive"][1]) < 0.1, "no trial-specific coupling expected"
    assert out["x_to_y"][1] > 0.2, "x_to_y should couple to the trial's own kinematics"


if __name__ == "__main__":
    run(globals())
