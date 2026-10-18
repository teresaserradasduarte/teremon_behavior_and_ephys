"""Real-data loader smoke test (skipped if the MI cache is not reachable).
Checks layout, label consistency, the max-reach alignment (forward x peaks at
ALIGN_FRAME) and prints the per-stratum trial counts."""
from _harness import run
import numpy as np
from gctools import config, io

SESSION = config.FIRST_PASS[0]


def _available():
    return config.mi_config.cache_path().is_file()


def test_load_first_pass_session():
    if not _available():
        print("  SKIP: MI cache not reachable"); return
    S = io.load_session(*SESSION)
    n_tr = len(S["trials"])
    assert S["kin"]["X"].shape == (n_tr, config.N_FRAMES)
    assert np.isfinite(S["kin"]["X"]).all()
    assert all(len(n["spikes"]) == n_tr for n in S["neurons"])
    assert set(S["trials"]["cat"]) <= set(config.REACH_CAT)
    assert abs(S["kin"]["t"][config.ALIGN_FRAME]) < 1e-12
    regs = {r: sum(n["region"] == r for n in S["neurons"]) for r in config.REGIONS}
    print(f"  {SESSION}: {n_tr} trials, neurons {regs}")
    print("  trials per stratum:", S["trials"]["cond"].value_counts().sort_index().to_dict())


def test_long_kinematics_match_short():
    """kin_source="long" (REACH_LONG): same trials, identical kinematics inside the
    201-frame window, real (finite) history before it."""
    if not (_available() and io.reach_long_path().exists()):
        print("  SKIP: MI cache or REACH_LONG not reachable"); return
    a = io.load_session(*SESSION, verbose=False)
    b = io.load_session(*SESSION, verbose=False, kin_source="long")
    assert np.array_equal(a["trials"]["trial_idx"], b["trials"]["trial_idx"])
    i0 = int(np.flatnonzero(np.isclose(b["kin"]["t"], a["kin"]["t"][0]))[0])
    assert np.array_equal(b["kin"]["X"][:, i0:i0 + config.N_FRAMES], a["kin"]["X"])
    pre = b["kin"]["X"][:, :i0]
    print(f"  long: {b['kin']['X'].shape[1]} frames ({b['kin']['t'][0]:.2f}..{b['kin']['t'][-1]:.2f} s); "
          f"finite before -375 ms: {np.isfinite(pre).mean():.3f}")
    assert np.isfinite(pre[:, -400:]).all()                       # >= 1 s of real history


def test_alignment_is_max_reach():
    if not _available():
        print("  SKIP"); return
    NEURONS, REACH_PX = io.load_mi()
    rp = REACH_PX[SESSION]
    fwd = np.nanargmax(np.nan_to_num(rp[:, 0, :], nan=-np.inf), axis=1)
    med = np.median(fwd)
    print(f"  median argmax(forward x) = frame idx {med:.0f} (ALIGN_FRAME = {config.ALIGN_FRAME})")
    assert abs(med - config.ALIGN_FRAME) <= 1


def test_neuron_table_and_sessions():
    if not _available():
        print("  SKIP"); return
    S = io.load_session(*SESSION, verbose=False)
    nt = io.neuron_table(S)
    assert len(nt) == len(S["neurons"]) and (nt["fr_win_hz"] >= 0).all()
    for r, g in nt.groupby("region"):
        print(f"  {r}: window rate pct 5/50/95 = "
              f"{np.percentile(g['fr_win_hz'], [5, 50, 95]).round(1).tolist()} Hz")
    ls = io.list_sessions()
    assert len(ls) == 29
    print(ls.sort_values("n_kept", ascending=False).head(5).to_string(index=False))


if __name__ == "__main__":
    run(globals())
