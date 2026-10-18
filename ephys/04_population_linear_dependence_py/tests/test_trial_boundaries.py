"""No lag crosses a trial boundary (CLAUDE.md rule 2), and full vs reduced models use
identical rows (rule 3).

Every sample is encoded as z = 1000*trial + bin, so any lagged cell can be decoded
back to (source trial, source bin) and checked against its row's own (trial, bin).
"""
from _harness import run
import numpy as np
from gctools import design as ds

N_TR, N_B = 7, 60
TRIAL = np.repeat(np.arange(N_TR)[:, None], N_B, axis=1)
TBIN  = np.repeat(np.arange(N_B)[None, :], N_TR, axis=0)
Z     = 1000.0 * TRIAL + TBIN                                  # (n_trials, n_bins)
IN_WIN = np.arange(N_B) >= 10


def _design(Z=Z):
    lags_own = np.arange(1, 6)
    lags_x   = np.arange(1, 13)
    blocks = dict(own=ds.lagged_block(Z, ds.lag_matrix, lags_own),
                  cross=ds.lagged_block(Z, ds.lag_matrix, lags_x))
    return ds.build_design(blocks, Z, IN_WIN), lags_own, lags_x


def test_lagged_cells_decode_to_own_trial_and_past():
    D, lags_own, lags_x = _design()
    Xc = D["X"][:, D["groups"]["cross"]]
    src_trial = np.floor(Xc / 1000).astype(int)
    src_bin   = (Xc - 1000 * src_trial).round().astype(int)
    assert np.all(src_trial == D["trial"][:, None]), "a lag reached into another trial"
    assert np.all(src_bin == D["tbin"][:, None] - lags_x[None, :]), "wrong lag"
    assert np.all(src_bin < D["tbin"][:, None]), "a lag is not in the past"
    print(f"  {D['X'].shape[0]} rows x {D['X'].shape[1]} cols: every cell from its own "
          f"trial, strictly earlier bin")


def test_rows_needing_pre_trial_history_are_dropped():
    D, _, lags_x = _design()
    first_ok = max(lags_x.max(), 10)             # longest lag vs window start
    assert D["tbin"].min() == first_ok
    rm = D["row_mask"]
    assert not rm[:, :first_ok].any() and rm[:, first_ok:].all()
    assert D["X"].shape[0] == N_TR * (N_B - first_ok)
    # NaN inside a trial (e.g. untracked kinematics) also removes the affected rows
    Z2 = Z.copy(); Z2[3, 30] = np.nan
    D2, _, _ = _design(Z2)
    lost = set(zip(*np.nonzero(D["row_mask"] & ~D2["row_mask"])))
    expect = {(3, 30)} | {(3, 30 + L) for L in range(1, 13) if 30 + L < N_B}
    assert lost == expect, sorted(lost)


def test_trial_order_irrelevant():
    D, _, _ = _design()
    perm = np.random.default_rng(0).permutation(N_TR)
    Dp, _, _ = _design(Z[perm])
    for new_i, old_i in enumerate(perm):
        a = D["X"][D["trial"] == old_i]
        b = Dp["X"][Dp["trial"] == new_i]
        assert np.array_equal(a, b)


def test_full_and_reduced_share_rows():
    D, _, _ = _design()
    R = ds.reduced(D, "cross")
    assert np.array_equal(R["y"], D["y"]) and np.array_equal(R["trial"], D["trial"])
    assert np.array_equal(R["tbin"], D["tbin"])
    assert np.array_equal(R["X"], D["X"][:, D["groups"]["own"]])
    assert list(R["groups"]) == ["own"]
    print(f"  full {D['X'].shape}  reduced {R['X'].shape}  -- identical rows and y")


if __name__ == "__main__":
    run(globals())
