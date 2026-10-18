"""Parameter-hashed cache: stable keys, roundtrip, load vs rerun. Uses a temp dir,
never the Dropbox tree."""
from _harness import run
import tempfile
import numpy as np
from gctools import cache


def test_hash_stable_and_sensitive():
    p = dict(bin_ms=2.5, window=(-375, 125), arr=np.arange(5.0), name="x")
    assert cache.param_hash(p) == cache.param_hash(dict(reversed(list(p.items()))))
    assert cache.param_hash(p) == cache.param_hash(dict(p, window=[-375, 125]))
    for k, v in [("bin_ms", 10.0), ("window", (-370, 125)), ("arr", np.arange(5.0) + 1e-9),
                 ("name", "y")]:
        assert cache.param_hash(dict(p, **{k: v})) != cache.param_hash(p), k
    big = np.zeros(10_000); big2 = big.copy(); big2[5_000] = 1   # repr would truncate these
    assert cache.param_hash(dict(a=big)) != cache.param_hash(dict(a=big2))


def test_roundtrip_load_and_rerun():
    calls = []
    def compute():
        calls.append(1); return dict(v=np.arange(3), n=len(calls))
    with tempfile.TemporaryDirectory() as d:
        a = cache.cached("t", dict(x=1), compute, cache_dir=d)
        b = cache.cached("t", dict(x=1), compute, cache_dir=d)
        assert len(calls) == 1 and b["n"] == 1 and np.array_equal(a["v"], b["v"])
        c = cache.cached("t", dict(x=1), compute, mode="rerun", cache_dir=d)
        assert len(calls) == 2 and c["n"] == 2
        cache.cached("t", dict(x=2), compute, cache_dir=d)
        assert len(calls) == 3


if __name__ == "__main__":
    run(globals())
