"""cache.py - parameter-hashed disk cache + Dropbox-safe saving.

Every parameter that changes the numbers goes into `params`; the file is
CACHE_DIR/<name>_<hash>.pkl with a <name>_<hash>.json sidecar listing the params,
so a cache can be identified by eye. Style/figure parameters stay OUT of params.

The output tree is Dropbox-synced: Dropbox briefly locks a file while uploading
and an overwrite then raises OSError. save_pickle / savefig_retry retry instead of
crashing (same idea as _savefig_retry in MI 12_MI_2D_kin_x_neuron.py).
"""
import hashlib, json, os, pickle, time
from pathlib import Path
import numpy as np


def _canon(obj):
    """JSON-able canonical form. Arrays are hashed by content (repr would truncate)."""
    if isinstance(obj, dict):
        return {str(k): _canon(v) for k, v in sorted(obj.items(), key=lambda kv: str(kv[0]))}
    if isinstance(obj, (list, tuple)):
        return [_canon(v) for v in obj]
    if isinstance(obj, np.ndarray):
        a = np.ascontiguousarray(obj)
        return {"__ndarray__": hashlib.sha1(a.tobytes()).hexdigest(),
                "shape": list(a.shape), "dtype": str(a.dtype)}
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, Path):
        return str(obj)
    return obj


def param_hash(params, n=10):
    """First n hex chars of SHA-1 over the canonical JSON of params."""
    s = json.dumps(_canon(params), sort_keys=True, default=repr)
    return hashlib.sha1(s.encode()).hexdigest()[:n]


def _long(path):
    """Windows refuses file paths over 260 characters unless given the \\\\?\\ prefix;
    the Dropbox output tree is deep, so long paths get it."""
    from .config import long_path, warn_long
    warn_long(path)
    return long_path(path)


def save_pickle(obj, path, tries=6):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for t in range(tries):
        try:
            with open(_long(path), "wb") as f:
                pickle.dump(obj, f, protocol=4)
            return path
        except OSError:
            if t == tries - 1: raise
            time.sleep(0.5 * (t + 1))


def load_pickle(path):
    with open(_long(path), "rb") as f:
        return pickle.load(f)


def savefig_retry(fig, path, tries=6, **kw):
    kw = {"dpi": 150, "bbox_inches": "tight", **kw}
    for t in range(tries):
        try:
            return fig.savefig(_long(path), **kw)
        except OSError:
            if t == tries - 1: raise
            time.sleep(0.5 * (t + 1))


def cached(name, params, compute_fn, mode="load", save=True, cache_dir=None):
    """Return compute_fn() through the cache.
    mode "load"  -> reuse CACHE_DIR/<name>_<hash>.pkl if it exists, else compute
         "rerun" -> always compute (and overwrite if save)."""
    if cache_dir is None:
        from . import config
        cache_dir = config.CACHE_DIR
    cache_dir = Path(cache_dir)
    h    = param_hash(params)
    path = cache_dir / f"{name}_{h}.pkl"
    if mode == "load" and path.is_file():
        print(f"[cache] loaded   {path.name}")
        return load_pickle(path)
    if mode not in ("load", "rerun"):
        raise ValueError(f"mode must be 'load' or 'rerun', got {mode!r}")
    t0  = time.time()
    out = compute_fn()
    print(f"[cache] computed {name}_{h} in {time.time()-t0:.1f}s")
    if save:
        save_pickle(out, path)
        with open(path.with_suffix(".json"), "w") as f:
            json.dump(_canon(params), f, indent=1, default=repr)
    return out
