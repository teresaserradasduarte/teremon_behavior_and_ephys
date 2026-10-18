"""nulls.py - trial-shuffle null (CLAUDE.md rule 4).

Neuron trial i is paired with the partner signal of trial perm[i], drawn from the
SAME condition stratum (SHUFFLE_WITHIN, e.g. spout side) and never i itself. This
keeps each signal's autocorrelation, its time-locking to the reach, and the
condition-level association; what the null destroys is trial-by-trial coupling
within a condition -- which is exactly what survives a significant result.

The permutations are drawn ONCE per session (make_perms) and reused for every
neuron and both directions, so all tests of a session share one null.
"""
import numpy as np


def perm_within(cond, rng):
    """Derangement within each stratum -> (n_trials,) partner index.
    Strata with a single trial cannot be shuffled and are left in place (the
    caller should drop them; io.list_sessions reports stratum sizes)."""
    cond = np.asarray(cond)
    perm = np.arange(cond.size)
    for c in np.unique(cond):
        idx = np.flatnonzero(cond == c)
        if idx.size < 2:
            continue
        while True:
            p = rng.permutation(idx.size)
            if not np.any(p == np.arange(idx.size)):
                break
        perm[idx] = idx[p]
    return perm


def make_perms(cond, n_shuffle, seed=0):
    """(n_shuffle, n_trials) int array of within-stratum derangements. Drawn sequentially
    from one seeded generator: make_perms(cond, n, s)[:k] == make_perms(cond, k, s), so a
    longer sequence EXTENDS a shorter one (used by the adaptive shuffles)."""
    rng = np.random.default_rng(seed)
    return np.stack([perm_within(cond, rng) for _ in range(n_shuffle)]) if n_shuffle \
        else np.zeros((0, len(cond)), dtype=int)


def perm_batch(cond, start, stop, seed=0):
    """Shuffles start..stop-1 of the session's sequence (identical to
    make_perms(cond, stop, seed)[start:stop]). The runs generate the whole sequence once
    per session; this is the single-batch convenience form."""
    return make_perms(cond, stop, seed)[start:stop]
