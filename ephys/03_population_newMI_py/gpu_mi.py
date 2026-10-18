"""gpu_mi.py - batched GPU (torch) mutual-information estimators + shuffle-null
Excess_MI. Falls back to CPU automatically if no CUDA device is available.

KSG (continuous-continuous) replicates sklearn.feature_selection.mutual_info_regression
(Kraskov, Chebyshev joint distance). Ross (continuous-discrete) replicates
mutual_info_classif (same-class kNN radius). Both operate on a BATCH of independent
MI problems at once: shape (B, N) where B = 1 real + N_SHUF shuffle nulls, optionally
times however many (start,stop) grid windows are chunked together, and N = n_trials.

Binning (plug-in) stays on CPU/numpy - it's the cheapest estimator already (O(N)
counting vs KSG/Ross's O(N^2) kNN search); porting it to GPU would add kernel-launch
overhead exceeding its CPU cost at N in the low hundreds.
"""
import os
# torch bundles its own OpenMP runtime (libiomp5md.dll), which collides with the
# MKL-linked OpenMP numpy/scipy already loaded in the kernel process (OMP Error #15).
# Must be set before torch touches any threaded op - setting it here, before import,
# covers both fresh processes and long-lived Jupyter/VS Code interactive kernels.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import torch

NAT2BIT = 1.0 / np.log(2)
DEVICE  = "cuda" if torch.cuda.is_available() else "cpu"


def zscore(x, axis=-1, eps=1e-12):
    """Z-score along axis. Guards zero-std with eps (returns zeros, not NaN/inf)."""
    x = np.asarray(x, dtype=np.float64)
    mu = x.mean(axis=axis, keepdims=True)
    sd = x.std(axis=axis, keepdims=True)
    return (x - mu) / np.maximum(sd, eps)


# ── KSG (continuous-continuous) ───────────────────────────────────────────────
def ksg_mi_batch(X, Y, k=3, device=DEVICE):
    """
    X, Y : (B, N) array-like, already z-scored. Independent MI problems along B.
    Returns (B,) numpy float32, MI in bits.

    Replicates the classic Kraskov estimator (KSG-1) as used by sklearn's
    mutual_info_regression: Chebyshev joint distance, strict (open-ball) marginal
    neighbor counts.
    """
    Xt = torch.as_tensor(np.asarray(X), dtype=torch.float32, device=device)
    Yt = torch.as_tensor(np.asarray(Y), dtype=torch.float32, device=device)
    if Yt.shape[0] == 1 and Xt.shape[0] > 1:
        Yt = Yt.expand(Xt.shape[0], -1)
    B, N = Xt.shape

    dX = (Xt.unsqueeze(2) - Xt.unsqueeze(1)).abs()   # (B,N,N)
    dY = (Yt.unsqueeze(2) - Yt.unsqueeze(1)).abs()   # (B,N,N)
    dZ = torch.maximum(dX, dY)

    eye = torch.eye(N, dtype=torch.bool, device=device).unsqueeze(0)
    dZ = dZ.masked_fill(eye, float("inf"))

    eps = torch.kthvalue(dZ, k, dim=2).values         # (B,N) - k-th NN dist (self excluded)

    dX_ = dX.masked_fill(eye, float("inf"))
    dY_ = dY.masked_fill(eye, float("inf"))
    n_x = (dX_ < eps.unsqueeze(2)).sum(dim=2).to(torch.float32)  # strict, self excluded
    n_y = (dY_ < eps.unsqueeze(2)).sum(dim=2).to(torch.float32)

    psi = torch.special.digamma
    mi_nats = psi(torch.tensor(float(k), device=device)) \
              - (psi(n_x + 1) + psi(n_y + 1)).mean(dim=1) \
              + psi(torch.tensor(float(N), device=device))
    mi_bits = (mi_nats * NAT2BIT).clamp(min=0)
    return mi_bits.cpu().numpy().astype(np.float32)


# ── Ross (continuous-discrete) ────────────────────────────────────────────────
def ross_mi_batch(X, y, k=3, device=DEVICE):
    """
    X : (B, N) z-scored float, array-like.
    y : (B, N) or (1, N) int class labels (raw counts - NOT z-scored).
    Returns (B,) numpy float32, MI in bits.

    Replicates sklearn's mutual_info_classif / _compute_mi_cd (Ross 2014) exactly:
      mi = psi(n_valid) + mean(psi(k_eff)) - mean(psi(class_size)) - mean(psi(m_all))
    where, per point i (only over points whose class_size > 1 - singleton classes are
    dropped entirely, including from n_valid):
      class_size(i) = # points sharing i's label (self included)
      k_eff(i)       = min(k, class_size(i) - 1)
      radius(i)      = distance to the k_eff(i)-th nearest SAME-CLASS neighbor (self excl.)
      m_all(i)       = # points (ANY class) with distance < radius(i), SELF INCLUDED
                       (sklearn's query_radius counts the self point at distance 0)
    """
    Xt = torch.as_tensor(np.asarray(X), dtype=torch.float32, device=device)
    yt = torch.as_tensor(np.asarray(y), dtype=torch.int64, device=device)
    if yt.shape[0] == 1 and Xt.shape[0] > 1:
        yt = yt.expand(Xt.shape[0], -1)
    B, N = Xt.shape

    dX = (Xt.unsqueeze(2) - Xt.unsqueeze(1)).abs()                     # (B,N,N)
    same_class = (yt.unsqueeze(2) == yt.unsqueeze(1))                   # (B,N,N)
    eye = torch.eye(N, dtype=torch.bool, device=device).unsqueeze(0)

    class_size = same_class.sum(dim=2).to(torch.float32)                # (B,N), self included
    valid = class_size > 1
    valid_f = valid.to(torch.float32)
    n_valid = valid_f.sum(dim=1).clamp(min=1.0)                         # (B,)

    m_excl_self = (class_size - 1.0).clamp(min=1.0)                     # class size - self
    k_eff = torch.clamp(torch.full_like(m_excl_self, float(k)), max=m_excl_self).clamp(min=1.0)

    dX_same = dX.masked_fill(~same_class | eye, float("inf"))
    k_cap = int(min(k, N - 1))
    topk_vals, _ = torch.topk(dX_same, k_cap, dim=2, largest=False)      # (B,N,k_cap)
    k_idx = (k_eff.long() - 1).clamp(min=0, max=k_cap - 1)               # (B,N)
    radius = torch.gather(topk_vals, 2, k_idx.unsqueeze(2)).squeeze(2)   # (B,N)
    radius = torch.where(torch.isfinite(radius), radius, torch.zeros_like(radius))

    dX_all = dX.masked_fill(eye, float("inf"))
    N_d = (dX_all < radius.unsqueeze(2)).sum(dim=2).to(torch.float32)   # (B,N), self excluded, strict
    m_all = N_d + 1.0                                                    # sklearn counts self too

    def masked_mean(v):
        return (v * valid_f).sum(dim=1) / n_valid

    psi = torch.special.digamma
    mi_nats = psi(n_valid) \
              + masked_mean(psi(k_eff)) \
              - masked_mean(psi(class_size)) \
              - masked_mean(psi(m_all))
    mi_bits = (mi_nats * NAT2BIT).clamp(min=0)
    return mi_bits.cpu().numpy().astype(np.float32)


# ── Binning (plug-in), CPU/numpy, vectorized over the batch (shuffle) dim ─────
def bin_mi_batch(kin_batch, count_batch, k=None, n_bins=16):
    """
    kin_batch, count_batch : (B, N) array-like (kin usually shuffled per row,
    count usually identical across rows - broadcasting handled by caller).
    Returns (B,) numpy float64, MI in bits. Plain per-row loop (cheap; Binning is
    O(N) per call, not worth vectorizing further or moving to GPU).
    """
    kin_batch = np.asarray(kin_batch)
    count_batch = np.asarray(count_batch)
    if count_batch.shape[0] == 1 and kin_batch.shape[0] > 1:
        count_batch = np.broadcast_to(count_batch, kin_batch.shape)

    B = kin_batch.shape[0]
    out = np.zeros(B, dtype=np.float64)
    for b in range(B):
        kin_m, c_w = kin_batch[b], count_batch[b]
        edges = np.unique(np.quantile(kin_m, np.linspace(0, 1, n_bins + 1)))
        kb = np.clip(np.digitize(kin_m, edges[1:-1]), 0, len(edges) - 2)
        cv = np.unique(c_w)
        cm_ = {v: idx for idx, v in enumerate(cv)}
        ci = np.array([cm_[v] for v in c_w])
        J = np.zeros((len(cv), len(edges) - 1))
        np.add.at(J, (ci, kb), 1.0)
        P = J / J.sum()
        pc = P.sum(1, keepdims=True); pk_ = P.sum(0, keepdims=True); nz = P > 0
        out[b] = float(np.sum(P[nz] * np.log2(P[nz] / (pc @ pk_ + 1e-300)[nz])))
    return out


# ── generic shuffle-null Excess_MI wrapper ────────────────────────────────────
def excess_mi(estimator_fn, kin_z, count_for_est, n_shuf, rng, k=3, **kwargs):
    """
    Builds the (1+n_shuf, N) batch (real row + n_shuf trial-permutations of kin_z),
    calls estimator_fn ONCE, and splits the result into observed vs null.

    kin_z         : (N,) already z-scored (or raw int kin for estimators that don't
                    z-score, e.g. Ross/Binning callers still pass whatever kin array
                    is appropriate for that estimator - only kin gets shuffled).
    count_for_est : (N,) count (z-scored for KSG, raw int for Ross/Binning).
    estimator_fn  : one of ksg_mi_batch, ross_mi_batch, bin_mi_batch - signature
                    (X_batch, Y_batch, k=..., **kwargs) -> (B,) bits.

    Returns (mi_obs, mi_excess, mi_null_mean, mi_null_std) as plain floats.
    """
    kin_z = np.asarray(kin_z)
    N = kin_z.shape[0]
    X_batch = np.empty((1 + n_shuf, N), dtype=kin_z.dtype)
    X_batch[0] = kin_z
    for j in range(n_shuf):
        X_batch[1 + j] = rng.permutation(kin_z)

    Y_batch = np.asarray(count_for_est).reshape(1, N)   # broadcast inside estimator_fn

    mi_all = estimator_fn(X_batch, Y_batch, k=k, **kwargs)
    mi_obs = float(mi_all[0])
    mi_null = mi_all[1:]
    mi_null_mean = float(mi_null.mean())
    mi_null_std = float(mi_null.std())
    mi_excess = mi_obs - mi_null_mean
    return mi_obs, mi_excess, mi_null_mean, mi_null_std


def excess_mi_windows(estimator_fn, kin_stack, count_windows, k=3, chunk_size=40, **kwargs):
    """
    Batches Excess_MI over MANY (start,stop) grid windows at once - the main GPU
    utilization lever (one estimator call handles chunk_size windows x (1+n_shuf) rows
    instead of one call per window).

    kin_stack     : (1+n_shuf, N) - real kin (row 0) + n_shuf trial-permutations of it.
                    IDENTICAL across all windows for a given neuron (kin doesn't depend
                    on the window; only count does), so build this ONCE per neuron.
    count_windows : (n_windows, N) - one row per (start,stop) window (already
                    preprocessed: z-scored per-row for KSG, raw int for Ross/Binning).
    Returns (mi_obs, mi_excess, mi_null_mean, mi_null_std), each (n_windows,) numpy.
    """
    kin_stack = np.asarray(kin_stack)
    count_windows = np.asarray(count_windows)
    n_shuf_plus1, N = kin_stack.shape
    n_windows = count_windows.shape[0]

    mi_obs = np.full(n_windows, np.nan)
    mi_excess = np.full(n_windows, np.nan)
    mi_null_mean = np.full(n_windows, np.nan)
    mi_null_std = np.full(n_windows, np.nan)

    for start in range(0, n_windows, chunk_size):
        stop = min(start + chunk_size, n_windows)
        w_block = stop - start
        X_chunk = np.tile(kin_stack, (w_block, 1))                          # (w_block*B, N)
        Y_chunk = np.repeat(count_windows[start:stop], n_shuf_plus1, axis=0)  # (w_block*B, N)

        mi_all = estimator_fn(X_chunk, Y_chunk, k=k, **kwargs).reshape(w_block, n_shuf_plus1)
        mi_obs[start:stop] = mi_all[:, 0]
        null = mi_all[:, 1:]
        mi_null_mean[start:stop] = null.mean(axis=1)
        mi_null_std[start:stop] = null.std(axis=1)
        mi_excess[start:stop] = mi_obs[start:stop] - mi_null_mean[start:stop]

    return mi_obs, mi_excess, mi_null_mean, mi_null_std


def build_kin_stack(kin_z, n_shuf, rng):
    """(1+n_shuf, N): real z-scored kin (row 0) + n_shuf trial-permutations of it."""
    kin_z = np.asarray(kin_z)
    N = kin_z.shape[0]
    stack = np.empty((1 + n_shuf, N), dtype=kin_z.dtype)
    stack[0] = kin_z
    for j in range(n_shuf):
        stack[1 + j] = rng.permutation(kin_z)
    return stack
