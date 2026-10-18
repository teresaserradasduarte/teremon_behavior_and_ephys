# %% ── validate gpu_mi.py against sklearn + analytic MI (run before using in 05) ──
import numpy as np
from sklearn.feature_selection import mutual_info_regression, mutual_info_classif
import gpu_mi

print(f"gpu_mi.DEVICE = {gpu_mi.DEVICE}")
rng = np.random.default_rng(0)

# %% ── KSG: Gaussian X, Y = rho*X + sqrt(1-rho^2)*noise -> known analytic MI ─────
print("\n=== KSG validation (Gaussian, known analytic MI) ===")
for N in (100, 300, 500):
    for rho in (0.3, 0.6, 0.85):
        X = rng.normal(size=N)
        noise = rng.normal(size=N)
        Y = rho * X + np.sqrt(1 - rho**2) * noise
        analytic = -0.5 * np.log2(1 - rho**2)

        Xz, Yz = gpu_mi.zscore(X), gpu_mi.zscore(Y)

        sk = mutual_info_regression(Xz.reshape(-1, 1), Yz + rng.normal(0, 1e-9, N),
                                     n_neighbors=3, random_state=0)[0] * gpu_mi.NAT2BIT
        gp = gpu_mi.ksg_mi_batch(Xz.reshape(1, -1), Yz.reshape(1, -1), k=3)[0]

        print(f"  N={N:4d}  rho={rho:.2f}  analytic={analytic:.3f}  "
              f"sklearn={sk:.3f}  gpu_mi={gp:.3f}  |gpu-sk|={abs(gp-sk):.3f}")

# %% ── Ross: continuous X, discrete Y classes with known separation ──────────────
print("\n=== Ross validation (continuous X, discrete Y) ===")
for N in (100, 300):
    for n_classes in (3, 6):
        y = rng.integers(0, n_classes, size=N)
        # X depends on y (class means spread out) + noise -> real dependency
        X = y.astype(float) * 1.5 + rng.normal(size=N)
        Xz = gpu_mi.zscore(X)

        sk = mutual_info_classif(Xz.reshape(-1, 1), y, n_neighbors=3,
                                  random_state=0)[0] * gpu_mi.NAT2BIT
        gp = gpu_mi.ross_mi_batch(Xz.reshape(1, -1), y.reshape(1, -1), k=3)[0]

        print(f"  N={N:4d}  n_classes={n_classes}  sklearn={sk:.3f}  gpu_mi={gp:.3f}  "
              f"|gpu-sk|={abs(gp-sk):.3f}")

# %% ── Ross: null case (X independent of Y) -> both should be near 0 ─────────────
print("\n=== Ross validation (null case, X independent of Y) ===")
N = 300
y = rng.integers(0, 4, size=N)
X = rng.normal(size=N)
Xz = gpu_mi.zscore(X)
sk = mutual_info_classif(Xz.reshape(-1, 1), y, n_neighbors=3, random_state=0)[0] * gpu_mi.NAT2BIT
gp = gpu_mi.ross_mi_batch(Xz.reshape(1, -1), y.reshape(1, -1), k=3)[0]
print(f"  N={N}  sklearn={sk:.3f}  gpu_mi={gp:.3f}  (expect both near 0)")

# %% ── batching sanity check: batched call == looped single calls ────────────────
print("\n=== Batch consistency check ===")
N = 150
X_batch = np.stack([rng.normal(size=N) for _ in range(5)])
Y_fixed = rng.normal(size=N).reshape(1, -1)
mi_batched = gpu_mi.ksg_mi_batch(X_batch, Y_fixed, k=3)
mi_looped = np.array([gpu_mi.ksg_mi_batch(X_batch[i:i+1], Y_fixed, k=3)[0] for i in range(5)])
print("  batched:", np.round(mi_batched, 4))
print("  looped :", np.round(mi_looped, 4))
print("  max abs diff:", np.max(np.abs(mi_batched - mi_looped)))

# %% ── excess_mi wrapper smoke test ────────────────────────────────────────────
print("\n=== excess_mi wrapper smoke test (KSG, null case) ===")
N = 200
kin = rng.normal(size=N)
count = rng.normal(size=N)  # independent -> expect mi_excess near 0
kin_z = gpu_mi.zscore(kin)
count_z = gpu_mi.zscore(count)
mi_obs, mi_excess, mi_null_mean, mi_null_std = gpu_mi.excess_mi(
    gpu_mi.ksg_mi_batch, kin_z, count_z, n_shuf=20, rng=rng, k=3)
print(f"  mi_obs={mi_obs:.4f}  mi_excess={mi_excess:.4f}  "
      f"mi_null_mean={mi_null_mean:.4f}  mi_null_std={mi_null_std:.4f}")
print("  (mi_excess should be small/near-zero-ish since kin/count are independent)")

# %% ── excess_mi_windows: chunked multi-window batching vs per-window loop ──────
print("\n=== excess_mi_windows consistency check (chunked vs per-window loop) ===")
N, n_windows, n_shuf = 120, 15, 10
kin = rng.normal(size=N)
kin_z = gpu_mi.zscore(kin)
kin_stack = gpu_mi.build_kin_stack(kin_z, n_shuf, np.random.default_rng(1))

count_windows_raw = rng.normal(size=(n_windows, N)) + 0.3 * kin[None, :]  # some real dependency
count_windows_z = gpu_mi.zscore(count_windows_raw, axis=1)

mi_obs_chunk, mi_excess_chunk, _, _ = gpu_mi.excess_mi_windows(
    gpu_mi.ksg_mi_batch, kin_stack, count_windows_z, k=3, chunk_size=4)

# reference: loop one window at a time via excess_mi (same shuffles reused via same kin_stack)
mi_obs_loop = np.full(n_windows, np.nan)
mi_excess_loop = np.full(n_windows, np.nan)
for w in range(n_windows):
    X_single = kin_stack
    Y_single = count_windows_z[w:w+1]
    mi_all = gpu_mi.ksg_mi_batch(X_single, Y_single, k=3)
    mi_obs_loop[w] = mi_all[0]
    mi_excess_loop[w] = mi_all[0] - mi_all[1:].mean()

print("  mi_obs    max abs diff:", np.max(np.abs(mi_obs_chunk - mi_obs_loop)))
print("  mi_excess max abs diff:", np.max(np.abs(mi_excess_chunk - mi_excess_loop)))
