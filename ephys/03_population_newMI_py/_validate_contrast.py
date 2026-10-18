# Do the candidate contrasts SEPARATE real maps from shuffled ones?
# Same test that killed the prominence rule: build single-shuffle maps (a real
# null realisation, not an average) and score them exactly like real maps.
# A usable statistic gives real >> shuffled. Anything else is noise-fitting.
import numpy as np, re, time
import mi_io, config, mi_rates
from joblib import Parallel, delayed

N_PER_REGION = 20
N_NULL_MAPS  = 3

src  = open("12_MI_2D_kin_x_neuron.py", encoding="utf-8").read()
head = src[:re.search(r"(?m)^payload, keep = \[\], \[\]", src).start()]
head = re.sub(r"(?m)^N_TEST      = .*$", f"N_TEST      = {N_PER_REGION}", head, count=1)
head = re.sub(r"(?m)^KIN_STEP  = .*$",   "KIN_STEP  = 8", head, count=1)
exec(compile(head, "12h", "exec"))

def _separable(E):
    A = np.nan_to_num(E - np.nanmean(E))
    return A.mean(axis=1, keepdims=True) + A.mean(axis=0, keepdims=True)

def score(E, sd):
    """-> dict of candidate statistics for ONE map."""
    out = {}
    # A: within-map MAD  (current CONTRAST="within")
    zw = _z_within(E)
    out["within_maxz"] = np.nanmax(zw)
    # B: residual after removing the separable part, scaled by the NULL sd
    R  = E - _separable(E)
    zr = R / sd
    out["resid_maxz"]  = np.nanmax(zr)
    # C: residual, scaled by the residual's own robust spread
    m  = np.nanmedian(R); mad = np.nanmedian(np.abs(R - m)) * 1.4826
    out["resid_selfz"] = np.nanmax((R - m) / mad) if mad > 0 else np.nan
    # D: plain null z (reference: known to pass everything)
    out["null_maxz"]   = np.nanmax(E / sd)
    return out

pay, keep2 = {}, []
for pi, i in enumerate(run_idx):
    n = NEURONS[i]; K = KIN[(n["mouse"], n["sess"])]
    cat = n["labels"]["idx_reach_cat"]
    m = np.isfinite(K).all(axis=(1,2)) & np.isin(cat, list(CAT_KEEP))
    if m.sum() < MIN_TRIALS: continue
    C = np.asarray(mi_rates.windowed(n["spikes"], win_starts, WIN_WIDTH,
                                     NEURAL_VAR), float)[m]
    pay[pi] = (K[m], C); keep2.append(pi)
print(f"{len(keep2)} neurons | grid {n_kin}x{n_wins} = {n_kin*n_wins} cells")

def one(K, C, seed, shuffle):
    rs = np.random.default_rng(seed)
    Ks = K[rs.permutation(K.shape[0])] if shuffle else K
    M  = np.full((K.shape[1], C.shape[1]), np.nan)
    SH = np.full((N_SHUF_2D, K.shape[1], C.shape[1]), np.nan)
    perms = [rs.permutation(K.shape[0]) for _ in range(N_SHUF_2D)]
    okw = [w for w in range(C.shape[1]) if np.unique(C[:, w]).size >= 2]
    for ki in range(K.shape[1]):
        X = Ks[:, ki, :]
        for w in okw:
            M[ki, w] = _ksg_md(X, C[:, w], K_KSG, seed + ki*97 + w)
            for s, p in enumerate(perms):
                SH[s, ki, w] = _ksg_md(X[p], C[:, w], K_KSG, seed + ki*97 + w)
    mu = np.nanmean(SH, axis=0); sd = _floor_sd(np.nanstd(SH, axis=0, ddof=1))
    return score(M - mu, sd)

t0 = time.time()
real = Parallel(n_jobs=-1, backend="loky")(
    delayed(one)(pay[pi][0], pay[pi][1], 100+pi, False) for pi in keep2)
print(f"real maps: {(time.time()-t0)/60:.1f} min")
t0 = time.time()
null = Parallel(n_jobs=-1, backend="loky")(
    delayed(one)(pay[pi][0], pay[pi][1], 9000+pi*13+s, True)
    for pi in keep2 for s in range(N_NULL_MAPS))
print(f"null maps: {(time.time()-t0)/60:.1f} min")

print(f"\n{'statistic':16s} {'REAL median':>12s} {'NULL median':>12s} "
      f"{'% real > null p95':>18s} {'verdict':>22s}")
print("-" * 84)
for k in ("null_maxz", "within_maxz", "resid_maxz", "resid_selfz"):
    r = np.array([d[k] for d in real]); r = r[np.isfinite(r)]
    u = np.array([d[k] for d in null]); u = u[np.isfinite(u)]
    if not r.size or not u.size: continue
    thr = np.quantile(u, 0.95)
    frac = 100 * np.mean(r > thr)
    v = ("USELESS" if frac < 10 else "weak" if frac < 40 else
         "usable" if frac < 75 else "GOOD")
    print(f"{k:16s} {np.median(r):12.2f} {np.median(u):12.2f} {frac:17.0f}% {v:>22s}")
print("\n(% real > null p95 is the detection rate at a 5% false-positive rate.)")
