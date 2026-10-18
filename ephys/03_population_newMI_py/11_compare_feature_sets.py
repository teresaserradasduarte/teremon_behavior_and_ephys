# %% ── 0. imports & parameters ────────────────────────────────────────────────
# Which KINEMATIC FEATURE SET should 09 use? Compare candidate sets on EXCESS MI
# (observed - that set's OWN shuffle null), which is the only fair comparison
# across different dimensionality: raw MI rises with d through KSG bias alone, so
# a bigger set would win by construction.
#
# Everything -- neurons, windows, estimator, nulls -- is inherited from
# 09_MI_multiD_kin.py by executing its parameter/helper cells (see cell 1), so
# this script cannot drift from the analysis it is advising. Only FEATURES varies.
import numpy as np
import matplotlib.pyplot as plt
import re, time, pickle
import mi_io, config, mi_rates
from joblib import Parallel, delayed

# ── the sets to compare ──────────────────────────────────────────────────────
# Any feature name 09 understands: x, y, z, PC1, PC2, speed, accel,
# dist_water, speed_water. Add/remove rows freely.
FEATURE_SETS = [
    ("x", "y", "z", "speed", "dist_water", "speed_water"),
    ("x", "y", "z", "speed", "dist_water"),
    ("x", "y", "z", "speed", "speed_water"),
    ("PC1", "PC2", "speed", "dist_water"),
    ("x", "y", "z", "speed"),
    ("x", "y", "z"),
    ("y",),
    ("PC1", "PC2"),
    ("dist_water", "speed_water"),
]

# A one-element tuple needs the trailing comma: ("y",) not ("y"). Without it you
# get the STRING "y", which iterates per CHARACTER -- harmless for a 1-letter name
# but ("speed") would silently become s+p+e+e+d. Normalise and say so.
_fixed = []
for _s in FEATURE_SETS:
    if isinstance(_s, str):
        print(f'NOTE: FEATURE_SETS entry ("{_s}") is a string, not a tuple '
              f'-- read as ("{_s}",). Add the trailing comma to silence this.')
        _fixed.append((_s,))
    else:
        _fixed.append(tuple(_s))
FEATURE_SETS = _fixed

# ── how much to run ──────────────────────────────────────────────────────────
N_TEST      = 500      # neurons PER REGION for a quick answer;
                       # None = THE WHOLE POPULATION (~2872 neurons).
                       # Cost scales with n_neurons x n_sets x (1 + N_SHUF).
N_SHUF      = 10       # null permutations per window per neuron. 5 is usually
                       # enough to rank sets; 10+ if two sets are within ~0.005.
N_JOBS      = -1

# ── locality diagnostic / k sweep ────────────────────────────────────────────
DO_LOCALITY = True     # measure how LOCAL the k-NN neighbourhoods actually are
K_SWEEP     = (2, 3, 4)   # k values to test; () to skip the sweep
K_SWEEP_SET = None     # which feature set to sweep k on; None = the highest-d set
                       # (the one where locality is worst and k matters most)

SAVE_TXT    = True     # table -> OUT_DIR/<OUT_SUBDIR>/
SAVE_FIGS   = True
SAVE_VARS   = True     # full landscapes -> VAR_DIR/<RUN_NAME>.pkl
SHOW_FIGS   = True
OUT_SUBDIR  = "compare_feature_sets"
RUN_NAME    = "featuresets100"


# %% ── 1. inherit 09's parameters, data and estimator ────────────────────────
# 09's cells 0-2 build: NEURONS, REACH_MM, PCA_FIT/KEY_POOL, WATER, the window
# grid, run_idx/regions_run, and the estimator helpers (_ksg_md, _landscape_one,
# _null_one_md). We re-exec them here with two edits:
#   FEATURES -> the UNION of every set below, so the PCA fit and the water loader
#               switch on iff some set actually needs them
#   N_TEST   -> this script's N_TEST
# Line-based regex, so it survives you editing either value in 09.
UNION = tuple(dict.fromkeys(f for s in FEATURE_SETS for f in s))   # ordered, unique

src  = open("09_MI_multiD_kin.py", encoding="utf-8").read()
head = src[:src.index("# %% ── 3.")]
head, n1 = re.subn(r"(?m)^FEATURES\s*=.*$", f"FEATURES    = {UNION!r}", head, count=1)
head, n2 = re.subn(r"(?m)^N_TEST\s*=.*$",   f"N_TEST      = {N_TEST!r}", head, count=1)
if not (n1 and n2):
    raise RuntimeError("could not find FEATURES / N_TEST in 09 -- has it been renamed?")

# 09's cell 0 assigns OUT_SUBDIR, RUN_NAME, SAVE_*, SHOW_FIGS and N_JOBS too, so
# exec'ing it CLOBBERS this script's own settings and the outputs land in 09's
# folder under 09's name. Stash ours, exec, then put them back.
_MINE = {k: globals()[k] for k in
         ("OUT_SUBDIR", "RUN_NAME", "SAVE_TXT", "SAVE_FIGS", "SAVE_VARS",
          "SHOW_FIGS", "N_JOBS")}
exec(compile(head, "09_head", "exec"))
globals().update(_MINE)

print(f"\ncomparing {len(FEATURE_SETS)} feature sets on {len(run_idx)} neurons "
      f"(S={N_SHUF} nulls)")
print(f"union of features: {list(UNION)}")


# %% ── 2. build every candidate column once, per session ─────────────────────
# Cheaper and safer than rebuilding inside the set loop: identical columns are
# then guaranteed to be shared between sets, so differences are only the set.
COLS = {}
for key in REACH_MM:
    rp  = REACH_MM[key]
    vel = savgol_filter(rp, SAVGOL_WIN, SAVGOL_ORD, deriv=1, delta=1/FPS, axis=2)
    spd = np.linalg.norm(vel, axis=1)
    acc = savgol_filter(spd, SAVGOL_WIN, SAVGOL_ORD, deriv=1, delta=1/FPS, axis=1)
    c = {"x": rp[:, 0, EP].mean(1), "y": rp[:, 1, EP].mean(1), "z": rp[:, 2, EP].mean(1),
         "speed": spd[:, EP].mean(1), "accel": acc[:, EP].mean(1)}
    if NEEDS_PCA:
        mu, comp, _ = PCA_FIT[KEY_POOL[key]]
        pcs = (_pca_input(rp) - mu) @ comp.T
        c["PC1"], c["PC2"] = pcs[:, 0], pcs[:, 1]
    if NEEDS_WATER:
        W = WATER[key]
        c["dist_water"]  = W["dist_water"][:, EP].mean(1)
        c["speed_water"] = W["speed_water"][:, EP].mean(1)
    COLS[key] = c

missing = {f for s in FEATURE_SETS for f in s} - set(COLS[list(COLS)[0]])
if missing: raise KeyError(f"features not available: {sorted(missing)}")
print(f"columns built for {len(COLS)} sessions: {sorted(COLS[list(COLS)[0]])}")


# %% ── 3. run every set ──────────────────────────────────────────────────────
RES = {}
t_all = time.time()
for feats in FEATURE_SETS:
    t0 = time.time()
    payload, keep = [], []
    for pi, i in enumerate(run_idx):
        n = NEURONS[i]; c = COLS[(n["mouse"], n["sess"])]
        F   = np.column_stack([c[f] for f in feats])
        cat = n["labels"]["idx_reach_cat"]
        msk = np.isfinite(F).all(1) & np.isin(cat, list(CAT_KEEP))
        if msk.sum() < MIN_TRIALS: payload.append(None); continue
        C = np.asarray(mi_rates.windowed(n["spikes"], win_starts, WIN_WIDTH,
                                         NEURAL_VAR), float)[msk]
        payload.append((F[msk], C)); keep.append(pi)

    obs = Parallel(n_jobs=N_JOBS, backend="loky")(
        delayed(_landscape_one)(payload[pi][0], payload[pi][1], K_KSG, 100+pi)
        for pi in keep)
    nul = Parallel(n_jobs=N_JOBS, backend="loky")(
        delayed(_null_one_md)(payload[pi][0], payload[pi][1], N_SHUF, 1000+pi, K_KSG)
        for pi in keep)

    O = np.full((len(run_idx), n_wins), np.nan); N = np.full_like(O, np.nan)
    for pi, v in zip(keep, obs): O[pi] = v
    for pi, v in zip(keep, nul): N[pi] = v
    RES[feats] = dict(obs=O, null=N, excess=O - N, keep=keep,
                      secs=time.time() - t0)
    print(f"  {'+'.join(feats):40s} done in {time.time()-t0:.0f}s")
print(f"total {time.time()-t_all:.0f}s")


# %% ── 3b. LOCALITY: how local are the k-NN neighbourhoods, really? ──────────
# KSG estimates a LOCAL density from the k nearest neighbours, and its derivation
# assumes the density is roughly constant inside that neighbourhood. With k fixed
# and d growing, the neighbourhood must swell to keep finding k points, until it
# is no longer local and the assumption silently fails.
#
# Measured here, per feature set, on the SAME preprocessed data the estimator sees
# (_prep: each column scaled to unit sd, then the Chebyshev/max norm):
#   span_frac  = 2 * r_k / (per-axis range), i.e. the fraction of the data range
#                covered by the neighbourhood cube. ~0.1 is local; ~0.5 is not.
#   frac_pts   = k / N, the fraction of the data the neighbourhood contains.
#   theory     = (k/N)^(1/(d+1)), the uniform-density prediction, for comparison.
# N is the median number of usable TRIALS per neuron -- the "total amount of data"
# the estimator has to work with, which is what sets how bad this gets.
def locality(X, y, k):
    """-> (span_frac, N, d+1) for one neuron's (X, y) at a given k."""
    r  = np.random.default_rng(0)
    Xs = _prep(np.asarray(X, float).reshape(len(y), -1), r)
    ys = _prep(np.asarray(y, float).reshape(-1, 1), r)
    J  = np.hstack([Xs, ys])
    nn = NearestNeighbors(metric="chebyshev", n_neighbors=k).fit(J)
    rk = np.median(nn.kneighbors()[0][:, -1])          # median k-th NN distance
    span = np.median(J.max(axis=0) - J.min(axis=0))    # median per-axis range
    return float(2*rk/span), J.shape[0], J.shape[1]

LOC = {}
if DO_LOCALITY:
    lines_loc = ["", "LOCALITY of the k-NN neighbourhoods (k=%d)" % K_KSG,
                 "  span_frac = fraction of the data RANGE the neighbourhood covers",
                 "              (~0.1 local, >0.4 the local-density assumption is gone)",
                 "  N         = median usable trials per neuron", ""]
    hdr_loc = (f"{'features':40s} {'d+1':>4s} {'N':>6s} {'span_frac':>10s} "
               f"{'k/N':>8s} {'theory':>8s}")
    lines_loc += [hdr_loc, "-" * len(hdr_loc)]
    for feats in FEATURE_SETS:
        vals, Ns = [], []
        for pi in RES[feats]["keep"][:60]:            # a sample is plenty
            i = run_idx[pi]; n = NEURONS[i]; c = COLS[(n["mouse"], n["sess"])]
            F   = np.column_stack([c[f] for f in feats])
            cat = n["labels"]["idx_reach_cat"]
            msk = np.isfinite(F).all(1) & np.isin(cat, list(CAT_KEEP))
            C = np.asarray(mi_rates.windowed(n["spikes"], win_starts, WIN_WIDTH,
                                             NEURAL_VAR), float)[msk]
            wi = C.shape[1] // 2                       # a mid landscape window
            if np.unique(C[:, wi]).size < 2: continue
            sf, N_, dp1 = locality(F[msk], C[:, wi], K_KSG)
            vals.append(sf); Ns.append(N_)
        sf = float(np.median(vals)); N_ = int(np.median(Ns)); d1 = len(feats) + 1
        LOC[feats] = dict(span_frac=sf, N=N_, dp1=d1)
        lines_loc.append(f"{'+'.join(feats):40s} {d1:4d} {N_:6d} {sf:10.3f} "
                         f"{K_KSG/N_:8.4f} {(K_KSG/N_)**(1/d1):8.3f}")
    print("\n".join(lines_loc))

# --- does changing k help? ---------------------------------------------------
# k trades two errors against each other:
#   larger k  -> smoother density estimate, LOWER variance, but a BIGGER
#                neighbourhood, so WORSE locality and more bias
#   smaller k -> better locality, but noisier
# There is no k that fixes high-d: k>=1 always, and the radius grows with d
# whatever k is. The sweep shows how steep that trade-off is on THIS data.
KSWEEP = {}
if K_SWEEP:
    ks_set = K_SWEEP_SET or max(FEATURE_SETS, key=len)
    print(f"\nk sweep on {'+'.join(ks_set)} (d={len(ks_set)}):")
    payload, keep = [], []
    for pi, i in enumerate(run_idx):
        n = NEURONS[i]; c = COLS[(n["mouse"], n["sess"])]
        F   = np.column_stack([c[f] for f in ks_set])
        cat = n["labels"]["idx_reach_cat"]
        msk = np.isfinite(F).all(1) & np.isin(cat, list(CAT_KEEP))
        if msk.sum() < MIN_TRIALS: payload.append(None); continue
        C = np.asarray(mi_rates.windowed(n["spikes"], win_starts, WIN_WIDTH,
                                         NEURAL_VAR), float)[msk]
        payload.append((F[msk], C)); keep.append(pi)

    hdr_k = (f"  {'k':>3s} {'span_frac':>10s} {'peak excess':>12s} "
             f"{'peak raw':>10s} {'peak null':>10s}")
    lines_k = ["", f"k SWEEP on {'+'.join(ks_set)} (d+1={len(ks_set)+1}, "
                   f"N~{LOC.get(ks_set, {}).get('N', '?')} trials)",
               hdr_k, "-" * len(hdr_k)]
    for k in K_SWEEP:
        obs = Parallel(n_jobs=N_JOBS, backend="loky")(
            delayed(_landscape_one)(payload[pi][0], payload[pi][1], k, 100+pi)
            for pi in keep)
        nul = Parallel(n_jobs=N_JOBS, backend="loky")(
            delayed(_null_one_md)(payload[pi][0], payload[pi][1], N_SHUF, 1000+pi, k)
            for pi in keep)
        O = np.nanmean(np.vstack(obs), axis=0); N = np.nanmean(np.vstack(nul), axis=0)
        j = int(np.nanargmax(O - N))
        sf = float(np.median([locality(payload[pi][0],
                                       payload[pi][1][:, payload[pi][1].shape[1]//2], k)[0]
                              for pi in keep[:40]]))
        KSWEEP[k] = dict(excess=float((O-N)[j]), obs=float(O[j]), null=float(N[j]),
                         span_frac=sf)
        lines_k.append(f"  {k:3d} {sf:10.3f} {(O-N)[j]:12.4f} {O[j]:10.4f} {N[j]:10.4f}")
        print(lines_k[-1], flush=True)
    bk = max(KSWEEP, key=lambda k: KSWEEP[k]["excess"])
    lines_k += ["", f"  best k = {bk} ({KSWEEP[bk]['excess']:.4f} bits); "
                    f"k={K_KSG} gives {KSWEEP.get(K_KSG, {}).get('excess', float('nan')):.4f}"]
else:
    lines_k = []


# %% ── 3c. LEAD/LAG peak distributions, CP vs CB side by side ────────────────
# The same quantity 09 prints at the end of its per-mouse cell, but with one ROW
# per feature set and CP / CB as separate COLUMN GROUPS, so the lead/lag result
# can be read off across feature choices at a glance.
#
# Peaks are found on the RAW landscape (not excess), exactly as 09 does: the
# d-dependent bias is near-flat across windows, so it shifts the height but not
# the TIME of the peak, and using raw keeps these numbers directly comparable to
# 09's own output. "n sig" applies BOTH criteria 09 applies -- a prominent peak
# (PEAK_PROM_RATIO / PEAK_RUNNERUP_MAX) AND the peak-window shuffle test.
DO_LEADLAG = True      # costs one extra shuffle test per neuron per set

lo_ll, hi_ll = PEAK_RANGE if PEAK_RANGE is not None else (win_centers[0], win_centers[-1])
mask_ll = (win_centers >= lo_ll) & (win_centers <= hi_ll)

LEADLAG = {}
if DO_LEADLAG:
    for feats in FEATURE_SETS:
        keep = RES[feats]["keep"]
        O    = RES[feats]["obs"]
        # rebuild the payload (cheap) so _sig_one can permute the real trials
        pay = {}
        for pi in keep:
            i = run_idx[pi]; n = NEURONS[i]; c = COLS[(n["mouse"], n["sess"])]
            F   = np.column_stack([c[f] for f in feats])
            cat = n["labels"]["idx_reach_cat"]
            msk = np.isfinite(F).all(1) & np.isin(cat, list(CAT_KEEP))
            C = np.asarray(mi_rates.windowed(n["spikes"], win_starts, WIN_WIDTH,
                                             NEURAL_VAR), float)[msk]
            pay[pi] = (F[msk], C)
        sig = Parallel(n_jobs=N_JOBS, backend="loky")(
            delayed(_sig_one)(pay[pi][0], pay[pi][1], O[pi], K_KSG, N_SHUF_LEADLAG,
                              300+pi, SMOOTH_WIN, PEAK_PROM_RATIO, PEAK_RUNNERUP_MAX)
            for pi in keep)
        sig = dict(zip(keep, sig))

        row = {}
        for reg in ("CP", "CB"):
            idx = [pi for pi in keep if regions_run[pi] == reg]
            pks = []
            for pi in idx:
                if not sig[pi]: continue
                r_ = O[pi]
                if not np.any(~np.isnan(r_)): continue
                pk = _find_peak(np.where(mask_ll, _smooth(r_, SMOOTH_WIN), np.nan))
                if pk is not None: pks.append(win_centers[pk[0]] * 1000)
            pks = np.asarray(pks)
            row[reg] = dict(n=len(idx), nsig=len(pks),
                            med=float(np.median(pks)) if pks.size else np.nan,
                            q1=float(np.percentile(pks, 25)) if pks.size else np.nan,
                            q3=float(np.percentile(pks, 75)) if pks.size else np.nan)
        row["lead"] = row["CB"]["med"] - row["CP"]["med"]
        LEADLAG[feats] = row
        print(f"  lead/lag {'+'.join(feats):38s} "
              f"CP {row['CP']['med']:6.0f} ms (n={row['CP']['nsig']})   "
              f"CB {row['CB']['med']:6.0f} ms (n={row['CB']['nsig']})   "
              f"CB-CP {row['lead']:+6.0f} ms", flush=True)

    # --- the table ---
    h1 = (f"{'':40s} {'--------------- CP ---------------':>34s}  "
          f"{'--------------- CB ---------------':>34s}  {'':>9s}")
    h2 = (f"{'feature set':40s} {'n':>5s} {'n sig':>6s} {'median':>8s} "
          f"{'IQR':>13s}  {'n':>5s} {'n sig':>6s} {'median':>8s} {'IQR':>13s}  "
          f"{'CB-CP':>9s}")
    lines_ll = ["", "LEAD / LAG: peak-MI time of the significantly-peaked neurons (ms "
                    "re reach)", f"peaks on the RAW landscape, smoothed win={SMOOTH_WIN}, "
                    f"prominence>={PEAK_PROM_RATIO}, shuffle test S={N_SHUF_LEADLAG}",
                "negative CB-CP = CB LEADS", "", h1, h2, "-" * len(h2)]
    for feats in FEATURE_SETS:
        r = LEADLAG[feats]
        cells = []
        for reg in ("CP", "CB"):
            d_ = r[reg]
            cells.append(f"{d_['n']:5d} {d_['nsig']:6d} {d_['med']:8.0f} "
                         f"[{d_['q1']:5.0f},{d_['q3']:5.0f}]")
        lines_ll.append(f"{'+'.join(feats):40s} {cells[0]}  {cells[1]}  "
                        f"{r['lead']:+9.0f}")
    print("\n".join(lines_ll))
else:
    lines_ll = []


# %% ── 4. the table ──────────────────────────────────────────────────────────
def _peak(E, sel):
    """(peak mean excess MI, its window centre in ms) over neurons `sel`."""
    if not sel: return (np.nan, np.nan)
    m = np.nanmean(E[sel], axis=0)
    j = int(np.nanargmax(m))
    return float(m[j]), float(win_centers[j] * 1000)

rows = []
for feats, R in RES.items():
    E, keep = R["excess"], R["keep"]
    cp = _peak(E, [p for p in keep if regions_run[p] == "CP"])
    cb = _peak(E, [p for p in keep if regions_run[p] == "CB"])
    both = _peak(E, keep)
    rows.append(("+".join(feats), len(feats), len(keep),
                 cp[0], cb[0], both[0], cp[1], cb[1]))

hdr = (f"{'features':40s} {'d':>2s} {'n':>5s} {'CP excess':>10s} {'CB excess':>10s} "
       f"{'BOTH':>8s} {'CP t':>6s} {'CB t':>6s}")
txt = [f"feature-set comparison  |  {len(run_idx)} neurons"
       f"{'' if N_TEST is None else f' ({N_TEST}/region)'}, S={N_SHUF} nulls, "
       f"KSG k={K_KSG}, {NEURAL_VAR}, mm",
       f"EXCESS MI = observed - that set's own per-window shuffle null (bits)", "",
       hdr, "-" * len(hdr)]
for r in sorted(rows, key=lambda r: -r[5]):
    txt.append(f"{r[0]:40s} {r[1]:2d} {r[2]:5d} {r[3]:10.4f} {r[4]:10.4f} "
               f"{r[5]:8.4f} {r[6]:6.0f} {r[7]:6.0f}")
best = max(rows, key=lambda r: r[5])
txt += ["", f"BEST: {best[0]}  ({best[5]:.4f} bits)"]

# pairwise: what does each set cost/buy against the best one?
txt += ["", "vs BEST (negative = worse):"]
for r in sorted(rows, key=lambda r: -r[5])[1:]:
    txt.append(f"  {r[0]:40s} {r[5]-best[5]:+.4f} bits  (d {r[1]} vs {best[1]})")

txt += (lines_loc if DO_LOCALITY else []) + lines_k + lines_ll
REPORT = "\n".join(txt)
print("\n" + REPORT)
if SAVE_TXT:
    p = config.fig_dir(OUT_SUBDIR) / f"{RUN_NAME}_table.txt"
    p.write_text(REPORT + "\n", encoding="utf-8")
    print(f"\ntable -> {p}")


# %% ── 5. figures ────────────────────────────────────────────────────────────
fdir  = config.fig_dir(OUT_SUBDIR)
order = [f for f, _ in sorted(RES.items(), key=lambda kv: -_peak(kv[1]["excess"],
                                                                kv[1]["keep"])[0])]
labels = ["+".join(f) for f in order]

# (a) peak excess MI per set, per region
fig, ax = plt.subplots(figsize=(max(8, 1.5*len(order)), 4.5))
w = 0.38
for k, reg in enumerate(("CP", "CB")):
    vals = [_peak(RES[f]["excess"], [p for p in RES[f]["keep"]
                                     if regions_run[p] == reg])[0] for f in order]
    ax.bar(np.arange(len(order)) + (k-0.5)*w, vals, w, color=CLR[reg], label=reg)
ax.set_xticks(np.arange(len(order)))
ax.set_xticklabels([l.replace("+", "\n+") for l in labels], fontsize=7)
ax.set_ylabel("peak mean EXCESS MI (bits)")
ax.set_title(f"Feature sets ranked  —  {len(run_idx)} neurons, S={N_SHUF} nulls")
ax.spines[["top", "right"]].set_visible(False)
ax.legend(fontsize=8)
plt.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / f"{RUN_NAME}_ranking.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)

# (b) the excess landscapes, so a set that wins on shape not just peak is visible
fig, axes = plt.subplots(1, 2, figsize=(13, 4), sharey=True)
cmap = plt.get_cmap("tab10")
for ax, reg in zip(axes, ("CP", "CB")):
    for i, f in enumerate(order):
        sel = [p for p in RES[f]["keep"] if regions_run[p] == reg]
        if not sel: continue
        ax.plot(win_centers, np.nanmean(RES[f]["excess"][sel], axis=0),
                color=cmap(i % 10), lw=1.8, label="+".join(f))
    ax.axvline(0, color="k", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.8)
    ax.set(xlabel="window center re reach (s)", title=reg)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("mean EXCESS MI (bits)")
axes[0].legend(fontsize=7)
plt.suptitle(f"Excess-MI landscape per feature set  —  KSG k={K_KSG}, {NEURAL_VAR}",
             fontsize=11)
plt.tight_layout()
if SAVE_FIGS: fig.savefig(fdir / f"{RUN_NAME}_landscapes.png", dpi=150, bbox_inches="tight")
if SHOW_FIGS: plt.show()
plt.close(fig)
# (c) the lead/lag table rendered as a figure -- same text as the report, so the
#     two can never disagree
if DO_LEADLAG and lines_ll:
    # ~0.075 in per char and ~0.17 in per line at fontsize 8.5 monospace, so the
    # canvas ends up the size of the text instead of mostly white space
    fig, ax = plt.subplots(figsize=(max(9, 0.075*max(len(l) for l in lines_ll)),
                                    0.17*len(lines_ll) + 0.3))
    ax.axis("off")
    ax.text(0.0, 1.0, "\n".join(lines_ll), family="monospace", fontsize=8.5,
            va="top", ha="left", transform=ax.transAxes)
    plt.tight_layout()
    if SAVE_FIGS:
        fig.savefig(fdir / f"{RUN_NAME}_leadlag_table.png", dpi=150,
                    bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

    # (d) the same numbers as a plot: median peak time per set, CP vs CB
    fig, ax = plt.subplots(figsize=(max(8, 1.3*len(FEATURE_SETS)), 4.5))
    xs = np.arange(len(FEATURE_SETS))
    for k, reg in enumerate(("CP", "CB")):
        med = [LEADLAG[f][reg]["med"] for f in FEATURE_SETS]
        q1  = [LEADLAG[f][reg]["q1"]  for f in FEATURE_SETS]
        q3  = [LEADLAG[f][reg]["q3"]  for f in FEATURE_SETS]
        err = np.abs(np.vstack([np.array(med)-np.array(q1),
                                np.array(q3)-np.array(med)]))
        ax.errorbar(xs + (k-0.5)*0.16, med, yerr=err, fmt="o", capsize=3,
                    color=CLR[reg], label=f"{reg} (median, IQR)")
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xticks(xs)
    ax.set_xticklabels(["+".join(f).replace("+", "\n+") for f in FEATURE_SETS],
                       fontsize=7)
    ax.set_ylabel("peak MI time re reach (ms)")
    ax.set_title("Lead/lag per feature set  —  significantly-peaked neurons")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8)
    plt.tight_layout()
    if SAVE_FIGS:
        fig.savefig(fdir / f"{RUN_NAME}_leadlag.png", dpi=150, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

print(f"figures -> {fdir}")

if SAVE_VARS:
    p = config.var_path(f"{RUN_NAME}.pkl")
    with open(p, "wb") as f:
        pickle.dump(dict(sets=[list(s) for s in FEATURE_SETS],
                         excess={"+".join(k): v["excess"] for k, v in RES.items()},
                         obs={"+".join(k): v["obs"] for k, v in RES.items()},
                         null={"+".join(k): v["null"] for k, v in RES.items()},
                         win_centers=win_centers, regions=regions_run,
                         uids=[uids[i] for i in run_idx],
                         locality={"+".join(k): v for k, v in LOC.items()},
                         k_sweep=KSWEEP,
                         leadlag={"+".join(k): v for k, v in LEADLAG.items()},
                         N_TEST=N_TEST, N_SHUF=N_SHUF, K=K_KSG,
                         NEURAL_VAR=NEURAL_VAR), f, protocol=4)
    print(f"saved -> {p}")

# %%
