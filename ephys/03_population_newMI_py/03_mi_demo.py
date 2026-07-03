# %% ── imports & constants ────────────────────────────────────────────────────
import numpy as np
import matplotlib.pyplot as plt
import mi_io, config
from sklearn.feature_selection import mutual_info_regression, mutual_info_classif

if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

CLR      = {"CP": np.array([106,176,134])/255, "CB": np.array([154,123,193])/255}
NAT2BIT  = 1 / np.log(2)
rng      = np.random.default_rng(0)

win        = (-0.10, 0.50)   # spike-count integration window (s)
raster_win = (-0.50, 1.00)   # raster display window (s)
COORD      = 1               # kinematic coordinate: 0=x, 1=y med-lat, 2=z
FRAMES     = (147, 167)      # endpoint window, MATLAB 1-indexed
FPS        = 400             # camera fps
CAT_KEEP   = {1, 2}         # reach categories to include
N_SHUF     = 200             # shuffle permutations
K_LIST     = [3, 5, 7]      # k values for KSG sensitivity check
SAVE_FIGS  = True          # save figures to disk (config.fig_dir("compare"))

# ── tiny helpers (reused across cells) ────────────────────────────────────────
def ksg_mi(kin, count, k=3):
    """KSG MI in bits. count treated as continuous (+ infinitesimal jitter)."""
    X = kin.reshape(-1, 1)
    return mutual_info_regression(X, count + rng.normal(0, 1e-9, count.size),
                                  n_neighbors=k, random_state=0)[0] * NAT2BIT

def bin_mi(kin, count, nb=16):
    """Plug-in (binning) MI in bits. Biased — for comparison only."""
    edges = np.unique(np.quantile(kin, np.linspace(0, 1, nb + 1)))
    kb = np.clip(np.digitize(kin, edges[1:-1]), 0, len(edges) - 2)
    cv = np.unique(count); cm = {v: i for i, v in enumerate(cv)}
    ci = np.array([cm[v] for v in count])
    J  = np.zeros((len(cv), len(edges)-1)); np.add.at(J, (ci, kb), 1.0); P = J / J.sum()
    pc = P.sum(1, keepdims=True); pk = P.sum(0, keepdims=True); nz = P > 0
    return float(np.sum(P[nz] * np.log2(P[nz] / (pc @ pk + 1e-300)[nz])))

def ross_mi(kin, count, k=3):
    """Ross MI in bits. count stays discrete (classif)."""
    X = kin.reshape(-1, 1)
    return mutual_info_classif(X, count.astype(int),
                               n_neighbors=k, random_state=0)[0] * NAT2BIT


# %% ── pick representative neurons ────────────────────────────────────────────
# For each neuron: compute mean count + KSG MI on cat-1&2 trials.
# Neurons with all-zero counts are skipped (MI=0, not interesting).
# This cell takes ~1-2 min (one MI call per neuron); result stored in neu_stats.
neu_stats = []   # (region, idx, mean_count, mi)
for i, n in enumerate(NEURONS):
    if i % 200 == 0: print(f"  scanning {i}/{len(NEURONS)} ...", flush=True)
    kin  = mi_io.endpoint_kin(n, REACH_PX, coord=COORD, frames=FRAMES)
    cat  = n["labels"]["idx_reach_cat"]
    mask = np.isfinite(kin) & np.isin(cat, list(CAT_KEEP))
    if mask.sum() < 20:
        neu_stats.append((n["region"], i, 0.0, 0.0)); continue
    c_m   = mi_io.counts(n, win)[mask]
    kin_m = kin[mask]
    mc    = c_m.mean()
    mi    = ksg_mi(kin_m, c_m, k=3) if np.unique(c_m).size >= 2 else 0.0
    neu_stats.append((n["region"], i, mc, mi))

neu_stats.sort(key=lambda r: r[2])   # sort by mean count for sparsest/median/densest picks
cp = [r for r in neu_stats if r[0] == "CP"]
cb = [r for r in neu_stats if r[0] == "CB"]

picks = []
if cp: picks += [("CP sparsest", cp[0][1]), ("CP median", cp[len(cp)//2][1])]
if cb: picks += [("CB median",   cb[len(cb)//2][1]), ("CB densest", cb[-1][1])]

# add highest-MI neuron per region (for k-sensitivity & shuffle tests)
cp_by_mi = sorted(cp, key=lambda r: r[3], reverse=True)
cb_by_mi = sorted(cb, key=lambda r: r[3], reverse=True)
if cp_by_mi: picks.append(("CP best MI", cp_by_mi[0][1]))
if cb_by_mi: picks.append(("CB best MI", cb_by_mi[0][1]))

print("\npicks:")
for label, idx in picks:
    r = neu_stats[idx]
    print(f"  {label:13s}  {NEURONS[idx]['uid']:35s}  mean_count={r[2]:.2f}  MI={r[3]:.4f}")


# %% ── demo figure: 5 panels per neuron ───────────────────────────────────────
for label, idx in picks:
    n   = NEURONS[idx]
    reg = n["region"]
    col = CLR[reg]

    # ── data preparation ──
    c   = mi_io.counts(n, win)
    kin = mi_io.endpoint_kin(n, REACH_PX, coord=COORD, frames=FRAMES)
    cat = n["labels"]["idx_reach_cat"]
    mask = np.isfinite(kin) & np.isin(cat, list(CAT_KEEP))

    kin_m = kin[mask]
    c_m   = c[mask]
    sort_i = np.argsort(kin_m)      # trial order ↑ kin (shared across panels)

    rp   = REACH_PX[(n["mouse"], n["sess"])]
    traj = rp[mask, COORD, :]       # (n_kept, 201)
    t_ms = (np.arange(201) - FRAMES[0] + 1) / FPS * 1000   # ms re endpoint-window start

    maxc = max(int(c_m.max()), 1)
    mb   = bin_mi(kin_m, c_m)
    mk   = ksg_mi(kin_m, c_m, k=3)
    mr   = ross_mi(kin_m, c_m, k=3)

    # ── figure ──
    fig, axes = plt.subplots(1, 5, figsize=(18, 4))
    fig.suptitle(f"{label}  |  {n['uid']}  (cat 1&2, n={mask.sum()})\n"
                 f"mean count={c_m.mean():.2f}  {100*np.mean(c_m==0):.0f}% empty  "
                 f"bin={mb:.3f}  KSG={mk:.3f}  Ross={mr:.3f} bits", fontsize=9)

    # [0] raster sorted by kin
    ax = axes[0]
    spk_trials = [n["spikes"][i] for i in np.where(mask)[0]]
    for rank, ti in enumerate(sort_i):
        spk_w = spk_trials[ti]
        spk_w = spk_w[(spk_w >= raster_win[0]) & (spk_w < raster_win[1])]
        ax.scatter(spk_w, np.full(len(spk_w), rank), s=1.5, c=[col], linewidths=0)
    ax.axvspan(win[0], win[1], alpha=0.12, color="k", label="count win")
    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.set(xlabel="time re reach (s)", ylabel="trial (↑ kin)",
           title="raster", xlim=raster_win)
    ax.legend(fontsize=7, loc="upper right")

    # [1] trajectory heatmap sorted by kin
    ax = axes[1]
    vmax = np.nanpercentile(np.abs(traj), 98)
    im   = ax.imshow(traj[sort_i], aspect="auto", origin="lower",
                     extent=[t_ms[0], t_ms[-1], 0, mask.sum()],
                     cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ep_end_ms = (FRAMES[1] - FRAMES[0]) / FPS * 1000
    ax.axvspan(0, ep_end_ms, alpha=0.25, color="k")
    plt.colorbar(im, ax=ax, label=f"coord {COORD} (px)")
    ax.set(xlabel="time re endpoint win (ms)", ylabel="trial (↑ kin)",
           title="trajectory heatmap")

    # [2] p(kin)
    ax = axes[2]
    ax.hist(kin_m, bins=25, color=col, edgecolor="white")
    ax.set(xlabel="endpoint y (px)", ylabel="# trials", title="p(kin)")

    # [3] p(count)
    ax = axes[3]
    ax.hist(c_m, bins=np.arange(-0.5, maxc + 1.5), color=col, edgecolor="white")
    ax.set(xlabel="spike count", ylabel="# trials", title="p(count)")

    # [4] count vs endpoint scatter
    ax = axes[4]
    ax.scatter(kin_m, c_m, s=8, color=col, alpha=0.35, linewidths=0)
    ax.set(xlabel="endpoint y (px)", ylabel="spike count",
           title="count vs endpoint", yticks=np.arange(0, maxc + 1))

    plt.tight_layout()
    if SAVE_FIGS:
        fig.savefig(config.fig_dir("compare") / f"mi_demo_{n['uid'].replace('|','_')}.png",
                    dpi=150, bbox_inches="tight")
    plt.show()
    print(f"{label:13s}  {n['uid']:35s}  bin={mb:.3f}  KSG={mk:.3f}  Ross={mr:.3f}")


# %% ── k-sensitivity: KSG MI for k = 3, 5, 7 ─────────────────────────────────
# result table: rows = neurons, cols = k values
k_results = np.full((len(picks), len(K_LIST)), np.nan)

for pi, (label, idx) in enumerate(picks):
    n   = NEURONS[idx]
    c   = mi_io.counts(n, win)
    kin = mi_io.endpoint_kin(n, REACH_PX, coord=COORD, frames=FRAMES)
    cat = n["labels"]["idx_reach_cat"]
    mask = np.isfinite(kin) & np.isin(cat, list(CAT_KEEP))
    kin_m = kin[mask]; c_m = c[mask]

    for ki, k in enumerate(K_LIST):
        k_results[pi, ki] = ksg_mi(kin_m, c_m, k=k)

    print(f"{label:13s}  " + "  ".join(f"k={k}: {k_results[pi,ki]:.4f}"
                                        for ki, k in enumerate(K_LIST)))

# figure
labels_short = [lbl for lbl, _ in picks]
x = np.arange(len(picks))
width = 0.25

fig, ax = plt.subplots(figsize=(8, 4))
for ki, k in enumerate(K_LIST):
    ax.bar(x + ki*width, k_results[:, ki], width, label=f"k={k}", alpha=0.85)

ax.set_xticks(x + width); ax.set_xticklabels(labels_short, rotation=15, ha="right")
ax.set(ylabel="MI (bits)", title=f"KSG MI sensitivity to k  |  win={win}  cat 1&2")
ax.legend()
ax.axhline(0, color="k", lw=0.5)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir("compare") / "mi_k_sensitivity.png", dpi=150, bbox_inches="tight")
plt.show()


# %% ── shuffle test: 200 permutations, p95 threshold ─────────────────────────
shuf_results = {}   # label -> {"real": float, "null": array}

for label, idx in picks:
    n   = NEURONS[idx]
    c   = mi_io.counts(n, win)
    kin = mi_io.endpoint_kin(n, REACH_PX, coord=COORD, frames=FRAMES)
    cat = n["labels"]["idx_reach_cat"]
    mask = np.isfinite(kin) & np.isin(cat, list(CAT_KEEP))
    kin_m = kin[mask]; c_m = c[mask]

    real_mi = ksg_mi(kin_m, c_m, k=3)

    null = np.array([ksg_mi(rng.permutation(kin_m), c_m, k=3) for _ in range(N_SHUF)])

    p95 = np.percentile(null, 95)
    sig = real_mi > p95
    shuf_results[label] = {"real": real_mi, "null": null, "p95": p95, "sig": sig}
    print(f"{label:13s}  real={real_mi:.4f}  p95={p95:.4f}  {'* sig' if sig else 'ns'}")

# figure: one panel per neuron
fig, axes = plt.subplots(1, len(picks), figsize=(4*len(picks), 3.8), sharey=False)
axes = np.atleast_1d(axes)

for ax, (label, idx) in zip(axes, picks):
    res = shuf_results[label]
    col = CLR[NEURONS[idx]["region"]]

    ax.hist(res["null"], bins=30, color="lightgray", edgecolor="white", label="shuffle null")
    ax.axvline(res["p95"],  color="gray",  lw=1.5, ls="--", label="p95 null")
    ax.axvline(res["real"], color=col,     lw=2,   ls="-",  label=f"real MI")
    sig_str = "* p<0.05" if res["sig"] else "ns"
    ax.set(xlabel="MI (bits)", ylabel="# shuffles" if ax is axes[0] else "",
           title=f"{label}\n{sig_str}  real={res['real']:.4f}  p95={res['p95']:.4f}")
    ax.legend(fontsize=7)

fig.suptitle(f"Shuffle test (n={N_SHUF})  |  KSG k=3  |  win={win}  cat 1&2", fontsize=9)
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir("compare") / "mi_shuffle_test.png", dpi=150, bbox_inches="tight")
plt.show()
