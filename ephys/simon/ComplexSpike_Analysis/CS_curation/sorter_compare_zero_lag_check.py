# -*- coding: utf-8 -*-
"""
sorter_compare_zero_lag_check.py — Do CS spikes also appear inside nearby units (same event, two clusters)?

Question (from the step0 v3 inspection of 5_FerreroRocher S6): half its CS have a nearby unit with a zero-lag
CCG peak > 2. If CS spikes (or their spikelets) are also landing in neighbouring clusters, those neighbours are
CS-contaminated and the "no SS pause" result is partly a sorting artefact.

For every extra-good CB PkC_cs unit (step0 role == "CS") and every non-CS unit within NEAR_UM, any label (v5):
  window       = overlap of the two units' step0 v3 stable spans (unit_table.csv span_start/span_end)
  p_zero       = fraction of CS spikes with a neighbour spike within ±ZERO_HALF
  p_zero_null  = 1 - exp(-fr_nb · 2·ZERO_HALF)      (Poisson chance for that neighbour rate)
  zero_excess  = p_zero - p_zero_null                (CS spikes "explained" beyond chance)
  frac_nb_cs   = fraction of the neighbour's spikes within ±ZERO_HALF of a CS spike
  post_ratio   = neighbour spikes in (ZERO_HALF, POST_END] after CS / expected from its rate
                 (> 1: excess after the CS = spikelet capture; < 1: a pause, as a real SS would show)
  peak_lag_ms  = lag of the max of the 0.1 ms CCG within ±10 ms
Per CS: union p_zero across all its neighbours vs the union chance, and the best (max zero_excess) neighbour.
Shared sinks: neighbours whose zero_excess > SINK_EXCESS for >= 2 CS.

Spike times are read (read-only, mmap) from the sorter folder; units and spans come from the step0 v3 batch.
Outputs (local only): CS_curation/_sorter_comparison/zero_lag_check_<BATCH_TAG>/  (first run: zero_lag_check/, v3batch input)
Batch tag: env STEP0_BATCH_TAG (default "v5batch"), shared with run_step0_batch.py and sorter_compare_step0v3.py.
GOLDEN RULE: Dropbox is READ-ONLY.
"""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

DROPBOX_ROOT = Path(r"D:\Learning Lab Dropbox")
CODE_DIR = Path(r"C:\Users\simon\Desktop\Code\EPhys\ComplexSpike_Analysis")
LOCAL = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation")
CMP = LOCAL / "_sorter_comparison"
BATCH_TAG = os.environ.get("STEP0_BATCH_TAG", "v5batch")
RUN_DIR = f"step0_explore_bhvspan_{BATCH_TAG}"
NEAR_UM = 100.0
ZERO_HALF = 0.5e-3
POST_END = 8e-3
SINK_EXCESS = 0.1
STRONG_EXCESS = 0.2          # "a large share of this CS's spikes is also in a neighbour"
INSPECT = [("5_FerreroRocher", "S6"), ("3_Toblerone", "S5")]
SORTER_COLOR = {"ks4": "#2a78d6", "ibl": "#eb6834"}
SORTER_NAME = {"ks4": "Kilosort4", "ibl": "IBL sorter"}


def safe_out_dir(p):
    p = Path(p).resolve()
    for forbidden in (DROPBOX_ROOT, CODE_DIR):
        try:
            p.relative_to(forbidden.resolve())
        except ValueError:
            continue
        raise RuntimeError(f"Refusing to write under {forbidden}: {p}")
    p.mkdir(parents=True, exist_ok=True)
    return p


OUT = safe_out_dir(CMP / f"zero_lag_check_{BATCH_TAG}")
safe_out_dir(OUT / "figs")


def count_within(ref, tgt, lo, hi):
    """per ref spike: number of tgt spikes with lag in [lo, hi]"""
    return np.searchsorted(tgt, ref + hi, "right") - np.searchsorted(tgt, ref + lo, "left")


def peak_lag(ref, tgt, bin_s=0.1e-3, win=10e-3):
    lo = np.searchsorted(tgt, ref - win, "left")
    hi = np.searchsorted(tgt, ref + win, "left")
    n = hi - lo
    if n.sum() == 0:
        return np.nan
    ri = np.repeat(np.arange(len(ref)), n)
    j = lo[ri] + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
    d = tgt[j] - ref[ri]
    c, e = np.histogram(d, np.arange(-win, win + bin_s / 2, bin_s))
    return float((e[np.argmax(c)] + bin_s / 2) * 1e3)


tab = pd.read_csv(CMP / "session_table.csv")
pair_rows, cs_rows = [], []
for row in tab.itertuples():
    ut_p = LOCAL / f"{row.mouse}_{row.R}" / RUN_DIR / "tables" / "unit_table.csv"
    if not ut_p.exists():
        continue
    ut = pd.read_csv(ut_p).set_index("phyID")
    cs_ids = ut.index[ut["role"] == "CS"]
    if len(cs_ids) == 0:
        continue
    nb_pool = ut[ut["role"] != "CS"]   # any label (v5: no label gate, same as step0 tier 2)
    sd = Path(row.sorter_dir)
    params = {}
    for line in (sd / "params.py").read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                params[k.strip()] = eval(v.strip(), {}, {})
            except Exception:
                pass
    SR = float(params["sample_rate"])
    st_all = np.asarray(np.load(sd / "spike_times.npy", mmap_mode="r")).ravel()
    clu = np.asarray(np.load(sd / "spike_clusters.npy", mmap_mode="r")).ravel().astype(np.int64)
    # neighbours of any CS
    need = set(cs_ids)
    near = {}
    for c in cs_ids:
        d = np.hypot(nb_pool["x"] - ut.at[c, "x"], nb_pool["y"] - ut.at[c, "y"])
        near[c] = d[d <= NEAR_UM]
        need |= set(near[c].index)
    sel = np.isin(clu, list(need))
    t_s, c_s = st_all[sel].astype(np.float64) / SR, clu[sel]
    o = np.argsort(c_s, kind="stable")
    t_s, c_s = t_s[o], c_s[o]
    u, s0, n = np.unique(c_s, return_index=True, return_counts=True)
    ST = {int(a): np.sort(t_s[b:b + k]) for a, b, k in zip(u, s0, n)}
    del st_all, clu

    for c in cs_ids:
        cs_full = ST.get(int(c), np.empty(0))
        any_hit, chance_rate, best = None, 0.0, None
        a0, a1 = ut.at[c, "span_start"], ut.at[c, "span_end"]
        cs_w = cs_full[(cs_full >= a0) & (cs_full < a1)]
        for nb, dist in near[c].items():
            w0, w1 = max(a0, ut.at[nb, "span_start"]), min(a1, ut.at[nb, "span_end"])
            if w1 - w0 < 60:
                continue
            cs = cs_full[(cs_full >= w0) & (cs_full < w1)]
            nbt = ST.get(int(nb), np.empty(0))
            nbt = nbt[(nbt >= w0) & (nbt < w1)]
            if len(cs) < 20 or len(nbt) < 20:
                continue
            fr_nb = len(nbt) / (w1 - w0)
            hit = count_within(cs, nbt, -ZERO_HALF, ZERO_HALF) > 0
            p0 = hit.mean()
            p0_null = 1 - np.exp(-fr_nb * 2 * ZERO_HALF)
            post = count_within(cs, nbt, ZERO_HALF, POST_END).sum() / (len(cs) * fr_nb * (POST_END - ZERO_HALF))
            frac_nb_cs = (count_within(nbt, cs, -ZERO_HALF, ZERO_HALF) > 0).mean()
            r = dict(mouse=row.mouse, S=row.S, sorter=row.sorter, cs=int(c), nb=int(nb), dist_um=float(dist),
                     nb_group=ut.at[nb, "group"], nb_extragood=bool(ut.at[nb, "extragood"]),
                     nb_c4=ut.at[nb, "c4_label"] if pd.notna(ut.at[nb, "c4_label"]) else "unlabelled",
                     nb_role=ut.at[nb, "role"] if pd.notna(ut.at[nb, "role"]) else "",
                     nb_fr=fr_nb, cs_fr=len(cs) / (w1 - w0), n_cs=len(cs), p_zero=p0, p_zero_null=p0_null,
                     zero_excess=p0 - p0_null, frac_nb_cs=frac_nb_cs, post_ratio=post,
                     peak_lag_ms=peak_lag(cs, nbt), win_s=w1 - w0)
            pair_rows.append(r)
            if best is None or r["zero_excess"] > best["zero_excess"]:
                best = r
            # union over neighbours, on the CS's own span
            if len(cs_w) >= 20:
                nbu = nbt  # restricted to overlap; good enough for a union estimate
                h = count_within(cs_w, nbu, -ZERO_HALF, ZERO_HALF) > 0
                any_hit = h if any_hit is None else (any_hit | h)
                chance_rate += fr_nb
        cs_rows.append(dict(mouse=row.mouse, S=row.S, sorter=row.sorter, cs=int(c), y_um=ut.at[c, "y"],
                            cs_fr=ut.at[c, "fr"], n_near=len(near[c]),
                            union_p_zero=float(any_hit.mean()) if any_hit is not None else np.nan,
                            union_null=1 - np.exp(-chance_rate * 2 * ZERO_HALF),
                            best_nb=best["nb"] if best else np.nan,
                            best_excess=best["zero_excess"] if best else np.nan,
                            best_frac_nb_cs=best["frac_nb_cs"] if best else np.nan,
                            best_nb_fr=best["nb_fr"] if best else np.nan,
                            best_nb_c4=best["nb_c4"] if best else "",
                            best_nb_group=best["nb_group"] if best else "",
                            best_post_ratio=best["post_ratio"] if best else np.nan,
                            best_peak_lag=best["peak_lag_ms"] if best else np.nan))
    print(f"{row.mouse} {row.S} ({row.sorter}): {len(cs_ids)} CS done", flush=True)

pairs = pd.DataFrame(pair_rows)
csd = pd.DataFrame(cs_rows)
csd["union_excess"] = csd["union_p_zero"] - csd["union_null"]
csd["strong_shared"] = csd["best_excess"] > STRONG_EXCESS
pairs.to_csv(OUT / "cs_neighbour_pairs.csv", index=False)
csd.to_csv(OUT / "cs_zero_lag.csv", index=False)

# step0 v3 pairing outcome per CS
v3 = pd.read_csv(CMP / f"step0_{BATCH_TAG}" / "cs_units.csv")[["mouse", "S", "cs", "paired_any", "unpaired_reason"]]
csd = csd.merge(v3, on=["mouse", "S", "cs"], how="left")

pd.set_option("display.width", 250)
print("\n--- per session: CS whose spikes also appear in a nearby unit (zero lag ±0.5 ms, beyond chance)")
ss = csd.groupby(["mouse", "S", "sorter"]).agg(
    n_CS=("cs", "size"), best_excess_med=("best_excess", "median"),
    frac_strong=("strong_shared", "mean"), union_excess_med=("union_excess", "median"),
    frac_paired=("paired_any", "mean")).round(2)
sinks = (pairs[pairs["zero_excess"] > SINK_EXCESS].groupby(["mouse", "S", "nb"]).size()
         .rename("n_cs").reset_index())
ss["n_sink_units"] = sinks[sinks["n_cs"] >= 2].groupby(["mouse", "S"]).size().reindex(ss.index.droplevel(2)).fillna(0).values
print(ss.to_string())
ss.reset_index().to_csv(OUT / "session_zero_lag.csv", index=False)

print("\n--- pooled by sorter")
print(csd.groupby("sorter").agg(n=("cs", "size"), best_excess_med=("best_excess", "median"),
                                frac_strong=("strong_shared", "mean")).round(3).to_string())
print("\n--- does sharing go with failing to pair? (all CS)")
print(csd.groupby("paired_any")[["best_excess", "union_excess"]].median().round(3).to_string())
print(pd.crosstab(csd["strong_shared"], csd["paired_any"], margins=True).to_string())
g = ss.reset_index()
p = mannwhitneyu(g.loc[g.sorter == "ks4", "frac_strong"], g.loc[g.sorter == "ibl", "frac_strong"]).pvalue
print(f"session frac_strong KS4 vs IBL: MW p={p:.3f}")

for m, s in INSPECT:
    c = csd[(csd.mouse == m) & (csd.S == s)]
    pp = pairs[(pairs.mouse == m) & (pairs.S == s)]
    print("\n" + "=" * 80 + f"\n{m} {s}: {len(c)} CS")
    print(f"  best-neighbour zero-lag excess: median {c.best_excess.median():.2f}, "
          f"CS with >{STRONG_EXCESS:.0%} of spikes shared: {c.strong_shared.sum()}")
    print(f"  union over all neighbours: median p_zero {c.union_p_zero.median():.2f} vs chance {c.union_null.median():.2f}")
    print("  best neighbour of strongly shared CS — c4:", c.loc[c.strong_shared, "best_nb_c4"].value_counts().to_dict(),
          "| group:", c.loc[c.strong_shared, "best_nb_group"].value_counts().to_dict())
    print(f"  best neighbour post-CS ratio (0.5–8 ms; <1 pause, >1 excess): median "
          f"{c.loc[c.strong_shared, 'best_post_ratio'].median():.2f}; peak lag ms median "
          f"{c.loc[c.strong_shared, 'best_peak_lag'].median():.2f}")
    sk = pp[pp.zero_excess > SINK_EXCESS].groupby("nb").agg(
        n_cs=("cs", "size"), nb_fr=("nb_fr", "first"), nb_c4=("nb_c4", "first"), nb_group=("nb_group", "first"),
        nb_eg=("nb_extragood", "first"), frac_nb_cs=("frac_nb_cs", "sum"), excess_med=("zero_excess", "median"))
    print("  neighbour units carrying zero-lag spikes of >=2 CS (frac_nb_cs summed over those CS):")
    print(sk[sk.n_cs >= 2].sort_values("n_cs", ascending=False).round(2).head(15).to_string())

# figure: per-session distribution of best-neighbour zero-lag excess
fig, ax = plt.subplots(figsize=(11, 4.5))
order = csd.groupby(["mouse", "S"]).ngroup()
labs = csd.drop_duplicates(["mouse", "S"])[["mouse", "S", "sorter"]]
labs = labs.assign(k=range(len(labs)))
key = {(r.mouse, r.S): r.k for r in labs.itertuples()}
rng = np.random.default_rng(0)
for (m, s, so), g in csd.groupby(["mouse", "S", "sorter"]):
    k = key[(m, s)]
    xx = k + rng.uniform(-0.25, 0.25, len(g))
    ax.scatter(xx[~g.paired_any.fillna(False).values], g.best_excess[~g.paired_any.fillna(False).values], s=12,
               facecolor="none", edgecolor=SORTER_COLOR[so], lw=1)
    ax.scatter(xx[g.paired_any.fillna(False).values], g.best_excess[g.paired_any.fillna(False).values], s=14,
               color=SORTER_COLOR[so], edgecolor="white", lw=0.5)
ax.axhline(STRONG_EXCESS, color="#999999", lw=0.8, ls="--")
ax.set_xticks(labs.k)
ax.set_xticklabels([f"{r.mouse.split('_')[1][:6]} {r.S}" for r in labs.itertuples()], rotation=60, fontsize=7)
for t, r in zip(ax.get_xticklabels(), labs.itertuples()):
    t.set_color(SORTER_COLOR[r.sorter])
ax.set_ylabel("share of CS spikes also in best nearby unit\n(±0.5 ms, minus chance)")
ax.set_title("Do CS spikes reappear in a nearby unit?  filled = CS paired by step0 v3, open = unpaired;  "
             "tick colour: blue Kilosort4, orange IBL", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "figs" / "fig8_cs_zero_lag_sharing.png", dpi=150)
print("\noutputs in", OUT)
