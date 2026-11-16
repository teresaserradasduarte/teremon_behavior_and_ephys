# -*- coding: utf-8 -*-
"""
sorter_compare_pass2_goodmua.py — second-pass SS-partner search for extra-good CS units, by sorter.

Pass 1 (sorter_compare_CS_SS.py) searched SS partners among extra-good CB non-CS units only.
Pass 2 adds tier-2 candidates: CB units whose Phy group is good or mua, that are NOT extra-good, and fire at
>= SS2_FR_MIN Hz in the window. Same definitions as step0_explore_session.py.
  - CS units: extra-good CB C4 PkC_cs (unchanged).
  - Tier 1: extra-good CB non-CS units (any C4 label, incl. unlabelled).
  - Tier 2: good/mua CB, not extra-good, FR >= 20 Hz in win_selected (no C4 label: C4 ran on extra-good only).
  - Candidate partner: pause_ratio < 0.7, asym > 0.2, no collision dip, jitter p < 0.05.
    The jitter test is run on the N_JITTER_TOP lowest-pause-ratio candidates of each CS.

Inputs: pass-1 outputs (session_table.csv, units_all.csv: CS ids, window in imec s, CS peak positions) and the
sorter folders (read-only). Window = bhv.win_selected mapped to imec (taken from units_all.csv).
Outputs: <OUT_DIR>/sessions/<mouse>_<S>/cs_ss_pairs_pass2.csv, pass2_cs_units.csv, pass2_session_summary.csv,
         pass2_stats.csv, figs/fig5_pass2_pairing.png

GOLDEN RULE: Dropbox is READ-ONLY.
"""

# %% 1. PARAMETERS
from pathlib import Path
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

DROPBOX_ROOT = Path(r"D:\Learning Lab Dropbox")
CODE_DIR = Path(r"C:\Users\simon\Desktop\Code\EPhys\ComplexSpike_Analysis")
OUT_DIR = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior\CS_curation\_sorter_comparison")
FORCE = False
ONLY = None                       # e.g. ("4_Milka", "S4") to test one session
import os
if os.environ.get("PASS2_ONLY"):  # e.g. PASS2_ONLY=4_Milka:S4 (process one session, then stop before stats)
    ONLY = tuple(os.environ["PASS2_ONLY"].split(":"))

MICE = ["1_CoteDor", "2_Lindt", "3_Toblerone", "4_Milka", "5_FerreroRocher"]
REGION, CS_LABEL = "CB", "PkC_cs"
SS2_GROUPS, SS2_FR_MIN = ("good", "mua"), 20.0
MIN_CS_SPIKES, MIN_SS_SPIKES = 20, 100
EDGE_PAD = 0.1
COARSE_BIN, COARSE_WIN = 1e-3, 50e-3
FINE_BIN, FINE_WIN = 0.1e-3, 10e-3
PAUSE_WIN, MIRROR_WIN, BASELINE_WIN = (0.0, 8e-3), (-8e-3, 0.0), (-50e-3, -10e-3)
FINE_BASELINE_WIN, DIP_HALF = (-10e-3, -2e-3), 0.5e-3
JITTER_WIN, JITTER_HALF, N_SURR, RNG_SEED = (1e-3, 8e-3), 25e-3, 500, 0
CAND_PAUSE_RATIO_MAX, CAND_JITTER_P_MAX, CAND_ASYM_MIN = 0.7, 0.05, 0.2
N_JITTER_TOP = 5
N_PERM = 10000
SORTER_COLOR = {"ks4": "#2a78d6", "ibl": "#eb6834"}
SORTER_NAME = {"ks4": "Kilosort4", "ibl": "IBL sorter"}
MOUSE_MARKER = dict(zip(MICE, ["o", "s", "^", "D", "v"]))


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


OUT_DIR = safe_out_dir(OUT_DIR)
FIG_DIR = safe_out_dir(OUT_DIR / "figs")
sessions = pd.read_csv(OUT_DIR / "session_table.csv")
units1 = pd.read_csv(OUT_DIR / "units_all.csv")


# %% 2. FUNCTIONS (same definitions as step0 / pass 1)
def _pair_diffs(ref, tgt, win, chunk=200_000):
    out = []
    for s in range(0, len(ref), chunk):
        r = ref[s:s + chunk]
        lo = np.searchsorted(tgt, r - win, "left")
        hi = np.searchsorted(tgt, r + win, "left")
        n = hi - lo
        tot = int(n.sum())
        if tot == 0:
            continue
        ri = np.repeat(np.arange(len(r)), n)
        j = lo[ri] + (np.arange(tot) - np.repeat(np.cumsum(n) - n, n))
        out.append(tgt[j] - r[ri])
    return np.concatenate(out) if out else np.empty(0)


def ccg(ref, tgt, bin_s, win_s):
    nb = int(round(win_s / bin_s))
    d = _pair_diffs(ref, tgt, nb * bin_s)
    idx = np.floor(d / bin_s + 1e-9).astype(np.int64) + nb
    idx = idx[(idx >= 0) & (idx < 2 * nb)]
    counts = np.bincount(idx, minlength=2 * nb).astype(float)
    return counts / (max(len(ref), 1) * bin_s), (np.arange(2 * nb) - nb) * bin_s


def _wm(left, bin_s, win):
    return (left >= win[0] - 1e-12) & (left + bin_s <= win[1] + 1e-12)


def pause_metrics(cs, ss):
    rc, lc = ccg(cs, ss, COARSE_BIN, COARSE_WIN)
    rf, lf = ccg(cs, ss, FINE_BIN, FINE_WIN)
    base = rc[_wm(lc, COARSE_BIN, BASELINE_WIN)].mean()
    m = dict(baseline_rate=base)
    m["pause_ratio"] = rc[_wm(lc, COARSE_BIN, PAUSE_WIN)].mean() / base if base > 0 else np.nan
    m["mirror_ratio"] = rc[_wm(lc, COARSE_BIN, MIRROR_WIN)].mean() / base if base > 0 else np.nan
    m["asym"] = m["mirror_ratio"] - m["pause_ratio"]
    pos = lc >= -1e-12
    ok = np.nonzero(rc[pos] >= 0.5 * base)[0]
    m["pause_dur50_ms"] = lc[pos][ok[0]] * 1e3 if len(ok) else np.nan
    fb = rf[_wm(lf, FINE_BIN, FINE_BASELINE_WIN)].mean()
    dn = rf[_wm(lf, FINE_BIN, (-DIP_HALF, 0.0))].mean()
    dp = rf[_wm(lf, FINE_BIN, (0.0, DIP_HALF))].mean()
    late = rc[_wm(lc, COARSE_BIN, (1e-3, PAUSE_WIN[1]))].mean()
    m["collision_dip"] = bool(fb > 0 and dn < 0.5 * fb and dp < 0.5 * fb and late >= 0.8 * base)
    m["zero_peak_ratio"] = rf[_wm(lf, FINE_BIN, (-0.5e-3, 0.5e-3))].mean() / fb if fb > 0 else np.nan
    return m


def jitter_pause_test(cs, ss):
    rng = np.random.default_rng(RNG_SEED)
    obs = int((np.searchsorted(ss, cs + JITTER_WIN[1]) - np.searchsorted(ss, cs + JITTER_WIN[0])).sum())
    flat = (cs[None, :] + rng.uniform(-JITTER_HALF, JITTER_HALF, (N_SURR, len(cs)))).ravel()
    c = np.searchsorted(ss, flat + JITTER_WIN[1]) - np.searchsorted(ss, flat + JITTER_WIN[0])
    surr = c.reshape(N_SURR, len(cs)).sum(1)
    return dict(jitter_p=(1 + np.sum(surr <= obs)) / (1 + N_SURR),
                jitter_deficit=1 - obs / surr.mean() if surr.mean() > 0 else np.nan)


# %% 3. PER SESSION
def process_session(row):
    sd = Path(row.sorter_dir)
    files = {p.name.lower(): p for p in sd.iterdir() if p.is_file()}
    u1 = units1[(units1["mouse"] == row.mouse) & (units1["S"] == row.S)].set_index("phyID")
    T0, T1 = float(u1["win_start_imec"].iloc[0]), float(u1["win_end_imec"].iloc[0])
    SR = float(u1["sample_rate"].iloc[0])

    grp = pd.read_csv(files["cluster_group.tsv"], sep="\t")
    grp = grp.set_index(grp.columns[0]).iloc[:, 0].rename("group")
    reg = pd.read_csv(files["cluster_channel_region_map.txt"], sep="\t").set_index("cluster_id")["brain_region"]
    info = pd.read_csv(files["cluster_info.tsv"], sep="\t").set_index("cluster_id")
    chan_pos = np.load(files["channel_positions.npy"])

    cs_ids = u1.index[u1["c4_label"] == CS_LABEL].tolist()
    tier1 = u1.index[u1["c4_label"] != CS_LABEL].tolist()
    eg_all = set(u1.index)  # extra-good CB units
    cand2 = [c for c in grp.index if grp[c] in SS2_GROUPS and reg.get(c) == REGION and c not in eg_all]

    st_all = np.asarray(np.load(files["spike_times.npy"], mmap_mode="r"), dtype=np.float64).ravel() / SR
    clu = np.asarray(np.load(files["spike_clusters.npy"], mmap_mode="r")).ravel().astype(np.int64)
    want = np.array(cs_ids + tier1 + cand2, dtype=np.int64)
    sel = np.isin(clu, want) & (st_all >= T0 - EDGE_PAD) & (st_all < T1 + EDGE_PAD)
    t_s, c_s = st_all[sel], clu[sel]
    del st_all, clu
    order = np.argsort(c_s, kind="stable")
    t_s, c_s = t_s[order], c_s[order]
    ids, starts, counts = np.unique(c_s, return_index=True, return_counts=True)
    ST_PAD = {int(c): t_s[s:s + n] for c, s, n in zip(ids, starts, counts)}
    ST = {c: t[(t >= T0) & (t < T1)] for c, t in ST_PAD.items()}

    T = T1 - T0
    tier2 = [c for c in cand2 if len(ST.get(c, [])) / T >= SS2_FR_MIN]

    def pos(c):
        if c in u1.index and pd.notna(u1.at[c, "x_um"]):
            return u1.at[c, "x_um"], u1.at[c, "depth_um"]
        ch = int(info.at[c, "ch"]) if c in info.index else None
        return (chan_pos[ch, 0], chan_pos[ch, 1]) if ch is not None and ch < len(chan_pos) else (np.nan, np.nan)

    P = {c: pos(c) for c in cs_ids + tier1 + tier2}
    rows = []
    for cs in cs_ids:
        if len(ST.get(cs, [])) < MIN_CS_SPIKES:
            continue
        cand = []
        for ss, tier in [(s, 1) for s in tier1] + [(s, 2) for s in tier2]:
            if len(ST_PAD.get(ss, [])) < MIN_SS_SPIKES:
                continue
            m = pause_metrics(ST[cs], ST_PAD[ss])
            lab = u1.at[ss, "c4_label"] if tier == 1 else "none(tier2)"
            m.update(cs=cs, ss=ss, tier=tier, ss_group=grp.get(ss, ""), ss_c4=lab, ss_fr=len(ST[ss]) / T,
                     n_cs=len(ST[cs]),
                     dist_um=float(np.hypot(P[cs][0] - P[ss][0], P[cs][1] - P[ss][1])))
            cand.append(m)
        cand.sort(key=lambda d: (np.nan_to_num(d["pause_ratio"], nan=9), -d["asym"]))
        for k, m in enumerate(cand):
            m["rank_pr"] = k + 1
            if k < N_JITTER_TOP:
                m.update(jitter_pause_test(ST[cs], ST_PAD[m["ss"]]))
        rows.extend(cand)
    pairs = pd.DataFrame(rows)
    if len(pairs):
        if "jitter_p" not in pairs:
            pairs["jitter_p"] = np.nan
        pairs["candidate"] = ((pairs["pause_ratio"] < CAND_PAUSE_RATIO_MAX) & (pairs["asym"] > CAND_ASYM_MIN)
                              & ~pairs["collision_dip"] & (pairs["jitter_p"].fillna(1.0) < CAND_JITTER_P_MAX))
    meta = dict(n_tier1=len(tier1), n_cand2_groupmatch=len(cand2), n_tier2=len(tier2))
    return pairs, meta


SES_DIR = safe_out_dir(OUT_DIR / "sessions")
metas = []
for row in sessions.itertuples():
    if ONLY and (row.mouse, row.S) != tuple(ONLY):
        continue
    d = safe_out_dir(SES_DIR / f"{row.mouse}_{row.S}")
    f = d / "cs_ss_pairs_pass2.csv"
    fm = d / "pass2_meta.csv"
    if f.exists() and fm.exists() and not FORCE:
        metas.append(pd.read_csv(fm).iloc[0].to_dict())
        print(f"cached  {row.mouse} {row.S}")
        continue
    t0 = time.time()
    try:
        pairs, meta = process_session(row)
    except Exception as e:
        print(f"FAILED  {row.mouse} {row.S}: {type(e).__name__}: {e}")
        continue
    meta.update(mouse=row.mouse, S=row.S, sorter=row.sorter)
    pairs.to_csv(f, index=False)
    pd.DataFrame([meta]).to_csv(fm, index=False)
    metas.append(meta)
    nc = pairs.loc[pairs["candidate"], "cs"].nunique() if len(pairs) else 0
    ncs = pairs["cs"].nunique() if len(pairs) else 0
    print(f"done    {row.mouse} {row.S} ({row.sorter}) {time.time() - t0:.0f} s | tier1 {meta['n_tier1']} | "
          f"tier2 {meta['n_tier2']} (of {meta['n_cand2_groupmatch']} good/mua CB) | CS paired {nc}/{ncs}")
meta = pd.DataFrame(metas)
if ONLY:
    raise SystemExit("ONLY set: single-session test done, summary skipped")


# %% 4. PER-CS AND PER-SESSION SUMMARY
cs_rows = []
for r in meta.itertuples():
    f = SES_DIR / f"{r.mouse}_{r.S}" / "cs_ss_pairs_pass2.csv"
    try:
        p = pd.read_csv(f)
    except pd.errors.EmptyDataError:
        continue
    for cs, g in p.groupby("cs"):
        c = g[g["candidate"]].sort_values("pause_ratio")
        c1 = c[c["tier"] == 1]
        b = c.iloc[0] if len(c) else None
        cs_rows.append(dict(mouse=r.mouse, S=r.S, sorter=r.sorter, cs=cs, n_cs=g["n_cs"].iloc[0],
                            paired_tier1=len(c1) > 0, paired_any=len(c) > 0,
                            paired_only_tier2=len(c) > 0 and len(c1) == 0,
                            n_candidates=len(c), best_ss=b["ss"] if b is not None else np.nan,
                            best_tier=b["tier"] if b is not None else np.nan,
                            best_group=b["ss_group"] if b is not None else "",
                            best_c4=b["ss_c4"] if b is not None else "",
                            best_pause_ratio=b["pause_ratio"] if b is not None else np.nan,
                            best_pause_dur50_ms=b["pause_dur50_ms"] if b is not None else np.nan,
                            best_jitter_deficit=b["jitter_deficit"] if b is not None else np.nan,
                            best_dist_um=b["dist_um"] if b is not None else np.nan,
                            best_ss_fr=b["ss_fr"] if b is not None else np.nan))
csu = pd.DataFrame(cs_rows)
# one-to-one: how many CS share their best SS with another CS in the same session
csu["best_ss_shared"] = csu.duplicated(["mouse", "S", "best_ss"], keep=False) & csu["paired_any"]
csu.to_csv(OUT_DIR / "pass2_cs_units.csv", index=False)

ss_ = (csu.groupby(["mouse", "S", "sorter"])
       .agg(n_CS=("cs", "size"), frac_paired_tier1=("paired_tier1", "mean"), frac_paired_any=("paired_any", "mean"),
            n_paired_any=("paired_any", "sum"), n_best_shared=("best_ss_shared", "sum"),
            best_pause_ratio_med=("best_pause_ratio", "median"), best_pause_dur_med=("best_pause_dur50_ms", "median"))
       .reset_index())
ss_ = meta.merge(ss_, on=["mouse", "S", "sorter"], how="left")
ss_.to_csv(OUT_DIR / "pass2_session_summary.csv", index=False)
with pd.option_context("display.width", 250):
    print("--- per session")
    print(ss_[["mouse", "S", "sorter", "n_tier1", "n_cand2_groupmatch", "n_tier2", "n_CS", "frac_paired_tier1",
               "frac_paired_any", "n_best_shared", "best_pause_ratio_med", "best_pause_dur_med"]].round(2).to_string(index=False))
    print("--- pooled CS units per sorter")
    print(csu.groupby("sorter").agg(n=("cs", "size"), paired_tier1=("paired_tier1", "mean"),
                                    paired_any=("paired_any", "mean"), only_tier2=("paired_only_tier2", "mean"),
                                    shared_best=("best_ss_shared", "mean")).round(3))
    pa = csu[csu["paired_any"]]
    print("--- best partner of paired CS: tier / group / C4")
    print(pd.crosstab(pa["sorter"], [pa["best_tier"], pa["best_group"]]))
    print(pd.crosstab(pa["sorter"], pa["best_c4"]))
    print("--- best-partner pause properties (paired CS)")
    print(pa.groupby("sorter")[["best_pause_ratio", "best_pause_dur50_ms", "best_jitter_deficit", "best_dist_um",
                                "best_ss_fr"]].median().round(2))
    print("--- tier-2 candidate pool per sorter (sum over sessions)")
    print(meta.groupby("sorter")[["n_tier1", "n_cand2_groupmatch", "n_tier2"]].sum())


# %% 5. STATS (session level)
def within_mouse_perm(df, col, n_perm=N_PERM, seed=0):
    d = df[["mouse", "sorter", col]].dropna()
    groups = [g for _, g in d.groupby("mouse") if g["sorter"].nunique() == 2]
    labs = [g["sorter"].values for g in groups]

    def stat(ls):
        return np.mean([np.median(g[col].values[l == "ks4"]) - np.median(g[col].values[l == "ibl"])
                        for g, l in zip(groups, ls)])

    obs = stat(labs)
    rng = np.random.default_rng(seed)
    null = np.array([stat([rng.permutation(l) for l in labs]) for _ in range(n_perm)])
    return obs, (1 + np.sum(np.abs(null) >= abs(obs) - 1e-12)) / (1 + n_perm), len(groups)


ss_["gain_any_minus_tier1"] = ss_["frac_paired_any"] - ss_["frac_paired_tier1"]
ss_["tier2_per_eg"] = ss_["n_tier2"] / ss_["n_tier1"].replace(0, np.nan)
st = []
for col in ["n_tier2", "tier2_per_eg", "frac_paired_tier1", "frac_paired_any", "gain_any_minus_tier1",
            "best_pause_ratio_med", "best_pause_dur_med"]:
    k, i = ss_.loc[ss_["sorter"] == "ks4", col].dropna(), ss_.loc[ss_["sorter"] == "ibl", col].dropna()
    if len(k) < 2 or len(i) < 2:
        continue
    d, p, nm = within_mouse_perm(ss_, col)
    st.append(dict(metric=col, ks4_median=k.median(), ibl_median=i.median(), n_ks4=len(k), n_ibl=len(i),
                   p_mannwhitney=mannwhitneyu(k, i).pvalue, within_mouse_diff=d, p_within_mouse=p, n_mice=nm))
st = pd.DataFrame(st)
st.to_csv(OUT_DIR / "pass2_stats.csv", index=False)
print("--- session-level stats (uncorrected, exploratory)")
print(st.round(4).to_string(index=False))


# %% 6. FIGURE
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "font.size": 9, "savefig.dpi": 150})
rng = np.random.default_rng(0)
fig, axs = plt.subplots(1, 4, figsize=(14, 4))
for ax, (col, t) in zip(axs[:3], [("frac_paired_tier1", "CS paired: extra-good SS only"),
                                  ("frac_paired_any", "CS paired: + good/mua (tier 2)"),
                                  ("n_tier2", "tier-2 SS candidates (good/mua, ≥20 Hz)")]):
    for xi, s in enumerate(["ks4", "ibl"]):
        sub = ss_[ss_["sorter"] == s].dropna(subset=[col])
        for mouse, g in sub.groupby("mouse"):
            ax.scatter(xi + rng.uniform(-0.18, 0.18, len(g)), g[col], s=28, marker=MOUSE_MARKER[mouse],
                       facecolor=SORTER_COLOR[s], edgecolor="white", linewidth=0.6, zorder=3)
        ax.plot([xi - 0.3, xi + 0.3], [sub[col].median()] * 2, color="#333333", lw=2, zorder=4)
    r = st.set_index("metric")
    ax.set_title(f"{t}\nMW p={r.at[col, 'p_mannwhitney']:.3f} | within-mouse p={r.at[col, 'p_within_mouse']:.3f}",
                 fontsize=8)
    ax.set_xticks([0, 1], [SORTER_NAME["ks4"], SORTER_NAME["ibl"]])
    ax.set_xlim(-0.6, 1.6)
ax = axs[3]
pa = csu[csu["paired_any"]]
for s in ("ks4", "ibl"):
    v = pa.loc[pa["sorter"] == s, "best_pause_dur50_ms"].dropna()
    ax.hist(v, bins=np.arange(0, 51, 2), histtype="step", lw=2, color=SORTER_COLOR[s], density=True,
            label=f"{SORTER_NAME[s]} (n={len(v)})")
ax.set_xlabel("pause duration to 50% recovery (ms)")
ax.set_ylabel("density")
ax.set_title("best partner pause duration (paired CS)", fontsize=8)
ax.legend(frameon=False, fontsize=8)
h = [plt.Line2D([], [], marker=MOUSE_MARKER[m], ls="", color="#666666", label=m) for m in MICE]
fig.legend(handles=h, loc="lower center", ncol=5, frameon=False)
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(FIG_DIR / "fig5_pass2_pairing.png")
plt.close(fig)
print(f"saved {FIG_DIR / 'fig5_pass2_pairing.png'}")
