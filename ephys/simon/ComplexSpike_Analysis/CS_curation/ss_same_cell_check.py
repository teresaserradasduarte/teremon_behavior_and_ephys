# -*- coding: utf-8 -*-
"""
ss_same_cell_check.py — are two SS units the same cell (e.g. soma + axon), i.e. detections of the same spikes?

Question (Simon, 2026-11-09): CS 734 pauses SS 925 (52 µm, 194 sp/s) and SS 938 (200 µm, 14 sp/s). Is 938 the same
Purkinje cell as 925, seen at another site (axon, reversed template)? Generalised to every pair of SS that are both
candidates of the same CS, in the 20 Hz (v7batch) or 10 Hz (v7batch_fr10) step0 run of this session.

Tests per SS pair (A = the faster unit, B = the other; inside the overlap of both units' stable spans):
  1. SS–SS cross-correlogram (B after A) vs both autocorrelograms, fine (0.05 ms) and coarse (0.5 ms).
       same spikes, two sites  → sharp peak at a fixed small lag (propagation delay), above chance
       same cell, other spikes → a refractory dip at 0 like A's own ACG (B fills in spikes A missed)
       two cells               → flat around 0 (only slow co-modulation, e.g. both paused by the same CS)
  2. Coincidences at the peak lag (± COINC): share of B spikes with an A spike, vs chance from A's rate.
  3. Cross-triggered raw signal (the decisive one): mean raw waveform on a channel strip covering both units for
       A spikes with / without a B spike, B spikes with / without an A spike, and random times.
     Per event, a "signature" score = least-squares amplitude of B's mean waveform on B's 3 largest channels at the
     expected time (A + lag, best within ± SHIFT), and the reverse on B events. 1 = full-size spike, 0 = none.
       same cell (soma + axon) → A spikes WITHOUT a detected B still carry a B-like deflection (score near 1)
       two cells               → scores of the "without" events ≈ random-time scores
Preprocessing = raw layer (int16 → µV, mean removed, common median reference over 384 channels, 300 Hz high-pass).
Read-only on Dropbox (sorter output + .meta are only read); the raw .bin is the LOCAL copy.
Output: CS_curation/<mouse>_<R#>/ss_same_cell/<A>_<B>/ (figure + summary.json) and
        CS_curation/<mouse>_<R#>/ss_same_cell/ss_same_cell.js (window.CS_SSPAIR, loaded by browser.html).
Nothing is merged; information only. Env SSPAIR_A + SSPAIR_B = run one pair only (JS not written).
"""
# %% 1. PARAMETERS + PATHS
from itertools import combinations
from pathlib import Path
import os as _os
import json
import os
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt

LOCAL_ROOT = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior")
MOUSE, SESSION_R = "1_CoteDor", "R5"
MOUSE, SESSION_R = (_os.environ["RAW_SESSION"].rsplit("_", 1) if _os.environ.get("RAW_SESSION") else (MOUSE, SESSION_R))  # batch: raw_all_sessions.py
SESS_DIR = LOCAL_ROOT / "CS_curation" / f"{MOUSE}_{SESSION_R}"
STEP0_RUNS = {"20 Hz": SESS_DIR / "step0_explore_bhvspan_v7batch", "10 Hz": SESS_DIR / "step0_explore_bhvspan_v7batch_fr10"}
RAW_BIN = (LOCAL_ROOT / "raw_copies" / "1_CoteDor" / "04082023_CoteDor_Cer_S5_g0_imec0"
           / "04082023_CoteDor_Cer_S5_g0_t0.imec0.ap.bin")
RAW_META = (Path(r"D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_data"
                 r"\20230801_ChocolateGroup\1_CoteDor\04082023_CoteDor_Cer_S5_g0\04082023_CoteDor_Cer_S5_g0_imec0")
            / "04082023_CoteDor_Cer_S5_g0_t0.imec0.ap.meta")             # read only
if _os.environ.get("RAW_BIN"):                       # batch: local copy + its local .meta (raw_all_sessions.py)
    RAW_BIN, RAW_META = Path(_os.environ["RAW_BIN"]), Path(_os.environ["RAW_META"])
OUT_ROOT = SESS_DIR / "ss_same_cell"

WIN = (-2e-3, 3e-3)       # s around the trigger spike
PAD = 7e-3
HP_HZ, HP_ORDER = 300.0, 3
COINC = 0.2e-3            # ± s around the CCG peak lag that counts as "the same event"
NO_PARTNER = 1.5e-3       # "without" = no partner spike within ± this
N_EV = 1500               # events per group
SIG_HALF = (0.5e-3, 1.0e-3)  # signature template: −0.5…+1 ms around the unit's own spike
SHIFT = 0.2e-3            # score = best fit within ± this around the expected time (same for every group)
HIST_BINS = np.linspace(-0.5, 1.8, 47)

# %% 2. LOAD
p0 = json.loads((STEP0_RUNS["20 Hz"] / "params_used.json").read_text())
SORTER_DIR = Path(p0["sorter_dir"])
meta = dict(l.strip().split("=", 1) for l in RAW_META.read_text().splitlines() if "=" in l)
FS = float(meta["imSampRate"]); N_CH_SAVED = int(meta["nSavedChans"])
_gain = float(meta["~imroTbl"].split(")(")[1].split()[3])
UV = float(meta["imAiRangeMax"]) / float(meta["imMaxInt"]) / _gain * 1e6
assert RAW_BIN.exists() and RAW_BIN.stat().st_size == int(meta["fileSizeBytes"]), "raw copy missing or incomplete"
raw = np.memmap(RAW_BIN, dtype=np.int16, mode="r").reshape(-1, N_CH_SAVED)
N_AP = 384
chan_pos = np.load(SORTER_DIR / "channel_positions.npy"); ych = chan_pos[:, 1]
chan_map = np.load(SORTER_DIR / "channel_map.npy").ravel()
spike_samp = np.asarray(np.load(SORTER_DIR / "spike_times.npy", mmap_mode="r")).ravel().astype(np.int64)
spike_clu = np.asarray(np.load(SORTER_DIR / "spike_clusters.npy", mmap_mode="r")).ravel()
ut = pd.read_csv(STEP0_RUNS["20 Hz"] / "tables" / "unit_table.csv").set_index("phyID")
SOS = butter(HP_ORDER, HP_HZ, "highpass", fs=FS, output="sos")
S0, S1, SP = int(round(WIN[0] * FS)), int(round(WIN[1] * FS)), int(round(PAD * FS))
OFF = np.arange(S0 - SP, S1 + SP); CROP = slice(SP, SP + (S1 - S0)); T_MS = np.arange(S0, S1) / FS * 1e3

# pairs: two SS both candidates of the same CS (either floor); which CS and which floor
PAIRS = {}
for floor, d in STEP0_RUNS.items():
    if not (d / "tables" / "cs_ss_pairs.csv").exists():
        continue
    c = pd.read_csv(d / "tables" / "cs_ss_pairs.csv"); c = c[c.candidate]
    for cs, ss in c.groupby("cs").ss:
        for a, b in combinations(sorted(set(ss)), 2):
            PAIRS.setdefault((int(a), int(b)), set()).add((int(cs), floor))
if os.environ.get("SSPAIR_A"):
    k = tuple(sorted((int(os.environ["SSPAIR_A"]), int(os.environ["SSPAIR_B"]))))
    PAIRS = {k: PAIRS.get(k, set())}
print(f"{len(PAIRS)} SS pairs: " + ", ".join(f"{a}–{b} (CS {sorted(PAIRS[(a, b)])})" for a, b in PAIRS))


def snippets(times_s, chs):
    s = np.round(np.asarray(times_s) * FS).astype(np.int64)
    s = s[(s + OFF[0] >= 0) & (s + OFF[-1] < raw.shape[0])]
    out = []
    for i in range(0, len(s), 128):
        x = raw[s[i:i + 128, None] + OFF[None, :], :N_AP].astype(np.float32)[:, :, chan_map]
        x -= x.mean(axis=1, keepdims=True)
        x -= np.median(x, axis=2, keepdims=True)
        out.append(sosfiltfilt(SOS, x[:, :, chs] * UV, axis=1)[:, CROP])
    return np.concatenate(out)


def ccg(ref, tgt, half, bw, self_=False):
    edges = np.arange(-half, half + bw / 2, bw)
    i0, i1 = np.searchsorted(tgt, ref - half), np.searchsorted(tgt, ref + half)
    d = np.concatenate([tgt[a:b] - r for r, a, b in zip(ref, i0, i1)]) if len(ref) else np.array([])
    if self_:
        d = d[d != 0]                                                           # drop each spike's own zero lag
    return edges[:-1] + bw / 2, np.histogram(d, edges)[0] / (len(ref) * bw)     # rate (sp/s) of tgt around ref


def has_near(src, tgt, h):
    i = np.searchsorted(tgt, src - h); ii = np.clip(i, 0, len(tgt) - 1)
    return (i < len(tgt)) & (tgt[ii] <= src + h)


def run_pair(u1, u2, cs_list):
    rng = np.random.default_rng(0)
    pick = lambda t, n: np.sort(rng.choice(t, n, replace=False)) if len(t) > n else t
    w = (max(ut.loc[u1, "span_start"], ut.loc[u2, "span_start"]), min(ut.loc[u1, "span_end"], ut.loc[u2, "span_end"]))
    tr = {u: (lambda t: t[(t >= w[0]) & (t < w[1])])(np.sort(spike_samp[spike_clu == u]) / FS) for u in (u1, u2)}
    A, B = sorted((u1, u2), key=lambda u: -len(tr[u]))                    # A = the faster unit
    ta, tb = tr[A], tr[B]; T = w[1] - w[0]
    # ---- CCG + coincidences
    xf, cf = ccg(ta, tb, 3e-3, 0.05e-3)
    xc, cc = ccg(ta, tb, 30e-3, 0.5e-3)
    _, aa = ccg(ta, ta, 30e-3, 0.5e-3, True); _, ab = ccg(tb, tb, 30e-3, 0.5e-3, True)
    far = np.abs(xc) >= 15e-3
    base = cc[far].mean()
    near = np.abs(xf) <= 1.5e-3
    lag = float(xf[near][np.argmax(cf[near])])
    peak_ratio = float(cf[near].max() / base)
    zero_ratio = float(cc[np.abs(xc) <= 1e-3].mean() / base)
    j = np.searchsorted(ta, tb - lag - COINC); jj = np.clip(j, 0, len(ta) - 1)
    b_with_a = (j < len(ta)) & (ta[jj] <= tb - lag + COINC)
    k = np.searchsorted(tb, ta + lag - COINC); kk = np.clip(k, 0, len(tb) - 1)
    a_with_b = (k < len(tb)) & (tb[kk] <= ta + lag + COINC)
    chance_b = 1 - np.exp(-len(ta) / T * 2 * COINC); chance_a = 1 - np.exp(-len(tb) / T * 2 * COINC)
    a_free = ~has_near(ta, tb, NO_PARTNER); b_free = ~has_near(tb, ta, NO_PARTNER)
    # ---- channel strip covering both units
    pa = np.ptp(snippets(pick(ta, 300), np.arange(len(chan_map))).mean(0), 0)
    pb = np.ptp(snippets(pick(tb, 300), np.arange(len(chan_map))).mean(0), 0)
    ya, yb = ych[np.argmax(pa)], ych[np.argmax(pb)]
    rows = np.unique(ych[(ych >= min(ya, yb) - 60) & (ych <= max(ya, yb) + 60)])
    if len(rows) > 24:                                                   # far apart: rows around each unit only
        rows = rows[(np.abs(rows - ya) <= 60) | (np.abs(rows - yb) <= 60)]
    strip = [int((lambda cand: cand[np.argmax(np.maximum(pa[cand] / pa.max(), pb[cand] / pb.max()))])(np.where(ych == yv)[0]))
             for yv in rows]
    ch_a = [int(c) for c in np.argsort(pa)[::-1][:3]]; ch_b = [int(c) for c in np.argsort(pb)[::-1][:3]]
    chs = sorted(set(strip) | set(ch_a) | set(ch_b), key=lambda c: (-ych[c], chan_pos[c, 0]))
    ia = [chs.index(c) for c in ch_a]; ib = [chs.index(c) for c in ch_b]
    # ---- cross-triggered snippets
    t_rand = np.sort(rng.uniform(w[0] + 0.01, w[1] - 0.01, N_EV))
    G = {"a_with": (ta[a_with_b], "A"), "a_without": (ta[a_free], "A"), "b_with": (tb[b_with_a], "B"),
         "b_without": (tb[b_free], "B"), "a_all": (ta, "A"), "b_all": (tb, "B"), "random": (t_rand, "R")}
    X = {g: snippets(pick(tt, N_EV), chs) if len(tt) else np.zeros((0, len(T_MS), len(chs))) for g, (tt, _) in G.items()}
    sel = (T_MS >= -SIG_HALF[0] * 1e3) & (T_MS <= SIG_HALF[1] * 1e3); idx = np.where(sel)[0]
    TA = X["a_all"].mean(0)[:, ia]; TB = X["b_all"].mean(0)[:, ib]; ns = int(round(SHIFT * FS)); i0 = int(np.argmin(np.abs(T_MS)))

    def score(x, ii, tmpl, centre_ms):
        if not len(x):
            return np.array([])
        c0 = int(np.argmin(np.abs(T_MS - centre_ms))) - i0
        tv = tmpl[sel]; nrm = (tv ** 2).sum(); best = np.full(len(x), -np.inf)
        for s in range(-ns, ns + 1):
            j2 = idx + c0 + s
            if j2.min() >= 0 and j2.max() < x.shape[1]:
                best = np.maximum(best, (x[:, j2][:, :, ii] * tv).sum((1, 2)) / nrm)
        return best
    L = lag * 1e3
    SC = {}
    for g, (_, kind) in G.items():
        if kind in ("A", "R"):
            SC[f"bsig_{g}"] = score(X[g], ib, TB, L)
        if kind in ("B", "R"):
            SC[f"asig_{g}"] = score(X[g], ia, TA, -L)
    med = {k_: (float(np.median(v)) if len(v) else None) for k_, v in SC.items()}
    # ---- verdict (information only)
    excess = b_with_a.mean() - chance_b
    if med["bsig_a_without"] is not None and med["bsig_a_without"] - med["bsig_random"] > 0.4:
        verdict = "same cell? B's waveform follows A spikes where B was not detected"
    elif excess > 0.2 and peak_ratio > 3:
        verdict = "same spikes detected twice (sharp CCG peak)"
    elif zero_ratio < 0.5:
        verdict = "possibly the same cell, complementary detection (refractory dip in the CCG)"
    else:
        verdict = "two different cells"
    # ---- figure
    od = OUT_ROOT / f"{A}_{B}"; od.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(17, 11)); gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1.15])
    ax = fig.add_subplot(gs[0, 0:2]); ax.bar(xf * 1e3, cf, width=0.05, color="k")
    ax.axhline(base, color="tab:orange", lw=1, ls="--"); ax.axvspan((lag - COINC) * 1e3, (lag + COINC) * 1e3, color="tab:red", alpha=.2)
    ax.set(xlabel=f"time of {B} spike after {A} spike (ms)", ylabel=f"{B} rate (sp/s)",
           title=f"CCG {A}→{B}, 0.05 ms bins: highest bin {L:+.2f} ms = {peak_ratio:.1f}× baseline (dashed)")
    ax = fig.add_subplot(gs[0, 2:4])
    ax.step(xc * 1e3, cc / base, where="mid", color="k", label=f"CCG {A}→{B}")
    ax.step(xc * 1e3, aa / aa[far].mean(), where="mid", color="tab:purple", alpha=.7, label=f"ACG {A}")
    ax.step(xc * 1e3, ab / ab[far].mean(), where="mid", color="tab:blue", alpha=.7, label=f"ACG {B}")
    ax.axhline(1, color="#999", lw=.8); ax.set(xlabel="lag (ms)", ylabel="rate / baseline", title="0.5 ms bins", ylim=(0, None)); ax.legend(fontsize=8)
    sc_ = max(np.ptp(TA, 0).max(), np.ptp(TB, 0).max())
    for i, (g, lab) in enumerate([("a_with", f"{A} with {B}"), ("a_without", f"{A} without {B}"),
                                  ("b_with", f"{B} with {A}"), ("b_without", f"{B} without {A}")]):
        ax = fig.add_subplot(gs[1, i]); m = X[g].mean(0) if len(X[g]) else np.zeros((len(T_MS), len(chs))); r = X["random"].mean(0)
        for c in range(len(chs)):
            col = "tab:purple" if c in ia else "tab:blue" if c in ib else "k"
            ax.plot(T_MS, r[:, c] / sc_ * 1.6 - c, color="#ccc", lw=.8); ax.plot(T_MS, m[:, c] / sc_ * 1.6 - c, color=col, lw=1)
        ax.axvline(0, color="#999", lw=.6, ls=":")
        ax.set_yticks(-np.arange(len(chs))); ax.set_yticklabels([f"{chs[c]} y{ych[chs[c]]:.0f}" for c in range(len(chs))], fontsize=6)
        ax.set_title(f"mean raw: {lab} (n={len(X[g])} of {len(G[g][0])})", fontsize=9); ax.set_xlabel("ms from trigger spike")
    for i, (lab, keys) in enumerate([(f"{B}'s waveform at {A}{L:+.2f} ms", ["bsig_a_with", "bsig_a_without", "bsig_a_all", "bsig_random"]),
                                     (f"{A}'s waveform at {B}{-L:+.2f} ms", ["asig_b_with", "asig_b_without", "asig_b_all", "asig_random"])]):
        ax = fig.add_subplot(gs[2, 2 * i:2 * i + 2])
        for k_, col in zip(keys, ["tab:green", "tab:red", "tab:gray", "k"]):
            if len(SC[k_]):
                ax.hist(SC[k_], HIST_BINS, histtype="step", density=True, color=col, lw=1.3, label=f"{k_.split('_', 1)[1]} (median {med[k_]:.2f})")
        ax.axvline(1, color="#999", ls=":"); ax.axvline(0, color="#999", ls=":")
        ax.set(xlabel="fitted amplitude (1 = the unit's own mean spike)", ylabel="density", title=lab); ax.legend(fontsize=8)
    fig.suptitle(f"{MOUSE} {SESSION_R}: same cell? SS {A} (A, {len(ta) / T:.0f} sp/s, ch {ch_a[0]}) vs SS {B} (B, {len(tb) / T:.0f} sp/s, "
                 f"ch {ch_b[0]}) — both candidates of CS {', '.join(str(c) for c, _ in sorted(cs_list))} — {w[0]:.0f}–{w[1]:.0f} s — {verdict}", fontsize=10)
    fig.tight_layout(); fig.savefig(od / f"same_cell_{A}_{B}.png", dpi=110); plt.close(fig)
    summ = dict(ss_a=A, ss_b=B, cs=[dict(cs=c, floor=f) for c, f in sorted(cs_list)], window_s=[float(w[0]), float(w[1])],
                n_a=int(len(ta)), n_b=int(len(tb)), fr_a=len(ta) / T, fr_b=len(tb) / T, peak_lag_ms=L, peak_ratio=peak_ratio,
                zero_1ms_ratio=zero_ratio, frac_b_with_a=float(b_with_a.mean()), chance_b_with_a=float(chance_b),
                frac_a_with_b=float(a_with_b.mean()), chance_a_with_b=float(chance_a), n_a_with=int(a_with_b.sum()),
                n_b_with=int(b_with_a.sum()), n_a_without=int(a_free.sum()), n_b_without=int(b_free.sum()),
                peak_ch_a=ch_a[0], peak_ch_b=ch_b[0], dist_um=float(np.hypot(*(chan_pos[ch_a[0]] - chan_pos[ch_b[0]]))),
                median_scores=med, verdict=verdict)
    (od / "summary.json").write_text(json.dumps(summ, indent=1))
    # ---- browser payload (µV rounded to 0.1)
    rnd = lambda a, n=1: [round(float(v), n) for v in a]
    js = dict(summ, ch=[int(c) for c in chs], ch_y=[float(ych[c]) for c in chs], ch_x=[float(chan_pos[c, 0]) for c in chs],
              ch_a=ch_a, ch_b=ch_b, fine_x=rnd(xf * 1e3, 3), fine=rnd(cf), coarse_x=rnd(xc * 1e3, 2), base=round(float(base), 2),
              coarse=rnd(cc / base, 3), acg_a=rnd(aa / aa[far].mean(), 3), acg_b=rnd(ab / ab[far].mean(), 3),
              means={g: [rnd(X[g].mean(0)[:, c]) for c in range(len(chs))] if len(X[g]) else None
                     for g in ("a_with", "a_without", "b_with", "b_without", "random")},
              n_read={g: int(len(X[g])) for g in X},
              hist={k_: np.histogram(v, HIST_BINS)[0].tolist() for k_, v in SC.items()})
    print(f"  {A}–{B}: lag {L:+.2f} ms ×{peak_ratio:.2f}, |lag|≤1 ms ×{zero_ratio:.2f}, B with A {b_with_a.mean():.1%} "
          f"(chance {chance_b:.1%}), B-sig without {med['bsig_a_without']} vs random {med['bsig_random']:.2f} → {verdict}")
    return js


# %% 3. RUN + SAVE
_t = time.time()
OUT = [run_pair(a, b, PAIRS[(a, b)]) for a, b in PAIRS]
if not os.environ.get("SSPAIR_A"):
    js = "window.CS_SSPAIR = " + json.dumps(dict(t_ms=[round(float(v), 4) for v in T_MS], hist_edges=HIST_BINS.round(3).tolist(),
                                                  coinc_ms=COINC * 1e3, no_partner_ms=NO_PARTNER * 1e3, shift_ms=SHIFT * 1e3,
                                                  sig_ms=[-SIG_HALF[0] * 1e3, SIG_HALF[1] * 1e3], n_ev=N_EV, hp_hz=HP_HZ,
                                                  created=time.strftime("%Y-%m-%d %H:%M"), pairs=OUT),
                                             separators=(",", ":")) + ";\n"
    OUT_ROOT.mkdir(parents=True, exist_ok=True)          # 2026-11-11: sessions with no SS pairs have no pair folder
    (OUT_ROOT / "ss_same_cell.js").write_text(js, encoding="utf-8")
    print(f"saved {OUT_ROOT / 'ss_same_cell.js'} ({len(js) / 1e6:.2f} MB)")
print(f"done ({time.time() - _t:.0f} s)")
