# -*- coding: utf-8 -*-
"""
pause_figure_data.py — data for the browser's "Pause figure" and "Ephaptic check" cards (spike times only, no raw file).

Inspired by Han, Chen et al. 2020 (Nat Neurosci 23:1399, Fig. 1 and 8c): a CS pauses its own PC's SS for ~10+ ms, and
ALSO briefly (~0–2 ms) suppresses the SS of neighbouring PCs (ephaptic, strong ≤ 25 µm, gone > 75 µm).

Per session (pairs of both step0 floors, v7batch = 20 Hz and v7batch_fr10 = 10 Hz; each pair uses its own pair window
= overlap of both stable spans, as in step0; baseline = SS rate −50…−10 ms before the CS, step0's BASELINE_WIN):
  * per exported pair (candidate in either floor | step0 rank ≤ 8 | #1/#2 among this CS's SS on pause, asym, jitter z or
    distance | combined rank ≤ 8 | the CS's control SS):
      raster   ≤ 300 random CS events (in time order), SS spikes −10…+25 ms (units of 0.1 ms)
      psth     −10…+25 ms, 0.5 ms bins, ÷ baseline, with a 95 % band (bootstrap over CS events, 200×)
      zoom     −2…+5 ms, 0.2 ms bins, ÷ baseline, with band (relative timing of CS and pause onset, their Fig. 1f–g)
      quarters the pair window in 4 equal parts, psth of each ÷ that part's own baseline (stability, like ba–bd)
  * per CS: control SS = the nearest SS that is not a candidate of this CS (either floor) and not CS-flagged
  * ephaptic, every CS×SS pair: r_0_2 / r_3_8 = SS rate 0–2 / 3–8 ms ÷ baseline, distance, partner (= candidate of
    this CS in either floor), flagged (cs_share_flag or cs_event_flag: the "SS" carries CS spikes or CS pieces → its
    0–2 ms is not a suppression measure); per CS and for the session: mean ± SEM normalised curves (−5…+20 ms) of the
    non-partner, non-flagged SS in the distance groups ≤25, 25–50, 50–75, >75 µm (their Fig. 1d–e).

Output (local only): CS_curation/<session>/pause_figure/pause_figure.js (window.CS_PAUSEFIG, loaded by browser.html)
                     CS_curation/<session>/pause_figure/ephaptic_pairs.csv
Dropbox: only the sorter's spike_times.npy / spike_clusters.npy / params.py are READ.
Usage: python pause_figure_data.py [session ...]     (default: every session with a step0 v7batch run)
"""
from pathlib import Path
import json
import re
import sys
import time

import numpy as np
import pandas as pd

LOCAL_ROOT = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior")
CUR = LOCAL_ROOT / "CS_curation"
RUNS = {"20 Hz": "step0_explore_bhvspan_v7batch", "10 Hz": "step0_explore_bhvspan_v7batch_fr10"}
DROPBOX = "D:\\Learning Lab Dropbox\\"

BIN, ZBIN = 0.5e-3, 0.2e-3
PSTH_WIN, ZOOM_WIN, GROUP_WIN = (-10e-3, 25e-3), (-2e-3, 5e-3), (-5e-3, 20e-3)
EARLY_WIN, PAUSE_WIN = (0.0, 2e-3), (3e-3, 8e-3)
N_RASTER, N_BOOT, TOP_N = 300, 200, 8
DIST_GROUPS = [(0, 25, "≤25 µm"), (25, 50, "25–50 µm"), (50, 75, "50–75 µm"), (75, 1e9, ">75 µm")]
RNG = np.random.default_rng(0)


def per_event_counts(cs, ss, edges):
    """(n_events, n_bins) SS spike counts in [cs + edges[k], cs + edges[k+1])."""
    idx = np.searchsorted(ss, cs[:, None] + edges[None, :], "left")
    return np.diff(idx, axis=1)


def load_session(sess):
    p0 = json.loads((CUR / sess / RUNS["20 Hz"] / "params_used.json").read_text())
    sd = Path(p0["sorter_dir"])
    sr = float(re.search(r"sample_rate\s*=\s*([0-9.eE+]+)", (sd / "params.py").read_text()).group(1))
    t = np.asarray(np.load(sd / "spike_times.npy", mmap_mode="r")).ravel().astype(np.float64) / sr
    clu = np.asarray(np.load(sd / "spike_clusters.npy", mmap_mode="r")).ravel().astype(np.int64)
    o = np.argsort(clu, kind="stable"); cu, start = np.unique(clu[o], return_index=True)
    stop = np.r_[start[1:], len(o)]
    ST = {int(c): np.sort(t[o[a:b]]) for c, a, b in zip(cu, start, stop)}
    return p0, ST


def load_pairs(sess):
    """Union of both floors' pair tables; the row of the 20 Hz run wins when a pair is in both (same window)."""
    out = []
    for floor, run in RUNS.items():
        f = CUR / sess / run / "tables" / "cs_ss_pairs.csv"
        if f.exists():
            out.append(pd.read_csv(f).assign(floor=floor))
    p = pd.concat(out, ignore_index=True)
    cand_any = p[p.candidate.astype(bool)].groupby(["cs", "ss"]).size().index
    p = p.drop_duplicates(["cs", "ss"], keep="first").reset_index(drop=True)
    p["partner"] = [k in cand_any for k in zip(p.cs, p.ss)]
    p["flagged"] = p.get("cs_share_flag", False).fillna(False).astype(bool) | p.get("cs_event_flag", False).fillna(False).astype(bool)
    return p


def export_set(p):
    """Pairs to export in full: same rule as the browser's top N + best-per-criterion rows (2026-11-16)."""
    keep = set()
    for cs, g in p.groupby("cs"):
        r = {}
        for k, asc in (("pause_ratio", True), ("asym", False), ("jitter_z", True), ("dist_um", True)):
            r[k] = g[k].rank(method="min", ascending=asc, na_option="bottom")
            keep |= {(cs, s) for s in g.ss[r[k] <= 2]}
        comb = pd.concat(r, axis=1).mean(axis=1)
        keep |= {(cs, s) for s in g.ss[comb.rank(method="first") <= TOP_N]}
        keep |= {(cs, s) for s in g.ss[(g["rank"] <= TOP_N) | g.partner]}
    return keep


def pair_curves(cs_t, ss_t, base_win, boot):
    """Normalised psth / zoom (+ bands) and the per-event matrices needed for quarters."""
    e_main = np.arange(base_win[0], PSTH_WIN[1] + BIN / 2, BIN)          # −50…+25 ms
    C = per_event_counts(cs_t, ss_t, e_main)
    left = e_main[:-1]
    bmask = (left >= base_win[0] - 1e-12) & (left + BIN <= base_win[1] + 1e-12)
    pmask = left >= PSTH_WIN[0] - 1e-12
    n = len(cs_t)
    base = C[:, bmask].sum() / (n * bmask.sum() * BIN)
    out = dict(n_cs=int(n), base_hz=round(float(base), 2))
    if base <= 0 or n == 0:
        return out, C, left, bmask
    out["psth"] = np.round(C[:, pmask].mean(0) / BIN / base, 3).tolist()
    e_z = np.arange(ZOOM_WIN[0], ZOOM_WIN[1] + ZBIN / 2, ZBIN)
    Z = per_event_counts(cs_t, ss_t, e_z)
    out["zoom"] = np.round(Z.mean(0) / ZBIN / base, 3).tolist()
    if boot:
        W = RNG.multinomial(n, np.full(n, 1 / n), size=N_BOOT).astype(np.float32) / n      # (boot, n) weights
        bb = (W @ C[:, bmask]).sum(1) / (bmask.sum() * BIN)
        bb[bb <= 0] = np.nan
        for key, M, w in (("psth", C[:, pmask], BIN), ("zoom", Z, ZBIN)):
            r = (W @ M) / w / bb[:, None]
            lo, hi = np.nanpercentile(r, [2.5, 97.5], axis=0)
            out[key + "_lo"], out[key + "_hi"] = np.round(lo, 3).tolist(), np.round(hi, 3).tolist()
    return out, C, left, bmask


def run(sess):
    t0 = time.time()
    p0, ST = load_session(sess)
    base_win = tuple(p0.get("BASELINE_WIN", [-0.05, -0.01]))
    P = load_pairs(sess)
    keep = export_set(P)
    # control SS per CS: nearest SS that is neither a partner nor flagged
    control = {}
    for cs, g in P.groupby("cs"):
        g = g[~g.partner & ~g.flagged].sort_values("dist_um")
        if len(g):
            control[int(cs)] = int(g.ss.iloc[0]); keep.add((cs, g.ss.iloc[0]))
    e_main = np.arange(base_win[0], PSTH_WIN[1] + BIN / 2, BIN)
    left = e_main[:-1]
    gmask = (left >= GROUP_WIN[0] - 1e-12) & (left + BIN <= GROUP_WIN[1] + 1e-12)
    emask = (left >= EARLY_WIN[0] - 1e-12) & (left + BIN <= EARLY_WIN[1] + 1e-12)
    qmask = (left >= PAUSE_WIN[0] - 1e-12) & (left + BIN <= PAUSE_WIN[1] + 1e-12)
    pmask = left >= PSTH_WIN[0] - 1e-12
    pairs_out, eph, curves = {}, [], []
    for r in P.itertuples():
        cs, ss = int(r.cs), int(r.ss)
        if cs not in ST or ss not in ST:
            continue
        w = (float(r.win_start), float(r.win_end))
        a, b = ST[cs], ST[ss]
        cs_t = a[(a >= w[0] - base_win[0]) & (a < w[1] - PSTH_WIN[1])]          # full −50…+25 ms inside the window
        ss_t = b[(b >= w[0]) & (b < w[1])]
        full = (cs, ss) in keep
        o, C, _, bmask = pair_curves(cs_t, ss_t, base_win, boot=full)
        if "psth" not in o:
            continue
        rate = C.mean(0) / BIN / (o["base_hz"])
        eph.append(dict(cs=cs, ss=ss, dist_um=r.dist_um, r_0_2=float(rate[emask].mean()), r_3_8=float(rate[qmask].mean()),
                        partner=bool(r.partner), flagged=bool(r.flagged), floor=r.floor, n_cs=o["n_cs"]))
        if not r.partner and not r.flagged:
            curves.append((cs, r.dist_um, rate[gmask]))
        if not full:
            continue
        # raster
        k = np.sort(RNG.choice(len(cs_t), min(N_RASTER, len(cs_t)), replace=False))
        i0 = np.searchsorted(ss_t, cs_t[k] + PSTH_WIN[0]); i1 = np.searchsorted(ss_t, cs_t[k] + PSTH_WIN[1])
        o["raster"] = [np.round((ss_t[x:y] - cs_t[j]) * 1e4).astype(int).tolist() for j, x, y in zip(k, i0, i1)]
        # quarters of the pair window, each ÷ its own baseline
        qs = []
        q_edges = np.linspace(w[0], w[1], 5)
        for qi in range(4):
            m = (cs_t >= q_edges[qi]) & (cs_t < q_edges[qi + 1])
            bq = C[m][:, bmask].sum() / (max(m.sum(), 1) * bmask.sum() * BIN)
            qs.append(dict(t=[round(q_edges[qi]), round(q_edges[qi + 1])], n_cs=int(m.sum()), base_hz=round(float(bq), 2),
                           psth=np.round(C[m][:, pmask].mean(0) / BIN / bq, 3).tolist() if m.sum() and bq > 0 else None))
        o["quarters"] = qs
        o["dist_um"] = round(float(r.dist_um), 1)
        pairs_out[f"{cs}_{ss}"] = o
    E = pd.DataFrame(eph)
    out_dir = CUR / sess / "pause_figure"; out_dir.mkdir(parents=True, exist_ok=True)
    E.round(4).to_csv(out_dir / "ephaptic_pairs.csv", index=False)

    def group_curves(rows):
        g = []
        for lo, hi, lab in DIST_GROUPS:
            Y = np.array([y for _, d, y in rows if lo < d <= hi or (lo == 0 and d == 0)])
            g.append(dict(label=lab, n=len(Y),
                          mean=np.round(Y.mean(0), 3).tolist() if len(Y) else None,
                          sem=np.round(Y.std(0, ddof=1) / np.sqrt(len(Y)), 3).tolist() if len(Y) > 1 else None))
        return g
    per_cs = {int(c): group_curves([x for x in curves if x[0] == c]) for c in P.cs.unique()}
    js = dict(created=time.strftime("%Y-%m-%d %H:%M"), session=sess,
              psth_t0_ms=PSTH_WIN[0] * 1e3, zoom_t0_ms=ZOOM_WIN[0] * 1e3, group_t0_ms=GROUP_WIN[0] * 1e3,
              bin_ms=BIN * 1e3, zbin_ms=ZBIN * 1e3, base_win_ms=[base_win[0] * 1e3, base_win[1] * 1e3],
              early_win_ms=[EARLY_WIN[0] * 1e3, EARLY_WIN[1] * 1e3], pause_win_ms=[PAUSE_WIN[0] * 1e3, PAUSE_WIN[1] * 1e3],
              n_boot=N_BOOT, n_raster=N_RASTER, groups=[g[2] for g in DIST_GROUPS],
              control=control, pairs=pairs_out,
              eph=dict(cs=E.cs.tolist(), ss=E.ss.tolist(), d=E.dist_um.round(1).tolist(), r02=E.r_0_2.round(3).tolist(),
                       r38=E.r_3_8.round(3).tolist(), partner=E.partner.astype(int).tolist(), flagged=E.flagged.astype(int).tolist()),
              eph_cs=per_cs, eph_session=group_curves(curves))
    txt = "window.CS_PAUSEFIG = " + json.dumps(js, separators=(",", ":")).replace("NaN", "null") + ";\n"
    (out_dir / "pause_figure.js").write_text(txt, encoding="utf-8")
    print(f"{sess}: {len(pairs_out)} pairs in full, {len(E)} ephaptic pairs, {len(txt) / 1e6:.1f} MB, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    sessions = sys.argv[1:] or sorted(d.parent.name for d in CUR.glob(f"*/{RUNS['20 Hz']}"))
    for s in sessions:
        run(s)
