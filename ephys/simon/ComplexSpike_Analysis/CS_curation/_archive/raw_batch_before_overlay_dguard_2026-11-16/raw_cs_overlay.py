# -*- coding: utf-8 -*-
"""
raw_cs_overlay.py — raw traces around CS events, overlaid (like a classic "CS + SS pause" panel: many sweeps aligned
to the CS, the CS waveform at 0 and no SS spikes for ~10 ms after it).

Per session with a raw file on this computer (local copy in raw_copies/, or a Dropbox .ap.bin that is already on the
device — never one that is online-only, so nothing is downloaded), per CS that has a partner SS (candidate in either
step0 floor; if none, its step0 #1 SS):
  * N_EV random CS events inside the pair window of its best partner (time order), −10…+30 ms
  * traces on the CS peak channel and on the peak channel of each partner SS (≤ MAX_SS, by step0 rank) and of the
    CS's control SS (pause_figure.js: nearest non-partner, non-flagged SS)
  * preprocessing as raw_waveforms_session.py: int16 → µV, per-channel mean removed, common median reference over all
    channels, 300–6000 Hz band-pass (Butterworth 3, zero-phase; the low-pass is for display only); stored at 15 kHz (every 2nd sample) for the browser
  * sorter spike times of the CS and of each SS in the same window (ticks under the traces)
Output (local only): CS_curation/<session>/raw_overlay/raw_overlay.js (window.CS_RAWOVERLAY, loaded by browser.html)
                     CS_curation/<session>/raw_overlay/figures/cs<CS>_ss<SS>.png (best partner)
Dropbox is only READ (no attrib, no pin). Usage: python raw_cs_overlay.py [session ...]
"""
from pathlib import Path
import json
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pause_figure_data import load_pairs, load_session          # noqa: E402  (pure functions, no side effects)
from raw_all_sessions import resolve_bin, is_local, meta_of      # noqa: E402

LOCAL_ROOT = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior")
CUR = LOCAL_ROOT / "CS_curation"
RAW_COPIES = LOCAL_ROOT / "raw_copies"
RUN20 = "step0_explore_bhvspan_v7batch"

WIN = (-10e-3, 30e-3)
PAD = 7e-3
HP_HZ, HP_ORDER = 300.0, 3
LP_HZ = 6000.0               # also low-pass for display (less high-frequency noise; SS/CS spikes ~1–3 kHz)
N_EV = 40
DEC = 2                      # keep every 2nd sample (30 → 15 kHz) for the browser
MAX_SS = 4
N_AP = 384
RNG = np.random.default_rng(0)


def raw_source(sess):
    b = resolve_bin(sess)
    local = sorted(RAW_COPIES.rglob(b.name))
    if local and local[0].stat().st_size == b.stat().st_size:
        return local[0], local[0].with_suffix(".meta")
    if is_local(b):
        return b, b.with_suffix(".meta")
    return None, None


def run(sess):
    t0 = time.time()
    src, meta_f = raw_source(sess)
    if src is None:
        print(f"{sess}: no raw file on this computer (online-only in Dropbox) — skipped", flush=True)
        return
    meta = dict(l.strip().split("=", 1) for l in meta_f.read_text().splitlines() if "=" in l)
    FS = float(meta["imSampRate"]); NCH = int(meta["nSavedChans"])
    gain = float(meta["~imroTbl"].split(")(")[1].split()[3])
    UV = float(meta["imAiRangeMax"]) / float(meta["imMaxInt"]) / gain * 1e6
    raw = np.memmap(src, dtype=np.int16, mode="r").reshape(-1, NCH)
    p0, ST = load_session(sess)
    sd = Path(p0["sorter_dir"])
    chan_map = np.load(sd / "channel_map.npy").ravel()
    units = pd.read_csv(CUR / sess / RUN20 / "tables" / "unit_table.csv").set_index("phyID")
    P = load_pairs(sess)
    pf_js = CUR / sess / "pause_figure" / "pause_figure.js"
    control = {}
    if pf_js.exists():
        t = pf_js.read_text(encoding="utf-8"); control = {int(k): v for k, v in json.loads(t[t.index("=") + 1:].rstrip().rstrip(";"))["control"].items()}

    sos = butter(HP_ORDER, [HP_HZ, LP_HZ], "bandpass", fs=FS, output="sos")
    s0, s1, sp = int(round(WIN[0] * FS)), int(round(WIN[1] * FS)), int(round(PAD * FS))
    off = np.arange(s0 - sp, s1 + sp)
    out = dict(created=time.strftime("%Y-%m-%d %H:%M"), session=sess, raw_file=src.name,
               raw_where="local copy" if str(src).startswith(str(RAW_COPIES)) else "Dropbox (on device, read only)",
               t0_ms=WIN[0] * 1e3, dt_ms=1e3 * DEC / FS, hp_hz=HP_HZ, lp_hz=LP_HZ, n_ev=N_EV, cs={})
    fig_dir = CUR / sess / "raw_overlay" / "figures"; fig_dir.mkdir(parents=True, exist_ok=True)
    for cs, g in P.groupby("cs"):
        cs = int(cs)
        part = g[g.partner].sort_values("rank")
        if not len(part):                                   # no partner: the step0 #1 SS (shown as "not a candidate")
            part = g[g.ss.isin(ST.keys())].sort_values("rank").head(1)
        if not len(part) or cs not in units.index:
            continue
        ss_list = [int(s) for s in part.ss.head(MAX_SS)]
        if control.get(cs) is not None and control[cs] not in ss_list:
            ss_list.append(int(control[cs]))
        best = part.iloc[0]
        w = (float(best.win_start), float(best.win_end))
        a = ST[cs]; ev = a[(a >= w[0] - WIN[0] + PAD) & (a < w[1] - WIN[1] - PAD)]
        if len(ev) < 5:
            continue
        ev = np.sort(RNG.choice(ev, min(N_EV, len(ev)), replace=False))
        chans = {int(units.at[cs, "peak_ch"])} | {int(units.at[s, "peak_ch"]) for s in ss_list if s in units.index}
        chans = sorted(chans)
        samp = np.round(ev * FS).astype(np.int64)
        x = raw[samp[:, None] + off[None, :], :N_AP].astype(np.float32)[:, :, chan_map]       # sorter channel order
        x -= x.mean(axis=1, keepdims=True)
        x -= np.median(x, axis=2, keepdims=True)
        x = sosfiltfilt(sos, x[:, :, chans] * UV, axis=1)[:, sp:sp + (s1 - s0), :]
        tr = {str(ch): np.round(x[:, ::DEC, k]).astype(int).tolist() for k, ch in enumerate(chans)}
        ticks = lambda u: [np.round((ST[u][(ST[u] >= e + WIN[0]) & (ST[u] < e + WIN[1])] - e) * 1e4).astype(int).tolist() for e in ev]
        ss_out = {}
        for s in ss_list:
            if s not in units.index or s not in ST:
                continue
            r = g[g.ss == s].iloc[0]
            ss_out[str(s)] = dict(ch=int(units.at[s, "peak_ch"]), partner=bool(r.partner), control=s == control.get(cs),
                                  dist_um=round(float(r.dist_um), 1), rank=int(r["rank"]), ticks=ticks(s))
        out["cs"][str(cs)] = dict(ch=int(units.at[cs, "peak_ch"]), win=[round(w[0]), round(w[1])], n_cs_win=int(
            ((a >= w[0]) & (a < w[1])).sum()), t_ev=np.round(ev, 3).tolist(), traces=tr, cs_ticks=ticks(cs), ss=ss_out)
        # static figure: CS channel and best partner's channel (like a classic CS + SS-pause sweep overlay)
        sb = ss_list[0]; chs = [int(units.at[cs, "peak_ch"]), int(units.at[sb, "peak_ch"])]
        tm = WIN[0] * 1e3 + np.arange(x.shape[1]) / FS * 1e3
        ncol = 2 if chs[0] != chs[1] else 1
        fig, ax = plt.subplots(2, ncol, figsize=(4.6 * ncol, 5.2), squeeze=False, gridspec_kw=dict(height_ratios=[2.2, 1]), sharex=True)
        for j, ch in enumerate(chs[:ncol]):
            k = chans.index(ch); A = ax[0, j]
            A.plot(tm, x[:, :, k].T, color="k", lw=0.35, alpha=0.18)
            A.plot(tm, x[:, :, k].mean(0), color="#b3261e", lw=1.0)
            A.axvspan(3, 8, color="#e8e8e8", zorder=0); A.axvline(0, color="#b3261e", lw=0.6, ls="--")
            R = ax[1, j]
            for i, e in enumerate(ev):
                for u, col in ((cs, "#b3261e"), (sb, "k")):
                    ts = (ST[u][(ST[u] >= e + WIN[0]) & (ST[u] < e + WIN[1])] - e) * 1e3
                    R.vlines(ts, i + 0.1, i + 0.9, color=col, lw=0.8)
            R.axvspan(3, 8, color="#e8e8e8", zorder=0); R.set_ylim(0, len(ev)); R.set_ylabel("CS event"); R.set_xlabel("time from CS (ms)")
            R.spines[["top", "right"]].set_visible(False)
            A.set_xlim(WIN[0] * 1e3, WIN[1] * 1e3); A.set_ylabel(f"µV ({HP_HZ:.0f}–{LP_HZ:.0f} Hz)")
            A.set_title(f"ch {ch} ({'CS ' + str(cs) if ch == chs[0] else ''}{' / ' if ch == chs[0] == chs[1] else ''}{'SS ' + str(sb) if ch == chs[1] else ''} peak)", fontsize=9)
            A.spines[["top", "right"]].set_visible(False)
        fig.suptitle(f"{sess}: {len(ev)} CS {cs} events, raw traces overlaid, red = mean (SS {sb}: step0 rank {int(best['rank'])}, {best.dist_um:.0f} µm)\n"
                     f"below: sorter spikes, red = CS {cs}, black = SS {sb}; grey = step0 pause window 3–8 ms", fontsize=9)
        fig.tight_layout(); fig.savefig(fig_dir / f"cs{cs}_ss{sb}.png", dpi=130); plt.close(fig)
    od = CUR / sess / "raw_overlay"
    txt = "window.CS_RAWOVERLAY = " + json.dumps(out, separators=(",", ":")) + ";\n"
    (od / "raw_overlay.js").write_text(txt, encoding="utf-8")
    print(f"{sess}: {len(out['cs'])} CS, {len(txt) / 1e6:.1f} MB, {out['raw_where']}, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    sessions = sys.argv[1:] or ["1_CoteDor_R5", "1_CoteDor_R1", "5_FerreroRocher_R4"]
    for s in sessions:
        run(s)
