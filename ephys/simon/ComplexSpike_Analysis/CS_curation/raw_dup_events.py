# -*- coding: utf-8 -*-
"""
raw_dup_events.py — SINGLE raw events of CS–CS pairs, one unit's spike after the other (add-on to the raw layer).

Why: the raw-layer means (raw_waveforms_session.py, section 7b) average each pair's events at the peak lag, which hides
HOW the two detections follow each other in one event. For a suspected spikelet pair (B = a later piece of A's event)
we want to see single events where both units were detected, in a row: A's spike, then B's spike at its own lag, on
both units' peak channels, next to events where only one of them fired.

For every CS–CS pair in step0's cs_cs_evidence.csv (same pairs as the raw layer):
  - the pair is ordered so that "first" is the unit that fires first at the CCG peak (peak lag >= 0 after swapping)
  - coincident events: first-unit spikes with a second-unit spike within peak lag ± COINC_DELTA (step0's definition)
  - first-only / second-only events (no partner at the peak lag); second-only snippets are placed so that its spike
    sits at the peak lag (same time frame as the coincident events)
  - per event: raw snippet WIN around the first unit's spike time (= template trough), on the channels where the two
    units are largest (CH_EACH per unit), plus the times of ALL spikes of either unit inside the window (so double
    detections and spikelet-like repeats show up as extra marks)
  - all coincident lags (histogram) and the fraction of coincident events with >= 2 spikes of the same unit in WIN
Preprocessing = raw layer: int16 → µV, per-channel mean removed, common median reference over all 384 channels,
300 Hz high-pass (Butterworth 3, zero-phase). Read-only on Dropbox (only the .meta file is read); the raw .bin is the
LOCAL copy. Output: CS_curation/<mouse>_<R#>/raw_waveforms/raw_dup_events.js (window.CS_RAW_DUPEV), loaded by
browser.html next to raw_data.js. Nothing is merged; information only.
"""
# %% 1. PARAMETERS + PATHS
from pathlib import Path
import os as _os
import json
import time

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

LOCAL_ROOT = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior")
MOUSE, SESSION_R = "1_CoteDor", "R5"
MOUSE, SESSION_R = (_os.environ["RAW_SESSION"].rsplit("_", 1) if _os.environ.get("RAW_SESSION") else (MOUSE, SESSION_R))  # batch: raw_all_sessions.py
STEP0_DIR = LOCAL_ROOT / "CS_curation" / f"{MOUSE}_{SESSION_R}" / "step0_explore_bhvspan_v7batch"
RAW_BIN = (LOCAL_ROOT / "raw_copies" / "1_CoteDor" / "04082023_CoteDor_Cer_S5_g0_imec0"
           / "04082023_CoteDor_Cer_S5_g0_t0.imec0.ap.bin")
RAW_META = (Path(r"D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_data"
                 r"\20230801_ChocolateGroup\1_CoteDor\04082023_CoteDor_Cer_S5_g0\04082023_CoteDor_Cer_S5_g0_imec0")
            / "04082023_CoteDor_Cer_S5_g0_t0.imec0.ap.meta")             # read only
if _os.environ.get("RAW_BIN"):                       # batch: local copy + its local .meta (raw_all_sessions.py)
    RAW_BIN, RAW_META = Path(_os.environ["RAW_BIN"]), Path(_os.environ["RAW_META"])
OUT_DIR = LOCAL_ROOT / "CS_curation" / f"{MOUSE}_{SESSION_R}" / "raw_waveforms"

WIN = (-3e-3, 10e-3)          # s around the FIRST unit's spike (wide enough for both spikes and the spikelets after)
PAD = 7e-3                    # s read on each side for the filter
HP_HZ, HP_ORDER = 300.0, 3
N_COINC = 40                  # coincident events stored per pair (random, then sorted by the measured lag)
N_SINGLE = 15                 # first-only and second-only events stored per pair
CH_EACH = 2                   # channels per unit (largest mean ptp), union → up to 4 channels
N_MEAN_CH = 200               # spikes per unit used to pick its largest channels
RNG = np.random.default_rng(0)

# %% 2. LOAD
p0 = json.loads((STEP0_DIR / "params_used.json").read_text())
SORTER_DIR = Path(p0["sorter_dir"])
COINC_D = float(p0.get("COINC_DELTA", 0.5e-3))
meta = dict(l.strip().split("=", 1) for l in RAW_META.read_text().splitlines() if "=" in l)
FS = float(meta["imSampRate"]); N_CH_SAVED = int(meta["nSavedChans"])
_gain = float(meta["~imroTbl"].split(")(")[1].split()[3])
UV = float(meta["imAiRangeMax"]) / float(meta["imMaxInt"]) / _gain * 1e6
assert RAW_BIN.exists() and RAW_BIN.stat().st_size == int(meta["fileSizeBytes"]), "raw copy missing or incomplete"
raw = np.memmap(RAW_BIN, dtype=np.int16, mode="r").reshape(-1, N_CH_SAVED)
N_AP = 384
chan_pos = np.load(SORTER_DIR / "channel_positions.npy")
chan_map = np.load(SORTER_DIR / "channel_map.npy").ravel()
spike_samp = np.asarray(np.load(SORTER_DIR / "spike_times.npy", mmap_mode="r")).ravel().astype(np.int64)
spike_clu = np.asarray(np.load(SORTER_DIR / "spike_clusters.npy", mmap_mode="r")).ravel()
_evf = STEP0_DIR / "tables" / "cs_cs_evidence.csv"                # empty file (just a newline) when step0 found no CS–CS pairs
ev0 = pd.read_csv(_evf) if _evf.stat().st_size > 5 else pd.DataFrame(columns=["cs_a", "cs_b"])
need = set(ev0["cs_a"]).union(ev0["cs_b"])
SAMP = {int(c): np.sort(spike_samp[spike_clu == c]) for c in need}
SOS = butter(HP_ORDER, HP_HZ, "highpass", fs=FS, output="sos")
S0, S1, SP = int(round(WIN[0] * FS)), int(round(WIN[1] * FS)), int(round(PAD * FS))
OFF = np.arange(S0 - SP, S1 + SP)
CROP = slice(SP, SP + (S1 - S0))
T_MS = np.arange(S0, S1) / FS * 1e3
print(f"raw {raw.shape}, {FS:.0f} Hz | {len(ev0)} CS–CS pairs | step0 {STEP0_DIR.name}")


def snippets(samples, chs):
    """(n,) samples → (n, L, len(chs)) µV, CAR over all channels, high-passed, cropped to WIN."""
    s = np.asarray(samples, dtype=np.int64)
    s = s[(s + OFF[0] >= 0) & (s + OFF[-1] < raw.shape[0])]
    out = []
    for b in range(0, len(s), 128):
        x = raw[s[b:b + 128, None] + OFF[None, :], :N_AP].astype(np.float32)[:, :, chan_map]
        x -= x.mean(axis=1, keepdims=True)
        x -= np.median(x, axis=2, keepdims=True)
        out.append(sosfiltfilt(SOS, x[:, :, chs] * UV, axis=1)[:, CROP])
    return np.concatenate(out) if out else np.zeros((0, len(T_MS), len(chs)))


def top_ch(c, w):
    s = SAMP[c][(SAMP[c] >= w[0] * FS) & (SAMP[c] < w[1] * FS)]
    s = s if len(s) <= N_MEAN_CH else RNG.choice(s, N_MEAN_CH, replace=False)
    m = snippets(np.sort(s), np.arange(len(chan_map))).mean(0)
    return list(np.argsort(np.ptp(m, 0))[::-1][:CH_EACH])


def pick(idx, n):
    return np.sort(RNG.choice(idx, n, replace=False)) if len(idx) > n else idx


def marks(t0, tt):
    """Spike times of train tt (s) inside WIN of reference t0 (s), in ms relative to t0."""
    i0, i1 = np.searchsorted(tt, t0 + WIN[0]), np.searchsorted(tt, t0 + WIN[1])
    return [round(float((v - t0) * 1e3), 3) for v in tt[i0:i1]]


# %% 3. PER PAIR
_t = time.time()
OUT = []
for r in ev0.itertuples():
    a, b, L = int(r.cs_a), int(r.cs_b), float(r.peak_lag_ms) * 1e-3
    first, second, lag = (a, b, L) if L >= 0 else (b, a, -L)          # "first" fires first at the peak
    w = (r.win_start, r.win_end)
    tf = SAMP[first] / FS; tf = tf[(tf >= w[0]) & (tf < w[1])]
    ts = SAMP[second] / FS; ts = ts[(ts >= w[0]) & (ts < w[1])]
    if len(tf) < 5 or len(ts) < 5:
        continue
    j = np.searchsorted(ts, tf + lag - COINC_D)
    jj = np.clip(j, 0, len(ts) - 1)
    co = (j < len(ts)) & (ts[jj] <= tf + lag + COINC_D)
    co_lag = (ts[jj[co]] - tf[co]) * 1e3                                 # measured lag of every coincident event
    k = np.searchsorted(tf, ts - lag - COINC_D); kk = np.clip(k, 0, len(tf) - 1)
    co_s = (k < len(tf)) & (tf[kk] <= ts - lag + COINC_D)
    # how often a coincident event contains a further spike of either unit (double detection / repeat)
    rep_f = np.array([len(marks(t, tf)) > 1 for t in tf[co]]) if co.any() else np.array([])
    rep_s = np.array([len(marks(t, ts)) > 1 for t in tf[co]]) if co.any() else np.array([])
    chs = sorted(set(top_ch(first, w)) | set(top_ch(second, w)), key=lambda q: (chan_pos[q, 1], chan_pos[q, 0]))
    sets = {}
    for name, refs, part in (("coincident", tf[co], None), ("first-only", tf[~co], None), ("second-only", ts[~co_s] - lag, None)):
        n = N_COINC if name == "coincident" else N_SINGLE
        sel = pick(np.arange(len(refs)), n)
        if not len(sel):
            continue
        ref = refs[sel]
        x = snippets(np.round(ref * FS).astype(np.int64), chs)
        mf = [marks(t, tf) for t in ref]
        msd = [marks(t, ts) for t in ref]
        order = np.argsort([min([v for v in m2 if abs(v - lag * 1e3) <= COINC_D * 1e3 + 1e-9] or [99]) for m2 in msd]) \
            if name == "coincident" else np.arange(len(ref))
        sets[name] = dict(n_total=int(len(refs)), t_s=[round(float(v), 4) for v in ref[order]],
                          y=[[[int(round(float(v))) for v in x[e, :, c]] for c in range(len(chs))] for e in order],
                          first=[mf[e] for e in order], second=[msd[e] for e in order])
    OUT.append(dict(cs_a=a, cs_b=b, first=first, second=second, lag_ms=round(lag * 1e3, 3), status=r.evidence_status,
                    coinc=round(float(r.coinc_frac_peak), 3), ch=[int(c) for c in chs],
                    x=[float(chan_pos[c, 0]) for c in chs], y=[float(chan_pos[c, 1]) for c in chs],
                    peak_first=int(top_ch(first, w)[0]), peak_second=int(top_ch(second, w)[0]),
                    n_first=int(len(tf)), n_second=int(len(ts)), n_coinc=int(co.sum()),
                    lag_hist=np.histogram(co_lag, np.arange(lag * 1e3 - COINC_D * 1e3, lag * 1e3 + COINC_D * 1e3 + 1e-9, 0.05))[0].tolist(),
                    lag_hist_lo=round(lag * 1e3 - COINC_D * 1e3, 3),
                    frac_repeat_first=float(rep_f.mean()) if len(rep_f) else None,
                    frac_repeat_second=float(rep_s.mean()) if len(rep_s) else None, sets=sets))
    print(f"  {first} → {second}: lag {lag * 1e3:+.2f} ms, {co.sum()} coincident of {len(tf)} / {len(ts)}, channels {chs} "
          f"({time.time() - _t:.0f} s)")

# %% 4. SAVE
js = "window.CS_RAW_DUPEV = " + json.dumps(dict(t_ms=[round(float(v), 4) for v in T_MS], win_ms=[WIN[0] * 1e3, WIN[1] * 1e3],
                                                 coinc_delta_ms=COINC_D * 1e3, step0_dir=STEP0_DIR.name,
                                                 created=time.strftime("%Y-%m-%d %H:%M"), pairs=OUT),
                                            separators=(",", ":")) + ";\n"
(OUT_DIR / "raw_dup_events.js").write_text(js, encoding="utf-8")
print(f"saved {OUT_DIR / 'raw_dup_events.js'} ({len(js) / 1e6:.1f} MB, {len(OUT)} pairs, {time.time() - _t:.0f} s)")
