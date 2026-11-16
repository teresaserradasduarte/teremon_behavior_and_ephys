# -*- coding: utf-8 -*-
"""
raw_waveforms_session.py — CS events seen in the RAW AP data (one session; first pass: CoteDor R5).

Why: sorter templates span only ±1.3 ms, so spikelets (≈1–8 ms after the CS main spike) and CS pieces absorbed by
other clusters cannot be judged from them. This script reads short snippets of the raw AP band around spike times
and builds, for every CS and every cluster its event spreads into:
  - mean waveforms −2…+8 ms (channels within ±MEAN_CH_UM of the unit's peak) + a peak-to-peak profile over the probe
  - conditional means: CS events with vs without a spike of each "capturing" cluster
  - capturing clusters' own spikes split into CS-locked vs other (e.g. R5 SS 938)
  - single CS events on the CS peak channel (with which clusters fired in the event window)
  - a per-event spikelet count from the raw trace (negative peaks after the main spike)
Nothing is merged, nothing changes step0's decisions. Read-only on Dropbox (spike/cluster npy files are only read);
the raw .bin is read from a LOCAL copy (RAW_BIN). Outputs go to CS_curation/<mouse>_<R#>/raw_waveforms/.

Preprocessing per snippet (reads −2−PAD … +8+PAD ms, all 384 channels): int16 → µV, per-channel mean removed,
common median reference across channels, 300 Hz high-pass (Butterworth order 3, zero-phase). Means are built from
the CAR'd snippets and filtered once (the filter is linear, so this equals averaging filtered snippets).
"""
# %% 1. PARAMETERS + PATHS
from pathlib import Path
import json
import time

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt, find_peaks

LOCAL_ROOT = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior")
MOUSE, SESSION_R = "1_CoteDor", "R5"
STEP0_DIR = LOCAL_ROOT / "CS_curation" / f"{MOUSE}_{SESSION_R}" / "step0_explore_bhvspan_v7batch"   # 2026-11-09: v7 (was v5batch)
RAW_BIN = (LOCAL_ROOT / "raw_copies" / "1_CoteDor" / "04082023_CoteDor_Cer_S5_g0_imec0"
           / "04082023_CoteDor_Cer_S5_g0_t0.imec0.ap.bin")
RAW_META = (Path(r"D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_data"
                 r"\20230801_ChocolateGroup\1_CoteDor\04082023_CoteDor_Cer_S5_g0\04082023_CoteDor_Cer_S5_g0_imec0")
            / "04082023_CoteDor_Cer_S5_g0_t0.imec0.ap.meta")             # read only
OUT_DIR = LOCAL_ROOT / "CS_curation" / f"{MOUSE}_{SESSION_R}" / "raw_waveforms"

WIN = (-2e-3, 8e-3)              # s around the spike time (spike time = template trough)
PAD = 7e-3                       # s read on each side for the filter, cropped afterwards
HP_HZ, HP_ORDER = 300.0, 3
NEIGH_UM = 250.0                 # clusters within this distance of a CS can be "capturing" clusters / SS candidates
CAPTURE_WIN = (-0.5e-3, 4e-3)    # s after the CS: a spike of another cluster here belongs to the CS event?
CAPTURE_MIN_EXCESS = 0.05        # fraction of CS events above Poisson chance
CAPTURE_MIN_N = 30               # and at least this many CS events with a spike of the cluster
N_MEAN = 500                     # random spikes per unit for its mean waveform
N_SINGLE = 40                    # single CS events stored per CS
N_COND_MIN = 20                  # min events on each side for a conditional mean
MEAN_CH_UM = 100.0               # stored channels: |y - y_peak| <= this
COND_N_CH = 8                    # conditional / locked means: top-ptp channels stored
SPK_WIN = (0.5e-3, 8e-3)         # spikelet search window after the CS main spike
SPK_THR = 3.5                    # × robust noise SD of the residuals (event − mean CS waveform), pre-CS window
SPK_PROM = 2.0                   # × noise SD: minimum prominence (so the tail of the main spike is not counted)
SPK_MIN_SEP = 0.4e-3             # s between detected spikelets
N_NULL = 1500                    # random-time snippets per CS: same detector → false-positive "spikelets"/event
NULL_EXCL = 0.02                 # s: random times at least this far from any spike of the CS
BATCH = 192
RNG = np.random.default_rng(0)

OUT_DIR.mkdir(parents=True, exist_ok=True)
params_used = {k: v for k, v in dict(globals()).items() if k.isupper() and isinstance(v, (int, float, str, tuple))}
params_used = {k: (str(v) if isinstance(v, Path) else v) for k, v in params_used.items()}

# %% 2. LOAD (step0 tables, sorter spikes, raw meta)
p0 = json.loads((STEP0_DIR / "params_used.json").read_text())
SORTER_DIR = Path(p0["sorter_dir"])
units = pd.read_csv(STEP0_DIR / "tables" / "unit_table.csv").set_index("phyID")
meta = dict(l.strip().split("=", 1) for l in RAW_META.read_text().splitlines() if "=" in l)
FS = float(meta["imSampRate"]); N_CH_SAVED = int(meta["nSavedChans"])
_gain = float(meta["~imroTbl"].split(")(")[1].split()[3])
UV = float(meta["imAiRangeMax"]) / float(meta["imMaxInt"]) / _gain * 1e6
assert RAW_BIN.exists() and RAW_BIN.stat().st_size == int(meta["fileSizeBytes"]), "raw copy missing or incomplete"
raw = np.memmap(RAW_BIN, dtype=np.int16, mode="r").reshape(-1, N_CH_SAVED)
N_AP = 384
chan_pos = np.load(SORTER_DIR / "channel_positions.npy")
chan_map = np.load(SORTER_DIR / "channel_map.npy").ravel()             # sorter channel -> raw channel
assert len(chan_map) == N_AP
print(f"raw {raw.shape}, {FS:.0f} Hz, {UV:.3f} µV/bit | sorter {SORTER_DIR.name}")

spike_samp = np.asarray(np.load(SORTER_DIR / "spike_times.npy", mmap_mode="r")).ravel().astype(np.int64)
spike_clu = np.asarray(np.load(SORTER_DIR / "spike_clusters.npy", mmap_mode="r")).ravel()
_o = np.argsort(spike_clu, kind="stable")
_ids, _st = np.unique(spike_clu[_o], return_index=True)
_en = np.r_[_st[1:], len(_o)]
SAMP = {int(c): np.sort(spike_samp[_o[a:b]]) for c, a, b in zip(_ids, _st, _en)}


def span(c):
    return float(units.at[c, "span_start"]), float(units.at[c, "span_end"])


def times_in(c, w):
    s = SAMP[c] / FS
    return s[(s >= w[0]) & (s < w[1])]


def dist(a, b):
    return float(np.hypot(units.at[a, "x"] - units.at[b, "x"], units.at[a, "y"] - units.at[b, "y"]))


# %% 3. UNIT SELECTION: CS, SS candidates near a CS, clusters that capture CS events
CS_IDS = units.index[units["role"] == "CS"].tolist()
SS_NEAR = [u for u in units.index[units["role"].isin(["SS_cand", "SS_cand2"])]
           if min(dist(u, c) for c in CS_IDS) <= NEIGH_UM]
cap_rows = []
for c in CS_IDS:
    for k in units.index:
        if k == c or k not in SAMP or dist(c, k) > NEIGH_UM:
            continue
        w = (max(span(c)[0], span(k)[0]), min(span(c)[1], span(k)[1]))
        if w[1] - w[0] < 60:
            continue
        cs, sk = times_in(c, w), times_in(k, (w[0] - 0.01, w[1] + 0.01))
        if len(cs) < 20 or len(sk) == 0:
            continue
        i0, i1 = np.searchsorted(sk, cs + CAPTURE_WIN[0]), np.searchsorted(sk, cs + CAPTURE_WIN[1])
        hit = i1 > i0
        rate = len(sk) / (w[1] - w[0] + 0.02)
        null = 1 - np.exp(-rate * (CAPTURE_WIN[1] - CAPTURE_WIN[0]))
        if hit.sum() < CAPTURE_MIN_N or hit.mean() - null < CAPTURE_MIN_EXCESS:
            continue
        lg = sk[i0[hit]] - cs[hit]                                       # first spike in the window
        kin = sk[(sk >= w[0]) & (sk < w[1])]
        j = np.searchsorted(cs, kin - CAPTURE_WIN[1])
        locked_k = (j < len(cs)) & (cs[np.clip(j, 0, len(cs) - 1)] <= kin - CAPTURE_WIN[0])
        cap_rows.append(dict(cs=c, unit=k, role=units.at[k, "role"] if isinstance(units.at[k, "role"], str) else "",
                             group=units.at[k, "group"], c4=units.at[k, "c4_label"], dist_um=dist(c, k),
                             fr=units.at[k, "fr"], win_start=w[0], win_end=w[1], n_cs=len(cs), n_hit=int(hit.sum()),
                             frac_cs_events=hit.mean(), chance=null, excess=hit.mean() - null,
                             frac_unit_locked=locked_k.mean(), lag_med_ms=np.median(lg) * 1e3,
                             lag_q25_ms=np.quantile(lg, .25) * 1e3, lag_q75_ms=np.quantile(lg, .75) * 1e3))
cap = pd.DataFrame(cap_rows).sort_values(["cs", "excess"], ascending=[True, False])
cap.to_csv(OUT_DIR / "cs_event_capture.csv", index=False)
CAP_UNITS = sorted(set(cap["unit"]))
ALL_UNITS = sorted(set(CS_IDS) | set(SS_NEAR) | set(CAP_UNITS))
print(f"{len(CS_IDS)} CS, {len(SS_NEAR)} SS candidates within {NEIGH_UM:.0f} µm, {len(cap)} CS×cluster capture pairs "
      f"({len(CAP_UNITS)} clusters) -> {len(ALL_UNITS)} units")

# %% 3b. PAIRING RE-CHECK (spike times only): pause measured AFTER the CS event window
# v5's pause window (0–8 ms) contains the CS event itself (main spike + spikelets, ≈ 0–4+ ms). A unit that absorbed
# CS pieces gets extra spikes there and looks like it pauses less (R5: 938 vs 925). step0 v7 uses 3–8 ms (Simon
# 2026-11-07; v6's 3–8 ms reached the post-pause rebound). Here: 0–8 (v5), 3–8 (= v7) and 8–16 ms (late) with the same
# baseline as step0; post_rank = rank by 3–8 ms. Information only, nothing re-decided.
PR_BASE = tuple(p0.get("BASELINE_WIN", (-0.05, -0.01)))
PR_WINS = {"pr_0_8": (0.0, 0.008), "pr_3_8": (0.003, 0.008), "pr_8_16": (0.008, 0.016)}
RC_CCG = (-0.02, 0.05, 0.0005)                    # exported CCG: lag range and bin (s)
RC_TOP = 6                                        # per CS: export CCGs for the top-N by step0 rank and by pr_3_8


def n_in(ref, tgt, lo, hi):
    return int((np.searchsorted(tgt, ref + hi) - np.searchsorted(tgt, ref + lo)).sum())


pairs0 = pd.read_csv(STEP0_DIR / "tables" / "cs_ss_pairs.csv")
rc_rows = []
for p in pairs0.itertuples():
    c, k = int(p.cs), int(p.ss)
    if c not in SAMP or k not in SAMP:
        continue
    w = (p.win_start, p.win_end)
    cs, sk = times_in(c, w), times_in(k, (w[0] - 0.1, w[1] + 0.1))
    if len(cs) < 20:
        continue
    base = n_in(cs, sk, *PR_BASE) / (len(cs) * (PR_BASE[1] - PR_BASE[0]))
    row = dict(cs=c, ss=k, step0_rank=int(p.rank), step0_candidate=bool(p.candidate), ss_tier=int(p.ss_tier),
               ss_group=p.ss_group, dist_um=p.dist_um, ss_fr=p.ss_fr, n_cs=len(cs), baseline=base,
               step0_pause_ratio=p.pause_ratio, cs_share_excess=p.cs_share_excess, jitter_p=p.jitter_p)
    for name, (lo, hi) in PR_WINS.items():
        row[name] = n_in(cs, sk, lo, hi) / (len(cs) * (hi - lo)) / base if base > 0 else np.nan
    rc_rows.append(row)
recheck = pd.DataFrame(rc_rows)
recheck["post_rank"] = recheck.groupby("cs")["pr_3_8"].rank(method="first")
recheck = recheck.sort_values(["cs", "post_rank"])
recheck.to_csv(OUT_DIR / "pairing_recheck.csv", index=False)
RC_CCGS = {}
_edges = np.arange(RC_CCG[0], RC_CCG[1] + 1e-12, RC_CCG[2])
for c, g in recheck.groupby("cs"):
    keep = set(g.nsmallest(RC_TOP, "step0_rank")["ss"]) | set(g.nsmallest(RC_TOP, "pr_3_8")["ss"])
    for r in g[g["ss"].isin(keep)].itertuples():
        cs = times_in(c, (r_w := pairs0.loc[(pairs0.cs == c) & (pairs0.ss == r.ss), ["win_start", "win_end"]].iloc[0].tolist()))
        sk = times_in(r.ss, (r_w[0] - 0.1, r_w[1] + 0.1))
        i0, i1 = np.searchsorted(sk, cs + RC_CCG[0]), np.searchsorted(sk, cs + RC_CCG[1])
        lg = np.concatenate([sk[a:b] - t for t, a, b in zip(cs, i0, i1)]) if len(cs) else np.array([])
        h = np.histogram(lg, _edges)[0] / (len(cs) * RC_CCG[2])
        RC_CCGS[(c, int(r.ss))] = h / r.baseline if r.baseline > 0 else h
_best = recheck[recheck.post_rank == 1].set_index("cs")
_s0 = recheck[recheck.step0_rank == 1].set_index("cs")
print("pairing re-check: CS where the best post-window pause (3–8 ms) is not step0's rank-1 SS:")
for c in _best.index:
    if c in _s0.index and _best.at[c, "ss"] != _s0.at[c, "ss"]:
        print(f"  CS {c}: step0 #1 {_s0.at[c, 'ss']} (pr0-8 {_s0.at[c, 'pr_0_8']:.2f}, pr3-8 {_s0.at[c, 'pr_3_8']:.2f})"
              f"  vs best post {_best.at[c, 'ss']} (pr0-8 {_best.at[c, 'pr_0_8']:.2f}, pr3-8 {_best.at[c, 'pr_3_8']:.2f})")

# %% 4. SNIPPET READER
SOS = butter(HP_ORDER, HP_HZ, "highpass", fs=FS, output="sos")
S0, S1, SP = int(round(WIN[0] * FS)), int(round(WIN[1] * FS)), int(round(PAD * FS))
OFF = np.arange(S0 - SP, S1 + SP)                                         # read offsets (samples)
CROP = slice(SP, SP + (S1 - S0))
T_MS = np.arange(S0, S1) / FS * 1e3
N_RAW = raw.shape[0]


def car_snippets(samples):
    """(n,) spike samples -> (n, L, 384) float32 µV, mean-removed per channel, common median referenced."""
    s = np.asarray(samples, dtype=np.int64)
    s = s[(s + OFF[0] >= 0) & (s + OFF[-1] < N_RAW)]
    x = raw[s[:, None] + OFF[None, :], :N_AP].astype(np.float32)
    x = x[:, :, chan_map]                                                 # sorter channel order
    x -= x.mean(axis=1, keepdims=True)
    x -= np.median(x, axis=2, keepdims=True)
    return x * UV


def run_batches(samples, per_batch):
    """Stream CAR'd snippets in time order; per_batch(x, idx) gets each batch and the indices into `samples`."""
    order = np.argsort(samples)
    for b in range(0, len(order), BATCH):
        ix = order[b:b + BATCH]
        per_batch(car_snippets(samples[ix]), ix)


def hp(x, axis):
    return sosfiltfilt(SOS, x, axis=axis)


def mean_waveform(samples):
    """Mean over snippets, all channels, filtered and cropped: (n_samp, 384) µV."""
    acc = np.zeros((len(OFF), N_AP))
    n = [0]

    def f(x, ix):
        acc[:] += x.sum(0); n[0] += len(x)
    run_batches(np.asarray(samples), f)
    return hp(acc / max(n[0], 1), 0)[CROP], n[0]


def stored_channels(m, peak):
    y = chan_pos[:, 1]
    ch = np.flatnonzero(np.abs(y - y[peak]) <= MEAN_CH_UM)
    return ch[np.lexsort((chan_pos[ch, 0], chan_pos[ch, 1]))]


def top_channels(m, n, extra=()):
    ptp = np.ptp(m, 0)
    ch = list(np.argsort(ptp)[::-1][:n])
    for e in extra:
        if e not in ch:
            ch.append(int(e))
    return sorted(ch, key=lambda c: (chan_pos[c, 1], chan_pos[c, 0]))


def pick(arr, n):
    return arr if len(arr) <= n else np.sort(RNG.choice(arr, n, replace=False))


# %% 5. MEAN WAVEFORMS (every selected unit; its own stable span)
_t = time.time()
MEAN = {}
for i, u in enumerate(ALL_UNITS):
    s0, s1 = span(u)
    smp = SAMP[u][(SAMP[u] >= s0 * FS) & (SAMP[u] < s1 * FS)]
    m, n = mean_waveform(pick(smp, N_MEAN))
    MEAN[u] = dict(mean=m, n=n, peak=int(np.argmax(np.ptp(m, 0))))
    if i % 20 == 0:
        print(f"  mean {i + 1}/{len(ALL_UNITS)} ({time.time() - _t:.0f} s)")
print(f"mean waveforms: {len(MEAN)} units in {time.time() - _t:.0f} s")

# %% 6. CS EVENTS: single events, spikelet detection, conditional means
_t = time.time()
CSR = {}
pre = (T_MS >= -2) & (T_MS <= -0.3)
spk_lo, spk_hi = int(round((SPK_WIN[0] - WIN[0]) * FS)), int(round((SPK_WIN[1] - WIN[0]) * FS))
for c in CS_IDS:
    w = span(c)
    cs_t = times_in(c, w)
    cs_s = np.round(cs_t * FS).astype(np.int64)
    pk = MEAN[c]["peak"]
    caps = cap[cap["cs"] == c]
    # which capturing clusters fired in each CS event
    member = {}
    for k in caps["unit"]:
        sk = SAMP[k] / FS
        member[int(k)] = np.searchsorted(sk, cs_t + CAPTURE_WIN[1]) > np.searchsorted(sk, cs_t + CAPTURE_WIN[0])
    # one pass over all events: peak-channel trace per event + conditional sums
    cond_ch = {int(k): top_channels(MEAN[c]["mean"], COND_N_CH, extra=[MEAN[int(k)]["peak"]]) for k in caps["unit"]}
    trace = np.zeros((len(cs_s), len(OFF)), np.float32)
    sums = {k: [np.zeros((len(OFF), len(ch))), np.zeros((len(OFF), len(ch)))] for k, ch in cond_ch.items()}

    def f(x, ix):
        trace[ix] = x[:, :, pk]
        for k, ch in cond_ch.items():
            mk = member[k][ix]
            sums[k][0] += x[mk][:, :, ch].sum(0); sums[k][1] += x[~mk][:, :, ch].sum(0)
    run_batches(cs_s, f)
    tr = hp(trace, 1)[:, CROP]
    # polarity of the main spike on this channel (molecular-layer CS can be positive-going); spikelets = same polarity
    mtr = tr.mean(0); z0 = np.flatnonzero(np.abs(T_MS) <= 0.5)
    pol = float(np.sign(mtr[z0[np.argmax(np.abs(mtr[z0]))]])) or -1.0
    res = tr[:, pre] - mtr[pre]                                          # residual noise (fixed waveform removed)
    sd = 1.4826 * np.median(np.abs(res - np.median(res)))

    def detect(trs):
        cnt, lg, am, eo = [], [], [], []
        for e in range(len(trs)):
            seg = -pol * trs[e, spk_lo:spk_hi]
            ip, pr = find_peaks(seg, height=SPK_THR * sd, prominence=SPK_PROM * sd, distance=max(1, int(SPK_MIN_SEP * FS)))
            cnt.append(len(ip))
            lg.extend((ip + spk_lo) / FS * 1e3 + WIN[0] * 1e3); am.extend(-pol * pr["peak_heights"]); eo.extend([e] * len(ip))
        return map(np.asarray, (cnt, lg, am, eo))
    counts, lags, amps, ev_of = detect(tr)
    # null: same channel, same detector, random times away from this CS's spikes
    allc = SAMP[c] / FS
    rt = np.sort(RNG.uniform(w[0] + 0.05, w[1] - 0.05, N_NULL * 3))
    j = np.searchsorted(allc, rt)
    far = (np.abs(rt - allc[np.clip(j - 1, 0, len(allc) - 1)]) > NULL_EXCL) & (np.abs(allc[np.clip(j, 0, len(allc) - 1)] - rt) > NULL_EXCL)
    rs = np.round(rt[far][:N_NULL] * FS).astype(np.int64)
    ntrace = np.zeros((len(rs), len(OFF)), np.float32)

    def fn(x, ix):
        ntrace[ix] = x[:, :, pk]
    run_batches(rs, fn)
    ntr = hp(ntrace, 1)[:, CROP]
    null_counts, null_lags, _, _ = detect(ntr)
    # raw spikelet presence near each capturing cluster's lag, when the cluster fired vs not
    agree = []
    for k, row in caps.set_index("unit").iterrows():
        lo, hi = row.lag_q25_ms - 0.3, row.lag_q75_ms + 0.3
        has = np.zeros(len(tr), bool)
        has[np.unique(ev_of[(lags >= lo) & (lags <= hi)])] = True
        mk = member[int(k)]
        agree.append(dict(cs=c, unit=int(k), lag_lo_ms=lo, lag_hi_ms=hi, n_fired=int(mk.sum()),
                          raw_spikelet_if_fired=has[mk].mean() if mk.any() else np.nan,
                          raw_spikelet_if_not=has[~mk].mean() if (~mk).any() else np.nan))
    cond = []
    for k, ch in cond_ch.items():
        nw, nwo = int(member[k].sum()), int((~member[k]).sum())
        if nw < N_COND_MIN or nwo < N_COND_MIN:
            continue
        cond.append(dict(unit=k, ch=ch, n_with=nw, n_without=nwo,
                         with_=hp(sums[k][0] / nw, 0)[CROP], without=hp(sums[k][1] / nwo, 0)[CROP]))
    sel = pick(np.arange(len(tr)), N_SINGLE)
    z = (T_MS >= -0.3) & (T_MS <= 0.3)
    CSR[c] = dict(peak=pk, sd=sd, pol=pol, null_counts=null_counts, null_lags=null_lags,
                  n_events=len(tr), counts=counts, lags=lags, amps=amps, agree=agree, cond=cond,
                  all_traces=tr.astype(np.float16), ev_t=cs_t, main_amp=pol * (pol * tr[:, z]).max(1), ev_of=ev_of,
                  members={k: v.copy() for k, v in member.items()},
                  single=tr[sel], single_members={k: member[k][sel] for k in member},
                  single_spk=[lags[ev_of == e] for e in sel])
    print(f"  CS {c}: {len(tr)} events, polarity {pol:+.0f}, noise SD {sd:.1f} µV, spikelets/event {counts.mean():.2f} "
          f"(null {null_counts.mean():.2f}; 0 spikelets in {np.mean(counts == 0):.0%}), {len(cond)} conditional means "
          f"({time.time() - _t:.0f} s)")

# %% 7. CAPTURING CLUSTERS: their CS-locked spikes vs the others
_t = time.time()
LOCKED = []
for k in CAP_UNITS:
    s0, s1 = span(k)
    sk = SAMP[k][(SAMP[k] >= s0 * FS) & (SAMP[k] < s1 * FS)] / FS
    lk = np.zeros(len(sk), bool)
    for c in cap.loc[cap["unit"] == k, "cs"]:
        cs_t = times_in(c, (s0, s1))
        j = np.searchsorted(cs_t, sk - CAPTURE_WIN[1])
        lk |= (j < len(cs_t)) & (cs_t[np.clip(j, 0, len(cs_t) - 1)] <= sk - CAPTURE_WIN[0])
    if lk.sum() < N_COND_MIN or (~lk).sum() < N_COND_MIN:
        continue
    ml, nl = mean_waveform(np.round(pick(sk[lk], N_MEAN) * FS).astype(np.int64))
    mo, no = mean_waveform(np.round(pick(sk[~lk], N_MEAN) * FS).astype(np.int64))
    ch = top_channels(MEAN[k]["mean"], COND_N_CH)
    LOCKED.append(dict(unit=k, ch=ch, n_locked=int(lk.sum()), n_other=int((~lk).sum()), frac_locked=lk.mean(),
                       locked=ml[:, ch], other=mo[:, ch]))
print(f"locked vs other means: {len(LOCKED)} clusters in {time.time() - _t:.0f} s")

# %% 7b. CS–CS PAIRS FROM step0 (duplicate evidence): raw means of coincident / a-only / b-only events
# All three in a's time frame (b-only events shifted by the pair's peak lag), on the top channels of both units, so a
# re-detection shows the same waveform in all three, and a split CS shows complementary pieces.
_t = time.time()
ev0 = pd.read_csv(STEP0_DIR / "tables" / "cs_cs_evidence.csv")
COINC_D = float(p0.get("COINC_DELTA", 0.5e-3))
DUPS = []
for r in ev0.itertuples():
    a, b, L = int(r.cs_a), int(r.cs_b), float(r.peak_lag_ms) * 1e-3
    if a not in MEAN or b not in MEAN:
        continue
    w = (r.win_start, r.win_end)
    ta, tb = times_in(a, w), times_in(b, w)
    ia = np.searchsorted(tb, ta + L - COINC_D); co_a = (ia < len(tb)) & (tb[np.clip(ia, 0, len(tb) - 1)] <= ta + L + COINC_D)
    ib = np.searchsorted(ta, tb - L - COINC_D); co_b = (ib < len(ta)) & (ta[np.clip(ib, 0, len(ta) - 1)] <= tb - L + COINC_D)
    sets = {"coincident": ta[co_a], f"{a}-only": ta[~co_a], f"{b}-only": tb[~co_b] - L}
    ch = sorted(set(top_channels(MEAN[a]["mean"], 6)) | set(top_channels(MEAN[b]["mean"], 6)),
                key=lambda q: (chan_pos[q, 1], chan_pos[q, 0]))
    means, ns = {}, {}
    for name, tt in sets.items():
        if len(tt) < N_COND_MIN:
            continue
        m, n = mean_waveform(np.round(pick(tt, N_MEAN) * FS).astype(np.int64))
        means[name], ns[name] = m[:, ch], n
    # whole-mean similarity on the union channels: a vs b shifted by the lag, full window and late part only (> 1.5 ms)
    shift = int(round(L * FS))
    mb = np.roll(MEAN[b]["mean"], shift, axis=0)[:, ch]; ma = MEAN[a]["mean"][:, ch]
    late = T_MS > 1.5
    cosf = lambda x, y: float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum() + 1e-12))
    DUPS.append(dict(cs_a=a, cs_b=b, lag_ms=L * 1e3, status=r.evidence_status, dup_candidate=bool(r.dup_candidate),
                     coinc=r.coinc_frac_peak, ch=ch, means=means, n=ns,
                     cos_full=cosf(ma, mb), cos_late=cosf(ma[late], mb[late])))
print(f"CS–CS raw comparisons: {len(DUPS)} pairs in {time.time() - _t:.0f} s")

# %% 8. SAVE: tables, npz, browser data (raw_data.js, loaded by browser.html from ../raw_waveforms/)
spk_rows = [dict(cs=c, n_events=r["n_events"], polarity=r["pol"], noise_sd_uv=r["sd"], spikelets_mean=r["counts"].mean(),
                 null_mean=r["null_counts"].mean(), excess_mean=r["counts"].mean() - r["null_counts"].mean(),
                 **{f"frac_{i}": np.mean(r["counts"] == i) for i in range(4)},
                 frac_4plus=np.mean(r["counts"] >= 4),
                 lag_mode_ms=(np.bincount(np.round(r["lags"] * 10).astype(int)).argmax() / 10) if len(r["lags"]) else np.nan)
            for c, r in CSR.items()]
pd.DataFrame(spk_rows).round(4).to_csv(OUT_DIR / "cs_raw_spikelets.csv", index=False)
pd.DataFrame([a for r in CSR.values() for a in r["agree"]]).round(4).to_csv(OUT_DIR / "cs_capture_vs_raw_spikelet.csv",
                                                                           index=False)
np.savez_compressed(OUT_DIR / "raw_means.npz", t_ms=T_MS, units=np.array(ALL_UNITS),
                    means=np.stack([MEAN[u]["mean"] for u in ALL_UNITS]).astype(np.float32),
                    n=np.array([MEAN[u]["n"] for u in ALL_UNITS]), chan_pos=chan_pos)
# every CS event's peak-channel trace + per-event info (for raw_waveforms_figures.py)
_ev = {}
for c, r in CSR.items():
    _ev[f"tr_{c}"] = r["all_traces"]; _ev[f"t_{c}"] = r["ev_t"]; _ev[f"amp_{c}"] = r["main_amp"]
    _ev[f"cnt_{c}"] = r["counts"]; _ev[f"lag_{c}"] = r["lags"]; _ev[f"lagev_{c}"] = r["ev_of"]
    for k, v in r["members"].items():
        _ev[f"mem_{c}_{k}"] = v
np.savez_compressed(OUT_DIR / "cs_events.npz", t_ms=T_MS, cs_ids=np.array(list(CSR)),
                    peaks=np.array([CSR[c]["peak"] for c in CSR]), sd=np.array([CSR[c]["sd"] for c in CSR]), pol=np.array([CSR[c]["pol"] for c in CSR]),
                    **{f"null_{c}": CSR[c]["null_counts"] for c in CSR}, **{f"nulllag_{c}": CSR[c]["null_lags"] for c in CSR}, **_ev)
import pickle
with open(OUT_DIR / "raw_extra.pkl", "wb") as fh:                         # conditional / locked / dup means for figures
    pickle.dump(dict(cond={c: r["cond"] for c, r in CSR.items()}, locked=LOCKED, dups=DUPS, rc_ccgs=RC_CCGS,
                     rc_edges=_edges, cap=cap, recheck=recheck), fh)


def r1(a):
    return np.round(np.asarray(a, float), 1).tolist()


ystep = np.unique(chan_pos[:, 1])
RAWD = dict(meta=dict(mouse=MOUSE, session=SESSION_R, step0_dir=STEP0_DIR.name, fs=FS, t_ms=np.round(T_MS, 4).tolist(), win_ms=[WIN[0] * 1e3, WIN[1] * 1e3],
                      uv_per_bit=UV, hp_hz=HP_HZ, car="median", n_mean=N_MEAN, spk_thr=SPK_THR, spk_prom=SPK_PROM,
                      spk_win_ms=[SPK_WIN[0] * 1e3, SPK_WIN[1] * 1e3], capture_win_ms=[CAPTURE_WIN[0] * 1e3,
                                                                                         CAPTURE_WIN[1] * 1e3],
                      created=time.strftime("%Y-%m-%d %H:%M"), raw_bin=RAW_BIN.name),
            units={}, cs={}, locked={}, capture=cap.round(4).to_dict("records"))
for u in ALL_UNITS:
    m = MEAN[u]["mean"]; ch = stored_channels(m, MEAN[u]["peak"])
    prof = np.ptp(m, 0)
    py = [float(y) for y in ystep]
    pv = [float(prof[chan_pos[:, 1] == y].max()) for y in ystep]
    RAWD["units"][str(u)] = dict(n=MEAN[u]["n"], peak=MEAN[u]["peak"], ch=[int(c) for c in ch],
                                 x=r1(chan_pos[ch, 0]), y=r1(chan_pos[ch, 1]), mean=[r1(m[:, c]) for c in ch],
                                 prof_y=py, prof_ptp=r1(pv))
for c, r in CSR.items():
    hc = np.bincount(np.minimum(r["counts"], 6), minlength=7)
    hl = np.histogram(r["lags"], np.arange(SPK_WIN[0] * 1e3, SPK_WIN[1] * 1e3 + 1e-9, 0.1))[0]
    RAWD["cs"][str(c)] = dict(
        peak=r["peak"], noise_sd=round(float(r["sd"]), 2), n_events=r["n_events"], count_hist=hc.tolist(), pol=r["pol"],
        null_mean=round(float(r["null_counts"].mean()), 3), null_count_hist=np.bincount(np.minimum(r["null_counts"], 6), minlength=7).tolist(),
        null_lag_hist=np.histogram(r["null_lags"], np.arange(SPK_WIN[0] * 1e3, SPK_WIN[1] * 1e3 + 1e-9, 0.1))[0].tolist(), n_null=len(r["null_counts"]),
        lag_hist=hl.tolist(), lag_edges=[SPK_WIN[0] * 1e3, 0.1], agree=pd.DataFrame(r["agree"]).round(3).to_dict("records"),
        cond=[dict(unit=d["unit"], ch=[int(x) for x in d["ch"]], y=r1(chan_pos[d["ch"], 1]), x=r1(chan_pos[d["ch"], 0]),
                   n_with=d["n_with"], n_without=d["n_without"],
                   with_=[r1(d["with_"][:, i]) for i in range(len(d["ch"]))],
                   without=[r1(d["without"][:, i]) for i in range(len(d["ch"]))]) for d in r["cond"]],
        single=[r1(t) for t in r["single"]], single_spk=[r1(x) for x in r["single_spk"]],
        single_members={str(k): v.astype(int).tolist() for k, v in r["single_members"].items()})
for d in LOCKED:
    RAWD["locked"][str(d["unit"])] = dict(ch=[int(x) for x in d["ch"]], y=r1(chan_pos[d["ch"], 1]),
                                         x=r1(chan_pos[d["ch"], 0]), n_locked=d["n_locked"], n_other=d["n_other"],
                                         frac_locked=round(float(d["frac_locked"]), 4),
                                         locked=[r1(d["locked"][:, i]) for i in range(len(d["ch"]))],
                                         other=[r1(d["other"][:, i]) for i in range(len(d["ch"]))])
RAWD["recheck"] = recheck.round(4).to_dict("records")
RAWD["rc_ccg"] = dict(lo_ms=RC_CCG[0] * 1e3, bin_ms=RC_CCG[2] * 1e3, base_ms=[PR_BASE[0] * 1e3, PR_BASE[1] * 1e3],
                      wins={k: [v[0] * 1e3, v[1] * 1e3] for k, v in PR_WINS.items()},
                      y={f"{c}_{k}": np.round(v, 3).tolist() for (c, k), v in RC_CCGS.items()})
RAWD["dups"] = [dict(cs_a=d["cs_a"], cs_b=d["cs_b"], lag_ms=round(d["lag_ms"], 3), status=d["status"],
                     dup_candidate=d["dup_candidate"], coinc=d["coinc"], cos_full=round(d["cos_full"], 3),
                     cos_late=round(d["cos_late"], 3), ch=[int(x) for x in d["ch"]], y=r1(chan_pos[d["ch"], 1]),
                     x=r1(chan_pos[d["ch"], 0]), n=d["n"],
                     means={k: [r1(v[:, i]) for i in range(len(d["ch"]))] for k, v in d["means"].items()})
                for d in DUPS]
# per CS: spikelet count and main-spike amplitude over the recording (60 s bins) + count vs amplitude
for c, r in CSR.items():
    tb = np.floor(r["ev_t"] / 60).astype(int)
    ub = np.unique(tb)
    RAWD["cs"][str(c)].update(
        time_min=[float(b) for b in ub], cnt_t=[round(float(r["counts"][tb == b].mean()), 3) for b in ub],
        amp_t=[round(float(np.median(r["main_amp"][tb == b])), 1) for b in ub],
        n_t=[int((tb == b).sum()) for b in ub],
        amp_by_cnt=[round(float(np.median(r["main_amp"][np.minimum(r["counts"], 6) == i])), 1)
                    if (np.minimum(r["counts"], 6) == i).sum() >= 5 else None for i in range(7)],
        spikelets_mean=round(float(r["counts"].mean()), 3))


def clean(o):
    """JSON-safe: NaN/inf -> None, numpy scalars -> python."""
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    return o


_js = "window.CS_RAW = " + json.dumps(clean(RAWD), separators=(",", ":"), allow_nan=False) + ";\n"
(OUT_DIR / "raw_data.js").write_text(_js, encoding="utf-8")
(OUT_DIR / "params_used.json").write_text(json.dumps(params_used, indent=1, default=str))
print(f"saved {OUT_DIR} (raw_data.js {len(_js) / 1e6:.1f} MB)")
print(pd.DataFrame(spk_rows).round(2).to_string(index=False))
