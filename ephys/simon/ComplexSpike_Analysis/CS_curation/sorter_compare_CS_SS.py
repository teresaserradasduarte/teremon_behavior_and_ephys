# -*- coding: utf-8 -*-
"""
sorter_compare_CS_SS.py — Does CS / SS detection (as classified by C4) differ between the two sorters?

Population: extra-good units in CB (brain_region == "CB"), C4 PkC_cs (CS) and PkC_ss (SS).
Sessions: every CB sorter folder that carries C4 output (cluster_predicted_cell_type.tsv). Each session
          has exactly ONE such folder, so the sorter is a between-session factor (no session was
          labelled with both sorters). Sorter type: "ks4" (catGT/kilosort4) or "ibl" (ibl_sorter_results*).
Window:   bhv.win_selected (behavior_fundamentals.mat, Harp clock) mapped to the imec clock with the
          ephys→Harp linear fit recovered from neurons_session.mat (st_raw_all / st_all), as in step0.
          Every spike-train metric uses only spikes inside this window.

Outputs (local only): C:\\Users\\simon\\Desktop\\Data\\ephys_and_behavior\\CS_curation\\_sorter_comparison\\
  sessions/<mouse>_<S#>/units.csv, cs_ss_pairs.csv, cs_cs_pairs.csv   (per-session cache)
  session_table.csv, units_all.csv, session_summary.csv, stats_session_level.csv, figs/*.png

Statistics: session-level (session = the unit of replication for a session-level factor).
  - Mann–Whitney U on per-session values (KS4 vs IBL).
  - Within-mouse permutation test: sorter labels are shuffled WITHIN each mouse; statistic = mean over
    mice of (median KS4 sessions − median IBL sessions). Controls for animal differences.

GOLDEN RULE: Dropbox is READ-ONLY. Data are opened read-only; outputs refuse any Dropbox path.
"""

# %% 1. PARAMETERS + PATHS
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scipy.io as sio
import h5py
from scipy.stats import mannwhitneyu

DROPBOX_ROOT = Path(r"D:\Learning Lab Dropbox")
TD_ROOT = DROPBOX_ROOT / r"Learning Lab Team Folder\Patlab protocols\Data\TD"
EPHYS_ROOT = TD_ROOT / r"ephys_data\20230801_ChocolateGroup"
MAT_ROOT = TD_ROOT / r"ephys_and_behavior\mat_files\20230801_ChocolateGroup"
LOCAL_ROOT = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior")
CODE_DIR = Path(r"C:\Users\simon\Desktop\Code\EPhys\ComplexSpike_Analysis")
OUT_DIR = LOCAL_ROOT / "CS_curation" / "_sorter_comparison"
FORCE = False                     # True: recompute sessions even if cached

MICE = ["1_CoteDor", "2_Lindt", "3_Toblerone", "4_Milka", "5_FerreroRocher"]
REGION = "CB"
CS_LABEL, SS_LABEL = "PkC_cs", "PkC_ss"

# spike-train metrics
STAB_BIN = 60.0                   # s; presence ratio / rate stability bins
DOUBLET_ISI = 5e-3                # CS ISIs below this = multiple detections of one CS (spikelets)
REFRAC_ISI = 1e-3                 # refractory-violation threshold
# CS→SS pause (same definitions as step0_explore_session.py)
COARSE_BIN, COARSE_WIN = 1e-3, 50e-3
FINE_BIN, FINE_WIN = 0.1e-3, 10e-3
PAUSE_WIN, MIRROR_WIN, BASELINE_WIN = (0.0, 8e-3), (-8e-3, 0.0), (-50e-3, -10e-3)
FINE_BASELINE_WIN, DIP_HALF = (-10e-3, -2e-3), 0.5e-3
JITTER_WIN, JITTER_HALF, N_SURR, RNG_SEED = (1e-3, 8e-3), 25e-3, 500, 0
CAND_PAUSE_RATIO_MAX, CAND_JITTER_P_MAX, CAND_ASYM_MIN = 0.7, 0.05, 0.2
N_JITTER_TOP = 3                  # jitter test only for the N lowest-pause-ratio SS per CS (speed)
EDGE_PAD = 0.1
# CS–CS duplicates
COINC_DELTA = 0.3e-3
DUP_COINC_MIN = 0.2               # step0 candidate flag
DUP_COINC_STRONG = 0.5            # "likely duplicate"
# waveform
FOOTPRINT_FRAC = 0.3              # channels with ptp > this × peak ptp
N_PERM = 10000

SORTER_COLOR = {"ks4": "#2a78d6", "ibl": "#eb6834"}   # dataviz reference palette slots 1–2
SORTER_NAME = {"ks4": "Kilosort4", "ibl": "IBL sorter"}
MOUSE_MARKER = dict(zip(MICE, ["o", "s", "^", "D", "v"]))


def safe_out_dir(p: Path) -> Path:
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
SES_DIR = safe_out_dir(OUT_DIR / "sessions")


# %% 2. SESSION TABLE (discovered: the CB sorter folder is the one with C4 output)
def discover_sessions():
    rows = []
    for mouse in MICE:
        for sess_dir in sorted((EPHYS_ROOT / mouse).iterdir()):
            if not sess_dir.is_dir():
                continue
            m = re.search(r"S(\d)", sess_dir.name.split("_", 1)[1])
            if m is None or not sess_dir.name.split("_")[1].startswith(mouse.split("_")[1][:4]):
                continue  # e.g. a Toblerone folder misplaced under CoteDor
            for imec_dir in sorted(p for p in sess_dir.glob("*_imec?") if p.is_dir()):
                for sd in [imec_dir / "catGT" / "kilosort4", *sorted(imec_dir.glob("ibl_sorter_results*"))]:
                    if sd.is_dir() and any(f.name.lower() == "cluster_predicted_cell_type.tsv" for f in sd.iterdir()):
                        rows.append(dict(mouse=mouse, S=f"S{m.group(1)}", R=f"R{m.group(1)}",
                                         imec=imec_dir.name.rsplit("_", 1)[1],
                                         sorter="ks4" if sd.name == "kilosort4" else "ibl",
                                         sorter_folder=sd.name, sorter_dir=str(sd)))
    df = pd.DataFrame(rows)
    dup = df.duplicated(["mouse", "S"], keep=False)
    assert not dup.any(), f"more than one C4 sorter folder for a session:\n{df[dup]}"
    return df


t0 = time.time()
sessions = discover_sessions()
sessions["session_index"] = sessions["S"].str[1:].astype(int)
print(f"{len(sessions)} sessions with C4 output ({time.time() - t0:.0f} s)")
print(pd.crosstab(sessions["mouse"], sessions["sorter"], margins=True))
print(sessions[["mouse", "S", "imec", "sorter", "sorter_folder"]].to_string(index=False))
sessions.to_csv(OUT_DIR / "session_table.csv", index=False)


# %% 3. FUNCTIONS
def _pair_diffs(ref, tgt, win, exclude_self=False, chunk=200_000):
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
        if exclude_self:
            keep = j != (ri + s)
            ri, j = ri[keep], j[keep]
        out.append(tgt[j] - r[ri])
    return np.concatenate(out) if out else np.empty(0)


def ccg(ref, tgt, bin_s, win_s):
    nb = int(round(win_s / bin_s))
    d = _pair_diffs(ref, tgt, nb * bin_s)
    idx = np.floor(d / bin_s + 1e-9).astype(np.int64) + nb
    idx = idx[(idx >= 0) & (idx < 2 * nb)]
    counts = np.bincount(idx, minlength=2 * nb).astype(float)
    left = (np.arange(2 * nb) - nb) * bin_s
    return counts / (max(len(ref), 1) * bin_s), left, counts


def _win_mask(left, bin_s, win):
    return (left >= win[0] - 1e-12) & (left + bin_s <= win[1] + 1e-12)


def pause_metrics(cs, ss):
    rc, lc, _ = ccg(cs, ss, COARSE_BIN, COARSE_WIN)
    rf, lf, _ = ccg(cs, ss, FINE_BIN, FINE_WIN)
    base = rc[_win_mask(lc, COARSE_BIN, BASELINE_WIN)].mean()
    m = dict(baseline_rate=base)
    m["pause_ratio"] = rc[_win_mask(lc, COARSE_BIN, PAUSE_WIN)].mean() / base if base > 0 else np.nan
    m["mirror_ratio"] = rc[_win_mask(lc, COARSE_BIN, MIRROR_WIN)].mean() / base if base > 0 else np.nan
    m["asym"] = m["mirror_ratio"] - m["pause_ratio"]
    pos = lc >= -1e-12
    ok = np.nonzero(rc[pos] >= 0.5 * base)[0]
    m["pause_dur50_ms"] = lc[pos][ok[0]] * 1e3 if len(ok) else np.nan
    fb = rf[_win_mask(lf, FINE_BIN, FINE_BASELINE_WIN)].mean()
    dn = rf[_win_mask(lf, FINE_BIN, (-DIP_HALF, 0.0))].mean()
    dp = rf[_win_mask(lf, FINE_BIN, (0.0, DIP_HALF))].mean()
    late = rc[_win_mask(lc, COARSE_BIN, (1e-3, PAUSE_WIN[1]))].mean()
    m["collision_dip"] = bool(fb > 0 and dn < 0.5 * fb and dp < 0.5 * fb and late >= 0.8 * base)
    m["zero_peak_ratio"] = rf[_win_mask(lf, FINE_BIN, (-0.5e-3, 0.5e-3))].mean() / fb if fb > 0 else np.nan
    return m


def jitter_pause_test(cs, ss):
    rng = np.random.default_rng(RNG_SEED)
    obs = int((np.searchsorted(ss, cs + JITTER_WIN[1]) - np.searchsorted(ss, cs + JITTER_WIN[0])).sum())
    flat = (cs[None, :] + rng.uniform(-JITTER_HALF, JITTER_HALF, (N_SURR, len(cs)))).ravel()
    c = np.searchsorted(ss, flat + JITTER_WIN[1]) - np.searchsorted(ss, flat + JITTER_WIN[0])
    surr = c.reshape(N_SURR, len(cs)).sum(1)
    return dict(jitter_p=(1 + np.sum(surr <= obs)) / (1 + N_SURR),
                jitter_deficit=1 - obs / surr.mean() if surr.mean() > 0 else np.nan)


def coincident_mask(ref, tgt, lag, delta=COINC_DELTA):
    return np.searchsorted(tgt, ref + lag + delta, "right") > np.searchsorted(tgt, ref + lag - delta, "left")


def cs_cs_coinc(a, b):
    _, lf, cf = ccg(a, b, FINE_BIN, FINE_WIN)
    lag = (lf + FINE_BIN / 2)[int(np.argmax(cf))]
    small, big, sl = (a, b, lag) if len(a) <= len(b) else (b, a, -lag)
    return dict(peak_lag_ms=lag * 1e3, coinc_frac_peak=float(coincident_mask(small, big, sl).mean()),
                coinc_frac_zero=float(coincident_mask(small, big, 0.0).mean()))


def train_metrics(st, amp, T0, T1):
    """Window-restricted spike-train metrics for one unit (st sorted, imec s)."""
    n, T = len(st), T1 - T0
    m = dict(n_spikes=n, fr=n / T)
    isi = np.diff(st)
    if len(isi) >= 2:
        m["cv_isi"] = isi.std() / isi.mean()
        m["cv2"] = float(np.mean(2 * np.abs(np.diff(isi)) / (isi[1:] + isi[:-1])))
        m["frac_isi_lt1ms"] = float(np.mean(isi < REFRAC_ISI))
        m["frac_isi_lt5ms"] = float(np.mean(isi < DOUBLET_ISI))
        m["median_isi_ms"] = float(np.median(isi) * 1e3)
    edges = np.arange(T0, T1 + 1e-9, STAB_BIN)
    if len(edges) >= 3:
        c = np.histogram(st, edges)[0] / STAB_BIN
        m["presence_ratio"] = float(np.mean(c > 0))
        m["fr_cv_60s"] = float(c.std() / c.mean()) if c.mean() > 0 else np.nan
        k = max(1, len(c) // 5)
        m["fr_last_first"] = float(c[-k:].mean() / c[:k].mean()) if c[:k].mean() > 0 else np.nan
    if n >= 20:
        q = np.percentile(amp, [25, 50, 75])
        m["amp_rcv"] = float((q[2] - q[0]) / q[1]) if q[1] > 0 else np.nan   # scale-free amplitude spread
        k = max(5, n // 5)
        m["amp_last_first"] = float(np.median(amp[-k:]) / np.median(amp[:k]))  # amplitude drift in window
    return m


def ephys_to_harp_fit(path, n_check=3):
    fits = []
    with h5py.File(path, "r") as h:
        raw, harp = h["neurons"]["st_raw_all"], h["neurons"]["st_all"]
        for i in np.linspace(0, raw.shape[0] - 1, n_check).astype(int):
            x, y = h[raw[i, 0]][()].ravel(), h[harp[i, 0]][()].ravel()
            if len(x) < 2 or len(x) != len(y):
                continue
            p = np.polyfit(x, y, 1)
            fits.append((p[0], p[1], float(np.abs(y - np.polyval(p, x)).max())))
    fits = np.array(fits)
    assert len(fits), "no usable unit for the clock fit"
    assert np.ptp(fits[:, 0]) < 1e-9 and np.ptp(fits[:, 1]) < 1e-4 and fits[:, 2].max() < 1e-4, f"fits disagree {fits}"
    return float(fits[0, 0]), float(fits[0, 1])


def read_win_selected(path):
    try:
        b = sio.loadmat(path, squeeze_me=True, struct_as_record=False, variable_names=["bhv"])["bhv"]
        w = b.win_selected
        return float(w.select_start_time), float(w.time_sess_end), np.atleast_1d(w.trials_vec)
    except NotImplementedError:
        with h5py.File(path, "r") as h:
            w = h["bhv"]["win_selected"]
            return (float(w["select_start_time"][()].ravel()[0]), float(w["time_sess_end"][()].ravel()[0]),
                    w["trials_vec"][()].ravel())


def read_label(files, name, col):
    p = files.get(name.lower())
    if p is None:
        return None
    df = pd.read_csv(p, sep="\t")
    return df.set_index(df.columns[0]).iloc[:, 0].rename(col)


# %% 4. PER-SESSION PROCESSING (cached)
def process_session(row):
    sd = Path(row.sorter_dir)
    files = {p.name.lower(): p for p in sd.iterdir() if p.is_file()}
    params = {}
    for line in files["params.py"].read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                params[k.strip()] = eval(v.strip(), {}, {})
            except Exception:
                params[k.strip()] = v.strip()
    SR = float(params["sample_rate"])

    # labels
    lab = pd.DataFrame(index=pd.Index(pd.read_csv(files["cluster_group.tsv"], sep="\t").iloc[:, 0], name="phyID"))
    lab = lab.join(read_label(files, "cluster_group.tsv", "group"))
    eg = read_label(files, "cluster_extragood.tsv", "eg")
    lab["extragood"] = lab.index.isin(eg[eg.astype(str).str.lower() == "y"].index)
    for f, c in [("cluster_predicted_cell_type.tsv", "c4_label"), ("cluster_confidence_ratio.tsv", "c4_conf"),
                 ("cluster_pred_probability.tsv", "c4_prob"), ("cluster_model_votes.tsv", "c4_votes")]:
        s = read_label(files, f, c)
        lab = lab.join(s, how="outer") if s is not None else lab
    v = lab["c4_votes"].dropna().astype(str).str.split("/", expand=True).astype(float)
    lab["c4_vote_ratio"] = (v[0] / v[1]).reindex(lab.index)
    reg = pd.read_csv(files["cluster_channel_region_map.txt"], sep="\t").set_index("cluster_id")
    lab["region"] = reg["brain_region"].reindex(lab.index)
    lab["extragood"] = lab["extragood"].fillna(False).astype(bool)
    pop = lab[lab["extragood"] & (lab["region"] == REGION)].copy()
    pop["c4_label"] = pop["c4_label"].fillna("unlabelled")

    # spikes
    st_all = np.asarray(np.load(files["spike_times.npy"], mmap_mode="r"), dtype=np.float64).ravel() / SR
    clu = np.asarray(np.load(files["spike_clusters.npy"], mmap_mode="r")).ravel().astype(np.int64)
    amps = np.asarray(np.load(files["amplitudes.npy"], mmap_mode="r")).ravel()
    T_REC = st_all[-1]

    # window
    bhv_path = MAT_ROOT / row.mouse / row.R / "behavior_fundamentals.mat"
    h0, h1, trials = read_win_selected(bhv_path)
    a, b = ephys_to_harp_fit(files["neurons_session.mat"])
    T0, T1 = (h0 - b) / a, (h1 - b) / a
    assert 0 <= T0 < T1 <= T_REC + 1, f"window [{T0:.1f},{T1:.1f}] outside recording 0–{T_REC:.1f}"
    T1 = min(T1, T_REC)

    ids = pop.index.values
    sel = np.isin(clu, ids)
    order = np.argsort(clu[sel], kind="stable")
    t_s, c_s, a_s = st_all[sel][order], clu[sel][order], amps[sel][order]
    uids, starts, counts = np.unique(c_s, return_index=True, return_counts=True)
    ST, ST_PAD, AMP = {}, {}, {}
    for c, s0, n in zip(uids, starts, counts):
        t, am = t_s[s0:s0 + n], a_s[s0:s0 + n]
        w = (t >= T0) & (t < T1)
        ST[int(c)], AMP[int(c)] = t[w], am[w]
        ST_PAD[int(c)] = t[(t >= T0 - EDGE_PAD) & (t < T1 + EDGE_PAD)]
    del st_all, clu, amps

    # templates (un-whitened, spike-count-weighted over window spikes) → waveform shape metrics
    templates_w = np.load(files["templates.npy"], mmap_mode="r")
    winv = np.load(files["whitening_mat_inv.npy"])
    chan_pos = np.load(files["channel_positions.npy"])
    spk_tmpl = np.asarray(np.load(files["spike_templates.npy"], mmap_mode="r")).ravel()[sel][order]

    rows = []
    for c in pop.index:
        st = ST.get(int(c), np.empty(0))
        m = dict(phyID=int(c))
        m.update(train_metrics(st, AMP.get(int(c), np.empty(0)), T0, T1))
        idx = np.nonzero(c_s == c)[0]
        tw = spk_tmpl[idx][(t_s[idx] >= T0) & (t_s[idx] < T1)]
        if len(tw):
            tid, tn = np.unique(tw, return_counts=True)
            tmpl = np.tensordot(tn / tn.sum(), np.asarray(templates_w[tid]), axes=1) @ winv
            ptp = tmpl.max(0) - tmpl.min(0)
            pk = int(np.argmax(ptp))
            wf = tmpl[:, pk]
            tr = int(np.argmin(wf))
            pk_after = tr + int(np.argmax(wf[tr:]))
            fp = ptp > FOOTPRINT_FRAC * ptp[pk]
            m.update(peak_ch=pk, depth_um=float(chan_pos[pk, 1]), x_um=float(chan_pos[pk, 0]),
                     trough_peak_ms=(pk_after - tr) / SR * 1e3,
                     n_ch_footprint=int(fp.sum()),
                     footprint_span_um=float(np.ptp(chan_pos[fp, 1])) if fp.sum() > 1 else 0.0,
                     n_templates_used=len(tid), wf_polarity=float(wf.max() / -wf.min()) if wf.min() < 0 else np.nan)
        rows.append(m)
    units = pop.join(pd.DataFrame(rows).set_index("phyID"))
    units.index.name = "phyID"

    # CS→SS pairing (SS pool = extra-good CB non-CS units, tier 1, as in step0)
    cs_ids = units.index[units["c4_label"] == CS_LABEL].tolist()
    ss_pool = units.index[units["c4_label"] != CS_LABEL].tolist()
    pr = []
    for cs in cs_ids:
        if len(ST[cs]) < 20:
            continue
        cand = []
        for ss in ss_pool:
            if len(ST_PAD.get(ss, [])) < 100:
                continue
            m = pause_metrics(ST[cs], ST_PAD[ss])
            m.update(cs=cs, ss=ss, ss_c4=units.at[ss, "c4_label"],
                     dist_um=float(np.hypot(units.at[cs, "x_um"] - units.at[ss, "x_um"],
                                            units.at[cs, "depth_um"] - units.at[ss, "depth_um"])))
            cand.append(m)
        cand = sorted(cand, key=lambda d: (np.nan_to_num(d["pause_ratio"], nan=9), -d["asym"]))
        for k, m in enumerate(cand):
            if k < N_JITTER_TOP:
                m.update(jitter_pause_test(ST[cs], ST_PAD[m["ss"]]))
            m["rank_pr"] = k + 1
        pr.extend(cand)
    pairs = pd.DataFrame(pr)
    if len(pairs):
        pairs["candidate"] = ((pairs["pause_ratio"] < CAND_PAUSE_RATIO_MAX) & (pairs["asym"] > CAND_ASYM_MIN)
                              & ~pairs["collision_dip"] & (pairs.get("jitter_p", 1.0).fillna(1.0) < CAND_JITTER_P_MAX))

    # CS–CS duplicate evidence
    dr = []
    for i, ca in enumerate(cs_ids):
        for cb in cs_ids[i + 1:]:
            if len(ST[ca]) < 20 or len(ST[cb]) < 20:
                continue
            m = cs_cs_coinc(ST[ca], ST[cb])
            m.update(cs_a=ca, cs_b=cb, dist_um=float(np.hypot(units.at[ca, "x_um"] - units.at[cb, "x_um"],
                                                              units.at[ca, "depth_um"] - units.at[cb, "depth_um"])))
            dr.append(m)
    dups = pd.DataFrame(dr)

    # per-unit pairing summary
    units["has_ss_partner"] = False
    units["n_ss_candidates"] = 0
    units["partner_c4"] = ""
    units["has_cs_partner"] = False
    if len(pairs):
        cnt = pairs.groupby("cs")["candidate"].sum()
        units.loc[cnt.index, "n_ss_candidates"] = cnt.values
        units.loc[cnt.index, "has_ss_partner"] = cnt.values > 0
        best = pairs[pairs["candidate"]].sort_values("pause_ratio").groupby("cs").first()
        units.loc[best.index, "partner_c4"] = best["ss_c4"].values
        units.loc[units.index.isin(pairs.loc[pairs["candidate"], "ss"]), "has_cs_partner"] = True
        bpr = pairs[pairs["rank_pr"] == 1].set_index("cs")
        units.loc[bpr.index, "best_pause_ratio"] = bpr["pause_ratio"].values
    units["max_coinc_cs"] = np.nan
    if len(dups):
        for c in cs_ids:
            v = dups.loc[(dups["cs_a"] == c) | (dups["cs_b"] == c), "coinc_frac_peak"]
            if len(v):
                units.at[c, "max_coinc_cs"] = v.max()

    for k, v in dict(mouse=row.mouse, S=row.S, R=row.R, sorter=row.sorter, sorter_folder=row.sorter_folder,
                     imec=row.imec, win_start_imec=T0, win_end_imec=T1, win_dur_s=T1 - T0, rec_dur_s=T_REC,
                     trials_first=int(np.min(trials)), trials_last=int(np.max(trials)), sample_rate=SR).items():
        units[k] = v
    return units.reset_index(), pairs, dups


all_units = []
for row in sessions.itertuples():
    d = safe_out_dir(SES_DIR / f"{row.mouse}_{row.S}")
    f_u, f_p, f_d = d / "units.csv", d / "cs_ss_pairs.csv", d / "cs_cs_pairs.csv"
    if f_u.exists() and not FORCE:
        all_units.append(pd.read_csv(f_u))
        print(f"cached  {row.mouse} {row.S}")
        continue
    t0 = time.time()
    try:
        u, p, dd = process_session(row)
    except Exception as e:  # report and continue; failed sessions are listed below
        print(f"FAILED  {row.mouse} {row.S} ({row.sorter}): {type(e).__name__}: {e}")
        continue
    u.to_csv(f_u, index=False)
    p.to_csv(f_p, index=False)
    dd.to_csv(f_d, index=False)
    all_units.append(u)
    ncs = (u["c4_label"] == CS_LABEL).sum()
    print(f"done    {row.mouse} {row.S} ({row.sorter}) in {time.time() - t0:.0f} s | win {u['win_dur_s'].iloc[0]:.0f} s "
          f"of {u['rec_dur_s'].iloc[0]:.0f} s | eg CB {len(u)} | CS {ncs} | SS {(u['c4_label'] == SS_LABEL).sum()}")

units = pd.concat(all_units, ignore_index=True)
units.to_csv(OUT_DIR / "units_all.csv", index=False)
missing = sessions.merge(units[["mouse", "S"]].drop_duplicates(), how="left", indicator=True)
print("sessions not processed:", missing.loc[missing["_merge"] == "left_only", ["mouse", "S"]].values.tolist())


# %% 5. SESSION SUMMARY (raw numbers)
def sess_summary(g):
    cs, ss = g[g["c4_label"] == CS_LABEL], g[g["c4_label"] == SS_LABEL]
    r = dict(sorter=g["sorter"].iloc[0], win_dur_min=g["win_dur_s"].iloc[0] / 60,
             n_eg_CB=len(g), n_CS=len(cs), n_SS=len(ss),
             n_MLI=(g["c4_label"] == "MLI").sum(), n_GoC=(g["c4_label"] == "GoC").sum(),
             n_MFB=(g["c4_label"] == "MFB").sum(), n_unlabelled=(g["c4_label"] == "unlabelled").sum())
    r["frac_CS_of_eg"] = r["n_CS"] / r["n_eg_CB"] if r["n_eg_CB"] else np.nan
    r["frac_SS_of_eg"] = r["n_SS"] / r["n_eg_CB"] if r["n_eg_CB"] else np.nan
    r["CS_SS_ratio"] = r["n_CS"] / r["n_SS"] if r["n_SS"] else np.nan
    r["frac_CS_paired"] = cs["has_ss_partner"].mean() if len(cs) else np.nan
    r["frac_SS_paired"] = ss["has_cs_partner"].mean() if len(ss) else np.nan
    r["n_CS_dup_cand"] = int((cs["max_coinc_cs"] > DUP_COINC_MIN).sum())
    r["frac_CS_dup_strong"] = (cs["max_coinc_cs"] > DUP_COINC_STRONG).mean() if len(cs) else np.nan
    for lab, sub in (("CS", cs), ("SS", ss)):
        for k in ("fr", "cv_isi", "cv2", "frac_isi_lt1ms", "frac_isi_lt5ms", "presence_ratio", "fr_cv_60s",
                  "amp_rcv", "trough_peak_ms", "footprint_span_um", "c4_conf", "c4_vote_ratio"):
            r[f"{lab}_{k}_med"] = sub[k].median() if len(sub) else np.nan
    return pd.Series(r)


summ = units.groupby(["mouse", "S"]).apply(sess_summary).reset_index()
summ.to_csv(OUT_DIR / "session_summary.csv", index=False)
with pd.option_context("display.width", 250, "display.max_columns", 40):
    print("--- raw numbers per session (extra-good CB units, C4 labels)")
    print(summ[["mouse", "S", "sorter", "win_dur_min", "n_eg_CB", "n_CS", "n_SS", "n_MLI", "n_GoC", "n_MFB",
                "n_unlabelled", "frac_CS_paired", "frac_SS_paired", "n_CS_dup_cand"]].round(2).to_string(index=False))
    print("--- totals per sorter")
    tot = summ.groupby("sorter")[["n_eg_CB", "n_CS", "n_SS", "n_MLI", "n_GoC", "n_MFB", "n_unlabelled",
                                  "n_CS_dup_cand", "win_dur_min"]].sum()
    tot.insert(0, "n_sessions", summ.groupby("sorter").size())
    print(tot.round(1).to_string())
    print("--- per mouse × sorter: sessions, CS, SS")
    print(summ.pivot_table(index="mouse", columns="sorter", values=["n_CS", "n_SS"], aggfunc=["count", "sum"]).to_string())


# %% 6. SESSION-LEVEL STATISTICS (Mann–Whitney + within-mouse permutation)
def within_mouse_perm(df, col, n_perm=N_PERM, seed=0):
    d = df[["mouse", "sorter", col]].dropna()
    groups = [g for _, g in d.groupby("mouse") if g["sorter"].nunique() == 2]
    if not groups:
        return np.nan, np.nan, 0

    def stat(labs):
        diffs = []
        for g, l in zip(groups, labs):
            v = g[col].values
            diffs.append(np.median(v[l == "ks4"]) - np.median(v[l == "ibl"]))
        return np.mean(diffs)

    labs = [g["sorter"].values for g in groups]
    obs = stat(labs)
    rng = np.random.default_rng(seed)
    null = np.array([stat([rng.permutation(l) for l in labs]) for _ in range(n_perm)])
    p = (1 + np.sum(np.abs(null) >= abs(obs) - 1e-12)) / (1 + n_perm)
    return obs, p, len(groups)


METRICS = (["n_eg_CB", "n_CS", "n_SS", "frac_CS_of_eg", "frac_SS_of_eg", "CS_SS_ratio",
            "frac_CS_paired", "frac_SS_paired", "frac_CS_dup_strong", "win_dur_min"]
           + [c for c in summ.columns if c.endswith("_med")])
st_rows = []
for col in METRICS:
    k, i = summ.loc[summ["sorter"] == "ks4", col].dropna(), summ.loc[summ["sorter"] == "ibl", col].dropna()
    if len(k) < 2 or len(i) < 2:
        continue
    p_mw = mannwhitneyu(k, i, alternative="two-sided").pvalue
    d_wm, p_wm, n_m = within_mouse_perm(summ, col)
    st_rows.append(dict(metric=col, ks4_median=k.median(), ibl_median=i.median(), n_ks4=len(k), n_ibl=len(i),
                        p_mannwhitney=p_mw, within_mouse_diff_ks4_minus_ibl=d_wm, p_within_mouse_perm=p_wm,
                        n_mice=n_m))
stats = pd.DataFrame(st_rows)
stats.to_csv(OUT_DIR / "stats_session_level.csv", index=False)
with pd.option_context("display.width", 250, "display.max_rows", 100):
    print("--- session-level comparison (KS4 vs IBL); uncorrected p, exploratory")
    print(stats.round(4).to_string(index=False))

# unit-level descriptive (pooled; NOT independent → descriptive only)
desc = (units[units["c4_label"].isin([CS_LABEL, SS_LABEL])]
        .groupby(["c4_label", "sorter"])[["fr", "cv_isi", "cv2", "frac_isi_lt1ms", "frac_isi_lt5ms", "presence_ratio",
                                          "fr_cv_60s", "fr_last_first", "amp_rcv", "amp_last_first", "trough_peak_ms",
                                          "n_ch_footprint", "footprint_span_um", "n_templates_used", "c4_conf",
                                          "c4_vote_ratio", "c4_prob"]]
        .median().T)
with pd.option_context("display.width", 250):
    print("--- unit-level medians (pooled units; descriptive)")
    print(desc.round(3).to_string())
    print("--- CS units: pairing / duplicate rates pooled")
    cs_u = units[units["c4_label"] == CS_LABEL]
    print(cs_u.groupby("sorter").agg(n=("phyID", "size"), paired=("has_ss_partner", "mean"),
                                     dup_cand=("max_coinc_cs", lambda v: (v > DUP_COINC_MIN).mean()),
                                     dup_strong=("max_coinc_cs", lambda v: (v > DUP_COINC_STRONG).mean())).round(3))
    print("    partner C4 label of paired CS:")
    print(pd.crosstab(cs_u.loc[cs_u["has_ss_partner"], "sorter"], cs_u.loc[cs_u["has_ss_partner"], "partner_c4"]))


# %% 7. FIGURES
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.2,
                     "font.size": 9, "savefig.dpi": 150, "figure.facecolor": "white"})
rng = np.random.default_rng(0)


def strip(ax, df, col, by_session=True, log=False):
    for xi, s in enumerate(["ks4", "ibl"]):
        sub = df[df["sorter"] == s].dropna(subset=[col])
        for mouse, g in sub.groupby("mouse"):
            ax.scatter(xi + rng.uniform(-0.18, 0.18, len(g)), g[col], s=26 if by_session else 9,
                       marker=MOUSE_MARKER[mouse], facecolor=SORTER_COLOR[s], edgecolor="white", linewidth=0.6,
                       alpha=0.9 if by_session else 0.55, zorder=3)
        if len(sub):
            med = sub[col].median()
            ax.plot([xi - 0.3, xi + 0.3], [med, med], color="#333333", lw=2, zorder=4)
    ax.set_xticks([0, 1], [SORTER_NAME["ks4"], SORTER_NAME["ibl"]])
    ax.set_xlim(-0.6, 1.6)
    if log:
        ax.set_yscale("log")


def mouse_legend(fig):
    h = [plt.Line2D([], [], marker=MOUSE_MARKER[m], ls="", color="#666666", label=m) for m in MICE]
    fig.legend(handles=h, loc="lower center", ncol=5, frameon=False)


def ptxt(col):
    r = stats.set_index("metric")
    if col not in r.index:
        return ""
    return f"MW p={r.at[col, 'p_mannwhitney']:.3f} | within-mouse p={r.at[col, 'p_within_mouse_perm']:.3f}"


# Fig 1: raw numbers per session
cols1 = [("n_eg_CB", "extra-good CB units"), ("n_CS", "CS units (PkC_cs)"), ("n_SS", "SS units (PkC_ss)"),
         ("frac_CS_of_eg", "CS / extra-good CB"), ("CS_SS_ratio", "CS / SS units"), ("win_dur_min", "window (min)"),
         ("frac_CS_paired", "CS with SS-pause partner"), ("frac_SS_paired", "SS with CS partner"),
         ("frac_CS_dup_strong", "CS with likely duplicate")]
fig, axs = plt.subplots(3, 3, figsize=(10, 9))
for ax, (c, t) in zip(axs.ravel(), cols1):
    strip(ax, summ, c)
    ax.set_title(f"{t}\n{ptxt(c)}", fontsize=8)
fig.suptitle("Per-session counts and rates (one point per session; bar = median)")
mouse_legend(fig)
fig.tight_layout(rect=(0, 0.04, 1, 0.97))
fig.savefig(FIG_DIR / "fig1_session_counts.png")
plt.close(fig)

# Fig 2/3: unit-level CS and SS statistics
ucols = [("fr", "firing rate (Hz)", True), ("cv_isi", "CV ISI", False), ("cv2", "CV2", False),
         ("frac_isi_lt1ms", "ISI < 1 ms (frac)", False), ("frac_isi_lt5ms", "ISI < 5 ms (frac)", False),
         ("presence_ratio", "presence ratio (60 s bins)", False), ("fr_cv_60s", "rate CV across 60 s bins", False),
         ("amp_rcv", "amplitude IQR/median", False), ("amp_last_first", "amplitude last/first 20%", False),
         ("trough_peak_ms", "trough→peak (ms)", False), ("footprint_span_um", "footprint span (µm)", False),
         ("c4_conf", "C4 confidence ratio", False)]
for lab, fname in ((CS_LABEL, "fig2_CS_unit_stats.png"), (SS_LABEL, "fig3_SS_unit_stats.png")):
    sub = units[units["c4_label"] == lab]
    fig, axs = plt.subplots(3, 4, figsize=(13, 9))
    for ax, (c, t, lg) in zip(axs.ravel(), ucols):
        strip(ax, sub, c, by_session=False, log=lg)
        sc = f"{'CS' if lab == CS_LABEL else 'SS'}_{c}_med"
        ax.set_title(f"{t}\nsession-median: {ptxt(sc)}", fontsize=7.5)
    n = sub.groupby("sorter").size()
    fig.suptitle(f"{lab} extra-good CB units in win_selected (one point per unit; KS4 n={n.get('ks4', 0)}, "
                 f"IBL n={n.get('ibl', 0)}). p-values are session-level.")
    mouse_legend(fig)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    fig.savefig(FIG_DIR / fname)
    plt.close(fig)

# Fig 4: CS rate & doublets per mouse (within-mouse view)
fig, axs = plt.subplots(1, 3, figsize=(12, 3.8))
for ax, (c, t) in zip(axs, [("n_CS", "CS units"), ("CS_fr_med", "median CS rate (Hz)"),
                            ("frac_CS_paired", "CS with SS partner")]):
    for mi, m in enumerate(MICE):
        g = summ[summ["mouse"] == m]
        for s, dx in (("ks4", -0.15), ("ibl", 0.15)):
            v = g.loc[g["sorter"] == s, c]
            ax.scatter(np.full(len(v), mi + dx), v, color=SORTER_COLOR[s], s=24, edgecolor="white", linewidth=0.6)
    ax.set_xticks(range(len(MICE)), [m.split("_", 1)[1] for m in MICE], rotation=30)
    ax.set_title(f"{t}\n{ptxt(c)}", fontsize=8)
h = [plt.Line2D([], [], marker="o", ls="", color=SORTER_COLOR[s], label=SORTER_NAME[s]) for s in ("ks4", "ibl")]
fig.legend(handles=h, loc="upper right", frameon=False)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig4_within_mouse.png")
plt.close(fig)
print(f"figures in {FIG_DIR}")
