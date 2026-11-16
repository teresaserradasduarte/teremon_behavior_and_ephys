# -*- coding: utf-8 -*-
"""
raw_waveforms_figures.py — static figures from raw_waveforms_session.py outputs (one session; first pass CoteDor R5).

Reads only local files in CS_curation/<mouse>_<R#>/raw_waveforms/ and writes PNGs to raw_waveforms/figures/:
  session_summary.png                 all CS: spikelets, capture, pairing re-check (pause in 0–8 vs 3–8 ms)
  cs_<id>_events.png                  every event of one CS as a heat map (sorted by spikelet count) + which
                                      capturing clusters fired in each event (shown in the browser's Raw tab)
  cs_<id>_summary.png                 multi-panel page per CS (mean, spikelets, time course, conditional means,
                                      re-check CCGs)
  dup_<a>_<b>.png                     CS–CS pairs: raw means of coincident / a-only / b-only events
"""
from pathlib import Path
import os as _os
import json
import pickle

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

LOCAL_ROOT = Path(r"C:\Users\simon\Desktop\Data\ephys_and_behavior")
MOUSE, SESSION_R = "1_CoteDor", "R5"
MOUSE, SESSION_R = (_os.environ["RAW_SESSION"].rsplit("_", 1) if _os.environ.get("RAW_SESSION") else (MOUSE, SESSION_R))  # batch: raw_all_sessions.py
RAW_DIR = LOCAL_ROOT / "CS_curation" / f"{MOUSE}_{SESSION_R}" / "raw_waveforms"
STEP0_DIR = LOCAL_ROOT / "CS_curation" / f"{MOUSE}_{SESSION_R}" / "step0_explore_bhvspan_v7batch"   # 2026-11-09: v7 (was v5batch)
FIG = RAW_DIR / "figures"
FIG.mkdir(exist_ok=True)
HUE = ["#3380cc", "#d1495b", "#2a9d55", "#e08a00", "#7a4fb5", "#00838f", "#8d6e63", "#c2185b"]
CB = "#9a7bc1"
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110})

ev = np.load(RAW_DIR / "cs_events.npz")
mz = np.load(RAW_DIR / "raw_means.npz")
X = pickle.load(open(RAW_DIR / "raw_extra.pkl", "rb"))
cap, recheck, dups, cond, locked = X["cap"], X["recheck"], X["dups"], X["cond"], X["locked"]
rc_ccgs, rc_edges = X["rc_ccgs"], X["rc_edges"]
spk = pd.read_csv(RAW_DIR / "cs_raw_spikelets.csv").set_index("cs")
units = pd.read_csv(STEP0_DIR / "tables" / "unit_table.csv").set_index("phyID")
T = ev["t_ms"]; chan_pos = mz["chan_pos"]
MEANS = dict(zip(mz["units"].tolist(), mz["means"]))
CS_IDS = ev["cs_ids"].tolist()
PEAK = dict(zip(CS_IDS, ev["peaks"].tolist())); SD = dict(zip(CS_IDS, ev["sd"].tolist()))
POL = dict(zip(CS_IDS, ev["pol"].tolist()))
PRM = json.loads((RAW_DIR / "params_used.json").read_text())
THR = PRM.get("SPK_THR", 3.5)


def evget(name, c):
    return ev[f"{name}_{c}"]


def members(c):
    return {int(k.split("_")[2]): ev[k] for k in ev.files if k.startswith(f"mem_{c}_")}


def caps_of(c):
    return cap[cap.cs == c].sort_values("excess", ascending=False)


# ---------------------------------------------------------------- session summary
def session_summary():
    fig = plt.figure(figsize=(18, 17))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[2.2, 1.3, 5], hspace=.45, wspace=.3, top=.95)
    order = sorted(CS_IDS, key=lambda c: units.at[c, "y"])
    lab = [f"{c} ({units.at[c, 'y']:.0f})" for c in order]
    # spikelet detections per CS (stacked fractions)
    ax = fig.add_subplot(gs[0, :])
    cols = ["frac_0", "frac_1", "frac_2", "frac_3", "frac_4plus"]
    bottom = np.zeros(len(order))
    for i, k in enumerate(cols):
        v = spk.loc[order, k].values
        ax.bar(range(len(order)), v, bottom=bottom, color=plt.cm.Purples(.2 + .18 * i), label=k.replace("frac_", "").replace("plus", "+") + " detections")
        bottom += v
    ax.set_xticks(range(len(order)), lab, fontsize=7, rotation=90); ax.set_xlim(-.6, len(order) - .4)
    ax.set_ylabel("fraction of events"); ax.legend(ncol=5, fontsize=7, loc="upper center", bbox_to_anchor=(.5, -.42), frameon=False)
    ax.set_title("Raw spikelet detections per CS event on the CS peak channel (x: CS id (depth µm), sorted by depth)")
    # excess over random times
    ax = fig.add_subplot(gs[1, :2])
    ex = spk.loc[order, "excess_mean"].values
    ax.bar(range(len(order)), ex, color=[CB if v > .2 else "0.7" for v in ex])
    ax.scatter(range(len(order)), spk.loc[order, "null_mean"].values, marker="_", s=60, c="k", label="random-time detections/event", zorder=4)
    ax.axhline(0, c="k", lw=.6); ax.set_xticks(range(len(order)), [str(c) for c in order], fontsize=6.5, rotation=90)
    ax.set_xlim(-.6, len(order) - .4); ax.set_ylabel("spikelets / event"); ax.legend(fontsize=7, frameon=False)
    ax.set_title("Excess spikelets per event (CS events − random times); purple > 0.2. The rest is background on the channel")
    # capture count per CS
    ax = fig.add_subplot(gs[1, 2])
    nc = cap.groupby("cs").size().reindex(order).fillna(0)
    top = cap.groupby("cs")["excess"].max().reindex(order).fillna(0)
    ax.scatter(nc, top, c=CB, s=22)
    for c in order:
        if nc[c] >= 16 or top[c] > .45:
            ax.annotate(str(c), (nc[c], top[c]), fontsize=7, xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("n clusters capturing part of the event"); ax.set_ylabel("largest capture excess")
    ax.set_title("How spread is each CS event over clusters?")
    # pairing re-check: pr 0–8 vs pr 3–8 for all pairs, highlight step0 #1 and best post
    ax = fig.add_subplot(gs[2, 0])
    r = recheck
    ax.scatter(r.pr_0_8, r.pr_3_8, s=4, c="0.75", label="all CS×SS pairs")
    s0 = r[r.step0_rank == 1]; bp = r[r.post_rank == 1]
    ax.scatter(s0.pr_0_8, s0.pr_3_8, s=18, c=CB, label="step0 rank 1")
    ax.scatter(bp.pr_0_8, bp.pr_3_8, s=18, facecolors="none", edgecolors="k", label="best pause 3–8 ms")
    for c in bp.cs:
        b = bp[bp.cs == c].iloc[0]
        if b.pr_3_8 < .75:
            ax.annotate(f"{c}→{b.ss}", (b.pr_0_8, b.pr_3_8), fontsize=6.5, xytext=(3, -6), textcoords="offset points")
    ax.plot([0, 2], [0, 2], c="0.6", lw=.7, ls=":"); ax.set_xlim(0, 1.6); ax.set_ylim(0, 1.6)
    ax.set_xlabel("pause ratio 0–8 ms (v5 window)"); ax.set_ylabel("pause ratio 3–8 ms (after the event)")
    ax.legend(fontsize=7); ax.set_title("Pairing re-check: below the diagonal = pause hidden by the CS event")
    # table of CS whose best post-window SS differs from step0 #1 (and pauses clearly)
    ax = fig.add_subplot(gs[2, 1:]); ax.axis("off")
    rows = []
    for c in order:
        a, b = s0[s0.cs == c], bp[bp.cs == c]
        if not len(a) or not len(b):
            continue
        a, b = a.iloc[0], b.iloc[0]
        flag = "DIFFERS" if a.ss != b.ss and b.pr_3_8 < .75 else ("" if a.ss == b.ss else "differs (weak)")
        rows.append([c, int(a.ss), f"{a.pr_0_8:.2f}", f"{a.pr_3_8:.2f}", int(b.ss), f"{b.pr_0_8:.2f}", f"{b.pr_3_8:.2f}",
                     f"{b.dist_um:.0f}", f"{b.ss_fr:.0f}", flag])
    tb = ax.table(cellText=rows, colLabels=["CS", "step0 #1", "pr0-8", "pr3-8", "best 3-8", "pr0-8", "pr3-8", "d µm", "FR", ""],
                  loc="upper center", cellLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(6.5); tb.scale(1, 1.05)
    for i, rw in enumerate(rows, start=1):
        if rw[-1] == "DIFFERS":
            for j in range(len(rw)):
                tb[i, j].set_facecolor("#f6d6d6")
    ax.set_title("step0 rank-1 SS vs SS with the deepest pause after the event (3–8 ms); red = clear disagreement")
    fig.suptitle(f"{MOUSE} {SESSION_R} — raw-data layer, session summary", fontsize=11, y=.97)
    fig.savefig(FIG / "session_summary.png", bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- per CS: event heat map + membership raster
def events_fig(c):
    tr = evget("tr", c).astype(np.float32); cnt = evget("cnt", c); amp = evget("amp", c)
    mem = members(c); cp = caps_of(c)
    ks = [int(k) for k in cp.unit if int(k) in mem][:10]
    o = np.lexsort((amp, cnt))                                           # by spikelet count, then main amplitude
    fig = plt.figure(figsize=(9, 6.2))
    gs = GridSpec(1, 3, width_ratios=[6, .25 * max(len(ks), 1), .6], figure=fig, wspace=.05)
    ax = fig.add_subplot(gs[0])
    v = np.percentile(np.abs(tr), 99)
    ax.imshow(tr[o], aspect="auto", cmap="RdBu_r", vmin=-v, vmax=v, extent=[T[0], T[-1], len(o), 0], interpolation="nearest")
    ax.axvline(0, c="k", lw=.5, alpha=.4)
    # boundaries between spikelet-count groups
    for b in np.flatnonzero(np.diff(cnt[o])) + 1:
        ax.axhline(b, c="k", lw=.6)
    for n in np.unique(cnt):
        ii = np.flatnonzero(cnt[o] == n)
        ax.text(T[-1] - .1, ii.mean(), f"{n} spk", ha="right", va="center", fontsize=7, color="k",
                bbox=dict(fc="w", ec="none", alpha=.7, pad=1))
    ax.set_xlabel("time after CS spike (ms)"); ax.set_ylabel("event (sorted by spikelet count, then amplitude)")
    ax.set_title(f"CS {c}: all {len(o)} events, peak ch {PEAK[c]} (blue = negative, ±{v:.0f} µV)")
    ax2 = fig.add_subplot(gs[1], sharey=ax)
    if ks:
        M = np.stack([mem[k][o] for k in ks], 1).astype(float)
        rgb = np.ones((len(o), len(ks), 3))
        for j, k in enumerate(ks):
            col = np.array(matplotlib.colors.to_rgb(HUE[j % 8]))
            rgb[M[:, j] > 0, j] = col
        ax2.imshow(rgb, aspect="auto", extent=[0, len(ks), len(o), 0], interpolation="nearest")
        ax2.set_xticks(np.arange(len(ks)) + .5, [str(k) for k in ks], rotation=90, fontsize=7)
    ax2.tick_params(labelleft=False); ax2.set_title("fired", fontsize=7)
    ax3 = fig.add_subplot(gs[2], sharey=ax)
    ax3.plot(amp[o], np.arange(len(o)) + .5, c="0.3", lw=.5); ax3.tick_params(labelleft=False)
    ax3.set_xlabel("main µV", fontsize=7); ax3.invert_xaxis()
    fig.savefig(FIG / f"cs_{c}_events.png", bbox_inches="tight", dpi=100); plt.close(fig)


# ---------------------------------------------------------------- per CS: summary page
def summary_fig(c):
    tr = evget("tr", c).astype(np.float32); cnt = evget("cnt", c); amp = evget("amp", c); et = evget("t", c)
    lags = evget("lag", c); cp = caps_of(c); mem = members(c)
    fig = plt.figure(figsize=(16, 11)); gs = GridSpec(3, 4, figure=fig, hspace=.45, wspace=.3)
    # mean ± SD on the peak channel, with spikelet lag histogram
    ax = fig.add_subplot(gs[0, 0:2])
    m, s = tr.mean(0), tr.std(0)
    ax.fill_between(T, m - s, m + s, color=CB, alpha=.2, lw=0); ax.plot(T, m, c=CB, lw=1.6, label="mean ± SD")
    for n, col in [(0, "0.6"), (1, "#3380cc"), (2, "#2a9d55"), (3, "#d1495b")]:
        if (cnt == n).sum() >= 10:
            ax.plot(T, tr[cnt == n].mean(0), c=col, lw=.9, label=f"{n} spikelets (n {(cnt == n).sum()})")
    ax.axhline(POL[c] * THR * SD[c], c="k", ls=":", lw=.7, label=f"spikelet threshold ({THR}× residual SD)")
    ax.set_xlabel("ms"); ax.set_ylabel("µV"); ax.legend(fontsize=6.5, ncol=2); ax.set_title(f"CS {c}: peak ch {PEAK[c]}, means by spikelet count")
    axh = ax.twinx(); be = np.arange(0.5, 8.01, .1); nl = ev[f"nulllag_{c}"]
    axh.hist(lags, be, color="k", alpha=.18)
    axh.hist(nl, be, weights=np.full(len(nl), len(cnt) / max(len(ev[f"null_{c}"]), 1)), histtype="step", color="k", lw=.6, alpha=.6); axh.set_yticks([]); axh.spines["right"].set_visible(False)
    for j, r in enumerate(cp.head(8).itertuples()):
        ax.axvline(r.lag_med_ms, c=HUE[j % 8], ls="--", lw=.8)
    # spikelet count histogram
    ax = fig.add_subplot(gs[0, 2])
    h = np.bincount(np.minimum(cnt, 6), minlength=7) / len(cnt); nc = ev[f"null_{c}"]
    hn = np.bincount(np.minimum(nc, 6), minlength=7) / max(len(nc), 1)
    ax.bar(np.arange(7) - .2, h, .4, color=CB, label=f"CS events (mean {cnt.mean():.2f})")
    ax.bar(np.arange(7) + .2, hn, .4, color="0.7", label=f"random times (mean {nc.mean():.2f})")
    ax.set_xticks(range(7), ["0", "1", "2", "3", "4", "5", "6+"]); ax.legend(fontsize=6.5)
    ax.set_xlabel("detections per event"); ax.set_ylabel("fraction")
    ax.set_title(f"spikelets: excess {cnt.mean() - nc.mean():+.2f}/event (noise SD {SD[c]:.0f} µV, pol {POL[c]:+.0f})")
    # count vs amplitude
    ax = fig.add_subplot(gs[0, 3])
    ax.boxplot([amp[np.minimum(cnt, 4) == i] if (np.minimum(cnt, 4) == i).any() else [np.nan] for i in range(5)], positions=range(5), widths=.6, showfliers=False)
    ax.set_xticks(range(5), ["0", "1", "2", "3", "4+"]); ax.set_xlabel("spikelets"); ax.set_ylabel("main spike µV")
    ax.set_title("main-spike amplitude vs spikelet count")
    # over time
    ax = fig.add_subplot(gs[1, 0:2])
    tb = np.floor(et / 60).astype(int); ub = np.unique(tb)
    ax.plot(ub, [cnt[tb == b].mean() for b in ub], c=CB, lw=1.4); ax.set_ylabel("spikelets / event", color=CB)
    ax2 = ax.twinx(); ax2.plot(ub, [np.median(amp[tb == b]) for b in ub], c="0.4", lw=1); ax2.set_ylabel("main µV (median)", color="0.4")
    ax2.spines["right"].set_visible(True); ax.set_xlabel("time (min)"); ax.set_title("over the recording (60 s bins)")
    # capture: which clusters, how often, at which lag
    ax = fig.add_subplot(gs[1, 2:])
    if len(cp):
        y = np.arange(len(cp))
        ax.barh(y, cp.frac_cs_events, color="0.85", label="CS events with a spike")
        ax.barh(y, cp.chance, color="0.6", label="chance")
        for j, r in enumerate(cp.itertuples()):
            ax.text(1.01, j, f"lag {r.lag_med_ms:.2f} ms · {r.dist_um:.0f} µm · {r.fr:.0f} Hz · {r.role if isinstance(r.role, str) else ''} {r.group}",
                    va="center", fontsize=6.5, transform=ax.get_yaxis_transform())
        ax.set_yticks(y, [str(u) for u in cp.unit]); ax.invert_yaxis(); ax.set_xlim(0, 1)
        ax.legend(fontsize=7, loc="lower right"); ax.set_xlabel("fraction of CS events")
    ax.set_title("clusters capturing part of the event (−0.5…+4 ms)")
    # conditional means (top 3 capturing clusters) on the CS peak channel
    cd = {d["unit"]: d for d in cond.get(c, [])}
    pk = PEAK[c]
    for j, r in enumerate([r for r in cp.itertuples() if int(r.unit) in cd][:3]):
        ax = fig.add_subplot(gs[2, j]); d = cd[int(r.unit)]
        i = d["ch"].index(pk) if pk in d["ch"] else int(np.argmax(np.ptp(d["without"], 0)))
        ax.plot(T, d["without"][:, i], c="0.5", lw=1.3, label=f"without (n {d['n_without']})")
        ax.plot(T, d["with_"][:, i], c=HUE[j], lw=1.5, label=f"with {r.unit} (n {d['n_with']})")
        ax.axvline(r.lag_med_ms, c=HUE[j], ls="--", lw=.8)
        ax.set_title(f"CS {c} events with / without cluster {int(r.unit)}"); ax.legend(fontsize=6.5); ax.set_xlabel("ms")
    # re-check CCGs: step0 #1 vs best post
    ax = fig.add_subplot(gs[2, 3])
    rc = recheck[recheck.cs == c]
    ctr = (rc_edges[:-1] + rc_edges[1:]) / 2 * 1e3
    picks = list(dict.fromkeys(rc.nsmallest(2, "step0_rank").ss.tolist() + rc.nsmallest(2, "pr_3_8").ss.tolist()))
    for j, k in enumerate(picks[:4]):
        y = rc_ccgs.get((c, int(k)))
        if y is None:
            continue
        rr = rc[rc.ss == k].iloc[0]
        ax.plot(ctr, np.convolve(y, np.ones(3) / 3, "same"), lw=1.1, c=HUE[j],
                label=f"SS {k}: #{rr.step0_rank} step0, pr {rr.pr_0_8:.2f} / {rr.pr_3_8:.2f}")
    ax.axvspan(0, 8, color="#e4e4e4", zorder=0); ax.axvspan(3, 8, color="#cfe0f1", alpha=.5, zorder=0)
    ax.axhline(1, c="0.6", lw=.6); ax.set_xlim(-20, 50); ax.set_ylim(0, None)
    ax.set_xlabel("ms after CS"); ax.set_ylabel("SS rate / baseline"); ax.legend(fontsize=6)
    ax.set_title("pairing re-check (grey 0–8, blue 3–8 ms)")
    fig.suptitle(f"{MOUSE} {SESSION_R} CS {c} (y {units.at[c, 'y']:.0f} µm) — raw data", fontsize=11)
    fig.savefig(FIG / f"cs_{c}_summary.png", bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------- CS–CS pairs
def dup_fig(d):
    if not d["means"]:
        return
    ch = d["ch"]; n = len(ch)
    fig, axs = plt.subplots(int(np.ceil(n / 4)), 4, figsize=(14, 2.1 * np.ceil(n / 4)), sharex=True, sharey=True, squeeze=False)
    cols = {"coincident": "k", f"{d['cs_a']}-only": CB, f"{d['cs_b']}-only": "#3380cc"}
    for i, q in enumerate(ch):
        ax = axs.ravel()[i]
        for k, m in d["means"].items():
            ax.plot(T, m[:, i], c=cols[k], lw=1.1, label=f"{k} (n {d['n'][k]})")
        ax.set_title(f"ch {q} · x {chan_pos[q, 0]:.0f} y {chan_pos[q, 1]:.0f}", fontsize=7)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    axs.ravel()[0].legend(fontsize=6.5)
    fig.suptitle(f"CS {d['cs_a']} / {d['cs_b']} ({d['status']}, coinc {d['coinc']:.2f}, lag {d['lag_ms']:.2f} ms): raw means in "
                 f"{d['cs_a']}'s time frame; cosine of unit means {d['cos_full']:.2f} (after 1.5 ms: {d['cos_late']:.2f})", fontsize=9)
    fig.supxlabel("ms"); fig.tight_layout()
    fig.savefig(FIG / f"dup_{d['cs_a']}_{d['cs_b']}.png", bbox_inches="tight"); plt.close(fig)


def conditional_fig(c):
    cp = caps_of(c); cd = {d["unit"]: d for d in cond.get(c, [])}
    ks = [int(k) for k in cp.unit if int(k) in cd]
    if not ks:
        return
    n = len(ks); nc = 4; nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(15, 2.6 * nr), sharex=True, squeeze=False)
    pk = PEAK[c]
    for ax, k in zip(axs.ravel(), ks):
        d = cd[k]; r = cp[cp.unit == k].iloc[0]
        i = d["ch"].index(pk) if pk in d["ch"] else 0
        ax.plot(T, d["without"][:, i], c="0.5", lw=1.2); ax.plot(T, d["with_"][:, i], c="#d1495b", lw=1.3)
        diff = d["with_"][:, i] - d["without"][:, i]
        ax.plot(T, diff - np.abs(d["without"][:, i]).max() * .6, c="k", lw=.8)
        ax.axvline(r.lag_med_ms, c="#d1495b", ls="--", lw=.7)
        ax.set_title(f"{k}: {r.frac_cs_events:.0%} of events (chance {r.chance:.0%}) · lag {r.lag_med_ms:.2f} ms · "
                     f"{r.dist_um:.0f} µm · {r.fr:.0f} Hz", fontsize=7)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(f"CS {c}, peak ch {pk}: events WITH (red) vs WITHOUT (grey) a spike of each capturing cluster; black (shifted down) = difference", fontsize=10)
    fig.supxlabel("ms after CS spike"); fig.tight_layout()
    fig.savefig(FIG / f"cs_{c}_conditional.png", bbox_inches="tight"); plt.close(fig)


def locked_fig():
    L = sorted([d for d in locked if d["frac_locked"] > .01], key=lambda d: -d["frac_locked"])
    if not L:
        return
    nc = 4; nr = int(np.ceil(len(L) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(15, 2.4 * nr), squeeze=False)
    for ax, d in zip(axs.ravel(), L):
        i = int(np.argmax(np.ptp(d["other"], 0)))
        ax.plot(T, d["other"][:, i], c="0.45", lw=1.2, label=f"other (n {d['n_other']})")
        ax.plot(T, d["locked"][:, i], c="#d1495b", lw=1.3, label=f"CS-locked (n {d['n_locked']})")
        u = d["unit"]; role = units.at[u, "role"] if isinstance(units.at[u, "role"], str) else ""
        ax.set_title(f"{u} {role} {units.at[u, 'group']}: {d['frac_locked']:.1%} of spikes CS-locked", fontsize=7.5)
        ax.legend(fontsize=6)
    for ax in axs.ravel()[len(L):]:
        ax.axis("off")
    fig.suptitle("Clusters capturing CS events: their CS-locked spikes vs their other spikes (own peak channel, own spike time). "
                 "Different shapes → the locked spikes are CS pieces, not the unit's own spikes.", fontsize=9)
    fig.supxlabel("ms after the unit's spike"); fig.tight_layout(rect=(0, 0, 1, 1 - 0.5 / (2.4 * nr)))
    fig.savefig(FIG / "capturing_clusters_locked_vs_other.png", bbox_inches="tight"); plt.close(fig)


session_summary(); locked_fig()
for c in CS_IDS:
    events_fig(c); summary_fig(c); conditional_fig(c)
for d in dups:
    dup_fig(d)
print("figures in", FIG, len(list(FIG.glob("*.png"))))
