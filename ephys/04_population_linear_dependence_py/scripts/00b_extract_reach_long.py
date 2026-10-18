# %% ── 0. imports & parameters ────────────────────────────────────────────────
# 00b -- REACH_LONG: kinematics of every analysed reach in a LONG window (default the
# spike window, -3..+3 s around max reach), cut from the continuous session tracking.
#
# How the current 201-frame reaches were made (MATLAB, behavior/paw_and_body_trajectories):
#   step1  per trial tt: ppDom_paw(:, :, tt) = cleaned/interpolated dominant-paw x,y,z (px);
#          peaks found TRIAL BY TRIAL (findpeaks on x); reach window = 150 frames before /
#          50 after the peak, borrowing frames from the previous / next trial at the edges
#          -> reaches.reach_mat, reaches.reach_trial (video trial), reaches.reach_range_mat
#          (frame-in-trial; row 151 = the max-reach frame).
#          session.ppDom_paw_sess = the trials concatenated (continuous video, 400 fps).
#   step3  keeps the reaches of the analysed trials: bhv.reach_inVec_idx -> bhv.reaches_inVec_px
#   s05    st_eg_neurons.mat: reach_px = bhv.reaches_inVec_px (what the MI loader reads).
# Here: global frame of reach i = (frames of trials before reach_trial(i)) + reach_range(151, i);
# REACH_LONG = ppDom_paw_sess around that frame, for the reaches with reach_inVec_idx.
#
# Validation (per session, printed + figure): the long window, cut back to -375..+125 ms,
# must reproduce reaches.reach_mat (all reaches) AND the MI reach_px (selected reaches)
# exactly. Only frames outside the session are NaN; frames are never altered. Timing caveats are
# FLAGGED per reach: n_frames_ts_irregular (timestamps not 2.5 ms apart, e.g. Lindt R3 trials
# 79/81 step 6.8 ms -- the current data has the same frames) and trial_boundary_s (below).
# CAVEAT, trial boundaries: video trials are concatenated exactly as step1 did (so the
# common window matches), but at each boundary the next trial's first frame has the SAME
# timestamp as the previous frame and the position jumps ~2 normal steps -> timing there is
# uncertain by ~1-2 frames (2.5-5 ms). Their positions per reach: trial_boundary_s.
# (ppDom_paw_sess starts at step1's start_trials, e.g. CoteDor R2 skips the preparation trial.)
#
# Saved as ONE file, like the MI reach_px: <VAR_DIR>/reach_long/reach_long_win-3000to3000.pkl
#   REACH_LONG[(mouse, sess)] = (n_reaches, 3, n_frames) px  -- same reaches and order as REACH_PX
#   INFO   = table, one row per reach: mouse, sess, trial_idx (= row in REACH_PX / io trial_idx),
#            video_trial, prev/next reach (s), timing flags
#   SESSIONS = per-session validation and sources;  t, fps, px2mm, ...
# Load anywhere with:  D = io.load_reach_long();  D["REACH_LONG"][(mouse, sess)], D["INFO"]
import os, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

_here = globals().get("__vsc_ipynb_file__") or globals().get("__file__")
_starts = ([Path(_here).resolve().parent] if _here else []) + [Path.cwd()]
_cands = [q for s in _starts for p in [s, *s.parents]
          for q in (p, p / "04_population_linear_dependence_py",
                    p / "ephys" / "04_population_linear_dependence_py")]
ROOT = next((p for p in _cands if (p / "gctools").is_dir()), None)
if ROOT is None:
    raise RuntimeError(f"gctools/ not found from {_starts}")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import importlib
import h5py
import scipy.io as sio
from gctools import config, io, cache, plotting
for _m in (config, io, cache, plotting):
    importlib.reload(_m)

WINDOW_S   = (-3.0, 3.0)        # re max reach; the spike window is config.SPIKE_SPAN_S
SESSIONS   = "all"              # "all" = every session in the MI data, or list of (mouse, sess)
OVERWRITE  = False              # False: sessions already in the saved file are kept, not recomputed
SAVE_VARS  = True
SAVE_FIGS  = True
SHOW_FIGS  = True
N_EXAMPLES = 5                  # random reaches per mouse in the sanity figure
SEED       = 0

FPS     = config.FPS
K       = np.arange(int(round(WINDOW_S[0] * FPS)), int(round(WINDOW_S[1] * FPS)) + 1)   # frames re max reach
T_LONG  = K / FPS                                                                       # s
SHORT   = slice(int(np.flatnonzero(K == -config.ALIGN_FRAME)[0]),
                int(np.flatnonzero(K == config.N_FRAMES - 1 - config.ALIGN_FRAME)[0]) + 1)  # -375..+125 ms
_, REACH_PX = io.load_mi()
SESS_LIST = sorted(REACH_PX) if SESSIONS == "all" else [tuple(s) for s in SESSIONS]
FIG_DIR = config.fig_dir("reach_long") if SAVE_FIGS else None
print(f"window {WINDOW_S} s = {K.size} frames; {len(SESS_LIST)} sessions")


# %% ── 1. extract (or load) every session ───────────────────────────────────
def _bhv_selection(mouse, sess):
    """bhv.reach_inVec_idx (step3): which reaches of reaches.reach_mat are analysed."""
    p = config.mi_config.bhv_fundamentals_mat(mouse, sess)
    try:
        b = sio.loadmat(str(p), squeeze_me=True, struct_as_record=False, variable_names=["bhv"])["bhv"]
        return np.asarray(b.reach_inVec_idx).astype(bool).ravel()
    except NotImplementedError:                                   # v7.3
        with h5py.File(p, "r") as f:
            return np.asarray(f["bhv/reach_inVec_idx"]).astype(bool).ravel()


def extract_session(mouse, sess):
    src = config.mi_config.behav_session_mat(mouse, sess)
    with h5py.File(src, "r") as f:
        cont = np.asarray(f["session/ppDom_paw_sess"], float)          # (3, n_frames) px
        siz  = np.asarray(f["trials/size_trial"]).ravel().astype(int)
        TS   = np.asarray(f["trials/sync/frames_timestamps"], float)  # (frame, trial)
        RM   = np.asarray(f["reaches/reach_mat"], float)               # (n_reach, 3, 201)
        rtr  = np.asarray(f["reaches/reach_trial"]).ravel().astype(int)          # video trial, 1-idx
        rmax = np.asarray(f["reaches/reach_range_mat"])[config.ALIGN_FRAME].astype(int)  # frame in trial, 1-idx
        calib = (float(f["session/video/calib/px2mm/px2mm_xz"][0, 0]),
                 float(f["session/video/calib/px2mm/px2mm_y"][0, 0]),
                 float(f["session/video/calib/px2mm/px2mm_xz"][0, 0]))
    # ppDom_paw_sess = trials start_trials..end concatenated (step1); find start_trials
    st = np.flatnonzero(np.cumsum(siz[::-1])[::-1] == cont.shape[1])
    if st.size != 1:
        raise ValueError(f"{mouse}/{sess}: cannot place ppDom_paw_sess ({cont.shape[1]} frames) in the trials")
    st = int(st[0])                                                     # 0-idx first trial in the trace
    if (rtr - 1 < st).any():
        raise ValueError(f"{mouse}/{sess}: reaches in trials before start_trials={st+1}")
    ts = np.concatenate([TS[:siz[t], t] for t in range(st, siz.size)])  # harp time of every frame
    first = np.concatenate([[0], np.cumsum(siz[st:])])                  # trace index of each trial's frame 1
    gframe = first[rtr - 1 - st] + rmax - 1                             # 0-idx global max-reach frame
    tb = first[1:-1]                                                    # trace index where a new trial starts

    def cut(g):
        """Frames g+K of the trace, concatenated exactly as step1 did (frames are NOT
        altered, so the common window reproduces reach_mat). Timing caveats are FLAGGED:
          trial boundaries: the next trial's first frame has the SAME timestamp as the
            previous frame (+-1-2 frames uncertainty) -> trial_boundary_s
          irregular timestamps (a step != 1 frame +- 0.5, other than those boundaries; e.g.
            whole long trials with 6.8 ms steps) -> n_frames_ts_irregular
        Only frames outside the session are NaN."""
        idx = g + K
        ok = (idx >= 0) & (idx < cont.shape[1])
        out = np.full((3, K.size), np.nan)
        out[:, ok] = cont[:, idx[ok]]
        step = np.diff(ts[idx[ok]])
        irregular = (np.abs(step - 1 / FPS) > 0.5 / FPS) & (np.abs(step) > 0.1 / FPS)
        kb = K[np.isin(idx, tb)] / FPS                                  # s re max reach
        return out, int((~ok).sum()), int(irregular.sum()), kb

    # validation 1: every reach of reach_mat, cut back to -375..+125 ms
    d_rm = max(np.nanmax(np.abs(cut(g)[0][:, SHORT] - RM[i])) for i, g in enumerate(gframe))

    sel = _bhv_selection(mouse, sess)
    if sel.size != gframe.size:
        raise ValueError(f"{mouse}/{sess}: reach_inVec_idx has {sel.size} entries, reach_mat {gframe.size}")
    gs = gframe[sel]
    cuts = [cut(g) for g in gs]
    X = np.stack([c[0] for c in cuts])                                   # (n_sel, 3, n_k) px
    # validation 2: selected reaches vs the MI reach_px (what every script uses now)
    rp = REACH_PX[(mouse, sess)]
    d_mi = (np.nanmax(np.abs(X[:, :, SHORT] - rp)) if rp.shape == X[:, :, SHORT].shape else np.inf)
    nan_mismatch = (int((np.isnan(X[:, :, SHORT]) != np.isnan(rp)).sum())
                    if rp.shape == X[:, :, SHORT].shape else -1)

    # other reaches nearby (any detected reach of the session, analysed or not)
    allg = np.sort(gframe)
    pos = np.searchsorted(allg, gs)
    prev_s = np.array([(g - allg[p - 1]) / FPS if p > 0 else np.inf for g, p in zip(gs, pos)])
    nxt = np.searchsorted(allg, gs, side="right")
    next_s = np.array([(allg[n] - g) / FPS if n < allg.size else np.inf for g, n in zip(gs, nxt)])

    return dict(mouse=mouse, sess=sess, window_s=WINDOW_S, t=T_LONG, fps=FPS, align_index=int(np.flatnonzero(K == 0)[0]),
                px=X.astype(np.float32), coords=("x", "y", "z"), px2mm=tuple(config.PX2MM), px2mm_session_calib=calib,
                global_frame=gs, video_trial=rtr[sel], reach_index=np.flatnonzero(sel),
                prev_reach_s=prev_s, next_reach_s=next_s,
                n_frames_outside_session=np.array([c[1] for c in cuts]),
                n_frames_ts_irregular=np.array([c[2] for c in cuts]),  # timestamps not 2.5 ms apart (frames kept)
                trial_boundary_s=[c[3] for c in cuts],    # per reach: where a new video trial starts
                start_trials=st + 1,
                validation=dict(max_abs_diff_vs_reach_mat_px=float(d_rm), max_abs_diff_vs_MI_reach_px=float(d_mi),
                                nan_pattern_mismatch_vs_MI=nan_mismatch, n_reach_mat=int(gframe.size),
                                n_selected=int(sel.sum()), n_MI=int(rp.shape[0])),
                source=dict(behav_mat=str(src), behav_mtime=src.stat().st_mtime,
                            bhv_mat=str(config.mi_config.bhv_fundamentals_mat(mouse, sess))))

PATH = io.reach_long_path(WINDOW_S)
OLD = (cache.load_pickle(PATH) if os.path.exists(cache._long(PATH)) and not OVERWRITE else None)
REACH_LONG, INFO, SESS_INFO, t0 = {}, [], {}, time.time()
for mouse, sess in SESS_LIST:
    key = (mouse, sess)
    if OLD is not None and key in OLD["REACH_LONG"]:                 # already extracted
        REACH_LONG[key] = OLD["REACH_LONG"][key]
        INFO.append(OLD["INFO"][(OLD["INFO"].mouse == mouse) & (OLD["INFO"].sess == sess)])
        SESS_INFO[key] = OLD["SESSIONS"][key]
        print(f"  {mouse}/{sess}: loaded from {PATH.name}")
        continue
    try:
        R = extract_session(mouse, sess)
    except (OSError, KeyError, ValueError) as e:
        print(f"  {mouse}/{sess}: FAILED -- {e}")
        continue
    REACH_LONG[key] = R["px"]
    n = R["px"].shape[0]
    INFO.append(pd.DataFrame(dict(
        mouse=mouse, sess=sess, trial_idx=np.arange(n), video_trial=R["video_trial"],
        reach_index=R["reach_index"], global_frame=R["global_frame"],
        prev_reach_s=R["prev_reach_s"], next_reach_s=R["next_reach_s"],
        n_frames_outside_session=R["n_frames_outside_session"],
        n_frames_ts_irregular=R["n_frames_ts_irregular"],
        n_trial_boundaries=[len(b) for b in R["trial_boundary_s"]],
        trial_boundary_s=R["trial_boundary_s"])))
    SESS_INFO[key] = {k: R[k] for k in ("validation", "start_trials", "px2mm_session_calib", "source")}
    v = R["validation"]
    print(f"  {mouse}/{sess}: {v['n_selected']}/{v['n_reach_mat']} reaches, max |diff| vs reach_mat "
          f"{v['max_abs_diff_vs_reach_mat_px']:.1e} px, vs MI {v['max_abs_diff_vs_MI_reach_px']:.1e} px  "
          f"({time.time()-t0:.0f}s)")
INFO = pd.concat(INFO, ignore_index=True)
if OLD is not None:                                                  # keep sessions not in SESS_LIST
    for key in OLD["REACH_LONG"]:
        if key not in REACH_LONG:
            REACH_LONG[key] = OLD["REACH_LONG"][key]; SESS_INFO[key] = OLD["SESSIONS"][key]
            INFO = pd.concat([INFO, OLD["INFO"][(OLD["INFO"].mouse == key[0]) & (OLD["INFO"].sess == key[1])]],
                             ignore_index=True)
if SAVE_VARS:
    cache.save_pickle(dict(REACH_LONG=REACH_LONG, INFO=INFO, SESSIONS=SESS_INFO, t=T_LONG, fps=FPS,
                           window_s=WINDOW_S, align_index=int(np.flatnonzero(K == 0)[0]),
                           coords=("x", "y", "z"), unit="px", px2mm=tuple(config.PX2MM)), PATH)
    print(f"saved {len(REACH_LONG)} sessions, {len(INFO)} reaches -> {PATH}")


# %% ── 2. validation table ───────────────────────────────────────────────────
# Exact reproduction = 0 (validated in float64 before float32 storage). Outside session = NaN.
rows = []
for (mouse, sess) in REACH_LONG:
    v = SESS_INFO[(mouse, sess)]["validation"]
    I = INFO[(INFO.mouse == mouse) & (INFO.sess == sess)]
    rows.append(dict(mouse=mouse, sess=sess, n_reaches=v["n_selected"], n_MI=v["n_MI"],
                     diff_vs_reach_mat_px=v["max_abs_diff_vs_reach_mat_px"],
                     diff_vs_MI_px=v["max_abs_diff_vs_MI_reach_px"],
                     nan_mismatch_vs_MI=v["nan_pattern_mismatch_vs_MI"],
                     reaches_hitting_session_edge=int((I.n_frames_outside_session > 0).sum()),
                     reaches_ts_irregular=int((I.n_frames_ts_irregular > 0).sum()),
                     trial_boundaries_per_window=I.n_trial_boundaries.mean().round(2),
                     short_window_crosses_boundary=int(sum(((np.asarray(b) >= -0.375) & (np.asarray(b) <= 0.125)).any()
                                                           for b in I.trial_boundary_s)),
                     prev_reach_lt_1s=np.mean(I.prev_reach_s < 1.0).round(2),
                     prev_reach_in_window=np.mean(I.prev_reach_s < -WINDOW_S[0]).round(2)))
VAL = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print(VAL.to_string(index=False))
ok = (VAL.diff_vs_reach_mat_px < 1e-6) & (VAL.diff_vs_MI_px < 1e-6) & (VAL.nan_mismatch_vs_MI == 0)
print(f"\n{ok.sum()}/{len(VAL)} sessions reproduce reach_mat AND the MI reach_px exactly"
      + ("" if ok.all() else f"  -- CHECK: {VAL.loc[~ok, ['mouse', 'sess']].values.tolist()}"))
print("prev_reach_in_window = fraction of reaches with another reach's max within the window before "
      "max reach (the long window then contains that reach too)")


# %% ── 3. sanity figure per mouse: long vs short window ──────────────────────
# Left: random reaches, long window (colour) with the current 201-frame reach on top (black
# dashed). Middle: mean over all reaches of the mouse, long vs short. Right: long - short
# mean in the common window (should be 0). Shaded = current window (-375..+125 ms).
rng = np.random.default_rng(SEED)
LBL = {0: "x (mm)", 1: "y = med-lat (mm)", 2: "z (mm)"}
for mouse in sorted({m for m, _ in REACH_LONG}):
    keys = [k for k in REACH_LONG if k[0] == mouse]
    L = np.concatenate([REACH_LONG[k].astype(float) for k in keys])           # (n, 3, n_k)
    S = np.concatenate([REACH_PX[k] for k in keys])                                   # (n, 3, 201)
    mm = np.asarray(config.PX2MM)[None, :, None]
    L, S = L * mm, S * mm
    ex = rng.choice(L.shape[0], min(N_EXAMPLES, L.shape[0]), replace=False)
    ts_ = T_LONG[SHORT]
    fig, axs = plt.subplots(3, 3, figsize=(16, 9), gridspec_kw=dict(width_ratios=[2, 2, 1]))
    for c in range(3):
        ax = axs[c, 0]
        for j, i in enumerate(ex):
            ax.plot(T_LONG, L[i, c], lw=0.8, color=f"C{j}", label=f"reach {i}" if c == 0 else None)
            ax.plot(ts_, S[i, c], color="k", ls="--", lw=1.2, label="current 201 frames" if (c == 0 and j == 0) else None)
        ax = axs[c, 1]
        ax.plot(T_LONG, np.nanmean(L[:, c], 0), color="tab:blue", lw=2, label="long window (mean)")
        ax.plot(ts_, np.nanmean(S[:, c], 0), color="k", ls="--", lw=1.5, label="current (mean)")
        ax = axs[c, 2]
        ax.plot(ts_ * 1000, np.nanmean(L[:, c, SHORT], 0) - np.nanmean(S[:, c], 0), color="tab:red")
        ax.axhline(0, color="k", lw=0.5)
        ax.set(xlabel="time re max reach (ms)", ylabel="long − short (mm)")
        for ax in axs[c, :2]:
            ax.axvspan(ts_[0], ts_[-1], color="0.9", zorder=-2)
            ax.axvline(0, color="k", ls=":", lw=0.8)
            ax.set(ylabel=LBL[c], xlim=WINDOW_S)
        axs[c, 0].set_xlabel("time re max reach (s)"); axs[c, 1].set_xlabel("time re max reach (s)")
    axs[0, 0].legend(fontsize=7, frameon=False, ncol=2)
    axs[0, 1].legend(fontsize=7, frameon=False)
    axs[0, 0].set_title(f"{N_EXAMPLES} random reaches"); axs[0, 1].set_title(f"mean of all {L.shape[0]} reaches")
    axs[0, 2].set_title("difference of the means\n(common window; should be 0)")
    fig.suptitle(f"{mouse}: REACH_LONG {WINDOW_S} s vs current reach_px  |  sessions "
                 f"{', '.join(s for _, s in keys)}  |  max |diff| "
                 f"{VAL[VAL.mouse == mouse].diff_vs_MI_px.max():.1e} px", fontsize=10)
    plotting.finish(fig, FIG_DIR, f"reach_long_vs_short_{mouse}", SAVE_FIGS, SHOW_FIGS)

# %%
