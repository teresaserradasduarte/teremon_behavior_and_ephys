# %% ── 0. imports & constants ─────────────────────────────────────────────────
# Build the WATER-REFERENCED behavioural features, per reach:
#   (the paw path itself is NOT stored: it is exactly REACH_PX * PX2MM, already in
#    the parsed cache, and duplicating it cost ~140 MB. 09 derives REACH_MM itself.)
#   water_mm      (n_reach, 3)       water position for that reach's TRIAL
#   dist_water    (n_reach, 201)     Euclidean |paw - water| in mm
#   speed_water   (n_reach, 201)     d(dist)/dt in mm/s (NEGATIVE = approaching)
#   lcr / cat_reach (n_reach,)       target 1/2/3 = LEFT/CENTER/RIGHT, reach category
#
# LCR CODE VERIFIED against bhv.left_inVec / center_inVec / right_inVec:
# code 1 -> left (830 matches, 1 stray), 2 -> center (1035, 0), 3 -> right (1242, 5).
#
# Two source trees have to be joined, and they are DIFFERENT group folders:
#   config.BEHAV_MAT_DIR/<mouse>/<sess>/session_reaching_data_paw.mat  (v7.3, h5py)
#       trials.water                -> (ntime, 3, ntrials) in MATLAB order
#       trials.water_location_trial -> (ntrials,) target identity 1/2/3
#   config.MAT_DIR/<mouse>/<sess>/behavior_fundamentals.mat            (v7, scipy)
#       bhv.reaches_inVec_px    -> (201, 3, n_reach); REACH_PX is its (2,1,0) transpose
#       bhv.reach_trials_inVec  -> (n_reach,) which TRIAL each reach belongs to
#       bhv.LCR_inVec_idx       -> (n_reach,) target identity per reach
#
# INDEX CONVENTION (measured, not assumed -- see cell 2): reach_trials_inVec
# indexes the water array DIRECTLY, with NO -1 offset. Checked on all 29 sessions
# by matching LCR_inVec_idx against water_location_trial: direct indexing scores
# 99.2-100%, the -1 offset only 88.4-97.9%.
import numpy as np
import matplotlib.pyplot as plt
import h5py, scipy.io as sio
import pickle, time, importlib
import config
importlib.reload(config)
from scipy.signal import savgol_filter

FPS         = 400.0            # camera frame rate
PX2MM       = np.asarray(config.PX2MM, float)   # (0.1, 0.07, 0.1)
SAVGOL_WIN  = 11               # same derivative filter as 09
SAVGOL_ORD  = 3
FRAMES      = (147, 167)       # endpoint window, MATLAB 1-indexed (same as 04/08/09)
EP          = slice(FRAMES[0]-1, FRAMES[1])

MIN_WATER_FRAMES = 20          # a trial needs this many valid water samples
MAX_WATER_MAD    = 15.0        # px; above this the "stable" water pos is not stable

OUT_PKL     = "behavior_water.pkl"
OUT_SUBDIR  = "behavior_water"
SHOW_FIGS   = True
SAVE_FIGS   = True
SAVE_VARS   = True

print(f"behaviour tree : {config.BEHAV_MAT_DIR}")
print(f"ephys tree     : {config.MAT_DIR}")
print(f"px -> mm       : {PX2MM}")


# %% ── 1. helpers ────────────────────────────────────────────────────────────
def water_median_px(behav_mat):
    """-> (water_med (ntrials,3), water_mad (ntrials,3), n_valid (ntrials,3),
            water_loc (ntrials,)) for one session, in PIXELS.

    trials.water is (ntime, 3, ntrials) in MATLAB and arrives from h5py as
    (ntrials, 3, ntime). It is NaN wherever the drop is absent (before delivery,
    after it is taken, or tracking failure), so a per-coordinate nanmedian over
    the whole trial recovers the stable location. Median (not mean) because the
    frames around delivery/collection are exactly the outliers we must ignore.
    NOTE the three coordinates do NOT share a validity mask -- on 1_CoteDor/R1
    30504 frames are valid in x but not y -- so each coordinate is reduced
    independently rather than through a common all-finite mask."""
    with h5py.File(behav_mat, "r") as f:
        t = f["trials"]
        w = t["water"][:]                                  # (ntrials, 3, ntime)
        loc = np.array(t["water_location_trial"]).ravel()
    n_valid = np.sum(~np.isnan(w), axis=2).astype(int)      # (ntrials, 3)
    with np.errstate(all="ignore"):
        med = np.nanmedian(w, axis=2)                       # (ntrials, 3)
        mad = np.nanmedian(np.abs(w - med[:, :, None]), axis=2)
    return med, mad, n_valid, loc


def load_bhv(fund_mat):
    """behavior_fundamentals.mat -> the few bhv fields we need."""
    b = sio.loadmat(fund_mat, squeeze_me=True, struct_as_record=False)["bhv"]
    return dict(
        px    = np.asarray(b.reaches_inVec_px, float).transpose(2, 1, 0),  # (n,3,201)
        trial = np.asarray(b.reach_trials_inVec).ravel().astype(int),
        lcr   = np.asarray(b.LCR_inVec_idx).ravel().astype(int),
        cat   = np.asarray(b.cat_reach_inVec).ravel().astype(int),
    )


def dist_and_speed(reach_mm, water_mm):
    """reach_mm (n,3,T), water_mm (n,3) -> dist (n,T) mm, speed (n,T) mm/s.

    Water is fixed within a trial, so the distance is just the norm of the paw
    position relative to it. The derivative is SIGNED: negative = closing on the
    water. Same Savitzky-Golay settings as 09 so the two scripts' derivatives are
    directly comparable."""
    d = np.linalg.norm(reach_mm - water_mm[:, :, None], axis=1)     # (n, T)
    # savgol's edge polyfit raises on any row containing NaN (and a reach whose
    # trial has no valid water is entirely NaN), so filter only the clean rows
    # and leave the rest NaN rather than letting one bad reach kill the session.
    v  = np.full_like(d, np.nan)
    ok = np.isfinite(d).all(axis=1)
    if ok.any():
        v[ok] = savgol_filter(d[ok], SAVGOL_WIN, SAVGOL_ORD, deriv=1,
                              delta=1/FPS, axis=1)
    return d, v


# %% ── 2. build every session ────────────────────────────────────────────────
# Walks the ephys tree (that is what defines which sessions exist for the MI
# analysis) and joins each one to its raw-behaviour twin.
SESSIONS, REPORT = {}, []
t0 = time.time()
# only real mouse folders: MAT_DIR also holds MutualInformation_py etc.
mice = sorted(p.name for p in config.MAT_DIR.iterdir()
              if p.is_dir() and any(a in p.name for a in config.ANIMALS))
for mouse in mice:
    for sdir in sorted((config.MAT_DIR / mouse).iterdir()):
        if not sdir.is_dir(): continue
        sess = sdir.name
        fund, behav = config.bhv_fundamentals_mat(mouse, sess), \
                      config.behav_session_mat(mouse, sess)
        if not (fund.is_file() and behav.is_file()):
            REPORT.append((mouse, sess, "MISSING FILE", np.nan, np.nan, np.nan)); continue

        B   = load_bhv(fund)
        med, mad, nval, loc = water_median_px(behav)
        n_reach = B["px"].shape[0]

        # --- trial -> water, with NO index offset (validated below) ---
        tr = B["trial"]
        in_range = (tr >= 0) & (tr < med.shape[0])
        w_px = np.full((n_reach, 3), np.nan)
        w_px[in_range] = med[tr[in_range]]

        # drop reaches whose trial has too few / too unstable water samples
        bad = np.zeros(n_reach, bool)
        bad[~in_range] = True
        bad[in_range] |= (nval[tr[in_range]].min(axis=1) < MIN_WATER_FRAMES)
        bad[in_range] |= (mad[tr[in_range]].max(axis=1) > MAX_WATER_MAD)
        w_px[bad] = np.nan

        # --- SANITY CHECK: does water identity track the LCR label? ---
        agree = np.nan
        if in_range.any():
            wl = loc[tr[in_range]]
            ok = ~np.isnan(wl)
            if ok.any(): agree = float(np.mean(wl[ok] == B["lcr"][in_range][ok]) * 100)

        reach_mm = B["px"] * PX2MM[None, :, None]
        w_mm     = w_px * PX2MM[None, :]
        dist, spd = dist_and_speed(reach_mm, w_mm)

        SESSIONS[(mouse, sess)] = dict(
            reach_mm=reach_mm, water_px=w_px, water_mm=w_mm,
            dist_water=dist, speed_water=spd,
            trial=tr, lcr=B["lcr"], cat_reach=B["cat"], water_loc=np.where(in_range, loc[np.clip(tr,0,med.shape[0]-1)], np.nan),
            n_valid_water=np.where(in_range[:,None], nval[np.clip(tr,0,med.shape[0]-1)], 0),
            water_mad=np.where(in_range[:,None], mad[np.clip(tr,0,med.shape[0]-1)], np.nan),
        )
        REPORT.append((mouse, sess, "ok", n_reach, agree, float(np.mean(bad)*100)))

print(f"\nbuilt {len(SESSIONS)} sessions in {time.time()-t0:.0f}s")
print(f"\n{'mouse/session':26s} {'status':12s} {'n reach':>8s} {'LCR agree':>10s} {'% dropped':>10s}")
for mouse, sess, st, n, ag, dr in REPORT:
    ns = f"{int(n)}" if np.isfinite(n) else "-"
    ags = f"{ag:9.1f}%" if np.isfinite(ag) else f"{'-':>10s}"
    drs = f"{dr:9.1f}%" if np.isfinite(dr) else f"{'-':>10s}"
    print(f"{mouse+'/'+sess:26s} {st:12s} {ns:>8s} {ags} {drs}")


# %% ── 3. cross-check against REACH_PX used by the MI scripts ────────────────
# The reach arrays must be IDENTICAL to what 04/08/09 analyse, or the water
# features would be attached to the wrong reaches.
import mi_io
if "REACH_PX" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())

n_match = n_diff = n_absent = 0
for key, S in SESSIONS.items():
    if key not in REACH_PX: n_absent += 1; continue
    a, b = REACH_PX[key], S["reach_mm"] / PX2MM[None, :, None]
    if a.shape == b.shape and np.allclose(a, b, equal_nan=True): n_match += 1
    else:
        n_diff += 1
        print(f"  MISMATCH {key}: REACH_PX {a.shape} vs reaches_inVec_px {b.shape}")
print(f"\nreach arrays: {n_match} identical, {n_diff} differ, "
      f"{n_absent} built here but absent from the parsed cache")
missing = [k for k in REACH_PX if k not in SESSIONS]
if missing: print(f"in REACH_PX but NOT built here: {missing}")


# %% ── 4. save ───────────────────────────────────────────────────────────────
if SAVE_VARS:
    # reach_mm stays in the in-memory SESSIONS (cells 3/5/7 use it) but is stripped
    # from the FILE: it is REACH_PX * PX2MM, which every consumer already has.
    _slim = {k: {kk: vv for kk, vv in S.items() if kk != "reach_mm"}
             for k, S in SESSIONS.items()}
    out = dict(sessions=_slim, PX2MM=tuple(PX2MM), FPS=FPS, FRAMES=FRAMES,
               SAVGOL_WIN=SAVGOL_WIN, SAVGOL_ORD=SAVGOL_ORD,
               MIN_WATER_FRAMES=MIN_WATER_FRAMES, MAX_WATER_MAD=MAX_WATER_MAD)
    p = config.var_path(OUT_PKL)
    with open(p, "wb") as f: pickle.dump(out, f, protocol=4)
    print(f"saved -> {p}  ({p.stat().st_size/1e6:.1f} MB)")


# %% ── 5. example figures ────────────────────────────────────────────────────
# Per example session: single reaches thin, the mean thick, water as a dashed
# horizontal line (it is fixed in time within a trial).
EX_KEYS    = [("1_CoteDor", "R1"), ("3_Toblerone", "R1"), ("5_FerreroRocher", "R1")]
N_EX_REACH = 25
fdir = config.fig_dir(OUT_SUBDIR)
tms  = (np.arange(201) - EP.start) / FPS * 1000
epw  = (EP.stop - EP.start) / FPS * 1000
rng  = np.random.default_rng(0)
COORD = ("paw x", "paw y (med-lat)", "paw z")

# Target identity 1/2/3 = LEFT/CENTER/RIGHT (verified against bhv.left_inVec etc.).
# Colours match the MATLAB reach figures: sand yellow / near-black brown / orange.
LOC_CLR  = {1: "#E5CE87", 2: "#3D352F", 3: "#F57F4F"}
LOC_NAME = {1: "left reaches", 2: "center reaches", 3: "right reaches"}
EX_CAT   = 1        # draw ONLY this reach category; other categories start from very
                    # different postures and swamp the panel with position spread.
                    # None -> use every category.

for key in EX_KEYS:
    if key not in SESSIONS: continue
    S  = SESSIONS[key]
    ok = np.isfinite(S["water_mm"]).all(axis=1)
    if EX_CAT is not None: ok = ok & (S["cat_reach"] == EX_CAT)
    if not ok.any(): continue
    locs = [l for l in (1, 2, 3) if np.any(ok & (S["water_loc"] == l))]

    fig, axes = plt.subplots(1, 5, figsize=(24, 4))
    for c in range(3):                       # the three coordinates, in mm
        ax = axes[c]
        for l in locs:
            m   = np.flatnonzero(ok & (S["water_loc"] == l))
            sel = rng.choice(m, min(N_EX_REACH, m.size), replace=False)
            ax.plot(tms, S["reach_mm"][sel, c, :].T, color=LOC_CLR[l], lw=0.5, alpha=0.25)
            ax.plot(tms, np.nanmean(S["reach_mm"][m, c, :], axis=0), color=LOC_CLR[l],
                    lw=2.5, label=f"{LOC_NAME[l]} (n={m.size})")
            # water is fixed in time within a trial -> a horizontal line per target
            ax.axhline(np.nanmedian(S["water_mm"][m, c]), color=LOC_CLR[l], lw=1.6,
                       ls="--", alpha=0.9)
        ax.axvspan(0, epw, color="k", alpha=0.12)
        ax.set(xlabel="time re endpoint-window start (ms)", ylabel=f"{COORD[c]} (mm)",
               title=f"{COORD[c]}   (dashed = water)")
        ax.legend(fontsize=7); ax.spines[["top","right"]].set_visible(False)

    ax = axes[3]                             # distance to water
    for l in locs:
        m   = np.flatnonzero(ok & (S["water_loc"] == l))
        sel = rng.choice(m, min(N_EX_REACH, m.size), replace=False)
        ax.plot(tms, S["dist_water"][sel].T, color=LOC_CLR[l], lw=0.5, alpha=0.25)
        ax.plot(tms, np.nanmean(S["dist_water"][m], axis=0), color=LOC_CLR[l], lw=2.5,
                label=f"{LOC_NAME[l]} (n={m.size})")
    ax.axhline(0, color="k", lw=1.4, ls="--", label="water (distance 0)")
    ax.axvspan(0, epw, color="k", alpha=0.12)
    ax.set(xlabel="time re endpoint-window start (ms)", ylabel="distance to water (mm)",
           title="distance to water")
    ax.legend(fontsize=7); ax.spines[["top","right"]].set_visible(False)

    ax = axes[4]                             # speed towards water
    for l in locs:
        m   = np.flatnonzero(ok & (S["water_loc"] == l))
        sel = rng.choice(m, min(N_EX_REACH, m.size), replace=False)
        ax.plot(tms, S["speed_water"][sel].T, color=LOC_CLR[l], lw=0.5, alpha=0.25)
        ax.plot(tms, np.nanmean(S["speed_water"][m], axis=0), color=LOC_CLR[l], lw=2.5,
                label=f"{LOC_NAME[l]} (n={m.size})")
    ax.axhline(0, color="k", lw=0.8)
    ax.axvspan(0, epw, color="k", alpha=0.12)
    # a handful of single trials hit +-3000 mm/s and flatten every mean curve, so
    # clip the VIEW (not the data) to a robust range around what is actually there
    _v = S["speed_water"][np.flatnonzero(ok)]
    _lim = np.nanpercentile(np.abs(_v), 99.0) * 1.15
    if np.isfinite(_lim) and _lim > 0: ax.set_ylim(-_lim, _lim)
    ax.set(xlabel="time re endpoint-window start (ms)", ylabel="d(distance)/dt (mm/s)",
           title="speed towards water (negative = approaching)")
    ax.legend(fontsize=7); ax.spines[["top","right"]].set_visible(False)

    _cat = f"cat_reach == {EX_CAT}" if EX_CAT is not None else "all categories"
    fig.suptitle(f"{key[0]} / {key[1]}  --  {int(ok.sum())} reaches ({_cat}) with a valid "
                 f"water position, split by target;  thin = single reaches "
                 f"(up to {N_EX_REACH} per target), thick = mean, dashed = water",
                 fontsize=12)
    plt.tight_layout()
    if SAVE_FIGS:
        fig.savefig(fdir / f"example_{key[0]}_{key[1]}.png", dpi=140, bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

print(f"figures -> {fdir}")


# %% ── 6. summary of the derived features ────────────────────────────────────
print(f"\n{'quantity':28s} {'median':>10s} {'IQR':>22s}")
for nm, arr in [("water x (mm)",   np.concatenate([S["water_mm"][:,0] for S in SESSIONS.values()])),
                ("water y (mm)",   np.concatenate([S["water_mm"][:,1] for S in SESSIONS.values()])),
                ("water z (mm)",   np.concatenate([S["water_mm"][:,2] for S in SESSIONS.values()])),
                ("dist @ endpoint (mm)",
                 np.concatenate([S["dist_water"][:,EP].mean(1) for S in SESSIONS.values()])),
                ("speed @ endpoint (mm/s)",
                 np.concatenate([S["speed_water"][:,EP].mean(1) for S in SESSIONS.values()]))]:
    a = arr[np.isfinite(arr)]
    q = np.percentile(a, [25, 75])
    print(f"{nm:28s} {np.median(a):10.2f}  [{q[0]:9.2f}, {q[1]:9.2f}]")

tot = sum(S["reach_mm"].shape[0] for S in SESSIONS.values())
val = sum(int(np.isfinite(S["water_mm"]).all(1).sum()) for S in SESSIONS.values())
print(f"\nreaches: {tot} total, {val} with a valid water position ({val/tot*100:.1f}%)")


# %% ── 7. what is in here? one ROW per variable ─────────────────────────────
# A table of names / shapes / dtypes -- double-click it in the Variables pane.
import pandas as pd

def vars_table(d):
    """dict -> DataFrame with one row per entry: name, shape, dtype, NaN %."""
    rows = []
    for k, v in d.items():
        a = np.asarray(v) if not isinstance(v, (dict, list, tuple)) else None
        if a is not None and a.dtype.kind in "fiu":
            nan = 100*float(np.mean(~np.isfinite(a))) if a.dtype.kind == "f" else 0.0
            rows.append((k, str(a.shape), str(a.dtype), a.size, round(nan, 1)))
        elif a is not None:
            rows.append((k, str(a.shape), str(a.dtype), a.size, np.nan))
        else:
            rows.append((k, f"{type(v).__name__}({len(v)})", "-", len(v), np.nan))
    return pd.DataFrame(rows, columns=["variable", "shape", "dtype", "size", "NaN %"])

def sess(mouse, session):
    """SESSIONS is keyed by the TUPLE (mouse, session), which is tedious to type.
        S = sess("1_CoteDor", "R1")   ==   SESSIONS[("1_CoteDor", "R1")]
    Accepts a partial mouse name too: sess("Lindt", "R3")."""
    if (mouse, session) in SESSIONS: return SESSIONS[(mouse, session)]
    hit = [k for k in SESSIONS if mouse in k[0] and k[1] == session]
    if len(hit) == 1: return SESSIONS[hit[0]]
    raise KeyError(f"{(mouse, session)} not found; candidates: {hit or list(SESSIONS)[:5]}")

VARS = vars_table(SESSIONS[list(SESSIONS)[0]])       # one session
print(f"variables per session (example {list(SESSIONS)[0]}):")
print(VARS.to_string(index=False))

# every session at a glance: n reaches, and how many have a water position
SESS_TABLE = pd.DataFrame([
    (m, s_, S["reach_mm"].shape[0], int(np.isfinite(S["water_mm"]).all(1).sum()))
    for (m, s_), S in SESSIONS.items()],
    columns=["mouse", "sess", "n_reach", "n_with_water"])
print(f"{len(SESS_TABLE)} sessions, {SESS_TABLE.n_reach.sum()} reaches "
      f"({SESS_TABLE.n_with_water.sum()} with water)")

# %%
