"""io.py - session loader. Reuses the MI project's parsed cache (mi_io.load on
mi_config.cache_path()); nothing is re-parsed from the .mat here.

One plain dict per session, IDENTICAL layout for real and simulated data:
  sess = dict(
    key=(mouse, sess), mouse, sess, source="real"|"sim", align_event="max_reach",
    spike_span=(t0, t1) s      # where spike times exist (history may not reach beyond)
    trials  = DataFrame        # one row per KEPT trial: trial_idx (original index),
                               #   LCR, PP, cat, hit, succ, cond (shuffle stratum label)
    kin     = dict(t=(n_frames,) s re align, fps, var, unit="mm", source="short"|"long",
                   X=(n_trials, n_frames))  # short: 201 frames, frame i at t=(i-ALIGN_FRAME)/FPS;
                                            # long: REACH_LONG, -3..+3 s (same reaches)
    neurons = [dict(uid, region, phyID, spikes=[array s re align] * n_trials)]
    truth   = None | dict      # simulations only
  )
t=0 is the MAX-REACH frame (see config.py), not reach onset.

    from gctools import io
    S = io.load_session("5_FerreroRocher", "R1"); io.describe(S)
"""
import numpy as np
import pandas as pd
from . import config

mi_io = config.mi_io

LABEL_MAP = {"idx_reach_LCR": "LCR", "idx_reach_PP": "PP", "idx_reach_cat": "cat",
             "idx_reach_hit": "hit", "idx_reach_succ": "succ"}

_MI = {}          # in-process memo of the MI cache: (neurons, reach_px)


# ---------------------------------------------------------------------------
# MI cache
# ---------------------------------------------------------------------------
def load_mi():
    """(NEURONS, REACH_PX) from the MI parsed cache, loaded once per process."""
    path = config.mi_config.cache_path()
    if str(path) not in _MI:
        _MI.clear()
        _MI[str(path)] = mi_io.load(path)
    return _MI[str(path)]

def reach_long_path(window_s=(-3.0, 3.0)):
    """The single REACH_LONG file written by scripts/00b_extract_reach_long.py."""
    w = f"win{int(round(window_s[0] * 1000))}to{int(round(window_s[1] * 1000))}"
    return config.var_path("reach_long", f"reach_long_{w}.pkl")


def load_reach_long(window_s=(-3.0, 3.0)):
    """REACH_LONG (made by scripts/00b_extract_reach_long.py), saved like the MI reach_px:
    dict(REACH_LONG={(mouse, sess): (n_reaches, 3, n_frames) px, float32 -- same reaches and
    order as REACH_PX; mm = px * px2mm}, INFO=table one row per reach (mouse, sess, trial_idx,
    prev/next reach, timing flags), SESSIONS=per-session validation, t (s re max reach), ...)."""
    from .cache import load_pickle
    return load_pickle(reach_long_path(window_s))


_RL = {}

def _reach_long():
    """REACH_LONG file, loaded once per process."""
    p = reach_long_path()
    if str(p) not in _RL:
        _RL.clear()
        _RL[str(p)] = load_reach_long()
    return _RL[str(p)]


def mi_stamp():
    """Identity of the MI cache (path + mtime) -- put into every real-data cache key so
    re-parsing on the MI side invalidates our caches."""
    p = config.mi_config.cache_path()
    return dict(mi_cache=str(p), mi_mtime=p.stat().st_mtime)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _labels_equal(a, b):
    return all(np.array_equal(a[k], b[k], equal_nan=True) for k in LABEL_MAP)

def cond_labels(trials, shuffle_within):
    """Stratum label per trial, e.g. "LCR=1" or "LCR=1|PP=0"; "all" if None."""
    if not shuffle_within:
        return np.full(len(trials), "all", dtype=object)
    cols = list(shuffle_within)
    return np.array(["|".join(f"{c}={int(v)}" if np.isfinite(v) else f"{c}=nan"
                              for c, v in zip(cols, row))
                     for row in trials[cols].to_numpy(float)], dtype=object)

def kin_time():
    """(N_FRAMES,) frame times in s re the align event."""
    return (np.arange(config.N_FRAMES) - config.ALIGN_FRAME) / config.FPS


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------
def load_session(mouse, sess, regions=None, kin_var=None, reach_cat=None,
                 shuffle_within="default", verbose=True, kin_source=None):
    """Session dict (see module docstring) for real data.
    Keeps trials with idx_reach_cat in reach_cat, finite kinematics IN THE 201-FRAME WINDOW,
    and a finite stratum label. Checks that trial labels are identical across all neurons.
    kin_source "short" (MI reach_px) | "long" (REACH_LONG: same trials kept, longer history;
    adds the REACH_LONG timing flags to `trials`). None = config.KIN_SOURCE."""
    regions   = list(regions or config.REGIONS)
    kin_var   = kin_var or config.KIN_VAR
    kin_source = kin_source or config.KIN_SOURCE
    reach_cat = tuple(reach_cat or config.REACH_CAT)
    if shuffle_within == "default":
        shuffle_within = config.SHUFFLE_WITHIN

    NEURONS, REACH_PX = load_mi()
    sel = [n for n in NEURONS if n["mouse"] == mouse and n["sess"] == sess
           and n["region"] in regions]
    if not sel:
        raise ValueError(f"no neurons for {mouse}/{sess} in regions {regions}")
    lab0 = sel[0]["labels"]
    bad = [n["uid"] for n in sel if not _labels_equal(n["labels"], lab0)]
    if bad:
        raise ValueError(f"trial labels differ across neurons of {mouse}/{sess}: {bad[:3]}")
    n_all = sel[0]["n_trials"]
    if any(n["n_trials"] != n_all for n in sel):
        raise ValueError(f"n_trials differs across neurons of {mouse}/{sess}")

    trials = pd.DataFrame({LABEL_MAP[k]: np.asarray(lab0[k], float) for k in LABEL_MAP})
    trials.insert(0, "trial_idx", np.arange(n_all))

    rp = REACH_PX[(mouse, sess)]                                 # (n, 3, 201) px
    c  = config.KIN_COORD[kin_var]
    X  = rp[:, c, :] * config.PX2MM[c]                           # mm
    kin_t = kin_time()
    if kin_source == "long":
        L = _reach_long()
        XL = np.asarray(L["REACH_LONG"][(mouse, sess)], float)[:, c, :] * config.PX2MM[c]
        if XL.shape[0] != n_all:
            raise ValueError(f"REACH_LONG has {XL.shape[0]} reaches for {mouse}/{sess}, MI {n_all}")
        info = L["INFO"][(L["INFO"].mouse == mouse) & (L["INFO"].sess == sess)].sort_values("trial_idx")
        for col in ("n_frames_ts_irregular", "n_trial_boundaries", "prev_reach_s", "next_reach_s"):
            trials[col] = info[col].to_numpy()
        kin_t = np.asarray(L["t"], float)
        i0 = int(L["align_index"]) - config.ALIGN_FRAME             # 201-frame window inside REACH_LONG
        XL[:, i0:i0 + config.N_FRAMES] = X                          # exact MI values there (file is float32)
    elif kin_source != "short":
        raise ValueError(f"kin_source must be 'short' or 'long', got {kin_source!r}")

    in_cat = np.isin(trials["cat"].to_numpy(), reach_cat)
    fin_x  = np.isfinite(X).all(axis=1)                           # the 201-frame window
    if kin_source == "long":
        X = XL
    fin_c  = (np.isfinite(trials[list(shuffle_within)].to_numpy(float)).all(axis=1)
              if shuffle_within else np.ones(n_all, bool))
    keep = in_cat & fin_x & fin_c
    if verbose:
        print(f"{mouse}/{sess}: {n_all} trials -> keep {keep.sum()}  "
              f"(cat {reach_cat}: {in_cat.sum()}, non-finite kin dropped: "
              f"{(in_cat & ~fin_x).sum()}, non-finite stratum dropped: "
              f"{(in_cat & fin_x & ~fin_c).sum()})")

    trials = trials[keep].reset_index(drop=True)
    trials["cond"] = cond_labels(trials, shuffle_within)
    ki = np.flatnonzero(keep)
    neurons = [dict(uid=n["uid"], region=n["region"], phyID=n["phyID"],
                    spikes=[np.asarray(n["spikes"][i], float) for i in ki])
               for n in sel]
    return dict(key=(mouse, sess), mouse=mouse, sess=sess, source="real",
                align_event=config.ALIGN_EVENT, spike_span=config.SPIKE_SPAN_S,
                reach_cat=reach_cat, shuffle_within=shuffle_within,
                trials=trials,
                kin=dict(t=kin_t, fps=config.FPS, var=kin_var, unit="mm", source=kin_source,
                         X=X[keep]),
                neurons=neurons, truth=None)


def list_sessions(regions=None, reach_cat=None, shuffle_within="default"):
    """DataFrame, one row per session: neurons per region, trials (all / kept), and
    trial counts per shuffle stratum (min_stratum = the smallest)."""
    regions   = list(regions or config.REGIONS)
    reach_cat = tuple(reach_cat or config.REACH_CAT)
    if shuffle_within == "default":
        shuffle_within = config.SHUFFLE_WITHIN
    NEURONS, REACH_PX = load_mi()
    rows = {}
    for n in NEURONS:
        if n["region"] not in regions: continue
        k = (n["mouse"], n["sess"])
        if k not in rows:
            tr = pd.DataFrame({LABEL_MAP[f]: np.asarray(n["labels"][f], float) for f in LABEL_MAP})
            keep = np.isin(tr["cat"], reach_cat) & np.isfinite(
                REACH_PX[k][:, config.KIN_COORD[config.KIN_VAR], :]).all(axis=1)
            cnt = pd.Series(cond_labels(tr[keep], shuffle_within)).value_counts()
            rows[k] = dict(mouse=k[0], sess=k[1], n_trials=len(tr), n_kept=int(keep.sum()),
                           **{r: 0 for r in regions},
                           **{f"n[{s}]": int(v) for s, v in sorted(cnt.items())},
                           min_stratum=int(cnt.min()) if len(cnt) else 0)
        rows[k][n["region"]] += 1
    df = pd.DataFrame(list(rows.values()))
    cnt_cols = [c for c in df.columns if c.startswith("n[")]
    df[cnt_cols] = df[cnt_cols].fillna(0).astype(int)
    return df


def neuron_table(sess, window_ms=None, pre_ms=(-1000.0, -375.0)):
    """Per-neuron firing statistics (rate-dependent decisions are made PER NEURON,
    never per region). One row per neuron:
      fr_win_hz   mean rate in the analysis window
      fr_pre_hz   mean rate in pre_ms (just before the window)
      frac_empty  fraction of trials with 0 spikes in the window
      n_spk_win   total spikes in the window, all trials (GLM feasibility)
      cv2         median CV2 of consecutive ISIs in the window (NaN if too few)"""
    w0, w1 = [v / 1000 for v in (window_ms or config.TRIAL_WINDOW_MS)]
    p0, p1 = [v / 1000 for v in pre_ms]
    rows = []
    for n in sess["neurons"]:
        cw = np.array([np.sum((s >= w0) & (s < w1)) for s in n["spikes"]])
        cp = np.array([np.sum((s >= p0) & (s < p1)) for s in n["spikes"]])
        cv2 = []
        for s in n["spikes"]:
            isi = np.diff(np.sort(s[(s >= w0) & (s < w1)]))
            if isi.size >= 2:
                cv2.extend(2 * np.abs(np.diff(isi)) / (isi[1:] + isi[:-1]))
        rows.append(dict(uid=n["uid"], region=n["region"],
                         fr_win_hz=cw.mean() / (w1 - w0),
                         fr_pre_hz=cp.mean() / (p1 - p0),
                         frac_empty=float(np.mean(cw == 0)),
                         n_spk_win=int(cw.sum()),
                         cv2=float(np.median(cv2)) if len(cv2) >= 20 else np.nan))
    return pd.DataFrame(rows)


def get_data(DATA_SOURCE="real", SIM_SCENARIO=None, session=None, sim_params=None, **load_kw):
    """The one entry point for every script's cell 0.
    DATA_SOURCE "real" -> load_session(*session)
                "sim"  -> simulate_session(template=load_session(*session), SIM_SCENARIO, **sim_params)
    Simulations always use the REAL kinematics of the template session."""
    session = session or config.FIRST_PASS[0]
    S = load_session(*session, **load_kw)
    if DATA_SOURCE == "real":
        return S
    if DATA_SOURCE == "sim":
        from . import sim
        return sim.simulate_session(S, SIM_SCENARIO, **(sim_params or {}))
    raise ValueError(f"DATA_SOURCE must be 'real' or 'sim', got {DATA_SOURCE!r}")


def describe(obj, name="sess", max_depth=3, max_items=12):
    """Print the structure of a session dict (or any nested object).
    Wraps mi_io.describe, and prints DataFrames as shape + columns + head."""
    if isinstance(obj, dict):
        print(f"{name}: dict({len(obj)} keys)")
        for k, v in obj.items():
            if isinstance(v, pd.DataFrame):
                print(f"  {k}: DataFrame {v.shape}  columns={list(v.columns)}")
                print("    " + v.head(3).to_string().replace("\n", "\n    "))
            else:
                mi_io.describe(v, max_depth=max_depth - 1, max_items=max_items,
                               _depth=1, _key=k)
        return
    mi_io.describe(obj, name=name, max_depth=max_depth, max_items=max_items)
