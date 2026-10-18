"""sim.py - simulator with KNOWN ground truth, on REAL kinematics of a template
session (so autocorrelation and task structure are realistic).

Discrete-time conditional-intensity model on a fine grid (dt_ms, default 0.5 ms):
    log lambda(t) = b0 + hist(t) + coupling * drive(t)   [drive: x_to_y, null_shared_drive]
    P(spike in [t, t+dt)) = 1 - exp(-lambda(t) dt)
  hist(t)  = own history, lags >= 1 step: absolute refractory (-20 for <= 1.5 ms)
             + relative recovery -A exp(-tau / tau_r), tau_r set by the target RATE
             (never by region): tau_r = clip(250 / rate_hz, 3, 25) ms.
  drive(t) = sum_{j>=1} k_x[j] x(t - j dt): causal gamma kernel with its MODE at
             lag_ms, applied to the session-z-scored kinematics (ZOH to the dt grid;
             the first tracked frame is held before tracking starts).
  b0 is calibrated per neuron so the empirical rate in the analysis window hits
  rate_hz (a few short calibration runs).

Scenarios:
  "null_independent"  history-only spikes; kinematics untouched     -> no GC either way
  "x_to_y"            every neuron driven by past kinematics         -> x->y only
  "y_to_x"            neuron 0 (the DRIVER; n_drivers > 1: neurons 0..n-1, contributions
                      summed) pushes the kinematics:                 -> y->x for the
                      X_sim(t_f) = X(t_f) + coupling * sd(X) * z(k * s)(t_f), where s   driver only;
                      = the driver's spike counts in [t_f - 2.5 ms, t_f) and k a causal  the others
                      gamma kernel with its mode at lag_ms. Only spikes strictly before   are null
                      a frame reach it. Other neurons are history-only (null both ways).
  "null_shared_drive" every neuron driven by the SPOUT-SIDE MEAN of   -> NO trial-specific
                      the filtered kinematics (identical for all trials   coupling either way:
                      of a side; same kernel and coupling as x_to_y):     the null for
                      task-locked rate that tracks the AVERAGE reach,     task-locked activity
                      not the trial's own. Kinematics untouched. A test
                      that ignores the task structure (chi2) may fire; the
                      within-side trial shuffle must not.
Others (bidirectional, common_latent): Milestone 4.
Output: a session dict with the same layout as io.load_session, plus `truth`.
"""
import numpy as np
import pandas as pd
from . import config
from .preprocess import zoh_resample, causal_filter, gamma_kernel, remove_ensemble_mean, bin_spikes
from .design import lag_matrix
from .nulls import perm_within          # re-exported (it lived here in M1)

SCENARIOS = ("null_independent", "x_to_y", "y_to_x", "null_shared_drive")


# ---------------------------------------------------------------------------
# templates
# ---------------------------------------------------------------------------
def synthetic_template(n_trials=300, kin="white", seed=0, n_cond=3, pre_ms=0.0):
    """Session dict with SYNTHETIC kinematics and no neurons (for unit tests that must
    not depend on the data drive). kin "white" = iid N(0,1) per frame, which makes the
    spike-triggered average recover a coupling kernel directly.
    pre_ms > 0: kinematics also exist pre_ms before -375 ms (like io kin_source="long")."""
    rng = np.random.default_rng(seed)
    lcr = rng.integers(1, n_cond + 1, n_trials).astype(float)
    trials = pd.DataFrame(dict(trial_idx=np.arange(n_trials), LCR=lcr,
                               PP=rng.integers(0, 2, n_trials).astype(float),
                               cat=np.ones(n_trials), hit=np.nan, succ=np.nan))
    trials["cond"] = [f"LCR={int(v)}" for v in lcr]
    n_pre = int(round(pre_ms * config.FPS / 1000.0))
    n_fr = config.N_FRAMES + n_pre
    if kin == "white":
        X = rng.normal(size=(n_trials, n_fr))
    elif kin == "ar1":
        X = causal_filter(rng.normal(size=(n_trials, n_fr)), 0.1 * 0.9 ** np.arange(60))
    else:
        raise ValueError(kin)
    t = (np.arange(n_fr) - config.ALIGN_FRAME - n_pre) / config.FPS
    return dict(key=("synthetic", kin), mouse="synthetic", sess=kin, source="synthetic",
                align_event=config.ALIGN_EVENT, spike_span=config.SPIKE_SPAN_S,
                reach_cat=(1,), shuffle_within=("LCR",), trials=trials,
                kin=dict(t=t, fps=config.FPS, var="synthetic", unit="a.u.", X=X),
                neurons=[], truth=None)


# ---------------------------------------------------------------------------
# kernels
# ---------------------------------------------------------------------------
def coupling_kernel(lag_ms, shape, dt_ms):
    """Causal gamma kernel with its MODE at lag_ms, strictly lagged (k[0] = 0: the
    mass of [j dt, (j+1) dt) is assigned to lag j+1). -> (k, lags_ms)."""
    if shape <= 1:
        raise ValueError("shape must be > 1 for the mode to be at lag_ms")
    k, _ = gamma_kernel(shape, lag_ms / (shape - 1), dt_ms)
    k = np.concatenate([[0.0], k])
    return k, np.arange(k.size) * dt_ms


def _tau_r(rate_hz):
    return float(np.clip(250.0 / rate_hz, 3.0, 25.0))


def history_kernel(rate_hz, dt_ms, abs_ref_ms=1.5, A=3.0, abs_val=-20.0):
    """The own-history filter the simulator uses, as an explicit array (for `truth`)."""
    tau_r = _tau_r(rate_hz)
    tau = np.arange(1, int(np.ceil((abs_ref_ms + 6 * tau_r) / dt_ms)) + 1) * dt_ms
    h = -A * np.exp(-tau / tau_r)
    h[tau <= abs_ref_ms + 1e-9] = abs_val
    return tau, h


# ---------------------------------------------------------------------------
# core generator
# ---------------------------------------------------------------------------
def _run(b0, drive, coupling, tau_r, dt_ms, rng, abs_ref_ms=1.5, A=3.0, abs_val=-20.0):
    """b0 (n_neu,), drive (n_tr, n_steps) or None -> spikes bool (n_neu, n_tr, n_steps).
    History is updated recursively and only from PAST steps (j >= 1)."""
    n_neu, (n_tr, n_st) = len(b0), drive.shape
    dt_s  = dt_ms / 1000.0
    decay = np.exp(-dt_ms / np.asarray(tau_r))[:, None]            # (n_neu, 1)
    n_abs = int(round(abs_ref_ms / dt_ms))
    g     = np.zeros((n_neu, n_tr))                                 # sum_j s[t-j] e^{-j dt/tau}
    since = np.full((n_neu, n_tr), 10 ** 6)                         # steps since last spike
    S     = np.zeros((n_neu, n_tr, n_st), dtype=bool)
    b0    = np.asarray(b0, float)[:, None]
    for s in range(n_st):
        eta = b0 - A * g + np.where(since <= n_abs, abs_val, 0.0)
        if coupling:
            eta = eta + coupling * drive[None, :, s]
        spk = rng.random((n_neu, n_tr)) < -np.expm1(-np.exp(eta) * dt_s)
        S[:, :, s] = spk
        g = decay * (g + spk)               # becomes history for step s+1 onward
        since = np.where(spk, 1, since + 1)
    return S


def spike_drive(spikes, t_frames, k, t_start):
    """Causal spike -> kinematics drive at each frame time t_f (y_to_x scenario).
    Spike counts in [t_f - FRAME, t_f) are filtered with k (k[0] = 0 from
    coupling_kernel); history counted from t_start. -> (n_trials, n_frames)."""
    dt = config.FRAME_MS / 1000.0
    n_pre = int(round((t_frames[0] - t_start) / dt))
    edges = t_frames[0] + np.arange(-n_pre - 1, len(t_frames)) * dt
    C = bin_spikes(spikes, edges).astype(float)
    return causal_filter(C, k, pre="zero")[:, n_pre:]


def simulate_session(template, scenario, rate_hz=(2.0, 6.0, 20.0, 60.0), n_per_rate=5,
                     n_neurons=None, driver_rate_hz=20.0, n_trials=None, coupling=1.0,
                     lag_ms=50.0, kernel_shape=5.0, seed=0, dt_ms=0.5, window_ms=None,
                     hist_pad_ms=300.0, calib_iter=3, n_drivers=1):
    """Simulated session on the template's REAL kinematics (see module docstring).
      rate_hz         tuple of target rates (Hz), each used for n_per_rate neurons, or
                      "match" = the template's real per-neuron window rates (io.neuron_table)
      n_neurons       keep a random subset of that many neurons (None = all)
      n_drivers       y_to_x only: how many neurons LEAD the kinematics in this session
                      (neurons 0..n_drivers-1, "leaders"/"drivers"); the rest are null.
                      Their contributions ADD UP (each = coupling, as if alone), so the
                      total added signal grows ~sqrt(n_drivers): see truth["added_sd_frac"].
      driver_rate_hz  y_to_x only: rate of every leader (None = keep its own rate)
      coupling        x_to_y: gain on the z-scored kinematic drive (log-rate units)
                      y_to_x: size of each leader's added term, in units of sd(kinematics)
      n_trials        None = all template trials, else a random subset of that size
    -> session dict (io layout) with region "SIM" and a `truth` dict."""
    if scenario not in SCENARIOS:
        raise NotImplementedError(f"scenario {scenario!r}: this milestone has {SCENARIOS}")
    window_ms = window_ms or config.TRIAL_WINDOW_MS
    rng = np.random.default_rng(seed)

    trials, X = template["trials"], np.asarray(template["kin"]["X"], float)
    if n_trials is not None and n_trials < len(trials):
        pick = np.sort(rng.choice(len(trials), n_trials, replace=False))
        trials, X = trials.iloc[pick].reset_index(drop=True), X[pick]

    if isinstance(rate_hz, str) and rate_hz == "match":
        from .io import neuron_table
        rates = np.clip(neuron_table(template, window_ms)["fr_win_hz"].to_numpy(), 0.5, None)
    else:
        rates = np.repeat(np.asarray(rate_hz, float), n_per_rate)
    if n_neurons is not None and n_neurons < rates.size:
        rates = rates[np.sort(rng.choice(rates.size, n_neurons, replace=False))]
    rates = rates.copy()
    n_drv = int(min(max(n_drivers, 1), rates.size)) if scenario == "y_to_x" else 0
    if scenario == "y_to_x" and driver_rate_hz is not None:
        rates[:n_drv] = driver_rate_hz
    tau_r = np.array([_tau_r(r) for r in rates])

    # fine time grid: [window start - pad, window end)
    t0 = (window_ms[0] - hist_pad_ms) / 1000.0
    n_st = int(round((window_ms[1] - window_ms[0] + hist_pad_ms) / dt_ms))
    t_st = t0 + np.arange(n_st) * dt_ms / 1000.0
    in_win = t_st >= window_ms[0] / 1000.0 - 1e-12

    # kinematic drive (z-scored per session, ZOH, first frame held before tracking)
    # normalised on the analysis-window frames (same with a long REACH_LONG template)
    kt = np.asarray(template["kin"]["t"], float)
    wfr = (kt >= window_ms[0] / 1000.0 - 1e-9) & (kt <= window_ms[1] / 1000.0 + 1e-9)
    mu, sd = float(np.mean(X[:, wfr])), float(np.std(X[:, wfr])) or 1.0
    xz = (X - mu) / sd
    xs = zoh_resample(template["kin"]["t"], xz, t_st, max_hold_s=config.FRAME_MS / 1000.0)
    lead = t_st < template["kin"]["t"][0]
    xs[:, lead] = xz[:, :1]
    if np.isnan(xs).any():
        raise ValueError("simulation grid extends beyond the tracked kinematics")
    k_x, k_t = coupling_kernel(lag_ms, kernel_shape, dt_ms)
    drive = causal_filter(xs, k_x, pre="hold")
    if scenario == "null_shared_drive":                 # same drive for every trial of a side
        cnd = trials["cond"].to_numpy()
        for cv in np.unique(cnd):
            drive[cnd == cv] = drive[cnd == cv].mean(axis=0)
    c = coupling if scenario in ("x_to_y", "null_shared_drive") else 0.0

    # calibrate b0 so the in-window rate matches the target
    b0 = np.log(rates)
    win_s = in_win.sum() * dt_ms / 1000.0
    for it in range(calib_iter):
        S = _run(b0, drive, c, tau_r, dt_ms, np.random.default_rng(seed + 1000 + it))
        emp = S[:, :, in_win].sum(axis=(1, 2)) / (S.shape[1] * win_s)
        b0 = b0 + np.log(rates / np.maximum(emp, 1e-3))
    S = _run(b0, drive, c, tau_r, dt_ms, rng)

    # spike times: uniform within the step in which the spike occurred
    neurons = []
    for i in range(len(rates)):
        spikes = []
        for tr in range(S.shape[1]):
            st = np.flatnonzero(S[i, tr])
            spikes.append(t_st[st] + rng.random(st.size) * dt_ms / 1000.0)
        neurons.append(dict(uid=f"SIM|{template['mouse']}|{template['sess']}|s{seed}|{i}",
                            region="SIM", phyID=i, spikes=spikes,
                            rate_target_hz=float(rates[i]), b0=float(b0[i])))

    # y_to_x: the leaders' past spikes are added to the kinematics. Each leader adds
    # coupling * sd(X) * z(its filtered spikes) -- the SAME amount it would alone -- so the
    # total added signal grows ~sqrt(n_drivers) (leaders' spikes are independent).
    expect_xy = np.full(len(rates), scenario == "x_to_y")
    expect_yx = np.zeros(len(rates), bool)
    X_out, k_y, k_yt, drivers, added = X, None, None, [], 0.0
    if scenario == "y_to_x":
        k_y, k_yt = coupling_kernel(lag_ms, kernel_shape, config.FRAME_MS)
        drivers = list(range(n_drv))
        add = np.zeros_like(X)
        # spikes exist only on the simulated span: frames outside it (long template) get no
        # added signal -- the analysis never uses them (its history pad is shorter)
        span = (kt >= t_st[0] - 1e-9) & (kt <= t_st[-1] + dt_ms / 1000.0 + 1e-9)
        w_in = wfr[span]
        for j in drivers:
            D = spike_drive(neurons[j]["spikes"], kt[span], k_y, t_st[0])
            add[:, span] += coupling * sd * (D - D[:, w_in].mean()) / (D[:, w_in].std() or 1.0)
        X_out = X + add
        expect_yx[drivers] = True
        added = float(add[:, wfr].std() / sd)

    truth = dict(scenario=scenario, rate_hz=rates, coupling=float(coupling), lag_ms=float(lag_ms),
                 kernel_shape=float(kernel_shape), kernel_t_ms=k_t, kernel=k_x,
                 kernel_y_t_ms=k_yt, kernel_y=k_y,
                 driver=drivers[0] if drivers else None, drivers=drivers,
                 added_sd_frac=added,                      # size of added signal / sd(med-lat)
                 shared_drive=(scenario == "null_shared_drive"),
                 hist_tau_r_ms=tau_r, kin_norm=(mu, sd), dt_ms=dt_ms, seed=seed,
                 expect={"x_to_y": expect_xy, "y_to_x": expect_yx})
    return dict(key=template["key"], mouse=template["mouse"], sess=template["sess"],
                source="sim", align_event=template["align_event"],
                spike_span=(float(t_st[0]), float(t_st[-1] + dt_ms / 1000.0)),
                reach_cat=template.get("reach_cat"),
                shuffle_within=template.get("shuffle_within"),
                trials=trials, kin=dict(template["kin"], X=X_out), neurons=neurons, truth=truth)


# ---------------------------------------------------------------------------
# diagnostics (sanity checks only; the tests of precedence are the models)
# ---------------------------------------------------------------------------
def sta_kin(Y, X, in_window, lags, pair=None, resid_cond=None):
    """Spike-triggered average of PAST kinematics, per lag.
      Y (n_trials, n_bins) counts of one neuron; X (n_trials, n_bins) binned kinematics
      pair       optional (n_trials,) partner index: neuron trial i vs kin trial pair[i]
      resid_cond optional condition labels: subtract the leave-one-out condition mean
                 trajectory from X first (removes the task-locked mean)
    -> (n_lags,) sum(y[t] x[t-j]) / sum(y[t]) - mean(x[t-j]) over the same rows."""
    X = np.asarray(X, float)
    if resid_cond is not None:
        X = remove_ensemble_mean(X, resid_cond)
    if pair is not None:
        X = X[pair]
    L = np.stack([lag_matrix(X[i], lags) for i in range(X.shape[0])])   # (n_tr, n_b, n_lags)
    ok = np.asarray(in_window, bool)[None, :, None] & np.isfinite(L)
    Yb = np.broadcast_to(np.asarray(Y, float)[:, :, None], L.shape)
    Lz = np.where(ok, L, 0.0)
    num = (Yb * Lz).sum(axis=(0, 1))
    den = np.where(ok, Yb, 0.0).sum(axis=(0, 1))
    base = Lz.sum(axis=(0, 1)) / np.maximum(ok.sum(axis=(0, 1)), 1)
    return num / np.maximum(den, 1e-12) - base
