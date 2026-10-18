"""runs.py - the run machinery shared by scripts 01, 02, ... (everything that is not
the model itself): which sessions, readable run names + run_info.txt, loading real or
simulated sessions with per-neuron inclusion, BH-FDR + direction classes, the CP vs CB
printout, the simulation truth check, saving tables, and the direction figures.

A "setting" is one combination of the values in `setting_cols` of the result table
RES (01: kin_smooth x order; 02: none). FDR is applied per setting x region x direction,
pooling neurons over sessions -- that is the CP vs CB comparison.
Result tables carry, per direction d in {"xy" (kin->neuron), "yx" (neuron->kin)}:
  pS_d (shuffle p), z_d, sig_d, q_d, and optional lag columns.
"""
import time
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu
from . import config, io, sim, stats, plotting
from .cache import save_pickle, _long

SCEN_NAME = {"null_independent": "noCoupling", "x_to_y": "kin2neuron", "y_to_x": "neuron2kin",
             "null_shared_drive": "sharedDrive"}


# ---------------------------------------------------------------------------
# which sessions, names, run info
# ---------------------------------------------------------------------------
def stag(mouse, sess):
    """Session tag for names: ("5_FerreroRocher", "R1") -> "FerreroRocher-R1"."""
    return f"{mouse.split('_', 1)[-1]}-{sess}"


def session_list(data_source, sessions, min_stratum, reach_cat, template):
    """Real: "first" | "all" (sessions whose smallest LCR stratum >= min_stratum) |
    list of (mouse, sess). Simulation: [template]."""
    if data_source != "real":
        return [tuple(template)]
    if sessions == "first":
        return [tuple(config.FIRST_PASS[0])]
    if sessions == "all":
        ls = io.list_sessions(reach_cat=reach_cat)
        ok = ls.min_stratum >= min_stratum
        out = [tuple(r) for r in ls.loc[ok, ["mouse", "sess"]].itertuples(index=False)]
        print(f"sessions passing MIN_STRATUM={min_stratum}: {len(out)}/{len(ls)}; skipped: "
              f"{[stag(*r) for r in ls.loc[~ok, ['mouse', 'sess']].itertuples(index=False)]}")
        return out
    return [tuple(s) for s in sessions]


def data_name(data_source, sess_list, sim_scenario, sim_params, sim_seeds):
    """First part of RUN_NAME: real_<session | N sessions> or
    sim-<scenario>_true<lag>ms_c<coupling>_<N>seed(s)."""
    if data_source == "real":
        return (f"real_{stag(*sess_list[0])}" if len(sess_list) == 1
                else f"real_{len(sess_list)}sessions")
    n = len(sim_seeds)
    k = int(sim_params.get("n_drivers", 1) or 1)
    return (f"sim-{SCEN_NAME.get(sim_scenario, sim_scenario)}_true{sim_params['lag_ms']:g}ms"
            f"_c{sim_params['coupling']:g}_{n}seed{'s' if n > 1 else ''}"
            + (f"_k{k}leaders" if sim_scenario == "y_to_x" and k > 1 else ""))


def common_info(data_source, sess_list, reach_cat, min_stratum, min_rate_hz, sim_scenario,
                sim_params, sim_seeds, template, n_shuffle, fdr_q, kin_var=None, kin_source=None):
    """run_info.txt lines shared by every script."""
    kin_var = kin_var or config.KIN_VAR
    kin_source = kin_source or config.KIN_SOURCE
    return [
        f"data: {data_source}" + (f", scenario {sim_scenario}, true lag {sim_params['lag_ms']:g} ms, "
                                  f"coupling {sim_params['coupling']:g}, seeds {list(sim_seeds)}, "
                                  f"template {template}"
                                  + (f", leaders per session (n_drivers) {sim_params.get('n_drivers', 1)}"
                                     if sim_scenario == "y_to_x" else "")
                                  if data_source == "sim" else ""),
        f"sessions ({len(sess_list)}): {', '.join(stag(*m) for m in sess_list)}",
        f"reach category: {list(reach_cat)};  min trials per spout side: {min_stratum};  "
        f"min neuron rate: {min_rate_hz:g} Hz",
        f"analysis window: {config.TRIAL_WINDOW_MS} ms re max reach;  lag range "
        f"{config.CROSS_LAG_RANGE_MS} ms;  kinematics: {kin_var} ({config.KIN_NAME[kin_var]}, mm), "
        f"source {kin_source!r}" + (" (REACH_LONG: real history before -375 ms)" if kin_source == "long" else ""),
        f"significance: trial shuffle within {config.SHUFFLE_WITHIN}, {n_shuffle} shuffles; "
        f"BH-FDR q < {fdr_q} per region x direction, pooled over sessions",
    ]


def write_run_info(text, fig_dir, lag_tag, out_subdir, run_name, save_figs, save_vars):
    if save_figs and fig_dir is not None:
        with open(_long(fig_dir / "run_info.txt"), "w", encoding="utf-8") as f:
            f.write(text + "\n")
    if save_vars:
        with open(_long(config.var_path(lag_tag, out_subdir, f"{run_name}_run_info.txt")),
                  "w", encoding="utf-8") as f:
            f.write(text + "\n")
    print(text)
    print(f"figures  {fig_dir}\ntables   {config.VAR_DIR / lag_tag / out_subdir}")


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
def load_data(data_source, sess_list, reach_cat, template, sim_scenario, sim_params,
              sim_seeds, min_rate_hz, kin_var=None, kin_source=None):
    """-> (SESSIONS_DATA [(tag, session dict)], NT_ALL per-neuron table).
    Sets S["_keep"] (rate >= min_rate_hz, per NEURON) and S["_rate"] on every session.
    kin_var: kinematic variable (None = config.KIN_VAR); simulations use it from the template.
    kin_source: "short" | "long" (None = config.KIN_SOURCE); simulations use the template's."""
    sessions, t0 = [], time.time()
    if data_source == "real":
        for ms in sess_list:
            sessions.append((stag(*ms), io.load_session(*ms, kin_var=kin_var, reach_cat=reach_cat, verbose=False,
                                                                kin_source=kin_source)))
    else:
        T = io.load_session(*template, kin_var=kin_var, reach_cat=reach_cat, kin_source=kin_source)
        for seed in sim_seeds:
            sessions.append((f"s{seed}", sim.simulate_session(T, sim_scenario, seed=seed, **sim_params)))
    print(f"{len(sessions)} session(s) ready in {time.time()-t0:.1f}s")
    nts = []
    for tag, S in sessions:
        nt = io.neuron_table(S)
        keep = nt["fr_win_hz"] >= min_rate_hz
        S["_keep"] = keep.to_numpy(); S["_rate"] = nt["fr_win_hz"].to_numpy()
        nts.append(nt.assign(tag=tag, keep=keep))
        tr = S.get("truth") or {}
        lead = (f", leaders {tr['drivers']} (added signal = {tr['added_sd_frac']:.2f} x sd(med-lat))"
                if tr.get("drivers") else "")
        print(f"  {tag:18s}: {len(S['trials']):4d} trials, neurons kept {keep.sum():3d}/{len(nt):3d} "
              f"{nt[keep].groupby('region').size().to_dict()}{lead}")
    nt_all = pd.concat(nts, ignore_index=True)
    print(f"TOTAL kept neurons: {nt_all[nt_all.keep].groupby('region').size().to_dict()}")
    return sessions, nt_all


def base_row(tag, S, i):
    """Identity columns of one neuron's result row (+ simulation truth if any)."""
    n = S["neurons"][i]
    row = dict(mouse=S["mouse"], sess=S["sess"], tag=tag, region=n["region"], uid=n["uid"],
               rate_hz=S["_rate"][i], neuron_idx=int(i))
    if S.get("truth"):
        e = S["truth"]["expect"]
        drv = S["truth"].get("drivers")
        if drv is None:
            drv = [] if S["truth"].get("driver") is None else [S["truth"]["driver"]]
        row.update(exp_xy=bool(e["x_to_y"][i]), exp_yx=bool(e["y_to_x"][i]), driver=(i in drv))
    return row


# ---------------------------------------------------------------------------
# settings, FDR, summaries
# ---------------------------------------------------------------------------
def settings_of(RES, setting_cols):
    """Unique setting tuples in table order ([()] if no setting columns)."""
    if not setting_cols:
        return [()]
    return list(dict.fromkeys(map(tuple, RES[setting_cols].itertuples(index=False))))


def select(RES, setting_cols, values):
    m = np.ones(len(RES), bool)
    for c, v in zip(setting_cols, values):
        m &= (RES[c] == v).to_numpy()
    return RES[m]


def fdr_and_classes(RES, setting_cols, q):
    """BH-FDR per setting x region x direction (neurons pooled over sessions) ->
    sig_d, q_d, and the direction class."""
    for _, g in RES.groupby(list(setting_cols) + ["region"]):
        for d in ("xy", "yx"):
            sig, qq = stats.bh_fdr(g[f"pS_{d}"].to_numpy(), q)
            RES.loc[g.index, f"sig_{d}"] = sig; RES.loc[g.index, f"q_{d}"] = qq
    RES[["sig_xy", "sig_yx"]] = RES[["sig_xy", "sig_yx"]].astype(bool)
    RES["direction"] = plotting.direction_class(RES.sig_xy, RES.sig_yx)
    return RES


def print_direction_summary(RES, setting_cols, q, n_sessions, label=lambda v: ""):
    print(f"\n=== DIRECTION CLASSES (BH-FDR q < {q} per direction, region pooled over "
          f"{n_sessions} session(s)) ===")
    for v in settings_of(RES, setting_cols):
        R = select(RES, setting_cols, v)
        tab = pd.crosstab(R.region, R.direction).reindex(columns=list(plotting.DIR_CLASSES), fill_value=0)
        pct = (100 * tab.div(tab.sum(1), axis=0)).round(1)
        print(f"\n[{label(v)}]  counts / %:" if label(v) else "\ncounts / %:")
        print(pd.concat([tab.assign(n=tab.sum(1)), pct.add_suffix(" %")], axis=1).to_string())
        if {"CP", "CB"} <= set(R.region):
            for d in ("yx", "xy"):
                a = [int(R[R.region == r][f"sig_{d}"].sum()) for r in ("CP", "CB")]
                n = [int((R.region == r).sum()) for r in ("CP", "CB")]
                _, pf = fisher_exact([[a[0], n[0] - a[0]], [a[1], n[1] - a[1]]])
                print(f"  {plotting.DIR_LABEL[d]} significant: CP {a[0]}/{n[0]} ({100*a[0]/n[0]:.0f}%) "
                      f"vs CB {a[1]}/{n[1]} ({100*a[1]/n[1]:.0f}%)  Fisher p = {pf:.2g}")


def _extend_task(null_obj, obs, perms, n0, cap, h, batch):
    """One neuron x direction: sequential p over the session's shuffle sequence."""
    return stats.sequential_p(obs, 0, n0, lambda a, b: null_obj.stats(perms[a:b]), cap, h, batch)


def adaptive_extend(RES, setting_cols, obs_cols, n0, q, cap_mode, h, batch, seed, cond_of,
                    null_of, cache_params, n_jobs, cache_name, cache_mode="load", save_cache=True,
                    label=lambda v: ""):
    """ADAPTIVE SHUFFLES. Neurons whose observed statistic beat all n0 shuffles (p at the
    floor 1/(n0+1)) AND are not already BH-significant get more shuffles -- the continuation of the SAME session sequence --
    in batches, until h shuffles reach them (p = h/L) or the cap (p = (g+1)/(cap+1));
    stats.sequential_p (Besag & Clifford). Cap per setting x region x direction =
    ceil(m / q) with m neurons (region pooled over sessions): beyond it, more shuffles
    cannot change the BH decision. Cached per session x setting.
      obs_cols      {"xy": column of the observed shuffle statistic, "yx": ...}
      cond_of(tag, v)        stratum labels of that session
      null_of(row, d)        models.VarNull / GlmNull for that neuron and direction
      cache_params(tag, v)   dict identifying the underlying fits (cache key)
    Updates pS_d and adds n_shuf_d, p_method_d ("fixed" | "sequential" | "cap")."""
    import time
    from joblib import Parallel, delayed
    from . import nulls
    from .cache import cached
    RES = RES.copy()
    floor = 1.0 / (n0 + 1)
    for d in ("xy", "yx"):
        RES[f"n_shuf_{d}"] = n0
        RES[f"p_method_{d}"] = "fixed"
    tasks = {}                                        # (tag, v) -> list of (index, d, obs, cap)
    print(f"\n=== ADAPTIVE SHUFFLES (h = {h}, batch {batch}; cap = ceil(m/q) per region x direction) ===")
    for v in settings_of(RES, setting_cols):
        R = select(RES, setting_cols, v)
        for reg, g in R.groupby("region"):
            for d in ("xy", "yx"):
                cap = int(np.ceil(len(g) / q)) if cap_mode == "auto" else int(cap_mode)
                # a floor neuron that is ALREADY BH-significant is settled: more shuffles only
                # lower its p, which never undoes its significance and can only relax BH for
                # the others. So only floor neurons that are NOT yet significant are extended.
                sig_now = stats.bh_fdr(g[f"pS_{d}"].to_numpy(), q)[0]
                at_all = np.isclose(g[f"pS_{d}"], floor)
                at = g[at_all & ~sig_now]
                if label(v):
                    print(f"  [{label(v)}] ", end="")
                print(f"  {reg} {plotting.DIR_LABEL[d]}: {int(at_all.sum())} at the floor of {len(g)} "
                      f"({int((at_all & sig_now).sum())} already significant = settled; "
                      f"{len(at)} to extend), cap {cap}")
                if cap <= n0:
                    continue
                for ix, row in at.iterrows():
                    tasks.setdefault((row.tag, v), []).append((ix, d, float(row[obs_cols[d]]), cap))
    n_tasks = sum(len(t) for t in tasks.values())
    if not n_tasks:
        print("  nothing to extend.")
        return RES
    worst = sum(cap - n0 for t in tasks.values() for (_, _, _, cap) in t)
    (tag0, v0), t0list = next(iter(tasks.items()))
    ix0, d0, _, _ = t0list[0]
    probe = null_of(RES.loc[ix0], d0)
    pp = nulls.make_perms(cond_of(tag0, v0), 60, seed)
    probe.stats(pp[:10])                               # setup
    tt = time.time(); probe.stats(pp[10:60]); per = (time.time() - tt) / 50
    print(f"  {n_tasks} neuron x direction(s) to extend; worst case {worst:,} extra shuffles "
          f"~ {worst * per / max(n_jobs if n_jobs > 0 else 8, 1) / 60:.0f} min on {n_jobs} jobs "
          f"({per*1000:.1f} ms per shuffle; most stop far earlier than the cap)", flush=True)
    t0 = time.time()
    for (tag, v), tl in tasks.items():
        max_cap = max(c for *_, c in tl)
        def _run(tl=tl, tag=tag, v=v, max_cap=max_cap):
            perms = nulls.make_perms(cond_of(tag, v), max_cap, seed)
            res = Parallel(n_jobs=n_jobs, backend="loky")(
                delayed(_extend_task)(null_of(RES.loc[ix], d), obs, perms, n0, cap, h, batch)
                for ix, d, obs, cap in tl)
            return {(RES.loc[ix, "uid"], d): r for (ix, d, _, _), r in zip(tl, res)}
        params = dict(fit=cache_params(tag, v), n0=n0, h=h, batch=batch, seed=seed,
                      tasks=[(RES.loc[ix, "uid"], d, obs, cap) for ix, d, obs, cap in tl])
        vtag = "_".join(map(str, v)) or "all"
        out = cached(f"{cache_name}_{tag}_{vtag}_adapt", params, _run, mode=cache_mode, save=save_cache)
        for ix, d, _, _ in tl:
            p, n_used, method = out[(RES.loc[ix, "uid"], d)]
            RES.loc[ix, [f"pS_{d}", f"n_shuf_{d}", f"p_method_{d}"]] = [p, n_used, method]
    for d in ("xy", "yx"):
        vc = RES.loc[RES[f"n_shuf_{d}"] > n0, f"p_method_{d}"].value_counts().to_dict()
        if vc:
            print(f"  {plotting.DIR_LABEL[d]}: {vc}  (sequential = resolved; cap = beat every shuffle up to the cap)")
    print(f"  adaptive shuffles done in {time.time()-t0:.0f}s")
    return RES


def print_floor_summary(RES, setting_cols, n_shuffle, q, label=lambda v: ""):
    """How much the shuffle p-value FLOOR (1/(N+1)) limits FDR, per setting x region x
    direction. A neuron 'at the floor' beat every shuffle: its true p could be far smaller.
    BH can pass floor neurons only if at least ceil(floor * n / q) neurons share a p that
    small; 'floor, not sig' counts neurons that beat every shuffle yet failed FDR -- the
    ones more shuffles (or a tail fit) would help."""
    floor = 1.0 / (n_shuffle + 1)
    print(f"\n=== SHUFFLE p-VALUE FLOOR ({n_shuffle} shuffles -> smallest p = 1/{n_shuffle + 1} "
          f"= {floor:.2g}) ===")
    adaptive = "n_shuf_xy" in RES
    rows = []
    for v in settings_of(RES, setting_cols):
        R = select(RES, setting_cols, v)
        for reg, g in R.groupby("region"):
            for d in ("xy", "yx"):
                n = len(g)
                row = {**({"setting": label(v)} if label(v) else {}), "region": reg,
                       "direction": plotting.DIR_LABEL[d], "n": n}
                if adaptive:          # floor = each neuron's OWN smallest p, 1/(n_used+1)
                    ext = g[f"n_shuf_{d}"] > n_shuffle
                    at = np.isclose(g[f"pS_{d}"], 1.0 / (g[f"n_shuf_{d}"] + 1))
                    row.update({f"extended (floor, not sig after {n_shuffle})": int(ext.sum()),
                                "resolved (sequential)": int((g[f"p_method_{d}"] == "sequential").sum()),
                                "beat all up to cap": int((g[f"p_method_{d}"] == "cap").sum() - 0),
                                "FDR sig": int(g[f"sig_{d}"].sum()),
                                "floor, not sig": int((at & ~g[f"sig_{d}"]).sum())})
                else:
                    at = np.isclose(g[f"pS_{d}"], floor)
                    row.update({"at floor": int(at.sum()),
                                "floor neurons BH needs": int(np.ceil(floor * n / q)),
                                "FDR sig": int(g[f"sig_{d}"].sum()),
                                "floor, not sig": int((at & ~g[f"sig_{d}"]).sum())})
                row["lone-neuron BH threshold"] = f"{q / n:.1e}"
                rows.append(row)
    T = pd.DataFrame(rows)
    print(T.to_string(index=False))
    if (T["floor, not sig"] > 0).any():
        print("  -> some neurons beat EVERY shuffle but failed FDR: more shuffles would resolve them.")
    else:
        print("  -> no neuron was lost to the floor: the number of shuffles is sufficient here.")
    return T


def _ci(lo, hi):
    """pi1 CI text: the 95% CI resamples whole SESSIONS (neurons of a session share the
    kinematics and the shuffles); needs >= 3 sessions."""
    return f"[{lo:.2f}, {hi:.2f}]" if np.isfinite(lo) else "[CI n/a: < 3 sessions]"


def print_effect_summary(RES, setting_cols, regions, effect_cols, label=lambda v: ""):
    """HOW STRONG and HOW MANY, per setting x direction x region (never across directions):
    pi1 = estimated fraction of neurons with an effect (Storey; includes effects too weak
    for FDR), and the median EXCESS effect (effect minus its shuffle mean) for FDR-significant
    neurons and for all neurons, with CP vs CB Mann-Whitney.
    effect_cols = {"xy": (column, unit label), "yx": (column, unit label)}."""
    print("\n=== EFFECT STRENGTH and FRACTION OF NEURONS WITH AN EFFECT (compare regions within a "
          "direction only) ===")
    for v in settings_of(RES, setting_cols):
        R = select(RES, setting_cols, v)
        if label(v):
            print(f"[{label(v)}]")
        for d in ("xy", "yx"):
            col, unit = effect_cols[d]
            print(f"  {plotting.DIR_LABEL[d]}   effect = {unit}")
            vals = {}
            for reg in regions:
                g = R[R.region == reg]
                if not len(g):
                    continue
                pi1, lo, hi = stats.storey_pi1(g[f"pS_{d}"], groups=g["tag"])   # CI: resample SESSIONS
                s, a = g[g[f"sig_{d}"]][col].dropna(), g[col].dropna()
                vals[reg] = (s, a)
                print(f"    {reg:4s} n={len(g):5d}  FDR sig {int(g[f'sig_{d}'].sum()):5d}  "
                      f"pi1 = {pi1:5.2f} {_ci(lo, hi)}  |  median effect: significant "
                      f"{s.median() if len(s) else np.nan:8.4f} (n={len(s)}), all {a.median():8.4f}")
            if {"CP", "CB"} <= set(vals):
                for k, nm in ((0, "significant"), (1, "all")):
                    x, y = vals["CP"][k], vals["CB"][k]
                    if len(x) >= 3 and len(y) >= 3:
                        print(f"    CP vs CB ({nm}): Mann-Whitney p = {mannwhitneyu(x, y).pvalue:.2g}")


def plot_effect_figures(RES, setting_cols, regions, effect_cols, fig_dir, fig_name, save, show,
                        title, label=lambda v: "", rate_flag_hz=None):
    """8_effect_strength: per direction (rows), A = p-value histograms with pi1, B = ECDF of
    the excess effect (solid = FDR-significant, dashed = all), C = effect vs -log10 p with
    each region's FDR cutoff. Simulations: truly coupled neurons circled in C.
    rate_flag_hz: mark neurons above this rate with 'x' in C (Poisson effect mis-scaled)."""
    import matplotlib.pyplot as plt
    for v in settings_of(RES, setting_cols):
        R = select(RES, setting_cols, v)
        fig, axs = plt.subplots(2, 3, figsize=(16, 8.5))
        for r_, d in enumerate(("xy", "yx")):
            col, unit = effect_cols[d]
            ax = axs[r_, 0]                                       # A: p-value histogram
            edges = np.linspace(0, 1, 21)
            for reg in regions:
                g = R[R.region == reg]
                if not len(g):
                    continue
                pi1, lo, hi = stats.storey_pi1(g[f"pS_{d}"], groups=g["tag"])   # CI: resample SESSIONS
                ax.hist(g[f"pS_{d}"].dropna(), bins=edges, density=True, histtype="step", lw=1.8,
                        color=config.CLR.get(reg, (0.4, 0.4, 0.4)),
                        label=f"{reg} (n={len(g)}): π₁ = {pi1:.2f} {_ci(lo, hi)}")
            ax.axhline(1, color="k", ls="--", lw=0.8, label="no neuron with an effect")
            ax.set(xlabel="shuffle p", ylabel="density", xlim=(0, 1),
                   title=f"A. p-value distribution — {plotting.DIR_LABEL[d]}")
            ax.legend(fontsize=7, frameon=False)
            ax = axs[r_, 1]                                       # B: ECDF of excess effect
            allv = R[col].dropna()
            lo_x, hi_x = (np.nanpercentile(allv, [1, 99.5]) if len(allv) else (0, 1))
            for reg in regions:
                g = R[R.region == reg]
                c = config.CLR.get(reg, (0.4, 0.4, 0.4))
                for sel, ls, nm in ((g[f"sig_{d}"], "-", "FDR-significant"),
                                    (np.ones(len(g), bool), "--", "all")):
                    x = np.sort(g[sel][col].dropna().to_numpy())
                    if x.size:
                        ax.step(x, np.arange(1, x.size + 1) / x.size, where="post", color=c, ls=ls,
                                lw=1.6 if ls == "-" else 1.0,
                                label=f"{reg} {nm} (n={x.size}, median {np.median(x):.3g})")
            ax.axvline(0, color="k", ls=":", lw=0.8)
            ax.set(xlim=(min(lo_x, 0), hi_x), xlabel=f"excess effect: {unit}",
                   ylabel="cumulative fraction of neurons",
                   title="B. how strong? (0 = chance level)")
            ax.legend(fontsize=6, frameon=False, loc="lower right")
            ax = axs[r_, 2]                                       # C: strength vs evidence
            for reg in regions:
                g = R[R.region == reg]
                if not len(g):
                    continue
                c = config.CLR.get(reg, (0.4, 0.4, 0.4))
                lp = -np.log10(g[f"pS_{d}"].clip(lower=1e-12))
                flag = rate_flag_hz and col.startswith("bps")       # only Poisson effects mis-scale
                hi_rate = (g.rate_hz > rate_flag_hz).to_numpy() if flag else np.zeros(len(g), bool)
                ax.scatter(lp[~hi_rate], g[col][~hi_rate], s=9, alpha=0.5, color=c, label=reg)
                if hi_rate.any():
                    ax.scatter(lp[hi_rate], g[col][hi_rate], s=14, marker="x", color=c,
                               label=f"{reg} > {rate_flag_hz:g} Hz (effect mis-scaled)")
                sig = g[f"sig_{d}"].to_numpy()
                if sig.any():
                    ax.axvline(-np.log10(g.loc[sig, f"pS_{d}"].max()), color=c, lw=1, ls="--")
                e = f"exp_{d}"
                if e in g and g[e].any():
                    ax.scatter(lp[g[e]], g[col][g[e]], s=60, facecolors="none", edgecolors="k",
                               lw=0.8, label="truly coupled (sim)")
            ax.axhline(0, color="k", ls=":", lw=0.8)
            ax.set(ylim=(min(lo_x, 0), hi_x), xlabel="-log10 shuffle p (evidence)",
                   ylabel=f"excess effect: {unit}",
                   title="C. strength vs evidence (dashed = FDR cutoff per region)")
            ax.legend(fontsize=6, frameon=False)
        fig.suptitle(f"{title}  |  effect strength & fraction with an effect"
                     f"{' — ' + label(v) if label(v) else ''}  (compare regions within a row only)",
                     fontsize=9)
        plotting.finish(fig, fig_dir, fig_name("8_effect_strength", v), save, show)


def plot_pi1(RES, setting_cols, regions, fig_dir, fig_name, save, show, title, label=lambda v: "",
             null_res=None, rate_flag_hz=100):
    """8c_fraction_with_effect: pi1 (Storey; estimated fraction of neurons with an effect,
    including effects too weak for FDR) with its 95% SESSION-bootstrap CI, per region (x)
    and direction (panels). Filled = all neurons; open = only neurons <= rate_flag_hz.
    Grey band = pi1 [95% CI] of a simulated NULL run (null_res: that run's neuron table),
    i.e. what pi1 reads when no neuron has an effect -- read real pi1 against it."""
    import matplotlib.pyplot as plt
    for v in settings_of(RES, setting_cols):
        R = select(RES, setting_cols, v)
        fig, axs = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
        for ax, d in zip(axs, ("xy", "yx")):
            if null_res is not None and f"pS_{d}" in null_res:
                N = null_res if not set(setting_cols) <= set(null_res) else (
                    select(null_res, setting_cols, v) if len(select(null_res, setting_cols, v)) else null_res)
                p0, lo0, hi0 = stats.storey_pi1(N[f"pS_{d}"], groups=N["tag"])
                ax.axhspan(lo0 if np.isfinite(lo0) else p0, hi0 if np.isfinite(hi0) else p0,
                           color="0.85", zorder=0, label=f"simulated null: {p0:.2f} {_ci(lo0, hi0)}")
                ax.axhline(p0, color="0.5", lw=1, zorder=0)
            for j, reg in enumerate(regions):
                g = R[R.region == reg]
                if not len(g):
                    continue
                c = config.CLR.get(reg, (0.4, 0.4, 0.4))
                for dx, sub, mfc, nm in ((-0.1, g, c, "all"),
                                         (0.1, g[g.rate_hz <= rate_flag_hz], "w", f"<= {rate_flag_hz:g} Hz")):
                    if len(sub) < 3:
                        continue
                    p, lo, hi = stats.storey_pi1(sub[f"pS_{d}"], groups=sub["tag"])
                    yerr = ([[p - lo], [hi - p]] if np.isfinite(lo) else None)
                    ax.errorbar(j + dx, p, yerr=yerr, fmt="o", color=c, mfc=mfc, ms=8, capsize=4, lw=2,
                                label=f"{reg} {nm} (n={len(sub)}): {p:.2f} {_ci(lo, hi)}")
            ax.set_xticks(range(len(regions)), regions)
            ax.set(xlim=(-0.6, len(regions) - 0.4), ylim=(-0.02, 1.02), title=plotting.DIR_LABEL[d],
                   ylabel="π₁ = estimated fraction of neurons with an effect" if d == "xy" else None)
            ax.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.1))
        fig.suptitle(f"{title}\nfraction of neurons with an effect, 95% CI resampling SESSIONS"
                     + (f"  |  {label(v)}" if label(v) else ""), fontsize=9)
        plotting.finish(fig, fig_dir, fig_name("8c_fraction_with_effect", v), save, show)


def plot_null_vs_observed(pairs, regions, units, fig_dir, name, save, show, title):
    """Observed effect vs the trial-shuffle null, per region (rows) x direction (columns).
      pairs[(region, d)] = list of (observed value, array of that neuron's shuffle values)
      units[d]           = x-axis label
    Grey = shuffle values pooled over the region's neurons (each neuron weighs 1); coloured
    line = observed values; dashed = medians. Log x axis. For intuition only: significance is
    always each neuron against ITS OWN shuffles (the pooled null mixes sessions with different
    trial counts)."""
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(len(regions), 2, figsize=(12, 4 * len(regions)), squeeze=False)
    for c_, d in enumerate(("xy", "yx")):
        vals = [np.atleast_1d(np.asarray(x, float)) for r in regions
                for o_, n_ in pairs.get((r, d), []) for x in (o_, n_)]
        allv = np.concatenate(vals) if vals else np.zeros(0)
        allv = allv[np.isfinite(allv) & (allv > 0)]
        edges = (np.logspace(np.log10(np.percentile(allv, 0.1)), np.log10(allv.max()), 50)
                 if allv.size else np.logspace(-4, 0, 50))
        for r_, reg in enumerate(regions):
            ax = axs[r_, c_]
            P = pairs.get((reg, d), [])
            if not P:
                ax.text(0.5, 0.5, "no neurons", ha="center", va="center", transform=ax.transAxes, color="0.4")
                ax.set(title=f"{reg} — {plotting.DIR_LABEL[d]}")
                continue
            o = np.array([v for v, _ in P], float)
            N = [np.asarray(n, float) for _, n in P if len(n)]
            nn = np.concatenate(N)
            w = np.concatenate([np.full(len(x), 1 / (len(x) * len(N))) for x in N])
            ax.hist(nn, bins=edges, weights=w, color="0.6", alpha=0.7,
                    label=f"shuffle null ({len(N)} neurons x {len(N[0])} shuffles)")
            ax.hist(o, bins=edges, weights=np.full(o.size, 1 / o.size), histtype="step", lw=2,
                    color=config.CLR.get(reg, (0.4, 0.4, 0.4)), label=f"observed ({reg}, n={o.size})")
            ax.axvline(np.median(nn), color="0.3", ls="--", lw=1, label=f"null median {np.median(nn):.1e}")
            ax.axvline(np.median(o), color=config.CLR.get(reg, (0.4, 0.4, 0.4)), ls="--", lw=1.4,
                       label=f"observed median {np.median(o):.1e}")
            ax.set_xscale("log")
            ax.set(xlim=(edges[0], edges[-1]), xlabel=units[d], ylabel="fraction of neurons",
                   title=f"{reg} — {plotting.DIR_LABEL[d]}")
            ax.legend(fontsize=7, frameon=False)
    fig.suptitle(f"{title}  |  observed vs trial-shuffle null (grey = chance level; each neuron "
                 f"weighted equally; compare within a panel)", fontsize=9)
    plotting.finish(fig, fig_dir, name, save, show)


def print_lag_summary(RES, setting_cols, regions, keys, prefix, label=lambda v: "",
                      unreliable=()):
    """Median [IQR] of lag columns lag_<key>_<d> among FDR-significant neurons, per region,
    with a CP vs CB Mann-Whitney test."""
    for v in settings_of(RES, setting_cols):
        R = select(RES, setting_cols, v)
        for d in ("xy", "yx"):
            for key in keys:
                col = f"lag_{key}_{d}"
                if col not in R:
                    continue
                parts, vals = [], {}
                for reg in regions:
                    x = R[(R.region == reg) & R[f"sig_{d}"]][col].dropna()
                    vals[reg] = x
                    if len(x):
                        parts.append(f"{reg} median {x.median():.0f} ms "
                                     f"[IQR {x.quantile(.25):.0f}-{x.quantile(.75):.0f}] n={len(x)}")
                mw = ""
                if all(len(vals.get(r, [])) >= 3 for r in ("CP", "CB")):
                    mw = f"  | CP vs CB Mann-Whitney p = {mannwhitneyu(vals['CP'], vals['CB']).pvalue:.2g}"
                tag = f"[{label(v)}] " if label(v) else ""
                print(f"  {prefix} {tag}{plotting.DIR_LABEL[d]}, {key}"
                      f"{' (UNRELIABLE, see note)' if key in unreliable else ''}: "
                      f"{'; '.join(parts) or 'no significant neurons'}{mw}")


def sim_truth_check(RES, setting_cols, label=lambda v: "", lag_keys=(), true_lag=None):
    """Detection vs simulation truth (p_shuffle < 0.05). Rows expected False are
    false-positive rates. -> list of (label, k, n, expected) for plotting.detection_bars."""
    print("\nsimulation truth check (p_shuffle < 0.05; expected-False rows are FALSE-POSITIVE rates):")
    det = []
    for v in settings_of(RES, setting_cols):
        R = select(RES, setting_cols, v)
        for d, e in (("xy", "exp_xy"), ("yx", "exp_yx")):
            for exp in (False, True):
                m = R[e] == exp
                k, n = int((R.loc[m, f"pS_{d}"] < 0.05).sum()), int(m.sum())
                if not n:
                    continue
                lo, hi = stats.binom_ci(k, n)
                lab = f"{label(v)} {d}".strip()
                det.append((f"{lab} {'TRUE' if exp else 'null'}", k, n, exp))
                lagtxt = ""
                if exp and lag_keys:
                    lagtxt = (f" | recovered lag (true {true_lag:g} ms): " +
                              ", ".join(f"{kk} median {R.loc[m, f'lag_{kk}_{d}'].median():.0f}"
                                        for kk in lag_keys if f"lag_{kk}_{d}" in R) + " ms")
                print(f"  {lab:24s} expected {str(exp):5s}: {k:4d}/{n:<4d} = {k/n:.3f} "
                      f"[{lo:.3f}, {hi:.3f}]{lagtxt}")
    return det


def save_tables(obj, RES, lag_tag, out_subdir, run_name):
    base = config.var_path(lag_tag, out_subdir, f"{run_name}.pkl")
    save_pickle(dict(obj, RES=RES), base)
    csv = base.with_name(f"{run_name}_neurons.csv")
    RES.to_csv(_long(csv), index=False)
    print(f"\nsaved -> {base}\n         {csv}   (one row per neuron x setting; open in Excel)")


# ---------------------------------------------------------------------------
# figures shared by every script
# ---------------------------------------------------------------------------
def plot_direction_figures(RES, setting_cols, regions, n_shuffle, q, fig_dir, fig_name,
                           save, show, title, label=lambda v: ""):
    """1_direction_classes (stacked bars per region) and 2_direction_map (per region)."""
    import matplotlib.pyplot as plt
    sets = settings_of(RES, setting_cols)
    fig, axs = plt.subplots(1, len(sets), figsize=(4.2 * len(sets) + 1.5, 4), squeeze=False)
    for ax, v in zip(axs[0], sets):
        plotting.direction_classes(ax, select(RES, setting_cols, v), regions,
                                   title=f"{label(v)}\nFDR q<{q} per direction".strip())
    fig.suptitle(f"{title} direction classes", fontsize=9)
    plotting.finish(fig, fig_dir, fig_name("1_direction_classes"), save, show)
    for v in sets:
        R = select(RES, setting_cols, v)
        fig, axs = plt.subplots(1, len(regions), figsize=(5.4 * len(regions), 5), squeeze=False)
        for ax, reg in zip(axs[0], regions):
            g = R[R.region == reg]
            plotting.direction_map(ax, g, n_shuffle, title=f"{reg}  (n={len(g)})",
                                   highlight=g["driver"].to_numpy() if "driver" in g else None)
        fig.suptitle(f"shuffle p per neuron{', ' + label(v) if label(v) else ''}  "
                     f"(solid = FDR threshold, dashed = p 0.05)", fontsize=9)
        plotting.finish(fig, fig_dir, fig_name("2_direction_map", v), save, show)


TRUE_CLASS_LABEL = {"kin→neuron": "neuron LAGS\n(kin→neuron)", "neuron→kin": "neuron LEADS\n(neuron→kin)",
                    "both": "both", "none": "no coupling"}


def sim_truth_table(R):
    """Confusion counts: true class (from the simulation) x detected FDR class."""
    true = plotting.direction_class(R["exp_xy"], R["exp_yx"])
    tab = pd.crosstab(pd.Series(true, name="true"), pd.Series(R["direction"].to_numpy(), name="detected"))
    rows = [c for c in plotting.DIR_CLASSES if c in tab.index]
    return tab.reindex(index=rows, columns=list(plotting.DIR_CLASSES), fill_value=0)


def plot_sim_truth(RES, setting_cols, lag_keys, true_lag, fig_dir, fig_name, save, show,
                   label=lambda v: ""):
    """Simulations only. Per setting: (A) where each neuron SHOULD land (true class, rows)
    vs where it DOES land (detected FDR class, columns), counts and row %; (B) recovered lag
    of the truly coupled neurons vs the true lag, per direction (FDR-significant filled,
    missed open)."""
    import matplotlib.pyplot as plt
    for v in settings_of(RES, setting_cols):
        R = select(RES, setting_cols, v)
        tab = sim_truth_table(R)
        print(f"\nwhere simulated neurons land [{label(v)}] (rows = truth, columns = detected, FDR):"
              if label(v) else "\nwhere simulated neurons land (rows = truth, columns = detected, FDR):")
        print(tab.to_string())
        fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), gridspec_kw=dict(width_ratios=[1.2, 1, 1]))
        ax = axs[0]
        pct = tab.div(tab.sum(1), axis=0).to_numpy()
        ax.imshow(pct, cmap="Greens", vmin=0, vmax=1, aspect="auto")
        for i in range(tab.shape[0]):
            for j in range(tab.shape[1]):
                ok = tab.index[i] == tab.columns[j]
                ax.text(j, i, f"{tab.iat[i, j]}\n{100*pct[i, j]:.0f}%", ha="center", va="center",
                        fontsize=8, fontweight="bold" if ok else "normal",
                        color="w" if pct[i, j] > 0.6 else "k")
        ax.set_xticks(range(tab.shape[1]), [TRUE_CLASS_LABEL[c] for c in tab.columns], fontsize=7)
        ax.set_yticks(range(tab.shape[0]), [TRUE_CLASS_LABEL[c] for c in tab.index], fontsize=7)
        ax.set(xlabel="detected (BH-FDR)", ylabel="simulated truth",
               title="where neurons SHOULD land (row) vs DO land (column)")
        for ax, d, e in ((axs[1], "xy", "exp_xy"), (axs[2], "yx", "exp_yx")):
            g = R[R[e]]
            if not len(g) or not lag_keys:
                ax.text(0.5, 0.5, "no truly coupled neurons\nin this direction", ha="center",
                        va="center", transform=ax.transAxes, fontsize=8, color="0.4")
                ax.set(title=plotting.DIR_LABEL[d]); ax.set_xticks([]); ax.set_yticks([])
                continue
            for k, key in enumerate(lag_keys):
                col = f"lag_{key}_{d}"
                if col not in g:
                    continue
                x = k + np.random.default_rng(k).uniform(-0.15, 0.15, len(g))
                sig = g[f"sig_{d}"].to_numpy()
                ax.scatter(x[sig], g[col][sig], s=18, color="tab:blue", label="detected" if k == 0 else None)
                ax.scatter(x[~sig], g[col][~sig], s=18, facecolors="none", edgecolors="tab:red",
                           label="missed" if k == 0 else None)
                ax.plot([k - 0.3, k + 0.3], [g[col].median()] * 2, color="k", lw=2)
            ax.axhline(true_lag, color="tab:green", ls="--", lw=1.2, label=f"true lag {true_lag:g} ms")
            ax.set_xticks(range(len(lag_keys)), lag_keys)
            ax.set(ylabel="recovered lag (ms)", ylim=(0, config.CROSS_LAG_RANGE_MS[1] + 5),
                   title=f"{plotting.DIR_LABEL[d]}: true vs recovered lag (n={len(g)})")
            ax.legend(fontsize=7, frameon=False)
        fig.suptitle(f"simulation truth check{' — ' + label(v) if label(v) else ''}  "
                     "(black bar = median recovered lag)", fontsize=9)
        plotting.finish(fig, fig_dir, fig_name("6b_sim_truth", v), save, show)


def plot_detection(det, scenario, fig_dir, name, save, show):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(1 + 0.6 * len(det), 4))
    plotting.detection_bars(ax, det, title=f"{scenario}: shuffle-test detection (red = true effect)")
    plotting.finish(fig, fig_dir, name, save, show)
