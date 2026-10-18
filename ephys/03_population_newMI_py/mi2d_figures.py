"""mi2d_figures.py - the analysis and figures of the TIME x TIME MI maps.

One function per analysis step, taken from 12_MI_2D_kin_x_neuron.py cells 4b-7:

    hlag_per_mouse             4b  h(lag) per mouse, robust peak, median lag, insets
    diagonal_vs_neural         5   per neuron: diagonal (lag) vs neural-time R2
    individual_landscapes      5b  one figure per example neuron (slow)
    diagonal_vs_neural_figure  6   the summary figure of 5
    per_mouse_maps             7   mean maps per mouse and for all mice
    timing_summary                 where the time went
    save_results                   the run's results as one pickle

Every function takes `ctx`, the run's shared inputs (a SimpleNamespace): the
maps EXC2 (neurons x kinematic windows x neural windows, excess MI), which
neurons are used (keep, uids, run_idx, regions_run), the window grids
(win_centers, kin_centers, n_kin, n_wins, WIN_STEP, KIN_WIN, KIN_STEP), the
colours CLR, the figure title `lab`, the output folder `fdir` and RUN_NAME, the
save/show switches, the map drawing settings (extent, ASPECT, CMAP), the TIMES
dict, and `zero_txt`, the name of t = 0 written on the axes ("reach endpoint",
"the cue (detection)", "the force onset").

Used by 18_MI_2D_pushpull.py. 12 still runs its own copy of these cells (to be
switched to these functions), so until then A CHANGE HERE MUST ALSO BE MADE IN
12'S CELLS 4b-7, or the reach and push/pull analyses stop being comparable.
"""
import time, pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import mannwhitneyu, wilcoxon
import config


def _savefig_retry(fig, path, tries=6, **kw):
    """The output tree is Dropbox-synced: Dropbox briefly locks a file while it
    uploads it, and overwriting one then raises OSError 22. Retry, do not crash."""
    for t in range(tries):
        try:
            return fig.savefig(path, **kw)
        except OSError:
            if t == tries - 1: raise
            time.sleep(0.5 * (t + 1))


def shade_partial(ax, lag_axis, lag_full, scale=1000.0):
    """grey where a diagonal average uses only SOME kinematic windows"""
    lo_, hi_ = lag_axis[lag_full][[0, -1]] * scale
    ax.axvspan(lag_axis[0] * scale, lo_, color="0.92", zorder=0, lw=0)
    ax.axvspan(hi_, lag_axis[-1] * scale, color="0.92", zorder=0, lw=0)


def hlag_per_mouse(ctx, LAG_BIN=None, PER_MOUSE_MIN_N=5):
    """12's cell 4b: h(lag) of the mean landscape per mouse -- z-scored curves,
    robust peak and median lag with paired CB-vs-CP insets, the full two-row
    figure and the presentation figure (PNG + PDF). Returns the per-mouse
    curves and timing measures."""
    EXC2 = ctx.EXC2
    keep = ctx.keep
    uids = ctx.uids
    run_idx = ctx.run_idx
    regions_run = ctx.regions_run
    win_centers = ctx.win_centers
    kin_centers = ctx.kin_centers
    n_kin = ctx.n_kin
    n_wins = ctx.n_wins
    WIN_STEP = ctx.WIN_STEP
    CLR = ctx.CLR
    lab = ctx.lab
    fdir = ctx.fdir
    RUN_NAME = ctx.RUN_NAME
    SAVE_FIGS = ctx.SAVE_FIGS
    SHOW_FIGS = ctx.SHOW_FIGS
    MED4B = PE4B = None          # set below when there is something to plot

    # 4b. cell --> h(lag): each region's mean map averaged ALONG each line parallel
    # to the dashed diagonal. Every cell on one such line has the same
    # lag = neural - kinematic, so each line becomes one point of the curve
    # (the same h(lag) as cell 5, here on the mean map).
    #   thick line   all neurons of the region (mean +- SEM over neurons)
    #   thin lines   one per mouse (the mean over that mouse's neurons), so it is
    #                visible whether the population curve is shared by the animals
    #                or carried by one of them. Mice with fewer than PER_MOUSE_MIN_N
    #                neurons in the region are left out.
    # Lags are binned at LAG_BIN. Near the corners of the map a line crosses only a
    # corner and averages few cells: those bins are shaded grey (fewer than half
    # the cells of the fullest bin).
    _bin4b = WIN_STEP if LAG_BIN is None else LAG_BIN
    _NN, _KK = np.meshgrid(win_centers, kin_centers)        # (n_kin, n_wins), seconds
    _LAG4B = _NN - _KK                                       # constant along each diagonal
    _ib4b  = np.round(_LAG4B / _bin4b).astype(int)
    _i04b  = _ib4b.min(); _ib4b = _ib4b - _i04b
    LAG4B  = (np.arange(_ib4b.max() + 1) + _i04b) * _bin4b    # bin centres (s)
    NCELL4B = np.bincount(_ib4b.ravel(), minlength=LAG4B.size)


    def _h_lag(M):
        """h(lag) of one map: mean over the cells sharing each binned lag"""
        ok = np.isfinite(M)
        n  = np.bincount(_ib4b[ok], minlength=LAG4B.size)
        sm = np.bincount(_ib4b[ok], M[ok], minlength=LAG4B.size)
        return np.where(n > 0, sm / np.maximum(n, 1), np.nan)


    _mouse_of = {pi: uids[run_idx[pi]].split("|")[1] for pi in keep}
    _MICE4B   = sorted(set(_mouse_of.values()))
    _MLS      = ["-", "--", ":", "-.", (0, (3, 1, 1, 1)), (0, (5, 2))]
    H4B = {}                # region -> dict(all=(mean, sem, n), mice={mouse: (mean, n)})
    for reg in ("CP", "CB"):
        idx = [pi for pi in keep if regions_run[pi] == reg]
        if not idx: continue
        P = np.stack([_h_lag(EXC2[pi]) for pi in idx])
        H4B[reg] = dict(all=(np.nanmean(P, axis=0),
                             np.nanstd(P, axis=0) / np.sqrt(np.maximum(np.isfinite(P).sum(0), 1)),
                             len(idx)),
                        mice={})
        for m in _MICE4B:
            sel = [j for j, pi in enumerate(idx) if _mouse_of[pi] == m]
            if len(sel) >= PER_MOUSE_MIN_N:
                H4B[reg]["mice"][m] = (np.nanmean(P[sel], axis=0), len(sel))


    def _shade_low(ax):
        for x, low in zip(LAG4B, NCELL4B < 0.5 * NCELL4B.max()):
            if low: ax.axvspan((x - _bin4b / 2) * 1000, (x + _bin4b / 2) * 1000,
                               color="0.92", lw=0, zorder=0)


    if H4B:
        fig, axes = plt.subplots(1, 3, figsize=(19, 5.2), sharey=True)
        for ax, reg in zip(axes[:2], ("CP", "CB")):
            if reg not in H4B: ax.set_visible(False); continue
            _shade_low(ax)
            for j, (m, (mn, n)) in enumerate(H4B[reg]["mice"].items()):
                ax.plot(LAG4B * 1000, mn, color=CLR[reg], lw=1.0, alpha=0.5,
                        ls=_MLS[j % len(_MLS)], label=f"{m} (n={n})")
            mn, se, n = H4B[reg]["all"]
            ax.plot(LAG4B * 1000, mn, color=CLR[reg], lw=3, label=f"all neurons (n={n})")
            ax.fill_between(LAG4B * 1000, mn - se, mn + se, color=CLR[reg], alpha=0.25, lw=0)
            ax.set(title=f"{reg}: h(lag), all neurons (thick) and each mouse (thin)",
                   xlabel="lag = neural - kinematic (ms)     < 0: spikes BEFORE the movement")
            ax.title.set_fontsize(10)
            ax.legend(fontsize=7, loc="upper right")
        ax = axes[2]; _shade_low(ax)
        for reg in ("CP", "CB"):
            if reg not in H4B: continue
            mn, se, n = H4B[reg]["all"]
            ax.plot(LAG4B * 1000, mn, color=CLR[reg], lw=3, label=f"{reg} (n={n})")
            ax.fill_between(LAG4B * 1000, mn - se, mn + se, color=CLR[reg], alpha=0.25, lw=0)
        ax.set(title="CP and CB, all neurons", xlabel="lag = neural - kinematic (ms)")
        ax.title.set_fontsize(10); ax.legend(fontsize=8)
        for ax in axes:
            ax.axvline(0, color="k", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.5)
            ax.spines[["top", "right"]].set_visible(False)
        axes[0].set_ylabel("mean excess MI (bits)")
        fig.suptitle(f"{lab}  —  h(lag): the mean landscape averaged along each line "
                     f"parallel to the diagonal (bins {_bin4b*1000:.0f} ms; grey = lags "
                     f"built from few cells, near the corners of the map)", fontsize=11)
        plt.tight_layout(rect=(0, 0, 1, 0.93))
        if SAVE_FIGS:
            _savefig_retry(fig, fdir / f"{RUN_NAME}_hlag_per_mouse.png", dpi=150,
                           bbox_inches="tight")
        if SHOW_FIGS: plt.show()
        plt.close(fig)

        # the same curves on ONE axis: every mouse a thin solid line in its region's
        # colour (no identity, no SEM), and each region's all-neuron curve thick.
        # The curves below are Z-SCORED (_norm): each curve minus its own mean, over
        #           its own SD. Uses
        #           every bin, so one noisy bin cannot set the scale -> the robust
        #           way to compare WHERE the curves peak. The z-score uses the
        #           well-covered lags only (not the grey corners, which average few
        #           cells and would otherwise set the scale), and only those lags
        #           are drawn.
        #
        # Each curve's MEDIAN LAG is marked (two-row figure, bottom): the lag that splits
        # the area under the curve in half -- half of the curve's information sits at
        # shorter lags, half at longer ones. Unlike the peak it uses the whole curve,
        # so one noisy bin cannot move it. The area is taken ABOVE the curve's own
        # baseline (its MEDIAN_BASE_PCT percentile over the well-covered lags):
        # without that, a curve riding on a high floor (one CP mouse never drops
        # below ~0.03 bits) would have its median pulled to the middle of the lag
        # range by the floor alone. Because of that baseline, the median lag is the
        # same whether it is measured on the curve in bits or z-scored.
        MEDIAN_BASE_PCT = 10
        _full = NCELL4B >= 0.5 * NCELL4B.max()

        def _median_lag(h):
            """lag (s) splitting the area of h above its baseline in half, over the
            well-covered lags; linear interpolation between bins"""
            v, x = h[_full], LAG4B[_full]
            c = np.clip(np.nan_to_num(v - np.nanpercentile(v, MEDIAN_BASE_PCT)), 0, None)
            if c.sum() <= 0: return np.nan
            cs = np.r_[0.0, np.cumsum(c) / c.sum()]
            return float(np.interp(0.5, cs, np.r_[x[0] - _bin4b / 2, x + _bin4b / 2]))

        MED4B = {reg: dict(all=_median_lag(H4B[reg]["all"][0]),
                           mice={m: _median_lag(mn) for m, (mn, _) in H4B[reg]["mice"].items()})
                 for reg in H4B}
        print(f"\nmedian lag of h(lag) (ms; splits the area above the curve's "
              f"{MEDIAN_BASE_PCT}th-percentile baseline in half, lags "
              f"{LAG4B[_full][0]*1000:.0f}..{LAG4B[_full][-1]*1000:.0f} ms):")
        for reg, d in MED4B.items():
            print(f"  {reg}: all neurons {d['all']*1000:+5.0f}   |   "
                  + "  ".join(f"{m.split('_')[-1]} {v*1000:+.0f}" for m, v in d["mice"].items()))

        def _norm(h, mode):
            if mode == "bits": return h
            v = h[_full]
            mu, sd = np.nanmean(v), np.nanstd(v)
            return (h - mu) / sd if sd > 0 else h * np.nan

        # ── the same z-scored curves, marking each one's PEAK and ROBUST PEAK ────
        #   peak         the lag where the curve (smoothed over PEAK_SMOOTH bins) is
        #                highest -- decided by ONE bin, so noise can move it, and it
        #                can only take values on the 25 ms lag grid.
        #   robust peak  the CENTRE OF MASS of the curve's peak region: every lag
        #                where the z-scored curve is above PEAKREG_Z, weighted by how
        #                far above it is. No single bin decides it, a noisy maximum
        #                barely moves it, and it is not stuck on the lag grid.
        #                PEAKREG_Z = 1 means "at least one SD above the curve's own
        #                mean". On a skewed curve it sits slightly towards the
        #                longer tail of the peak, which is part of what it measures.
        # Only the well-covered lags are used.
        PEAK_SMOOTH = 3
        PEAKREG_Z    = 1.0

        def _smooth(v, w=PEAK_SMOOTH):
            if w <= 1: return v.copy()
            k = np.ones(w) / w
            num = np.convolve(np.nan_to_num(v), k, mode="same")
            den = np.convolve(np.isfinite(v).astype(float), k, mode="same")
            return np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)

        def _peaks(h):
            """-> (peak lag, robust peak lag) in s, on the well-covered lags"""
            x  = LAG4B[_full]
            ip = int(np.nanargmax(_smooth(h[_full])))
            z  = _norm(h, "z")[_full]
            w  = np.clip(np.nan_to_num(z - PEAKREG_Z), 0, None)
            return x[ip], (float(np.sum(w * x) / w.sum()) if w.sum() > 0 else np.nan)

        PE4B = {reg: dict(all=_peaks(H4B[reg]["all"][0]),
                          mice={m: _peaks(mn) for m, (mn, _) in H4B[reg]["mice"].items()})
                for reg in H4B}
        print(f"\npeak / robust peak of h(lag) (ms; peak on the curve smoothed over "
              f"{PEAK_SMOOTH} bins; robust peak = centre of mass of the curve above "
              f"z = {PEAKREG_Z:g}; lags {LAG4B[_full][0]*1000:.0f}..{LAG4B[_full][-1]*1000:.0f} ms):")
        for reg, d in PE4B.items():
            print(f"  {reg}: all neurons peak {d['all'][0]*1000:+5.0f}, robust peak "
                  f"{d['all'][1]*1000:+5.0f}   |   "
                  + "  ".join(f"{m.split('_')[-1]} {pk*1000:+.0f}/{el*1000:+.0f}"
                              for m, (pk, el) in d["mice"].items()))

        # one figure each, drawn like the median figure: thin dashed line per mouse,
        # thick dashed line for all neurons
        # ── per-mouse timing measures: CB vs CP, paired ──────────────────────────
        # One point per mouse per region, a box per region, and a line joining each
        # mouse's CB and CP points. Only mice with enough neurons in BOTH regions are
        # paired. The points are jittered sideways (the same offset for a mouse in
        # both regions, so its line stays straight) to keep close values apart.
        # Paired t-test on the per-mouse values, CB vs CP.
        #   robust peak  centre of mass of the curve above z = PEAKREG_Z
        #   median lag   splits the area above the curve's baseline in half
        # Drawn as insets of the z-scored curves (two-row figure: robust peak on top,
        # median lag at the bottom; and the robust peak alone, for presenting).
        from scipy.stats import ttest_rel
        BOX_JITTER  = 0.07
        INSET_SEP   = 0.55     # inset: distance between the CB and CP positions (1 = far apart)
        INSET_BOX_W = 0.22     # inset: box width, in the same units
        _pair = [m for m in _MICE4B if all(m in PE4B.get(r, {}).get("mice", {})
                                           and m in MED4B.get(r, {}).get("mice", {})
                                           for r in ("CP", "CB"))]
        _order = ("CB", "CP")                                    # left to right
        _jit = np.random.default_rng(0).uniform(-BOX_JITTER, BOX_JITTER, len(_pair))
        TIMING4B = {"robust peak": {r: np.array([PE4B[r]["mice"][m][1] for m in _pair]) * 1000
                                    for r in _order},
                    "median lag":  {r: np.array([MED4B[r]["mice"][m] for m in _pair]) * 1000
                                    for r in _order}}

        def _stars(p_):
            """significance marker for a p value"""
            return ("***" if p_ < 0.001 else "**" if p_ < 0.01 else "*" if p_ < 0.05
                    else "n.s.")

        def _paired_box(ax, vals, name, labels=True, compact=False, stars=False):
            """box + jittered points + a line per mouse, CB then CP; paired t-test.
            compact: small version for an inset (no names, one-line test, small
            fonts). stars: (compact only) mark the test with a bracket and stars
            (* p<0.05, ** p<0.01, *** p<0.001, n.s.) instead of writing t and p.
            -> (t, p)"""
            # the inset (compact) is drawn square, with the two regions closer
            # together and narrower boxes; the jitter shrinks with the spacing
            _x1 = INSET_SEP if compact else 1.0
            _jt = _jit * (_x1 / 1.0) * (INSET_BOX_W / 0.45 if compact else 1.0)
            bp = ax.boxplot([vals[r] for r in _order], positions=[0, _x1],
                            widths=INSET_BOX_W if compact else 0.45,
                            patch_artist=True, showfliers=False, zorder=1)
            for patch, r in zip(bp["boxes"], _order):
                patch.set(facecolor=CLR[r], alpha=0.35, edgecolor="none")
            for k_ in ("whiskers", "caps", "medians"):
                for ln in bp[k_]: ln.set(color="0.3", lw=1.2)
            for j in range(len(_pair)):
                ax.plot([0 + _jt[j], _x1 + _jt[j]], [vals["CB"][j], vals["CP"][j]],
                        color="0.55", lw=1.0 if compact else 1.2, zorder=2)
            for x_, r in zip((0, _x1), _order):
                ax.plot(x_ + _jt, vals[r], "o", ms=5 if compact else 8, color=CLR[r],
                        mec="w", mew=0.8 if compact else 1.0, zorder=3)
            if labels and not compact:
                # mouse names beside the CP points, pushed apart where points are
                # close, with a thin leader back to each point
                _all_y = np.r_[vals["CP"], vals["CB"]]
                _gap = 0.05 * (np.nanmax(_all_y) - np.nanmin(_all_y) + 1e-9)
                _o = np.argsort(vals["CP"]); _ly = vals["CP"][_o].astype(float).copy()
                for i_ in range(1, _ly.size): _ly[i_] = max(_ly[i_], _ly[i_ - 1] + _gap)
                _ly -= (_ly.mean() - vals["CP"][_o].mean())
                for y_lab, j in zip(_ly, _o):
                    ax.plot([1 + _jit[j] + 0.05, 1.22], [vals["CP"][j], y_lab], color="0.75",
                            lw=0.6, zorder=1)
                    ax.text(1.25, y_lab, _pair[j].split("_")[-1], va="center", fontsize=7.5,
                            color="0.35")
            d = vals["CB"] - vals["CP"]
            t_, p_ = ttest_rel(vals["CB"], vals["CP"]) if len(_pair) >= 2 else (np.nan, np.nan)
            if compact:
                ax.axhline(0, color="k", lw=0.6, ls=":")
                ax.set(xticks=[0, _x1], xticklabels=list(_order),
                       xlim=(-0.6 * _x1, 1.6 * _x1))
                ax.set_box_aspect(1)                     # square inset
                ax.tick_params(labelsize=7)
                if stars:
                    _top = np.nanmax(np.r_[vals["CB"], vals["CP"]])
                    _rng = _top - np.nanmin(np.r_[vals["CB"], vals["CP"]])
                    _yb = _top + 0.10 * _rng
                    ax.plot([0, 0, _x1, _x1], [_yb - 0.04 * _rng, _yb, _yb, _yb - 0.04 * _rng],
                            color="k", lw=0.8)
                    ax.text(_x1 / 2, _yb, _stars(p_), ha="center", va="bottom", fontsize=10)
                    ax.set_ylim(top=_yb + 0.22 * _rng)
                    ax.set_title(name, fontsize=7.5)
                else:
                    ax.set_title(f"{name}, {len(_pair)} mice" + chr(10)
                                 + f"paired t({len(_pair)-1}) = {t_:.2f}, p = {p_:.3g}",
                                 fontsize=7)
                ax.spines[["top", "right"]].set_visible(False)
                return t_, p_
            ax.text(0.5, 1.0, f"paired t-test, CB vs CP: t({len(_pair)-1}) = {t_:.2f}, "
                              f"p = {p_:.3g}" + chr(10)
                              + f"CB - CP: mean {np.mean(d):+.0f} ms, CB earlier in "
                                f"{int((d < 0).sum())}/{len(d)} mice",
                    transform=ax.transAxes, ha="center", va="bottom", fontsize=8)
            ax.axhline(0, color="k", lw=0.8, ls=":")
            ax.set(xticks=[0, 1], xticklabels=list(_order), xlim=(-0.5, 1.75 if labels else 1.5))
            ax.set_title(f"{name} per mouse (n={len(_pair)} mice)", fontsize=10, pad=34)
            ax.spines[["top", "right"]].set_visible(False)
            print(f"  {name:12s} CB - CP (ms): "
                  + ", ".join(f"{m.split('_')[-1]} {x:+.0f}" for m, x in zip(_pair, d))
                  + f"   | mean {np.mean(d):+.0f}, median {np.median(d):+.0f}; CB earlier in "
                    f"{int((d < 0).sum())}/{len(d)}; paired t({len(_pair)-1}) = {t_:.2f}, "
                    f"p = {p_:.3g}")
            return t_, p_

        if _pair:
            _ylab = "lag (ms)" + chr(10) + "< 0: spikes BEFORE the movement"
            print("\nper mouse, CB vs CP:")
            # the z-scored curves with one timing measure marked on every curve
            # (thin dashed per mouse, thick dashed for all neurons) and its per-mouse
            # box plot as an inset -- no mouse identities anywhere
            _VAL = {"robust peak": lambda r, m: PE4B[r]["mice"][m][1],
                    "median lag":  lambda r, m: MED4B[r]["mice"][m]}
            _ALL = {"robust peak": lambda r: PE4B[r]["all"][1],
                    "median lag":  lambda r: MED4B[r]["all"]}
            _DEF = {"robust peak": f"centre of mass of the curve above z = {PEAKREG_Z:g}",
                    "median lag":  f"splits the area above the curve's {MEDIAN_BASE_PCT}th-"
                                   f"percentile baseline in half"}

            def _zplot(ax, name, xlabel=True, present=False):
                """present: for a talk -- no legend (region and neuron count written
                inside the axes, each all-neuron value written next to its line), short
                title and axis labels, no definitions"""
                xs = LAG4B[_full] * 1000
                for reg in ("CP", "CB"):
                    if reg not in H4B: continue
                    for m, (mn, n) in H4B[reg]["mice"].items():
                        ax.plot(xs, _norm(mn, "z")[_full], color=CLR[reg], lw=1.6, alpha=0.45)
                        v = _VAL[name](reg, m)
                        if np.isfinite(v):
                            ax.axvline(v * 1000, color=CLR[reg], lw=1.2, alpha=0.45, ls="--")
                    mn, _, n = H4B[reg]["all"]
                    ax.plot(xs, _norm(mn, "z")[_full], color=CLR[reg], lw=3.2,
                            label=(f"{reg} (n={n} neurons, {len(H4B[reg]['mice'])} mice)"
                                   if present else
                                   f"{reg}: all neurons (n={n})" + chr(10)
                                   + f"thin = each of {len(H4B[reg]['mice'])} mice"))
                for reg in ("CP", "CB"):
                    v = _ALL[name](reg)
                    if np.isfinite(v):
                        ax.axvline(v * 1000, color=CLR[reg], lw=2.6, ls="--",
                                   label=(f"{reg} {name}: {v*1000:+.0f} ms" if present else
                                          f"{reg} {name} (all neurons): {v*1000:+.0f} ms"))
                        if present:
                            # the value next to its line, at the top: the earlier of
                            # the two lines gets its label on its left, the later one
                            # on its right, so they never overlap
                            _other = _ALL[name]("CB" if reg == "CP" else "CP")
                            _left = np.isfinite(_other) and v < _other
                            ax.text(v * 1000 + (-8 if _left else 8), 0.98,
                                    f"{v*1000:+.0f} ms", transform=ax.get_xaxis_transform(),
                                    ha="right" if _left else "left", va="top", fontsize=9,
                                    color=CLR[reg], fontweight="bold")
                ax.axvline(0, color="k", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.5)
                ax.set(ylabel="MI per lag (z-scored)" if present else "h(lag), z-scored",
                       title=(name if present else
                              f"vertical lines = {name.upper()}  ({_DEF[name]})"))
                # The inset goes in whichever UPPER corner the curves leave clearer:
                # where the peak falls depends on the lag range (in 12 the curves are
                # low on the right, in 18 on the left), so a fixed corner would sit on
                # top of them in one script or the other.
                ax.autoscale_view()
                (_x0, _x1), (_y0, _y1) = ax.get_xlim(), ax.get_ylim()
                _xf = (xs - _x0) / (_x1 - _x0)
                _Z = np.stack([_norm(mn, "z")[_full] for reg in H4B
                               for mn in [H4B[reg]["all"][0]]
                                         + [v[0] for v in H4B[reg]["mice"].values()]])
                _zf = (_Z - _y0) / (_y1 - _y0)

                def _clutter(a_, b_):
                    m_ = (_xf >= a_) & (_xf <= b_)
                    return float(np.nanmax(_zf[:, m_])) if m_.any() else 0.0

                _INSET = {"right": [0.70, 0.50, 0.26, 0.40], "left": [0.08, 0.50, 0.26, 0.40]}
                _side = "right" if _clutter(0.68, 0.98) <= _clutter(0.06, 0.36) else "left"
                if present:
                    # which colour is which region, and how many neurons, inside the
                    # axes, in the upper corner the inset does NOT use
                    _tx, _ha = (0.02, "left") if _side == "right" else (0.98, "right")
                    for i_, reg in enumerate(("CP", "CB")):
                        if reg in H4B:
                            ax.text(_tx, 0.97 - 0.07 * i_, f"{reg}  (n = {H4B[reg]['all'][2]} neurons)",
                                    transform=ax.transAxes, ha=_ha, va="top", fontsize=10,
                                    color=CLR[reg], fontweight="bold")
                if xlabel:
                    ax.set_xlabel("lag = neural - kinematic (ms)" if present else
                                  "lag = neural - kinematic (ms)     < 0: spikes BEFORE the movement")
                ax.title.set_fontsize(10)
                ax.spines[["top", "right"]].set_visible(False)
                if not present:
                    ax.legend(fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5),
                              frameon=False)
                _paired_box(ax.inset_axes(_INSET[_side]), TIMING4B[name],
                            f"{name} (ms)", compact=True, stars=present)

            # the robust peak alone, with its inset -- the PRESENTATION figure: no
            # legend (region and neuron count inside the axes, each all-neuron value
            # next to its line), short labels, a significance star instead of t and
            # p in the inset, and no explanatory text. Saved as PDF as well, with the
            # text kept as text (fonttype 42) so it stays editable in Illustrator.
            fig, ax = plt.subplots(figsize=(9.0, 5.6))
            _zplot(ax, "robust peak", present=True)
            ax.set_title("")
            plt.tight_layout()
            if SAVE_FIGS:
                _stem = fdir / f"{RUN_NAME}_hlag_zscored_robustpeak_inset"
                _savefig_retry(fig, f"{_stem}.png", dpi=150, bbox_inches="tight")
                with plt.rc_context({"pdf.fonttype": 42}):
                    _savefig_retry(fig, f"{_stem}.pdf", bbox_inches="tight")
            if SHOW_FIGS: plt.show()
            plt.close(fig)

            # robust peak on top, median lag at the bottom, each with its inset,
            # with the full labels and legend
            fig, axes = plt.subplots(2, 1, figsize=(10.5, 10.4), sharex=True)
            _zplot(axes[0], "robust peak", xlabel=False)
            _zplot(axes[1], "median lag")
            fig.suptitle(lab + f"   (h(lag) z-scored per curve: each mouse thin, all "
                               f"neurons thick; lags {LAG4B[_full][0]*1000:.0f}.."
                               f"{LAG4B[_full][-1]*1000:.0f} ms)", fontsize=9)
            plt.tight_layout()
            if SAVE_FIGS:
                _savefig_retry(fig, fdir / f"{RUN_NAME}_hlag_zscored_robustpeak_median_inset.png",
                               dpi=150, bbox_inches="tight")
            if SHOW_FIGS: plt.show()
            plt.close(fig)

            # the per-mouse numbers and the paired t-tests behind the two insets
            for name in ("robust peak", "median lag"):
                v_ = TIMING4B[name]; d = v_["CB"] - v_["CP"]
                t_, p_ = ttest_rel(v_["CB"], v_["CP"]) if len(_pair) >= 2 else (np.nan, np.nan)
                print(f"  {name:12s} CB - CP (ms): "
                      + ", ".join(f"{m.split('_')[-1]} {x:+.0f}" for m, x in zip(_pair, d))
                      + f"   | mean {np.mean(d):+.0f}, median {np.median(d):+.0f}; CB earlier "
                        f"in {int((d < 0).sum())}/{len(d)}; paired t({len(_pair)-1}) = "
                        f"{t_:.2f}, p = {p_:.3g}")

    return dict(H4B=H4B, LAG4B=LAG4B, MED4B=MED4B, PE4B=PE4B)


def diagonal_vs_neural(ctx, LAG_BIN=None, CV_GAP=None):
    """12's cell 5: per neuron, does the map vary along the DIAGONAL (lag) or
    along the NEURAL axis? Cross-validated R2 of both descriptions. Returns
    the per-neuron results, used by the later figures and saved."""
    EXC2 = ctx.EXC2
    keep = ctx.keep
    regions_run = ctx.regions_run
    win_centers = ctx.win_centers
    kin_centers = ctx.kin_centers
    n_kin = ctx.n_kin
    n_wins = ctx.n_wins
    WIN_STEP = ctx.WIN_STEP
    KIN_WIN = ctx.KIN_WIN
    KIN_STEP = ctx.KIN_STEP

    # 5. cell --> per neuron: does the map vary along the DIAGONAL (lag) or along the
    # NEURAL axis? Both descriptions are built on one half of the kinematic windows
    # and scored on the other half, row-centred. See the block in cell 0.
    _lag_bin = WIN_STEP if LAG_BIN is None else LAG_BIN
    _gap     = KIN_WIN // KIN_STEP if CV_GAP is None else CV_GAP
    if n_kin // 2 - _gap < 2:
        raise ValueError(f"only {n_kin} kinematic windows: with CV_GAP={_gap} each fold "
                         f"trains on {n_kin//2 - _gap}. Widen KIN_RANGE or lower CV_GAP.")

    LAG      = win_centers[None, :] - kin_centers[:, None]    # neural - kinematic (s)
    LAG_IDX  = np.round(LAG / _lag_bin).astype(int)
    lag_axis = np.arange(LAG_IDX.min(), LAG_IDX.max() + 1) * _lag_bin
    LAG_IDX -= LAG_IDX.min()
    n_lag    = lag_axis.size
    # Lags that EVERY kinematic window reaches. Outside this range a diagonal average
    # is built from only some kinematic windows (the corners of the map), so the ends
    # of h(lag) are thinner -- drawn shaded in the figures.
    lag_full = ((lag_axis >= LAG.min(axis=1).max() - _lag_bin / 2) &
                (lag_axis <= LAG.max(axis=1).min() + _lag_bin / 2))


    def _rowc(M):
        """subtract each kinematic row's own mean, ignoring NaN. This is what makes
        the comparison independent of how strong the map is at that kinematic time."""
        ok = np.isfinite(M); ct = ok.sum(axis=1, keepdims=True)
        mu = np.where(ct > 0,
                      np.where(ok, M, 0).sum(axis=1, keepdims=True) / np.maximum(ct, 1),
                      np.nan)
        return M - mu


    def _lag_profile(E, rows=None):
        """h(lag): mean of the map along each diagonal (optionally only some rows)."""
        Er = E if rows is None else E[rows]
        Lr = LAG_IDX if rows is None else LAG_IDX[rows]
        ok = np.isfinite(Er)
        sm = np.bincount(Lr[ok], Er[ok], n_lag)
        ct = np.bincount(Lr[ok], minlength=n_lag)
        return np.where(ct > 0, sm / np.maximum(ct, 1), np.nan)


    def _neu_profile(E, rows=None):
        """f(neural): mean of the map down each column (optionally only some rows)."""
        Er = E if rows is None else E[rows]
        ok = np.isfinite(Er)
        ct = ok.sum(axis=0)
        return np.where(ct > 0, np.where(ok, Er, 0).sum(axis=0) / np.maximum(ct, 1), np.nan)


    def _cv_axes(E):
        """-> (R2_neural, R2_lag) on HELD-OUT kinematic windows.
        Fold 1 builds both profiles from the first half of the kinematic windows and
        scores them on the second; fold 2 the reverse. _gap windows next to the split
        are left out of the building half. Both are scored on exactly the same cells,
        and data and predictions alike are ROW-CENTRED first (see _rowc), so neither
        is asked to predict how strong the map is at that kinematic time."""
        h = n_kin // 2
        folds = [(np.arange(0, h - _gap), np.arange(h, n_kin)),
                 (np.arange(h + _gap, n_kin), np.arange(0, h))]
        sse_n = sse_l = sst = 0.0
        for tr, te in folds:
            pn = _rowc(np.broadcast_to(_neu_profile(E, tr), (te.size, n_wins)).copy())
            pl = _rowc(_lag_profile(E, tr)[LAG_IDX[te]])
            y  = _rowc(E[te])
            ok = np.isfinite(y) & np.isfinite(pn) & np.isfinite(pl)
            if ok.sum() < 5: continue
            y = y[ok]
            sst   += np.sum((y - y.mean()) ** 2)
            sse_n += np.sum((y - pn[ok]) ** 2)
            sse_l += np.sum((y - pl[ok]) ** 2)
        if sst <= 0: return np.nan, np.nan
        return 1 - sse_n / sst, 1 - sse_l / sst


    R2N, R2L = {}, {}
    for pi in keep:
        R2N[pi], R2L[pi] = _cv_axes(EXC2[pi])
    DR2      = {pi: R2L[pi] - R2N[pi] for pi in keep}
    STRUCT   = {pi: bool(np.nanmax([R2N[pi], R2L[pi]]) > 0) for pi in keep}
    PROF_NEU = {pi: _neu_profile(EXC2[pi]) for pi in keep}      # bits, per neural window
    PROF_LAG = {pi: _lag_profile(EXC2[pi]) for pi in keep}      # bits, per lag bin

    print(f"\nlag bin {_lag_bin*1000:.0f} ms, {n_lag} lags "
          f"({lag_axis[0]*1000:.0f}..{lag_axis[-1]*1000:.0f} ms; every kinematic window "
          f"reaches {lag_axis[lag_full][0]*1000:.0f}..{lag_axis[lag_full][-1]*1000:.0f} ms)")
    print(f"CV: 2 blocked folds over {n_kin} kinematic windows, gap {_gap} "
          f"-> each fold builds on {n_kin//2 - _gap}; scored row-centred")


    def _med_iqr(v):
        q = np.percentile(v, [25, 75])
        return f"{np.median(v):+.3f} [{q[0]:+.2f},{q[1]:+.2f}]"


    for sub_name, sub_ok in (("ALL neurons  <- the headline", lambda pi: True),
                             ("only those whose shape repeats (secondary, power-limited)",
                              lambda pi: STRUCT[pi])):
        print(f"\n{sub_name}")
        print(f"{'region':6s} {'n':>5s} {'R2_neural':>10s} {'R2_lag':>8s} "
              f"{'dR2 median [IQR]':>24s} {'% diagonal':>11s} {'dR2 vs 0':>10s}")
        print("-" * 80)
        _d = {}
        for reg in ("CP", "CB"):
            idx = [pi for pi in keep if regions_run[pi] == reg and sub_ok(pi)]
            d = np.array([DR2[pi] for pi in idx]); m = np.isfinite(d); d = d[m]
            _d[reg] = d
            if d.size == 0:
                print(f"{reg:6s} {0:5d}"); continue
            rn = np.array([R2N[pi] for pi in idx])[m]; rl = np.array([R2L[pi] for pi in idx])[m]
            p0 = wilcoxon(d).pvalue if d.size >= 10 else np.nan
            print(f"{reg:6s} {d.size:5d} {np.median(rn):+10.3f} {np.median(rl):+8.3f} "
                  f"{_med_iqr(d):>24s} {100*np.mean(d > 0):10.0f}% {p0:10.1e}")
        if _d.get("CP", np.array([])).size and _d.get("CB", np.array([])).size:
            _u, _p = mannwhitneyu(_d["CP"], _d["CB"], alternative="two-sided")
            print(f"  CB - CP median dR2: {np.median(_d['CB']) - np.median(_d['CP']):+.3f}"
                  f"   (Mann-Whitney p = {_p:.1e})")
    print("\n  dR2 < 0: MI is set by WHEN the spikes are (vertical stripes)")
    print("  dR2 > 0: MI is set by the LAG between spikes and movement (diagonal)")
    print("  GRAIN OF SALT: a ridge whose slope is not 1 fits neither description, and")
    print("  a noisy map fits neither either, so both land near dR2 = 0 for different")
    print("  reasons. Read this as a population tendency, not a per-neuron verdict.")

    return dict(R2N=R2N, R2L=R2L, DR2=DR2, STRUCT=STRUCT, PROF_NEU=PROF_NEU, PROF_LAG=PROF_LAG, lag_axis=lag_axis, lag_full=lag_full, _lag_bin=_lag_bin, _gap=_gap)


def individual_landscapes(ctx, res, out_dir_name, N_EXAMPLES=40, EXAMPLE_MODE="spread"):
    """12's cell 5b: one figure per example neuron (its map, f(neural), h(lag)).
    Slow; nothing else depends on it."""
    EXC2 = ctx.EXC2
    keep = ctx.keep
    uids = ctx.uids
    run_idx = ctx.run_idx
    regions_run = ctx.regions_run
    win_centers = ctx.win_centers
    kin_centers = ctx.kin_centers
    CLR = ctx.CLR
    SAVE_FIGS = ctx.SAVE_FIGS
    SHOW_FIGS = ctx.SHOW_FIGS
    extent = ctx.extent
    ASPECT = ctx.ASPECT
    CMAP = ctx.CMAP
    zero_txt = ctx.zero_txt
    TIMES = ctx.TIMES
    R2N = res["R2N"]
    R2L = res["R2L"]
    DR2 = res["DR2"]
    STRUCT = res["STRUCT"]
    PROF_NEU = res["PROF_NEU"]
    PROF_LAG = res["PROF_LAG"]
    lag_axis = res["lag_axis"]
    lag_full = res["lag_full"]
    def _shade_partial(ax, scale=1000.0):
        return shade_partial(ax, lag_axis, lag_full, scale)

    # One figure per neuron: its map, f(neural) underneath, h(lag) below that, the
    # kinematic marginal on the right. These are for LOOKING at the neuron; the dR2
    # numbers are printed small at the bottom, deliberately not in the title, because
    # on a single map they are weak (a ridge whose slope is not 1, or simply a noisy
    # map, fits neither description).
    # Ordered by dR2 so a run of them sweeps from diagonal-like to neural-time-like.
    EX_DIR = config.fig_dir(out_dir_name, "individual_landscapes")


    _rng_ex = np.random.default_rng(0)
    examples = []
    for reg in ("CP", "CB"):
        idx = [pi for pi in keep if regions_run[pi] == reg]
        if not idx: continue
        order = sorted(idx, key=lambda pi: -DR2[pi] if np.isfinite(DR2[pi]) else np.inf)
        if N_EXAMPLES is None:
            sel = order                                     # every neuron
        else:
            k_ = min(N_EXAMPLES, len(order))
            if EXAMPLE_MODE == "top":      sel = order[:k_]
            elif EXAMPLE_MODE == "random": sel = list(_rng_ex.choice(order, k_, replace=False))
            else:                          sel = [order[int(round(t))]
                                                  for t in np.linspace(0, len(order)-1, k_)]
        for pi in sel:
            examples.append((reg, pi, order.index(pi) + 1, len(order)))
    print(f"drawing {len(examples)} individual landscapes"
          f"{' (ALL neurons)' if N_EXAMPLES is None else ''}")

    _t_ex = time.time()
    for reg, pi, rank, ntot in examples:
        M   = EXC2[pi]
        uid = uids[run_idx[pi]]
        col = CLR[reg]
        fig = plt.figure(figsize=(11, 11))
        gs  = fig.add_gridspec(3, 2, width_ratios=[3, 1], height_ratios=[3, 1.1, 1.1],
                               wspace=0.5, hspace=0.35)
        axm = fig.add_subplot(gs[0, 0])
        axr = fig.add_subplot(gs[0, 1], sharey=axm)     # kinematic marginal (viewing)
        axb = fig.add_subplot(gs[1, 0], sharex=axm)     # f(neural), aligned with the map
        axl = fig.add_subplot(gs[2, 0])                 # h(lag), its own time axis

        lim = np.nanpercentile(np.abs(M), 99) if np.isfinite(M).any() else 1.0
        im = axm.imshow(M, origin="lower", aspect=ASPECT, extent=extent, cmap=CMAP,
                        vmin=min(0, np.nanpercentile(M, 1)), vmax=lim)
        d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
        axm.plot(d, d, color="w", lw=1.0, ls="--", alpha=0.7)
        axm.axvline(0, color="w", lw=0.7, ls=":", alpha=0.6)
        axm.axhline(0, color="w", lw=0.7, ls=":", alpha=0.6)
        axm.set_ylabel(f"KINEMATIC window centre re {zero_txt} (s)")
        axm.set_xlabel(f"NEURAL window centre re {zero_txt} (s)")

        # kinematic marginal: the map averaged ACROSS each row (over neural windows).
        _ok = np.isfinite(M)
        g_kin = np.where(_ok.any(axis=1),
                         np.where(_ok, M, 0).sum(axis=1) / np.maximum(_ok.sum(axis=1), 1),
                         np.nan)
        axr.plot(g_kin, kin_centers, color=col, lw=1.8)
        axr.axvline(0, color="k", lw=0.6); axr.axhline(0, color="k", lw=0.6, ls=":")
        axr.set_xlabel("mean over neural\nwindows (bits)")
        axr.set_title("kinematic marginal\n(viewing only)", fontsize=8)
        plt.setp(axr.get_yticklabels(), visible=False)
        fig.colorbar(im, cax=axr.inset_axes([1.12, 0, 0.1, 1]), label="MI (bits)")

        axb.plot(win_centers, PROF_NEU[pi], color=col, lw=1.8)
        axb.axhline(0, color="k", lw=0.6); axb.axvline(0, color="k", lw=0.6, ls=":")
        axb.set(xlabel=f"NEURAL window centre re {zero_txt} (s)", ylabel="f(neural)\n(bits)",
                title="f(neural): averaged down each COLUMN")
        axb.title.set_fontsize(9)

        _shade_partial(axl, scale=1.0)
        axl.plot(lag_axis, PROF_LAG[pi], color=col, lw=1.8)
        axl.axhline(0, color="k", lw=0.6); axl.axvline(0, color="k", lw=0.6, ls=":")
        axl.set(xlabel="lag = neural - kinematic (s)     < 0: spikes BEFORE the movement",
                ylabel="h(lag)\n(bits)",
                title="h(lag): averaged along each DIAGONAL  "
                      "(grey: only some kinematic windows reach this lag)")
        axl.title.set_fontsize(9)
        for a_ in (axb, axl, axr): a_.spines[["top", "right"]].set_visible(False)

        fig.suptitle(f"{uid}  |  {reg}", fontsize=11)
        # secondary, small and grey: a single map does not settle this
        fig.text(0.5, 0.055,
                 f"held-out R2  neural-time {R2N[pi]:+.2f} / lag {R2L[pi]:+.2f}"
                 f"   ->  dR2 {DR2[pi]:+.2f}   (order {rank}/{ntot})"
                 f"{'' if STRUCT[pi] else '   neither repeats across the kinematic halves'}",
                 ha="center", va="top", fontsize=7.5, color="0.45")
        if SAVE_FIGS:
            _savefig_retry(fig, EX_DIR / f"{uid.replace('|','_')}.png", dpi=130,
                           bbox_inches="tight")
        if SHOW_FIGS: plt.show()
        plt.close(fig)

    TIMES["individual figures"] = time.time() - _t_ex
    print(f"{len(examples)} individual landscapes -> {EX_DIR} "
          f"({TIMES['individual figures']:.0f}s)")


def diagonal_vs_neural_figure(ctx, res, FIG_SUBSET="all"):
    """12's cell 6: the diagonal-vs-neural-time summary figure."""
    keep = ctx.keep
    regions_run = ctx.regions_run
    win_centers = ctx.win_centers
    n_kin = ctx.n_kin
    n_wins = ctx.n_wins
    CLR = ctx.CLR
    lab = ctx.lab
    fdir = ctx.fdir
    RUN_NAME = ctx.RUN_NAME
    SAVE_FIGS = ctx.SAVE_FIGS
    SHOW_FIGS = ctx.SHOW_FIGS
    zero_txt = ctx.zero_txt
    R2N = res["R2N"]
    R2L = res["R2L"]
    DR2 = res["DR2"]
    STRUCT = res["STRUCT"]
    PROF_NEU = res["PROF_NEU"]
    PROF_LAG = res["PROF_LAG"]
    lag_axis = res["lag_axis"]
    lag_full = res["lag_full"]
    _lag_bin = res["_lag_bin"]
    _gap = res["_gap"]
    def _shade_partial(ax, scale=1000.0):
        return shade_partial(ax, lag_axis, lag_full, scale)

    # Top    : the per-neuron comparison (dR2), its two ingredients, its distribution
    # Bottom : WHAT the two profiles look like, averaged per region
    _sub_ok  = (lambda pi: True) if FIG_SUBSET == "all" else (lambda pi: STRUCT[pi])
    _sub_lab = ("all neurons" if FIG_SUBSET == "all"
                else "only neurons whose shape repeats")
    _lf = plt.figure(figsize=(18, 10))
    gs  = _lf.add_gridspec(2, 3, hspace=0.45, wspace=0.28)

    # A: dR2 per neuron, as a proportion of each region
    ax = _lf.add_subplot(gs[0, 0])
    _bins = np.linspace(-1.5, 1.5, 41)
    for reg in ("CP", "CB"):
        d = np.array([DR2[pi] for pi in keep if regions_run[pi] == reg and _sub_ok(pi)])
        d = d[np.isfinite(d)]
        if not d.size: continue
        ax.hist(np.clip(d, _bins[0], _bins[-1]), bins=_bins,
                weights=np.full(d.size, 1 / d.size), color=CLR[reg], alpha=0.55,
                edgecolor="white", label=f"{reg} (n={d.size})")
        ax.axvline(np.median(d), color=CLR[reg], lw=1.8, ls="--",
                   label=f"{reg} median {np.median(d):+.2f}")
    ax.axvline(0, color="k", lw=1.0, ls=":")
    ax.set(xlabel="dR2 = R2_lag - R2_neural   (beyond ±1.5 piled in the end bins)",
           ylabel="proportion of neurons",
           title=f"diagonal vs neural time -- {_sub_lab}\n"
                 f"< 0 neural-time-like   |   > 0 diagonal-like")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=7)

    # B: the two ingredients, every neuron
    ax = _lf.add_subplot(gs[0, 1])
    _lo, _hi = -1.0, 1.0
    for reg in ("CP", "CB"):
        idx = [pi for pi in keep if regions_run[pi] == reg]
        x = np.clip([R2N[pi] for pi in idx], _lo, _hi)
        y = np.clip([R2L[pi] for pi in idx], _lo, _hi)
        ax.plot(x, y, "o", ms=3, alpha=0.3, color=CLR[reg], label=f"{reg} ({len(idx)})")
    ax.plot([_lo, _hi], [_lo, _hi], color="k", lw=1.0, ls="--")
    ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
    ax.set(xlim=(_lo, _hi), ylim=(_lo, _hi), aspect="equal",
           xlabel="held-out R2, neural-time  f(neural)",
           ylabel="held-out R2, lag  h(lag)",
           title="the two ingredients, every neuron\nabove the dashed line = diagonal wins")
    ax.legend(fontsize=7, loc="lower right")

    # C: cumulative dR2
    ax = _lf.add_subplot(gs[0, 2])
    for reg in ("CP", "CB"):
        d = np.sort([DR2[pi] for pi in keep if regions_run[pi] == reg and _sub_ok(pi)
                     and np.isfinite(DR2[pi])])
        if d.size:
            ax.plot(d, np.arange(1, d.size + 1) / d.size, color=CLR[reg], lw=1.8,
                    label=f"{reg} (n={d.size})")
    ax.axvline(0, color="k", lw=1.0, ls=":")
    ax.set(xlim=(-2, 1.5), xlabel="dR2", ylabel="cumulative proportion of neurons",
           title="same, cumulative\ncurve further RIGHT = more diagonal-like")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=7)

    # D: mean h(lag) per region
    ax = _lf.add_subplot(gs[1, 0])
    _shade_partial(ax)
    for reg in ("CP", "CB"):
        idx = [pi for pi in keep if regions_run[pi] == reg and _sub_ok(pi)]
        if not idx: continue
        Pm = np.stack([PROF_LAG[pi] for pi in idx])
        mn = np.nanmean(Pm, axis=0)
        se = np.nanstd(Pm, axis=0) / np.sqrt(np.maximum(np.isfinite(Pm).sum(0), 1))
        ax.plot(lag_axis * 1000, mn, color=CLR[reg], lw=2, label=f"{reg} (n={len(idx)})")
        ax.fill_between(lag_axis * 1000, mn - se, mn + se, color=CLR[reg], alpha=0.2)
    ax.axvline(0, color="k", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.5)
    ax.set(xlabel="lag = neural - kinematic (ms)     < 0: spikes BEFORE the movement",
           ylabel="mean excess MI ± SEM (bits)",
           title="h(lag): the map averaged along each DIAGONAL\n"
                 "(grey: only some kinematic windows reach this lag)")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

    # E: mean f(neu) per region
    ax = _lf.add_subplot(gs[1, 1])
    for reg in ("CP", "CB"):
        idx = [pi for pi in keep if regions_run[pi] == reg and _sub_ok(pi)]
        if not idx: continue
        Pm = np.stack([PROF_NEU[pi] for pi in idx])
        mn = np.nanmean(Pm, axis=0)
        se = np.nanstd(Pm, axis=0) / np.sqrt(np.maximum(np.isfinite(Pm).sum(0), 1))
        ax.plot(win_centers * 1000, mn, color=CLR[reg], lw=2, label=f"{reg} (n={len(idx)})")
        ax.fill_between(win_centers * 1000, mn - se, mn + se, color=CLR[reg], alpha=0.2)
    ax.axvline(0, color="k", lw=0.8, ls=":"); ax.axhline(0, color="k", lw=0.5)
    ax.set(xlabel=f"NEURAL window centre re {zero_txt} (ms)",
           ylabel="mean excess MI ± SEM (bits)",
           title="f(neural): the map averaged down each COLUMN\n"
                 "(whatever kinematic window is read)")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(fontsize=8)

    # F: how to read it, caveats included
    ax = _lf.add_subplot(gs[1, 2]); ax.axis("off")
    ax.text(0, 1, "\n".join([
        "Two 1-D summaries of each neuron's 2-D map:",
        "",
        "  f(neural)  MI depends on WHEN the spikes are,",
        "             whatever kinematics is read",
        "             -> vertical stripes",
        "  h(lag)     MI depends on neural - kinematic time",
        "             -> diagonal stripes",
        "",
        "Each is built from half of the kinematic windows and",
        "scored (R2) on the OTHER half, row by row, so how",
        "STRONG the map is at a kinematic time never enters,",
        "and h(lag)'s extra flexibility buys nothing on noise.",
        "",
        "  dR2 = R2_lag - R2_neural",
        "",
        "WITH A GRAIN OF SALT: a ridge whose slope is not 1",
        "fits neither summary, and so does a noisy map, for",
        "completely different reasons. This is a population",
        "tendency, not a per-neuron verdict.",
        "",
        f"grid {n_kin} kin x {n_wins} neural windows, lag bin "
        f"{_lag_bin*1000:.0f} ms, CV gap {_gap}",
    ]), va="top", ha="left", fontsize=9, family="monospace", transform=ax.transAxes)

    _lf.suptitle(f"Diagonal vs neural time  —  {lab}", fontsize=12)
    if SAVE_FIGS:
        _savefig_retry(_lf, fdir / f"{RUN_NAME}_diagonal_vs_neural.png", dpi=150,
                       bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(_lf)


    print(f"\nfigures -> {fdir}")


def per_mouse_maps(ctx, PER_MOUSE_MIN_N=5, PER_MOUSE_FILES=True, save_pdf_all_mice=False):
    """12's cell 7: the mean maps (CP, CB, CB - CP) per mouse and for all mice, as
    one grid (one colour scale) and one file per mouse (scaled to itself);
    optionally the ALL MICE file as PDF too."""
    EXC2 = ctx.EXC2
    keep = ctx.keep
    uids = ctx.uids
    run_idx = ctx.run_idx
    regions_run = ctx.regions_run
    lab = ctx.lab
    fdir = ctx.fdir
    RUN_NAME = ctx.RUN_NAME
    SAVE_FIGS = ctx.SAVE_FIGS
    SHOW_FIGS = ctx.SHOW_FIGS
    extent = ctx.extent
    ASPECT = ctx.ASPECT
    CMAP = ctx.CMAP
    zero_txt = ctx.zero_txt
    TIMES = ctx.TIMES

    # 7. cell --> the same three maps as cell 4 (CP, CB, CB - CP), one row per mouse,
    # plus an ALL MICE row at the bottom for reference. Saves one grid figure and,
    # if PER_MOUSE_FILES, one figure per mouse as well.

    MOUSE_OF = {pi: uids[run_idx[pi]].split("|")[1] for pi in keep}
    MICE_RUN = sorted({MOUSE_OF[pi] for pi in keep})
    ROWS = [(m, [pi for pi in keep if MOUSE_OF[pi] == m]) for m in MICE_RUN]
    ROWS.append(("ALL MICE", list(keep)))

    def _mean_maps(idx):
        """-> {'CP': map or None, 'CB': ..., 'n': {'CP': n, 'CB': n}}"""
        out, ns = {}, {}
        for reg in ("CP", "CB"):
            sel = [pi for pi in idx if regions_run[pi] == reg]
            ns[reg] = len(sel)
            out[reg] = (np.nanmean(EXC2[sel], axis=0)
                        if len(sel) >= PER_MOUSE_MIN_N else None)
        out["n"] = ns
        return out

    MAPS = {name: _mean_maps(idx) for name, idx in ROWS}
    print(f"\nper-mouse mean landscapes ({len(MICE_RUN)} mice)")
    for name, _ in ROWS:
        print(f"  {name:12s} CP n={MAPS[name]['n']['CP']:4d}   CB n={MAPS[name]['n']['CB']:4d}")

    # One colour scale for CP/CB across ALL panels, and one for every difference
    # panel, so the rows can be compared with each other rather than each being
    # stretched to its own range.
    _v = [m[r].ravel() for m in MAPS.values()
          for r in ("CP", "CB") if m[r] is not None]
    if not _v:
        raise ValueError(f"no mouse has >= PER_MOUSE_MIN_N ({PER_MOUSE_MIN_N}) neurons "
                         f"in either region -- lower it, or run more neurons")
    _v = np.concatenate(_v)
    VMAX = np.nanpercentile(_v, 99); VMIN = min(0.0, np.nanpercentile(_v, 1))
    _d = [m["CB"] - m["CP"] for m in MAPS.values()
          if m["CP"] is not None and m["CB"] is not None]
    DLIM = np.nanpercentile(np.abs(np.concatenate([d.ravel() for d in _d])), 99) if _d else 1.0


    def _panel(ax, M, title, diff=False, vmin=None, vmax=None, dlim=None):
        """vmin/vmax/dlim default to the figure-wide scales (VMIN/VMAX/DLIM) that the
        GRID uses, where every mouse must share one scale to be comparable. The
        single-mouse files pass their own, so each is scaled to itself."""
        if M is None:
            ax.text(0.5, 0.5, "too few neurons", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9, color="0.5")
            ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
            return None
        if diff:
            _d = DLIM if dlim is None else dlim
            im = ax.imshow(M, origin="lower", aspect=ASPECT, extent=extent, cmap="RdBu_r",
                           norm=TwoSlopeNorm(vcenter=0, vmin=-_d, vmax=_d))
            c = "k"
        else:
            im = ax.imshow(M, origin="lower", aspect=ASPECT, extent=extent, cmap=CMAP,
                           vmin=VMIN if vmin is None else vmin,
                           vmax=VMAX if vmax is None else vmax)
            c = "w"
        d = np.linspace(max(extent[0], extent[2]), min(extent[1], extent[3]), 2)
        ax.plot(d, d, color=c, lw=1.0, ls="--", alpha=0.7)
        ax.axvline(0, color=c, lw=0.7, ls=":", alpha=0.7)
        ax.axhline(0, color=c, lw=0.7, ls=":", alpha=0.7)
        ax.set_title(title, fontsize=9)
        return im


    _t_pm = time.time()
    _H = 3.6 * len(ROWS)
    fig, axes = plt.subplots(len(ROWS), 3, figsize=(15, _H), squeeze=False)
    for r, (name, _) in enumerate(ROWS):
        M = MAPS[name]
        D = (M["CB"] - M["CP"]) if (M["CP"] is not None and M["CB"] is not None) else None
        for c, (key, Mi, ttl, isd) in enumerate((
                ("CP", M["CP"], f"CP  (n={M['n']['CP']})", False),
                ("CB", M["CB"], f"CB  (n={M['n']['CB']})", False),
                ("D",  D,       "CB - CP", True))):
            ax = axes[r, c]
            im = _panel(ax, Mi, ttl, diff=isd)
            if im is not None:
                plt.colorbar(im, ax=ax, label=("Δ mean MI (bits)" if isd else "mean MI (bits)"),
                             fraction=0.046, pad=0.02)
            if c == 0:
                ax.set_ylabel(f"{name}\nKINEMATIC centre (s)", fontsize=9)
            if r == len(ROWS) - 1:
                ax.set_xlabel(f"NEURAL window centre re {zero_txt} (s)")
    fig.suptitle(f"Mean 2D landscape per mouse  —  {lab}\n"
                 f"dashed = simultaneity; ONE colour scale for every mouse, so the "
                 f"rows compare (the single-mouse files are scaled to themselves)",
                 fontsize=12)
    # leave room for the two-line suptitle: the figure grows with the number of
    # mice, so the margin has to be a fraction of the actual height, not a constant
    plt.tight_layout(rect=(0, 0, 1, 1 - 0.75 / _H))
    if SAVE_FIGS:
        _savefig_retry(fig, fdir / f"{RUN_NAME}_mean2D_per_mouse.png", dpi=150,
                       bbox_inches="tight")
    if SHOW_FIGS: plt.show()
    plt.close(fig)

    if PER_MOUSE_FILES:
        for name, _ in ROWS:
            M = MAPS[name]
            D = (M["CB"] - M["CP"]) if (M["CP"] is not None and M["CB"] is not None) else None
            # scale each single-mouse figure to ITS OWN maps: CP and CB share one
            # scale so the regions stay comparable within the animal, the difference
            # gets its own. (The grid above keeps one scale across all mice.)
            _mv = [m.ravel() for m in (M["CP"], M["CB"]) if m is not None]
            if _mv:
                _mv = np.concatenate(_mv)
                _vmax = np.nanpercentile(_mv, 99); _vmin = min(0.0, np.nanpercentile(_mv, 1))
            else:
                _vmin, _vmax = VMIN, VMAX
            _dlim = np.nanpercentile(np.abs(D), 99) if D is not None else DLIM
            f1, ax1 = plt.subplots(1, 3, figsize=(15, 4.2))
            for c, (Mi, ttl, isd) in enumerate((
                    (M["CP"], f"CP  (n={M['n']['CP']})", False),
                    (M["CB"], f"CB  (n={M['n']['CB']})", False),
                    (D,       "CB - CP", True))):
                im = _panel(ax1[c], Mi, ttl, diff=isd,
                            vmin=_vmin, vmax=_vmax, dlim=_dlim)
                if im is not None:
                    plt.colorbar(im, ax=ax1[c],
                                 label=("Δ mean MI (bits)" if isd else "mean MI (bits)"),
                                 fraction=0.046, pad=0.02)
                ax1[c].set_xlabel(f"NEURAL window centre re {zero_txt} (s)")
            ax1[0].set_ylabel(f"KINEMATIC window centre re {zero_txt} (s)")
            f1.suptitle(f"{name}  —  {lab}\ncolour scale fitted to this "
                        f"mouse (the per-mouse GRID shares one scale across mice)",
                        fontsize=11)
            plt.tight_layout()
            if SAVE_FIGS:
                _savefig_retry(f1, fdir / f"{RUN_NAME}_mean2D_{name.replace(' ', '_')}.png",
                               dpi=150, bbox_inches="tight")
                if name == "ALL MICE" and save_pdf_all_mice:
                    with plt.rc_context({"pdf.fonttype": 42}):
                        _savefig_retry(f1, fdir / f"{RUN_NAME}_mean2D_ALL_MICE.pdf",
                                       bbox_inches="tight")
            if SHOW_FIGS: plt.show()
            plt.close(f1)

    TIMES["per-mouse figures"] = time.time() - _t_pm
    print(f"per-mouse landscapes -> {fdir} ({TIMES['per-mouse figures']:.0f}s)")


def timing_summary(ctx, N_SHUF_2D, n_population):
    """where the time went (TIMES), and the extrapolation to the full population"""
    keep = ctx.keep
    n_kin = ctx.n_kin
    n_wins = ctx.n_wins
    TIMES = ctx.TIMES

    # ── where the time went ─────────────────────────────────────────────────────
    def _fmt(v):
        return f"{v:.1f}s" if v < 90 else f"{v/60:.1f} min"

    _tot = sum(TIMES.values())
    print(f"\n{'step':32s} {'time':>10s} {'% of total':>11s}")
    print("-" * 56)
    for _k, _v in sorted(TIMES.items(), key=lambda kv: -kv[1]):
        print(f"{_k:32s} {_fmt(_v):>10s} {100*_v/max(_tot, 1e-9):10.1f}%")
    print("-" * 56)
    print(f"{'TOTAL (timed steps)':32s} {_fmt(_tot):>10s}")
    if "2D landscapes (MI + nulls)" in TIMES:
        _per = TIMES["2D landscapes (MI + nulls)"] / max(len(keep), 1)
        print(f"\n{_per:.2f}s per neuron over {len(keep)} neurons x {n_kin*n_wins} cells "
              f"x {1+N_SHUF_2D} (obs+nulls)")
        print(f"scaling to the full population (~{n_population} neurons) would take "
              f"~{_per*n_population/60:.0f} min")


def save_results(path, ctx, res, MI2, NU2, NSD2, DO_NULL, **meta):
    """the run's results as one pickle (12's end of cell 7): the maps, their
    nulls, the per-neuron diagonal-vs-neural results, the axes, and `meta`
    (the parameters and the selection of the run, written as given)."""
    with open(path, "wb") as f:
        pickle.dump(dict(mi=MI2, null=NU2 if DO_NULL else None,
                         null_std=NSD2 if DO_NULL else None,
                         # diagonal vs neural time, per neuron (cross-validated)
                         R2_neural=res["R2N"], R2_lag=res["R2L"], dR2=res["DR2"],
                         structured=res["STRUCT"],
                         prof_neural=res["PROF_NEU"], prof_lag=res["PROF_LAG"],
                         lag_axis=res["lag_axis"], lag_full=res["lag_full"],
                         LAG_BIN=res["_lag_bin"], CV_GAP=res["_gap"],
                         kin_centers=ctx.kin_centers, win_centers=ctx.win_centers,
                         uids=[ctx.uids[i] for i in ctx.run_idx],
                         regions=ctx.regions_run, **meta), f, protocol=4)
    print(f"saved -> {path}")
