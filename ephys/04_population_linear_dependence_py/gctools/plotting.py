"""plotting.py - figure helpers in the MI project's style (top/right spines off,
tight_layout, dpi 150, Dropbox-safe save, show/close). Region colours from
config.CLR (CP green, CB purple, SIM grey).

Result tables (one row per neuron) carry, per direction d in {"xy", "yx"}:
  gc_d, F_d, pF_d, pS_d (shuffle p), z_d (null-standardised GC), sig_d (BH-FDR)
"xy" = kinematics -> neuron (x->y), "yx" = neuron -> kinematics (y->x).
Directions are compared only through shuffle p / null-z (rule 5).
"""
import numpy as np
import matplotlib.pyplot as plt
from . import config
from .cache import savefig_retry
from .stats import binom_ci

DIR_LABEL = {"xy": "kin → neuron (x→y)", "yx": "neuron → kin (y→x)"}


def finish(fig, save_dir=None, name=None, save=True, show=True):
    for ax in fig.axes:
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    fig.tight_layout()
    if save and save_dir is not None:
        savefig_retry(fig, save_dir / f"{name}.png")
    if show:
        plt.show()
    plt.close(fig)


def _colors(df, highlight=None):
    c = [config.CLR.get(r, (0.4, 0.4, 0.4)) for r in df["region"]]
    if highlight is not None:
        c = [(0.85, 0.2, 0.2) if h else ci for ci, h in zip(c, highlight)]
    return c


# direction classes (from BH-FDR significance per direction)
DIR_CLASSES   = ("kin→neuron", "neuron→kin", "both", "none")
DIR_CLASS_CLR = {"kin→neuron": "tab:blue", "neuron→kin": "tab:red",
                 "both": "0.15", "none": "0.8"}


def direction_class(sig_xy, sig_yx):
    """Per-neuron class from FDR significance of x->y (kin->neuron) and y->x."""
    sig_xy, sig_yx = np.asarray(sig_xy, bool), np.asarray(sig_yx, bool)
    return np.where(sig_xy & sig_yx, "both",
                    np.where(sig_xy, "kin→neuron", np.where(sig_yx, "neuron→kin", "none")))


def direction_map(ax, df, n_shuffle, title="", highlight=None):
    """-log10 shuffle p: x->y (horizontal) vs y->x (vertical), one point per neuron,
    coloured by direction class. Solid lines = the FDR threshold actually reached in
    this set (largest significant p); dashed = p 0.05. Quadrants are labelled."""
    floor = 1.0 / (n_shuffle + 1)
    lx = -np.log10(np.clip(df["pS_xy"].to_numpy(float), floor, 1))
    ly = -np.log10(np.clip(df["pS_yx"].to_numpy(float), floor, 1))
    cls = df["direction"].to_numpy()
    for c in DIR_CLASSES:
        m = cls == c
        if m.any():
            ax.scatter(lx[m], ly[m], s=16, color=DIR_CLASS_CLR[c], alpha=0.85,
                       label=f"{c} ({m.sum()})", edgecolors="none")
    if highlight is not None and np.any(highlight):
        h = np.asarray(highlight, bool)
        ax.scatter(lx[h], ly[h], s=60, facecolors="none", edgecolors="k", lw=1,
                   label="sim driver")
    m = -np.log10(floor) * 1.05
    for d, f in (("xy", ax.axvline), ("yx", ax.axhline)):
        f(-np.log10(0.05), color="0.6", ls="--", lw=0.6)
        sig = df[f"sig_{d}"].to_numpy(bool)
        if sig.any():
            f(-np.log10(df.loc[sig, f"pS_{d}"].max()), color="k", lw=0.8)
    kw = dict(fontsize=7, color="0.35", ha="center", va="center")
    ax.text(0.8 * m, 0.03 * m, "kin→neuron only", **kw)
    ax.text(0.2 * m, 0.97 * m, "neuron→kin only", **kw)
    ax.text(0.8 * m, 0.97 * m, "both", **kw)
    ax.text(0.2 * m, 0.03 * m, "none", **kw)
    ax.set(xlim=(0, m), ylim=(0, m), title=title,
           xlabel="-log10 p_shuffle   kin → neuron (x→y)",
           ylabel="-log10 p_shuffle   neuron → kin (y→x)")
    ax.legend(fontsize=6, loc="center right", frameon=False)


def direction_classes(ax, df, regions, title=""):
    """Stacked bars: % of neurons per direction class, one bar per region."""
    for i, reg in enumerate(regions):
        g = df[df["region"] == reg]
        n, bottom = len(g), 0.0
        if n == 0:
            continue
        for c in DIR_CLASSES:
            f = float(np.mean(g["direction"] == c))
            ax.bar(i, f, bottom=bottom, color=DIR_CLASS_CLR[c], width=0.6,
                   label=c if i == 0 else None)
            if f >= 0.04:
                ax.text(i, bottom + f / 2, f"{100*f:.0f}%", ha="center", va="center",
                        fontsize=7, color="w" if c in ("kin→neuron", "neuron→kin", "both") else "k")
            bottom += f
        ax.text(i, 1.02, f"n={n}", ha="center", fontsize=7)
    ax.set_xticks(range(len(regions)), regions)
    ax.set(ylim=(0, 1.08), ylabel="fraction of neurons", title=title)
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False)


def lag_hist(ax, df, d, key, regions, bins):
    """Histogram of the preview lag (key "peak" | "centroid") for neurons FDR-significant
    in direction d, per region, with medians."""
    for reg in regions:
        g = df[(df["region"] == reg) & df[f"sig_{d}"]]
        v = g[f"lag_{key}_{d}"].dropna().to_numpy(float)
        if v.size == 0:
            continue
        c = config.CLR.get(reg, (0.4, 0.4, 0.4))
        ax.hist(v, bins=bins, color=c, alpha=0.55, label=f"{reg} n={v.size}, median {np.median(v):.0f} ms")
        ax.axvline(np.median(v), color=c, lw=1.5, ls="--")
    ax.set(xlabel=f"{key} lag (ms)  {DIR_LABEL[d]}", ylabel="# significant neurons")
    ax.legend(fontsize=6, frameon=False)


def z_compare(ax, a, b, df, la, lb, title="", highlight=None):
    """Null-standardised GC of the same neurons under two settings (e.g. raw vs matched)."""
    ax.scatter(a, b, s=14, c=_colors(df, highlight), alpha=0.8)
    lo = np.nanmin([np.nanmin(a), np.nanmin(b), -3]); hi = np.nanmax([np.nanmax(a), np.nanmax(b), 3])
    ax.plot([lo, hi], [lo, hi], color="0.6", lw=0.6)
    ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
    ax.set(xlabel=la, ylabel=lb, title=title)


def detection_bars(ax, rows, title=""):
    """rows: list of (label, k, n, expected_true). Bars = detection rate with 95% CI;
    dashed line at 0.05 (the nominal rate for null cases)."""
    for i, (lab, k, n, exp) in enumerate(rows):
        if n == 0:
            continue
        r = k / n; lo, hi = binom_ci(k, n)
        ax.bar(i, r, color="tab:red" if exp else "0.6", alpha=0.8)
        ax.errorbar(i, r, yerr=[[r - lo], [hi - r]], color="k", lw=1, capsize=3)
        ax.text(i, hi + 0.02, f"{k}/{n}", ha="center", fontsize=6)
    ax.axhline(0.05, color="k", ls="--", lw=0.7)
    ax.set_xticks(range(len(rows)), [r[0] for r in rows], rotation=60, ha="right", fontsize=7)
    ax.set(ylim=(0, 1.08), ylabel="fraction p_shuffle < 0.05", title=title)
