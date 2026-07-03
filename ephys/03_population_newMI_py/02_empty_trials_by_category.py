# %% setup (loads cache; skips if NEURONS already in the interactive kernel)
import numpy as np, pandas as pd, matplotlib.pyplot as plt, mi_io, config
if "NEURONS" not in globals():
    NEURONS, REACH_PX = mi_io.load(config.cache_path())
CLR = {"CP": np.array([106, 176, 134]) / 255, "CB": np.array([154, 123, 193]) / 255}

win, bin_size = (-0.05, 0.05), 0.0025          # <<< window that DEFINES "empty"; vary me
SAVE_FIGS = False
CATS = {"hit": "idx_reach_hit", "push/pull": "idx_reach_PP",
        "L/C/R": "idx_reach_LCR", "reach cat": "idx_reach_cat"}

# %% build long table — vectorized (one array per neuron, no Python trial loop)
chunks = []
for n in NEURONS:
    c = mi_io.counts(n, win, bin_size)           # (n_trials,) int array
    chunk = {"region": n["region"], "uid": n["uid"], "empty": c == 0}
    for name, lf in CATS.items():
        chunk[name] = n["labels"][lf]
    chunks.append(pd.DataFrame(chunk))
df = pd.concat(chunks, ignore_index=True)

print("overall empty trials by region:")
g = df.groupby("region")["empty"]
print(pd.DataFrame({
    "n_empty": g.sum().astype(int),
    "n_total": g.count().astype(int),
    "pct":     (g.mean() * 100).round(1),
}))

# %% absolute + % empty per category level, per region (pooled over trials)
for name in CATS:
    sub = df.dropna(subset=[name])
    g = sub.groupby(["region", name])["empty"]
    t = pd.DataFrame({
        "n_empty": g.sum().astype(int),
        "n_total": g.count().astype(int),
        "pct":     (g.mean() * 100).round(1),
    }).reset_index()
    print(f"\nempty by {name} (pooled):")
    print(t.pivot(index=name, columns="region").to_string())

# %% within-neuron mean % empty (so high-trial neurons don't dominate)
def per_neuron_rate(name):
    g = df.dropna(subset=[name]).groupby(["region", "uid", name])["empty"].mean()
    return (g.groupby(["region", name]).mean() * 100).round(1).unstack(0)

def per_neuron_counts(name):
    sub = df.dropna(subset=[name])
    pooled = sub.groupby(["region", name])["empty"].agg(
        n_empty="sum", n_total="count"
    ).astype(int)
    within_pct = per_neuron_rate(name).stack().rename("within_neu_pct")
    return pooled.join(within_pct).reset_index()

for name in CATS:
    print(f"\nwithin-neuron mean % empty by {name}:")
    print(per_neuron_counts(name).to_string(index=False))

# %% plot: grouped bars of within-neuron % empty per category
fig, axes = plt.subplots(1, len(CATS), figsize=(4 * len(CATS), 3.4), sharey=True)
axes = np.atleast_1d(axes)
for ax, name in zip(axes, CATS):
    t = per_neuron_rate(name)
    t.plot.bar(ax=ax, color=[CLR.get(c, "gray") for c in t.columns], legend=(ax is axes[0]))
    ax.set(title=name, xlabel="", ylabel="% empty (within-neuron)")
plt.tight_layout()
if SAVE_FIGS:
    fig.savefig(config.fig_dir("compare") / "empty_by_category.png", dpi=150, bbox_inches="tight")
plt.show()
