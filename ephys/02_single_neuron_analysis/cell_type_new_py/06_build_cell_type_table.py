# %% ── 0. parameters ─────────────────────────────────────────────────────────
# STEP 06 -- the master table: one row per CB / CP neuron, with every label it
# ever had. This is the ONLY step that writes cell_types.csv / cell_types.pkl,
# the files other projects read (cell_types.load_cell_types()).
#
#   cell_type = manual_cell_type   if you set one in the GUI (07)
#             = auto_cell_type     otherwise (04 for CB, 05 for CP)
# The old s04 label is carried along unchanged (old_cell_type); no .mat is touched.
# Run this again after every GUI session, or after changing a threshold in 04/05:
# the manual labels live in overrides.csv and are re-applied each time.
# A manual label whose neuron is no longer in the table is reported, not dropped.
import collections, time
import matplotlib.pyplot as plt
import numpy as np
import ct_config as C, ct_io as io

SHOW_FIGS = True
SAVE_FIGS = True
SAVE_VARS = True       # writes cell_types.csv and cell_types.pkl

COLS = ["uid", "region", "mouse", "sess", "phyID", "cell_type", "label_origin",
        "old_cell_type", "c4_status", "c4_pred", "auto_cell_type", "auto_source", "auto_conf",
        "flags", "partner", "manual_cell_type", "manual_note", "manual_time"]


# %% ── 1. merge ──────────────────────────────────────────────────────────────
FEAT = io.load_features()
AUTO = {r["uid"]: r for r in io.read_csv(C.CB_AUTO_CSV) + io.read_csv(C.CP_AUTO_CSV)}
for name, p in (("04", C.CB_AUTO_CSV), ("05", C.CP_AUTO_CSV)):
    if not p.is_file(): print(f"!! {p.name} not found (step {name} not run): those neurons stay unlabelled")
MAN = {}
for r in io.read_csv(C.OVERRIDES_CSV):                    # last line per uid wins
    MAN[r["uid"]] = r
MAN = {u: r for u, r in MAN.items() if r.get("label", "")}   # empty label = override cleared

rows = []
for f in FEAT["rows"]:
    a, m = AUTO.get(f["uid"], {}), MAN.get(f["uid"])
    auto = a.get("auto_cell_type", "") or C.UNLAB
    r = dict(uid=f["uid"], region=f["region"], mouse=f["mouse"], sess=f["sess"], phyID=f["phyID"],
             old_cell_type=a.get("old_cell_type", ""), c4_status=a.get("c4_status", ""),
             c4_pred=a.get("c4_pred", ""), auto_cell_type=auto, auto_source=a.get("auto_source", ""),
             auto_conf=a.get("auto_conf", ""), flags=a.get("flags", ""), partner=a.get("partner", ""),
             manual_cell_type=m["label"] if m else "", manual_note=m.get("note", "") if m else "",
             manual_time=m.get("time", "") if m else "")
    r["cell_type"] = r["manual_cell_type"] or auto
    r["label_origin"] = "manual" if m else ("auto" if auto != C.UNLAB else "")
    rows.append(r)

known = {r["uid"] for r in rows}
orphan = [u for u in MAN if u not in known]
if orphan: print(f"!! {len(orphan)} manual labels have no neuron in features.pkl, e.g. {orphan[:3]}")
for reg in ("CB", "CP"):
    sel = [r for r in rows if r["region"] == reg]
    print(f"\n{reg}: {len(sel)} neurons, {sum(r['label_origin'] == 'manual' for r in sel)} manual")
    for col in ("old_cell_type", "auto_cell_type", "cell_type"):
        c = collections.Counter(r[col] or "-" for r in sel)
        print(f"  {col:15s}", "  ".join(f"{k} {v}" for k, v in c.most_common()))

if SAVE_VARS:
    io.write_csv(C.MASTER_CSV, rows, COLS)
    io.save_pkl(C.MASTER_PKL, dict(
        columns=COLS, rows=rows, cell_type={r["uid"]: r["cell_type"] for r in rows},
        made=time.strftime("%Y-%m-%d %H:%M"),
        note="cell_type = manual_cell_type if set, else auto_cell_type; old_cell_type = s04"))
    print(f"\nsaved -> {C.MASTER_CSV}\n         {C.MASTER_PKL}")


# %% ── 2. figure: the final labels ───────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(15, 4.6))
for ax, reg in zip(axes, ("CB", "CP")):
    sel = [r for r in rows if r["region"] == reg]
    types = C.TYPES[reg]; x = np.arange(len(types)); bottom = np.zeros(len(types))
    for origin, alpha, hatch in (("auto", 1.0, None), ("manual", 1.0, "//"), ("", 0.4, None)):
        h = np.array([sum(r["cell_type"] == t and r["label_origin"] == origin for r in sel) for t in types], float)
        ax.bar(x, h, bottom=bottom, color=[C.COLORS[t] for t in types], alpha=alpha, hatch=hatch,
               edgecolor="w", label={"auto": "automatic", "manual": "manual (hatched)", "": None}[origin])
        bottom += h
    for xi, b in zip(x, bottom): ax.text(xi, b, int(b), ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x, types); ax.set(ylabel="neurons", title=f"{reg}: final cell_type (n={len(sel)})")
    ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
if SAVE_FIGS: fig.savefig(C.fig_dir("06_table") / "final_cell_types.png", dpi=150)
if SHOW_FIGS: plt.show()
plt.close(fig)
# %%
