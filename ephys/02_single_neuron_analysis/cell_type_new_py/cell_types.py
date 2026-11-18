"""cell_types.py - read the cell types from any project. Standard library only
(no pandas / numpy), so it works in every environment.

    import sys; sys.path.append(r"...\\ephys\\02_single_neuron_analysis\\cell_type_new_py")
    import cell_types
    CT = cell_types.load_cell_types("CB")          # {uid: "PkC_ss", ...}
    CT = cell_types.load_cell_types("CB", column="old_cell_type")   # s04's labels

uid is the MI project's key: "CB|1_CoteDor|R5|43"  (region|mouse|session|phyID).

Columns of the table (cell_types.csv, written by 06_build_cell_type_table.py):
    cell_type        THE label to use: the manual one if there is one, else the automatic one
    old_cell_type    s04_CB_cellType.m's label (CB only), unchanged
    c4_pred          C4's best guess (CB), also for neurons C4 did not report
    auto_cell_type   the automatic label of this pipeline (04 for CB, 05 for CP)
    auto_source      c4 | c4_relaxed | cs_pause | cs_pair | criteria
    auto_conf        high | low
    flags            ';'-separated things to check
    manual_cell_type / manual_note / manual_time     from the review GUI (07)
    partner          uid of the paired CS / SS unit (CB)
"""
import csv, os
from pathlib import Path

_ROOT = Path(os.environ.get(
    "MI_ROOT_DIR", r"D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD"))
TABLE = Path(os.environ.get("CELL_TYPES_CSV", str(
    _ROOT / "ephys_and_behavior" / "mat_files" / "20230801_ChocolateGroup"
    / "cell_types_py" / "cell_types.csv")))
UNLAB = "unlabelled"


def load_table(path=None):
    """-> list of dicts, one per neuron (all values are strings)."""
    p = Path(path or TABLE)
    if not p.is_file():
        raise FileNotFoundError(f"{p} not found -- run 06_build_cell_type_table.py")
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_cell_types(region=None, column="cell_type", min_conf=None, path=None):
    """-> {uid: label}.
    region    "CB" / "CP" / None (both)
    column    "cell_type" (final) | "auto_cell_type" | "old_cell_type" | "c4_pred"
    min_conf  "high" -> automatic labels of low confidence become "unlabelled"
              (manual labels are always kept)"""
    out = {}
    for r in load_table(path):
        if region and r["region"] != region: continue
        lab = r.get(column, "") or UNLAB
        if (min_conf == "high" and column == "cell_type" and not r.get("manual_cell_type")
                and r.get("auto_conf") != "high"):
            lab = UNLAB
        out[r["uid"]] = lab
    return out
