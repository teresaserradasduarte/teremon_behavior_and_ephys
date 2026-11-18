"""ct_config.py - paths and constants of the cell-type pipeline.

Steps (run in this order; each one only READS what the previous ones wrote):
  01_extract_features.py   per neuron: waveform + spike-train features      -> features.pkl
  02_CB_c4_review.py       cerebellum: the C4 output, incl. units C4 dropped -> cb_c4.csv
  03_CB_cs_ss_pairs.py     cerebellum: CS-triggered pause of nearby units    -> cb_pairs.pkl
  04_CB_relabel.py         cerebellum: C4 + pause rules -> automatic label   -> cb_auto.csv
  05_CP_classify.py        striatum: bombcell criteria + embedding/clusters  -> cp_auto.csv
  06_build_cell_type_table.py  old + automatic + manual -> the master table  -> cell_types.csv/.pkl
  07_review_gui.py         interactive check / relabel (writes overrides.csv only)
After relabelling in 07, run 06 again: it is the ONLY step that writes the master table.
Other projects read the result with   cell_types.load_cell_types()   (no pandas needed).

Nothing here writes into the old .mat files.
"""
import os
from pathlib import Path

ROOT_DIR = Path(os.environ.get(
    "MI_ROOT_DIR",
    r"D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD"))
GROUP   = "20230801_ChocolateGroup"
MICE    = ["1_CoteDor", "2_Lindt", "3_Toblerone", "4_Milka", "5_FerreroRocher"]
N_SESS  = 7

MAT_DIR   = ROOT_DIR / "ephys_and_behavior" / "mat_files" / GROUP
EPHYS_DIR = ROOT_DIR / "ephys_data" / GROUP
PATHS_MAT = EPHYS_DIR / "teremon_paths.mat"            # files_curation_CB / _BG (as s04 uses)
OLD_CB_LABELS_MAT = MAT_DIR / "eg_neurons_CB_cellType.mat"   # s04's labels; read only

# outputs: tables next to the .mat files, figures in out_files (as the MI project)
CT_DIR  = MAT_DIR / "cell_types_py"
FIG_DIR = ROOT_DIR / "ephys_and_behavior" / "out_files" / GROUP / "group" / "cell_types_py"

FEATURES_PKL  = CT_DIR / "features.pkl"
CB_C4_CSV     = CT_DIR / "cb_c4.csv"
CB_PAIRS_PKL  = CT_DIR / "cb_pairs.pkl"
CB_AUTO_CSV   = CT_DIR / "cb_auto.csv"
CP_AUTO_CSV   = CT_DIR / "cp_auto.csv"
CP_EMBED_PKL  = CT_DIR / "cp_embedding.pkl"
OVERRIDES_CSV = CT_DIR / "overrides.csv"       # manual labels (07); never deleted by code
MASTER_CSV    = CT_DIR / "cell_types.csv"      # THE table other scripts read
MASTER_PKL    = CT_DIR / "cell_types.pkl"      # same content, plain dicts/lists

# which eg_neurons_<KIND>.mat holds which region, and the brain_region kept from it
KINDS = {"CB": "CB", "BG": "CP"}               # file kind -> brain_region classified
FS = 30000.0                                   # Neuropixels AP sampling rate (Hz)

UNLAB = "unlabelled"
CB_TYPES = ["PkC_ss", "PkC_cs", "MLI", "GoC", "MFB", UNLAB]
CP_TYPES = ["MSN", "FSI", "TAN", "UIN", UNLAB]
TYPES = {"CB": CB_TYPES, "CP": CP_TYPES}

# C4's colours (npyx/c4/plots_functions.py) + striatum
COLORS = {"PkC_ss": "#1c78b5", "PkC_cs": "#000000", "MLI": "#e0559f", "MFB": "#d62529",
          "GrC": "#8f67a9", "GoC": "#38ae3e", UNLAB: "#969696",
          "MSN": "#d1495b", "FSI": "#00798c", "TAN": "#edae49", "UIN": "#66a182"}


def lp(p):
    """Windows long-path form. The C4 files sit deeper than 260 characters."""
    p = os.path.abspath(str(p))
    return p if os.name != "nt" or p.startswith("\\\\?\\") else "\\\\?\\" + p


def uid(region, mouse, sess, phy):
    """Same key as the MI project (mi_io): 'CB|1_CoteDor|R5|43'."""
    return f"{region}|{mouse}|{sess}|{int(phy)}"


def fig_dir(*parts):
    d = FIG_DIR.joinpath(*parts); d.mkdir(parents=True, exist_ok=True); return d


def out_path(p):
    Path(p).parent.mkdir(parents=True, exist_ok=True); return Path(p)
