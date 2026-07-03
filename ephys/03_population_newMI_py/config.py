"""config.py - paths and project constants. EDIT THE ROOTS for your machine.
DATA and OUTPUTS live under a data root, SEPARATE from this code folder
(mirrors the MATLAB rootdir / group / out_files layout). Per-machine overrides
via environment variables MI_ROOT_DIR and MI_DATA_MAT (no need to edit the file).
"""
#%%
import os
from pathlib import Path

#%% ---- EDIT THESE ---------------------------------------------------------
ROOT_DIR = Path(os.environ.get(
    "MI_ROOT_DIR",
    r"D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD"))
GROUP   = "20230801_ChocolateGroup"
SETUP   = "headfixed_dynamicTarget"
ANIMALS = ["CoteDor", "Lindt", "Toblerone", "Milka", "FerreroRocher"]

# region code -> struct variable name inside the .mat.
# NOTE: we need the "_sts" structs (they carry st_reach spike times). This is a
# DIFFERENT file/variant from the MATLAB "all_<REG>_neurons" (spikes removed there).
REGION_STRUCT = {r: f"all_{r}_sts" for r in ["CP", "CB", "CTX", "DCN", "BS", "MB", "P"]}
REGIONS = ["CP", "CB"]                      # regions to load by default

# ---- derived paths (mirror the MATLAB layout) ---------------------------
MAT_DIR = ROOT_DIR / "ephys_and_behavior" / "mat_files" / GROUP
OUT_DIR = ROOT_DIR / "ephys_and_behavior" / "out_files" / GROUP / "group" / "MutualInformation_NEW"

# .mat holding all_<REG>_sts (spike times + behaviour). Override with env MI_DATA_MAT.
DATA_MAT = Path(os.environ.get("MI_DATA_MAT", str(MAT_DIR / "st_eg_neurons.mat")))

# VAR_DIR: intermediate variables / cache (next to the .mat, not in out_files)
VAR_DIR = MAT_DIR / "MutualInformation_py"

# CACHE_DIR: parsed NEURONS pickle (inside VAR_DIR)
CACHE_DIR = VAR_DIR / "python_cache"

def struct_names(regions=None):
    return [REGION_STRUCT[r] for r in (regions or REGIONS)]

def cache_path(tag=None):
    tag = tag or "-".join(REGIONS)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f"parsed_{tag}.pkl"

def fig_dir(*parts):
    """fig_dir("ksg", "landscapes_individual") → OUT_DIR/ksg/landscapes_individual/"""
    d = OUT_DIR.joinpath(*parts)
    d.mkdir(parents=True, exist_ok=True)
    return d

def var_path(filename):
    """var_path("landscape_ksg.pkl") → VAR_DIR/landscape_ksg.pkl"""
    VAR_DIR.mkdir(parents=True, exist_ok=True)
    return VAR_DIR / filename

# %%
