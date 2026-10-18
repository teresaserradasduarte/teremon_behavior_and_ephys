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
BEHAV_GROUP = "20230511_ChocolateGroup"   # raw behaviour tree (water tracking) -- a
                                          # DIFFERENT group folder from GROUP above,
                                          # same mouse/session names underneath.
ANIMALS = ["CoteDor", "Lindt", "Toblerone", "Milka", "FerreroRocher"]

# region code -> struct variable name inside the .mat.
# NOTE: we need the "_sts" structs (they carry st_reach spike times). This is a
# DIFFERENT file/variant from the MATLAB "all_<REG>_neurons" (spikes removed there).
REGION_STRUCT = {r: f"all_{r}_sts" for r in ["CP", "CB", "CTX", "DCN", "BS", "MB", "P"]}
REGIONS = ["CP", "CB"]                      # regions to load by default

# ---- derived paths (mirror the MATLAB layout) ---------------------------
MAT_DIR = ROOT_DIR / "ephys_and_behavior" / "mat_files" / GROUP
# raw behaviour tree: <BEHAV_MAT_DIR>/<mouse>/<sess>/session_reaching_data_paw.mat
BEHAV_MAT_DIR = (ROOT_DIR / "behavior_data" / "analyzed_data" / "mat_files"
                 / BEHAV_GROUP / SETUP)
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

# px -> mm for (x, y, z). Verified exactly against bhv.reaches_inVec_mm /
# bhv.reaches_inVec_px, which agree to 0 std on every session.
PX2MM = (0.1, 0.07, 0.1)

# ── how the reaches were CUT in MATLAB (not an analysis choice) ──────────────
# step1_find_reaches...m keeps reach_goBack_frames=150 frames before the
# max-reach frame (max_reach=151, MATLAB 1-based) plus 50 after, so every reach
# record is 201 frames with the REACH ENDPOINT at 0-based index 150 -- and
# st_reach t=0 is that same frame. Anything indexed off these already-cut
# reaches inherits that zero; it is a property of the data, NOT a knob.
# Verified against reaches.tm_w in the behaviour mats, which runs -375..+125 ms
# with zero at this index.
# (A window placed AROUND the endpoint, like 09's FRAMES=(147,167), is a
# different thing: that one IS an analysis choice.)
REACH_N_FRAMES       = 201
REACH_ENDPOINT_FRAME = 150

# Raw behaviour tree: holds the LOAD CELL stream (behavior.inputs.paws_all,
# (n, 2) = absolute harp timestamp + force in mV/V at ~1 kHz) and the trial
# initiation times in SESSION SECONDS (behavior.init.timeof). The load cell is
# the signal that defines a push/pull, so it is the ground truth for that event.
BEHAV_RAW_DIR = ROOT_DIR / "behavior_data" / "raw_data" / BEHAV_GROUP / SETUP

def behav_raw_session_mat(mouse, sess):
    """Raw behavior_session.mat (load cells + initiation times) for one session."""
    return BEHAV_RAW_DIR / mouse / sess / "behavior_session.mat"

def behav_session_mat(mouse, sess):
    """Raw behaviour .mat (v7.3) holding trials.water for one session."""
    return BEHAV_MAT_DIR / mouse / sess / "session_reaching_data_paw.mat"

def bhv_fundamentals_mat(mouse, sess):
    """behavior_fundamentals.mat holding bhv.reach_trials_inVec etc."""
    return MAT_DIR / mouse / sess / "behavior_fundamentals.mat"

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
