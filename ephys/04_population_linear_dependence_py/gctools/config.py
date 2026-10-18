"""config.py - paths and shared analysis constants for the linear-dependence
(Granger-style) project. Mirrors the MI project's config.py.

DATA and OUTPUTS live under the MI data root, SEPARATE from this code folder:
  VAR_DIR   = <MAT_DIR>/LinearDependence_py            pkls   (next to MutualInformation_py)
  CACHE_DIR = <VAR_DIR>/python_cache                   parameter-hashed caches
  OUT_DIR   = <ROOT>/ephys_and_behavior/out_files/<GROUP>/group/LinearDependence   figures, tables
Roots come from the MI config (so its env overrides MI_ROOT_DIR / MI_DATA_MAT apply
here too). The first folder level under VAR_DIR and OUT_DIR is the LAG TAG (e.g.
"lag150_win-375to125") so short- and long-lag runs coexist and never overwrite.

The MI modules are loaded BY FILE PATH under the aliases "mi_config" / "mi_io"
(the folder name 03_... is not importable, and a bare "config" would clash).
Nothing in the MI project is modified.
"""
import sys
import importlib.util
from pathlib import Path

# ---- MI project (read-only reuse) -----------------------------------------
PROJECT_DIR = Path(__file__).resolve().parents[1]
MI_DIR      = PROJECT_DIR.parent / "03_population_newMI_py"

def _load_by_path(alias, path):
    if alias in sys.modules:
        return sys.modules[alias]
    spec = importlib.util.spec_from_file_location(alias, str(path))
    mod  = importlib.util.module_from_spec(spec)
    sys.modules[alias] = mod
    spec.loader.exec_module(mod)
    return mod

mi_config = _load_by_path("mi_config", MI_DIR / "config.py")
mi_io     = _load_by_path("mi_io",     MI_DIR / "mi_io.py")

# ---- derived paths (mirror the MI layout) ---------------------------------
ROOT_DIR  = mi_config.ROOT_DIR
GROUP     = mi_config.GROUP
MAT_DIR   = mi_config.MAT_DIR
VAR_DIR   = MAT_DIR / "LinearDependence_py"
CACHE_DIR = VAR_DIR / "python_cache"
OUT_DIR   = ROOT_DIR / "ephys_and_behavior" / "out_files" / GROUP / "group" / "LinearDependence"

# ---- recording / alignment constants --------------------------------------
# st_reach t=0 = the MAX-REACH frame (NOT reach onset): step1_find_reaches...m sets
# reach_goBack_frames=150, max_reach=151 (MATLAB 1-idx), and st_reach = st - timestamp
# of that frame. reach_px holds 201 frames = 150 before + max + 50 after.
FPS         = 400.0
FRAME_MS    = 1000.0 / FPS          # 2.5 ms
ALIGN_FRAME = 150                   # python 0-idx of the max-reach frame
N_FRAMES    = 201                   # tracked frames per reach -> t in [-375, +125] ms
SPIKE_SPAN_S = (-3.0, 3.0)          # st_reach covers ~[-3.09, 3.09] s; use the safe inner span
PX2MM       = tuple(mi_config.PX2MM)   # (x, y med-lat, z)
KIN_COORD   = {"x": 0, "y": 1, "z": 2}
KIN_NAME    = {"x": "x", "y": "med-lat", "z": "z"}   # axis labels / run info

# ---- shared analysis defaults (set ONCE; scripts may override in cell 0) ---
REGIONS            = ["CP", "CB"]
ALIGN_EVENT        = "max_reach"
TRIAL_WINDOW_MS    = (-375.0, 125.0)    # the span with tracked kinematics
CROSS_LAG_RANGE_MS = (0.0, 150.0)       # lags (lo, hi] ms; min lag is always 1 bin.
                                        # FIRST PASS ONLY -- longer lags to be run later.
KIN_VAR            = "y"                # kinematic variable: "x" | "y" (med-lat) | "z" (position, mm);
                                        # scripts can override (KIN_VAR in cell 0)
KIN_SOURCE         = "short"            # "short" = the MI reach_px (201 frames, -375..+125 ms);
                                        # "long" = REACH_LONG (scripts/00b; -3..+3 s, same reaches):
                                        # real kinematic HISTORY before -375 ms. Window unchanged.
SHUFFLE_WITHIN     = ("LCR",)           # trial-shuffle strata; None = plain permutation
REACH_CAT          = (1,)
FIRST_PASS         = [("5_FerreroRocher", "R1"), ("1_CoteDor", "R6")]

CLR = {"CP": (106/255, 176/255, 134/255), "CB": (154/255, 123/255, 193/255),
       "SIM": (0.45, 0.45, 0.45)}


def lag_tag(cross_lag_range_ms=None, trial_window_ms=None):
    """lag_tag((0,150), (-375,125)) -> "lag150_win-375to125". First folder level of
    every output, and part of result filenames, so runs with different lag ranges /
    windows live side by side."""
    lo, hi = cross_lag_range_ms or CROSS_LAG_RANGE_MS
    w0, w1 = trial_window_ms or TRIAL_WINDOW_MS
    lag = f"lag{lo:g}-{hi:g}" if lo > 0 else f"lag{hi:g}"
    return f"{lag}_win{w0:g}to{w1:g}"

LAG_TAG = lag_tag()


MAX_PATH = 259      # Windows limit when LongPathsEnabled is off (it is on this machine);
                    # Explorer / Dropbox may not handle longer paths, so we warn.

def long_path(path):
    """Path string usable by Python beyond MAX_PATH (\\\\?\\ prefix on Windows)."""
    s = str(Path(path).resolve())
    if sys.platform == "win32" and len(s) >= 240 and not s.startswith("\\\\?\\"):
        s = "\\\\?\\" + s
    return s

def warn_long(path):
    if len(str(path)) > MAX_PATH:
        print(f"WARNING: path is {len(str(path))} characters (> {MAX_PATH}); Python can write it, "
              f"but Windows Explorer / Dropbox may not handle it -- shorten the run name:\n  {path}")

def var_path(*parts):
    """var_path(LAG_TAG, "fit_x.pkl") -> VAR_DIR/LAG_TAG/fit_x.pkl (folders created)."""
    p = VAR_DIR.joinpath(*parts)
    import os
    os.makedirs(long_path(p.parent), exist_ok=True)
    return p

def fig_dir(*parts):
    """fig_dir(LAG_TAG, "00_core_checks") -> OUT_DIR/LAG_TAG/00_core_checks/ (created)."""
    d = OUT_DIR.joinpath(*parts)
    import os
    os.makedirs(long_path(d), exist_ok=True)
    return d
