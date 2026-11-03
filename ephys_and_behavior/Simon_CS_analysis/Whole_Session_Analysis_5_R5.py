# %%
# ██╗███╗   ███╗██████╗  ██████╗ ██████╗ ████████╗ █████╗ ████████╗██╗ ██████╗ ███╗   ██╗
# ██║████╗ ████║██╔══██╗██╔═══██╗██╔══██╗╚══██╔══╝██╔══██╗╚══██╔══╝██║██╔═══██╗████╗  ██║
# ██║██╔████╔██║██████╔╝██║   ██║██████╔╝   ██║   ███████║   ██║   ██║██║   ██║██╔██╗ ██║
# ██║██║╚██╔╝██║██╔═══╝ ██║   ██║██╔══██╗   ██║   ██╔══██║   ██║   ██║██║   ██║██║╚██╗██║
# ██║██║ ╚═╝ ██║██║     ╚██████╔╝██║  ██║   ██║   ██║  ██║   ██║   ██║╚██████╔╝██║ ╚████║
# ╚═╝╚═╝     ╚═╝╚═╝      ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝
                                                                                       
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.signal import convolve
import matplotlib.patches as patches
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter, argrelextrema
import random
import math
from matplotlib.colors import LinearSegmentedColormap
import pickle
from scipy.signal import find_peaks, correlate
from scipy import signal
from matplotlib.widgets import Cursor
%matplotlib qt
from sklearn.manifold import TSNE
from sklearn.mixture import GaussianMixture
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.widgets import Button


# CONFIGURATION
file_path = r"C:\Users\simon\Desktop\Data\ephys_and_behavior"  # Your file path
save_path = r"C:\Users\simon\Desktop\Data\ephys_and_behavior\1_CoteDor\R5"
mat_file = f'{save_path}\\session_reaching_data_paw.mat'
target_mouse = "1_CoteDor"
target_sess = "R5"
output_file = f'{save_path}\\kinematic_data.parquet'




#  ██╗ ██╗     ███████╗██╗   ██╗███╗   ██╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗
# ████████╗    ██╔════╝██║   ██║████╗  ██║██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║
# ╚██╔═██╔╝    █████╗  ██║   ██║██╔██╗ ██║██║        ██║   ██║██║   ██║██╔██╗ ██║
# ████████╗    ██╔══╝  ██║   ██║██║╚██╗██║██║        ██║   ██║██║   ██║██║╚██╗██║
# ╚██╔═██╔╝    ██║     ╚██████╔╝██║ ╚████║╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║
#  ╚═╝ ╚═╝     ╚═╝      ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝
                                                                               
                                                                               
                                                                               
                                                                               
                                                                               
                                                                               
                                                                               
# --- RECURSIVE DECODER ---
def recursive_resolve(f, item):
    """
    Recursively follows HDF5 references until it finds actual data (numbers/strings).
    Handles: References, Arrays of References, and Nested Cells.
    """
    # 1. unwrapping: array([ref]) -> ref
    if isinstance(item, np.ndarray) and item.size == 1 and item.dtype == 'object':
        item = item.item()

    # 2. Dereferencing: Ref -> Data
    if isinstance(item, h5py.h5r.Reference):
        try:
            # Jump to the data location
            target = f[item]
            # Read the value
            data = target[:]
            # RECURSE: The result might be *another* reference
            return recursive_resolve(f, data)
        except Exception:
            return [] # Dead reference or empty

    # 3. Handling Arrays of Objects (e.g. string lists or jagged arrays)
    if isinstance(item, np.ndarray):
        # If numeric, we are done! Return the array.
        if item.dtype != 'object':
            return item
        # If it's an array of objects/refs, recurse on every element
        return [recursive_resolve(f, x) for x in item.flatten()]

    # 4. Base Case (Strings, Ints, Floats)
    return item

def decode_matlab_ref(f, ref):
    """
    Smart decoder:
    - If ref points to a string (1D ints, small-ish): Returns str
    - If ref points to a matrix (2D/3D ints/floats): Returns numpy array or list
    """
    try:
        if not ref:
            return None
        
        # Follow the reference
        ds = f[ref]
        val = ds[:]
        
        # --- CASE 1: MULTIDIMENSIONAL DATA (Definitely not a string) ---
        # If it has 2+ dimensions (e.g. 3x2000), it's data. Keep it as an array.
        if val.ndim > 1:
            return val  # Return the raw numpy array

        # --- CASE 2: 1D DATA ---
        # It could be a string OR a 1D vector of data.
        # Heuristic: If it's huge (>1000 elements), it's likely data, not a label.
        if val.size > 1000:
            return val

        # --- CASE 3: STRINGS ---
        # If it's small, 1D, and integer-like, it's likely a MATLAB char array.
        if np.issubdtype(val.dtype, np.integer):
            try:
                # Try to decode as text
                return "".join(chr(c) for c in val.flatten() if c != 0)
            except:
                # If decoding fails (weird chars), return original data
                return val
        
        # Fallback: Just return the value
        return val
        
    except Exception:
        return None
    
def force_resolve(f, item):
    """
    Forces HDF5 items (Datasets/Groups) into memory as numpy arrays or dicts.
    Then runs recursive_resolve to handle any internal references.
    """
    # If it's a Group (struct), convert to dict
    if isinstance(item, h5py.Group):
        return {k: force_resolve(f, item[k]) for k in item.keys()}
    
    # If it's a Dataset (array), read strictly into memory
    if isinstance(item, h5py.Dataset):
        data = item[()] # Load into memory
        return recursive_resolve(f, data) # Resolve internal refs
    
    # Otherwise, it's already a value or ref, let recursive_resolve handle it
    return recursive_resolve(f, item)

def decode_ascii_cell(nested_list):
    try:
        # Flatten [[80], [107], [67]] into [80, 107, 67]
        flat_list = [item[0] for item in nested_list]
        # Convert to characters and join into a string
        return "".join(chr(int(c)) for c in flat_list if c != 0)
    except Exception:
        return "Unknown"

data_dict = {}
keys_to_extract = ['reaches', 'session', 'trials']

with h5py.File(mat_file, 'r') as f:
    print(f"Opening {mat_file}...")
    for key in keys_to_extract:
        if key in f:
            print(f"Extracting '{key}'...")
            # Use force_resolve instead of just recursive_resolve
            data_dict[key] = force_resolve(f, f[key])
        else:
            print(f"Warning: Key '{key}' not found.")

# Reconstruct timestamps (NEW)
if "trials" in data_dict and "session" in data_dict:
    try:
        ts_matrix = data_dict["trials"]["sync"]["frames_timestamps"]
        # Ensure it's a numpy array for easy referencing
        trial_lengths = np.array(data_dict["trials"]["size_trial"]).flatten()
        
        valid_ts = []
        # Flatten based on valid lengths
        for i, length in enumerate(trial_lengths):
            n_frames = int(length)
            # Row-based indexing: (trial_idx, 0:n_frames)
            valid_ts.append(ts_matrix[0:n_frames, i])
            
        reconstructed = np.concatenate(valid_ts)
        data_dict["session"]["reconstructed_timestamps"] = reconstructed
        print(f"Reconstructed timestamps shape: {reconstructed.shape}")
        
    except Exception as e:
        print(f"Failed to reconstruct timestamps: {e}")

# Save as Dictionary (safest for varying lengths)
pkl_path = f'{save_path}\\extracted_data.pkl'

try:
    with open(pkl_path, 'wb') as pfile:
        pickle.dump(data_dict, pfile)
    print(f"Successfully saved to: {pkl_path}")
    print(f"Keys: {list(data_dict.keys())}")
except TypeError as e:
    print(f"STILL FAILED: {e}") 
    # Debugging: check types in dict
    for k, v in data_dict.items():
        print(f"{k}: {type(v)}")
        if isinstance(v, (list, np.ndarray)) and len(v) > 0:
             print(f"  First element type: {type(v[0])}")
# %%

# ██╗  ██╗███████╗██╗   ██╗███████╗    ███████╗██╗   ██╗███╗   ███╗███╗   ███╗ █████╗ ██████╗ ██╗   ██╗    
# ██║ ██╔╝██╔════╝╚██╗ ██╔╝██╔════╝    ██╔════╝██║   ██║████╗ ████║████╗ ████║██╔══██╗██╔══██╗╚██╗ ██╔╝    
# █████╔╝ █████╗   ╚████╔╝ ███████╗    ███████╗██║   ██║██╔████╔██║██╔████╔██║███████║██████╔╝ ╚████╔╝     
# ██╔═██╗ ██╔══╝    ╚██╔╝  ╚════██║    ╚════██║██║   ██║██║╚██╔╝██║██║╚██╔╝██║██╔══██║██╔══██╗  ╚██╔╝      
# ██║  ██╗███████╗   ██║   ███████║    ███████║╚██████╔╝██║ ╚═╝ ██║██║ ╚═╝ ██║██║  ██║██║  ██║   ██║       
# ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚══════╝    ╚══════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝       
                                                                                                         

keys = list(data_dict["session"].keys())
print("Extracted Keys:", keys)
#%%

plt.plot(data_dict["session"]["reconstructed_timestamps"], data_dict["session"]["ppDom_paw_sess"][0])
# %%

from scipy.signal import argrelextrema

ONSET_THRESHOLD_PCT = 0.05
ONSET_STABILITY_FRAMES = 10
# --- 2. SPEED-BASED DETECTION ---
def get_onset_idx_speed(trajectory_3d):
    # 1. Calculate Euclidean Speed
    vel = np.diff(trajectory_3d, axis=0)
    speed = np.linalg.norm(vel, axis=1) # No smoothing as requested
    
    # 2. Find Peak (10-150 range)
    # Ensure window is large enough
    search_start = 10
    search_end = min(150, len(speed))
    
    roi = speed[search_start:search_end]
    if len(roi) == 0: return 0
    
    peak_idx = np.argmax(roi) + search_start
    peak_val = speed[peak_idx]
    
    # 3. Valley Search
    threshold = peak_val * ONSET_THRESHOLD_PCT
    pre_peak_speed = speed[0:peak_idx]
    valley_indices = argrelextrema(pre_peak_speed, np.less, order=1)[0]
    
    onset_idx = 0 # Default fallback
    
    if len(valley_indices) > 0:
        found_deep_valley = False
        # Search Right -> Left
        for v_idx in reversed(valley_indices):
            if speed[v_idx] < threshold:
                # Check 2: Stability (Was it low for n frames before/at this point?)
                check_start = max(0, v_idx - ONSET_STABILITY_FRAMES)
                if np.median(speed[check_start : v_idx+1]) < threshold:
                    onset_idx = v_idx
                    found_deep_valley = True
                    break
        
        if not found_deep_valley:
            onset_idx = valley_indices[-1] # Fallback: closest valley
            
    return onset_idx

# --- 3. PC1-BASED DETECTION ---
def get_onset_idx_pc1(trajectory_3d, global_pca):
    # 1. Project to Global PC1
    # trajectory_3d is already cleaned/interpolated ideally
    pc1_trace = global_pca.transform(trajectory_3d).flatten()
    
    # 2. Calculate PC1 Speed
    speed_pc1 = np.abs(np.diff(pc1_trace)) # No smoothing
    
    # 3. Find Peak
    search_start = 10
    search_end = min(150, len(speed_pc1))
    
    roi = speed_pc1[search_start:search_end]
    if len(roi) == 0: return 0
    
    peak_idx = np.argmax(roi) + search_start
    peak_val = speed_pc1[peak_idx]
    
    # 4. Valley Search
    threshold = peak_val * ONSET_THRESHOLD_PCT
    pre_peak_speed = speed_pc1[0:peak_idx]
    valley_indices = argrelextrema(pre_peak_speed, np.less, order=1)[0]
    
    onset_idx = 0
    
    if len(valley_indices) > 0:
        found_deep_valley = False
        for v_idx in reversed(valley_indices):
            if speed_pc1[v_idx] < threshold:
                check_start = max(0, v_idx - ONSET_STABILITY_FRAMES)
                if np.median(speed_pc1[check_start : v_idx+1]) < threshold:
                    onset_idx = v_idx
                    found_deep_valley = True
                    break
        if not found_deep_valley:
            onset_idx = valley_indices[-1]

    return onset_idx

# --- 4. ARBITRATION LOGIC (THE MASTER FUNCTION) ---
def get_arbitrated_onset(trajectory_3d, global_pca):
    """
    Decides between Speed-based and PC1-based onset.
    """
    # Get candidates
    idx_speed = get_onset_idx_speed(trajectory_3d)
    idx_pc1 = get_onset_idx_pc1(trajectory_3d, global_pca)
    
    # Check discrepancy
    if abs(idx_speed - idx_pc1) > 5:
        # Get PC1 Position trace (not speed)
        pc1_pos = global_pca.transform(trajectory_3d).flatten()
        
        # Safety check for indices out of bounds
        max_idx = len(pc1_pos) - 1
        val_speed = pc1_pos[min(idx_speed, max_idx)]
        val_pc1 = pc1_pos[min(idx_pc1, max_idx)]
        
        # Pick the index corresponding to the "Lower" (Start) position
        # Assuming lower PC1 means "retracted" or "start position"
        # If PC1 is aligned with forward reach, start might be low or high depending on axis alignment
        # Assuming standard "reach starts from low X/PC1"
        if val_pc1 < val_speed:
            return idx_pc1
        else:
            return idx_speed
    else:
        # Agreement (or close enough) -> Default to Speed
        return idx_speed

                                                                                      
                                                                                      
                                                                                      

                             #%%                                                         
                                                                                      
                                                                                      
                                                                                      
# ██████╗ ███████╗ █████╗ ██╗  ██╗    ██╗    ██╗██╗███╗   ██╗██████╗  ██████╗ ██╗    ██╗
# ██╔══██╗██╔════╝██╔══██╗██║ ██╔╝    ██║    ██║██║████╗  ██║██╔══██╗██╔═══██╗██║    ██║
# ██████╔╝█████╗  ███████║█████╔╝     ██║ █╗ ██║██║██╔██╗ ██║██║  ██║██║   ██║██║ █╗ ██║
# ██╔═══╝ ██╔══╝  ██╔══██║██╔═██╗     ██║███╗██║██║██║╚██╗██║██║  ██║██║   ██║██║███╗██║
# ██║     ███████╗██║  ██║██║  ██╗    ╚███╔███╔╝██║██║ ╚████║██████╔╝╚██████╔╝╚███╔███╔╝
# ╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝     ╚══╝╚══╝ ╚═╝╚═╝  ╚═══╝╚═════╝  ╚═════╝  ╚══╝╚══╝ 
                                                                                      
# Main Analysis Logic
all_reaches_metadata = []
grooming_metadata = []
small_mvmt_metadata = []

# Define mapping from your data keys to logical names
# (Assuming 'ppDom_paw_sess' is the Dominant paw)
paw_mapping = {
    'Dom':  {'data_key': 'ppDom_paw_sess',  'opponent': 'NonDom'},
    'NonDom': {'data_key': 'ppNonDom_paw_sess', 'opponent': 'Dom'}
}

for paw_role in ['Dom', 'NonDom']:
    current_paw_conf = paw_mapping[paw_role]
    opponent_role = current_paw_conf['opponent']
    
    # Check if keys exist
    if current_paw_conf['data_key'] not in data_dict['session']:
        print(f"Skipping {paw_role}: Key {current_paw_conf['data_key']} not found.")
        continue
        
    print(f"Analyzing {paw_role} Paw...")

    # 1. Prepare Data
    # Shape is likely (3, N_samples). We transpose to (N_samples, 3) for PCA
    X_paw = data_dict['session'][current_paw_conf['data_key']].T 
    X_opponent = data_dict['session'][paw_mapping[opponent_role]['data_key']].T

    # Handle NaNs (PCA cannot handle NaNs)
    df_paw = pd.DataFrame(X_paw, columns=['x', 'y', 'z']).interpolate(method='linear', limit_direction='both')
    X_clean = df_paw.values
    
    # 2. PCA for main movement direction
    pca = PCA(n_components=1)
    pc1 = pca.fit_transform(X_clean).flatten()
    print(f"  Variance explained by PC1: {pca.explained_variance_ratio_[0]:.4f}")

    # Align PC1 with X direction (assuming forward reach is +X)
    # If correlation is negative, flip PC1 so positive values = forward reach
    if np.corrcoef(X_clean[:, 0], pc1)[0, 1] < 0:
        pc1 = -pc1
    
    # 3. Peak Detection
    # 'prominence' and 'width' tuned for PCA units. Adjust if necessary.
    peaks, _ = find_peaks(pc1, prominence=20, distance=30, width=10)
    print(f"  Found {len(peaks)} peaks.")

    for peak in peaks:
        # Define window around peak for correlation check
        window_size = 50
        start_w = max(0, peak - window_size)
        end_w = min(len(pc1), peak + window_size)
        
        # Get Z-coordinates for correlation (as per your snippet)
        # Z is index 2
        paw_z_segment = X_paw[start_w:end_w, 2]
        opponent_z_segment = X_opponent[start_w:end_w, 2]

        # --- Correlation Check ---
        # Skip if segment is too short or all NaNs
        if len(paw_z_segment) < 10 or np.isnan(paw_z_segment).all():
            continue

        # Pearson Correlation
        # Create temp DF to handle NaNs gracefully
        df_corr = pd.DataFrame({'a': paw_z_segment, 'b': opponent_z_segment})
        pearson_r = df_corr['a'].corr(df_corr['b'])
        if np.isnan(pearson_r): pearson_r = 0

        # Correlation Lag (using cross-correlation)
        # Fill NaNs with 0 for signal.correlate
        a_filled = np.nan_to_num(paw_z_segment)
        b_filled = np.nan_to_num(opponent_z_segment)
        correlation = signal.correlate(a_filled, b_filled, mode="full")
        lag = (np.argmax(correlation) - (len(a_filled) - 1)) if len(correlation) > 0 else 0
        
        # --- Distance Check ---
        # Distance between paws in X axis
        dist_array = np.abs(X_paw[start_w:end_w, 0] - X_opponent[start_w:end_w, 0])
        distance_between_paws = np.nanmedian(dist_array) if len(dist_array) > 0 else 0
        
        # --- Logic Conditions ---
        is_uncorrelated = (pearson_r < 0.6)
        is_opponent_stationary = (pearson_r >= 0.6 and distance_between_paws > 30)
        is_far_enough = (distance_between_paws > 30)
        
        # --- ONSET DETECTION LOGIC ---
        w_start = peak - 150
        w_end = peak + 50
        
        # Safe slicing relative to session end
        # X_clean is defined earlier
        traj_3d_segment = X_clean[max(0, w_start):min(len(X_clean), w_end)]
        
        # Calculate Onsets
        if len(traj_3d_segment) > 10: # Minimum length check
            onset_speed = get_onset_idx_speed(traj_3d_segment)
            onset_pc1 = get_onset_idx_pc1(traj_3d_segment, pca)
            onset_arb = get_arbitrated_onset(traj_3d_segment, pca)
        else:
            onset_speed, onset_pc1, onset_arb = 0, 0, 0
        

        onset_global = w_start + onset_arb
        
        # Calculate Amplitude (X displacement from onset to peak)
        mov_amp = X_paw[peak, 0] - X_paw[onset_global, 0]
        
        # Metadata dictionary
        meta = {
            'paw': paw_role,
            'peak_idx': peak,
            'peak_timestamp': data_dict['session']['reconstructed_timestamps'][peak] if 'reconstructed_timestamps' in data_dict['session'] and peak < len(data_dict['session']['reconstructed_timestamps']) else None,
            'onset_idx': onset_global,
            'onset_timestamp': data_dict['session']['reconstructed_timestamps'][onset_global] if 'reconstructed_timestamps' in data_dict['session'] and onset_global < len(data_dict['session']['reconstructed_timestamps']) else None,
            'amplitude': mov_amp,
            'pearson_r': pearson_r,
            'distance': distance_between_paws,
            'lag': lag
        }

        # --- Classification ---
        if mov_amp > 60: # Large movement threshold
            if (is_uncorrelated or is_opponent_stationary) and is_far_enough:
                all_reaches_metadata.append(meta)
            else:
                grooming_metadata.append(meta)
        else:
            small_mvmt_metadata.append(meta)

print(f"\nSummary:")
print(f"  Valid Reaches: {len(all_reaches_metadata)}")
print(f"  Grooming/Synced: {len(grooming_metadata)}")
print(f"  Small Movements: {len(small_mvmt_metadata)}")
# %%

# --- PLOTTING FUNCTION ---
def plot_reaches_grid(metadata_list, data_storage, session_name="Session", batch_size=36):
    """
    Plots reaches in grid batches.
    Args:
        metadata_list: List of dicts containing 'peak_idx', 'onset_idx', etc.
        data_storage: The main data_dict containing 'session' -> 'ppDom_paw_sess' etc.
        session_name: Label for the plot title.
        batch_size: Number of plots per figure.
    """
    
    if not metadata_list:
        print("No metadata to plot.")
        return

    num_batches = (len(metadata_list) + batch_size - 1) // batch_size
    print(f"Plotting {len(metadata_list)} reaches in {num_batches} batches...")

    # Define layout
    n_cols = 6
    n_rows = 6  # 5*4 = 20 plots per page

    for batch_i in range(num_batches):
        start = batch_i * batch_size
        end = min(start + batch_size, len(metadata_list))
        batch_meta = metadata_list[start:end]
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 10), constrained_layout=True)
        axes = axes.flatten()
        
        for i, meta in enumerate(batch_meta):
            ax = axes[i]
            
            # 1. Identify Data Source
            paw_role = meta['paw'] # 'Dom' or 'NonDom'
            
            # Map roles to keys in data_dict
            if paw_role == 'Dom':
                key_main = 'ppDom_paw_sess'
                key_opp  = 'ppNonDom_paw_sess'
            else:
                key_main = 'ppNonDom_paw_sess'
                key_opp  = 'ppDom_paw_sess'
            
            # 2. Define Window indices (e.g. onset - 50 frames to peak + 50 frames)
            # We use the peak_idx from metadata to center the view
            peak_idx = meta['peak_idx']
            onset_idx = meta.get('onset_idx', peak_idx - 20)
            
            w_start = max(0, onset_idx - 150)
            w_end = min(data_storage['session'][key_main].shape[1], peak_idx + 50)
            
            # 3. Slice Data (X-coordinate is index 0)
            # Data shape is (3, n_frames), so we want row 0, cols w_start:w_end
            main_trace = data_storage['session'][key_main][0, w_start:w_end]
            opp_trace  = data_storage['session'][key_opp][0, w_start:w_end]
            
            # Time axis relative to window start
            time_axis = np.arange(w_start, w_end)
            
            # 4. Plot
            # Main Paw (Left Axis)
            ax.plot(time_axis, main_trace, 'b-', label=f'{paw_role} (X)', linewidth=1.5)
            ax.set_ylabel(f'{paw_role} Pos', color='b', fontsize=8)
            ax.tick_params(axis='y', labelcolor='b', labelsize=8)
            
            # Opponent Paw (Right Axis)
            ax2 = ax.twinx()
            ax2.plot(time_axis, opp_trace, 'orange', alpha=0.7, label='Opp (X)', linewidth=1)
            ax2.tick_params(axis='y', labelsize=6)
            
            # 5. Markers
            # Mark Peak
            if w_start <= peak_idx < w_end:
                 ax.plot(peak_idx, data_storage['session'][key_main][0, peak_idx], 'ro', markersize=4)
            
            # Mark Onset
            if w_start <= onset_idx < w_end:
                ax.axvline(x=onset_idx, color='green', linestyle='--', alpha=0.8, label='Onset')

            # 6. Title / Info
            # Amplitude, Correlation, etc.
            info = f"Amp:{meta.get('amplitude',0):.0f} | r:{meta.get('pearson_r',0):.2f}"
            ax.set_title(f"Idx: {peak_idx}\n{info}", fontsize=9)
            
            # Legend only on first plot
            if i == 0:
                lines, labels = ax.get_legend_handles_labels()
                lines2, labels2 = ax2.get_legend_handles_labels()
                ax.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=6)

        # Hide unused subplots
        for j in range(len(batch_meta), len(axes)):
            axes[j].axis('off')

        fig.suptitle(f"Reaches Batch {batch_i+1}/{num_batches} - {session_name}", fontsize=14)
        plt.show()

#%%
# Example Usage:
plot_reaches_grid(all_reaches_metadata, data_dict, session_name=target_sess)
# %%^
                                                                                                                         
                                                                                                                         
                                                                                                                         
                                                                                                                         
                                                                                                                         
                                                                                                                         
                                                                                                                         
# ██╗███╗   ██╗████████╗███████╗██████╗  █████╗  ██████╗████████╗██╗██╗   ██╗███████╗    ██████╗ ██╗      ██████╗ ████████╗
# ██║████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔══██╗██╔════╝╚══██╔══╝██║██║   ██║██╔════╝    ██╔══██╗██║     ██╔═══██╗╚══██╔══╝
# ██║██╔██╗ ██║   ██║   █████╗  ██████╔╝███████║██║        ██║   ██║██║   ██║█████╗      ██████╔╝██║     ██║   ██║   ██║   
# ██║██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██╔══██║██║        ██║   ██║╚██╗ ██╔╝██╔══╝      ██╔═══╝ ██║     ██║   ██║   ██║   
# ██║██║ ╚████║   ██║   ███████╗██║  ██║██║  ██║╚██████╗   ██║   ██║ ╚████╔╝ ███████╗    ██║     ███████╗╚██████╔╝   ██║   
# ╚═╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝   ╚═╝   ╚═╝  ╚═══╝  ╚══════╝    ╚═╝     ╚══════╝ ╚═════╝    ╚═╝   
                                                                                                                         
# --- INTERACTIVE PEAK INSPECTION PLOT ---
# Ensure you are using %matplotlib qt for interactivity

def setup_click_inspector(data_dict, timestamps):
    """
    Sets up a main timeline plot. Clicking a peak opens a detailed inspection window.
    Peaks are color-coded by classification (Accepted, Grooming, Rejected).
    """
    print("Setting up Interactive Inspector...")
    
    # --- 1. PREPARE DATA (Dominant Paw Focus) ---
    dom_key = 'ppDom_paw_sess'
    ndom_key = 'ppNonDom_paw_sess' 
    
    if dom_key not in data_dict['session']:
        print(f"Error: {dom_key} not found.")
        return

    # Get Raw Data
    X_dom = data_dict['session'][dom_key].T
    # Interpolate for PCA
    df_dom = pd.DataFrame(X_dom, columns=['x', 'y', 'z']).interpolate(method='linear', limit_direction='both')
    X_clean = df_dom.values
    
    # Calculate PC1
    pca = PCA(n_components=1)
    pc1 = pca.fit_transform(X_clean).flatten()
    
    # Align direction
    if np.corrcoef(X_clean[:, 0], pc1)[0, 1] < 0:
        pc1 = -pc1
        
    # Find Peaks
    peaks, props = find_peaks(pc1, prominence=20, distance=30, width=10)
    
    # Prepare Time Vector
    min_len = min(len(timestamps), len(pc1))
    t = timestamps[:min_len]
    trace_pc1 = pc1[:min_len]
    trace_dom_x = data_dict['session'][dom_key][0, :min_len]
    trace_ndom_x = data_dict['session'][ndom_key][0, :min_len]
    
    # --- 2. CLASSIFY PEAKS FOR PLOTTING ---
    # Create lists for scatter plots
    cat_indices = {'Accepted': [], 'Grooming': [], 'Rejected': []}
    
    accepted_indices = set(m['peak_idx'] for m in all_reaches_metadata if m['paw'] == 'Dom')
    grooming_indices = set(m['peak_idx'] for m in grooming_metadata if m['paw'] == 'Dom')
    
    for p in peaks:
        if p >= min_len: continue
        
        if p in accepted_indices:
            cat_indices['Accepted'].append(p)
        elif p in grooming_indices:
            cat_indices['Grooming'].append(p)
        else:
            cat_indices['Rejected'].append(p)

    # --- 3. MAIN TIMELINE PLOT ---
    fig_main, ax_main = plt.subplots(figsize=(18, 6))
    plt.subplots_adjust(bottom=0.2)
    
    # Plot Full Session Trace
    ax_main.plot(t, trace_pc1, 'k-', lw=0.5, alpha=0.4, label='Dom PC1')
    ax_main.set_ylabel('PC1 (Arb)', color='k')
    ax_main.hlines(xmin=min(t), xmax=max(t),y = 20, colors='gray', linestyles='--', alpha=0.5, linewidth=1) 
    # NEW: Plot Non-Dominant X for global context (secondary axis)
    ax_ndom = ax_main.twinx()
    ax_ndom.plot(t, trace_ndom_x, 'orange', lw=0.5, alpha=0.3, label='NDom X')
    ax_ndom.set_ylabel('NDom X (mm)', color='orange')
    ax_ndom.tick_params(axis='y', labelcolor='orange')

    # FIX INTERACTIVITY: Ensure ax_main is on top of ax_ndom so clicks reach the scatter
    ax_main.set_zorder(ax_ndom.get_zorder() + 1) 
    ax_main.patch.set_visible(False) # Make background transparent so ax_ndom is visible

    # Scatter Plots by Category
    colors = {'Accepted': 'green', 'Grooming': 'orange', 'Rejected': 'red'}
    markers = {'Accepted': 'o', 'Grooming': 's', 'Rejected': 'x'}
    scatters = {}

    # NEW: Plot Water Delivery Timestamps
    if "trials" in data_dict and "sync" in data_dict["trials"]:
        try:
            # Check structure specifically for timestamp_delivered
            sync_data = data_dict["trials"]["sync"]
            water_timestamps = None
            
            # recursive_resolve or force_resolve might have left it as dict or array
            if isinstance(sync_data, dict) and "timestamp_delivered" in sync_data:
                 water_timestamps = sync_data["timestamp_delivered"]
            
            if water_timestamps is not None:
                # Handle if it's a nested array or list of arrays
                if hasattr(water_timestamps, 'flatten'):
                     water_timestamps = water_timestamps.flatten()
                
                # Filter valid timestamps (remove NaNs or zeros if applicable)
                if len(water_timestamps) > 0:
                     valid_wt = [wt for wt in water_timestamps if wt > 0 and wt < t[-1]]
                     
                     if len(valid_wt) > 0:
                        # Plot as vertical lines
                        ax_main.vlines(valid_wt, ymin=trace_pc1.min(), ymax=trace_pc1.max(), 
                                      colors='cyan', linestyles=':', alpha=0.6, linewidth=1, zorder=1, label='Water')
            
        except Exception as e:
            print(f"Could not plot water timestamps: {e}")
    
    for cat, idx_list in cat_indices.items():
        if not idx_list: continue
        
        idx_array = np.array(idx_list)
        # Scatter is clickable (picker=True)
        sc = ax_main.scatter(t[idx_array], trace_pc1[idx_array], 
                             color=colors[cat], marker=markers[cat], s=60, 
                             picker=8, zorder=10, label=f'{cat} ({len(idx_list)})')
        scatters[cat] = (sc, idx_array) # Store object and original indices
    
    ax_main.legend(loc='upper right')
    ax_main.set_title("MAIN TIMELINE: Click Peaks to Inspect\n(Green=Reach, Orange=Grooming/Sync, Red=Small/Rejected)", fontsize=14, fontweight='bold')
    ax_main.set_xlabel("Time (s)")

    # --- 4. DETAIL WINDOW FUNCTION ---
    def open_detail_window(peak_idx, category):
        """Spawns a new figure for the selected peak."""
        # Find Metadata
        meta = {}
        target_list = all_reaches_metadata if category == 'Accepted' else (grooming_metadata if category == 'Grooming' else small_mvmt_metadata)
        
        for m in target_list:
            if m['peak_idx'] == peak_idx and m['paw'] == 'Dom':
                meta = m
                break
        
        # Window Range
        window = 100
        w_start = max(0, peak_idx - window)
        w_end = min(min_len, peak_idx + window)
        
        # Slice Data
        t_slice = t[w_start:w_end]
        pc1_slice = trace_pc1[w_start:w_end]
        dom_slice = trace_dom_x[w_start:w_end]
        ndom_slice = trace_ndom_x[w_start:w_end]
        
        # Create Popup Figure
        fig_detail, (ax_d1, ax_d2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        fig_detail.suptitle(f"Frame {peak_idx} ({t[peak_idx]:.2f}s) - {category}", fontsize=14, color=colors.get(category, 'black'))
        
        # Plot 1: PC1 & Onset
        ax_d1.plot(t_slice, pc1_slice, 'k-', alpha=0.7, label='PC1')
        ax_d1.plot(t[peak_idx], trace_pc1[peak_idx], color=colors.get(category,'red'), marker='o', label='Peak')
        ax_d1.set_ylabel('PC1')
        ax_d1.legend(loc='upper left')
        ax_d1.grid(True, alpha=0.3)
        
        # Mark Onset if available
        if 'onset_idx' in meta:
            onset_t = t[meta['onset_idx']]
            ax_d1.axvline(onset_t, color='blue', linestyle='--', label='Onset')
            ax_d1.text(onset_t, min(pc1_slice), "Onset", color='blue', fontsize=8)

        # Plot 2: Raw Coordinates (Dom vs NDom)
        ax_d2.plot(t_slice, dom_slice, 'b-', label='Dom X (Main)')
        ax_d2.plot(t_slice, ndom_slice, 'orange', label='NDom X (Opposite)')
        ax_d2.set_ylabel('X Position (mm)')
        ax_d2.set_xlabel('Time (s)')
        
        # Add a visual bar for amplitude
        if 'onset_idx' in meta:
             onset_frame = meta['onset_idx']
             if onset_frame >= w_start and onset_frame < w_end:
                 peak_val_x = trace_dom_x[peak_idx]
                 onset_val_x = trace_dom_x[onset_frame]
                 ax_d2.plot([t[onset_frame], t[peak_idx]], [onset_val_x, peak_val_x], 'k--', alpha=0.5, label='Amplitude')

        ax_d2.legend(loc='upper left')
        ax_d2.grid(True, alpha=0.3)
        
        # Text Info Box
        info_text = (
            f"Metrics:\n"
            f"Amplitude: {meta.get('amplitude', 0):.1f}\n"
            f"Correlation (r): {meta.get('pearson_r', 0):.2f}\n"
            f"Lag: {meta.get('lag', 0)}\n"
            f"Paw Dist: {meta.get('distance', 0):.1f}"
        )
        ax_d2.text(0.95, 0.95, info_text, transform=ax_d2.transAxes, 
                   fontsize=10, verticalalignment='top', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

        fig_detail.show()

    # --- 4. CLICK HANDLER ---
    def on_pick(event):
        # Determine which scatter plot was clicked
        clicked_sc = event.artist
        category = None
        indices = None
        
        for cat, (sc_obj, idx_list) in scatters.items():
            if sc_obj == clicked_sc:
                category = cat
                indices = idx_list
                break
        
        if category is None: return
            
        # Get index of the point clicked in the specific scatter array
        ind = event.ind[0]
        # Map back to original frame index using the stored index list
        peak_idx = indices[ind]
        
        print(f"Clicked {category} Peak at Frame {peak_idx}")
        open_detail_window(peak_idx, category)

    fig_main.canvas.mpl_connect('pick_event', on_pick)
    plt.show()

# Run
if "reconstructed_timestamps" in data_dict["session"]:
    setup_click_inspector(data_dict, data_dict["session"]["reconstructed_timestamps"])
else:
    print("Timestamps needed.")
# %%


# ████████╗███████╗███╗   ██╗███████╗    
# ╚══██╔══╝██╔════╝████╗  ██║██╔════╝    
#    ██║   ███████╗██╔██╗ ██║█████╗      
#    ██║   ╚════██║██║╚██╗██║██╔══╝      
#    ██║   ███████║██║ ╚████║███████╗    
#    ╚═╝   ╚══════╝╚═╝  ╚═══╝╚══════╝    
                                       
# --- 1. DATA PREPARATION (Run Once) ---
# Check if windows are already extracted
if 'X_windows_3d' not in locals():
    print("Extracting reach windows...")
    window_before = 150
    window_after = 50
    min_length = window_before + window_after
    
    valid_reaches_list = []      # List of (3, 200) arrays
    valid_metadata_list = []     # List of corresponding metadata dicts
    
    # Use Dominant Paw (X, Y, Z) for clustering
    dom_key = 'ppDom_paw_sess' 
    if dom_key in data_dict['session']:
        raw_data = data_dict['session'][dom_key] # Shape (3, n_frames)
        
        # We need to iterate carefully to keep alignment
        for m in all_reaches_metadata:
            if m['paw'] == 'Dom':
                idx = m['peak_idx']
                # Check bounds
                if idx - window_before >= 0 and idx + window_after < raw_data.shape[1]:
                    # Extract (3, 200) window
                    w = raw_data[:, idx-window_before : idx+window_after]
                    if not np.isnan(w).any():
                        valid_reaches_list.append(w)
                        valid_metadata_list.append(m)
    
    if len(valid_reaches_list) > 10:
        X_windows_3d = np.array(valid_reaches_list) # Shape (N, 3, 200)
        # Flatten for ML: (N, 600)
        X_windows_flat = X_windows_3d.reshape(len(X_windows_3d), -1)
        print(f"Data prepared: {X_windows_flat.shape}")
    else:
        print("Not enough valid reaches found.")
        X_windows_flat = None
        X_windows_3d = None

# --- 2. TSNE COMPUTATION (Run Once) ---
if 'tsne_result' not in locals() and X_windows_flat is not None:
    print("Computing t-SNE (this takes time)...")
    p_val = min(30, len(X_windows_flat) - 1)
    tsne = TSNE(n_components=2, perplexity=p_val, random_state=42, init='pca', learning_rate='auto')
    tsne_result = tsne.fit_transform(X_windows_flat)
    print("t-SNE done.")
elif 'tsne_result' in locals():
    print("Using existing t-SNE result.")

# --- 3. GMM CLUSTERING & INTERACTIVE PLOT ---
if 'tsne_result' in locals() and X_windows_3d is not None:
    print("Fitting GMM...")
    gmm = GaussianMixture(n_components=15, random_state=42)
    cluster_labels = gmm.fit_predict(tsne_result)
    
    # --- INTERACTIVE CLASS DEFINITION ---
    class InteractiveGMMPlot:
        def __init__(self, tsne_data, windows_3d, labels, metadata):
            self.tsne = tsne_data
            self.windows = windows_3d
            self.labels = labels
            self.metadata = metadata
            self.n_clusters = len(np.unique(labels))
            
            # Setup Main Figure (Scatter)
            self.fig_main, self.ax = plt.subplots(figsize=(10, 8))
            self.scatter = self.ax.scatter(self.tsne[:, 0], self.tsne[:, 1], 
                                           c=self.labels, cmap='nipy_spectral', 
                                           alpha=0.6, picker=5)
            self.fig_main.colorbar(self.scatter, label='Cluster ID')
            self.ax.set_title(f"t-SNE + GMM ({self.n_clusters} clusters)\nClick a point to see cluster average")
            
            # Connect Event
            self.cid = self.fig_main.canvas.mpl_connect('pick_event', self.on_pick)
            
            # Setup Detail Figure (4 Rows: X, Y, Z, Onset Hist)
            self.fig_avg, self.axes_avg = plt.subplots(4, 1, figsize=(7, 12), sharex=True, 
                                                       gridspec_kw={'height_ratios': [1, 1, 1, 0.6]})
            self.fig_avg.suptitle("Cluster Averages")
            plt.show()

        def on_pick(self, event):
            if event.artist != self.scatter: return
            ind = event.ind[0]
            selected_cluster = self.labels[ind]
            self.update_average_plots(selected_cluster)

        def update_average_plots(self, cluster_id):
            # Get indices for this cluster
            indices = np.where(self.labels == cluster_id)[0]
            cluster_data = self.windows[indices] # (N, 3, Time)
            n_points = len(indices)

            mean_traj = np.mean(cluster_data, axis=0) # (3, Time)
            std_traj = np.std(cluster_data, axis=0)
            
            t_axis = np.arange(-150, 50)
            dims = ['X', 'Y', 'Z']
            colors = ['black', 'green', 'magenta']
            
            # --- Plot Kinematics (Rows 0-2) ---
            for i in range(3):
                ax = self.axes_avg[i]
                ax.clear()
                # Mean + Shaded Error
                ax.plot(t_axis, mean_traj[i], color=colors[i], lw=2)
                ax.fill_between(t_axis, mean_traj[i]-std_traj[i], mean_traj[i]+std_traj[i], 
                                color=colors[i], alpha=0.02)
                
                # Faint Individual Traces
                subset_idx = np.random.choice(len(cluster_data), min(len(cluster_data), 50), replace=False)
                for _idx in subset_idx:
                    ax.plot(t_axis, cluster_data[_idx, i, :], color=colors[i], lw=0.5, alpha=0.5)
                
                ax.axvline(0, color='k', linestyle='--', alpha=0.5)
                ax.set_ylabel(f'{dims[i]} Position')
                ax.grid(True, alpha=0.3)
            
            # --- Plot Onset Histogram (Row 3) ---
            ax_hist = self.axes_avg[3]
            ax_hist.clear()
            
            # Collect relative onsets (Onset - Peak)
            onsets_rel = []
            for global_idx in indices:
                # Assuming metadata list aligns 1:1 with windows array
                if global_idx < len(self.metadata):
                    m = self.metadata[global_idx]
                    # Calculate relative onset
                    if 'onset_idx' in m and 'peak_idx' in m:
                        rel = m['onset_idx'] - m['peak_idx']
                        onsets_rel.append(rel)
            
            if len(onsets_rel) > 0:
                sns.histplot(onsets_rel, ax=ax_hist, kde=True, color='purple', bins=20, element="step")
                med_onset = np.median(onsets_rel)
                ax_hist.axvline(med_onset, color='purple', linestyle='--', lw=2, label=f'Med: {med_onset:.0f}')
                ax_hist.legend(loc='upper right', fontsize='small')
            
            ax_hist.set_ylabel('Count')
            ax_hist.set_xlabel('Frames from Peak (0)')
            ax_hist.set_xlim(-150, 50)
            ax_hist.grid(True, alpha=0.3)

            self.fig_avg.suptitle(f"Cluster {cluster_id} (n={n_points})")
            self.fig_avg.canvas.draw()

    # Launch App, passing aligned metadata
    app = InteractiveGMMPlot(tsne_result, X_windows_3d, cluster_labels, all_reaches_metadata)
else:
    print("Missing data required for clustering.")
# %%
# ███╗   ██╗███████╗██╗   ██╗██████╗  █████╗ ██╗          █████╗  ██████╗████████╗██╗██╗   ██╗██╗████████╗██╗   ██╗
# ████╗  ██║██╔════╝██║   ██║██╔══██╗██╔══██╗██║         ██╔══██╗██╔════╝╚══██╔══╝██║██║   ██║██║╚══██╔══╝╚██╗ ██╔╝
# ██╔██╗ ██║█████╗  ██║   ██║██████╔╝███████║██║         ███████║██║        ██║   ██║██║   ██║██║   ██║    ╚████╔╝ 
# ██║╚██╗██║██╔══╝  ██║   ██║██╔══██╗██╔══██║██║         ██╔══██║██║        ██║   ██║╚██╗ ██╔╝██║   ██║     ╚██╔╝  
# ██║ ╚████║███████╗╚██████╔╝██║  ██║██║  ██║███████╗    ██║  ██║╚██████╗   ██║   ██║ ╚████╔╝ ██║   ██║      ██║   
# ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝    ╚═╝  ╚═╝ ╚═════╝   ╚═╝   ╚═╝  ╚═══╝  ╚═╝   ╚═╝      ╚═╝   
                                                                                                                  
neural_mat_file = f'{save_path}\\eg_neurons_CB.mat'

data_neurons_dict = {}
keys_to_extract = ['eg_neurons']

with h5py.File(neural_mat_file, 'r') as f:
    print(f"Available keys in .mat file: {list(f.keys())}")
    print(f"Opening {neural_mat_file}...")
    for key in keys_to_extract:
        if key in f:
            print(f"Extracting '{key}'...")
            # Use force_resolve instead of just recursive_resolve
            data_neurons_dict[key] = force_resolve(f, f[key])
        else:
            print(f"Warning: Key '{key}' not found.")

# Save as Dictionary (safest for varying lengths)
pkl_path = f'{save_path}\\extracted_data_full.pkl'

try:
    with open(pkl_path, 'wb') as pfile:
        pickle.dump(data_neurons_dict, pfile)
    print(f"Successfully saved to: {pkl_path}")
    print(f"Keys: {list(data_neurons_dict.keys())}")
except TypeError as e:
    print(f"STILL FAILED: {e}") 
    # Debugging: check types in dict
    for k, v in data_neurons_dict.items():
        print(f"{k}: {type(v)}")
        if isinstance(v, (list, np.ndarray)) and len(v) > 0:
             print(f"  First element type: {type(v[0])}")

# %%
                                                    
# ██╗  ██╗███████╗██╗   ██╗    ███╗   ██╗███████╗██╗   ██╗██████╗  ██████╗ ███╗   ██╗███████╗    
# ██║ ██╔╝██╔════╝╚██╗ ██╔╝    ████╗  ██║██╔════╝██║   ██║██╔══██╗██╔═══██╗████╗  ██║██╔════╝    
# █████╔╝ █████╗   ╚████╔╝     ██╔██╗ ██║█████╗  ██║   ██║██████╔╝██║   ██║██╔██╗ ██║███████╗    
# ██╔═██╗ ██╔══╝    ╚██╔╝      ██║╚██╗██║██╔══╝  ██║   ██║██╔══██╗██║   ██║██║╚██╗██║╚════██║    
# ██║  ██╗███████╗   ██║       ██║ ╚████║███████╗╚██████╔╝██║  ██║╚██████╔╝██║ ╚████║███████║    
# ╚═╝  ╚═╝╚══════╝   ╚═╝       ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝    
                                                                                               
for key in data_neurons_dict:
    print(f"{key}: {type(data_neurons_dict[key])}")
    if isinstance(data_neurons_dict[key], dict):
        for subkey in data_neurons_dict[key]:
            print(f"  {subkey}: {type(data_neurons_dict[key][subkey])}")
data_neurons_dict = data_neurons_dict["eg_neurons"]  # Focus on the relevant sub-dictionary


# Apply to the cell_type column. Convert to Pandas Series first if it's an array.
data_neurons_dict['cell_type'] = pd.Series(data_neurons_dict['cell_type']).apply(decode_ascii_cell).values

# %%
# --- ENRICH METADATA WITH SPIKE TIMES (PkC_cs) ---
print("Enriching metadata with PkC_cs spike times...")

# 1. Identify Target Neurons
cell_types = data_neurons_dict['cell_type']
target_type = "PkC_cs"
# Indices of neurons that match the cell type
neuron_indices = [i for i, ct in enumerate(cell_types) if ct == target_type]
print(f"Found {len(neuron_indices)} neurons of type {target_type}")

# 2. Extract Spikes for Each Reach
# Ensure timestamps are available
timestamps = data_dict["session"].get("reconstructed_timestamps")

if timestamps is not None and len(neuron_indices) > 0:
    for meta in all_reaches_metadata:
        # Get reach window in seconds
        peak_idx = meta['peak_idx']
        onset_idx = meta.get('onset_idx', peak_idx - 20)
        
        # Define a window (e.g., onset - 0.5s to peak + 0.5s)
        # Check bounds
        t_start_idx = peak_idx - 150 # similar to plotting window
        t_end_idx = peak_idx + 50
        
        t_start = timestamps[t_start_idx]
        t_end = timestamps[t_end_idx]
        
        # Dictionary to store spikes for this reach: {neuron_idx: [spike_times]}
        meta['spikes'] = {}
        
        for n_idx in neuron_indices:
            # cell_spikes might be an array or list of arrays
            all_spikes = data_neurons_dict['st'][n_idx]
            
            # Helper to flatten if needed (sometimes loaded as [[t1], [t2]])
            if hasattr(all_spikes, 'flatten'):
                all_spikes = all_spikes.flatten()
            
            # Filter spikes in window
            # Assuming spikes are in seconds. If timestamps are in seconds, this matches.
            # If spikes were frames, we'd compare to indices. But usually 'st' = spike times (s).
            spikes_in_window = all_spikes[(all_spikes >= t_start) & (all_spikes <= t_end)]
            
            # Store relative to onset time? Or absolute? 
            # Storing absolute for now, can calculate relative later.
            meta['spikes'][n_idx] = spikes_in_window

    print("Spike times added to 'all_reaches_metadata'.")
else:
    print("Skipping spike enrichment: Missing timestamps or no neurons found.")
# %%
                          
#  ██████╗███████╗    ████████╗███████╗███╗   ██╗███████╗    
# ██╔════╝██╔════╝    ╚══██╔══╝██╔════╝████╗  ██║██╔════╝    
# ██║     ███████╗       ██║   ███████╗██╔██╗ ██║█████╗      
# ██║     ╚════██║       ██║   ╚════██║██║╚██╗██║██╔══╝      
# ╚██████╗███████║       ██║   ███████║██║ ╚████║███████╗    
#  ╚═════╝╚══════╝       ╚═╝   ╚══════╝╚═╝  ╚═══╝╚══════╝    


# --- 3. GMM CLUSTERING & INTERACTIVE PLOT WITH SPIKES ---
if 'tsne_result' in locals() and X_windows_3d is not None:
    print("Fitting GMM...")
    gmm = GaussianMixture(n_components=15, random_state=42)
    cluster_labels = gmm.fit_predict(tsne_result)
    
    class InteractiveGMMPlot_Spikes:
        def __init__(self, tsne_data, windows_3d, labels, metadata, timestamps):
            self.tsne = tsne_data
            self.windows = windows_3d
            self.labels = labels
            self.metadata = metadata
            self.timestamps = timestamps
            self.n_clusters = len(np.unique(labels))
            
            # Setup Main Figure (Scatter)
            self.fig_main, self.ax = plt.subplots(figsize=(10, 8))
            self.scatter = self.ax.scatter(self.tsne[:, 0], self.tsne[:, 1], 
                                           c=self.labels, cmap='nipy_spectral', 
                                           alpha=0.6, picker=5)
            self.fig_main.colorbar(self.scatter, label='Cluster ID')
            self.ax.set_title(f"t-SNE + GMM ({self.n_clusters} clusters)\nClick a point to see cluster average + Spikes")
            
            self.cid = self.fig_main.canvas.mpl_connect('pick_event', self.on_pick)
            
            # Setup Detail Figure (4 Rows: X, Y, Z, Spike Hist)
            self.fig_avg, self.axes_avg = plt.subplots(4, 1, figsize=(7, 12), sharex=True, 
                                                       gridspec_kw={'height_ratios': [1, 1, 1, 0.6]})
            self.fig_avg.suptitle("Cluster Averages")
            plt.show()

        def on_pick(self, event):
            if event.artist != self.scatter: return
            ind = event.ind[0]
            selected_cluster = self.labels[ind]
            self.update_average_plots(selected_cluster)

        def update_average_plots(self, cluster_id):
            indices = np.where(self.labels == cluster_id)[0]
            cluster_data = self.windows[indices] # (N, 3, Time)
            n_points = len(indices)

            mean_traj = np.mean(cluster_data, axis=0) # (3, Time)
            std_traj = np.std(cluster_data, axis=0)
            
            t_axis = np.arange(-150, 50)
            dims = ['X', 'Y', 'Z']
            colors = ['black', 'green', 'magenta']
            
            # --- Collect Spike Data ---
            all_spike_rel_frames = []
            spike_markers = {0:[], 1:[], 2:[]} # Lists of (frame, val) for plotting

            for i, global_idx in enumerate(indices):
                if global_idx < len(self.metadata):
                    m = self.metadata[global_idx]
                    peak_idx = m['peak_idx']
                    
                    if 'spikes' in m:
                        for neuron_id, spike_times in m['spikes'].items():
                            for st in spike_times:
                                # Conversion: Timestamp -> Frame Index
                                # We find the index in the global timestamp array
                                # np.searchsorted finds the insertion point to maintain order
                                spike_global_idx = np.searchsorted(self.timestamps, st)
                                
                                # Calculate relative frame (0 is peak)
                                rel_frame = spike_global_idx - peak_idx
                                
                                # Check if inside our window (-150 to +50)
                                if -150 <= rel_frame < 50:
                                    all_spike_rel_frames.append(rel_frame)
                                    
                                    # Get kinematic value at this specific frame for checking
                                    # array index 0 = -150 relative
                                    arr_idx = int(rel_frame + 150)
                                    if 0 <= arr_idx < 200:
                                        for dim in range(3):
                                            # Position of this specific trial at this frame
                                            pos = cluster_data[i, dim, arr_idx]
                                            spike_markers[dim].append((rel_frame, pos))

            # --- Plot Kinematics (Rows 0-2) ---
            for i in range(3):
                ax = self.axes_avg[i]
                ax.clear()
                # Mean + Std
                ax.plot(t_axis, mean_traj[i], color=colors[i], lw=2)
                ax.fill_between(t_axis, mean_traj[i]-std_traj[i], mean_traj[i]+std_traj[i], 
                                color=colors[i], alpha=0.02)
                
                # Plot Red Circles for Spikes
                if len(spike_markers[i]) > 0:
                    frames, positions = zip(*spike_markers[i])
                    ax.scatter(frames, positions, color='red', s=15, alpha=0.2, zorder=5, label='CS Spike')
                                # Faint Individual Traces
                subset_idx = np.random.choice(len(cluster_data), min(len(cluster_data), 50), replace=False)
                
                for _idx in subset_idx:
                    ax.plot(t_axis, cluster_data[_idx, i, :], color=colors[i], lw=0.5, alpha=0.5)

                ax.axvline(0, color='k', linestyle='--', alpha=0.2)
                ax.set_ylabel(f'{dims[i]} Position')
                ax.grid(True, alpha=0.3)
                if i==0 and len(spike_markers[0]) > 0: 
                    ax.legend(loc='upper right', fontsize='x-small')
            
            # --- Plot Spike Histogram (Row 3) ---
            ax_hist = self.axes_avg[3]
            ax_hist.clear()
            
            if len(all_spike_rel_frames) > 0:
                sns.histplot(all_spike_rel_frames, ax=ax_hist, kde=False, color='red', bins=50, element="step", alpha=0.5)
                med_spike = np.median(all_spike_rel_frames)
                ax_hist.axvline(med_spike, color='red', linestyle='--', lw=2, label=f'Med: {med_spike:.0f}')
                ax_hist.legend(loc='upper right', fontsize='small')
            else:

                ax_hist.text(0.5, 0.5, "No Spikes Found", ha='center', transform=ax_hist.transAxes)

            ax_hist.set_ylabel('Spike Count')
            ax_hist.set_xlabel('Frames from Peak (0)')
            ax_hist.set_xlim(-150, 50)
            ax_hist.grid(True, alpha=0.3)

            self.fig_avg.suptitle(f"Cluster {cluster_id} (n={n_points})")
            self.fig_avg.canvas.draw()

    # Pass timestamps for conversion
    ts_ref = data_dict['session']['reconstructed_timestamps']
    app = InteractiveGMMPlot_Spikes(tsne_result, X_windows_3d, cluster_labels, all_reaches_metadata, ts_ref)
else:
    print("Missing data required for clustering.")
# %%

#  ██████╗███████╗    ███╗   ██╗███████╗██╗   ██╗██████╗  ██████╗ ███╗   ██╗
# ██╔════╝██╔════╝    ████╗  ██║██╔════╝██║   ██║██╔══██╗██╔═══██╗████╗  ██║
# ██║     ███████╗    ██╔██╗ ██║█████╗  ██║   ██║██████╔╝██║   ██║██╔██╗ ██║
# ██║     ╚════██║    ██║╚██╗██║██╔══╝  ██║   ██║██╔══██╗██║   ██║██║╚██╗██║
# ╚██████╗███████║    ██║ ╚████║███████╗╚██████╔╝██║  ██║╚██████╔╝██║ ╚████║
#  ╚═════╝╚══════╝    ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝

if 'tsne_result' in locals() and X_windows_3d is not None:
    print("Fitting GMM...")
    gmm = GaussianMixture(n_components=15, random_state=42)
    cluster_labels = gmm.fit_predict(tsne_result)
    
    from matplotlib.widgets import TextBox

    class InteractiveGMMPlot_NeuronSelect:
        def __init__(self, tsne_data, windows_3d, labels, metadata, timestamps):
            self.tsne = tsne_data
            self.windows = windows_3d
            self.labels = labels
            self.metadata = metadata
            self.timestamps = timestamps
            self.n_clusters = len(np.unique(labels))
            
            # Identify all available neurons across metadata
            self.available_neurons = set()
            for m in self.metadata:
                if 'spikes' in m:
                    self.available_neurons.update(m['spikes'].keys())
            self.available_neurons = sorted(list(self.available_neurons))
            self.current_neuron = self.available_neurons[0] if self.available_neurons else None
            print(f"Available Neurons: {self.available_neurons}")
            
            # Track selection state
            self.selected_cluster = 0

            # Setup Main Figure (Scatter)
            self.fig_main, self.ax = plt.subplots(figsize=(10, 8))
            plt.subplots_adjust(bottom=0.2) # Make room for widget
            
            self.scatter = self.ax.scatter(self.tsne[:, 0], self.tsne[:, 1], 
                                           c=self.labels, cmap='nipy_spectral', 
                                           alpha=0.6, picker=5)
            self.fig_main.colorbar(self.scatter, label='Cluster ID')
            self.ax.set_title(f"t-SNE + GMM ({self.n_clusters} clusters)\nClick a point to select cluster")
            
            # Text Box for Neuron Selection
            axbox = plt.axes([0.4, 0.05, 0.2, 0.075])
            self.text_box = TextBox(axbox, 'Neuron ID:', initial=str(self.current_neuron))
            self.text_box.on_submit(self.submit_neuron)

            self.cid = self.fig_main.canvas.mpl_connect('pick_event', self.on_pick)
            
            # Setup Detail Figure (4 Rows: X, Y, Z, Spike Hist)
            self.fig_avg, self.axes_avg = plt.subplots(4, 1, figsize=(7, 12), sharex=True, 
                                                       gridspec_kw={'height_ratios': [1, 1, 1, 0.6]})
            self.fig_avg.suptitle("Cluster Averages")
            plt.show()

        def submit_neuron(self, text):
            try:
                n_id = int(text)
                # Check directly or allow partial match if needed, but strict is safer
                if n_id in self.available_neurons:
                    self.current_neuron = n_id
                    print(f"Selected Neuron: {n_id}")
                    # Update plots with new neuron for current cluster
                    self.update_average_plots(self.selected_cluster)
                else:
                    print(f"Neuron {n_id} not found in data.")
            except ValueError:
                print("Invalid ID")

        def on_pick(self, event):
            if event.artist != self.scatter: return
            ind = event.ind[0]
            # Update selected cluster
            self.selected_cluster = self.labels[ind]
            self.update_average_plots(self.selected_cluster)

        def update_average_plots(self, cluster_id):
            indices = np.where(self.labels == cluster_id)[0]
            if len(indices) == 0: return

            cluster_data = self.windows[indices] # (N, 3, Time)
            n_points = len(indices)

            mean_traj = np.mean(cluster_data, axis=0) # (3, Time)
            std_traj = np.std(cluster_data, axis=0)
            
            t_axis = np.arange(-150, 50)
            dims = ['X', 'Y', 'Z']
            colors = ['black', 'green', 'magenta']
            
            # --- Collect Spike Data for SELECTED NEURON ---
            all_spike_rel_frames = []
            spike_markers = {0:[], 1:[], 2:[]} # Lists of (frame, val) for plotting
            
            # Statistics Counters
            n_reaches_with_spikes = 0
            total_spikes_found = 0

            for i, global_idx in enumerate(indices):
                if global_idx < len(self.metadata):
                    m = self.metadata[global_idx]
                    peak_idx = m['peak_idx']
                    Peak_time = self.timestamps[peak_idx] if peak_idx < len(self.timestamps) else float('inf')

                    has_spike_in_window = False

                    if 'spikes' in m and self.current_neuron in m['spikes']:
                        spike_times = m['spikes'][self.current_neuron]
                        
                        # Handle potential single value or array
                        if isinstance(spike_times, (float, np.floating)) or not hasattr(spike_times, '__iter__'):
                             if not np.isnan(spike_times): spike_times = [spike_times]
                             else: spike_times = []
                        
                        for st in spike_times:
                            # 1. Find closest frame index in global timestamp array
                            spike_global_idx = np.searchsorted(self.timestamps, st)
                            
                            # Calculate relative frame (0 is peak)
                            rel_frame = spike_global_idx - peak_idx
                            
                            # Check window (-150 to +50)
                            if -150 <= rel_frame < 50:
                                all_spike_rel_frames.append(rel_frame)
                                has_spike_in_window = True
                                
                                # Get specific trial position at that frame
                                arr_idx = int(rel_frame + 150)
                                if 0 <= arr_idx < 200:
                                    for dim in range(3):
                                        pos = cluster_data[i, dim, arr_idx]
                                        spike_markers[dim].append((rel_frame, pos))
                    
                    if has_spike_in_window:
                        n_reaches_with_spikes += 1

            total_spikes_found = len(all_spike_rel_frames)
            proportion = (n_reaches_with_spikes / n_points) * 100 if n_points > 0 else 0

            # --- Plot Kinematics (Rows 0-2) ---
            for i in range(3):
                ax = self.axes_avg[i]
                ax.clear()
                # Mean + Std
                ax.plot(t_axis, mean_traj[i], color=colors[i], lw=2)
                ax.fill_between(t_axis, mean_traj[i]-std_traj[i], mean_traj[i]+std_traj[i], 
                                color=colors[i], alpha=0.02)
                
                # Plot Red Circles for Spikes (Specific Neuron)
                if len(spike_markers[i]) > 0:
                    frames, positions = zip(*spike_markers[i])
                    ax.scatter(frames, positions, color='red', s=20, alpha=0.8, zorder=5, label=f'N{self.current_neuron}')

                ax.axvline(0, color='k', linestyle='--', alpha=0.5)
                ax.set_ylabel(f'{dims[i]} Position')
                ax.grid(True, alpha=0.3)
                if i==0 and len(spike_markers[0]) > 0: 
                    ax.legend(loc='upper right', fontsize='x-small')
            
            # --- Plot Spike Histogram (Row 3) ---
            ax_hist = self.axes_avg[3]
            ax_hist.clear()
            
            title_stats = f"Count: {total_spikes_found} | Reaches with Spikes: {n_reaches_with_spikes}/{n_points} ({proportion:.1f}%)"
            
            if len(all_spike_rel_frames) > 0:
                sns.histplot(all_spike_rel_frames, ax=ax_hist, kde=False, color='red', bins=50, element="step", alpha=0.6)
                med_spike = np.median(all_spike_rel_frames)
                ax_hist.axvline(med_spike, color='red', linestyle='--', lw=2, label=f'Med: {med_spike:.0f}')
                ax_hist.legend(loc='upper right', fontsize='small')
            else:
                ax_hist.text(0.5, 0.5, f"No Spikes for Neuron {self.current_neuron}", ha='center', transform=ax_hist.transAxes)

            ax_hist.set_title(title_stats, fontsize=10, color='darkred')
            ax_hist.set_ylabel('Spike Count')
            ax_hist.set_xlabel('Frames from Peak (0)')
            ax_hist.set_xlim(-150, 50)
            ax_hist.grid(True, alpha=0.3)

            self.fig_avg.suptitle(f"Cluster {cluster_id} - Neuron {self.current_neuron}")
            self.fig_avg.canvas.draw()

    # Launch App
    ts_ref = data_dict['session']['reconstructed_timestamps']
    app = InteractiveGMMPlot_NeuronSelect(tsne_result, X_windows_3d, cluster_labels, all_reaches_metadata, ts_ref)
else:
    print("Missing data required for clustering.")
# %%
