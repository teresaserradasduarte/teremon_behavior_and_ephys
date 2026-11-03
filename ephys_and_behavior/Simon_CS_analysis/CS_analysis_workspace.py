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

%matplotlib qt

# CONFIGURATION
file_path = r"C:\Users\simon\Desktop\Data\ephys_and_behavior"  # Your file path
save_path = r"C:\Users\simon\Desktop\Data\ephys_and_behavior\1_CoteDor\R5"
mat_file = f'{file_path}\\st_CB_DCN.mat'
target_mouse = "1_CoteDor"
target_sess = "R5"
output_file = f'{file_path}\\filtered_data.parquet'


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
    
def process_huge_mat_file(mat_file,target_mouse,target_sess,output_file):
    print(f"Opening {mat_file}...")
    
    with h5py.File(mat_file, 'r') as f:
        # 1. INSPECT STRUCTURE
        # MATLAB variables are usually at the top level. 
        # If your data is in a struct, e.g., 'my_data.mouse', access it like f['my_data']['mouse']
        # Let's assume columns are top-level variables or inside a main struct key.
        
        # Uncomment this line to see the keys in your file:
        print("Keys found:", list(f.keys()))
        
        # NOTE: Adjust 'data_struct' below to the actual key name of your struct if it exists.
        # If the columns are at the root, use `group = f`
        group = f["all_CB_sts"]  # or f['data_struct_name']
        
        # 2. LOAD FILTER COLUMNS ONLY
        print("Loading filter columns...")
        print("Keys found in all_CB_neurons:", list(group.keys()))
        # MATLAB saves arrays as (N, 1). We flatten to (N,) for easy handling.
        # Warning: Dereferencing strings in a loop can be slow. 
        # If 'mouse' and 'sess' are simple numeric codes, this is instant.
        # If they are strings, we need to decode them carefully.
        
        # This part depends heavily on how MATLAB saved the strings.
        # Assuming they are cell arrays of strings:
        mouse_refs = group['mouse'][:] 
        sess_refs = group['sess'][:]
        

        if mouse_refs.shape[0] == 1 and mouse_refs.shape[1] > 1:
            mouse_refs = mouse_refs[0]
            sess_refs = sess_refs[0]
        print(f"Total rows to check: {len(mouse_refs)}")

        # 3. FIND MATCHING INDICES
        valid_indices = []
        
        # Iterate and decode strings to check condition
        # (We iterate because loading 67GB of strings into RAM might crash it)
        for i in range(len(mouse_refs)):
            # Decode mouse string
            ref = mouse_refs[i] # Get reference from cell array
            mouse_val = "".join(chr(c) for c in f[ref][:].flatten())
            
            if mouse_val == target_mouse:
                # Decode sess string only if mouse matches (optimization)
                ref_sess = sess_refs[i]
                sess_val = "".join(chr(c) for c in f[ref_sess][:].flatten())
                
                if sess_val == target_sess:
                    valid_indices.append(i)

        valid_indices = np.array(valid_indices)
        print(f"Found {len(valid_indices)} matching rows.")

        if len(valid_indices) == 0:
            print("No matching data found.")
            return

# ... (Previous code: Setup and Matching Indices is fine) ...

        # --- 2. EXTRACT & RESOLVE ---
        print("Extracting and resolving columns...")
        data_dict = {}
        sorted_idx = np.sort(valid_indices)
        
        # DEFINITION: Which columns MUST be flattened to 1D?
        # Add any other columns that are simple lists of numbers/times here.
        FLATTEN_COLS = ['st_init', 'st_reach', 'DCnD_init', 'idx_reach_hit', 'idx_reach_succ']

        for col_name in group.keys():
            if col_name in ['mouse', 'sess']: continue
            
            ds = group[col_name]
            if not isinstance(ds, h5py.Dataset): continue

            try:
                # Get Subset
                # Force 1D slice if it's a simple vector (N,1) or (1,N)
                if ds.ndim == 2 and ds.shape[0] == 1:
                    raw_subset = ds[0, sorted_idx]
                elif ds.ndim == 2 and ds.shape[1] == 1:
                    raw_subset = ds[sorted_idx, 0]
                else:
                    raw_subset = ds[sorted_idx]

                cleaned_column = []
                for item in raw_subset:
                    # 1. Recursive Resolve (Dig out the data)
                    val = recursive_resolve(f, item)
                    
                    # 2. Handle Numpy Arrays -> Lists
                    if isinstance(val, np.ndarray):
                        # SPECIAL HANDLING: FLATTEN SPECIFIC COLUMNS
                        if col_name in FLATTEN_COLS:
                            # Force 1D: turn [[0.1, 0.2]] into [0.1, 0.2]
                            val = val.flatten().tolist()
                        else:
                            # Keep structure (e.g. for reach_px)
                            val = val.tolist()
                            
                    elif val is None:
                        val = []
                        
                    # 3. Handle Scalars (Single spike case)
                    # If a neuron has 1 spike, 'val' might be a float, not a list.
                    # We must wrap it in a list so the column type is consistent.
                    if col_name in FLATTEN_COLS and not isinstance(val, list):
                        val = [val]
                        
                    cleaned_column.append(val)
                    
                data_dict[col_name] = cleaned_column

            except Exception as e:
                print(f"Skipping {col_name}: {e}")

        # --- 3. FINAL SAVE ---
        print("Saving to Parquet...")
        df = pd.DataFrame(data_dict)
        
        try:
            # We explicitly tell PyArrow to handle the schema inference
            df.to_parquet(output_file, engine='pyarrow', index=False)
            print("Success!")
        except Exception as e:
            print(f"Parquet save failed: {e}")
            print("Attempting to save as Pickle (slower, but handles complex types)...")
            df.to_pickle(output_file.replace('.parquet', '.pkl'))
            print("Saved as .pkl")
# %%
#process_huge_mat_file(mat_file, target_mouse, target_sess, output_file)
# %%

# ██████╗  █████╗ ████████╗ █████╗ 
# ██╔══██╗██╔══██╗╚══██╔══╝██╔══██╗
# ██║  ██║███████║   ██║   ███████║
# ██║  ██║██╔══██║   ██║   ██╔══██║
# ██████╔╝██║  ██║   ██║   ██║  ██║
# ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝
                                 
parquet_data = r"C:\Users\simon\Desktop\Data\ephys_and_behavior\filtered_data.pkl"
df_neurons = pd.read_pickle(parquet_data)
print(df_neurons.columns)

def decode_ascii_cell(nested_list):
    try:
        # Flatten [[80], [107], [67]] into [80, 107, 67]
        flat_list = [item[0] for item in nested_list]
        # Convert to characters and join into a string
        return "".join(chr(int(c)) for c in flat_list if c != 0)
    except Exception:
        return "Unknown"

# Apply to the cell_type column
df_neurons['cell_type'] = df_neurons['cell_type'].apply(decode_ascii_cell)
#%%
df_pkc = df_neurons[df_neurons['cell_type'] == "PkC_cs"].copy()
df_pkc = df_pkc.reset_index(drop=True)
print(f"Original count: {len(df_neurons)}")
print(f"Filtered count: {len(df_pkc)}")
# %%
print(df_pkc.columns)
len(df_pkc["st_reach"][0])
# %%
# accessin y coordinate of first neuron df_neurons["reach_px"][0][0][1]
# accessin the firing rate of the first neuron during the first reach df_neurons["FR_reach"][0][0]


# ███████╗██╗   ██╗███╗   ███╗███╗   ███╗ █████╗ ██████╗ ██╗   ██╗
# ██╔════╝██║   ██║████╗ ████║████╗ ████║██╔══██╗██╔══██╗╚██╗ ██╔╝
# ███████╗██║   ██║██╔████╔██║██╔████╔██║███████║██████╔╝ ╚████╔╝ 
# ╚════██║██║   ██║██║╚██╔╝██║██║╚██╔╝██║██╔══██║██╔══██╗  ╚██╔╝  
# ███████║╚██████╔╝██║ ╚═╝ ██║██║ ╚═╝ ██║██║  ██║██║  ██║   ██║   
# ╚══════╝ ╚═════╝ ╚═╝     ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝   
                                               
                                                                
# % Reach Category
# % 1 -> balistic from resting bar
# % 2 -> lifted paw
# % 3 -> grooming / groomung-like
# % 4 -> lifted paw: drink/reach
# --- 1. EXTRACT SESSION DATA ---
# We use the first neuron's row because reach metadata is identical for the whole session
row0 = df_pkc.iloc[0]
r_cat = row0['idx_reach_cat']
r_succ = row0['idx_reach_succ']
r_lcr = row0['idx_reach_LCR']

# Definitions
cat_labels = {
    1: "Ballistic (Resting)",
    2: "Lifted Paw",
    3: "Grooming-like",
    4: "Lifted: Drink/Reach"
}
loc_labels = {1: "Left", 2: "Center", 3: "Right"}

# --- 2. TALLY COUNTS ---
# Stats 1: Category Only
stats_cat = {k: {'Success': 0, 'Fail': 0} for k in cat_labels}

# Stats 2: Category + Location
# Key: (Category, Location)
stats_detail = {} 

grand_total = 0

if r_cat is not None:
    for cat, succ, lcr in zip(r_cat, r_succ, r_lcr):
        # Unwrap values (handle lists or scalars)
        c_val = cat[0] if isinstance(cat, list) else cat
        s_val = succ[0] if isinstance(succ, list) else succ
        l_val = lcr[0] if isinstance(lcr, list) else lcr
        
        # Only count known categories (1-4)
        if c_val in cat_labels:
            grand_total += 1
            
            # Update Table 1 Stats
            if s_val == 1:
                stats_cat[c_val]['Success'] += 1
            else:
                stats_cat[c_val]['Fail'] += 1
                
            # Update Table 2 Stats (if location is valid)
            if l_val in loc_labels:
                key = (c_val, l_val)
                if key not in stats_detail:
                    stats_detail[key] = {'Success': 0, 'Fail': 0}
                
                if s_val == 1:
                    stats_detail[key]['Success'] += 1
                else:
                    stats_detail[key]['Fail'] += 1

# --- 3. BUILD TABLE 1 (General Summary) ---
data_cat = []
for cat_id, counts in stats_cat.items():
    n_succ = counts['Success']
    n_fail = counts['Fail']
    n_total = n_succ + n_fail
    
    succ_rate = (n_succ / n_total) if n_total > 0 else 0.0
    prevalence = (n_total / grand_total) if grand_total > 0 else 0.0
    
    data_cat.append({
        "Category": f"{cat_id}: {cat_labels[cat_id]}",
        "Total": n_total,
        "Success": n_succ,
        "Fail": n_fail,
        "Succ Rate": f"{succ_rate:.1%}",
        "% of All": f"{prevalence:.1%}"
    })

df_summary_cat = pd.DataFrame(data_cat)

# --- 4. BUILD TABLE 2 (Detailed Breakdown) ---
data_detail = []
for c_code in range(1, 5):
    for l_code in range(1, 4):
        key = (c_code, l_code)
        if key in stats_detail:
            counts = stats_detail[key]
            n_succ = counts['Success']
            n_fail = counts['Fail']
            n_total = n_succ + n_fail
            
            succ_rate = (n_succ / n_total) if n_total > 0 else 0.0
            prevalence = (n_total / grand_total) if grand_total > 0 else 0.0
            
            data_detail.append({
                "Category": cat_labels[c_code],
                "Location": loc_labels[l_code],
                "Total": n_total,
                "Success": n_succ,
                "Fail": n_fail,
                "Succ Rate": f"{succ_rate:.1%}",
                "% of All": f"{prevalence:.1%}"
            })

df_summary_detail = pd.DataFrame(data_detail)

# --- 5. DISPLAY ---
print(f"Total Reaches Analyzed: {grand_total}")
print("\n" + "="*40)
print("TABLE 1: SUMMARY BY CATEGORY (Aggregated)")
print("="*40)
print(df_summary_cat.to_string(index=False))

print("\n" + "="*40)
print("TABLE 2: DETAILED BREAKDOWN (Category + Location)")
print("="*40)
# Use to_string to ensure all rows print
print(df_summary_detail.to_string(index=False))

output_lines = []
output_lines.append(f"REACH STATISTICS SUMMARY")
output_lines.append(f"Total Reaches Analyzed: {grand_total}")
output_lines.append("\n" + "="*50)
output_lines.append("TABLE 1: SUMMARY BY CATEGORY (Aggregated)")
output_lines.append("="*50)
output_lines.append(df_summary_cat.to_string(index=False))
output_lines.append("\n" + "="*50)
output_lines.append("TABLE 2: DETAILED BREAKDOWN (Category + Location)")
output_lines.append("="*50)
output_lines.append(df_summary_detail.to_string(index=False))

# Join into single string
full_report = "\n".join(output_lines)

# --- 5. PRINT AND SAVE ---
print(full_report)

# Save to file
filename = file_path + "\\reach_summary_stats.txt"
with open(filename, "w") as f:
    f.write(full_report)

print(f"\n[INFO] Successfully saved stats to '{filename}'")


#%%
# init_valid PP  or invPP 
plt.plot(df_pkc["reach_px"][0][0][1])
# %%
# CONFIGURATION
# --- CONFIGURATION ---

# ██████╗ ██╗███████╗████████╗     ██████╗███████╗     ██████╗ █████╗ ████████╗     ██╗
# ██╔══██╗██║██╔════╝╚══██╔══╝    ██╔════╝██╔════╝    ██╔════╝██╔══██╗╚══██╔══╝    ███║
# ██║  ██║██║███████╗   ██║       ██║     ███████╗    ██║     ███████║   ██║       ╚██║
# ██║  ██║██║╚════██║   ██║       ██║     ╚════██║    ██║     ██╔══██║   ██║        ██║
# ██████╔╝██║███████║   ██║       ╚██████╗███████║    ╚██████╗██║  ██║   ██║        ██║
# ╚═════╝ ╚═╝╚══════╝   ╚═╝        ╚═════╝╚══════╝     ╚═════╝╚═╝  ╚═╝   ╚═╝        ╚═╝
                                                                                     
FPS = 400
CENTER_FRAME = 150

# Convert Frames to Seconds
CENTER_SEC = CENTER_FRAME / FPS  # 0.375s
WINDOW_SEC = [-150 / FPS, 50 / FPS]  # [-0.375s, +0.125s]

TARGET_CAT = 1
TARGET_SUCC = 1 


reaches_succ = df_pkc.iloc[0]['idx_reach_succ']
reaches_cat = df_pkc.iloc[0]['idx_reach_cat']

total_matching_reaches = 0

# 2. Iterate through the reaches of this SINGLE session
if reaches_succ is not None and len(reaches_succ) > 0:
    for succ, cat in zip(reaches_succ, reaches_cat):
        
        # Check Success (Handles 1 vs [1])
        is_succ = (succ == 1) or (isinstance(succ, list) and succ[0] == 1)
        
        # Check Category (Handles 1 vs [1])
        is_cat = (cat == 1) or (isinstance(cat, list) and cat[0] == 1)
        
        if is_succ and is_cat:
            total_matching_reaches += 1

print(f"Total Unique Reaches (Succ | Cat 1): {total_matching_reaches}")

all_aligned_spikes = []

# Iterate through each neuron
for idx, row in df_pkc.iterrows():
    
    reaches_spikes = row['st_reach']
    reaches_succ = row['idx_reach_succ']
    reaches_cat = row['idx_reach_cat']
    
    if reaches_spikes is None or len(reaches_spikes) == 0:
        continue
        
    for spikes, succ, cat in zip(reaches_spikes, reaches_succ, reaches_cat):
        
        # --- CHECK CONDITIONS ---
        # Handle formats: 1, [1], or 1.0
        is_succ = (succ == TARGET_SUCC) or (isinstance(succ, list) and succ[0] == TARGET_SUCC)
        is_cat = (cat == TARGET_CAT) or (isinstance(cat, list) and cat[0] == TARGET_CAT)
        
        if is_succ and is_cat:
            if spikes is None or len(spikes) == 0:
                continue

            spikes_arr = np.array(spikes)
            
            # --- ALIGNMENT (SECONDS) ---
            # 1. Align: Subtract the time of the event (0.375s)
            # This makes the event happen at t=0
            aligned_spikes = spikes_arr - CENTER_SEC
            
            # 2. Filter: Keep spikes within the relative window [-0.375, +0.125]
            mask = (aligned_spikes >= WINDOW_SEC[0]) & (aligned_spikes <= WINDOW_SEC[1])
            
            valid_spikes = aligned_spikes[mask]
            
            all_aligned_spikes.extend(valid_spikes)

all_aligned_spikes = np.array(all_aligned_spikes)
print(f"Collected {len(all_aligned_spikes)} aligned spikes.")

# --- PLOTTING ---
plt.figure(figsize=(10, 6))

# Bin width of 10ms (0.01s)
bin_width = 0.01
bins = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + bin_width, bin_width)

sns.histplot(all_aligned_spikes, bins=bins, kde=False, color='teal', label='Spikes')

plt.axvline(0, color='red', linestyle='--', label='Reach Center (Frame 150)')
plt.xlabel('Time relative to Reach Center (s)')
plt.ylabel('Spike Count')
plt.title(f'PSTH: Successful Reaches (Cat {TARGET_CAT} - N = {total_matching_reaches})\nNeuron: PkC_cs')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

#%% 

# ██████╗ ██╗███████╗████████╗     ██████╗███████╗     ██████╗ ██╗   ██╗███████╗██████╗      ██████╗ █████╗ ████████╗
# ██╔══██╗██║██╔════╝╚══██╔══╝    ██╔════╝██╔════╝    ██╔═══██╗██║   ██║██╔════╝██╔══██╗    ██╔════╝██╔══██╗╚══██╔══╝
# ██║  ██║██║███████╗   ██║       ██║     ███████╗    ██║   ██║██║   ██║█████╗  ██████╔╝    ██║     ███████║   ██║   
# ██║  ██║██║╚════██║   ██║       ██║     ╚════██║    ██║   ██║╚██╗ ██╔╝██╔══╝  ██╔══██╗    ██║     ██╔══██║   ██║   
# ██████╔╝██║███████║   ██║       ╚██████╗███████║    ╚██████╔╝ ╚████╔╝ ███████╗██║  ██║    ╚██████╗██║  ██║   ██║   
# ╚═════╝ ╚═╝╚══════╝   ╚═╝        ╚═════╝╚══════╝     ╚═════╝   ╚═══╝  ╚══════╝╚═╝  ╚═╝     ╚═════╝╚═╝  ╚═╝   ╚═╝   


# --- CONFIGURATION ---
FPS = 400
CENTER_FRAME = 150
CENTER_SEC = CENTER_FRAME / FPS  # 0.375s
WINDOW_SEC = [-1,0.5]  # [-0.375s, +0.125s]
MIN_REACHES = 5  # Threshold to plot
# --- CONFIGURATION ---

# Binning Setup (10ms bins)
bin_width = 0.01
bins = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + bin_width, bin_width)

cat_labels = {
    1: "Ballistic (Resting)",
    2: "Lifted Paw",
    3: "Grooming-like",
    4: "Lifted: Drink/Reach"
}
num_neurons = len(df_pkc)

# --- PLOTTING SETUP ---
# 4 Rows (Categories), 2 Columns (Success/Fail)
fig, axes = plt.subplots(4, 2, figsize=(12, 16), sharex=True, sharey='row')

# Column Titles
axes[0, 0].set_title("SUCCESSFUL", fontsize=14, fontweight='bold', color='cornflowerblue')
axes[0, 1].set_title("FAILED", fontsize=14, fontweight='bold', color='indianred')

print(f"Processing 8-Panel Plot (Min Reaches = {MIN_REACHES})...")
# --- MASTER LOOP ---
for cat_idx, target_cat in enumerate([1, 2, 3, 4]):
    
    for col_idx, target_succ in enumerate([1, 0]):
        
        ax = axes[cat_idx, col_idx]
        color = 'cornflowerblue' if target_succ == 1 else 'indianred'
        
        # --- 1. COUNT REACHES ---
        row0 = df_pkc.iloc[0]
        r_succ = row0['idx_reach_succ']
        r_cat = row0['idx_reach_cat']
        
        reach_count = 0
        if r_succ is not None:
            for s, c in zip(r_succ, r_cat):
                s_val = s[0] if isinstance(s, list) else s
                c_val = c[0] if isinstance(c, list) else c
                if s_val == target_succ and c_val == target_cat:
                    reach_count += 1
        
        # --- 2. COLLECT DATA ---
        all_aligned_spikes = []
        
        if reach_count >= MIN_REACHES:
            for _, row in df_pkc.iterrows():
                n_spikes = row['st_reach']
                n_succ = row['idx_reach_succ']
                n_cat = row['idx_reach_cat']
                
                if n_spikes is None: continue
                
                for spikes, s, c in zip(n_spikes, n_succ, n_cat):
                    s_val = s[0] if isinstance(s, list) else s
                    c_val = c[0] if isinstance(c, list) else c
                    
                    if s_val == target_succ and c_val == target_cat:
                        if spikes is not None:
                            # Robust Flattening
                            flat_reach = np.array(spikes).ravel().tolist()
                            for t in flat_reach:
                                aligned_t = t - CENTER_SEC
                                if WINDOW_SEC[0] <= aligned_t <= WINDOW_SEC[1]:
                                    all_aligned_spikes.append(aligned_t)

            # Plot Histogram
            if len(all_aligned_spikes) > 0:
                data_to_plot = np.array(all_aligned_spikes)
                total_seconds = reach_count * num_neurons * bin_width
                weights = np.full(len(data_to_plot), 1 / total_seconds)
                
                ax.hist(data_to_plot, bins=bins, weights=weights, 
                        color=color, alpha=0.7, edgecolor='white', linewidth=0.5)
                
                # N label top-right
                ax.text(0.95, 0.95, f"N={reach_count}", 
                        transform=ax.transAxes, ha='right', va='top', fontsize=10, 
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9))
            else:
                ax.text(0.5, 0.5, f"N={reach_count}\n(Silent)", 
                        transform=ax.transAxes, ha='center', va='center', color='gray')

        else:
            # --- HANDLING INSUFFICIENT DATA ---
            # transform=ax.transAxes ensures (0.5, 0.5) is always the visual center of the box
            ax.text(0.5, 0.5, f"Insufficient Data\n(N={reach_count})", 
                    transform=ax.transAxes, ha='center', va='center', 
                    color='gray', fontstyle='italic', fontsize=12)

        # --- 3. UNIVERSAL FORMATTING ---
        # Apply limits even if empty to maintain the "box" size
        ax.set_xlim(WINDOW_SEC)
        ax.set_ylim([0,6])

        ax.axvline(0, color='black', linestyle='--', alpha=0.5)
        ax.grid(True, alpha=0.2)
        
        # Row Labels
        if col_idx == 0:
            ax.set_ylabel(f"{cat_labels[target_cat]}\nFiring Rate (Hz)", fontsize=11, fontweight='bold')
        
        # X-Axis Labels
        if cat_idx == 3:
            ax.set_xlabel("Time relative to Reach (s)")

plt.suptitle(f"Population CS Activity (N={num_neurons} cells)", fontsize=16, y=0.99)
plt.tight_layout()
plt.show()
plt.savefig(save_path + "\\population_reach_activity.png", dpi=300)
# %%

#  ██████╗███████╗     █████╗ ███╗   ██╗██████╗     ██████╗ ███████╗ █████╗  ██████╗██╗  ██╗
# ██╔════╝██╔════╝    ██╔══██╗████╗  ██║██╔══██╗    ██╔══██╗██╔════╝██╔══██╗██╔════╝██║  ██║
# ██║     ███████╗    ███████║██╔██╗ ██║██║  ██║    ██████╔╝█████╗  ███████║██║     ███████║
# ██║     ╚════██║    ██╔══██║██║╚██╗██║██║  ██║    ██╔══██╗██╔══╝  ██╔══██║██║     ██╔══██║
# ╚██████╗███████║    ██║  ██║██║ ╚████║██████╔╝    ██║  ██║███████╗██║  ██║╚██████╗██║  ██║
#  ╚═════╝╚══════╝    ╚═╝  ╚═╝╚═╝  ╚═══╝╚═════╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝
                  
# --- CONFIGURATION ---
FPS = 400
CENTER_FRAME = 150
CENTER_SEC = CENTER_FRAME / FPS

# Trajectory Time Axis (Assuming 300 frames per reach)
FRAMES_PER_REACH = 201
time_axis = (np.arange(FRAMES_PER_REACH) - CENTER_FRAME) / FPS

# Spike Window
WINDOW_SEC = [-150/FPS, 50/FPS] # [-0.375s, +0.125s]

TARGET_CAT = 1
TARGET_SUCC = 1

# Mapping: Value -> (Name, Spikes, Trajectories, Count)
loc_data = {
    1: {'name': 'Left',   'spikes': [], 'trajectories': [], 'reach_count': 0},
    2: {'name': 'Center', 'spikes': [], 'trajectories': [], 'reach_count': 0},
    3: {'name': 'Right',  'spikes': [], 'trajectories': [], 'reach_count': 0}
}

# --- STEP 1: COLLECT TRAJECTORIES & COUNT (Session Level) ---
# We use row0 because trajectory info is the same for all neurons
row0 = df_pkc.iloc[0]
r_succ = row0['idx_reach_succ']
r_cat = row0['idx_reach_cat']
r_lcr = row0['idx_reach_LCR']
r_px = row0['reach_px']  # <--- Added Trajectory Column

if r_succ is not None:
    for succ, cat, lcr, px in zip(r_succ, r_cat, r_lcr, r_px):
        # Unwrap values
        s_val = succ[0] if isinstance(succ, list) else succ
        c_val = cat[0] if isinstance(cat, list) else cat
        l_val = lcr[0] if isinstance(lcr, list) else lcr
        
        # Check Filters
        if s_val == TARGET_SUCC and c_val == TARGET_CAT:
            if l_val in loc_data:
                # 1. Increment Count
                loc_data[l_val]['reach_count'] += 1
                
                # 2. Store Trajectory (X-dimension only)
                # px is [x, y, z], so we take index 0
                x_traj = np.array(px[0], dtype=float).flatten()
                
                # Safety check: ensure it matches our expected frame count
                if len(x_traj) == FRAMES_PER_REACH:
                    loc_data[l_val]['trajectories'].append(x_traj)

# --- STEP 2: COLLECT SPIKES (Neuron Level) ---
# (This part is exactly your previous logic)
for idx, row in df_pkc.iterrows():
    r_spikes = row['st_reach']
    r_succ = row['idx_reach_succ']
    r_cat = row['idx_reach_cat']
    r_lcr = row['idx_reach_LCR']
    
    if r_spikes is None or len(r_spikes) == 0:
        continue

    for spikes, succ, cat, lcr in zip(r_spikes, r_succ, r_cat, r_lcr):
        s_val = succ[0] if isinstance(succ, list) else succ
        c_val = cat[0] if isinstance(cat, list) else cat
        l_val = lcr[0] if isinstance(lcr, list) else lcr
        
        if s_val == TARGET_SUCC and c_val == TARGET_CAT and l_val in loc_data:
            if spikes is not None and len(spikes) > 0:
                spikes_arr = np.array(spikes)
                aligned = spikes_arr - CENTER_SEC
                mask = (aligned >= WINDOW_SEC[0]) & (aligned <= WINDOW_SEC[1])
                loc_data[l_val]['spikes'].extend(aligned[mask])

# --- PLOTTING (2 ROWS x 3 COLS) ---
# --- PLOTTING (NORMALIZED) ---
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)

# 1. Define Normalization Constants
num_neurons = len(df_pkc)  # Total number of cells recorded
bin_width = 0.01 
bins = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + bin_width, bin_width)
colors = {1: '#4c72b0', 2: '#dd8452', 3: '#55a868'}

for col_idx, lcr_code in enumerate([1, 2, 3]):
    data = loc_data[lcr_code]
    count = data['reach_count']
    name = data['name']
    
    # Select subplots
    ax_top = axes[0, col_idx]
    ax_bot = axes[1, col_idx]
    
    # -----------------------------
    # TOP ROW: TRAJECTORIES (Same as before)
    # -----------------------------
    trajs = np.array(data['trajectories'])
    if len(trajs) > 0:
        for traj in trajs:
            ax_top.plot(time_axis, traj, color=colors[lcr_code], alpha=0.3)
        mean_traj = np.mean(trajs, axis=0)
        std_traj = np.std(trajs, axis=0)
        ax_top.plot(time_axis, mean_traj, color=colors[lcr_code], lw=3)
        ax_top.fill_between(time_axis, mean_traj - std_traj, mean_traj + std_traj, 
                            color=colors[lcr_code], alpha=0.3)
    
    ax_top.set_title(f"{name}\n(n={count} reaches)", fontsize=14, fontweight='bold')
    ax_top.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax_top.grid(True, alpha=0.2)
    if col_idx == 0: 
        ax_top.set_ylabel('X Position (px)')
        ax_top.set_ylim(np.min(trajs)-10, np.max(trajs)+10)
    # -----------------------------
    # BOTTOM ROW: PROBABILITY HISTOGRAM
    # -----------------------------
    spikes_list = data['spikes']
    
    if len(spikes_list) > 0 and count > 0:
# --- NUCLEAR FLATTENER ---
        # We manually iterate to handle mixed types (scalars, arrays, nested arrays)
        clean_spikes = []
        for item in spikes_list:
            # 1. Force every item to be a numpy array
            # 2. Reshape(-1) forces it to be 1D, even if it was scalar or 2D matrix
            try:
                clean_spikes.extend(np.array(item).reshape(-1))
            except Exception:
                continue # Skip corrupt entries
                
        # Convert the final clean list to a Float Array
        spikes = np.array(clean_spikes, dtype=float)
        
        # --- CALCULATE WEIGHTS ---
        # Only proceed if we actually extracted spikes
        if len(spikes) > 0:
            norm_factor = count*num_neurons*bin_width#count * num_neurons * bin_width
            # Create weights array matching the EXACT length of spikes
            weights = np.full(len(spikes), 1 / norm_factor)
            
            # --- PLOT ---
            # We explicitly pass x and weights
            ax_bot.hist(spikes, bins=bins, weights=weights, 
                        color=colors[lcr_code], alpha=0.8)
    
    ax_bot.axvline(0, color='red', linestyle='--', alpha=0.5)
    ax_bot.grid(True, alpha=0.2)
    ax_bot.set_xlabel('Time (s)')
    
    # Update Y-Label to reflect probability
    if col_idx == 0: 
        ax_bot.set_ylabel('Firing Rate (Hz)')

# Final Formatting
plt.xlim(WINDOW_SEC[0], WINDOW_SEC[1])
plt.suptitle(f"Kinematics vs CS Probability (Succ | Cat {TARGET_CAT})", fontsize=16, y=0.98)
plt.tight_layout()
plt.show()

#%%

# ███████╗██╗███╗   ██╗ ██████╗ ██╗     ███████╗    ███╗   ██╗███████╗██╗   ██╗██████╗  ██████╗ ███╗   ██╗    
# ██╔════╝██║████╗  ██║██╔════╝ ██║     ██╔════╝    ████╗  ██║██╔════╝██║   ██║██╔══██╗██╔═══██╗████╗  ██║    
# ███████╗██║██╔██╗ ██║██║  ███╗██║     █████╗      ██╔██╗ ██║█████╗  ██║   ██║██████╔╝██║   ██║██╔██╗ ██║    
# ╚════██║██║██║╚██╗██║██║   ██║██║     ██╔══╝      ██║╚██╗██║██╔══╝  ██║   ██║██╔══██╗██║   ██║██║╚██╗██║    
# ███████║██║██║ ╚████║╚██████╔╝███████╗███████╗    ██║ ╚████║███████╗╚██████╔╝██║  ██║╚██████╔╝██║ ╚████║    
# ╚══════╝╚═╝╚═╝  ╚═══╝ ╚═════╝ ╚══════╝╚══════╝    ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝    
                                                                                                            
                                                                                                            
                                                                                                            
# --- CONFIGURATION ---
FPS = 400
CENTER_FRAME = 150
CENTER_SEC = CENTER_FRAME / FPS  # 0.375s
WINDOW_SEC = [-150/FPS, 50/FPS]  # [-0.375s, +0.125s]

# Exponential Decay Kernel Settings
TAU_MS = 20  # Time constant in ms
KERNEL_SIZE_MS = 100 # How long the tail lasts in ms

# Visualization Settings
COLORS = {1: '#4c72b0', 2: '#dd8452', 3: '#55a868'} # Blue, Orange, Green
LOC_NAMES = {1: 'Left', 2: 'Center', 3: 'Right'}

def generate_kernel(fps, tau_ms, kernel_size_ms):
    dt = 1000 / fps
    size = int(kernel_size_ms / dt)
    t = np.arange(size) * dt
    return np.exp(-t / tau_ms)

def inspect_neuron_sorted(df, neuron_idx, target_cat=1, target_succ=1):
    
    # 1. GET NEURON DATA
    row = df.iloc[neuron_idx]
    r_spikes = row['st_reach']
    r_succ = row['idx_reach_succ']
    r_cat = row['idx_reach_cat']
    r_lcr = row['idx_reach_LCR'] # We need location for sorting
    
    # 2. COLLECT AND TAG DATA
    # We store tuples: (spikes_array, location_code)
    collected_reaches = []
    
    if r_spikes is not None:
        for spikes, s, c, lcr in zip(r_spikes, r_succ, r_cat, r_lcr):
            # Unwrap
            s_val = s[0] if isinstance(s, list) else s
            c_val = c[0] if isinstance(c, list) else c
            l_val = lcr[0] if isinstance(lcr, list) else lcr
            
            if s_val == target_succ and c_val == target_cat:
                if spikes is not None:
                    # Robust flatten & align
                    flat_spikes = np.array(spikes).ravel().tolist()
                    aligned = np.array(flat_spikes) - CENTER_SEC
                    
                    # Store only valid spikes in window
                    valid_mask = (aligned >= WINDOW_SEC[0]) & (aligned <= WINDOW_SEC[1])
                    collected_reaches.append({
                        'spikes': aligned[valid_mask],
                        'loc': l_val
                    })

    n_reaches = len(collected_reaches)
    print(f"Neuron {neuron_idx}: Found {n_reaches} matching reaches.")
    
    if n_reaches == 0:
        print("No reaches found.")
        return

    # 3. SORT DATA (Left -> Center -> Right)
    # This groups the waterfall plot by location
    collected_reaches.sort(key=lambda x: x['loc'])

    # --- PLOTTING ---
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(10, 12), sharex=True, 
                                         gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05})
    
    dt_sec = 1 / FPS
    t_axis = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + dt_sec, dt_sec)
    kernel = generate_kernel(FPS, TAU_MS, KERNEL_SIZE_MS)
    y_offset_step = 1.2
    
    # Storage for bottom histograms
    spikes_by_loc = {1: [], 2: [], 3: []}
    counts_by_loc = {1: 0, 2: 0, 3: 0}

    # --------------------------
    # TOP PANEL: WATERFALL
    # --------------------------
    for i, trial_data in enumerate(collected_reaches):
        spikes = trial_data['spikes']
        loc = trial_data['loc']
        
        # Aggregate for bottom plot
        if loc in spikes_by_loc:
            spikes_by_loc[loc].extend(spikes)
            counts_by_loc[loc] += 1
        
        # 1. Convolve Trace
        counts, _ = np.histogram(spikes, bins=np.append(t_axis, t_axis[-1] + dt_sec))
        trace = convolve(counts, kernel, mode='full')[:len(t_axis)]
        
        # 2. Plot Trace
        y_base = i * y_offset_step
        ax_top.plot(t_axis, trace + y_base, color='black', lw=0.8, alpha=0.6)
        
        # 3. Red Dots
        spike_indices = np.where(counts > 0)[0]
        if len(spike_indices) > 0:
            ax_top.scatter(t_axis[spike_indices], trace[spike_indices] + y_base, 
                           color='red', s=10, zorder=5)
            
        # 4. COLOR BAR (The visual strip on the left)
        # We draw a rectangle at the very left edge of the plot window
        # x position: slightly shifted left of start time
        bar_width = (WINDOW_SEC[1] - WINDOW_SEC[0]) * 0.02 # 2% of width
        bar_x = WINDOW_SEC[0] - bar_width
        
        # Draw a colored rectangle for this specific trial row
        rect = patches.Rectangle((WINDOW_SEC[0], y_base), width=bar_width, height=y_offset_step, 
                                 facecolor=COLORS.get(loc, 'gray'), edgecolor=None, clip_on=False)
        ax_top.add_patch(rect)

    ax_top.set_title(f"Neuron {neuron_idx} | Cat {target_cat} | {'Success' if target_succ else 'Fail'}\n(Exponential Conv. $\\tau={TAU_MS}ms$)", fontsize=14)
    ax_top.set_ylabel(f"Trials # {n_reaches}", fontsize=12)
    ax_top.set_yticks([])
    ax_top.spines['left'].set_visible(False)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    
    # --------------------------
    # BOTTOM PANEL: OVERLAID PSTH
    # --------------------------
    bin_width = 0.01
    bins = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + bin_width, bin_width)
    
    # Iterate 1, 2, 3 so the legend order is consistent
    for loc_code in [1, 2, 3]:
        data = np.array(spikes_by_loc[loc_code])
        count = counts_by_loc[loc_code]
        
        if count > 0 and len(data) > 0:
            # Normalization (Hz per location)
            weights = np.full(len(data), 1 / (count))
            
            ax_bot.hist(data, bins=bins, weights=weights, 
                        color=COLORS[loc_code], alpha=0.3, 
                        label=f"{LOC_NAMES[loc_code]} (N={count})",
                        edgecolor=COLORS[loc_code], linewidth=1)
            
            # Add a smoothed line for clarity (optional, but helps with overlapping histograms)
            # You can comment this out if you strictly want bars only
            h_counts, h_edges = np.histogram(data, bins=bins, weights=weights)
            centers = (h_edges[:-1] + h_edges[1:]) / 2
            smooth = gaussian_filter1d(h_counts, sigma=1.5)
            ax_bot.plot(centers, smooth, color=COLORS[loc_code], lw=2)

    ax_bot.set_ylabel("Firing Rate (Hz)", fontsize=12)
    ax_bot.set_xlabel("Time relative to Reach (s)", fontsize=12)
    ax_bot.axvline(0, color='black', linestyle='--', alpha=0.5)
    ax_bot.legend(loc='upper right', frameon=False)
    ax_bot.grid(True, alpha=0.2)
    ax_bot.set_xlim(WINDOW_SEC) # Ensure alignment with top

    plt.tight_layout()
    plt.show()

# --- EXECUTE ---
inspect_neuron_sorted(df_pkc, neuron_idx=0, target_cat=1, target_succ=1)          


#%% 

# ██╗      ██████╗  ██████╗ ██████╗     ████████╗██╗  ██╗██████╗  ██████╗ ██╗   ██╗ ██████╗ ██╗  ██╗    ███╗   ██╗███████╗██╗   ██╗██████╗  ██████╗ ███╗   ██╗███████╗
# ██║     ██╔═══██╗██╔═══██╗██╔══██╗    ╚══██╔══╝██║  ██║██╔══██╗██╔═══██╗██║   ██║██╔════╝ ██║  ██║    ████╗  ██║██╔════╝██║   ██║██╔══██╗██╔═══██╗████╗  ██║██╔════╝
# ██║     ██║   ██║██║   ██║██████╔╝       ██║   ███████║██████╔╝██║   ██║██║   ██║██║  ███╗███████║    ██╔██╗ ██║█████╗  ██║   ██║██████╔╝██║   ██║██╔██╗ ██║███████╗
# ██║     ██║   ██║██║   ██║██╔═══╝        ██║   ██╔══██║██╔══██╗██║   ██║██║   ██║██║   ██║██╔══██║    ██║╚██╗██║██╔══╝  ██║   ██║██╔══██╗██║   ██║██║╚██╗██║╚════██║
# ███████╗╚██████╔╝╚██████╔╝██║            ██║   ██║  ██║██║  ██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║    ██║ ╚████║███████╗╚██████╔╝██║  ██║╚██████╔╝██║ ╚████║███████║
# ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝            ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝    ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝
                                                                                                                                                                    
# --- HELPER: Plot Logic for ONE Column (Top + Bottom) ---
def plot_single_neuron(ax_top, ax_bot, row, neuron_idx, target_cat, target_succ):
    
    # 1. EXTRACT DATA
    r_spikes = row['st_reach']
    r_succ = row['idx_reach_succ']
    r_cat = row['idx_reach_cat']
    r_lcr = row['idx_reach_LCR']
    
    collected_reaches = []
    
    if r_spikes is not None:
        for spikes, s, c, lcr in zip(r_spikes, r_succ, r_cat, r_lcr):
            s_val = s[0] if isinstance(s, list) else s
            c_val = c[0] if isinstance(c, list) else c
            l_val = lcr[0] if isinstance(lcr, list) else lcr
            
            if s_val == target_succ and c_val == target_cat:
                if spikes is not None:
                    flat_spikes = np.array(spikes).ravel().tolist()
                    aligned = np.array(flat_spikes) - CENTER_SEC
                    valid_mask = (aligned >= WINDOW_SEC[0]) & (aligned <= WINDOW_SEC[1])
                    collected_reaches.append({
                        'spikes': aligned[valid_mask],
                        'loc': l_val
                    })

    n_reaches = len(collected_reaches)
    
    # Handle Empty Neuron
    if n_reaches == 0:
        ax_top.text(0.5, 0.5, "No Data", ha='center', va='center')
        ax_bot.text(0.5, 0.5, "No Data", ha='center', va='center')
        return

    # 2. SORT (L -> C -> R)
    collected_reaches.sort(key=lambda x: x['loc'])

    # 3. PREPARE PLOTTING
    dt_sec = 1 / FPS
    t_axis = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + dt_sec, dt_sec)
    kernel = generate_kernel(FPS, TAU_MS, KERNEL_SIZE_MS)
    y_offset_step = 1.2
    
    spikes_by_loc = {1: [], 2: [], 3: []}
    counts_by_loc = {1: 0, 2: 0, 3: 0}

    # 4. PLOT TOP (Waterfall)
    for i, trial_data in enumerate(collected_reaches):
        spikes = trial_data['spikes']
        loc = trial_data['loc']
        
        # Aggregate
        if loc in spikes_by_loc:
            spikes_by_loc[loc].extend(spikes)
            counts_by_loc[loc] += 1
        
        # Convolve
        counts, _ = np.histogram(spikes, bins=np.append(t_axis, t_axis[-1] + dt_sec))
        trace = convolve(counts, kernel, mode='full')[:len(t_axis)]
        
        # Plot
        y_base = i * y_offset_step
        ax_top.plot(t_axis, trace + y_base, color='black', lw=0.5, alpha=0.6)
        
        # Red dots
        spike_indices = np.where(counts > 0)[0]
        if len(spike_indices) > 0:
            ax_top.scatter(t_axis[spike_indices], trace[spike_indices] + y_base, 
                           color='red', s=5, zorder=5) # Smaller size for dense plots
            
        # Color Bar
        bar_width = (WINDOW_SEC[1] - WINDOW_SEC[0]) * 0.04
        rect = patches.Rectangle((WINDOW_SEC[0], y_base), width=bar_width, height=y_offset_step, 
                                 facecolor=COLORS.get(loc, 'gray'), clip_on=False)
        ax_top.add_patch(rect)

    ax_top.set_title(f"Neuron {neuron_idx}", fontsize=10, fontweight='bold')
    ax_top.set_yticks([])
    ax_top.spines['left'].set_visible(False)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.set_xlim(WINDOW_SEC)

    # 5. PLOT BOTTOM (Histograms)
    bin_width = 0.01
    bins = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + bin_width, bin_width)
    
    for loc_code in [1, 2, 3]:
        data = np.array(spikes_by_loc[loc_code])
        count = counts_by_loc[loc_code]
        
        if count > 0 and len(data) > 0:
            # Metric: As requested (1/count)
            weights = np.full(len(data), 1 / (count))
        
            
            ax_bot.hist(data, bins=bins, weights=weights, 
                        color=COLORS[loc_code], alpha=0.3, 
                        label=f"{LOC_NAMES[loc_code]} (N={count})",
                        edgecolor=COLORS[loc_code], linewidth=1)
            
            # Smooth line
            h_counts, h_edges = np.histogram(data, bins=bins, weights=weights)
            centers = (h_edges[:-1] + h_edges[1:]) / 2
            smooth = gaussian_filter1d(h_counts, sigma=1.5)
            ax_bot.plot(centers, smooth, color=COLORS[loc_code], lw=1.5)

    ax_bot.axvline(0, color='black', linestyle='--', alpha=0.3)
    ax_bot.grid(True, alpha=0.2)
    ax_bot.set_xlim(WINDOW_SEC)
    ax_bot.legend(loc='upper right', frameon=False)


# --- MAIN EXECUTION LOOP ---
def plot_all_neurons_batched(df, target_cat=1, target_succ=1):
    
    BATCH_SIZE = 8
    total_neurons = len(df)
    
    # Loop in chunks of 5
    for start_idx in range(0, total_neurons, BATCH_SIZE):
        end_idx = min(start_idx + BATCH_SIZE, total_neurons)
        
        # Create 10-column figure
        # 2 Rows: Top (Waterfall) and Bottom (Hist)
        # Height ratio 3:1
        fig, axes = plt.subplots(2, BATCH_SIZE, figsize=(20, 8), 
                                 sharex='col', 
                                 gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.1, 'wspace': 0.3})
        
        print(f"Plotting Neurons {start_idx} to {end_idx - 1}...")
        
        # Iterate through the columns (0 to 4)
        for col_idx in range(BATCH_SIZE):
            neuron_idx = start_idx + col_idx
            
            ax_top = axes[0, col_idx]
            ax_bot = axes[1, col_idx]
            
            if neuron_idx < total_neurons:
                # Plot the neuron
                row = df.iloc[neuron_idx]
                plot_single_neuron(ax_top, ax_bot, row, neuron_idx, target_cat, target_succ)
                
                # Labels
                if col_idx == 0:
                    ax_top.set_ylabel("Trials", fontsize=10)
                    ax_bot.set_ylabel("Rate", fontsize=10)
                
                ax_bot.set_xlabel("Time (s)", fontsize=9)
                
            else:
                # Turn off empty subplots
                ax_top.axis('off')
                ax_bot.axis('off')
        
        plt.suptitle(f"Cat {target_cat} | Succ {target_succ} | Neurons {start_idx}-{end_idx-1}", y=0.98, fontsize=14)
        plt.show()
        plt.savefig(save_path + f"\\Neuron_Batch_{start_idx}_to_{end_idx-1}_Cat{target_cat}_Succ{target_succ}.png", dpi=300)

# Run the Plotter
plot_all_neurons_batched(df_pkc, target_cat=1, target_succ=0)
# %%
# ██████╗  ██████╗ █████╗ 
# ██╔══██╗██╔════╝██╔══██╗
# ██████╔╝██║     ███████║
# ██╔═══╝ ██║     ██╔══██║
# ██║     ╚██████╗██║  ██║
# ╚═╝      ╚═════╝╚═╝  ╚═╝
                                                                                              
# --- CONFIGURATION ---
FPS = 400
CENTER_FRAME = 150
# Time axis for plotting
time_axis = (np.arange(201) - CENTER_FRAME) / FPS

TARGET_CAT = 1
TARGET_SUCC = 1

# Storage: 1=Left, 2=Center, 3=Right
# We now store the full (3, 300) matrices instead of just X
loc_data = {
    1: {'name': 'Left',   '3d_trajs': []},
    2: {'name': 'Center', '3d_trajs': []},
    3: {'name': 'Right',  '3d_trajs': []}
}

# --- STEP 1: COLLECT 3D DATA ---
row0 = df_pkc.iloc[0]
r_px = row0['reach_px']      # List of [x, y, z]
r_succ = row0['idx_reach_succ']
r_cat = row0['idx_reach_cat']
r_lcr = row0['idx_reach_LCR']

if r_px is not None:
    for px, succ, cat, lcr in zip(r_px, r_succ, r_cat, r_lcr):
        # Unwrap
        s_val = succ[0] if isinstance(succ, list) else succ
        c_val = cat[0] if isinstance(cat, list) else cat
        l_val = lcr[0] if isinstance(lcr, list) else lcr
        
        if s_val == TARGET_SUCC and c_val == TARGET_CAT:
            if l_val in loc_data:
                # px is [x_array, y_array, z_array]
                # We define shape (300, 3) for PCA
                # 1. Stack to (3, 300)
                coords = np.vstack([px[0], px[1], px[2]])
                # 2. Transpose to (300, 3) -> (Samples, Features)
                coords_t = coords.T 
                
                # Check length (assume 300 frames)
                if coords_t.shape[0] == 201:
                    loc_data[l_val]['3d_trajs'].append(coords_t)

# --- STEP 2: COMPUTE PCAs ---

# A. Prepare Global Data (Stack ALL reaches from ALL locations)
all_trajs_list = []
for k in loc_data:
    all_trajs_list.extend(loc_data[k]['3d_trajs'])

# Stack into one giant matrix (N_total_frames, 3)
# This lets PCA find the global "average" movement axis
global_points = np.vstack(all_trajs_list)

# Fit Global PCA
pca_global = PCA(n_components=1) # We only need the main axis (PC1)
pca_global.fit(global_points)

# --- STEP 3: PLOTTING ---
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
colors = {1: '#4c72b0', 2: '#dd8452', 3: '#55a868'}

for col_idx, lcr_code in enumerate([1, 2, 3]):
    
    data = loc_data[lcr_code]
    name = data['name']
    trajs_3d = data['3d_trajs'] # List of (201, 3) arrays
    
    if len(trajs_3d) == 0: continue
    
    # ---------------------------------------
    # ROW 0: GLOBAL PCA (Same Basis for all)
    # ---------------------------------------
    ax_global = axes[0, col_idx]
    
    pc1_traces_global = []
    for traj in trajs_3d:
        # Project this reach onto the Global PC1
        # transform returns (201, 1), we flatten to (201,)
        proj = pca_global.transform(traj).flatten()
        pc1_traces_global.append(proj)
        
        # Plot faint individual line
        ax_global.plot(time_axis, proj, color=colors[lcr_code], alpha=0.05, lw=1)

    # Plot Mean Global
    pc1_matrix_g = np.array(pc1_traces_global)
    mean_g = np.mean(pc1_matrix_g, axis=0)
    std_g = np.std(pc1_matrix_g, axis=0)
    
    ax_global.plot(time_axis, mean_g, color='black', lw=2, linestyle='--', label='Mean')
    # Standard deviation shading
    ax_global.fill_between(time_axis, mean_g - std_g, mean_g + std_g, 
                           color='black', alpha=0.1)

    ax_global.set_title(f"{name} (Global PCA)", fontsize=12, fontweight='bold')
    ax_global.axvline(0, color='red', linestyle=':')
    ax_global.grid(True, alpha=0.2)
    if col_idx == 0: ax_global.set_ylabel('Position along Global PC1')

    # ---------------------------------------
    # ROW 1: LOCAL PCA (Specific Basis)
    # ---------------------------------------
    ax_local = axes[1, col_idx]
    
    # Fit PCA ONLY on this location's data
    local_points = np.vstack(trajs_3d)
    pca_local = PCA(n_components=1)
    pca_local.fit(local_points)
    
    pc1_traces_local = []
    for traj in trajs_3d:
        # Project onto Local PC1
        proj = pca_local.transform(traj).flatten()
        
        # Flip sign check: 
        # PCA direction is arbitrary (Left vs Right). 
        # We align them so they generally go "up" (positive) for consistency.
        if np.mean(proj) < 0: proj = -proj
            
        pc1_traces_local.append(proj)
        ax_local.plot(time_axis, proj, color=colors[lcr_code], alpha=0.05, lw=1)
        
    # Plot Mean Local
    pc1_matrix_l = np.array(pc1_traces_local)
    mean_l = np.mean(pc1_matrix_l, axis=0)
    std_l = np.std(pc1_matrix_l, axis=0) # Add std calculation here

    ax_local.plot(time_axis, mean_l, color=colors[lcr_code], lw=3, label='Mean')
    ax_local.fill_between(time_axis, mean_l - std_l, mean_l + std_l, 
                           color=colors[lcr_code], alpha=0.2)

    ax_local.set_title(f"{name} (Local PCA)", fontsize=12, fontweight='bold')
    ax_local.axvline(0, color='red', linestyle=':')
    ax_local.grid(True, alpha=0.2)
    ax_local.set_xlabel('Time (s)')
    if col_idx == 0: ax_local.set_ylabel('Position along Local PC1')

plt.suptitle(f"3D Trajectory Compressed via PCA (Global vs Local)", fontsize=16, y=0.95)
plt.tight_layout()
plt.show()
# %%

#  █████╗ ██╗     ██╗     ██╗ ██████╗ ███╗   ██╗    ████████╗ ██████╗     ██╗███╗   ██╗██╗████████╗██╗ █████╗ ████████╗██╗ ██████╗ ███╗   ██╗    ██████╗ ███████╗ █████╗  ██████╗██╗  ██╗
# ██╔══██╗██║     ██║     ██║██╔════╝ ████╗  ██║    ╚══██╔══╝██╔═══██╗    ██║████╗  ██║██║╚══██╔══╝██║██╔══██╗╚══██╔══╝██║██╔═══██╗████╗  ██║    ██╔══██╗██╔════╝██╔══██╗██╔════╝██║  ██║
# ███████║██║     ██║     ██║██║  ███╗██╔██╗ ██║       ██║   ██║   ██║    ██║██╔██╗ ██║██║   ██║   ██║███████║   ██║   ██║██║   ██║██╔██╗ ██║    ██████╔╝█████╗  ███████║██║     ███████║
# ██╔══██║██║     ██║     ██║██║   ██║██║╚██╗██║       ██║   ██║   ██║    ██║██║╚██╗██║██║   ██║   ██║██╔══██║   ██║   ██║██║   ██║██║╚██╗██║    ██╔══██╗██╔══╝  ██╔══██║██║     ██╔══██║
# ██║  ██║███████╗███████╗██║╚██████╔╝██║ ╚████║       ██║   ╚██████╔╝    ██║██║ ╚████║██║   ██║   ██║██║  ██║   ██║   ██║╚██████╔╝██║ ╚████║    ██║  ██║███████╗██║  ██║╚██████╗██║  ██║
# ╚═╝  ╚═╝╚══════╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝       ╚═╝    ╚═════╝     ╚═╝╚═╝  ╚═══╝╚═╝   ╚═╝   ╚═╝╚═╝  ╚═╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝    ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝


# --- CONFIGURATION ---
FPS = 400
SMOOTH_WINDOW = 21 
POLY_ORDER = 3
ONSET_THRESHOLD_PCT = 0.10 # 10% of peak speed


def get_onset_frame(trajectory_3d, fps):
    """
    Detects movement onset using Iterative Valley Search.
    1. Finds all valleys before the peak.
    2. Walks backwards through valleys to find one < 10% of peak.
    3. If none found, returns the valley closest to the peak.
    """
    # 1. Calculate Velocity & Speed
    vel = np.diff(trajectory_3d, axis=0) * fps
    speed = np.linalg.norm(vel, axis=1)
    
    # 2. Smooth
    if len(speed) > SMOOTH_WINDOW:
        speed_smooth = speed#savgol_filter(speed, window_length=SMOOTH_WINDOW, polyorder=POLY_ORDER)
    else:
        speed_smooth = speed 
    
    # 3. Find Peak (Restrict search to first 250 frames)

    roi = speed_smooth[10:150]
    if len(roi) == 0: return None
    
    peak_idx = np.argmax(roi)
    peak_val = roi[peak_idx]
    
    # Define the "Deep Valley" threshold
    threshold = peak_val * ONSET_THRESHOLD_PCT
    
    # 4. Find ALL Valleys (Local Minima) before the peak
    # We slice speed array up to the peak index
    pre_peak_speed = speed_smooth[0:peak_idx]
    
    # argrelextrema returns indices of local minima
    valley_indices = argrelextrema(pre_peak_speed, np.less, order=1)[0]
    
    onset_idx = None
    if len(valley_indices) > 0:
        # We have valleys. Let's search them from RIGHT to LEFT (Closest to peak -> Furthest)
        # valley_indices is sorted ascending (0, 10, 50...), so we reverse it.
        
        # Default Fallback: The valley closest to the peak (the last one in the list)
        fallback_onset = valley_indices[-1]
        
        found_deep_valley = False
        
        # Iterate backwards: [Valley_3, Valley_2, Valley_1]
        for v_idx in reversed(valley_indices):
            if speed_smooth[v_idx] < threshold:
                # FOUND IT: This valley is quiet enough to be a "true start"
                onset_idx = v_idx
                found_deep_valley = True
                break
        
        # If we checked them all and none were deep enough, use the fallback
        if not found_deep_valley:
            onset_idx = fallback_onset
            
    else:
        # NO VALLEYS FOUND (Monotonic rise from index 0)
        # The only logical options are index 0 or the Peak.
        # Usually index 0 is safer if it's a pure ramp-up.
        onset_idx = 0 

    return onset_idx, speed_smooth


def fit_global_pca(df):
    """
    Fits a PCA on ALL reaches in the session to find the 
    Global Main Axis of movement.
    """
    # 1. Collect all trajectories
    all_points = []
    
    # We only need row 0 since reach data is session-wide
    row0 = df.iloc[0]
    r_px = row0['reach_px']
    
    if r_px is not None:
        for px in r_px:
            # px is [x_array, y_array, z_array]
            # Transpose to (N_frames, 3)
            traj = np.array([px[0], px[1], px[2]]).T
            all_points.append(traj)
            
    # 2. Stack into one giant (Total_Frames, 3) matrix
    if len(all_points) > 0:
        global_data = np.vstack(all_points)
        
        # 3. Fit PCA
        pca = PCA(n_components=1)
        pca.fit(global_data)
        
        # OPTIONAL: Force polarity
        # We want positive PC1 to mean "Forward/Outward". 
        # We check the mean end position of reaches. 
        # If the average reach ends at a negative PC1 value, we flip the component.
        # (This is a heuristic, but helpful for visualization)
        test_proj = pca.transform(global_data[-100:]) # Check end of data
        if np.mean(test_proj) < 0:
            pca.components_ *= -1
            
        return pca
    else:
        return None
    

def get_onset_frame_global_pc1(trajectory_3d, fps, global_pca):
    """
    Detects onset using projection onto the GLOBAL PCA axis.
    """
    SMOOTH_WINDOW = 21 
    POLY_ORDER = 3
    ONSET_THRESHOLD_PCT = 0.10
    
    # 1. GLOBAL PROJECTION
    # Transform (N,3) -> (N,1) -> (N,) using the session-wide basis
    pc1_trace = global_pca.transform(trajectory_3d).flatten()
    
    # 2. PC1 SPEED (Absolute change along the main axis)
    speed_pc1 = np.abs(np.diff(pc1_trace) * fps)
    
    # 3. SMOOTHING
    if len(speed_pc1) > SMOOTH_WINDOW:
        speed_smooth = speed_pc1#savgol_filter(speed_pc1, window_length=SMOOTH_WINDOW, polyorder=POLY_ORDER)
    else:
        speed_smooth = speed_pc1 
    
    # 4. FIND PEAK (Search 0-250)

    roi = speed_smooth[10:150]
    if len(roi) == 0: return None
    
    peak_idx = np.argmax(roi)
    peak_val = roi[peak_idx]
    
    threshold = peak_val * ONSET_THRESHOLD_PCT
    
    # 5. VALLEY SEARCH (Same Logic as before)
    pre_peak_speed = speed_smooth[0:peak_idx]
    valley_indices = argrelextrema(pre_peak_speed, np.less, order=1)[0]
    
    onset_idx = None
    if len(valley_indices) > 0:
        fallback_onset = valley_indices[-1]
        found_deep_valley = False
        
        for v_idx in reversed(valley_indices):
            if speed_smooth[v_idx] < threshold:
                onset_idx = v_idx
                found_deep_valley = True
                break
        
        if not found_deep_valley:
            onset_idx = fallback_onset
    else:
        onset_idx = 0 

    return onset_idx,speed_smooth

# --- PLOTTING LOOP ---
row0 = df_pkc.iloc[0]
r_px = row0['reach_px'] # List of [x, y, z]
r_cat = row0['idx_reach_cat']

fig, axes = plt.subplots(4, 5, figsize=(20, 8)) # No sharey, scales differ
axes_flat = axes.flatten()

plotted_count = 0
target_cat = 2
global_pca = fit_global_pca(df_pkc)
if r_px is not None:
    for px, cat in zip(r_px, r_cat):
        c_val = cat[0] if isinstance(cat, list) else cat
        
        if c_val == target_cat and plotted_count < 20:
            
            # 1. Prepare Trajectory (300, 3)
            traj = np.array([px[0], px[1], px[2]]).T
            # 2. Compute Onset (Speed Based)
            onset_idx_pc, speed_prof_pc = get_onset_frame_global_pc1(traj, FPS, global_pca)
            onset_idx, speed_prof = get_onset_frame(traj, FPS)
            # 3. Compute PC1 (Geometry Based)
            pca = PCA(n_components=1)
            pc1 = pca.fit_transform(traj).flatten()
            
            # Align PC1 polarity: make it go "up" for easier visual comparison
            if np.mean(pc1[-50:]) < np.mean(pc1[:50]):
                pc1 = -pc1
            
            if onset_idx is not None:
                ax = axes_flat[plotted_count]
                
                # --- Plot Speed (Left Axis) ---
                ax.plot(speed_prof, color='black', alpha=0.8, label='Speed', lw=1.5)
                ax.set_ylabel("Speed (px/s)", color='black', fontsize=8)
                ax.tick_params(axis='y', labelcolor='black')
                
                # --- Plot PC1 (Right Axis) ---
                ax2 = ax.twinx()
                ax2.plot(pc1, color='cornflowerblue', alpha=0.8, label='PC1', lw=2)
                ax2.set_ylabel("PC1 Position", color='cornflowerblue', fontsize=8)
                ax2.tick_params(axis='y', labelcolor='cornflowerblue')
                
                # --- Markers ---
                # Onset (Red Dot)
                ax.plot(onset_idx, speed_prof[onset_idx], 'ro', markersize=5, zorder=10, label='Onset')
                ax.plot(onset_idx_pc, speed_prof_pc[onset_idx_pc], 'bo', markersize=5, zorder=10, label='Onset')

                # Peak Speed (Green X)
                peak_idx = np.argmax(speed_prof[10:150])
                ax.plot(peak_idx+10, speed_prof[peak_idx+10], 'gx', zorder=10)
                
                # Old Center (Gray Line)
                ax.axvline(150, color='gray', linestyle='--', alpha=0.5)
                ax.axvline(onset_idx, color='red', linestyle=':', alpha=0.5)

                ax.set_title(f"Reach {plotted_count+1} | Onset: {onset_idx}")
                plotted_count += 1

plt.suptitle(f"Onset Validation: Speed vs PC1 (Threshold={ONSET_THRESHOLD_PCT*100}%)", fontsize=16)
plt.tight_layout()
plt.show()
# %%

# ███╗   ██╗███████╗██╗   ██╗██████╗  ██████╗ ███╗   ██╗███████╗     █████╗ ██╗     ██╗ ██████╗ ███╗   ██╗███████╗██████╗     ██╗███╗   ██╗██╗████████╗
# ████╗  ██║██╔════╝██║   ██║██╔══██╗██╔═══██╗████╗  ██║██╔════╝    ██╔══██╗██║     ██║██╔════╝ ████╗  ██║██╔════╝██╔══██╗    ██║████╗  ██║██║╚══██╔══╝
# ██╔██╗ ██║█████╗  ██║   ██║██████╔╝██║   ██║██╔██╗ ██║███████╗    ███████║██║     ██║██║  ███╗██╔██╗ ██║█████╗  ██║  ██║    ██║██╔██╗ ██║██║   ██║   
# ██║╚██╗██║██╔══╝  ██║   ██║██╔══██╗██║   ██║██║╚██╗██║╚════██║    ██╔══██║██║     ██║██║   ██║██║╚██╗██║██╔══╝  ██║  ██║    ██║██║╚██╗██║██║   ██║   
# ██║ ╚████║███████╗╚██████╔╝██║  ██║╚██████╔╝██║ ╚████║███████║    ██║  ██║███████╗██║╚██████╔╝██║ ╚████║███████╗██████╔╝    ██║██║ ╚████║██║   ██║   
# ╚═╝  ╚═══╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝    ╚═╝  ╚═╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝╚═════╝     ╚═╝╚═╝  ╚═══╝╚═╝   ╚═╝   



WINDOW_SEC = [-0.125, +0.375]

# --- 2. SPEED-BASED DETECTION ---
def get_onset_idx_speed(trajectory_3d, fps):
    # 1. Calculate Euclidean Speed
    vel = np.diff(trajectory_3d, axis=0) * fps
    speed = np.linalg.norm(vel, axis=1) # No smoothing as requested
    
    # 2. Find Peak (10-150 range)
    search_start, search_end = 10, 150
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
                onset_idx = v_idx
                found_deep_valley = True
                break
        
        if not found_deep_valley:
            onset_idx = valley_indices[-1] # Fallback: closest valley
            
    return onset_idx

# --- 3. PC1-BASED DETECTION ---
def get_onset_idx_pc1(trajectory_3d, fps, global_pca):
    # 1. Project to Global PC1
    pc1_trace = global_pca.transform(trajectory_3d).flatten()
    
    # 2. Calculate PC1 Speed
    speed_pc1 = np.abs(np.diff(pc1_trace) * fps) # No smoothing
    
    # 3. Find Peak
    search_start, search_end = 10, 150
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
                onset_idx = v_idx
                found_deep_valley = True
                break
        if not found_deep_valley:
            onset_idx = valley_indices[-1]

    return onset_idx

# --- 4. ARBITRATION LOGIC (THE MASTER FUNCTION) ---
def get_arbitrated_onset(trajectory_3d, fps, global_pca):
    """
    Decides between Speed-based and PC1-based onset.
    Rule: Use Speed by default. 
          If diff > 5 frames, pick the one with lower PC1 position (more retracted).
    """
    # Get candidates
    idx_speed = get_onset_idx_speed(trajectory_3d, fps)
    idx_pc1 = get_onset_idx_pc1(trajectory_3d, fps, global_pca)
    
    # Check discrepancy
    if abs(idx_speed - idx_pc1) > 5:
        # Get PC1 Position trace (not speed)
        pc1_pos = global_pca.transform(trajectory_3d).flatten()
        
        # Safety check for indices out of bounds
        max_idx = len(pc1_pos) - 1
        val_speed = pc1_pos[min(idx_speed, max_idx)]
        val_pc1 = pc1_pos[min(idx_pc1, max_idx)]
        
        # Pick the index corresponding to the "Lower" (Start) position
        if val_pc1 < val_speed:
            return idx_pc1
        else:
            return idx_speed
    else:
        # Agreement (or close enough) -> Default to Speed
        return idx_speed
    

def plot_single_neuron_arbitrated(ax_top, ax_bot, row, neuron_idx, target_cat, target_succ, global_pca):
    
    r_spikes = row['st_reach']
    r_px = row['reach_px']
    r_succ = row['idx_reach_succ']
    r_cat = row['idx_reach_cat']
    r_lcr = row['idx_reach_LCR']
    
    collected_reaches = []
    
    if r_spikes and r_px and global_pca:
        for spikes, px, s, c, lcr in zip(r_spikes, r_px, r_succ, r_cat, r_lcr):
            s_val = s[0] if isinstance(s, list) else s
            c_val = c[0] if isinstance(c, list) else c
            l_val = lcr[0] if isinstance(lcr, list) else lcr
            
            if s_val == target_succ and c_val == target_cat:
                if spikes is not None and px is not None:
                    
                    # --- ARBITRATED ALIGNMENT ---
                    traj = np.array([px[0], px[1], px[2]]).T
                    onset_frame = get_arbitrated_onset(traj, FPS, global_pca)
                    
                    onset_sec = onset_frame / FPS
                    
                    # Shift Spikes
                    flat_spikes = np.array(spikes).ravel().tolist()
                    aligned = np.array(flat_spikes) - onset_sec
                    
                    valid_mask = (aligned >= WINDOW_SEC[0]) & (aligned <= WINDOW_SEC[1])
                    collected_reaches.append({
                        'spikes': aligned[valid_mask],
                        'loc': l_val
                    })

    n_reaches = len(collected_reaches)
    if n_reaches == 0:
        ax_top.text(0.5, 0.5, "No Data", ha='center', va='center')
        ax_bot.text(0.5, 0.5, "No Data", ha='center', va='center')
        return

    # Sort & Plot
    collected_reaches.sort(key=lambda x: x['loc'])
    
    dt_sec = 1 / FPS
    t_axis = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + dt_sec, dt_sec)
    kernel = generate_kernel(FPS, TAU_MS, KERNEL_SIZE_MS)
    y_offset_step = 1.2
    
    spikes_by_loc = {1: [], 2: [], 3: []}
    counts_by_loc = {1: 0, 2: 0, 3: 0}

    # Top Panel (Waterfall)
    for i, trial_data in enumerate(collected_reaches):
        spikes = trial_data['spikes']
        loc = trial_data['loc']
        
        if loc in spikes_by_loc:
            spikes_by_loc[loc].extend(spikes)
            counts_by_loc[loc] += 1
        
        counts, _ = np.histogram(spikes, bins=np.append(t_axis, t_axis[-1] + dt_sec))
        trace = convolve(counts, kernel, mode='full')[:len(t_axis)]
        
        y_base = i * y_offset_step
        ax_top.plot(t_axis, trace + y_base, color='black', lw=0.5, alpha=0.6)
        
        spike_indices = np.where(counts > 0)[0]
        if len(spike_indices) > 0:
            ax_top.scatter(t_axis[spike_indices], trace[spike_indices] + y_base, 
                           color='red', s=5, zorder=5)
            
        bar_width = (WINDOW_SEC[1] - WINDOW_SEC[0]) * 0.04
        rect = patches.Rectangle((WINDOW_SEC[0], y_base), width=bar_width, height=y_offset_step, 
                                 facecolor=COLORS.get(loc, 'gray'), clip_on=False)
        ax_top.add_patch(rect)

    ax_top.set_title(f"Neuron {neuron_idx}", fontsize=10, fontweight='bold')
    ax_top.axis('off')
    ax_top.set_xlim(WINDOW_SEC)

    # Bottom Panel (Histogram)
    bin_width = 0.01
    bins = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + bin_width, bin_width)
    for loc_code in [1, 2, 3]:
        data = np.array(spikes_by_loc[loc_code])
        count = counts_by_loc[loc_code]
        if count > 0 and len(data) > 0:
            weights = np.full(len(data), 1 / count)
            ax_bot.hist(data, bins=bins, weights=weights, color=COLORS[loc_code], alpha=0.3, edgecolor=None)
            h_counts, h_edges = np.histogram(data, bins=bins, weights=weights)
            centers = (h_edges[:-1] + h_edges[1:]) / 2
            smooth = gaussian_filter1d(h_counts, sigma=1.5)
            ax_bot.plot(centers, smooth, color=COLORS[loc_code], lw=1.5)

    ax_bot.axvline(0, color='red', linestyle='--', alpha=0.5)
    ax_bot.grid(True, alpha=0.2)
    ax_bot.set_xlim(WINDOW_SEC)

def plot_all_neurons_arbitrated(df, target_cat=1, target_succ=1):
    print("Fitting Global PCA...")
    global_pca = fit_global_pca(df)
    
    BATCH_SIZE = 5
    total_neurons = len(df)
    
    for start_idx in range(0, total_neurons, BATCH_SIZE):
        end_idx = min(start_idx + BATCH_SIZE, total_neurons)
        
        fig, axes = plt.subplots(2, 5, figsize=(20, 8), 
                                 sharex='col', 
                                 gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.1, 'wspace': 0.3})
        
        print(f"Plotting Neurons {start_idx}-{end_idx-1} (Arbitrated Onset)...")
        
        for col_idx in range(BATCH_SIZE):
            neuron_idx = start_idx + col_idx
            ax_top = axes[0, col_idx]
            ax_bot = axes[1, col_idx]
            
            if neuron_idx < total_neurons:
                row = df.iloc[neuron_idx]
                plot_single_neuron_arbitrated(ax_top, ax_bot, row, neuron_idx, target_cat, target_succ, global_pca)
                if col_idx == 0: ax_bot.set_ylabel("Rate (Hz)")
                ax_bot.set_xlabel("Time from Onset (s)")
            else:
                ax_top.axis('off')
                ax_bot.axis('off')
        
        succ_lbl = "Success" if target_succ else "Fail"
        plt.suptitle(f"Arbitrated Onset Aligned | Cat {target_cat} | {succ_lbl}", y=0.98, fontsize=14)
        plt.show()
        plt.savefig(save_path + f"\\Neuron_Batch_{start_idx}_to_{end_idx-1}_Cat{target_cat}_Succ{target_succ}_Aligned_Onset.png", dpi=300)

# Run
plot_all_neurons_batched(df_pkc, target_cat=1, target_succ=1)
# %%

# ██████╗ ██╗███████╗████████╗     █████╗ ██╗     ██╗ ██████╗ ███╗   ██╗███████╗██████╗      ██████╗ █████╗ ████████╗
# ██╔══██╗██║██╔════╝╚══██╔══╝    ██╔══██╗██║     ██║██╔════╝ ████╗  ██║██╔════╝██╔══██╗    ██╔════╝██╔══██╗╚══██╔══╝
# ██║  ██║██║███████╗   ██║       ███████║██║     ██║██║  ███╗██╔██╗ ██║█████╗  ██║  ██║    ██║     ███████║   ██║   
# ██║  ██║██║╚════██║   ██║       ██╔══██║██║     ██║██║   ██║██║╚██╗██║██╔══╝  ██║  ██║    ██║     ██╔══██║   ██║   
# ██████╔╝██║███████║   ██║       ██║  ██║███████╗██║╚██████╔╝██║ ╚████║███████╗██████╔╝    ╚██████╗██║  ██║   ██║   
# ╚═════╝ ╚═╝╚══════╝   ╚═╝       ╚═╝  ╚═╝╚══════╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝╚═════╝      ╚═════╝╚═╝  ╚═╝   ╚═╝   
                                                                                                                   
# --- CONFIGURATION ---
FPS = 400
# NEW WINDOW: Aligning to Onset [-0.2s before start, +0.5s after start]
WINDOW_SEC = [-0.5, 1] 
MIN_REACHES = 5
ONSET_THRESHOLD_PCT = 0.10

# Binning Setup (10ms bins)
bin_width = 0.01
bins = np.arange(WINDOW_SEC[0], WINDOW_SEC[1] + bin_width, bin_width)

print("Fitting Global PCA...")
global_pca = fit_global_pca(df_pkc)

fig, axes = plt.subplots(4, 2, figsize=(12, 16), sharex=True, sharey='row')
axes[0, 0].set_title("SUCCESSFUL", fontsize=14, fontweight='bold', color='cornflowerblue')
axes[0, 1].set_title("FAILED", fontsize=14, fontweight='bold', color='indianred')

print(f"Processing Onset-Aligned Plot (Min Reaches = {MIN_REACHES})...")

for cat_idx, target_cat in enumerate([1, 2, 3, 4]):
    for col_idx, target_succ in enumerate([1, 0]):
        
        ax = axes[cat_idx, col_idx]
        color = 'cornflowerblue' if target_succ == 1 else 'indianred'
        
        # --- A. COUNT REACHES ---
        row0 = df_pkc.iloc[0]
        r_succ = row0['idx_reach_succ']
        r_cat = row0['idx_reach_cat']
        
        reach_count = 0
        if r_succ is not None:
            for s, c in zip(r_succ, r_cat):
                s_val = s[0] if isinstance(s, list) else s
                c_val = c[0] if isinstance(c, list) else c
                if s_val == target_succ and c_val == target_cat:
                    reach_count += 1
        
        # --- B. COLLECT AND ALIGN DATA ---
        all_aligned_spikes = []
        
        if reach_count >= MIN_REACHES:
            for _, row in df_pkc.iterrows():
                # Need Trajectories (reach_px) for alignment
                r_spikes = row['st_reach']
                r_px = row['reach_px']
                n_succ = row['idx_reach_succ']
                n_cat = row['idx_reach_cat']
                
                if r_spikes is None or r_px is None: continue
                
                for spikes, px, s, c in zip(r_spikes, r_px, n_succ, n_cat):
                    s_val = s[0] if isinstance(s, list) else s
                    c_val = c[0] if isinstance(c, list) else c
                    
                    if s_val == target_succ and c_val == target_cat:
                        if spikes is not None and px is not None:
                            
                            # 1. Detect Onset
                            traj = np.array([px[0], px[1], px[2]]).T
                            onset_frame = get_arbitrated_onset(traj, FPS, global_pca)
                            onset_sec = onset_frame / FPS
                            
                            # 2. Shift Spikes
                            flat_reach = np.array(spikes).ravel().tolist()
                            for t in flat_reach:
                                aligned_t = t - onset_sec # Shift: 0 is now onset
                                if WINDOW_SEC[0] <= aligned_t <= WINDOW_SEC[1]:
                                    all_aligned_spikes.append(aligned_t)

            # Plot Histogram
            if len(all_aligned_spikes) > 0:
                data_to_plot = np.array(all_aligned_spikes)
                total_seconds = reach_count * num_neurons * bin_width
                weights = np.full(len(data_to_plot), 1 / total_seconds)
                
                ax.hist(data_to_plot, bins=bins, weights=weights, 
                        color=color, alpha=0.7, edgecolor='white', linewidth=0.5)
                
                ax.text(0.95, 0.95, f"N={reach_count}", 
                        transform=ax.transAxes, ha='right', va='top', fontsize=10, 
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9))
            else:
                ax.text(0.5, 0.5, f"N={reach_count}\n(Silent)", 
                        transform=ax.transAxes, ha='center', va='center', color='gray')

        else:
            ax.text(0.5, 0.5, f"Insufficient Data\n(N={reach_count})", 
                    transform=ax.transAxes, ha='center', va='center', 
                    color='gray', fontstyle='italic', fontsize=12)

        # --- C. FORMATTING ---
        ax.set_xlim(WINDOW_SEC)
        ax.set_ylim([0, 6])
        ax.axvline(0, color='red', linestyle='--', alpha=0.5, label='Onset') # Red line for Onset
        ax.grid(True, alpha=0.2)
        
        if col_idx == 0:
            ax.set_ylabel(f"{cat_labels[target_cat]}\nFiring Rate (Hz)", fontsize=11, fontweight='bold')
        if cat_idx == 3:
            ax.set_xlabel("Time relative to Onset (s)")

plt.suptitle(f"Population CS Activity - Onset Aligned (N={num_neurons} cells)", fontsize=16, y=0.99)
plt.tight_layout()
plt.show()

# Optional: Save
# plt.savefig(save_path + "\\population_reach_activity_onset.png", dpi=300)
# %%


#  ██████╗███████╗    ██╗███╗   ██╗    ██████╗ ███████╗ █████╗  ██████╗██╗  ██╗    
# ██╔════╝██╔════╝    ██║████╗  ██║    ██╔══██╗██╔════╝██╔══██╗██╔════╝██║  ██║    
# ██║     ███████╗    ██║██╔██╗ ██║    ██████╔╝█████╗  ███████║██║     ███████║    
# ██║     ╚════██║    ██║██║╚██╗██║    ██╔══██╗██╔══╝  ██╔══██║██║     ██╔══██║    
# ╚██████╗███████║    ██║██║ ╚████║    ██║  ██║███████╗██║  ██║╚██████╗██║  ██║    
#  ╚═════╝╚══════╝    ╚═╝╚═╝  ╚═══╝    ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝    
                                                                                 
# --- CONFIGURATION ---
FPS = 400
CENTER_FRAME = 150
# Create time axis for trajectory (matches your 201 frames logic)
# [-0.375s to +0.125s approximately]
traj_time_axis = (np.arange(201) - CENTER_FRAME) / FPS 

# Selection Criteria
TARGET_CAT = 1
TARGET_SUCC = 0

def analyze_single_reach(df, specific_idx=None):
    """
    Analyzes one specific reach event across the entire population.
    If specific_idx is None, picks a random one matching criteria.
    """
    
    # 1. FIND A VALID REACH INDEX
    # We use the first neuron's metadata as the "master" list
    row0 = df.iloc[0]
    r_succ = row0['idx_reach_succ']
    r_cat = row0['idx_reach_cat']
    r_px = row0['reach_px'] # Trajectory data
    
    valid_indices = []
    
    # Identify all valid reaches
    for i, (s, c) in enumerate(zip(r_succ, r_cat)):
        s_val = s[0] if isinstance(s, list) else s
        c_val = c[0] if isinstance(c, list) else c
        
        if s_val == TARGET_SUCC and c_val == TARGET_CAT:
            valid_indices.append(i)
            
    if not valid_indices:
        print("No matches found!")
        return

    # Randomly select if not specified
    if specific_idx is None:
        target_idx = random.choice(valid_indices)
    else:
        target_idx = specific_idx
        
    print(f"Analyzing Reach Index: {target_idx}")
    
    # 2. EXTRACT TRAJECTORY (Kinematics)
    # px is usually [x, y, z] arrays
    raw_traj = r_px[target_idx] 
    # Let's take the X-component (index 0) as primary, or compute speed
    traj_x = np.array(raw_traj[0]).flatten()
    
    # 3. COLLECT SPIKES (Population Activity)
    active_neurons = 0
    total_spikes = 0
    spike_times_flat = []
    neural_raster = [] # List of (neuron_idx, spike_time)
    
    print(f"Scanning {len(df)} neurons...")
    
    for neuron_idx, row in df.iterrows():
        # Get spikes for THIS specific reach
        # Note: st_reach is a list of lists. We access by index.
        spikes = row['st_reach'][target_idx]
        
        if spikes is not None and len(spikes) > 0:
            # Flatten & Align
            s_arr = np.array(spikes).ravel()
            
            # Align to center (0.375s)
            aligned_s = s_arr - (150/400)
            
            # Filter to relevant window only (matching trajectory)
            valid_s = aligned_s[(aligned_s >= traj_time_axis[0]) & (aligned_s <= traj_time_axis[-1])]
            
            if len(valid_s) > 0:
                active_neurons += 1
                total_spikes += len(valid_s)
                spike_times_flat.extend(valid_s)
                
                # For raster plot
                for t in valid_s:
                    neural_raster.append((neuron_idx, t))

    # --- REPORTING ---
    participation_rate = (active_neurons / len(df)) * 100
    print("-" * 30)
    print(f"REACH DATA SUMMARY (ID: {target_idx})")
    print("-" * 30)
    print(f"Total Spikes:       {total_spikes}")
    print(f"Active Neurons:     {active_neurons} / {len(df)} ({participation_rate:.1f}%)")
    print(f"Avg Spikes/Active:  {total_spikes/active_neurons:.2f}" if active_neurons > 0 else "N/A")
    
    # --- PLOTTING ---
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True, 
                             gridspec_kw={'height_ratios': [2, 3, 1], 'hspace': 0.1})
    
    # Plot 1: Kinematics (Trajectory)
    ax_kin = axes[0]
    ax_kin.plot(traj_time_axis, traj_x, color='black', lw=2, label='X-Position')
    ax_kin.set_ylabel("Position (px)")
    ax_kin.set_title(f"Reach #{target_idx} | Cat {TARGET_CAT} | Succ {TARGET_SUCC}")
    ax_kin.axvline(0, color='red', linestyle='--', alpha=0.5, label='Center')
    ax_kin.legend(loc='upper right')
    ax_kin.grid(True, alpha=0.3)
    
    # Plot 2: Raster Plot (Who fired and when?)
    ax_rast = axes[1]
    if neural_raster:
        y_pts, x_pts = zip(*neural_raster)
        # Using scatter for raster
        ax_rast.scatter(x_pts, y_pts, s=15, color='black', alpha=0.7, marker='|')
    
    ax_rast.set_ylabel("Neuron ID")
    ax_rast.set_ylim(-1, len(df))
    ax_rast.grid(True, alpha=0.1)
    
    # Add text box with stats
    stats_text = (f"Total Spikes: {total_spikes}\n"
                  f"Active Cells: {active_neurons}")
    ax_rast.text(0.02, 0.95, stats_text, transform=ax_rast.transAxes, 
                 va='top', bbox=dict(facecolor='white', alpha=0.8))

    # Plot 3: Population Histogram (Summed activity)
    ax_hist = axes[2]
    if spike_times_flat:
        bins = np.arange(traj_time_axis[0], traj_time_axis[-1], 0.01) # 10ms bins
        ax_hist.hist(spike_times_flat, bins=bins, color='teal', alpha=0.7)
    
    ax_hist.set_ylabel("Spike Count")
    ax_hist.set_xlabel("Time relative to Reach Center (s)")
    ax_hist.axvline(0, color='red', linestyle='--')
    ax_hist.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

# --- RUN IT ---
# Will pick a random reach each time you run this cell
analyze_single_reach(df_pkc) 
# %%

# ██████╗  █████╗ ████████╗ ██████╗██╗  ██╗    
# ██╔══██╗██╔══██╗╚══██╔══╝██╔════╝██║  ██║    
# ██████╔╝███████║   ██║   ██║     ███████║    
# ██╔══██╗██╔══██║   ██║   ██║     ██╔══██║    
# ██████╔╝██║  ██║   ██║   ╚██████╗██║  ██║    
# ╚═════╝ ╚═╝  ╚═╝   ╚═╝    ╚═════╝╚═╝  ╚═╝    
                                             

def analyze_reaches_batched_pca_with_marks(df, target_cat=1, target_succ=1, batch_size=8):
    
    print("Step 1: collecting all trajectories for Global PCA...")
    row0 = df.iloc[0]
    r_succ = row0['idx_reach_succ']
    r_cat = row0['idx_reach_cat']
    r_px = row0['reach_px']
    
    valid_indices = []
    all_3d_points = []
    
    # --- 1. COLLECT DATA for PCA ---
    for i, (s, c) in enumerate(zip(r_succ, r_cat)):
        s_val = s[0] if isinstance(s, list) else s
        c_val = c[0] if isinstance(c, list) else c
        
        if s_val == target_succ and c_val == target_cat:
            if i < len(r_px) and r_px[i] is not None:
                # px is [x_array, y_array, z_array]
                try:
                    # Stack to (N_samples, 3)
                    # transpose because we want (Time, Dimensions)
                    traj_3d = np.array([r_px[i][0], r_px[i][1], r_px[i][2]]).T 
                    
                    # Ensure it's valid numeric data
                    if traj_3d.shape[1] == 3 and traj_3d.shape[0] > 10:
                        valid_indices.append(i)
                        all_3d_points.append(traj_3d)
                except Exception:
                    continue

    if not all_3d_points:
        print("No valid 3D data found.")
        return

    # --- 2. FIT GLOBAL PCA ---
    # Stack all timepoints from all reaches to find the "Main Axis of Movement"
    # shape: (Total_Frames_All_Reaches, 3)
    grand_matrix = np.vstack(all_3d_points)
    
    pca = PCA(n_components=1)
    pca.fit(grand_matrix)
    print(f"PCA fit complete. Explained Variance: {pca.explained_variance_ratio_[0]:.2%}")
    
    # --- 3. BATCH PLOTTING ---
    total_reaches = len(valid_indices)
    FPS = 400
    CENTER_FRAME = 150
    t_axis = (np.arange(201) - CENTER_FRAME) / FPS 

    print(f"Plotting {total_reaches} reaches in batches of {batch_size}...")

    for start_i in range(0, total_reaches, batch_size):
        end_i = min(start_i + batch_size, total_reaches)
        batch_indices = valid_indices[start_i:end_i]
        
        fig, axes = plt.subplots(3, batch_size, figsize=(24, 10), sharex=True, 
                                 gridspec_kw={'height_ratios': [1.5, 3, 1], 'hspace': 0.1, 'wspace': 0.3})
        
        if batch_size == 1: axes = axes.reshape(3, 1)

        for col_idx, reach_idx in enumerate(batch_indices):
            # Select column axes
            if batch_size > 1:
                ax_traj, ax_rast, ax_hist = axes[0, col_idx], axes[1, col_idx], axes[2, col_idx]
            else:
                ax_traj, ax_rast, ax_hist = axes[0], axes[1], axes[2]

            # --- A. PROJECT TRAJECTORY ---
            # Get 3D data
            raw_px = row0['reach_px'][reach_idx]
            traj_3d = np.array([raw_px[0], raw_px[1], raw_px[2]]).T
            
            # Project onto PC1
            pc1_proj = pca.transform(traj_3d).flatten()
            
            # Use data length to slice time axis
            path_len = len(pc1_proj)
            if path_len >= 201:
                y_plot = pc1_proj[:201]
                x_plot = t_axis
            else:
                y_plot = pc1_proj
                x_plot = t_axis[:path_len]

            # Plot Trajectory Line
            ax_traj.plot(x_plot, y_plot, color='#8c564b', lw=2, alpha=0.8, label='PC1 Proj') 
            
            # --- B. COLLECT SPIKES ---
            all_spikes_in_window = []
            raster_x, raster_y = [], []
            
            for neuron_id, row in df.iterrows():
                spikes = row['st_reach'][reach_idx]
                if spikes is not None:
                    # Flatten standard list/array
                    s_arr = np.array(spikes).ravel() 
                    # Align to center (0s)
                    aligned_s = s_arr - (CENTER_FRAME/FPS)
                    
                    # Filter for plot window
                    valid = aligned_s[(aligned_s >= t_axis[0]) & (aligned_s <= t_axis[-1])]
                    
                    if len(valid) > 0:
                        all_spikes_in_window.extend(valid)
                        raster_x.extend(valid)
                        raster_y.extend([neuron_id] * len(valid))

            # --- C. PLOT SPIKES ON TRAJECTORY ---
            # For every spike, find its Y-position on the PC1 curve
            if len(all_spikes_in_window) > 0:
                # Interpolate PC1 value at spike times
                # np.interp requires sorted x-values (t_axis is sorted)
                y_spikes = np.interp(all_spikes_in_window, x_plot, y_plot)
                
                # Plot red dots ON the curve
                ax_traj.scatter(all_spikes_in_window, y_spikes, 
                                c='red', s=10, alpha=0.6, zorder=5, label='CS Event')

            ax_traj.axvline(0, color='gray', linestyle='--', alpha=0.3)
            ax_traj.set_title(f"Idx {reach_idx}", fontsize=10)
            ax_traj.grid(True, alpha=0.2)
            
            # --- D. STANDARD PLOTS ---
            # Raster
            if raster_x:
                ax_rast.scatter(raster_x, raster_y, s=3, c='black', alpha=0.6, marker='|')
            ax_rast.set_ylim(-1, len(df))
            ax_rast.axvline(0, color='gray', linestyle='--', alpha=0.3)
            ax_rast.grid(True, alpha=0.1)
            
            # Histogram
            if all_spikes_in_window:
                ax_hist.hist(all_spikes_in_window, bins=np.linspace(t_axis[0], t_axis[-1], 20), 
                             color='teal', alpha=0.7)
            ax_hist.axvline(0, color='gray', linestyle='--', alpha=0.3)
            ax_hist.grid(True, alpha=0.2)
            ax_hist.set_xlabel("Time (s)")

            # Labels (First col only)
            if col_idx == 0:
                ax_traj.set_ylabel("PC1 (a.u.)")
                ax_rast.set_ylabel("Neuron ID")
                ax_hist.set_ylabel("Spike Count")
            else:
                ax_traj.set_yticks([])
                ax_rast.set_yticks([])
                ax_hist.set_yticks([])
        
        # Hide empty columns
        for j in range(len(batch_indices), batch_size):
            if batch_size > 1:
                axes[0, j].axis('off')
                axes[1, j].axis('off')
                axes[2, j].axis('off')

        plt.suptitle(f"Global PC1 & CS Events (Batch {start_i}-{end_i})", y=0.98)
        #plt.tight_layout()
        plt.subplots_adjust(top=0.92)
        plt.show()

# --- EXECUTE ---
analyze_reaches_batched_pca_with_marks(df_pkc, target_cat=1, target_succ=1, batch_size=8)
# %%

# ████████╗███████╗██████╗  ██████╗██╗██╗     ███████╗    
# ╚══██╔══╝██╔════╝██╔══██╗██╔════╝██║██║     ██╔════╝    
#    ██║   █████╗  ██████╔╝██║     ██║██║     █████╗      
#    ██║   ██╔══╝  ██╔══██╗██║     ██║██║     ██╔══╝      
#    ██║   ███████╗██║  ██║╚██████╗██║███████╗███████╗    
#    ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝╚═╝╚══════╝╚══════╝    
                                                        



def analyze_reach_types_per_cell(df, target_cat=1, target_succ=1, neurons_per_page=15):
    
    print("Step 1: Fitting Global PCA...")
    valid_reaches_for_pca = []
    row0 = df.iloc[0]
    r_px = row0['reach_px']
    for px in r_px:
        if px is not None:
            try:
                traj = np.array([px[0], px[1], px[2]]).T
                if traj.shape[0] > 100: valid_reaches_for_pca.append(traj)
            except: continue
    
    if not valid_reaches_for_pca: return
    grand_matrix = np.vstack(valid_reaches_for_pca)
    global_pca = PCA(n_components=1)
    global_pca.fit(grand_matrix)

    print("Step 2: Defining Terciles...")
    reach_meta = {} 
    r_succ = row0['idx_reach_succ']
    r_cat = row0['idx_reach_cat']
    valid_reach_indices = []

    for i, (s, c) in enumerate(zip(r_succ, r_cat)):
        s_val = s[0] if isinstance(s, list) else s
        c_val = c[0] if isinstance(c, list) else c
        
        if s_val == target_succ and c_val == target_cat:
            if i < len(r_px) and r_px[i] is not None:
                traj_3d = np.array([r_px[i][0], r_px[i][1], r_px[i][2]]).T
                onset = get_arbitrated_onset(traj_3d, FPS, global_pca)
                offset = 150 # Center frame
                
                if onset < offset:
                    duration = offset - onset
                    t_len = duration / 3.0
                    reach_meta[i] = (onset, onset+t_len, onset+(2*t_len), offset)
                    valid_reach_indices.append(i)

    print("Step 3: Classifying Reaches...")
    cell_results = []
    
    # Categories of Reach Activity
    LABELS = ['Silent', 'T1', 'T2', 'T3', 'T1+T2', 'T1+T3', 'T2+T3', 'All']
    # Color map for the bars
    COLORS = {
        'Silent': '#f0f0f0', # Light Grey
        'T1': '#ff9999',     # Red
        'T2': '#66b3ff',     # Blue
        'T3': '#99ff99',     # Green
        'T1+T2': '#ba8cd7',  # Purple
        'T1+T3': '#ffcc00',  # Yellow
        'T2+T3': '#40e0d0',  # Turquoise
        'All': '#696969'     # Dark Grey
    }

    for neuron_idx, row in df.iterrows():
        spikes_list = row['st_reach']
        
        # Initialize counters for this neuron
        type_counts = {lbl: 0 for lbl in LABELS}
        
        if spikes_list is not None:
            for r_idx in valid_reach_indices:
                b1, b2, b3, b4 = reach_meta[r_idx]
                spikes = spikes_list[r_idx]
                
                if spikes is None or len(spikes) == 0:
                    type_counts['Silent'] += 1
                    continue
                
                sf = np.array(spikes).ravel() * FPS 
                
                # Check presence (>0 spikes) in each bin
                has_t1 = np.any((sf >= b1) & (sf < b2))
                has_t2 = np.any((sf >= b2) & (sf < b3))
                has_t3 = np.any((sf >= b3) & (sf <= b4))
                
                # Determine Label
                if not (has_t1 or has_t2 or has_t3):
                    tag = 'Silent' # Spikes existed but outside the Onset-Offset window
                elif has_t1 and has_t2 and has_t3:
                    tag = 'All'
                elif has_t1 and has_t2: tag = 'T1+T2'
                elif has_t1 and has_t3: tag = 'T1+T3'
                elif has_t2 and has_t3: tag = 'T2+T3'
                elif has_t1: tag = 'T1'
                elif has_t2: tag = 'T2'
                elif has_t3: tag = 'T3'
                else: tag = 'Silent'
                
                type_counts[tag] += 1
        
        # Convert to list for plotting (maintain order)
        counts_ordered = [type_counts[l] for l in LABELS]
        total_reaches = sum(counts_ordered)
        
        # Calculate "% Active" (Non-Silent)
        active_count = total_reaches - type_counts['Silent']
        active_pct = (active_count / total_reaches * 100) if total_reaches > 0 else 0
        
        cell_results.append({
            'id': neuron_idx,
            'counts': counts_ordered,
            'total': total_reaches,
            'active_pct': active_pct
        })

    # 4. PLOTTING
    total_cells = len(cell_results)
    
    for start_i in range(0, total_cells, neurons_per_page):
        end_i = min(start_i + neurons_per_page, total_cells)
        batch = cell_results[start_i:end_i]
        
        rows, cols = 5, 4
        fig, axes = plt.subplots(rows, cols, figsize=(14, 18))
        axes = axes.flatten()
        
        for i, cell in enumerate(batch):
            ax = axes[i]
            
            # Remove 'Silent' from plot to focus on activity types?
            # Or keep it to show participation rate? Let's keep it but maybe separate.
            # User likely wants to know logic of ACTIVITY.
            # Let's plot only the Active Types (indices 1 to 7)
            
            x_labels = LABELS[1:] # Skip Silent
            y_values = cell['counts'][1:] 
            bar_colors = [COLORS[l] for l in x_labels]
            
            bars = ax.bar(x_labels, y_values, color=bar_colors, edgecolor='black', alpha=0.7)
            
            # Annotate
            for bar, val in zip(bars, y_values):
                if val > 0:
                    ax.text(bar.get_x() + bar.get_width()/2, val, str(val), 
                            ha='center', va='bottom', fontsize=9)
            
            # Info Title
            silent_count = cell['counts'][0]
            ax.set_title(f"Cell {cell['id']} (Active: {int(cell['active_pct'])}% | Silent: {silent_count})", 
                         fontsize=10, fontweight='bold')
            
            ax.tick_params(axis='x', rotation=45, labelsize=8)
            ax.set_ylim(bottom=0)
            if max(y_values) > 0: ax.set_ylim(0, max(y_values)*1.2)
            
        for j in range(len(batch), len(axes)): axes[j].axis('off')
            
        plt.suptitle(f"Reach Type Counts (Batch {start_i}-{end_i-1})", fontsize=16, y=0.99)
        plt.tight_layout()
        plt.show()

# --- EXECUTE ---
analyze_reach_types_per_cell(df_pkc, target_cat=1, target_succ=0, neurons_per_page=20)
# %%

# ██╗    ██╗     ██████╗ ██████╗     ██╗    ██╗ ██████╗      ██████╗███████╗
# ██║    ██║    ██╔═══██╗██╔══██╗    ██║    ██║██╔═══██╗    ██╔════╝██╔════╝
# ██║ █╗ ██║    ██║   ██║██████╔╝    ██║ █╗ ██║██║   ██║    ██║     ███████╗
# ██║███╗██║    ██║   ██║██╔══██╗    ██║███╗██║██║   ██║    ██║     ╚════██║
# ╚███╔███╔╝    ╚██████╔╝██║  ██║    ╚███╔███╔╝╚██████╔╝    ╚██████╗███████║
#  ╚══╝╚══╝      ╚═════╝ ╚═╝  ╚═╝     ╚══╝╚══╝  ╚═════╝      ╚═════╝╚══════╝
                                                                          

def compute_cs_probability_map_combined(df, target_cat=1, target_succ=1, bin_size_ms=50):
    
    print("Step 1: Fitting Global PCA for Onset Detection...")
    valid_reaches_for_pca = []
    row0 = df.iloc[0]
    r_px = row0['reach_px']
    for px in r_px:
        if px is not None:
            try:
                traj = np.array([px[0], px[1], px[2]]).T
                if traj.shape[0] > 100: valid_reaches_for_pca.append(traj)
            except: continue
                
    if not valid_reaches_for_pca: return
    grand_matrix = np.vstack(valid_reaches_for_pca)
    global_pca = PCA(n_components=1)
    global_pca.fit(grand_matrix)
    
    # --- BINNING SETUP ---
    FPS = 400
    bin_sec = bin_size_ms / 1000.0
    t_start, t_end = -2, 2 
    bins = np.arange(t_start, t_end + bin_sec, bin_sec)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    n_bins = len(bins) - 1

    print(f"Step 2: Processing Reaches (Bin Size: {bin_size_ms}ms)...")
    
    reach_onsets = {} 
    r_succ = row0['idx_reach_succ']
    r_cat = row0['idx_reach_cat']
    valid_reach_indices = []

    for i, (s, c) in enumerate(zip(r_succ, r_cat)):
        s_val = s[0] if isinstance(s, list) else s
        c_val = c[0] if isinstance(c, list) else c
        
        if s_val == target_succ and c_val == target_cat:
            if i < len(r_px) and r_px[i] is not None:
                traj_3d = np.array([r_px[i][0], r_px[i][1], r_px[i][2]]).T
                onset = get_arbitrated_onset(traj_3d, FPS, global_pca)
                onset_sec = onset / FPS
                reach_onsets[i] = onset_sec
                valid_reach_indices.append(i)

    print(f"Valid Reaches: {len(valid_reach_indices)}")

    # Matrix: (Neurons x Bins)
    pop_prob_matrix = np.zeros((len(df), n_bins))
    
    for n_idx, (_, row) in enumerate(df.iterrows()):
        spikes_list = row['st_reach']
        if spikes_list is None: continue
        
        neuron_trial_matrix = np.zeros((len(valid_reach_indices), n_bins))
        
        for k, r_idx in enumerate(valid_reach_indices):
            onset_t = reach_onsets[r_idx]
            spikes = spikes_list[r_idx]
            
            if spikes is not None:
                s_arr = np.array(spikes).ravel()
                s_rel = s_arr - onset_t
                
                # Digitize is the fastest way to bin
                bin_indices = np.digitize(s_rel, bins) - 1
                valid_b = bin_indices[(bin_indices >= 0) & (bin_indices < n_bins)]
                
                # Binary flag for probability calculation
                neuron_trial_matrix[k, valid_b] = 1
        
        # Mean across trials = Probability
        if len(valid_reach_indices) > 0:
            pop_prob_matrix[n_idx, :] = np.mean(neuron_trial_matrix, axis=0)

    # --- 3. SORTING ---
    peak_times = np.argmax(pop_prob_matrix, axis=1)
    sorted_indices = np.argsort(peak_times)
    sorted_matrix = pop_prob_matrix[sorted_indices, :]

    # --- 4. PLOTTING COMBINED ---
    fig, (ax_heat, ax_peth) = plt.subplots(2, 1, figsize=(12, 12), 
                                           gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.1},
                                           sharex=True)
    
    # --- A. HEATMAP ---
    cmap = LinearSegmentedColormap.from_list("custom", ["white", "cornflowerblue", "firebrick"])
    
    # We use imshow for better control, but aspect='auto' is key
    im = ax_heat.imshow(sorted_matrix, aspect='auto', cmap=cmap, 
                        interpolation='nearest', vmin=0, vmax=0.4,
                        extent=[t_start, t_end, len(df), 0]) # Y-axis inverted (0 at top)
    
    # Add Colorbar nicely
    cbar = plt.colorbar(im, ax=ax_heat, pad=0.02)
    cbar.set_label(f'P(CS) per {bin_size_ms}ms', rotation=270, labelpad=15)
    
    ax_heat.set_ylabel("Neurons (Sorted by Peak Time)")
    ax_heat.set_title(f"CS Probability Map (N={len(df)} cells)")
    ax_heat.axvline(0, color='black', linestyle='--')

    # --- B. PETH (Population Mean) ---
    mean_prob = np.mean(pop_prob_matrix, axis=0)
    sem_prob = np.std(pop_prob_matrix, axis=0) / np.sqrt(len(df))
    
    # Plot
    ax_peth.plot(bin_centers, mean_prob, color='black', lw=2, label='Population Mean')
    ax_peth.fill_between(bin_centers, mean_prob - sem_prob, mean_prob + sem_prob, 
                         color='black', alpha=0.2)
    
    ax_peth.set_ylabel("Avg Probability")
    ax_peth.set_xlabel("Time from Movement Onset (s)")
    ax_peth.axvline(0, color='red', linestyle='--')
    ax_peth.grid(True, alpha=0.3)
    ax_peth.set_ylim(bottom=0)
    
    # Peak Annotation
    peak_idx = np.argmax(mean_prob)
    ax_peth.plot(bin_centers[peak_idx], mean_prob[peak_idx], 'ro')
    ax_peth.text(bin_centers[peak_idx], mean_prob[peak_idx] + 0.005, 
                 f"Peak: {mean_prob[peak_idx]:.1%}", ha='center', va='bottom', fontsize=9)
    
    plt.suptitle(f"Population CS Dynamics (Start-Aligned)", y=0.95, fontsize=16)
    plt.show()
    
    return pop_prob_matrix

# --- EXECUTE ---
probs = compute_cs_probability_map_combined(df_pkc, bin_size_ms=50)
# %%

# ████████╗██████╗  █████╗      ██╗
# ╚══██╔══╝██╔══██╗██╔══██╗     ██║
#    ██║   ██████╔╝███████║     ██║
#    ██║   ██╔══██╗██╔══██║██   ██║
#    ██║   ██║  ██║██║  ██║╚█████╔╝
#    ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚════╝ 
                                  
def analyze_trajectory_by_cs_presence(df, target_cat=1, target_succ=1, window_sec=[-0.3, 0.3]):
    """
    Plots the average PC1 trajectory for reaches WITH vs WITHOUT a Complex Spike.
    Two modes:
    1. Population: Any CS from ANY neuron in the window.
    2. Per Neuron: Specific CS from specific neuron.
    """
    
    print("Step 1: Preparing Trajectories & Onsets...")
    
    # --- 1. SHARED SETUP: TRAJECTORIES & ONSETS ---
    row0 = df.iloc[0]
    r_px = row0['reach_px']
    r_succ = row0['idx_reach_succ']
    r_cat = row0['idx_reach_cat']

    # Helper: Fit PCA if needed (or reuse global logic)
    valid_data_for_pca = []
    for px in r_px:
        if px is not None:
            try:
                t = np.array([px[0], px[1], px[2]]).T
                if len(t) > 50: valid_data_for_pca.append(t)
            except: continue
    
    if not valid_data_for_pca: return
    
    # Fit PCA
    grand_matrix = np.vstack(valid_data_for_pca)
    pca = PCA(n_components=1)
    pca.fit(grand_matrix)
    
    # Collect Processed Reaches
    # Structure: {reach_idx: {'pc1': array, 'time': array_relative_to_onset}}
    processed_reaches = {}
    
    # Common Time Axis for Interpolation (e.g. -0.5s to +0.5s)
    common_time = np.linspace(-0.5, 0.5, 400) # 400 points
    
    valid_indices = []
    
    for i, (s, c) in enumerate(zip(r_succ, r_cat)):
        s_val = s[0] if isinstance(s, list) else s
        c_val = c[0] if isinstance(c, list) else c
        
        if s_val == target_succ and c_val == target_cat:
            if i < len(r_px) and r_px[i] is not None:
                # 3D Trajectory
                traj_3d = np.array([r_px[i][0], r_px[i][1], r_px[i][2]]).T
                
                # Computed Onset
                onset_frame = get_arbitrated_onset(traj_3d, FPS, pca)
                onset_sec = onset_frame / FPS
                
                # PC1 Projection
                pc1 = pca.transform(traj_3d).flatten()
                
                # Create Time Axis relative to Onset
                # Frame 0 is time 0 in original recording? typically yes.
                # So t = (frame - onset) / FPS
                t_frames = np.arange(len(pc1))
                t_rel = (t_frames - onset_frame) / FPS
                
                # Interpolate to Common Time Axis (Data Standardization)
                # We align purely by time.
                pc1_interp = np.interp(common_time, t_rel, pc1, left=np.nan, right=np.nan)
                
                processed_reaches[i] = {
                    'pc1_interp': pc1_interp,
                    'onset_sec': onset_sec,
                    'orig_pc1': pc1,
                    'orig_t': t_rel
                }
                valid_indices.append(i)

    print(f"Processed {len(valid_indices)} reaches.")

    # --- HELPER PLOT FUNCTION ---
    def plot_comparison(ax, group_with, group_without, title_str=""):
        # Plot Individual Traces (Faint)
        # Without CS (Blue)
        for r_idx in group_without:
            trace = processed_reaches[r_idx]['pc1_interp']
            ax.plot(common_time, trace, color='tab:blue', alpha=0.2, lw=1)
            
        # With CS (Red)
        for r_idx in group_with:
            trace = processed_reaches[r_idx]['pc1_interp']
            ax.plot(common_time, trace, color='tab:red', alpha=0.2, lw=1)
            
        # Calculate Mean & Std (ignoring NaNs from interpolation edges)
        if len(group_without) > 0:
            stack_no = np.vstack([processed_reaches[i]['pc1_interp'] for i in group_without])
            mean_no = np.nanmean(stack_no, axis=0)
            std_no = np.nanstd(stack_no, axis=0)
            ax.plot(common_time, mean_no, color='tab:blue', lw=3, label=f'No CS (N={len(group_without)})')
            ax.fill_between(common_time, mean_no-std_no, mean_no+std_no, color='tab:blue', alpha=0.01)
            
        if len(group_with) > 0:
            stack_yes = np.vstack([processed_reaches[i]['pc1_interp'] for i in group_with])
            mean_yes = np.nanmean(stack_yes, axis=0)
            std_yes = np.nanstd(stack_yes, axis=0)
            ax.plot(common_time, mean_yes, color='tab:red', lw=3, label=f'With CS (N={len(group_with)})')
            ax.fill_between(common_time, mean_yes-std_yes, mean_yes+std_yes, color='tab:red', alpha=0.01)
            
        ax.axvline(0, color='black', linestyle='--', label='Onset')
        ax.set_title(title_str)
        ax.set_xlabel('Time from Onset (s)')
        ax.set_ylabel('PC1 Projection (a.u.)')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # --- 2. POPULATION LEVEL ANALYSIS ---
    # Condition: ANY CS in the window from ANY neuron
    print("Running Population Analysis...")
    pop_with_cs = set()
    pop_without_cs = set()
    
    for r_idx in valid_indices:
        onset_sec = processed_reaches[r_idx]['onset_sec']
        has_cs = False
        
        # Check ALL neurons
        for _, row in df.iterrows():
            spikes = row['st_reach'][r_idx]
            if spikes is not None:
                s_arr = np.array(spikes).ravel()
                # Relative check: Spike Time - Onset Time
                # Window is [-0.3, 0.3] relative to onset
                # So raw spike must be in [Onset-0.3, Onset+0.3]
                if np.any((s_arr >= onset_sec + window_sec[0]) & (s_arr <= onset_sec + window_sec[1])):
                    has_cs = True
                    break # Found one, that's enough
        
        if has_cs:
            pop_with_cs.add(r_idx)
        else:
            pop_without_cs.add(r_idx)
            
    # Plot Population
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_comparison(ax, pop_with_cs, pop_without_cs, 
                    title_str=f"Population: Reaches with ANY CS vs NONE ({window_sec}s window)")
    plt.show()


    # --- 3. SINGLE NEURON ANALYSIS (Example: Top 3 Most Active) ---
    print("Running Single Neuron Analysis (Top 3)...")
    
    # Compute total spikes to find most active cells
    spike_counts = []
    for idx, row in df.iterrows():
        total = 0
        if row['st_reach'] is not None:
            for s in row['st_reach']: 
                if s is not None: total += len(s)
        spike_counts.append((idx, total))
        
    top_neurons = sorted(spike_counts, key=lambda x: x[1], reverse=True)[:3]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for ax, (n_idx, total) in zip(axes, top_neurons):
        group_yes = set()
        group_no = set()
        
        row = df.iloc[n_idx]
        spikes_list = row['st_reach']
        
        for r_idx in valid_indices:
            onset_sec = processed_reaches[r_idx]['onset_sec']
            has_cs = False
            
            if spikes_list is not None:
                spikes = spikes_list[r_idx]
                if spikes is not None:
                    s_arr = np.array(spikes).ravel()
                    if np.any((s_arr >= onset_sec + window_sec[0]) & (s_arr <= onset_sec + window_sec[1])):
                        has_cs = True
            
            if has_cs: group_yes.add(r_idx)
            else: group_no.add(r_idx)
            
        plot_comparison(ax, group_yes, group_no, title_str=f"Neuron {n_idx}")
        
    plt.suptitle(f"Single Neuron Effect on Trajectory ({window_sec}s window)", y=0.98)
    plt.tight_layout()
    plt.show()

# --- EXECUTE ---
analyze_trajectory_by_cs_presence(df_pkc, target_cat=1, target_succ=1, window_sec=[-0.3, 0.3])
# %%
