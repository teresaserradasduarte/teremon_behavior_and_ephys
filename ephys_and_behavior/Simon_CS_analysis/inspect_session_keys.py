import h5py
import numpy as np

file_path = r"C:\Users\simon\Desktop\Data\ephys_and_behavior"
save_path = r"C:\Users\simon\Desktop\Data\ephys_and_behavior\1_CoteDor\R5"
mat_file = f'{save_path}\\session_reaching_data_paw.mat'

def recursive_resolve(f, item):
    if isinstance(item, h5py.Group):
        return {k: recursive_resolve(f, item[k]) for k in item.keys()}
    if isinstance(item, h5py.Dataset):
        return "Dataset" # just peek at keys for now
    return item

with h5py.File(mat_file, 'r') as f:
    if 'session' in f:
        print("Keys in 'session':")
        session_group = f['session']
        print(list(session_group.keys()))
    else:
        print("'session' key not found in MAT file.")
