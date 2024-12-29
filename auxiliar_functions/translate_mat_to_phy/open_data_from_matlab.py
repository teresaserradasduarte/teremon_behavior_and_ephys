#%% Import from matlab
import pathlib
import scipy.io
import pandas as pd

#%% behavior
data_path = pathlib.WindowsPath(r"E:\26082023_Lindt_StrCer_S4_g0\26082023_Lindt_StrCer_S4_g0_imec0\catGT\kilosort4\simplified_data")
behav_data = scipy.io.loadmat(data_path / 'behav_info.mat')

# Extract the struct
behav_struct = behav_data['behav_info']

# behaviorally relevant window (subsample of recording session)
sess_shortening_flag = behav_struct['was_sess_shortened'][0,0]
sess_last_trial = behav_struct['last_trial_considered'][0,0]
sess_time_win = behav_struct['time_win'][0,0]

# all behavior logs
behav_logs = behav_struct['logs'][0,0]

#%% neural data
neurons_data = scipy.io.loadmat(data_path / 'struct_spikes.mat')

# Extract the table
neurons_struct = neurons_data['struct_spikes']

# convert to pandas maybe

