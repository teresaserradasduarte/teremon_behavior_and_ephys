#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec, rcParams
import pathlib   
temp = pathlib.PosixPath   
pathlib.PosixPath = pathlib.WindowsPath

#%% 
# outputs saved to results_dir
results_dir = pathlib.WindowsPath(r"E:\18082023_Milka_StrCer_S4_g0\18082023_Milka_StrCer_S4_g0_imec0\ibl_sorter_results")
contam_pct = pd.read_csv(results_dir / 'cluster_ContamPct.tsv', sep='\t')['ContamPct'].values
chan_map =  np.load(results_dir / 'channel_map.npy')
chan_pos =  np.load(results_dir / 'channel_positions.npy')
templates =  np.load(results_dir / 'templates.npy')
chan_best = (templates**2).sum(axis=1).argmax(axis=-1)
chan_best = chan_map[chan_best]
amplitudes = np.load(results_dir / 'amplitudes.npy')
st = np.load(results_dir / 'spike_times.npy')
clu = np.load(results_dir / 'spike_clusters.npy')
tem = np.load(results_dir / 'spike_templates.npy')
pc_feat = np.load(results_dir / 'pc_features.npy')
pc_feat_ind = np.load(results_dir / 'pc_feature_ind.npy')
drift_time = np.load(results_dir / 'drift.times.npy')
drift_um = np.load(results_dir / 'drift.um.npy')
drift_depth = np.load(results_dir / 'drift_depths.um.npy')
firing_rates = np.unique(clu, return_counts=True)[1] * 30000 / st.max()
all_clu = np.unique(clu)

## get the depth of spikes by fiding the mean pos of the pc
sp_feat_ind = pc_feat_ind[tem,:];
ycoords = chan_pos[:,1];
sp_feat_y = ycoords[sp_feat_ind]; 
sp_depth = np.sum(sp_feat_y*pc_feat[:,1,:]**2,axis=1) / np.sum(pc_feat[:,1,:]**2,axis=1);

#%% plot figure
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
gray = .5 * np.ones(3)

fig = plt.figure(figsize=(15,8), dpi=100)
grid = gridspec.GridSpec(2, 4, figure=fig, hspace=0.5, wspace=0.5)

ax = fig.add_subplot(grid[0,0])
ax.plot(drift_time, drift_um);
ax.set_xlabel('time (sec)')
ax.set_ylabel('drift (um)')

ax = fig.add_subplot(grid[0,1:3])
t0 = 0
t1 = np.nonzero(st > 30000*5)[0][0]
ax.scatter(st/30000., chan_best[tem], s=0.00005, color='k', alpha=0.2)
#ax.set_xlim([0, 5])
ax.set_ylim([chan_map.max(), 0])
ax.invert_yaxis()
ax.set_xlabel('time (sec.)')
ax.set_ylabel('channel')
ax.set_title('spikes from units')

ax = fig.add_subplot(grid[1,0])
nb=ax.hist(firing_rates, 50, color=gray)
ax.set_xlim([0, 100])
ax.set_xlabel('firing rate (Hz)')
ax.set_ylabel('# of units')
ax.set_title('Total number of clusters: {}'.format(len(all_clu)))

ax = fig.add_subplot(grid[1,1])
nb=ax.hist(amplitudes, 100, color=gray)
ax.set_xlim([0, 50])
ax.set_xlabel('amplitude')
ax.set_ylabel('# of spikes')
ax.set_title('Total number of spikes: {}'.format(len(amplitudes)))

ax = fig.add_subplot(grid[1,2])
nb=ax.hist(np.minimum(100, contam_pct), np.arange(0,105,5), color=gray)
ax.plot([10, 10], [0, nb[0].max()], 'k--')
ax.set_xlabel('% contamination')
ax.set_ylabel('# of units')
ax.set_title('< 10% = good units')

ax = fig.add_subplot(grid[:,3])
nb=ax.hist(sp_depth, 70, color=gray, orientation='horizontal')
ax.set_xlabel('# of spikes')
ax.set_ylabel('depth of probe (um)')

plt.savefig(results_dir / 'ibsorter_output_overview.png')
plt.show()
