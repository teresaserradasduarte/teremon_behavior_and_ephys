#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec, rcParams
import pathlib   
temp = pathlib.PosixPath   
pathlib.PosixPath = pathlib.WindowsPath

#%% CHANGE DIRECTORY AND STAGE (pre / post behavior alignment)
# outputs saved to results_dir
results_dir = pathlib.WindowsPath(r"E:\18082023_Milka_StrCer_S4_g0\18082023_Milka_StrCer_S4_g0_imec0\ibl_sorter_results")
stage = 'post'

#%% Genrate overview figure
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

# Clusters types --> After curation
extragood_file = pd.read_csv(results_dir / 'cluster_Extragood.tsv', sep='\t')
clustergroup_file = pd.read_csv(results_dir / 'cluster_group.tsv', sep='\t')

extragood_clu = extragood_file[extragood_file['Extragood']=='y']['cluster_id'].values
good_clu = clustergroup_file[clustergroup_file['group']=='good']['cluster_id'].values
mua_clu = clustergroup_file[clustergroup_file['group']=='mua']['cluster_id'].values

idx_eg_sp = np.where(np.isin(clu, extragood_clu))
idx_g_sp  = np.where(np.isin(clu, good_clu))
idx_m_sp  = np.where(np.isin(clu, mua_clu))

idx_eg_clu = np.where(np.isin(all_clu, extragood_clu))
idx_g_clu  = np.where(np.isin(all_clu, good_clu))
idx_m_clu  = np.where(np.isin(all_clu, mua_clu))



#%% plot figure
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
gray = .5 * np.ones(3)

fig = plt.figure(figsize=(15,8), dpi=100)
grid = gridspec.GridSpec(2, 4, figure=fig, hspace=0.5, wspace=0.5)

ax = fig.add_subplot(grid[0,0])
ax.plot(drift_time, drift_um);
ax.set_xlabel('time (sec.)')
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
nb_all, bins_all, patches_all = ax.hist(firing_rates, bins=50, color=gray, alpha = 0.25, label='all')
nb_eg = ax.hist(firing_rates[idx_eg_clu], bins=bins_all, color='b', alpha = 0.25, label='extra good')
nb_g = ax.hist(firing_rates[idx_g_clu], bins=bins_all, color='g', alpha = 0.25, label='good')
nb_m = ax.hist(firing_rates[idx_m_clu], bins=bins_all, color='r', alpha = 0.10, label='mua')
#ax.set_xlim([0, 150])
ax.set_ylim([0, 80])
ax.set_xlabel('firing rate (Hz)')
ax.set_ylabel('# of units')
ax.set_title('Total number of clusters: {}'.format(len(all_clu)))

ax = fig.add_subplot(grid[1,1])
nb_all, bins_all, patches_all = ax.hist(amplitudes, bins=100, color=gray, alpha = 0.25, label='all')
nb_eg = ax.hist(amplitudes[idx_eg_sp], bins=bins_all, color='b', alpha = 0.25, label='extra good')
nb_g = ax.hist(amplitudes[idx_g_sp], bins=bins_all, color='g', alpha = 0.25, label='good')
nb_m = ax.hist(amplitudes[idx_m_sp], bins=bins_all, color='r', alpha = 0.2, label='mua')
ax.set_xlim([7, 35])
#ax.set_ylim([0, 80])
ax.set_xlabel('amplitude')
ax.set_ylabel('# of spikes')
ax.set_title('Total number of spikes: {}'.format(len(amplitudes)))

ax = fig.add_subplot(grid[1,2])
ax.text(0.5, 0.8, f'# extra good units: {len(extragood_clu)}', ha='center', va='center',fontsize=14)
ax.text(0.5, 0.6, f'# good units: {len(good_clu)}', ha='center', va='center',fontsize=14)
ax.text(0.5, 0.4, f'# mua units: {len(mua_clu)}', ha='center', va='center',fontsize=14)
ax.text(0.5, 0.2, f'# found units: {len(all_clu)}', ha='center', va='center',fontsize=14)
ax.axis('off')

ax = fig.add_subplot(grid[:,3])
nb_all, bins_all, patches_all = ax.hist(sp_depth, 70, color=gray, alpha = 0.25, label='all', orientation='horizontal')
nb_eg = ax.hist(sp_depth[idx_eg_sp], bins=bins_all, color='b', alpha = 0.25, label='extra good', orientation='horizontal')
nb_g = ax.hist(sp_depth[idx_g_sp], bins=bins_all, color='g', alpha = 0.25, label='good', orientation='horizontal')
nb_m = ax.hist(sp_depth[idx_m_sp], bins=bins_all, color='m', alpha = 0.25, label='mua', orientation='horizontal')
ax.set_xlabel('# of spikes')
ax.set_ylabel('depth of probe (um)')
ax.set_ylim([0, 4000])
plt.savefig(str(results_dir / ('output_ibl_soerter_classification_' + stage + '.png')))
plt.show()

# %% plot figure
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
gray = .5 * np.ones(3)

fig = plt.figure(figsize=(8,8), dpi=100)
grid = gridspec.GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.5)

ax = fig.add_subplot(grid[0,0])
ax.plot(drift_time, drift_um);
ax.set_xlabel('time (sec.)')
ax.set_ylabel('drift (um)')


ax = fig.add_subplot(grid[1,0])
ax.text(0.5, 0.9, f'# extra good units: {len(extragood_clu)}', ha='center', va='center',fontsize=14)
ax.text(0.5, 0.7, f'# good units: {len(good_clu)}', ha='center', va='center',fontsize=14)
ax.text(0.5, 0.5, f'# mua units: {len(mua_clu)}', ha='center', va='center',fontsize=14)
ax.text(0.5, 0.3, f'# found units: {len(all_clu)}', ha='center', va='center',fontsize=14)
ax.text(0.5, 0.1, f'# spikes total: {len(amplitudes)}', ha='center', va='center',fontsize=14)
ax.axis('off')

ax = fig.add_subplot(grid[:,1])
nb_all, bins_all, patches_all = ax.hist(sp_depth, 70, color=gray, alpha = .25, label='all', orientation='horizontal')
#ax.hist(sp_depth[:], bins=bins_all, color=gray, alpha = .25, label='all', orientation='horizontal')
ax.hist(sp_depth[idx_eg_sp], bins=bins_all, color='b', alpha = 0.25, label='extra good', orientation='horizontal')
ax.hist(sp_depth[idx_g_sp], bins=bins_all, color='g', alpha = 0.25, label='good', orientation='horizontal')
ax.hist(sp_depth[idx_m_sp], bins=bins_all, color='m', alpha = 0.25, label='mua', orientation='horizontal')
ax.set_xlabel('# of spikes')
ax.set_ylabel('depth of probe (um)')
ax.set_ylim([0, 4000])
plt.savefig(str(results_dir / ('simple_output_ibl_sorter_classification_' + stage + '.png')))
plt.show()

# %%
