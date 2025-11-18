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
results_dir = pathlib.WindowsPath(r"E:\18082023_Milka_StrCer_S4_g0\18082023_Milka_StrCer_S4_g0_imec0\catGT\kilosort4")
stage = 'pre'


#%% Genrate overview figure
ops = np.load(results_dir / 'ops.npy', allow_pickle=True).item()
contam_pct = pd.read_csv(results_dir / 'cluster_ContamPct.tsv', sep='\t')['ContamPct'].values
chan_map =  np.load(results_dir / 'channel_map.npy')
templates =  np.load(results_dir / 'templates.npy')
chan_best = (templates**2).sum(axis=1).argmax(axis=-1)
chan_best = chan_map[chan_best]
amplitudes = np.load(results_dir / 'amplitudes.npy')
st = np.load(results_dir / 'spike_times.npy')
clu = np.load(results_dir / 'spike_clusters.npy')
tem = np.load(results_dir / 'spike_templates.npy')
pos = np.load(results_dir / 'spike_positions.npy')
firing_rates = np.unique(clu, return_counts=True)[1] * 30000 / st.max()
dshift = ops['dshift']
all_clu = np.unique(clu)

# Clusters types --> After curation
extragood_file = pd.read_csv(results_dir / 'cluster_Extragood.tsv', sep='\t')
clustergroup_file = pd.read_csv(results_dir / 'cluster_group.tsv', sep='\t')

extragood_clu = extragood_file[extragood_file['Extragood']=='y']['cluster_id'].values
good_clu = clustergroup_file[clustergroup_file['group']=='good']['cluster_id'].values
mua_clu = clustergroup_file[clustergroup_file['group']=='mua']['cluster_id'].values

idx_eg_sp = ba = np.where(np.isin(clu, extragood_clu))
idx_g_sp = ba = np.where(np.isin(clu, good_clu))
idx_m_sp = ba = np.where(np.isin(clu, mua_clu))

idx_eg_clu = ba = np.where(np.isin(all_clu, extragood_clu))
idx_g_clu = ba = np.where(np.isin(all_clu, good_clu))
idx_m_clu = ba = np.where(np.isin(all_clu, mua_clu))



#%% plot figure
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
gray = .5 * np.ones(3)

fig = plt.figure(figsize=(15,8), dpi=100)
grid = gridspec.GridSpec(2, 4, figure=fig, hspace=0.5, wspace=0.5)

ax = fig.add_subplot(grid[0,0])
ax.plot(np.arange(0, ops['Nbatches'])*2, dshift);
ax.set_xlabel('time (sec.)')
ax.set_ylabel('drift (um)')

ax = fig.add_subplot(grid[0,1:3])
t0 = 0
t1 = np.nonzero(st > ops['fs']*5)[0][0]
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
#nb_m = ax.hist(firing_rates[idx_m_clu], bins=bins_all, color='r', alpha = 0.10, label='mua')
#ax.set_xlim([0, 150])
ax.set_ylim([0, 80])
ax.set_xlabel('firing rate (Hz)')
ax.set_ylabel('# of units')
ax.set_title('Total number of clusters: {}'.format(len(all_clu)))

ax = fig.add_subplot(grid[1,1])
nb_all, bins_all, patches_all = ax.hist(amplitudes, bins=100, color=gray, alpha = 0.25, label='all')
nb_eg = ax.hist(amplitudes[idx_eg_sp], bins=bins_all, color='b', alpha = 0.25, label='extra good')
nb_eg = ax.hist(amplitudes[idx_g_sp], bins=bins_all, color='g', alpha = 0.25, label='good')
nb_eg = ax.hist(amplitudes[idx_m_sp], bins=bins_all, color='r', alpha = 0.2, label='mua')
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
nb_all, bins_all, patches_all = ax.hist(pos[:,1], 70, color=gray, alpha = 0.25, label='all', orientation='horizontal')
nb_eg = ax.hist(pos[:,1][idx_eg_sp], bins=bins_all, color='b', alpha = 0.25, label='extra good', orientation='horizontal')
nb_g = ax.hist(pos[:,1][idx_g_sp], bins=bins_all, color='g', alpha = 0.25, label='good', orientation='horizontal')
nb_m = ax.hist(pos[:,1][idx_m_sp], bins=bins_all, color='m', alpha = 0.25, label='mua', orientation='horizontal')
ax.set_xlabel('# of spikes')
ax.set_ylabel('depth of probe (um)')
ax.set_ylim([0, 4000])
plt.savefig(str(results_dir / ('output_phy_curated_' + stage + '.png')))
plt.show()


#%% plot figure - less stuff
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
gray = .5 * np.ones(3)

fig = plt.figure(figsize=(8,8), dpi=100)
grid = gridspec.GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.5)

ax = fig.add_subplot(grid[0,0])
ax.plot(np.arange(0, ops['Nbatches'])*2, dshift);
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
nb_all, bins_all, patches_all = ax.hist(pos[:,1], 70, color=gray, alpha = 0.25, label='all', orientation='horizontal')
nb_eg = ax.hist(pos[:,1][idx_eg_sp], bins=bins_all, color='b', alpha = 0.25, label='extra good', orientation='horizontal')
nb_g = ax.hist(pos[:,1][idx_g_sp], bins=bins_all, color='g', alpha = 0.25, label='good', orientation='horizontal')
nb_m = ax.hist(pos[:,1][idx_m_sp], bins=bins_all, color='m', alpha = 0.25, label='mua', orientation='horizontal')
ax.set_xlabel('# of spikes')
ax.set_ylabel('depth of probe (um)')
ax.set_ylim([0, 4000])
plt.savefig(str(results_dir / ('simple_output_phy_curated_' + stage + '.png')))
plt.show()
