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
results_dir = pathlib.WindowsPath(r"E:\20082023_Ferrero_StrCer_S6_g0\20082023_Ferrero_StrCer_S6_g0_imec0\catGT\kilosort4")
ops_file = results_dir / 'ops.npy'
ops = np.load(ops_file, allow_pickle=True).item()
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
ax.set_title('Total number of clusters: {}'.format(len(amplitudes)))

ax = fig.add_subplot(grid[1,2])
nb=ax.hist(np.minimum(100, contam_pct), np.arange(0,105,5), color=gray)
ax.plot([10, 10], [0, nb[0].max()], 'k--')
ax.set_xlabel('% contamination')
ax.set_ylabel('# of units')
ax.set_title('< 10% = good units')

ax = fig.add_subplot(grid[:,3])
nb=ax.hist(pos[:,1], 70, color=gray, orientation='horizontal')
ax.set_xlabel('# of spikes')
ax.set_ylabel('depth of probe (um)')

plt.savefig(results_dir / 'ks_output_overview.png')
plt.show()
