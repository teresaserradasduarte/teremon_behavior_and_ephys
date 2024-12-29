clear; close all; clc

%% Load behavior
behavior_root ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
ephys_root = 'E:'; %group_ephys = '20230801_ChocolateGroup';
group_setup_behav = strcat('20230511_ChocolateGroup',filesep,'headfixed_dynamicTarget');
mouse = '4_Milka';
paw_pref = 'right';
session = 'R4';
imec_id = '0';
ephys_sess = '18082023_Milka_StrCer_S4_g0';
sorter_folder = 'catGT\kilosort4';

% Load behavior
behavior_path = strcat(behavior_root,filesep,group_setup_behav,filesep,mouse,filesep,session);
load(strcat(behavior_path,filesep,'behavior_session.mat'));

% Load neurons
path_neurons  = strcat(ephys_root,filesep,ephys_sess,filesep,ephys_sess,'_imec',num2str(imec_id),filesep,sorter_folder);
load(strcat(path_neurons,filesep,'neurons_session.mat'));

save_path = strcat(path_neurons,filesep,'CS_neurons');
if ~exist("save_path","dir"), mkdir(save_path); end

%% Load neuron of interest
idx_CSs = [neurons.potentialCS] == 1;
CS_neurons = neurons(idx_CSs);
nr_cs_uns = length(CS_neurons);

cs_un = 86;
ss_un = 89;
binsize_CC = 0.001;
duration_CC = 0.3;
% Cross correlogram -------------------------------
[ccg, ccg_t] = CCGBz([double(neurons(ss_un).st); double(neurons(cs_un).st)], [ones(size(neurons(ss_un).st, 1), 1); ...
    ones(size(neurons(cs_un).st, 1), 1) * 2], 'binSize', binsize_CC, 'duration', duration_CC, 'norm', 'rate'); %function

%%

figure
ff=tiledlayout(2,2);
nexttile
%subplot(221)
bar(ccg_t.*1000, ccg(:,1,1),'BarWidth', 1,...
    'facecolor',figProp.spk_clr_hist,'edgecolor','none','facealpha',.9)
xlabel('time lag (ms)'); ylabel('spike rate (sp/s)')
title('Auto-correlogram: SS with respect to SS')
set(gca,figProp.axeOpt{:})

nexttile
%subplot(222)
bar(ccg_t.*1000, ccg(:,1,2),'BarWidth', 1,...
    'facecolor',figProp.spk_clr_hist,'edgecolor','none','facealpha',.9)
xlabel('time lag (ms)'); ylabel('spike rate (sp/s)')
title('Cross-correlogram: CS with respect to SS')

set(gca,figProp.axeOpt{:})
nexttile
%subplot(223)
bar(ccg_t.*1000, ccg(:,2,1),'BarWidth', 1,...
    'facecolor',figProp.spk_clr_hist,'edgecolor','none','facealpha',.9)
set(gca,figProp.axeOpt{:})
xlabel('time lag (ms)'); ylabel('spike rate (sp/s)')
title('Cross-correlogram: SS with respect to CS')

nexttile
%subplot(224)
bar(ccg_t.*1000, ccg(:,2,2),'BarWidth', 1,...
    'facecolor',figProp.spk_clr_hist,'edgecolor','none','facealpha',.9)
set(gca,figProp.axeOpt{:})
xlabel('time lag (ms)'); ylabel('spike rate (sp/s)')
title('Auto-correlogram: CS with respect to CS')

title(ff, sprintf('%s%i%s%i', 'SS phyID = ',neurons(ss_un).phyID, ' | CS phyID = ' ,neurons(cs_un).phyID))

set(gcf,'position',[2269 129 1177 818],'color','w')

saveas(gcf, strcat(save_path,filesep,'SS',num2str(neurons(ss_un).phyID),'_CS',num2str(neurons(cs_un).phyID),'_pair.png'),'png')







