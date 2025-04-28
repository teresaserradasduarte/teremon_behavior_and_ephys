%% Figure neuron types
clear; close all; clc

%% %% Load data
% Group and individual
group = '20230801_ChocolateGroup';
group_behav = '20230511_ChocolateGroup';
setup = 'headfixed_dynamicTarget';
animals = {...
    'CoteDor',...
    'Lindt',...
    'Toblerone',...
    'Milka',...
    'FerreroRocher'};
animal_idx = 4;
mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});

%% Session & ephys source
sess = 'R4';
region = 'CB';
ephys_sess = '18082023_Milka_StrCer_S4_g0';
sorter_folder = 'catGT\kilosort4';
ephys_local_folder = 1;

%% paths
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
behavior_dir = fullfile(rootdir,"behavior_data","raw_data",group_behav,setup,mouse,sess);
bhv_dir = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess);
mat_dir = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess,region);
save_out = fullfile(rootdir,"ephys_and_behavior","out_files",group,mouse,sess,region,'figure_cellTypes');
if ~exist(save_out,"dir"), mkdir(save_out); end

% Load mat files
load(fullfile(bhv_dir,"behavior_fundamentals.mat"));
load(fullfile(behavior_dir,"behavior_session.mat"));
load(fullfile(mat_dir,"eg_neurons.mat"));

% Figure
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};
figProp.axeOpt = axeOpt;
push_clr = behavior.colors.push_clr;
pull_clr =  behavior.colors.pull_clr;
pp_clr = [63 130 109]./256;
right_clr = behavior.colors.right_color;
center_clr = behavior.colors.center_color;
left_clr = behavior.colors.left_color;
rcl_clr = [44 123 182]./256;

% Window
spk_bins = eg_neu_FR_params.bin_edges;
tm_before = eg_neu_FR_params.win_interest(1);
tm_after = eg_neu_FR_params.win_interest(2);

%% Select neurons
un = 6;

%% Push and pull mean
FR_push = eg_neurons(un).FR_init(:,bhv.push_idx);
FR_pull = eg_neurons(un).FR_init(:,bhv.pull_idx);

psth_push = mean(FR_push,2,'omitnan');
psth_pull = mean(FR_pull,2,'omitnan');
psth_push_sem = std(FR_push,[],2)./length(bhv.push_idx);
psth_pull_sem = std(FR_pull,[],2)./length(bhv.pull_idx);


% Left / center / right
left_reaches_idx = find(ismember(bhv.reach_trials_inVec,bhv.left_idx) & bhv.hit_inVec);
center_reaches_idx = find(ismember(bhv.reach_trials_inVec,bhv.center_idx) & bhv.hit_inVec);
right_reaches_idx = find(ismember(bhv.reach_trials_inVec,bhv.right_idx) & bhv.hit_inVec);

FR_left = eg_neurons(un).FR_reach(:,left_reaches_idx);
FR_center = eg_neurons(un).FR_reach(:,center_reaches_idx);
FR_right = eg_neurons(un).FR_reach(:,right_reaches_idx);

psth_left = mean(FR_left,2,'omitnan');
psth_center = mean(FR_center,2,'omitnan');
psth_right = mean(FR_right,2,'omitnan');
psth_left_sem = std(FR_left,[],2)./length(left_reaches_idx);
psth_center_sem = std(FR_center,[],2)./length(center_reaches_idx);
psth_right_sem = std(FR_right,[],2)./length(right_reaches_idx);

% Raster 
line_x = [tm_before;tm_before]+0.08;
% init
pull_bounds = behavior.init.pull_bounds;
push_bounds = behavior.init.push_bounds;
%reach
left_bounds =  behavior.reach.left_bounds;
right_bounds = behavior.reach.right_bounds;
center_bounds = behavior.reach.center_bounds;

trials_w_reach = bhv.reach_trials_inVec(bhv.hit_inVec);
spike_times_init = eg_neurons(un).st_init(bhv.init_invPush_invPull_idx==1)';
trials_idx_scat_init = cell2mat(arrayfun(@(i) i * ones(1, numel(spike_times_init{i})),...
    1:length(spike_times_init), 'UniformOutput', false));
reach_in_init = bhv.reach_times_inVec(bhv.hit_inVec)-bhv.pp_times(trials_w_reach);
nr_trials = length(spike_times_init);
trials_vec = 1:nr_trials;

spike_times_reach = eg_neurons(un).st_reach(bhv.hit_inVec)';
trials_idx_scat_reach = cell2mat(arrayfun(@(i) i * ones(1, numel(spike_times_reach{i})),...
     1:length(spike_times_reach), 'UniformOutput', false));
init_in_reach = bhv.pp_times(trials_w_reach)-bhv.reach_times_inVec(bhv.hit_inVec);


%% FIGURE WITHH LOTS OF STUFF 
tm_before = -2;
tm_after = 2;
tt = tiledlayout(3,7);

% DEPTH -----------------------------------------------------------
axx0 = nexttile(1);
axx0.Layout.TileSpan = [3,1];

[depth_ocupacy,depth]=histcounts(depths_all,50);
histogram('BinEdges', depth ,'BinCounts',  depth_ocupacy,...
    'facecolor',figProp.FR_clrFace ,'facealpha',0.6,'edgecolor','none');
set(gca,figProp.axeOpt{:},'XDir','reverse')
camroll(-90)
h = gca;
xlim([0 3500])
h.YAxis.Visible = 'off';
xline(eg_neurons(un).depth,'-',sprintf('%s%i','Unit ID: ',eg_neurons(un).phyID)...
    ,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle',...
    'color','k','linewidth',3)
xlabel('Depth of neuropixel probe (\mum)');

% INIT PSTH --------------------------------------------------------
axx1 = nexttile;
axx1.Layout.TileSpan = [1,2];

shadedErrorBar(spk_bins,psth_push,psth_push_sem,'lineProps',{'Color',push_clr,'LineWidth',2},'transparent',1,'patchSaturation',0.1); hold on
shadedErrorBar(spk_bins,psth_pull,psth_pull_sem,'lineProps',{'Color',pull_clr,'LineWidth',2},'transparent',1,'patchSaturation',0.1); hold on
ylabel('spike rate (sp/s)'); xlabel('time from trial init (s)');
xlim([tm_before tm_after])
xline(0,'--','Color',[pp_clr 0.3],'linewidth',2);
set(gca,axeOpt{:})
title('trial initiation (push/pull)')
%title(sprintf('%s%s','Region: ',region));

% REACH PSTH --------------------------------------------------------
axx2 = nexttile;
axx2.Layout.TileSpan = [1,2];

shadedErrorBar(spk_bins,psth_right,psth_right_sem,'lineProps',{'Color',right_clr,'LineWidth',2},'transparent',1,'patchSaturation',0.1); hold on
shadedErrorBar(spk_bins,psth_left,psth_left_sem,'lineProps',{'Color',left_clr,'LineWidth',2},'transparent',1,'patchSaturation',0.1); hold on
shadedErrorBar(spk_bins,psth_center,psth_center_sem,'lineProps',{'Color',center_clr,'LineWidth',2},'transparent',1,'patchSaturation',0.1); hold on
ylabel('spike rate (sp/s)'); xlabel('time from reach (s)');
xlim([tm_before tm_after])
xline(0,'--','Color',[rcl_clr 0.3],'linewidth',2);
set(gca,axeOpt{:})
title('water collection (reach)')
%title(sprintf('%s%s','Cell unit: ',num2str(neurons(un).phyID)));


% AUTO-CORRELOGRAM -------------------------------------------------
axxx1 = nexttile;
axxx1.Layout.TileSpan = [1,2];
bar(eg_neurons(un).CCG_bins*1000, eg_neurons(un).CCG,'BarWidth', 1,...
    'facecolor',figProp.spk_clr_hist,'edgecolor','none','facealpha',.9)
set(gca,figProp.axeOpt{:})
xlabel('time (ms)'); ylabel('spike rate (sp/s)')
title('auto-correlogram')

% RASTER INIT -----------------------------------------------------
axx2 = nexttile;
axx2.Layout.TileSpan = [2,2];
spk_size = 2;
event_size = 5;
patch_wd = 10;
scatter(cell2mat(spike_times_init),trials_idx_scat_init,spk_size,'k.'), hold on
scatter(reach_in_init,trials_vec(trials_w_reach),event_size,'filled','MarkerFaceColor',rcl_clr)
scatter(zeros(size(trials_vec)),trials_vec,event_size,'filled','MarkerFaceColor',pp_clr)
plot(line_x,pull_bounds,'linewidth',patch_wd,'color',pull_clr)
plot(line_x,push_bounds,'linewidth',patch_wd,'color',push_clr)
hold off
axis tight;
xlabel('time from trial init (s)'); ylabel('trials in session')
axis([tm_before tm_after 0 nr_trials])
set(gca,axeOpt{:})
%set(gcf,'Position',[2582,215,560,771],'color','w')

% RASTER REACH -----------------------------------------------------
axx2 = nexttile;
axx2.Layout.TileSpan = [2,2];
spk_size = 2;
event_size = 5;
patch_wd = 10;
scatter(cell2mat(spike_times_reach),trials_idx_scat_reach,spk_size,'k.'), hold on
scatter(init_in_reach,1:length(spike_times_reach),event_size,'filled','MarkerFaceColor',pp_clr)
scatter(zeros(size(trials_vec)),trials_vec,event_size,'filled','MarkerFaceColor',rcl_clr)
hold off
axis tight;
xlabel('time from trial init (s)'); ylabel('trials in session')
axis([tm_before tm_after 0 nr_trials])
set(gca,axeOpt{:})
%set(gcf,'Position',[2582,215,560,771],'color','w')

yticks_to_label = round(linspace(1, length(spike_times_reach), 10));
ytick_labels = trials_w_reach(yticks_to_label);
ylim([1 length(spike_times_reach)])
yticks(yticks_to_label);
yticklabels(string(ytick_labels));
ylabel('trials in session')
hold on
plot(line_x,left_bounds,'linewidth',patch_wd,'color',left_clr)
plot(line_x,right_bounds,'linewidth',patch_wd,'color',right_clr)
plot(line_x,center_bounds,'linewidth',patch_wd,'color',center_clr)
hold off

% --------------------------------------------------------------
% WAVWFORMS
axxx3 = nexttile;
axxx3.Layout.TileSpan = [2,2];
sample_rate = 30000;
time_wf = (0:1/sample_rate:(size(eg_neurons(un).templateWaveforms,1)-1)/sample_rate)*1000;
% plot(time_wf,zero2nan(eg_neurons(un).templateWaveforms(:,:)),'linewidth',1.5,'Color',[.8 .8 .8]);
% hold on
for i = 1:length(eg_neurons(un).templateID)
    plot(time_wf,zero2nan(eg_neurons(un).templateWaveforms(:,eg_neurons(un).templatePeakCh(i))),...
        'linewidth',2,'Color',[.1 .1 .1 eg_neurons(un).templateWeight(i)]);
end
hold off
set(gca,figProp.axeOpt{:})
xlabel('time (ms)'); ylabel('template waveform (pseudo-volts)')
str = 'Straight Line Plot from 1 to 10';
dim_ch=[0.81 0.12 0.1 0.1];
title('template waveform')
annotation('textbox',dim_ch,'String',...
    sprintf('%s%s','channel = ',num2str(eg_neurons(un).templatePeakCh)),...
    'edgecolor','w','fontsize',10.5)

    title(tt,sprintf('%s%s%s%s%s%s%s',...
        'region: ',eg_neurons(un).reg,...
        ' | phyID: ', num2str(eg_neurons(un).phyID),...
        ' | idx: ', num2str(eg_neurons(un).unitIdx),...
        ' | idxEG: ',num2str(eg_neurons(un).egIdx)),...
        'Interpreter', 'none','fontsize',12,'fontweight','normal')
set(gcf,'Position',[2164 377 1366 557],'color','w')

 saveas(gcf,fullfile(save_out,['neu',num2str(eg_neurons(un).phyID),'.png']),'png');
% print(gcf, fullfile(save_out, ['neu',num2str(eg_neurons(un).phyID),'.pdf']), '-dpdf', '-painters');



%% CELL TYPE FIGURE
un = 6

%select_uns = [6, 18, 19, 23, 24, 38];

%for un = select_uns
figure()
subplot(211)
bar(eg_neurons(un).CCG_bins*1000, eg_neurons(un).CCG,'BarWidth', 1,...
    'facecolor',figProp.spk_clr_hist,'edgecolor','none','facealpha',.9)
set(gca,figProp.axeOpt{:})
xlabel('time (ms)'); ylabel('spike rate (sp/s)')

title(sprintf('%s%s%s%s%s%s%s',...
    'region: ',eg_neurons(un).reg,...
    ' | phyID: ', num2str(eg_neurons(un).phyID),...
    ' | idxEG: ',num2str(eg_neurons(un).egIdx)),...
    'Interpreter', 'none','fontsize',10,'fontweight','normal');

inset_pos = [0.75 0.8 0.15 0.1];  % Top right corner
inset_ax = axes('Position', inset_pos);
plot(time_wf, zero2nan(eg_neurons(un).templateWaveforms(:, eg_neurons(un).templatePeakCh(1))), ...
    'LineWidth', 2, 'Color', figProp.spk_clr_hist)
axis tight
box off
axis off
set(inset_ax, 'xtick', [], 'ytick', [], 'color', 'none')  % Optional: clean look

subplot(212)

spike_times_reach = eg_neurons(un).st_reach(bhv.hit_inVec)';
trials_idx_scat_reach = cell2mat(arrayfun(@(i) i * ones(1, numel(spike_times_reach{i})),...
     1:length(spike_times_reach), 'UniformOutput', false));

spk_size = .01;
event_size = 2;
patch_wd = 10;
scatter(cell2mat(spike_times_reach),trials_idx_scat_reach,spk_size,'k.'), hold on
scatter(zeros(size(trials_vec)),trials_vec,event_size,'filled','MarkerFaceColor',rcl_clr)
hold off
axis tight;
xlabel('time from trial init (s)'); ylabel('trials in session')
axis([tm_before tm_after 0 nr_trials])
set(gca,axeOpt{:})
%set(gcf,'Position',[2582,215,560,771],'color','w')
hold off
axis tight;
xlabel('time from reach (s)'); ylabel('trials in session')
axis([tm_before tm_after 0 nr_trials])
set(gca,axeOpt{:})
%set(gcf,'Position',[2582,215,560,771],'color','w')

yticks_to_label = round(linspace(1, length(spike_times_reach), 3));
ytick_labels = trials_w_reach(yticks_to_label);
ylim([1 length(spike_times_reach)])
yticks(yticks_to_label);
yticklabels(string(ytick_labels));
ylabel('trials in session')

%set(gcf,'position',[2518         176         729         704],'color','none');
set(gcf,'position',[2518         588         262         292],'color','none');

%saveas(gcf,fullfile(save_out,['cellType_neu',num2str(eg_neurons(un).phyID),'.png']),'BackgroundColor', 'none','png');
print(gcf,fullfile(save_out,['cellType_neu',num2str(eg_neurons(un).phyID),'.png']), '-dpng', '-r300');
print(gcf, fullfile(save_out, ['cellType_neu',num2str(eg_neurons(un).phyID),'.pdf']), '-dpdf', '-painters');
%end

%% Depth fig
%select_uns=[6, 23, 27];
figure()
[depth_ocupacy,depth]=histcounts(depths_all,50);
histogram('BinEdges', depth ,'BinCounts',  depth_ocupacy,...
    'facecolor',figProp.FR_clrFace ,'facealpha',0.6,'edgecolor','none');
set(gca,figProp.axeOpt{:},'XDir','reverse')
camroll(-90)
h = gca;
xlim([0 3500])
h.YAxis.Visible = 'off';
for un = select_uns
xline(eg_neurons(un).depth,'-',sprintf('%s%i','Unit ID: ',eg_neurons(un).phyID)...
    ,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle',...
    'color','k','linewidth',3)
end
xlabel('Depth of neuropixel probe (\mum)');
set(gcf,'position',[2286         272         135         561],'color','w');

saveas(gcf,fullfile(save_out,'depth_cellType_neu.png'),'png');
print(gcf, fullfile(save_out, 'depth_cellType_neu.pdf'), '-dpdf', '-painters');

