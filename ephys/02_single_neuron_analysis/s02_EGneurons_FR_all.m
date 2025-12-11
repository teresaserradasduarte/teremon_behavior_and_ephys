%% NEURONS STEP 2 
% GET EG neurons, calculate instantaneous FR
% align to push/pull and reach from DLC
% neuron type?? sort to visualize?
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
animal_idx = 5;
mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});

%% Session & ephys source
sess = 'R6';
imec_id = 1;
ephys_sess = '20082023_Ferrero_StrCer_S6_g0';
sorter_folder = 'catGT\kilosort4';
%sorter_folder = 'ibl_sorter_results';
ephys_local_folder = false;
flag_plot_rasters = false;

%% paths - load
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';    
%if ephys_local_folder==1, ephys_root = 'E:\';
%if ephys_local_folder==1, ephys_root = fullfile('Y:\shared-paton\teremon\ephys_curated\20230801_ChocolateGroup\',mouse);
if ephys_local_folder==1, ephys_root = fullfile('F:\curated_ephys_data\20230801_ChocolateGroup\',mouse);
else, ephys_root = fullfile(rootdir,"ephys_data",group,mouse); end
ephys_dir = fullfile(ephys_root,ephys_sess,strcat(ephys_sess,'_imec',num2str(imec_id)),sorter_folder);

bhv_dir = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess);

%% Load mat files
load(fullfile(bhv_dir,"behavior_fundamentals.mat"));
load(fullfile(ephys_dir,"neurons_session.mat"));

% Check brain region
region = neurons(1).meta.region;

% Save output
save_mat = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess);
save_out = fullfile(rootdir,"ephys_and_behavior","out_files",group,mouse,sess,region);
if ~exist(save_mat,"dir"), mkdir(save_mat); end
if ~exist(save_out,"dir"), mkdir(save_out); end

% Figure
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};
figProp.axeOpt = axeOpt;

%% Create neuron struct for extra-good neurons only
neu_all = struct2table(neurons);
eg_idx = find(neu_all.extraGood == 1);
n_eg = length(eg_idx);
n_good = length(find(neu_all.quality == 2));
n_mua = length(find(neu_all.quality == 1));


eg_neurons = struct();
for i=1:n_eg
    eg_neurons(i).st = neurons(eg_idx(i)).st_all;
    eg_neurons(i).phyID = neurons(eg_idx(i)).phyID;
    eg_neurons(i).unitIdx = neurons(eg_idx(i)).unitIdx;
    eg_neurons(i).depth = neurons(eg_idx(i)).depth;
    eg_neurons(i).egIdx = i;
    eg_neurons(i).CCG = neurons(eg_idx(i)).CCG;
    eg_neurons(i).CCG_bins = neurons(eg_idx(i)).CCG_bins;
    eg_neurons(i).templateWeight = neurons(eg_idx(i)).templateWeight;
    eg_neurons(i).templatePeakCh = neurons(eg_idx(i)).templatePeakCh;
    eg_neurons(i).templateWaveforms = neurons(eg_idx(i)).templateWaveforms;
    eg_neurons(i).templateID = neurons(eg_idx(i)).templateID;
end
eg_neurons(1).meta.mouse = neurons(1).meta.recordedAnimal;
eg_neurons(1).meta.session = neurons(1).meta.recordedSession;
eg_neurons(1).meta.region = neurons(1).meta.region;
eg_neurons(1).meta.n_good = n_good;
eg_neurons(1).meta.n_mua = n_mua;



%% Sanity check - reaches keep their identities
% Hit + location
left_hit_idx = find(ismember(bhv.reach_trials_inVec,bhv.left_idx) & bhv.hit_inVec==1);
center_hit_idx = find(ismember(bhv.reach_trials_inVec,bhv.center_idx) & bhv.hit_inVec==1);
right_hit_idx = find(ismember(bhv.reach_trials_inVec,bhv.right_idx) & bhv.hit_inVec==1);
hitIdx = find(bhv.hit_inVec==1);
hLRC_idx = zeros(length(hitIdx),1);
hLRC_idx(ismember(hitIdx,left_hit_idx))=1;
hLRC_idx(ismember(hitIdx,center_hit_idx))=2;
hLRC_idx(ismember(hitIdx,right_hit_idx))=3;

figure
plot(squeeze(bhv.reaches_inVec_px(:,1,left_hit_idx)),'color',bhv.figProp.clr_left); hold on
plot(squeeze(bhv.reaches_inVec_px(:,1,center_hit_idx)),'color',bhv.figProp.clr_center); hold on
plot(squeeze(bhv.reaches_inVec_px(:,1,right_hit_idx)),'color',bhv.figProp.clr_right); hold off
%plot(squeeze(reaches_inVec_px(:,1,cat_reach_inVec==1)))

%% Align to events params 
win_interest = [-3 3];
bin_width = .002;
sig_pk = .02;
k_gaus = gausskernel('sig',sig_pk,'binwidth',bin_width);

win_padded = win_interest + k_gaus.paddx;
bin_edges_padded = win_padded(1):bin_width:win_padded(2);
n_bins_padded = length(bin_edges_padded);
bin_edges = bin_edges_padded(1)-k_gaus.paddx(1):bin_width:bin_edges_padded(end)-k_gaus.paddx(2)+bin_width;
nr_bins = length(bin_edges);

eg_neu_FR_params.win_interest = win_interest;
eg_neu_FR_params.bin_width = bin_width;
eg_neu_FR_params.sig_pk = sig_pk;
eg_neu_FR_params.k_gaus = k_gaus;
eg_neu_FR_params.bin_edges_padded = bin_edges_padded;
eg_neu_FR_params.win_padded = win_padded;
eg_neu_FR_params.bin_edges = bin_edges;
eg_neu_FR_params.nr_bins = nr_bins;

figProp.spk_size = 6;
figProp.rcl_clr = [44 123 182]./256;
if ~isfield(bhv.figProp, 'clr_push')
    bhv.figProp.clr_push=[0.5156    0.6250    0.4844];
    bhv.figProp.clr_pull=[ 0.1016    0.3672    0.3867];
    save(fullfile(bhv_dir,"behavior_fundamentals.mat"),'bhv');
end


%% Align to events
event_tm_reach = bhv.reach_times_inVec;
event_tm_init = bhv.pp_times;
n_reach_inVec = length(event_tm_reach);
nr_pp_events = length(event_tm_init);
FR_reach = zeros(nr_bins,n_reach_inVec,n_eg);
FR_init = zeros(nr_bins,nr_pp_events,n_eg);

save_out_reach = fullfile(save_out,'reach_align');
if ~exist(save_out_reach,'dir'), mkdir(save_out_reach); end
save_out_init = fullfile(save_out,'init_align');
if ~exist(save_out_init,'dir'), mkdir(save_out_init); end

tic
fprintf('Calculating FR alinged to init and reach... \n')
for i = 1:n_eg
    % ALIGN TO REACH ----------------------------------------
     [eg_neurons,FR_reach_raw] = calculate_instant_FR_reach(eg_neurons,i,event_tm_reach,win_padded,bin_width);
     % smoothing individual trials
     FR_reach(:,:,i) = conv2(k_gaus.pdf,1,FR_reach_raw,"valid");
     eg_neurons(i).FR_reach = FR_reach(:,:,i);


    % ALIGN TO PUSH/PULL ------------------------------------
    [eg_neurons,FR_init_raw] = calculate_instant_FR_init(eg_neurons,i,event_tm_init,win_padded,bin_width);
    % smoothing individual trials
    FR_init(:,:,i) = conv2(k_gaus.pdf,1,FR_init_raw,"valid");
    eg_neurons(i).FR_init = FR_init(:,:,i);

    % FOR PLOTING - TRIALS IDENTITIES
    eg_neurons(i).reg = region;   
    pp_idx=bhv.valPushPull_invalPushPull_idx(bhv.init_invPush_invPull_idx==1);
    colors_pp=[bhv.figProp.clr_push;bhv.figProp.clr_pull];
    [~, ~, color_idx_pp] = unique(pp_idx); % convert to 1, 2, 3
    hLRC_idx(hLRC_idx==0)=2;
    colors_lcr=[bhv.figProp.clr_left;bhv.figProp.clr_center;bhv.figProp.clr_right];
    [~, ~, color_idx_lcr] = unique(hLRC_idx); % convert to 1, 2, 3

    % Plot
    if flag_plot_rasters == true
        figure(1)
        trialsIdx = find(bhv.hit_inVec==1);
        figProp.name = 'Hit reaches';
        figProp.ticks_y = bhv.reach_trials_inVec(trialsIdx);
        figProp.y_name = 'trial idx across session';
        plot_reach_raster_and_psth(eg_neurons,i,bin_edges,win_interest,figProp,trialsIdx); hold on
        scatter(ones(length(hLRC_idx),1)*-3, (1:length(hLRC_idx))', 36, colors_lcr(hLRC_idx, :), 'filled'); hold off
        saveas(gcf,fullfile(save_out_reach,['neu',num2str(eg_neurons(i).phyID),'_reach.png']),'png')

        figure(2)
        trialsIdx = find(bhv.init_invPush_invPull_idx==1); 
        figProp.name = 'Push-Pull valid';
        figProp.ticks_y = trialsIdx;
        figProp.y_name = 'trial idx across session';
        plot_init_raster_and_psth(eg_neurons,i,bin_edges,win_interest,figProp,trialsIdx); hold on
        scatter(ones(length(pp_idx),1)*-3, (1:length(pp_idx))', 36, colors_pp(color_idx_pp, :), 'filled');  hold off
        saveas(gcf,fullfile(save_out_init,['neu',num2str(eg_neurons(i).phyID),'_init.png']),'png');
    end
end
toc
fprintf('done \n')

%% Save
fprintf('Saving... \n')
tic
save(fullfile(save_mat,strcat('eg_neurons_',region,'.mat')),'eg_neurons',...
    'figProp','eg_neu_FR_params','-v7.3');
toc
fprintf('done!! \n')






