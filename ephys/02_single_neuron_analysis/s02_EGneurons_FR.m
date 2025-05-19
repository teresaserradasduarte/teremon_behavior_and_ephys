%% GET EG neurons, calculate instantaneous FR
% align to push/pull and reach from DLC
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
animal_idx = 3;
mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});

%% Session & ephys source
sess = 'R6';
imec_id = 1;
ephys_sess = '28082023_Toblerone_StrCer_S6_g0';
sorter_folder = 'catGT\kilosort4';
ephys_local_folder = false;
flag_plot_rasters = true;
flag_save_behavior_funda = true;

%% Depth cortex 
depth_cortex_lim = 1; %1000; %65.1348;
depthLim_phyID = nan;


%% paths
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
if ephys_local_folder==1, ephys_root = 'E:\';
else, ephys_root = fullfile(rootdir,"ephys_data",group,mouse); end
behavior_dir = fullfile(rootdir,"behavior_data","raw_data",group_behav,setup,mouse,sess);
bhv_dir = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess);
reaching_dir = fullfile(rootdir,"behavior_data","analyzed_data","mat_files",group_behav,setup,mouse,sess);
ephys_dir = fullfile(ephys_root,ephys_sess,strcat(ephys_sess,'_imec',num2str(imec_id)),sorter_folder);

% Load mat files
load(fullfile(behavior_dir,"behavior_session.mat"));
load(fullfile(reaching_dir,"session_reaching_data_paw.mat"));
load(fullfile(ephys_dir,"neurons_session.mat"));

% Check brain region
paw_pref = mouse_info.paw_pref;
if (strcmp(paw_pref,'R') || strcmp(paw_pref,'right'))
    if imec_id == 0, region = 'CB';
    elseif imec_id == 1, region = 'BG';
    end
elseif  (strcmp(paw_pref,'L') || strcmp(paw_pref,'left'))
    if imec_id == 0, region = 'BG';
    elseif imec_id == 1, region = 'CB';
    end
end

% Save output
save_mat = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess,region);
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
depths_all = neu_all.depth;

if isnan(depth_cortex_lim)
    limit_depth_neuron = find(neu_all.phyID == depthLim_phyID);
    depth_lim = neurons(limit_depth_neuron).depth;
else
    depth_lim = depth_cortex_lim;
end

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
    if strcmp(region,"BG")
        if (eg_neurons(i).depth <= depth_lim)
            eg_neurons(i).reg = 'BG ';
            eg_neurons(i).isBG = 1;
        elseif (eg_neurons(i).depth > depth_lim)
            eg_neurons(i).reg = 'CTX';
            eg_neurons(i).isBG = 0;
        end
    elseif strcmp(region,"CB")
        if (eg_neurons(i).depth <= depth_lim)
            eg_neurons(i).reg = 'CB_DCN ';
            eg_neurons(i).isDCN = 1;
        elseif (eg_neurons(i).depth > depth_lim)
            eg_neurons(i).reg = 'CB_CTX';
            eg_neurons(i).isDCN = 0;    
        end
    end
end

neu_eg = struct2table(eg_neurons);
if strcmp(region,"BG")
    idx_BG = find(neu_eg.depth <= depth_lim);
    idx_CT = find(neu_eg.depth > depth_lim);
elseif strcmp(region,"CB")
    idx_CB = 1:n_eg;
    idx_CB_DCN = find(neu_eg.depth <= depth_lim);
    idx_CB_CTX = find(neu_eg.depth > depth_lim);
end
clear neu_all neurons;

%% Get event time from behavior
% Win of good behavior + ehpys
behav_start = behavior.behavior_duration.time_start;
behav_stop = behavior.behavior_duration.time_end + behav_start;
nr_trials = behavior.behavior_duration.trial_end;
trials_vec = 1:nr_trials;

% Trial identity defined by logs
push_idx = behavior.init.idx_trial_push(behavior.init.idx_trial_push<=nr_trials);
pull_idx = behavior.init.idx_trial_pull(behavior.init.idx_trial_pull<=nr_trials);
left_idx = behavior.reach.left_idx(behavior.reach.left_idx<=nr_trials);
right_idx = behavior.reach.right_idx(behavior.reach.right_idx<=nr_trials);
center_idx = behavior.reach.center_idx(behavior.reach.center_idx<=nr_trials);

% Convert L/R to Dom/non-Dom
if strcmp(mouse_info.paw_pref,'R')
    dom_idx = right_idx;
    nondom_idx = left_idx;
elseif strcmp(mouse_info.paw_pref,'L')
    dom_idx = left_idx;
    nondom_idx = right_idx;
end

% Event times from logs - trial init
initiation_times_all = behavior.inputs.read_log(behavior.logs.trial_init_ind,1);
inva_push_times_all = behavior.inputs.read_log(behavior.logs.inval_push_ind,1);
inva_pull_times_all = behavior.inputs.read_log(behavior.logs.inval_pull_ind,1);
time_sess_end = behavior.behavior_duration.time_start+behavior.behavior_duration.time_end+10;
initiation_times = initiation_times_all(initiation_times_all<time_sess_end);
inva_push_times = inva_push_times_all(inva_push_times_all<time_sess_end);
inva_pull_times = inva_pull_times_all(inva_pull_times_all<time_sess_end);
pp_times = cat(1,initiation_times,inva_push_times,inva_pull_times);
init_invPush_invPull_idx = cat(1,ones(length(initiation_times),1),...
    ones(length(inva_push_times),1)*2,ones(length(inva_pull_times),1)*3);
nr_pp_events = length(pp_times);
valPushPull_invalPushPull_idx = init_invPush_invPull_idx;
valPushPull_invalPushPull_idx(push_idx)=0;

% event times from reaches
reach_times = reaches.reach_timestamps_mat(:,reaches.reach_params.reach_interval.max_reach);
n_reach = length(reach_times);

% reach identities
reach_trials = reaches.reach_trial-1;
is_success = reaches.success_reach==1;
is_hit = reaches.hit_reach==1;
is_purpose = reaches.purpose_reach==1;
cat_reach = reaches.cat_reach;

% Idx of reaches within win
reach_inVec_idx = ismember(reach_trials,trials_vec);
reach_times_inVec = reach_times(reach_inVec_idx);
n_reach_inVec = length(reach_times_inVec);
reach_trials_inVec = reach_trials(reach_inVec_idx);
cat_reach_inVec = cat_reach(reach_inVec_idx);
reaches_inVec_px = reaches.reach_mat(:,:,reach_inVec_idx);
hit_inVec = is_hit(reach_inVec_idx);

% Save structure with fundamental behavior info
if flag_save_behavior_funda == true
    bhv.push_idx = push_idx;
    bhv.pull_idx = pull_idx;
    bhv.left_idx = left_idx;
    bhv.center_idx = center_idx;
    bhv.right_idx = right_idx;
    bhv.dom_idx = dom_idx;
    bhv.nondom_idx = nondom_idx;
    bhv.paw_pref = mouse_info.paw_pref;
    bhv.pp_times = pp_times;
    bhv.init_invPush_invPull_idx = init_invPush_invPull_idx;
    bhv.valPushPull_invalPushPull_idx = valPushPull_invalPushPull_idx;
    bhv.reach_inVec_idx = reach_inVec_idx;
    bhv.reach_times_inVec = reach_times_inVec;
    bhv.reach_trials_inVec = reach_trials_inVec;
    bhv.cat_reach_inVec = cat_reach_inVec;
    bhv.reaches_inVec_px = reaches_inVec_px;
    bhv.hit_inVec = hit_inVec;
    bhv.figProp = figProp;
    save(fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess,'behavior_fundamentals.mat'),'bhv');
end



%% Sanity check - reaches keep their identities
% Hit + location
left_hit_idx = find(ismember(bhv.reach_trials_inVec,bhv.left_idx) & bhv.hit_inVec==1);
center_hit_idx = find(ismember(bhv.reach_trials_inVec,bhv.center_idx) & bhv.hit_inVec==1);
right_hit_idx = find(ismember(bhv.reach_trials_inVec,bhv.right_idx) & bhv.hit_inVec==1);
hitIdx = find(hit_inVec==1);
hLRC_idx = zeros(length(hitIdx),1);
hLRC_idx(ismember(hitIdx,left_hit_idx))=1;
hLRC_idx(ismember(hitIdx,center_hit_idx))=2;
hLRC_idx(ismember(hitIdx,right_hit_idx))=3;

figure
plot(squeeze(reaches_inVec_px(:,1,left_hit_idx)),'color',behavior.colors.left_color); hold on
plot(squeeze(reaches_inVec_px(:,1,center_hit_idx)),'color',behavior.colors.center_color); hold on
plot(squeeze(reaches_inVec_px(:,1,right_hit_idx)),'color',behavior.colors.right_color); hold off
%plot(squeeze(reaches_inVec_px(:,1,cat_reach_inVec==1)))

%% Align to events
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

%% Align to events
event_tm_reach = reach_times_inVec;
event_tm_init = pp_times;
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

    % Plot
    if flag_plot_rasters == true
        figure(1)
        trialsIdx = find(hit_inVec==1);
        figProp.name = 'Hit reaches';
        figProp.ticks_y = reach_trials_inVec(trialsIdx);
        figProp.y_name = 'trial idx across session';
        plot_reach_raster_and_psth(eg_neurons,i,bin_edges,win_interest,figProp,trialsIdx); hold on
        colors_lcr=[behavior.colors.left_color;behavior.colors.center_color;behavior.colors.right_color];
        [~, ~, color_idx_lcr] = unique(hLRC_idx); % convert to 1, 2, 3
        scatter(ones(length(hLRC_idx),1)*-2.9, (1:length(hLRC_idx))', 36, colors_lcr(color_idx_lcr, :), 'filled'); hold off      
        saveas(gcf,fullfile(save_out_reach,['neu',num2str(eg_neurons(i).phyID),'_reach.png']),'png')


        figure(2)
        trialsIdx = find(init_invPush_invPull_idx==1); 
        pp_idx=valPushPull_invalPushPull_idx(init_invPush_invPull_idx==1);
        figProp.name = 'Push-Pull valid';
        figProp.ticks_y = trialsIdx;
        figProp.y_name = 'trial idx across session';
        plot_init_raster_and_psth(eg_neurons,i,bin_edges,win_interest,figProp,trialsIdx); hold on
        colors_pp=[behavior.colors.push_clr;behavior.colors.pull_clr];
        [~, ~, color_idx_pp] = unique(pp_idx); % convert to 1, 2, 3
        scatter(ones(length(pp_idx),1)*-2.9, (1:length(pp_idx))', 36, colors_pp(color_idx_pp, :), 'filled');  hold off
        saveas(gcf,fullfile(save_out_init,['neu',num2str(eg_neurons(i).phyID),'_init.png']),'png');
    end
end
toc
fprintf('done \n')

%% Save
fprintf('Saving... \n')
tic
if strcmp(region,"BG")
save(fullfile(save_mat,'eg_neurons.mat'),'eg_neurons','depths_all',...
    'FR_init','FR_reach','figProp','eg_neu_FR_params',...
    'depth_lim','idx_BG','idx_CT','-v7.3');
elseif strcmp(region,"CB")
    save(fullfile(save_mat,'eg_neurons.mat'),'eg_neurons','depths_all',...
        'FR_init','FR_reach','figProp','eg_neu_FR_params',...
        'depth_lim','idx_CB','idx_CB_DCN','idx_CB_CTX','-v7.3');
end
toc
fprintf('done!! \n')






