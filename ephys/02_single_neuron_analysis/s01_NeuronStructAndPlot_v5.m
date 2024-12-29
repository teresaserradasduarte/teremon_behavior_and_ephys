%% Get neuron struct
clear; close all; clc

%% Manage paths
person = 'teresa';
%person = 'simon';
if strcmp(person,'teresa')
    behavior_root ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
    ephys_root = 'E:'; %group_ephys = '20230801_ChocolateGroup';
    %ephys_root = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_data\20230801_ChocolateGroup';
    %ephys_root = 'G:\ePhys\'; %group_ephys = '20230801_ChocolateGroup';
elseif strcmp(person,'simon')
    behavior_root ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
    ephys_root = 'C:\Users\SimonZ\Documents\Data\EPhys\DTATA_2CURATE'; %group_ephys = '20230801_ChocolateGroup';
end
group_setup_behav = strcat('20230511_ChocolateGroup',filesep,'headfixed_dynamicTarget');

% IDs and Definitions
animals = {...
    'CoteDor',...
    'Lindt',...
    'Toblerone',...
    'Milka',...
    'FerreroRocher'};
animal_idx = 5;
mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});
%ephys_root = strcat(ephys_root,filesep,mouse);

paw_pref = 'right';
session = 'R6';
ephys_sess = '20082023_Ferrero_StrCer_S6_g0';
imec_id = 1;

%catGT_folder = 'catGT_KS_DSRemoved';
sorter_folder = 'catGT\kilosort4';
%sorter_folder = 'ibl_sorter_results_driftAdapt';
output_folder_name = 'neurons_overview_post_postPredictions2';

% To run / save
show_aux_plots = 1;
save_mat_flag = 1;
plot_neuron_fig = 1;
use_inferred_reach_times = 0;

% Probe side: BG or CB
if strcmp(paw_pref,'right')
    if imec_id == 0, region = 'CB';
    elseif imec_id == 1, region = 'BG';
    end
elseif  strcmp(paw_pref,'left')
    if imec_id == 0, region = 'BG';
    elseif imec_id == 1, region = 'CB';
    end
end

%% Path to load & save ------------------------------------------------
behavior_path = strcat(behavior_root,filesep,group_setup_behav,filesep,mouse,filesep,session);
ephys_path = strcat(ephys_root,filesep,ephys_sess,filesep,ephys_sess,'_imec',num2str(imec_id));
ephys_LFP_path = ephys_path;
myKsDir = strcat(ephys_path,filesep,sorter_folder);
% path to save fig
out_ephys_folder = strcat(myKsDir,filesep,output_folder_name);
if ~exist(out_ephys_folder,'dir'), mkdir(out_ephys_folder); end
% load behavior data
load(strcat(behavior_path,filesep,'behavior_session.mat'));

% load kilosort output with spikes and metadata
tic
fprintf('Loading data...\n')
sp = loadKSdir_TD(myKsDir);
meta_ap = readSpikeGLXmeta(strcat(ephys_path,filesep,dir(fullfile(ephys_path,'*.ap.meta')).name));
meta_lf = readSpikeGLXmeta(strcat(ephys_path,filesep,dir(fullfile(ephys_path,'*.lf.meta')).name));
toc

% save all session spikes and metadata
neurons_params.sp = sp;
neurons_params.sp = meta_ap;
neurons_params.sp = meta_lf;

%% Sync
nChansInFile = 385;  % neuropixels phase3a, from spikeGLX
syncChanIndex = 385;
tic
syncDat = extractSyncChannel(ephys_LFP_path, nChansInFile, syncChanIndex);

lfpFs = meta_lf.imSampRate;
eventTimes = spikeGLXdigitalParse(syncDat, lfpFs);
eventTimes_water_giv = eventTimes{7}{1};

eventTime_ind = find(diff([0,syncDat])~=0);
eventTime_ind_up = find(diff([0,syncDat])>1);
eventTime_ind_down = find(diff([0,syncDat])<-1);

eventTimes_realTrial_ephys = eventTimes_water_giv(2:end);
eventTimes_realTrial_harp = behavior.sync.time_newTrial_log(2:end);

if numel(eventTimes_realTrial_ephys)<numel(eventTimes_realTrial_harp)
    if round(sum(diff(eventTimes_realTrial_harp(2:end))-diff(eventTimes_realTrial_ephys))) == 0
        eventTimes_realTrial_harp = behavior.sync.time_newTrial_log(3:end);
        flag_remove_1st_behavior_frame = 1;
        behavior.sync.flag_remove_1st_behavior_frame = flag_remove_1st_behavior_frame;
        save(strcat(behavior_path,filesep,'behavior_session.mat'),'behavior');
    end
    fprintf('Removed the first trial from behavior! \n')
end

f=fit(eventTimes_realTrial_ephys,eventTimes_realTrial_harp,'poly1');
tm_bhv2ephys =  @(eventTimes_realTrial_ephys)  f.p2+f.p1*eventTimes_realTrial_ephys;
toc

if show_aux_plots
    figure
    plot(syncDat)
    hold on
    plot(eventTime_ind,ones(size(eventTime_ind)).*64,'*')
    hold off

    figure
    subplot(411), plot(diff(eventTimes_water_giv));
    subplot(412), plot(diff(eventTime_ind));
    subplot(413), plot(diff(eventTimes_realTrial_harp));
    subplot(414), plot(diff(behavior.sync.time_newTrial_log));

    if  numel(eventTimes_realTrial_harp)==numel(behavior.sync.time_newTrial_log(2:end))
        figure
        plot(eventTimes_water_giv,behavior.sync.time_newTrial_log,'-'); hold on
        plot(eventTimes_water_giv,behavior.sync.time_newTrial_log,'.');  hold off
        ylabel('trial time on harp'); xlabel('trial time on ephys');
    end

    figure
    plot(eventTimes_realTrial_ephys,eventTimes_realTrial_harp,'-'); hold on
    plot(eventTimes_realTrial_ephys,eventTimes_realTrial_harp,'.');  hold off
    ylabel('trial time on harp'); xlabel('trial time on ephys');
    axis square


    figure, plot(f,eventTimes_realTrial_ephys,eventTimes_realTrial_harp); shg
    ylabel('trial time on harp'); xlabel('trial time on ephys');
    title(sprintf('%s%.5f','slop = ',f.p1));
end

%%

neurons=[];
% Match the sample rate
Fs = sp.sample_rate;

% Find the unique units and filter according to quality
QualityFilter = [1,2]; %Quality Filter (1 = MUA, 2 = Good)
uniqueUnits = sp.cids;
ExtraGood = sp.cegs;
filterQualityUnits = ismember(sp.cgs, QualityFilter);
uniqueUnits = uniqueUnits(filterQualityUnits);
filteredQualities = sp.cgs(filterQualityUnits);

% Cell type: is complex spike?
if isfield(sp,'potential_CS')
    pCS = sp.potential_CS;
end

% Find total number of units
nOfUnits = numel(uniqueUnits);

% Find the duration of the behavior session with engagement
behav_start = behavior.behavior_duration.time_start;
behav_stop = behavior.behavior_duration.time_end + behav_start;
behavior_end_timelog = behavior.behavior_duration.time_end;
all_st = tm_bhv2ephys(sp.st);
st_behaviorEpochBool = all_st>behav_start & all_st<behav_stop;



%% Populate structure with unit's data
temp_neurons = struct();

fprintf('Generating neurons struct...\n')
tic
for iUnit = 1:nOfUnits
    %iUnit = 1;

    % find the spikes for the cluster and time
    cluSpikeBool = sp.clu == uniqueUnits(iUnit);
    thisSpikeBool = st_behaviorEpochBool & cluSpikeBool;


    % spike times (raw and in behavior time)
    temp_neurons(iUnit).st_raw_all = sp.st(cluSpikeBool);
    temp_neurons(iUnit).st_raw = sp.st(thisSpikeBool);
    temp_neurons(iUnit).st = all_st(thisSpikeBool);

    %unit ID used in Phy
    temp_neurons(iUnit).phyID = uniqueUnits(iUnit);

    %unit ID in the session
    temp_neurons(iUnit).unitIdx = iUnit;

    % templates waveforms
    templates_clu = unique(sp.spikeTemplates(thisSpikeBool));
    weighth_template = nan(1,length(templates_clu));
    neuron_template_waveform = nan(size(sp.templateWaveforms,2:3));
    temp_id = nan(1,length(templates_clu));
    neuron_template_PeakChannel = nan(1,length(templates_clu));
    for j = 1:length(templates_clu)
        clu_template = sp.spikeTemplates == templates_clu(j);
        clu_template_spikes = thisSpikeBool & clu_template;
        weighth_template(1,j) = sum(clu_template_spikes)./sum(thisSpikeBool);
        temp_id(1,j)=find(sp.template_phyIDs==templates_clu(j));
        neuron_template_waveform(:,:,j)  = squeeze(sp.templateWaveforms(temp_id(j),:,:));
        neuron_template_PeakChannel(1,j) = sp.templatePeakChn(temp_id(j));
    end
    %weighth_template(weighth_template>1)=1; % THIS IS AN EASY FIX, SOLVE IT
    temp_neurons(iUnit).templateID = temp_id;
    temp_neurons(iUnit).templateWeight = weighth_template;
    temp_neurons(iUnit).templateWaveforms = neuron_template_waveform;
    temp_neurons(iUnit).templatePeakCh =neuron_template_PeakChannel;

    %unit amplitudes
    temp_neurons(iUnit).amplitudes = sp.tempScalingAmps(thisSpikeBool);

    %unit quality
    temp_neurons(iUnit).quality = filteredQualities(iUnit);
    temp_neurons(iUnit).extraGood = ExtraGood(iUnit);

    % if relevant, is this unit a CS
    if ~isempty(pCS)
        temp_neurons(iUnit).potentialCS = pCS(iUnit);
    end

    %Depth from kilo
    temp_neurons(iUnit).depth =nanmedian(sp.spikeDepths(thisSpikeBool));

    %recorded animal
    temp_neurons(iUnit).meta.recordedAnimal = mouse;

    %recorded session
    temp_neurons(iUnit).meta.recordedSessionName = ephys_sess;
    temp_neurons(iUnit).meta.recordedSession = session;
    temp_neurons(iUnit).meta.imecID = imec_id;
    temp_neurons(iUnit).meta.region = region;
    temp_neurons(iUnit).meta.curator = person;
    temp_neurons(iUnit).meta.use_inferred_reach_times = use_inferred_reach_times;

end
neurons = [neurons, temp_neurons];
toc

%% ACTIVITY ALIGN to INIT & REACH + OVERALL SESSION

% trials
nr_trials = behavior.behavior_duration.trial_end;
trials_vec = 1:nr_trials;

% events
initiation_times = behavior.inputs.read_log(behavior.logs.trial_init_ind(trials_vec),1);
% reach time already corrected for180 offset related to water collection detection
if use_inferred_reach_times == 1
    %reach_time_log = behavior.reach.timeof.reach_time_inferred(trials_vec);
    reach_time_log = behavior.reach.timeof.reach_time_inferred_approxlog_corrected(trials_vec);
    reach_times = reach_time_log +...
        behavior.behavior_duration.time_start +...
        behavior.inputs.timelog(1); % because times are related to session ready in the rwd collect correction
else    
    reach_times = behavior.reach.timeof.reach_time_corrected(trials_vec) + behavior.behavior_duration.time_start;
end

% idx
push_idx = behavior.init.idx_trial_push(behavior.init.idx_trial_push<=nr_trials);
pull_idx = behavior.init.idx_trial_pull(behavior.init.idx_trial_pull<=nr_trials);
left_idx = behavior.reach.left_idx(behavior.reach.left_idx<=nr_trials);
right_idx = behavior.reach.right_idx(behavior.reach.right_idx<=nr_trials);
center_idx = behavior.reach.center_idx(behavior.reach.center_idx<=nr_trials);

% gamma for PSTH
peak_x = .05;
bin_width = .002;
k = gammakernel('peakx',peak_x,'binwidth',bin_width);

% parameters around events
tm_before = -3; % sec
tm_after = 2; % sec
bin_edges = tm_before+k.paddx(1):bin_width:tm_after+k.paddx(2);
spk_bins = tm_before:bin_width:tm_after;

% params for mean session FR
nbins_FR = 50;

% Params for cross-correlogram
binsize_CC = 0.001;
duration_CC  = 0.3;

% Params for isi
isi_bins = [0:0.0005:0.1];
isi_bins_cs = [0:0.02:5];

neurons_params.events.tm_before = tm_before;
neurons_params.events.tm_after = tm_after;
neurons_params.events.bin_width = bin_width;
neurons_params.events.bin_edges = bin_edges;
neurons_params.gamma.k = k;
neurons_params.gamma.spk_bins = spk_bins;
neurons_params.crosscorrelogram.binsize_CC = binsize_CC;
neurons_params.crosscorrelogram.duration_CC = duration_CC;
neurons_params.isi.isi_bins = isi_bins;
neurons_params.isi.isi_bins_cs = isi_bins_cs;

%% FIGURE PROPERTIES
% axis and colors
figProp.axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0.05,'ticklength',[1,1]*.01,'fontsize',10};

%clrss_in = get(gca,'ColorOrder');
push_clr = behavior.colors.push_clr;
pull_clr =  behavior.colors.pull_clr;
pp_clr = [63 130 109]./256;
right_color = behavior.colors.right_color;
center_color = behavior.colors.center_color;
left_color = behavior.colors.left_color;
rcl_clr = [44 123 182]./256;
lw=2;

spk_size = 2;
event_size = 5;
patch_wd = 10;

% Line plot
line_x = [tm_before;tm_before]+0.08;
% init
pull_bounds = behavior.init.pull_bounds;
push_bounds = behavior.init.push_bounds;
%reach
left_bounds =  behavior.reach.left_bounds;
right_bounds = behavior.reach.right_bounds;
center_bounds = behavior.reach.center_bounds;

axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10};
figProp.axeOpt =axeOpt;
figProp.lw = 2;
figProp.transitionBarClr = [.8 .8 .8 .5];
figProp.lw_ppT = 3;
figProp.lw_lcrT = 1;
figProp.pa_pp_alpha = .8;
figProp.pa_lcr_alpha = .8;
figProp.spk_clr =[.2 .2 .2 .1];
figProp.spk_clr_hist =[.2 .2 .2];
figProp.FR_clr =[.8 .8 .9];
figProp.FR_clrFace =[.8 .8 .9];
figProp.valid_clr = [0 0 0];
figProp.invalid_clr = [.7 .7 .7];
figProp.valid_clr_fc = [.8 .8 .8];
figProp.invalid_clr_fc = [.8 .8 .8];


%% Run through units
fprintf('Align activity to behavior...\n')
tic
for un = 1: length(neurons)
    %un = 90;
    if ~isempty(neurons(un).st_raw)

        %un=1;
        % ------------------------------------------------------
        % ACTIVITY ALIGN TO INIT / REACH -----------------------
        % ------------------------------------------------------
        %tic
        neurons(un).st_init = cell(nr_trials,1);
        neurons(un).spk_trials_init = cell(nr_trials,1);
        neurons(un).st_reach = cell(nr_trials,1);
        neurons(un).spk_trials_reach = cell(nr_trials,1);
        neurons(un).reach_in_init = nan(1,nr_trials);
        neurons(un).init_in_reach = nan(1,nr_trials);
        for tt = 1:nr_trials
            % init: push/pull
            spike_init_flags = ...
                neurons(un).st > initiation_times(tt) + bin_edges(1) & ...
                neurons(un).st <= initiation_times(tt) + bin_edges(end);
            neurons(un).st_init{tt} = ...
                neurons(un).st(spike_init_flags) - initiation_times(tt);
            neurons(un).reach_in_init(tt) = reach_times(tt) - initiation_times(tt);
            neurons(un).spk_trials_init{tt} = ones(size(neurons(un).st_init{tt}))*tt;
            % reach_ left/center/right
            spike_reach_flags = ...
                neurons(un).st > reach_times(tt) + bin_edges(1) & ...
                neurons(un).st <= reach_times(tt) + bin_edges(end);
            neurons(un).st_reach{tt} = ...
                neurons(un).st(spike_reach_flags) - reach_times(tt);
            neurons(un).init_in_reach(tt) = initiation_times(tt) - reach_times(tt);
            neurons(un).spk_trials_reach{tt} = ones(size(neurons(un).st_reach{tt}))*tt;
        end

        % SPK COUNTS ----------------------------------------------
        % Divide in 2 ms bins and count spikes
        bin_counts_init = nan(length(bin_edges)-1,nr_trials);
        bin_counts_reach = nan(length(bin_edges)-1,nr_trials);
        for tt = 1:nr_trials
            bin_counts_init(:,tt) = histcounts(neurons(un).st_init{tt},bin_edges);
            bin_counts_reach(:,tt) = histcounts(neurons(un).st_reach{tt},bin_edges);
        end
        % Spike rates: convolve with gamma function
        neurons(un).spk_rates_init = conv2(k.pdf,1,bin_counts_init./bin_width,"valid");
        neurons(un).spk_rates_reach = conv2(k.pdf,1,bin_counts_reach./bin_width,"valid");

        %  PSTHs mean and sem -------------------------------
        % Push and pull PSTH
        psth_push = mean(neurons(un).spk_rates_init(:,push_idx),2);
        psth_pull = mean(neurons(un).spk_rates_init(:,pull_idx),2);
        if length(push_idx)<2, psth_push=nan(size(psth_push)); end
        if length(pull_idx)<2, psth_pull=nan(psth_pull); end
        psth_push_sem = std(neurons(un).spk_rates_init(:,push_idx),[],2)./length(push_idx);
        psth_pull_sem = std(neurons(un).spk_rates_init(:,pull_idx),[],2)./length(pull_idx);

        % Right / Left / Center PSTH
        psth_right = mean(neurons(un).spk_rates_reach(:,right_idx),2);
        psth_left = mean(neurons(un).spk_rates_reach(:,left_idx),2);
        psth_center = mean(neurons(un).spk_rates_reach(:,center_idx),2);
        psth_right_sem = std(neurons(un).spk_rates_reach(:,right_idx),[],2)./length(right_idx);
        psth_left_sem = std(neurons(un).spk_rates_reach(:,left_idx),[],2)./length(left_idx);
        psth_center_sem = std(neurons(un).spk_rates_reach(:,center_idx),[],2)./length(center_idx);
        %toc

        % -------------------------------------------------
        % SESSION OVERVIEW
        % -------------------------------------------------
        % FR all session ----------------------------------
        [FR_sess,neurons(un).edges_sess] = histcounts((neurons(un).st-behav_start),nbins_FR);
        bin_size_sec =diff(neurons(un).edges_sess(1:2));
        neurons(un).FR_sess_mean = FR_sess./bin_size_sec;

        % Cross correlogram -------------------------------
        [ccg, ccg_t] = CCGBz([double(neurons(un).st); double(neurons(un).st)], [ones(size(neurons(un).st, 1), 1); ...
            ones(size(neurons(un).st, 1), 1) * 2], 'binSize', binsize_CC, 'duration', duration_CC, 'norm', 'rate'); %function
        neurons(un).CCG = ccg(:,1,1);
        neurons(un).CCG_bins = ccg_t;

        % ISI  --------------------------------------------
        ISI = diff(neurons(un).st);
        [neurons(un).isiProba, neurons(un).edgesISI] = histcounts(ISI,isi_bins);
        if sum(pCS) ~= 0
            if neurons(un).potentialCS == 1
                [neurons(un).isiProba, neurons(un).edgesISI] = histcounts(ISI,isi_bins_cs);
            end
        end


        % BEHAVIOR ALIGNED --------------------------------
        inval_trials_tmp = sort([behavior.init.timeof.invalPull_time;behavior.init.timeof.invalPush_time]);
        inval_trials = inval_trials_tmp(inval_trials_tmp>0 & inval_trials_tmp<behavior_end_timelog);
        [val_init,ed] = histcounts(behavior.init.timeof.init_time(trials_vec),neurons(un).edges_sess);
        [inval_init,ed2] = histcounts(inval_trials,neurons(un).edges_sess);


        %% FIGURE
        if plot_neuron_fig == 1
            if neurons(un).quality==2
                %define layout
                fig=figure();
                tt = tiledlayout(5,7);
                title(tt,sprintf('%s%s%s%s%s%s%s%s%s%s%s%s',...
                    'mouse: ',neurons(un).meta.recordedAnimal,...
                    ' | session: ',neurons(un).meta.recordedSession,...
                    ' | probe: imec', num2str(neurons(un).meta.imecID),...
                    ' | region: ',neurons(un).meta.region,...
                    ' | neuron phyID: ', num2str(neurons(un).phyID),...
                    ' | neuron idx: ', num2str(neurons(un).unitIdx)),...
                    'Interpreter', 'none','fontsize',12,'fontweight','normal')

                % --------------------------------------------------------------
                % NEURON DEPTH
                axx0 = nexttile;
                axx0.Layout.TileSpan = [3,1];

                [depth_ocupacy,depth]=histcounts(sp.spikeDepths,50);
                histogram('BinEdges', depth ,'BinCounts',  depth_ocupacy,...
                    'facecolor',figProp.FR_clrFace ,'facealpha',0.6,'edgecolor','none');
                set(gca,figProp.axeOpt{:},'XDir','reverse')
                camroll(-90)
                h = gca;
                xlim([0 3500])
                h.YAxis.Visible = 'off';
                xline(neurons(un).depth,'-',sprintf('%s%i','Unit ID: ',neurons(un).phyID)...
                    ,'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle',...
                    'color','k','linewidth',3)
                xlabel('Depth of neuropixel probe (\mum)');

                % --------------------------------------------------------------
                % PSTH ALIGN TO INIT
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

                % --------------------------------------------------------------
                % PSTH ALIGN TO REACH
                axx2 = nexttile;
                axx2.Layout.TileSpan = [1,2];

                shadedErrorBar(spk_bins,psth_right,psth_right_sem,'lineProps',{'Color',right_color,'LineWidth',2},'transparent',1,'patchSaturation',0.1); hold on
                shadedErrorBar(spk_bins,psth_left,psth_left_sem,'lineProps',{'Color',left_color,'LineWidth',2},'transparent',1,'patchSaturation',0.1); hold on
                shadedErrorBar(spk_bins,psth_center,psth_center_sem,'lineProps',{'Color',center_color,'LineWidth',2},'transparent',1,'patchSaturation',0.1); hold on
                ylabel('spike rate (sp/s)'); xlabel('time from reach (s)');
                xlim([tm_before tm_after])
                xline(0,'--','Color',[rcl_clr 0.3],'linewidth',2);
                set(gca,axeOpt{:})
                title('water collection (reach)')
                %title(sprintf('%s%s','Cell unit: ',num2str(neurons(un).phyID)));

                % --------------------------------------------------------------
                % CROSS-CORRELOGRAM
                axxx1 = nexttile;
                axxx1.Layout.TileSpan = [1,2];
                bar(ccg_t.*1000, ccg(:,1,1),'BarWidth', 1,...
                    'facecolor',figProp.spk_clr_hist,'edgecolor','none','facealpha',.9)
                set(gca,figProp.axeOpt{:})
                xlabel('time (ms)'); ylabel('spike rate (sp/s)')
                title('cross-correlogram')

                % --------------------------------------------------------------
                % SCATTER SPIKE TIMES ALIGNED TO INIT
                axx3 = nexttile;
                axx3.Layout.TileSpan = [2,2];

                scatter(cell2mat(neurons(un).st_init),cell2mat(neurons(un).spk_trials_init),spk_size,'k.'), hold on
                scatter(neurons(un).reach_in_init,trials_vec,event_size,'filled','MarkerFaceColor',rcl_clr)
                scatter(zeros(size(trials_vec)),trials_vec,event_size,'filled','MarkerFaceColor',pp_clr)
                plot(line_x,pull_bounds,'linewidth',patch_wd,'color',pull_clr)
                plot(line_x,push_bounds,'linewidth',patch_wd,'color',push_clr)
                hold off
                axis tight;
                xlabel('time from trial init (s)'); ylabel('trials in session')
                axis([tm_before tm_after 0 nr_trials])
                set(gca,axeOpt{:})
                set(gcf,'Position',[2582,215,560,771],'color','w')

                % --------------------------------------------------------------
                % SCATTER SPIKE TIMES ALIGNED TO REACH
                axx4 = nexttile;
                axx4.Layout.TileSpan = [2,2];

                scatter(cell2mat(neurons(un).st_reach),cell2mat(neurons(un).spk_trials_reach),spk_size,'k.'), hold on
                scatter(neurons(un).init_in_reach,trials_vec,event_size,'filled','MarkerFaceColor',pp_clr)
                scatter(zeros(size(trials_vec)),trials_vec,event_size,'filled','MarkerFaceColor',rcl_clr)
                plot(repmat(line_x,[1 size(right_bounds,2)]),right_bounds,'linewidth',patch_wd,'color',right_color)
                plot(repmat(line_x,[1 size(left_bounds,2)]),left_bounds,'linewidth',patch_wd,'color',left_color)
                plot(repmat(line_x,[1 size(center_bounds,2)]),center_bounds,'linewidth',patch_wd,'color',center_color)
                hold off
                axis tight;
                axis([tm_before tm_after 0 nr_trials])
                xlabel('time from water collection (s)'); ylabel('trials in session')
                set(gca,axeOpt{:})

                % --------------------------------------------------------------
                % ISI
                axxx2 = nexttile;
                axxx2.Layout.TileSpan = [2,2];

                histogram('BinEdges', neurons(un).edgesISI.*1000 ,'BinCounts',  neurons(un).isiProba,...
                    'facecolor',figProp.spk_clr_hist,'edgecolor','none','facealpha',.9)
                set(gca,axeOpt{:})
                xlabel('interspike interval (ms)'); ylabel ('# of spikes')
                %title('interspike interval')

                if neurons(un).quality==2
                    cell_qual = 'good';
                elseif neurons(un).quality==1
                    cell_qual = 'mua';
                else
                    cell_qual = 'noise';
                end
                dim_inf=[0.84 0.31 0.4 0.4];
                annotation('textbox',dim_inf,'String',{...
                    "quality = " + cell_qual,...
                    "extra good = " + num2str(neurons(un).extraGood),...
                    "# spikes = " + num2str(length(neurons(un).st)),...
                    },'edgecolor','w','fontsize',10.5, 'FitBoxToText','on')
                %
                if sum(pCS) ~= 0
                    if neurons(un).potentialCS == 1
                        annotation('textbox',[0.84 0.25 0.4 0.4],'String',...
                            "potential CS!",...
                            'edgecolor','w','fontsize',10.5, 'FitBoxToText','on')
                    end
                end

                % --------------------------------------------------------------
                % BEHAVIOR THROUGHOUT SESSION
                ax1 = nexttile;

                fr_plot = 0;
                maxBev =  max([val_init,inval_init])+1;
                maxFR = ceil(max(neurons(un).FR_sess_mean))+1;
                max_yLim = max([maxBev,maxFR])+2;
                plen = 2;
                patchY.extend = 1;
                patchY.showPatch = 1;

                % FR
                if fr_plot
                    yyaxis left
                    histogram('BinEdges', neurons(un).edges_sess./60 ,'BinCounts',  neurons(un).FR_sess_mean,...
                        'displaystyle','stairs','linewidth',1.5, 'edgecolor',figProp.FR_clr); hold on
                    histogram('BinEdges', neurons(un).edges_sess./60 ,'BinCounts',  neurons(un).FR_sess_mean,...
                        'facecolor',figProp.FR_clrFace ,'facealpha',0.3,'edgecolor','none')
                    ylim([0 maxFR+plen*2]);
                    set(gca, 'YColor', figProp.FR_clr);
                    if (maxFR>maxBev)
                        patchY.pp = [maxFR+plen maxFR+(plen*2) maxFR+(plen*2) maxFR+plen];
                        patchY.rlc = [maxFR maxFR+plen maxFR+plen maxFR];
                        landmarks_plt = behavior_landmarks(behavior, figProp,patchY);
                    end
                    ylabel('mean firing rate (sp/s)','color', figProp.FR_clr)

                    % behavior
                    yyaxis right
                end
                p1 = histogram('BinEdges', ed./60 ,'BinCounts',  val_init,...
                    'displaystyle','stairs','linewidth',1.5, 'edgecolor',figProp.valid_clr,'LineStyle','-');
                hold on
                histogram('BinEdges', ed./60 ,'BinCounts',  val_init,...
                    'facecolor',figProp.valid_clr_fc ,'facealpha',0.3,'edgecolor','none')
                p2 = histogram('BinEdges', ed2./60 ,'BinCounts',  inval_init,...
                    'displaystyle','stairs','linewidth',1.5, 'edgecolor',figProp.invalid_clr,'LineStyle','--');
                histogram('BinEdges', ed2./60 ,'BinCounts',  inval_init,...
                    'facecolor',figProp.invalid_clr_fc ,'facealpha',0.3,'edgecolor','none')
                ylim([0 maxBev+1]);
                xlim([0 behavior_end_timelog/60])
                ylabel('trial counts (3 min^{-1})','Color', figProp.spk_clr)
                set(gca, 'YColor', figProp.valid_clr);
                xlabel('time across session (min)');
                ylim([0 maxBev+plen*2]);
                if (maxFR<maxBev  || fr_plot==0)
                    patchY.pp = [maxBev+plen maxBev+(plen*2) maxBev+(plen*2) maxBev+plen];
                    patchY.rlc = [maxBev maxBev+plen maxBev+plen maxBev];
                    landmarks_plt = behavior_landmarks(behavior, figProp,patchY);
                end
                %     xlim([0 behavior.session_duration/60])
                %     h2 = gca;
                h2.XAxis.Visible = 'off';
                title('Behavior and cell activity across the session')
                lgd1 = legend([p1, p2], 'valid counts', 'invalid counts',...
                    'location','east','edgeColor', 'w');

                % legend
                axleg = nexttile(26);
                lgd2=legend(axleg,landmarks_plt,...
                    'push','pull','left','right','center',...
                    'edgeColor','w','fontsize',10.5);
                lgd2.Layout.Tile = 26;
                set(axleg,'Visible', 'off');
                ax1.Layout.TileSpan = [1, 4];



                % --------------------------------------------------------------
                % WAVWFORMS
                axxx3 = nexttile;
                axxx3.Layout.TileSpan = [2,2];

                time_wf = (0:1/sp.sample_rate:(size(neurons(un).templateWaveforms,1)-1)/sp.sample_rate)*1000;
                plot(time_wf,zero2nan(neurons(un).templateWaveforms(:,:)),'linewidth',1.5,'Color',[.8 .8 .8]);
                hold on
                for i = 1:length(neurons(un).templateID)
                    plot(time_wf,zero2nan(neurons(un).templateWaveforms(:,neurons(un).templatePeakCh(i))),...
                        'linewidth',2,'Color',[.1 .1 .1 neurons(un).templateWeight(i)]);
                end
                hold off
                set(gca,figProp.axeOpt{:})
                xlabel('time (ms)'); ylabel('template waveform (pseudo-volts)')
                str = 'Straight Line Plot from 1 to 10';
                dim_ch=[0.85 0.12 0.1 0.1];
                title('template waveform')
                annotation('textbox',dim_ch,'String',...
                    sprintf('%s%s','channel = ',num2str(neurons(un).templatePeakCh)),...
                    'edgecolor','w','fontsize',10.5)

                % --------------------------------------------------------------
                % AMPLITUDE & FIRING RATE
                ax2 = nexttile;
                ax2.Layout.TileSpan = [1,4];
                linkaxes([ax1, ax2], 'x');
                lp = 0; % if we want to display the patches

                % amplitude
                AmpLimMax = ceil(max(neurons(un).amplitudes));
                AmpLimMin = floor(min(neurons(un).amplitudes))-5;
                %colororder({figProp.FR_clr,figProp.spk_clr(1:3)})
                yyaxis left
                amp_plt = plot((neurons(un).st-behav_start)./60,neurons(un).amplitudes,'.','Color', figProp.spk_clr);
                hold on
                ylim([0 AmpLimMax+lp]);
                ylabel('amplitude (pseudo-volts)','Color', figProp.spk_clr)
                set(gca,figProp.axeOpt{:})
                set(gca,'YColor',figProp.spk_clr)

                % firing rate
                yyaxis right
                histogram('BinEdges', neurons(un).edges_sess./60 ,'BinCounts',  neurons(un).FR_sess_mean,...
                    'displaystyle','stairs','linewidth',1.5, 'edgecolor',figProp.FR_clr)
                histogram('BinEdges', neurons(un).edges_sess./60 ,'BinCounts',  neurons(un).FR_sess_mean,...
                    'facecolor',figProp.FR_clrFace ,'facealpha',0.6,'edgecolor','none')
                maxFR = ceil(max(neurons(un).FR_sess_mean));
                ylim([0 maxFR+(lp*2)]);
                set(gca, 'YColor', figProp.spk_clr);
                ylabel('mean firing rate (sp/s)','color', figProp.FR_clr)

                maxY = max([maxFR AmpLimMax]);
                patchY.extend = 1;
                patchY.showPatch = 0;
                landmarks_plt2 = behavior_landmarks(behavior, figProp,patchY);
                xlim([0 behavior_end_timelog./60])
                xlabel('time across session (min)');
                set(gca, 'YColor', figProp.FR_clr);


                % --------------------------------------------------------------
                % AMPLITUDE HISTOGRAM
                ax3 =nexttile;
                nbins = 30;
                histogram(neurons(un).amplitudes,nbins,'FaceColor',...
                    figProp.spk_clr_hist,'facealpha',.5);
                xlim([0 AmpLimMax+6])
                hax = gca;
                hax.YDir = 'reverse';
                camroll(90)
                set(gca,figProp.axeOpt{:})
                ylabel('density'); xlabel('amplitude');


                set(gcf,'Position',[2078 49 1583 949],'color','w')
                %drawnow;
                % Save figure
                saveas(gcf,strcat(out_ephys_folder,filesep,'neuron',num2str(neurons(un).phyID),'.png'),'png');
            end
        end
    end
end
toc

%% SAVE MAT
if save_mat_flag==1
    fprintf('Saving data...\n')
    save(strcat(myKsDir,filesep,'neurons_session'),'neurons',...
        'neurons_params','figProp','-v7.3');
       %save(strcat(myKsDir,filesep,'neuronStruct_allVar'),'-v7.3');
end
fprintf('Done!!')



