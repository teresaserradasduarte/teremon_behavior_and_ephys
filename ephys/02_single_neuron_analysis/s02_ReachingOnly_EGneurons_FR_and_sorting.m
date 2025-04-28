% Get extra-good neurons, calculate instantaneous FR, order neurons
clear; close all; clc

%% Load data
% Group and individual
group = '20230801_ChocolateGroup';
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
imec_id = 1;
ephys_sess = '18082023_Milka_StrCer_S4_g0';
sorter_folder = 'catGT\kilosort4';
ephys_local_folder = 1;

%% Depth cortex 
depth_cortex_lim = nan;
%depth_cortex_lim = 2000;
depthLim_phyID = 363;

%% paths
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
if ephys_local_folder==1, ephys_root = 'E:\';
else, ephys_root = fullfile(rootdir,"ephys_data",group,mouse); end
behavior_dir = fullfile(rootdir,"behavior_data","raw_data",group,setup,mouse,sess);
reaching_dir = fullfile(rootdir,"behavior_data","analyzed_data","mat_files",group,setup,mouse,sess);
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
save_mat = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,region);
save_out = fullfile(rootdir,"ephys_and_behavior","out_files",group,mouse,region);
if ~exist(save_mat,"dir"), mkdir(save_mat); end
if ~exist(save_out,"dir"), mkdir(save_out); end

% Figure
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};


%% Create neuron struct for extra-good neurons only
neu_all = struct2table(neurons);
eg_idx = find(neu_all.extraGood == 1);
n_eg = length(eg_idx);

if strcmp(region,"BG")
    if isnan(depth_cortex_lim)
        limit_depth_neuron = find(neu_all.phyID == depthLim_phyID);
        depth_lim = neurons(limit_depth_neuron).depth;
    else
        depth_lim = depth_cortex_lim;
    end
end

eg_neurons = struct();
for i=1:n_eg
    eg_neurons(i).st = neurons(eg_idx(i)).st_all;
    eg_neurons(i).phyID = neurons(eg_idx(i)).phyID;
    eg_neurons(i).unitIdx = neurons(eg_idx(i)).unitIdx;
    eg_neurons(i).depth = neurons(eg_idx(i)).depth;
    eg_neurons(i).egIdx = i;
    if strcmp(region,"BG")
        if (eg_neurons(i).depth <= depth_lim)
            eg_neurons(i).reg = 'BG ';
            eg_neurons(i).isBG = 1;
        elseif (eg_neurons(i).depth > depth_lim)
            eg_neurons(i).reg = 'CTX';
            eg_neurons(i).isBG = 0;
        end
    elseif strcmp(region,"CB")
        eg_neurons(i).reg = 'CB';
    end
end
if strcmp(region,"BG")
    idx_BG = find(struct2table(eg_neurons).depth <= depth_lim);
    idx_CT = find(struct2table(eg_neurons).depth > depth_lim);
elseif strcmp(region,"CB")
    idx_CB = 1:n_eg;
end

%% Save
save(fullfile(save_mat,'eg_neurons.mat'),'eg_neurons','idx_BG','idx_CT','BG_reach_order','-v7.3');
clear neu_all neurons;

%% Get event time from behavior
% Win of good behavior + ehpys
behav_start = behavior.behavior_duration.time_start;
behav_stop = behavior.behavior_duration.time_end + behav_start;
nr_trials = behavior.behavior_duration.trial_end;
trials_vec = 1:nr_trials;

% Trial identity defined by logs
push_idx = behavior.init.idx_trial_push;
pull_idx = behavior.init.idx_trial_pull;
left_idx = behavior.reach.left_idx;
right_idx = behavior.reach.right_idx;
center_idx = behavior.reach.center_idx;

% Event times from logs - trial init
initiation_times = behavior.inputs.read_log(behavior.logs.trial_init_ind,1);
inva_push_times = behavior.inputs.read_log(behavior.logs.inval_push_ind,1);
inva_pull_times = behavior.inputs.read_log(behavior.logs.inval_pull_ind,1);

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

figProp.spk_size = 6;
figProp.rcl_clr = [44 123 182]./256;

%% Align to reaches
event_time = reach_times_inVec;
FR_reach = zeros(nr_bins,n_reach_inVec,n_eg);

for i = 1:n_eg
    FR_reach_raw = zeros(n_bins_padded,n_reach_inVec);
    for j = 1:n_reach_inVec
        spike_reach_flags = ...
            eg_neurons(i).st > event_time(j) + win_padded(1) & ...
            eg_neurons(i).st <= event_time(j) + win_padded(2);

        eg_neurons(i).st_reach{j} = ...
            eg_neurons(i).st(spike_reach_flags) - event_time(j);
        eg_neurons(i).st_reachIdx{j} = ones(size(eg_neurons(i).st_reach{j}))*j;
        reach_spikes = eg_neurons(i).st_reach{j};

        if ~isempty(reach_spikes)
            % Find the spike before interval fo initial ISI
            first_spk_idx=find(spike_reach_flags,1,'first');
            if first_spk_idx == 1
                FR_reach_raw(bin_edges_padded<=reach_spikes(1),j) = 0;
            else
                start_instant_isi = eg_neurons(i).st(first_spk_idx) - ...
                    eg_neurons(i).st(first_spk_idx-1);
                FR_reach_raw(bin_edges_padded<=reach_spikes(1),j) = 1/start_instant_isi;
            end
            % Find spikes within the interval
            for sp=2:length(reach_spikes)
                with_isi_flag = (bin_edges_padded > reach_spikes(sp-1)) & ...
                    bin_edges_padded <= reach_spikes(sp);
                inst_isi = reach_spikes(sp)-reach_spikes(sp-1);
                FR_reach_raw(with_isi_flag,j) = 1/inst_isi;
            end

            % find the spike after the interval for final ISI
            last_spk_idx=find(spike_reach_flags,1,'last');
            finish_instant_isi = eg_neurons(i).st(last_spk_idx+1) - ...
                eg_neurons(i).st(last_spk_idx);
            FR_reach_raw(bin_edges_padded>reach_spikes(end),j) = 1/finish_instant_isi;
        end
    end
    % smoothing individual trials & store
    FR_reach(:,:,i) = conv2(k_gaus.pdf,1,FR_reach_raw,"valid");
    eg_neurons(i).FR_reach = FR_reach(:,:,i);

    % Plot
    figProp.name = 'All reaches';
    plot_raster_and_psth(eg_neurons,i,bin_edges,win_interest,figProp)

    saveas(gcf,fullfile(save_out,['neu',num2str(eg_neurons(i).phyID),'_reach_lessSmooth.png']),'png')
end


%% Order neurons by PC1 and 2 coefficients angular position
% find PCs using center cat1 reaches
% HEYYYYY THIS IS VERY SGTUPID!!!! ONLY THOSE IN RIGHT LEFT CENTER INX ARE
% WITHIN THE TRIALS OF INTEREST, THAT CONDITION IS ENOUGH....
left_cat1_idx = find(ismember(reach_trials_inVec,left_idx) & cat_reach_inVec==1);
center_cat1_idx = find(ismember(reach_trials_inVec,center_idx) & cat_reach_inVec==1);
right_cat1_idx = find(ismember(reach_trials_inVec,right_idx) & cat_reach_inVec==1);

BG_rCat1L_FR_mean = squeeze(mean(FR_reach(:,left_cat1_idx,idx_BG),2,'omitnan'));
BG_rCat1C_FR_mean = squeeze(mean(FR_reach(:,center_cat1_idx,idx_BG),2,'omitnan'));
BG_rCat1R_FR_mean = squeeze(mean(FR_reach(:,right_cat1_idx,idx_BG),2,'omitnan'));

BG_rLCR_mat = cat(1,BG_rCat1L_FR_mean,BG_rCat1C_FR_mean,BG_rCat1R_FR_mean);
[zs_BG_rLCR_FR, BG_FR_mu, BG_FR_sigma] = zscore(BG_rLCR_mat);
[coeff_rLCR_BG, score_rLCR_BG, latent_rLCR_BG, tsquare_rLCR_BG, explained_rLCR_BG, mus_rLCR_BG] = pca(zs_BG_rLCR_FR);

conditions_reach = ['left  ';'center';'right '];
nr_conditions_reach = size(conditions_reach,1);
n_bg_neu = length(idx_BG);
n_ct_neu = length(idx_CT);
score_rLCR_BG_rs = nan(nr_bins,nr_conditions_reach,n_bg_neu);

for j=1:nr_conditions_reach
    score_rLCR_BG_rs(:,j,:) = score_rLCR_BG((nr_bins*(j-1))+1:nr_bins*j,:);
end
% Check scores
lw=2;
figure
ff=tiledlayout(1,3);
title(ff,'Region: BG | All full reaches')
for pc = 1:3
nexttile
plot(bin_edges,squeeze(score_rLCR_BG_rs(:,1,pc)),'linewidth',lw,'Color',behavior.colors.left_color); hold on
plot(bin_edges,squeeze(score_rLCR_BG_rs(:,2,pc)),'linewidth',lw,'Color',behavior.colors.center_color); 
plot(bin_edges,squeeze(score_rLCR_BG_rs(:,3,pc)),'linewidth',lw,'Color',behavior.colors.right_color); 
xlim(win_interest);
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8 .5],'LineWidth',2);
xlabel('time from reach max (s)'); ylabel(strcat('PC',num2str(pc))); 
end
legend(conditions_reach,'box','off')

set(gcf,'position',[2128 489 1637 412])
saveas(gcf,fullfile(save_out,'pca_scores_reachesCat1_BG.png'),'png');


%% Transform coeficients to polar coordinates
r_coef_BG_reach = sqrt(coeff_rLCR_BG(:,1).^2 + coeff_rLCR_BG(:,2).^2);    
theta_coef_BG_reach = atan2(coeff_rLCR_BG(:,1), coeff_rLCR_BG(:,2));     
[theta_coef_BG_reach_sorted,BG_reach_order] = sort(theta_coef_BG_reach);
cmap = parula(length(BG_reach_order));

figure
subplot(121)
scatter(coeff_rLCR_BG(:,1),coeff_rLCR_BG(:,2),20,'filled');
xlabel('PC1 coefficients'); ylabel('PC2 coefficients')
set(gca,axeOpt{:})
title('Position of PC1 vs PC2 coeffs')
axis square

subplot(122)
for i=1:n_bg_neu
    polarplot(theta_coef_BG_reach_sorted(i), r_coef_BG_reach(BG_reach_order(i)),...
        'o','MarkerSize', 8, 'MarkerFaceColor',cmap(i,:),'markerEdgeColor','w'); hold on
end
title('Points in Polar Coordinates');
title('Angular position of PC1 vs PC2 coeffs')
saveas(gcf,fullfile(save_out,'PC1coeffs_vs_PC2coeffs_BG.png'),'png');


%% Save
save(fullfile(save_mat,'eg_neurons.mat'),'eg_neurons','idx_BG','idx_CT','BG_reach_order','-v7.3');


%% Order cells by theta position, different contrasts
zscore_cat1L_BG = zscore(BG_rCat1L_FR_mean);
zscore_cat1C_BG = zscore(BG_rCat1C_FR_mean);
zscore_cat1R_BG = zscore(BG_rCat1R_FR_mean);

% All full reaches, different directions
figure
subplot(131)
imagesc(bin_edges,1:n_bg_neu,zscore_cat1L_BG(:,BG_reach_order)')
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
colormap("hot")
axis xy;
set(gca,axeOpt{:},'TickDir','out')
xlabel('time  from reach max(s)'); ylabel('cells in BG');
title({'Mean FR in BG';'Full reaches to the left'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

subplot(132)
imagesc(bin_edges,1:n_bg_neu,zscore_cat1C_BG(:,BG_reach_order)')
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
colormap("hot")
axis xy;
set(gca,axeOpt{:},'TickDir','out')
xlabel('time  from reach max(s)'); ylabel('cells in BG');
title({'Mean FR in BG';'Full reaches to the center'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

subplot(133)
imagesc(bin_edges,1:n_bg_neu,zscore_cat1R_BG(:,BG_reach_order)')
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
colormap("hot")
axis xy;
set(gca,axeOpt{:},'TickDir','out')
xlabel('time  from reach max(s)'); ylabel('cells in BG');
title({'Mean FR in BG';'Full reaches to the right'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

set(gcf,'position',[1972 212 1870 644])
saveas(gcf,fullfile(save_out,'tile_LCR_full_reaches.png'),'png');

%% Lifted paw reaches
left_cat2_idx = find(ismember(reach_trials_inVec,left_idx) & cat_reach_inVec==2);
center_cat2_idx = find(ismember(reach_trials_inVec,center_idx) & cat_reach_inVec==2);
right_cat2_idx = find(ismember(reach_trials_inVec,right_idx) & cat_reach_inVec==2);

BG_rCat2L_FR_mean = squeeze(mean(FR_reach(:,left_cat2_idx,idx_BG),2,'omitnan'));
BG_rCat2C_FR_mean = squeeze(mean(FR_reach(:,center_cat2_idx,idx_BG),2,'omitnan'));
BG_rCat2R_FR_mean = squeeze(mean(FR_reach(:,right_cat2_idx,idx_BG),2,'omitnan'));

zscore_cat2L_BG = zscore(BG_rCat21L_FR_mean);
zscore_cat2C_BG = zscore(BG_rCat2C_FR_mean);
zscore_cat2R_BG = zscore(BG_rCat2R_FR_mean);

% All full reaches, different directions
figure
subplot(131)
imagesc(bin_edges,1:n_bg_neu,zscore_cat2L_BG(:,BG_reach_order)')
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
colormap("hot")
axis xy;
set(gca,axeOpt{:},'TickDir','out')
xlabel('time  from reach max(s)'); ylabel('cells in BG');
title({'Mean FR in BG';'Lifted paw reaches to the left'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

subplot(132)
imagesc(bin_edges,1:n_bg_neu,zscore_cat2C_BG(:,BG_reach_order)')
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
colormap("hot")
axis xy;
set(gca,axeOpt{:},'TickDir','out')
xlabel('time  from reach max(s)'); ylabel('cells in BG');
title({'Mean FR in BG';'Lifted paw reaches to the center'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

subplot(133)
imagesc(bin_edges,1:n_bg_neu,zscore_cat2R_BG(:,BG_reach_order)')
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
colormap("hot")
axis xy;
set(gca,axeOpt{:},'TickDir','out')
xlabel('time  from reach max(s)'); ylabel('cells in BG');
title({'Mean FR in BG';'Lifted paw to the right'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

set(gcf,'position',[1972 212 1870 644])
saveas(gcf,fullfile(save_out,'tile_LCR_full_reaches.png'),'png');

%% Hit vs miss


%% Success vs fail

%% Repeat for CTX


%%




