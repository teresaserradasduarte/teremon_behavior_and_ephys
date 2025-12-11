%% Simple preliminary ugly single session population analysis
clear; close all; clc

%% Load
% ehpys
ephys_root = 'E:\'; %group_ephys = '20230801_ChocolateGroup';
ephys_sess = '26082023_Lindt_StrCer_S4_g0';

% Behavior
behavior_root ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
group_setup_behav = strcat('20230511_ChocolateGroup',filesep,'headfixed_dynamicTarget');
mouse = '2_Lindt';
session = 'R4';
imec = '1';

% Load behavior
behavior_path = strcat(behavior_root,filesep,group_setup_behav,filesep,mouse,filesep,session);
load(strcat(behavior_path,filesep,'behavior_session.mat'));

% Load ehpys
imec_path = strcat(ephys_root,filesep,ephys_sess,filesep,ephys_sess,'_imec',imec,filesep,'catGT',filesep,'kilosort4');
neurons_imec = load(strcat(imec_path,filesep,'neurons_session.mat'));
out_folder = strcat(imec_path,filesep,'sess_out_figures');
if ~exist('out_folder','dir'), mkdir(out_folder); end

%% Load neuron table
neu = struct2table(neurons_imec.neurons);
eg_idx = find(neu.extraGood == 1);

% Dims
%nr_trials = behavior.behavior_duration.trial_end;
nr_trials = size(neu.spk_rates_init{1,1},2);
nr_time_bins = size(neu.spk_rates_init{1,1},1);
nr_eg_neurons = length(eg_idx);

init_rates = zeros(nr_eg_neurons,nr_time_bins,nr_trials);
reach_rates = zeros(nr_eg_neurons,nr_time_bins,nr_trials);

for i = 1:length(eg_idx)
    init_rates(i,:,:) = neu{eg_idx,"spk_rates_init"}{i,1};
    reach_rates(i,:,:) = neu{eg_idx,"spk_rates_reach"}{i,1};
end
time_bins = neurons_imec.neurons_params.gamma.spk_bins;


%%  trial-concatenated PCA
pca_mat_init = init_rates(:,:)';
pca_mat_reach = reach_rates(:,:)';

% z-score
[zXinit, zmu_in, zsigma_in] = zscore(pca_mat_init);
[zXreach, zmu_re, zsigma_re] = zscore(pca_mat_reach);

% PCA
[coeff_init, score_init, latent_init, tsquare_init, explained_init, mus_init] = pca(zXinit);
[coeff_reach, score_reach, latent_reach, tsquare_reach, explained_reach, mus_reach] = pca(zXreach);
L = 10;
coeff_L_init = coeff_init(:,1:L);
coeff_L_reach = coeff_reach(:,1:L);
explained_L_init = explained_init(1:L);


score_init_rs = zeros(nr_time_bins,nr_trials,nr_eg_neurons);
score_reach_rs = zeros(nr_time_bins,nr_trials,nr_eg_neurons);
for j=1:nr_trials
    score_init_rs(:,j,:) = score_init((nr_time_bins*(j-1))+1:nr_time_bins*j,:);
    score_reach_rs(:,j,:) = score_reach((nr_time_bins*(j-1))+1:nr_time_bins*j,:);
end

% Condition average + std
push_idx = behavior.init.idx_trial_push(behavior.init.idx_trial_push<=nr_trials);
pull_idx = behavior.init.idx_trial_pull(behavior.init.idx_trial_pull<=nr_trials);
left_idx = behavior.reach.left_idx(behavior.reach.left_idx<=nr_trials);
right_idx = behavior.reach.right_idx(behavior.reach.right_idx<=nr_trials);
center_idx = behavior.reach.center_idx(behavior.reach.center_idx<=nr_trials);

push_init_scores_mean = squeeze(mean(score_init_rs(:,push_idx,1:L),2));
pull_init_scores_mean = squeeze(mean(score_init_rs(:,pull_idx,1:L),2));
push_init_scores_sd = squeeze(std(score_init_rs(:,push_idx,1:L),[],2));%./length(push_idx);
pull_init_scores_sd = squeeze(std(score_init_rs(:,pull_idx,1:L),[],2));%./length(pull_idx);

left_reach_scores_mean = squeeze(mean(score_reach_rs(:,left_idx,1:L),2));
right_reach_scores_mean = squeeze(mean(score_reach_rs(:,right_idx,1:L),2));
center_reach_scores_mean = squeeze(mean(score_reach_rs(:,center_idx,1:L),2));
left_reach_scores_sd = squeeze(std(score_reach_rs(:,left_idx,1:L),[],2));%./length(left_idx);
right_reach_scores_sd = squeeze(std(score_reach_rs(:,right_idx,1:L),[],2));%./length(right_idx);
center_reach_scores_sd = squeeze(std(score_reach_rs(:,center_idx,1:L),[],2));%./length(center_idx);

% Save
% init
trial_concat_pca.init.pca_mat_init = pca_mat_init;
trial_concat_pca.init.zXinit = zXinit;
trial_concat_pca.init.zmu_in = zmu_in;
trial_concat_pca.init.zsigma_in = zsigma_in;
trial_concat_pca.init.coeff_init = coeff_init;
trial_concat_pca.init.score_init = score_init;
trial_concat_pca.init.latent_init = latent_init;
trial_concat_pca.init.tsquare_init = tsquare_init;
trial_concat_pca.init.explained_init = explained_init;
trial_concat_pca.init.mus_init = mus_init;
trial_concat_pca.init.score_init_rs = score_init_rs;
trial_concat_pca.init.push_init_scores_mean = push_init_scores_mean;
trial_concat_pca.init.pull_init_scores_mean = pull_init_scores_mean;
trial_concat_pca.init.push_init_scores_sd = push_init_scores_sd;
trial_concat_pca.init.pull_init_scores_sd = pull_init_scores_sd;
trial_concat_pca.init.left_reach_scores_mean = left_reach_scores_mean;
trial_concat_pca.init.right_reach_scores_mean = right_reach_scores_mean;
trial_concat_pca.init.center_reach_scores_mean = center_reach_scores_mean;
trial_concat_pca.init.left_reach_scores_sd = left_reach_scores_sd;
trial_concat_pca.init.right_reach_scores_sd = right_reach_scores_sd;
trial_concat_pca.init.center_reach_scores_sd = center_reach_scores_sd;
% reach
trial_concat_pca.reach.pca_mat_reach = pca_mat_reach;
trial_concat_pca.reach.zXreach = zXreach;
trial_concat_pca.reach.zmu_re = zmu_re;
trial_concat_pca.reach.zsigma_re = zsigma_re;
trial_concat_pca.reach.coeff_reach = coeff_reach;
trial_concat_pca.reach.latent_reach = latent_reach;
trial_concat_pca.reach.tsquare_reach = tsquare_reach;
trial_concat_pca.reach.explained_reach = explained_reach;
trial_concat_pca.reach.mus_reach = mus_reach;
trial_concat_pca.reach.score_reach_rs = score_reach_rs;
trial_concat_pca.reach.push_init_scores_mean = push_init_scores_mean;
trial_concat_pca.reach.pull_init_scores_mean = pull_init_scores_mean;
trial_concat_pca.reach.push_init_scores_sd = push_init_scores_sd;
trial_concat_pca.reach.pull_init_scores_sd = pull_init_scores_sd;
trial_concat_pca.reach.left_reach_scores_mean = left_reach_scores_mean;
trial_concat_pca.reach.right_reach_scores_mean = right_reach_scores_mean;
trial_concat_pca.reach.center_reach_scores_mean = center_reach_scores_mean;
trial_concat_pca.reach.left_reach_scores_sd = left_reach_scores_sd;
trial_concat_pca.reach.right_reach_scores_sd = right_reach_scores_sd;
trial_concat_pca.reach.center_reach_scores_sd = center_reach_scores_sd;



%% Figure
axeOpt=neurons_imec.figProp.axeOpt;
lw = 1.5;
patchSat = 0.1;

figure
subplot(321)
shadedErrorBar(time_bins,push_init_scores_mean(:,1),push_init_scores_sd(:,1),'lineProps',{'Color',behavior.colors.push_clr,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
hold on
shadedErrorBar(time_bins,pull_init_scores_mean(:,1),pull_init_scores_sd(:,1),'lineProps',{'Color',behavior.colors.pull_clr,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat);
set(gca,axeOpt{:})
xlabel('time from init (s)'); ylabel('PC1'); 

subplot(322)
shadedErrorBar(time_bins,left_reach_scores_mean(:,1),left_reach_scores_sd(:,1),'lineProps',{'Color',behavior.colors.left_color,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
hold on
shadedErrorBar(time_bins,right_reach_scores_mean(:,1),right_reach_scores_sd(:,1),'lineProps',{'Color',behavior.colors.right_color,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
shadedErrorBar(time_bins,center_reach_scores_mean(:,1),center_reach_scores_sd(:,1),'lineProps',{'Color',behavior.colors.center_color,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
set(gca,axeOpt{:})
xlabel('time from collection (s)'); ylabel('PC1'); 

subplot(323)
shadedErrorBar(time_bins,push_init_scores_mean(:,2),push_init_scores_sd(:,2),'lineProps',{'Color',behavior.colors.push_clr,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
hold on
shadedErrorBar(time_bins,pull_init_scores_mean(:,2),pull_init_scores_sd(:,2),'lineProps',{'Color',behavior.colors.pull_clr,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat);
set(gca,axeOpt{:})
xlabel('time from init (s)'); ylabel('PC2'); 

subplot(324)
shadedErrorBar(time_bins,left_reach_scores_mean(:,2),left_reach_scores_sd(:,2),'lineProps',{'Color',behavior.colors.left_color,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
hold on
shadedErrorBar(time_bins,right_reach_scores_mean(:,2),right_reach_scores_sd(:,2),'lineProps',{'Color',behavior.colors.right_color,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
shadedErrorBar(time_bins,center_reach_scores_mean(:,2),center_reach_scores_sd(:,2),'lineProps',{'Color',behavior.colors.center_color,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
set(gca,axeOpt{:})
xlabel('time from collection (s)'); ylabel('PC2'); 

subplot(325)
shadedErrorBar(time_bins,push_init_scores_mean(:,3),push_init_scores_sd(:,3),'lineProps',{'Color',behavior.colors.push_clr,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
hold on
shadedErrorBar(time_bins,pull_init_scores_mean(:,3),pull_init_scores_sd(:,3),'lineProps',{'Color',behavior.colors.pull_clr,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat);
set(gca,axeOpt{:})
xlabel('time from init (s)'); ylabel('PC3'); 

subplot(326)
shadedErrorBar(time_bins,left_reach_scores_mean(:,3),left_reach_scores_sd(:,3),'lineProps',{'Color',behavior.colors.left_color,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
hold on
shadedErrorBar(time_bins,right_reach_scores_mean(:,3),right_reach_scores_sd(:,3),'lineProps',{'Color',behavior.colors.right_color,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
shadedErrorBar(time_bins,center_reach_scores_mean(:,3),center_reach_scores_sd(:,3),'lineProps',{'Color',behavior.colors.center_color,'LineWidth',lw},'transparent',1,'patchSaturation',patchSat); 
set(gca,axeOpt{:})
xlabel('time from collection (s)'); ylabel('PC3');
set(gcf,'color','w')
saveas(gcf,strcat(out_folder,filesep,'trial_concat_pca.png'),'png');

figure
plot(cumsum(explained_init),'linewidth',lw); hold on
plot(cumsum(explained_reach),'linewidth',lw);
legend('init','reach','box','off','location','southeast')
set(gca,axeOpt{:})
xlabel('PCs'); ylabel('variance explained trial concatenated');
set(gcf,'color','w')
saveas(gcf,strcat(out_folder,filesep,'variance_explained_trial_concat.png'),'png');

%% Save
save(strcat(imec_path,filesep,'neurons_session.mat'),'trial_concat_pca',...
    '-append');


