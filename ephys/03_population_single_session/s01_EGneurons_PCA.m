% Get extra-good neurons, calculate instantaneous FR, order neurons
clear; close all; clc

%% Load data
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

% Session & ephys source
sess = 'R4';
reg = 'CB';
sorter_folder = 'catGT\kilosort4';
ephys_local_folder = 1;

if (strcmp(reg,'BG') || strcmp(reg,'CT')), region = 'BG';
elseif strcmp(reg,'CB'), region = 'CB'; end


%% paths
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
behavior_dir = fullfile(rootdir,"behavior_data","raw_data",group_behav,setup,mouse,sess);
reaching_dir = fullfile(rootdir,"behavior_data","analyzed_data","mat_files",group_behav,setup,mouse,sess);
bhv_dir = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess);
ephys_eg_dir = fullfile(rootdir,"ephys_and_behavior/","mat_files/",group,mouse,sess,region);

% Load mat files
load(fullfile(behavior_dir,"behavior_session.mat"));
load(fullfile(bhv_dir,"behavior_fundamentals.mat"));
load(fullfile(reaching_dir,"session_reaching_data_paw.mat"));
load(fullfile(ephys_eg_dir,"eg_neurons.mat"));

% Check name of idx
if strcmp(reg,'BG'), idx_reg = idx_BG;
elseif  strcmp(reg,'CT'), idx_reg = idx_CT;
elseif strcmp(reg,'CB'), idx_reg = 1:length(eg_neurons);
end

% % Save output
save_mat = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess,reg);
save_out = fullfile(rootdir,"ephys_and_behavior","out_files",group,mouse,sess,reg,'PCA_scores');
if ~exist(save_out,"dir"), mkdir(save_out); end
if ~exist(save_mat,"dir"), mkdir(save_mat); end

%% SMOOTHING PARAMS
bin_width = eg_neu_FR_params.bin_width;
sig_pk = .02;
k_pc = gausskernel('sig',sig_pk,'binwidth',bin_width);
bin_edges = eg_neu_FR_params.bin_edges;
nr_bins = eg_neu_FR_params.nr_bins;
bin_edges_sm = bin_edges(1)-k_pc.paddx(1):bin_width:bin_edges(end)-k_pc.paddx(2)+bin_width;
nr_bins_sm = length(bin_edges_sm);

% FigureS
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};

clr_push = behavior.colors.push_clr;
clr_pull = behavior.colors.pull_clr;

clr_left = behavior.colors.left_color;
clr_center = behavior.colors.center_color;
clr_right = behavior.colors.right_color;

%% PCA init
push_idx = bhv.push_idx;
pull_idx = bhv.pull_idx;
push_all_idx = cat(1,push_idx,find(bhv.init_invPush_invPull_idx==3));
pull_all_idx = cat(1,pull_idx,find(bhv.init_invPush_invPull_idx==2));
n_neu = length(idx_reg);

push_idx_forPCA = push_idx;
pull_idx_forPCA = pull_idx;

FR_mean_push = squeeze(mean(FR_init(:,push_idx_forPCA,idx_reg),2,'omitnan'));
FR_mean_pull = squeeze(mean(FR_init(:,pull_idx_forPCA,idx_reg),2,'omitnan'));

FR_pp_mat = cat(1,FR_mean_push,FR_mean_pull);
[zs_FR_pp, FR_pp_mu, FR_pp_sigma] = zscore(FR_pp_mat);
[coeff_pp, score_pp, latent_pp, tsquare_pp, explained_pp, mus_pp] = pca(zs_FR_pp);

conditions_pp = ['push';'pull'];
nr_conditions_pp = size(conditions_pp,1);
score_pp_rs = nan(nr_bins,nr_conditions_pp,n_neu);
for j=1:nr_conditions_pp
    score_pp_rs(:,j,:) = score_pp((nr_bins*(j-1))+1:nr_bins*j,:);
end

% smothed win
score_pp_sm = nan(nr_bins_sm,nr_conditions_pp,n_neu);
for i = 1:nr_conditions_pp
    score_pp_sm(:,i,:) = conv2(k_pc.pdf,1,squeeze(score_pp_rs(:,i,:)),"valid");
end

% Order of cells
r_coeff_pp = sqrt(coeff_pp(:,1).^2 + coeff_pp(:,2).^2);    
theta_coef_pp = atan2(coeff_pp(:,2), coeff_pp(:,1));     
[theta_coeffPP_sorted, neu_order_pp] = sort(theta_coef_pp);
cmap = parula(length(neu_order_pp));


%% Check scores
win_interest = [-3 3];
lw=2;
figure
ff=tiledlayout(3,1);
title(ff,sprintf('%s%s%s','Region: ',reg,' | Trial init'));
for pc = 1:3
nexttile
plot(bin_edges,squeeze(score_pp_rs(:,1,pc)),'linewidth',lw,'Color',clr_push); hold on
plot(bin_edges,squeeze(score_pp_rs(:,2,pc)),'linewidth',lw,'Color',clr_pull); 
xlim(win_interest);
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8 .5],'LineWidth',2);
xlabel('time from trial init (s)'); ylabel(strcat('PC',num2str(pc))); 
end
legend(conditions_pp,'box','off')

set(gcf,'position',[2661 175 533 758])
saveas(gcf,fullfile(save_out,'pca_scores_push_pull_val.png'),'png');

% Smoothed scores
win_interest = [-3 3];
lw=2;
figure
ff=tiledlayout(3,1);
title(ff,sprintf('%s%s%s','Region: ',reg,' | Trial init'));
for pc = 1:3
nexttile
plot(bin_edges_sm,squeeze(score_pp_sm(:,1,pc)),'linewidth',lw,'Color',clr_push); hold on
plot(bin_edges_sm,squeeze(score_pp_sm(:,2,pc)),'linewidth',lw,'Color',clr_pull); 
xlim(win_interest);
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8 .5],'LineWidth',2);
xlabel('time from trial init (s)'); ylabel(strcat('PC',num2str(pc))); 
end
legend(conditions_pp,'box','off')

set(gcf,'position',[2661 175 533 758])
saveas(gcf,fullfile(save_out,'pca_scores_push_pull_val_smooth.png'),'png');
print(gcf, fullfile(save_out, 'pca_scores_push_pull_val_smooth.pdf'), '-dpdf', '-painters');


%% Check sorting order
%for i=1:n_neu
figure
subplot(121)
scatter(coeff_pp(:,1),coeff_pp(:,2),50,cmap,'filled');
xlabel('PC1 coefficients'); ylabel('PC2 coefficients')
set(gca,axeOpt{:})
title('Position of PC1 vs PC2 coeffs')
axis square
hold on

for i=1:n_neu
subplot(122)
    polarplot(theta_coeffPP_sorted(i), r_coeff_pp(neu_order_pp(i)),...
        'o','MarkerSize', 8, 'MarkerFaceColor',cmap(i,:),'markerEdgeColor','w'); hold on
end
colorbar
%title('Points in Polar Coordinates');
title('Angular position of PC1 vs PC2 coeffs')
set(gcf,'Position',[2762         134        1027         734])
saveas(gcf,fullfile(save_out,'PC1coeffs_vs_PC2coeffs_pushPull.png'),'png');


% ------------------------------------------------------------------------------
%% PCA reach: full reach -------------------------------------------------------
% Sanity check - reaches keep their identities
figure
plot(squeeze(bhv.reaches_inVec_px(:,1,bhv.cat_reach_inVec==1)))

%%
left_r1_idx = find(ismember(bhv.reach_trials_inVec,bhv.left_idx) & bhv.cat_reach_inVec==1);
center_r1_idx = find(ismember(bhv.reach_trials_inVec,bhv.center_idx) & bhv.cat_reach_inVec==1);
right_r1_idx = find(ismember(bhv.reach_trials_inVec,bhv.right_idx) & bhv.cat_reach_inVec==1);

r1L_FR_mean = squeeze(mean(FR_reach(:,left_r1_idx,idx_reg),2,'omitnan'));
r1C_FR_mean = squeeze(mean(FR_reach(:,center_r1_idx,idx_reg),2,'omitnan'));
r1R_FR_mean = squeeze(mean(FR_reach(:,right_r1_idx,idx_reg),2,'omitnan'));

FR_rLCR_mat = cat(1,r1L_FR_mean,r1C_FR_mean,r1R_FR_mean);
[zs_FR_rLCR, FR_rLCR_mu, FR_rLCR_sigma] = zscore(FR_rLCR_mat);
[coeff_rLCR, score_rLCR, latent_rLCR, tsquare_rLCR, explained_rLCR, mus_rLCR] = pca(zs_FR_rLCR);

conditions_reach = ['left  ';'center';'right '];
nr_conditions_reach = size(conditions_reach,1);


score_rLCR_rs = nan(nr_bins,nr_conditions_reach,n_neu);
for j=1:nr_conditions_reach
    score_rLCR_rs(:,j,:) = score_rLCR((nr_bins*(j-1))+1:nr_bins*j,:);
end

score_rLCR_sm = nan(nr_bins_sm,nr_conditions_reach,n_neu);
for i = 1:nr_conditions_reach
    score_rLCR_sm(:,i,:) = conv2(k_pc.pdf,1,squeeze(score_rLCR_rs(:,i,:)),"valid");
end

% Order of cells
r_coeff_rLCR = sqrt(coeff_rLCR(:,1).^2 + coeff_rLCR(:,2).^2);    
theta_coef_rLCR = atan2(coeff_rLCR(:,2), coeff_rLCR(:,1));     
[theta_coeffrLCR_sorted, neu_order_rLCR] = sort(theta_coef_rLCR);
cmap_rLCR = parula(length(neu_order_rLCR));

%------------------------------------------------
%% Check scores
win_interest = [-3 3];
lw=2;
figure
ff=tiledlayout(3,1);
title(ff,sprintf('%s%s%s','Region: ',reg,' | Full reaches'));
for pc = 1:3
nexttile
plot(bin_edges,squeeze(score_rLCR_rs(:,1,pc)),'linewidth',lw,'Color',clr_left); hold on
plot(bin_edges,squeeze(score_rLCR_rs(:,2,pc)),'linewidth',lw,'Color',clr_center); 
plot(bin_edges,squeeze(score_rLCR_rs(:,3,pc)),'linewidth',lw,'Color',clr_right); 
xlim(win_interest);
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8 .5],'LineWidth',2);
xlabel('time from reach endpoint (s)'); ylabel(strcat('PC',num2str(pc))); 
end
legend(conditions_reach,'box','off')

set(gcf,'position',[2661 175 533 758])
saveas(gcf,fullfile(save_out,'pca_scores_fullReaches.png'),'png');

% Smoothed scores
win_interest = [-3 3];
lw=2;
figure
ff=tiledlayout(3,1);
title(ff,sprintf('%s%s%s','Region: ',reg,' | Full reaches'));
for pc = 1:3
nexttile
plot(bin_edges_sm,squeeze(score_rLCR_sm(:,1,pc)),'linewidth',lw,'Color',clr_left); hold on
plot(bin_edges_sm,squeeze(score_rLCR_sm(:,2,pc)),'linewidth',lw,'Color',clr_center); 
plot(bin_edges_sm,squeeze(score_rLCR_sm(:,3,pc)),'linewidth',lw,'Color',clr_right); 

xlim(win_interest);
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8 .5],'LineWidth',2);
xlabel('time from reach endpoint (s)'); ylabel(strcat('PC',num2str(pc))); 
end
legend(conditions_reach,'box','off')

set(gcf,'position',[2661 175 533 758])
saveas(gcf,fullfile(save_out,'pca_scores_fullReaches_smooth.png'),'png');


%% Check sorting order
%for i=1:n_neu
figure
subplot(121)
scatter(coeff_rLCR(:,1),coeff_rLCR(:,2),50,cmap,'filled');
xlabel('PC1 coefficients'); ylabel('PC2 coefficients')
set(gca,axeOpt{:})
title('Position of PC1 vs PC2 coeffs')
axis square
hold on

for i=1:n_neu
subplot(122)
    polarplot(theta_coeffrLCR_sorted(i), r_coeff_rLCR(neu_order_rLCR(i)),...
        'o','MarkerSize', 8, 'MarkerFaceColor',cmap(i,:),'markerEdgeColor','w'); hold on
end
colorbar
%title('Points in Polar Coordinates');
title('Angular position of PC1 vs PC2 coeffs')
set(gcf,'Position',[2762         134        1027         734])
saveas(gcf,fullfile(save_out,'PC1coeffs_vs_PC2_fullReaches.png'),'png');

% ------------------------------------------------------------------------------
%% PCA reach: hit reach -------------------------------------------------------
%Sanity check - reaches keep their identities
figure
plot(squeeze(bhv.reaches_inVec_px(:,1,bhv.hit_inVec==1)))
%%
left_hit_idx = find(ismember(bhv.reach_trials_inVec,bhv.left_idx) & bhv.hit_inVec==1);
center_hit_idx = find(ismember(bhv.reach_trials_inVec,bhv.center_idx) & bhv.hit_inVec==1);
right_hit_idx = find(ismember(bhv.reach_trials_inVec,bhv.right_idx) & bhv.hit_inVec==1);

hL_FR_mean = squeeze(mean(FR_reach(:,left_hit_idx,idx_reg),2,'omitnan'));
hC_FR_mean = squeeze(mean(FR_reach(:,center_hit_idx,idx_reg),2,'omitnan'));
hR_FR_mean = squeeze(mean(FR_reach(:,right_hit_idx,idx_reg),2,'omitnan'));

FR_hLCR_mat = cat(1,hL_FR_mean,hC_FR_mean,hR_FR_mean);
[zs_FR_hLCR, FR_hLCR_mu, FR_hLCR_sigma] = zscore(FR_hLCR_mat);
[coeff_hLCR, score_hLCR, latent_hLCR, tsquare_hLCR, explained_hLCR, mus_hLCR] = pca(zs_FR_hLCR);

score_hLCR_rs = nan(nr_bins,nr_conditions_reach,n_neu);
for j=1:nr_conditions_reach
    score_hLCR_rs(:,j,:) = score_hLCR((nr_bins*(j-1))+1:nr_bins*j,:);
end

score_hLCR_sm = nan(nr_bins_sm,nr_conditions_reach,n_neu);
for i = 1:nr_conditions_reach
    score_hLCR_sm(:,i,:) = conv2(k_pc.pdf,1,squeeze(score_hLCR_rs(:,i,:)),"valid");
end

% Order of cells
r_coeff_hLCR = sqrt(coeff_hLCR(:,1).^2 + coeff_hLCR(:,2).^2);    
theta_coef_hLCR = atan2(coeff_hLCR(:,2), coeff_hLCR(:,1));     
[theta_coeffhLCR_sorted, neu_order_hLCR] = sort(theta_coef_hLCR);
cmap_hLCR = parula(length(neu_order_hLCR));


%------------------------------------------------
%% Check scores
win_interest = [-3 3];
lw=2;
figure
ff=tiledlayout(3,1);
title(ff,sprintf('%s%s%s','Region: ',reg,' | Hit reaches'));
for pc = 1:3
nexttile
plot(bin_edges,squeeze(score_hLCR_rs(:,1,pc)),'linewidth',lw,'Color',clr_left); hold on
plot(bin_edges,squeeze(score_hLCR_rs(:,2,pc)),'linewidth',lw,'Color',clr_center); 
plot(bin_edges,squeeze(score_hLCR_rs(:,3,pc)),'linewidth',lw,'Color',clr_right); 
xlim(win_interest);
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8 .5],'LineWidth',2);
xlabel('time from reach endpoint (s)'); ylabel(strcat('PC',num2str(pc))); 
end
legend(conditions_reach,'box','off')

set(gcf,'position',[2661 175 533 758])
saveas(gcf,fullfile(save_out,'pca_scores_hitReaches.png'),'png');

% Smoothed scores
win_interest = [-3 3];
lw=2;
figure
ff=tiledlayout(3,1);
title(ff,sprintf('%s%s%s','Region: ',reg,' | Hit reaches'));
for pc = 1:3
nexttile
plot(bin_edges_sm,squeeze(score_hLCR_sm(:,1,pc)),'linewidth',lw,'Color',clr_left); hold on
plot(bin_edges_sm,squeeze(score_hLCR_sm(:,2,pc)),'linewidth',lw,'Color',clr_center); 
plot(bin_edges_sm,squeeze(score_hLCR_sm(:,3,pc)),'linewidth',lw,'Color',clr_right); 

xlim(win_interest);
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8 .5],'LineWidth',2);
xlabel('time from reach endpoint (s)'); ylabel(strcat('PC',num2str(pc))); 
end
legend(conditions_reach,'box','off')

set(gcf,'position',[2661 175 533 758])
saveas(gcf,fullfile(save_out,'pca_scores_hitReaches_smooth.png'),'png');
print(gcf, fullfile(save_out, 'pca_scores_hitReaches_smooth.pdf'), '-dpdf', '-painters');


%% Check sorting order
%for i=1:n_neu
figure
subplot(121)
scatter(coeff_hLCR(:,1),coeff_hLCR(:,2),50,cmap,'filled');
xlabel('PC1 coefficients'); ylabel('PC2 coefficients')
set(gca,axeOpt{:})
title('Position of PC1 vs PC2 coeffs')
axis square
hold on

for i=1:n_neu
subplot(122)
    polarplot(theta_coeffhLCR_sorted(i), r_coeff_hLCR(neu_order_hLCR(i)),...
        'o','MarkerSize', 8, 'MarkerFaceColor',cmap(i,:),'markerEdgeColor','w'); hold on
end
colorbar
%title('Points in Polar Coordinates');
title('Angular position of PC1 vs PC2 coeffs')
set(gcf,'Position',[2762         134        1027         734])
saveas(gcf,fullfile(save_out,'PC1coeffs_vs_PC2_fullReaches.png'),'png');


%% Save
save(fullfile(save_mat,['pca_',reg,'.mat']),'bin_edges',...
    'bin_edges_sm','k_pc',...
    'FR_pp_mat','zs_FR_pp','conditions_pp',...
    'coeff_pp','score_pp','explained_pp',...
    'score_pp_rs','score_pp_sm','neu_order_pp',...
    'FR_rLCR_mat','zs_FR_rLCR','conditions_reach',...
    'coeff_rLCR','score_rLCR','explained_rLCR',...
    'score_rLCR_rs','score_rLCR_sm','neu_order_rLCR',...
    'FR_hLCR_mat','zs_FR_hLCR',...
    'coeff_hLCR','score_hLCR','explained_hLCR',...
    'score_hLCR_rs','score_hLCR_sm','neu_order_hLCR');
fprintf('Done!!')



