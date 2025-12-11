%% Simple preliminary ugly single session population analysis
clear; close all; clc

%% Load
% ehpys
ephys_root = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_data\20230801_ChocolateGroup'; %group_ephys = '20230801_ChocolateGroup';
behavior_root ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
group_setup_behav = strcat('20230511_ChocolateGroup',filesep,'headfixed_dynamicTarget');
mouse = '4_Milka';
session = 'R4';
ephys_sess = '18082023_Milka_StrCer_S4_g0';

imec = '1';

% Load behavior
behavior_path = strcat(behavior_root,filesep,group_setup_behav,filesep,mouse,filesep,session);
load(strcat(behavior_path,filesep,'behavior_session.mat'));

% Load ehpys
imec_path = strcat(ephys_root,filesep,mouse,filesep,ephys_sess,filesep,ephys_sess,'_imec',imec,filesep,'catGT',filesep,'kilosort4');
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

%%  trial-averaged PCA

% Condition trials
push_idx = behavior.init.idx_trial_push(behavior.init.idx_trial_push<=nr_trials);
pull_idx = behavior.init.idx_trial_pull(behavior.init.idx_trial_pull<=nr_trials);
left_idx = behavior.reach.left_idx(behavior.reach.left_idx<=nr_trials);
right_idx = behavior.reach.right_idx(behavior.reach.right_idx<=nr_trials);
center_idx = behavior.reach.center_idx(behavior.reach.center_idx<=nr_trials);

% Mean by condition
% push / pull
push_init_mean = mean(init_rates(:,:,push_idx),3);
pull_init_mean = mean(init_rates(:,:,pull_idx),3);
% left / center / right
left_reach_mean = mean(reach_rates(:,:,left_idx),3);
center_reach_mean = mean(reach_rates(:,:,center_idx),3);
right_reach_mean = mean(reach_rates(:,:,right_idx),3);

% PCA input matrix
pca_init_mat = cat(2,push_init_mean,pull_init_mean)';
pca_reach_mat = cat(2,left_reach_mean,center_reach_mean,right_reach_mean)';

% z-score
[zXinit, zmu_in, zsigma_in] = zscore(pca_init_mat);
[zXreach, zmu_re, zsigma_re] = zscore(pca_reach_mat);

% PCA
[coeff_init, score_init, latent_init, tsquare_init, explained_init, mus_init] = pca(zXinit);
[coeff_reach, score_reach, latent_reach, tsquare_reach, explained_reach, mus_reach] = pca(zXreach);
L = 10;
coeff_L_init = coeff_init(:,1:L);
coeff_L_reach = coeff_reach(:,1:L);
explained_L_init = explained_init(1:L);

%% Divide in conditions
nr_conditions_init = 2;
nr_conditions_reach = 3;

score_init_rs = zeros(nr_time_bins,nr_conditions_init,nr_eg_neurons);
score_reach_rs = zeros(nr_time_bins,nr_conditions_reach,nr_eg_neurons);
for j=1:nr_conditions_init
    score_init_rs(:,j,:) = score_init((nr_time_bins*(j-1))+1:nr_time_bins*j,:);
end
for j=1:nr_conditions_reach
    score_reach_rs(:,j,:) = score_reach((nr_time_bins*(j-1))+1:nr_time_bins*j,:);
end

conditions_init = ['push';'pull'];
conditions_reach = ['left  ';'center';'right '];

%% Smoothing

% time win
time_bins = neurons_imec.neurons_params.gamma.spk_bins;
[~,init_idx]=min(abs(time_bins-0));

% gamma for PSTH
bin_width = .002;
peak_x = .05;
k_pc = gammakernel('peakx',peak_x,'binwidth',bin_width);

% smothed win
st_bins = time_bins(1)-k_pc.paddx(1):bin_width:time_bins(end)-k_pc.paddx(2)+bin_width;
nr_bins = length(st_bins);

score_init_sm = zeros(nr_bins,nr_conditions_init,nr_eg_neurons);
for i = 1:nr_conditions_init
    score_init_sm(:,i,:) = conv2(k_pc.pdf,1,squeeze(score_init_rs(:,i,:)),"valid");
end

score_reach_sm = zeros(nr_bins,nr_conditions_reach,nr_eg_neurons);
for i = 1:nr_conditions_reach
    score_reach_sm(:,i,:) = conv2(k_pc.pdf,1,squeeze(score_reach_rs(:,i,:)),"valid");
end

       

%% Save
% init
trial_avg_pca.init.pca_mat_init = pca_init_mat;
trial_avg_pca.init.zXinit = zXinit;
trial_avg_pca.init.zmu_in = zmu_in;
trial_avg_pca.init.zsigma_in = zsigma_in;
trial_avg_pca.init.coeff_init = coeff_init;
trial_avg_pca.init.score_init = score_init;
trial_avg_pca.init.latent_init = latent_init;
trial_avg_pca.init.tsquare_init = tsquare_init;
trial_avg_pca.init.explained_init = explained_init;
trial_avg_pca.init.mus_init = mus_init;
trial_avg_pca.init.score_init_rs = score_init_rs;
trial_avg_pca.init.conditions_init = conditions_init;

% reach
trial_avg_pca.reach.pca_mat_reach = pca_reach_mat;
trial_avg_pca.reach.zXreach = zXreach;
trial_avg_pca.reach.zmu_re = zmu_re;
trial_avg_pca.reach.zsigma_re = zsigma_re;
trial_avg_pca.reach.coeff_reach = coeff_reach;
trial_avg_pca.reach.latent_reach = latent_reach;
trial_avg_pca.reach.tsquare_reach = tsquare_reach;
trial_avg_pca.reach.explained_reach = explained_reach;
trial_avg_pca.reach.mus_reach = mus_reach;
trial_avg_pca.reach.score_reach_rs = score_reach_rs;
trial_avg_pca.reach.conditions_reach = conditions_reach;


%% Figure
axeOpt=neurons_imec.figProp.axeOpt;
lw = 2;

figure
subplot(321)
plot(time_bins,squeeze(score_init_rs(:,1,1)),'linewidth',lw,'Color',behavior.colors.push_clr); hold on
plot(time_bins,squeeze(score_init_rs(:,2,1)),'linewidth',lw,'Color',behavior.colors.pull_clr); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from init (s)'); ylabel('PC1'); 

subplot(322)
plot(time_bins,squeeze(score_reach_rs(:,1,1)),'linewidth',lw,'Color',behavior.colors.left_color); hold on
plot(time_bins,squeeze(score_reach_rs(:,2,1)),'linewidth',lw,'Color',behavior.colors.center_color); 
plot(time_bins,squeeze(score_reach_rs(:,3,1)),'linewidth',lw,'Color',behavior.colors.right_color); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from collection (s)'); ylabel('PC1'); 

subplot(323)
plot(time_bins,squeeze(score_init_rs(:,1,2)),'linewidth',lw,'Color',behavior.colors.push_clr); hold on
plot(time_bins,squeeze(score_init_rs(:,2,2)),'linewidth',lw,'Color',behavior.colors.pull_clr); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from init (s)'); ylabel('PC2'); 

subplot(324)
plot(time_bins,squeeze(score_reach_rs(:,1,2)),'linewidth',lw,'Color',behavior.colors.left_color); hold on
plot(time_bins,squeeze(score_reach_rs(:,2,2)),'linewidth',lw,'Color',behavior.colors.center_color); 
plot(time_bins,squeeze(score_reach_rs(:,3,2)),'linewidth',lw,'Color',behavior.colors.right_color); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from collection (s)'); ylabel('PC2'); 

subplot(325)
plot(time_bins,squeeze(score_init_rs(:,1,3)),'linewidth',lw,'Color',behavior.colors.push_clr); hold on
plot(time_bins,squeeze(score_init_rs(:,2,3)),'linewidth',lw,'Color',behavior.colors.pull_clr); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from init (s)'); ylabel('PC3'); 

subplot(326)
plot(time_bins,squeeze(score_reach_rs(:,1,3)),'linewidth',lw,'Color',behavior.colors.left_color); hold on
plot(time_bins,squeeze(score_reach_rs(:,2,3)),'linewidth',lw,'Color',behavior.colors.center_color); 
plot(time_bins,squeeze(score_reach_rs(:,3,3)),'linewidth',lw,'Color',behavior.colors.right_color); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from collection (s)'); ylabel('PC3');
set(gcf,'color','w')
saveas(gcf,strcat(out_folder,filesep,'trial_avg_pca.png'),'png');

figure
subplot(1,4,[1,2])
plot(cumsum(explained_init),'linewidth',lw); hold on
plot(cumsum(explained_reach),'linewidth',lw);
legend('init','reach','box','off','location','southeast')
set(gca,axeOpt{:})
xlabel('PCs'); ylabel('variance explained');

subplot(1,4,3)
imagesc(coeff_init(:,1:10))
xlabel('PCs'); ylabel('coefficients (loadings)');
title('init')

subplot(1,4,4)
imagesc(coeff_reach(:,1:10))
xlabel('PCs'); %ylabel('coefficients (loadings)');
title('reach')
set(gcf,'color','w','position',[2438 469 1108 476])

saveas(gcf,strcat(out_folder,filesep,'varianceExplained_coeffs_trial_avg.png'),'png');


%% SMOOTHED
figure
subplot(321)
plot(st_bins,squeeze(score_init_sm(:,1,1)),'linewidth',lw,'Color',behavior.colors.push_clr); hold on
plot(st_bins,squeeze(score_init_sm(:,2,1)),'linewidth',lw,'Color',behavior.colors.pull_clr); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from init (s)'); ylabel('PC1'); 

subplot(322)
plot(st_bins,squeeze(score_reach_sm(:,1,1)),'linewidth',lw,'Color',behavior.colors.left_color); hold on
plot(st_bins,squeeze(score_reach_sm(:,2,1)),'linewidth',lw,'Color',behavior.colors.center_color); 
plot(st_bins,squeeze(score_reach_sm(:,3,1)),'linewidth',lw,'Color',behavior.colors.right_color); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from collection (s)'); ylabel('PC1'); 

subplot(323)
plot(st_bins,squeeze(score_init_sm(:,1,2)),'linewidth',lw,'Color',behavior.colors.push_clr); hold on
plot(st_bins,squeeze(score_init_sm(:,2,2)),'linewidth',lw,'Color',behavior.colors.pull_clr); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from init (s)'); ylabel('PC2'); 

subplot(324)
plot(st_bins,squeeze(score_reach_sm(:,1,2)),'linewidth',lw,'Color',behavior.colors.left_color); hold on
plot(st_bins,squeeze(score_reach_sm(:,2,2)),'linewidth',lw,'Color',behavior.colors.center_color); 
plot(st_bins,squeeze(score_reach_sm(:,3,2)),'linewidth',lw,'Color',behavior.colors.right_color); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from collection (s)'); ylabel('PC2'); 

subplot(325)
plot(st_bins,squeeze(score_init_sm(:,1,3)),'linewidth',lw,'Color',behavior.colors.push_clr); hold on
plot(st_bins,squeeze(score_init_sm(:,2,3)),'linewidth',lw,'Color',behavior.colors.pull_clr); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from init (s)'); ylabel('PC3'); 

subplot(326)
plot(st_bins,squeeze(score_reach_sm(:,1,3)),'linewidth',lw,'Color',behavior.colors.left_color); hold on
plot(st_bins,squeeze(score_reach_sm(:,2,3)),'linewidth',lw,'Color',behavior.colors.center_color); 
plot(st_bins,squeeze(score_reach_sm(:,3,3)),'linewidth',lw,'Color',behavior.colors.right_color); 
xlim([-2.5 1.5]);
set(gca,axeOpt{:})
xlabel('time from collection (s)'); ylabel('PC3');
set(gcf,'color','w')
saveas(gcf,strcat(out_folder,filesep,'trial_avg_pca_smoothed.png'),'png');


%% plot 3d
figure
subplot(121)
transp = 0.3;
transp_p = 0.05;
scatter3(score_init_rs(:,1,1)',score_init_rs(:,1,2)',score_init_rs(:,1,3)',15,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp); hold on
scatter3(score_init_rs(1,1,1)',score_init_rs(1,1,2)',score_init_rs(1,1,3)',100,'o','MarkerEdgeColor',behavior.colors.push_clr,'LineWidth',2);
scatter3(score_init_rs(init_idx,1,1)',score_init_rs(init_idx,1,2)',score_init_rs(init_idx,1,3)',100,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',1);
scatter3(score_init_rs(:,2,1)',score_init_rs(:,2,2)',score_init_rs(:,2,3)',15,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp)
scatter3(score_init_rs(init_idx,2,1)',score_init_rs(init_idx,2,2)',score_init_rs(init_idx,2,3)',100,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',1);
scatter3(score_init_rs(1,2,1)',score_init_rs(1,2,2)',score_init_rs(1,2,3)',100,'o','MarkerEdgeColor',behavior.colors.pull_clr,'LineWidth',2);

extra_points = [-2 1];
yLi = get(gca,'YLim')+extra_points;
zLi = get(gca,'ZLim')+extra_points;
xLi = get(gca,'XLim')+extra_points;
oneMati = ones(size(score_init_rs(:,1,1)));

scatter3(score_init_rs(:,1,1)', oneMati .* yLi(2), score_init_rs(:,1,3)',10,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp_p);
scatter3(score_init_rs(init_idx,1,1)', oneMati(init_idx) .* yLi(2), score_init_rs(init_idx,1,3)',80,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp);
scatter3(score_init_rs(:,1,1)',  score_init_rs(:,1,2)', oneMati .* zLi(1),10,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp_p);
scatter3(score_init_rs(init_idx,1,1)',  score_init_rs(init_idx,1,2)', oneMati(init_idx) .* zLi(1),80,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp);
%scatter3(oneMat .* xL(2),  score_init_rs(:,1,2)',score_init_rs(:,1,3)',10,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp_p);
hold on;

scatter3(score_init_rs(:,2,1)', oneMati .* yLi(2), score_init_rs(:,2,3)',10,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp_p);
scatter3(score_init_rs(init_idx,2,1)', oneMati(init_idx) .* yLi(2), score_init_rs(init_idx,2,3)',80,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp);
scatter3(score_init_rs(:,2,1)',  score_init_rs(:,2,2)', oneMati .* zLi(1),10,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp_p);
scatter3(score_init_rs(init_idx,2,1)',  score_init_rs(init_idx,2,2)', oneMati(init_idx) .* zLi(1),80,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp);
%scatter3(oneMat .* xL(2),  score_init_rs(:,2,2)',score_init_rs(:,2,3)',10,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp_p);
hold off;
view(-40,20)
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
axis([xLi yLi zLi])
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1.5,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top',...
    'GridAlpha',0.05);




subplot(122)
transp = 0.3;
transp_p = 0.05;
scatter3(score_reach_rs(:,1,1)',score_reach_rs(:,1,2)',score_reach_rs(:,1,3)',15,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp); hold on
scatter3(score_reach_rs(1,1,1)',score_reach_rs(1,1,2)',score_reach_rs(1,1,3)',100,'o','MarkerEdgeColor',behavior.colors.left_color,'LineWidth',2);
scatter3(score_reach_rs(init_idx,1,1)',score_reach_rs(init_idx,1,2)',score_reach_rs(init_idx,1,3)',100,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',1);

scatter3(score_reach_rs(:,2,1)',score_reach_rs(:,2,2)',score_reach_rs(:,2,3)',15,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp)
scatter3(score_reach_rs(init_idx,2,1)',score_reach_rs(init_idx,2,2)',score_reach_rs(init_idx,2,3)',100,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',1);
scatter3(score_reach_rs(1,2,1)',score_reach_rs(1,2,2)',score_reach_rs(1,2,3)',100,'o','MarkerEdgeColor',behavior.colors.center_color,'LineWidth',2);

scatter3(score_reach_rs(:,3,1)',score_reach_rs(:,3,2)',score_reach_rs(:,3,3)',15,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp)
scatter3(score_reach_rs(init_idx,3,1)',score_reach_rs(init_idx,3,2)',score_reach_rs(init_idx,3,3)',100,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',1);
scatter3(score_reach_rs(1,3,1)',score_reach_rs(1,3,2)',score_reach_rs(1,3,3)',100,'o','MarkerEdgeColor',behavior.colors.right_color,'LineWidth',2);

extra_points = [-2 1];
yLr = get(gca,'YLim')+extra_points;
zLr = get(gca,'ZLim')+extra_points;
xLr = get(gca,'XLim')+extra_points;
oneMati = ones(size(score_reach_rs(:,1,1)));

scatter3(score_reach_rs(:,1,1)', oneMati .* yLr(2), score_reach_rs(:,1,3)',10,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_rs(init_idx,1,1)', oneMati(init_idx) .* yLr(2), score_reach_rs(init_idx,1,3)',80,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp);
scatter3(score_reach_rs(:,1,1)',  score_reach_rs(:,1,2)', oneMati .* zLr(1),10,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_rs(init_idx,1,1)',  score_reach_rs(init_idx,1,2)', oneMati(init_idx) .* zLr(1),80,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp);
%scatter3(oneMat .* xL(2),  score_reach_rs(:,1,2)',score_reach_rs(:,1,3)',10,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp_p);

scatter3(score_reach_rs(:,2,1)', oneMati .* yLr(2), score_reach_rs(:,2,3)',10,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_rs(init_idx,2,1)', oneMati(init_idx) .* yLr(2), score_reach_rs(init_idx,2,3)',80,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp);
scatter3(score_reach_rs(:,2,1)',  score_reach_rs(:,2,2)', oneMati .* zLr(1),10,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_rs(init_idx,2,1)',  score_reach_rs(init_idx,2,2)', oneMati(init_idx) .* zLr(1),80,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp);

scatter3(score_reach_rs(:,3,1)', oneMati .* yLr(2), score_reach_rs(:,3,3)',10,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_rs(init_idx,3,1)', oneMati(init_idx) .* yLr(2), score_reach_rs(init_idx,3,3)',80,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp);
scatter3(score_reach_rs(:,2,1)',  score_reach_rs(:,2,2)', oneMati .* zLr(1),10,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_rs(init_idx,3,1)',  score_reach_rs(init_idx,3,2)', oneMati(init_idx) .* zLr(1),80,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp);

%scatter3(oneMat .* xL(2),  score_reach_rs(:,2,2)',score_reach_rs(:,2,3)',10,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp_p);
hold off;
view(60, 26)
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
axis([xLr yLr zLr])
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1.5,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top',...
    'GridAlpha',0.05);
set(gcf,'Position',[2229 233 1361 687])
saveas(gcf,strcat(out_folder,filesep,'3d_scores.png'),'png');


%% plot 3d smoothed
figure
subplot(121)
transp = 0.3;
transp_p = 0.05;
scatter3(score_init_sm(:,1,1)',score_init_sm(:,1,2)',score_init_sm(:,1,3)',15,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp); hold on
scatter3(score_init_sm(1,1,1)',score_init_sm(1,1,2)',score_init_sm(1,1,3)',100,'o','MarkerEdgeColor',behavior.colors.push_clr,'LineWidth',2);
scatter3(score_init_sm(init_idx,1,1)',score_init_sm(init_idx,1,2)',score_init_sm(init_idx,1,3)',100,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',1);
scatter3(score_init_sm(:,2,1)',score_init_sm(:,2,2)',score_init_sm(:,2,3)',15,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp)
scatter3(score_init_sm(init_idx,2,1)',score_init_sm(init_idx,2,2)',score_init_sm(init_idx,2,3)',100,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',1);
scatter3(score_init_sm(1,2,1)',score_init_sm(1,2,2)',score_init_sm(1,2,3)',100,'o','MarkerEdgeColor',behavior.colors.pull_clr,'LineWidth',2);

extra_points = [-2 1];
yLi = get(gca,'YLim')+extra_points;
zLi = get(gca,'ZLim')+extra_points;
xLi = get(gca,'XLim')+extra_points;
oneMati = ones(size(score_init_sm(:,1,1)));

scatter3(score_init_sm(:,1,1)', oneMati .* yLi(2), score_init_sm(:,1,3)',10,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp_p);
scatter3(score_init_sm(init_idx,1,1)', oneMati(init_idx) .* yLi(2), score_init_sm(init_idx,1,3)',80,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp);
scatter3(score_init_sm(:,1,1)',  score_init_sm(:,1,2)', oneMati .* zLi(1),10,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp_p);
scatter3(score_init_sm(init_idx,1,1)',  score_init_sm(init_idx,1,2)', oneMati(init_idx) .* zLi(1),80,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp);
%scatter3(oneMat .* xL(2),  score_init_sm(:,1,2)',score_init_sm(:,1,3)',10,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp_p);
hold on;

scatter3(score_init_sm(:,2,1)', oneMati .* yLi(2), score_init_sm(:,2,3)',10,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp_p);
scatter3(score_init_sm(init_idx,2,1)', oneMati(init_idx) .* yLi(2), score_init_sm(init_idx,2,3)',80,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp);
scatter3(score_init_sm(:,2,1)',  score_init_sm(:,2,2)', oneMati .* zLi(1),10,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp_p);
scatter3(score_init_sm(init_idx,2,1)',  score_init_sm(init_idx,2,2)', oneMati(init_idx) .* zLi(1),80,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp);
%scatter3(oneMat .* xL(2),  score_init_sm(:,2,2)',score_init_sm(:,2,3)',10,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp_p);
hold off;
view(-40,20)
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
axis([xLi yLi zLi])
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1.5,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top',...
    'GridAlpha',0.05);




subplot(122)
transp = 0.3;
transp_p = 0.05;
scatter3(score_reach_sm(:,1,1)',score_reach_sm(:,1,2)',score_reach_sm(:,1,3)',15,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp); hold on
scatter3(score_reach_sm(1,1,1)',score_reach_sm(1,1,2)',score_reach_sm(1,1,3)',100,'o','MarkerEdgeColor',behavior.colors.left_color,'LineWidth',2);
scatter3(score_reach_sm(init_idx,1,1)',score_reach_sm(init_idx,1,2)',score_reach_sm(init_idx,1,3)',100,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',1);

scatter3(score_reach_sm(:,2,1)',score_reach_sm(:,2,2)',score_reach_sm(:,2,3)',15,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp)
scatter3(score_reach_sm(init_idx,2,1)',score_reach_sm(init_idx,2,2)',score_reach_sm(init_idx,2,3)',100,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',1);
scatter3(score_reach_sm(1,2,1)',score_reach_sm(1,2,2)',score_reach_sm(1,2,3)',100,'o','MarkerEdgeColor',behavior.colors.center_color,'LineWidth',2);

scatter3(score_reach_sm(:,3,1)',score_reach_sm(:,3,2)',score_reach_sm(:,3,3)',15,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp)
scatter3(score_reach_sm(init_idx,3,1)',score_reach_sm(init_idx,3,2)',score_reach_sm(init_idx,3,3)',100,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',1);
scatter3(score_reach_sm(1,3,1)',score_reach_sm(1,3,2)',score_reach_sm(1,3,3)',100,'o','MarkerEdgeColor',behavior.colors.right_color,'LineWidth',2);

extra_points = [-2 1];
yLr = get(gca,'YLim')+extra_points;
zLr = get(gca,'ZLim')+extra_points;
xLr = get(gca,'XLim')+extra_points;
oneMati = ones(size(score_reach_sm(:,1,1)));

scatter3(score_reach_sm(:,1,1)', oneMati .* yLr(2), score_reach_sm(:,1,3)',10,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_sm(init_idx,1,1)', oneMati(init_idx) .* yLr(2), score_reach_sm(init_idx,1,3)',80,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp);
scatter3(score_reach_sm(:,1,1)',  score_reach_sm(:,1,2)', oneMati .* zLr(1),10,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_sm(init_idx,1,1)',  score_reach_sm(init_idx,1,2)', oneMati(init_idx) .* zLr(1),80,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp);
%scatter3(oneMat .* xL(2),  score_reach_sm(:,1,2)',score_reach_sm(:,1,3)',10,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp_p);

scatter3(score_reach_sm(:,2,1)', oneMati .* yLr(2), score_reach_sm(:,2,3)',10,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_sm(init_idx,2,1)', oneMati(init_idx) .* yLr(2), score_reach_sm(init_idx,2,3)',80,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp);
scatter3(score_reach_sm(:,2,1)',  score_reach_sm(:,2,2)', oneMati .* zLr(1),10,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_sm(init_idx,2,1)',  score_reach_sm(init_idx,2,2)', oneMati(init_idx) .* zLr(1),80,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp);

scatter3(score_reach_sm(:,3,1)', oneMati .* yLr(2), score_reach_sm(:,3,3)',10,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_sm(init_idx,3,1)', oneMati(init_idx) .* yLr(2), score_reach_sm(init_idx,3,3)',80,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp);
scatter3(score_reach_sm(:,2,1)',  score_reach_sm(:,2,2)', oneMati .* zLr(1),10,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp_p);
scatter3(score_reach_sm(init_idx,3,1)',  score_reach_sm(init_idx,3,2)', oneMati(init_idx) .* zLr(1),80,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp);

%scatter3(oneMat .* xL(2),  score_reach_sm(:,2,2)',score_reach_sm(:,2,3)',10,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp_p);
hold off;
view(60, 26)
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
axis([xLr yLr zLr])
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1.5,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top',...
    'GridAlpha',0.05);
set(gcf,'Position',[2229 233 1361 687])
saveas(gcf,strcat(out_folder,filesep,'3d_scores_smooth.png'),'png');

%%
figure(10)
fps = 0.1;
transp = 0.3;
scatter3(score_init_sm(:,1,1)',score_init_sm(:,1,2)',score_init_sm(:,1,3)',15,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',transp);
hold on
scatter3(score_init_sm(1,1,1)',score_init_sm(1,1,2)',score_init_sm(1,1,3)',100,'o','MarkerEdgeColor',behavior.colors.push_clr,'LineWidth',2);
scatter3(score_init_sm(init_idx,1,1)',score_init_sm(init_idx,1,2)',score_init_sm(init_idx,1,3)',100,'o','filled','MarkerFaceColor',behavior.colors.push_clr,'MarkerFaceAlpha',1);
scatter3(score_init_sm(:,2,1)',score_init_sm(:,2,2)',score_init_sm(:,2,3)',15,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',transp)
scatter3(score_init_sm(init_idx,2,1)',score_init_sm(init_idx,2,2)',score_init_sm(init_idx,2,3)',100,'o','filled','MarkerFaceColor',behavior.colors.pull_clr,'MarkerFaceAlpha',1);
scatter3(score_init_sm(1,2,1)',score_init_sm(1,2,2)',score_init_sm(1,2,3)',100,'o','MarkerEdgeColor',behavior.colors.pull_clr,'LineWidth',2);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1.5,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top');
set(gcf,'Position',[2550 188 725 575],'color','w'); grid off
shg
axis([xLi-extra_points yLi-extra_points zLi-extra_points])

%rotate_3Dplot(fps)
rotate_and_save_3Dplot(fps, strcat(out_folder,filesep,'rotatingScores_init_smooth.gif'))

%%
figure(11)
fps = 0.1;
transp = 0.3;
scatter3(score_reach_sm(:,1,1)',score_reach_sm(:,1,2)',score_reach_sm(:,1,3)',15,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',transp); hold on
scatter3(score_reach_sm(1,1,1)',score_reach_sm(1,1,2)',score_reach_sm(1,1,3)',100,'o','MarkerEdgeColor',behavior.colors.left_color,'LineWidth',2);
scatter3(score_reach_sm(init_idx,1,1)',score_reach_sm(init_idx,1,2)',score_reach_sm(init_idx,1,3)',100,'o','filled','MarkerFaceColor',behavior.colors.left_color,'MarkerFaceAlpha',1);
scatter3(score_reach_sm(:,2,1)',score_reach_sm(:,2,2)',score_reach_sm(:,2,3)',15,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',transp)
scatter3(score_reach_sm(init_idx,2,1)',score_reach_sm(init_idx,2,2)',score_reach_sm(init_idx,2,3)',100,'o','filled','MarkerFaceColor',behavior.colors.center_color,'MarkerFaceAlpha',1);
scatter3(score_reach_sm(1,2,1)',score_reach_sm(1,2,2)',score_reach_sm(1,2,3)',100,'o','MarkerEdgeColor',behavior.colors.center_color,'LineWidth',2);
scatter3(score_reach_sm(:,3,1)',score_reach_sm(:,3,2)',score_reach_sm(:,3,3)',15,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',transp)
scatter3(score_reach_sm(init_idx,3,1)',score_reach_sm(init_idx,3,2)',score_reach_sm(init_idx,3,3)',100,'o','filled','MarkerFaceColor',behavior.colors.right_color,'MarkerFaceAlpha',1);
scatter3(score_reach_sm(1,3,1)',score_reach_sm(1,3,2)',score_reach_sm(1,3,3)',100,'o','MarkerEdgeColor',behavior.colors.right_color,'LineWidth',2);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1.5,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top');
set(gcf,'Position',[2550 188 725 575],'color','w'); grid off
shg
axis([xLr-extra_points yLr-extra_points zLr-extra_points])

%rotate_3Dplot(fps)
rotate_and_save_3Dplot(fps, strcat(out_folder,filesep,'rotatingScores_reach_smooth.gif'))

%% Save
save(strcat(imec_path,filesep,'neurons_session.mat'),'trial_avg_pca',...
    '-append');


