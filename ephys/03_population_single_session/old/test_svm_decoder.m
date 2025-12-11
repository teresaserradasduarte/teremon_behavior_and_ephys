%% Linear SVM with Monte-Carlo cross-valitaton: TESTING
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

% Folder to save
out_folder = strcat(imec_path,filesep,'sess_out_figures',filesep,'svm_win200_pp_L_R_C');
if ~exist('out_folder','dir'), mkdir(out_folder); end

% figs to plot
plot_rasters_flag = 0;

%% Load neuron table
neu = struct2table(neurons_imec.neurons);
eg_idx = find(neu.extraGood == 1);
neu_eg = neu(eg_idx,:);

% Dims
%nr_trials = behavior.behavior_duration.trial_end;
nr_trials =size(neu_eg.st_init{1,1},1);
nr_eg_neurons = length(eg_idx);

% Condition trials - behavior
push_idx = behavior.init.idx_trial_push(behavior.init.idx_trial_push<=nr_trials);
pull_idx = behavior.init.idx_trial_pull(behavior.init.idx_trial_pull<=nr_trials);
left_idx = behavior.reach.left_idx(behavior.reach.left_idx<=nr_trials);
right_idx = behavior.reach.right_idx(behavior.reach.right_idx<=nr_trials);
center_idx = behavior.reach.center_idx(behavior.reach.center_idx<=nr_trials);


%% SPK COUNTS for decoder!
tm_before = neurons_imec.neurons_params.events.tm_before + neurons_imec.neurons_params.gamma.k.paddx(1);
tm_after = neurons_imec.neurons_params.events.tm_after + neurons_imec.neurons_params.gamma.k.paddx(2); % sec
bin_width_d = 0.2;
bin_edges_d = tm_before:bin_width_d:tm_after;
n_bins = length(bin_edges_d)-1;

bin_counts_init = nan(length(bin_edges_d)-1,nr_trials,nr_eg_neurons);
bin_counts_reach = nan(length(bin_edges_d)-1,nr_trials,nr_eg_neurons);

for n = 1:nr_eg_neurons
    sts_init = neu_eg.st_init{n,1};
    sts_reach = neu_eg.st_reach{n,1};

    for tt = 1:nr_trials
        bin_counts_init(:,tt,n) = histcounts(sts_init{tt},bin_edges_d);
        bin_counts_reach(:,tt,n) = histcounts(sts_reach{tt},bin_edges_d);
    end

    if plot_rasters_flag
        figure
        subplot(121)
        imagesc(flipud(squeeze(bin_counts_init(:,:,n)'))); colormap(flipud(gray));
        ax = gca; xticklabels(ax,bin_edges_d(ax.XTick));
        axis square
        xlabel('')
        title('spike counts arount init')
        subplot(122)
        imagesc(flipud(squeeze(bin_counts_reach(:,:,n)'))); colormap(flipud(gray));
        ax = gca; xticklabels(ax,bin_edges_d(ax.XTick));
        axis square
        title('spike counts arount reach')
        set(gcf,'position',[2454 259 1226 677])
    end
end




% --------------------------------------------------------------------
%% SVM DECODER OF PUSH/PULL

% Parameters
ratio_train_val = 0.8;
ncv = 100;
nfold = 10;
Cvec = [0.001, 0.01, 0.1, 1, 10, 100, 1000];


bac_pp = nan(n_bins,2);
bac_right = nan(n_bins,2);
bac_left = nan(n_bins,2);
bac_center = nan(n_bins,2);

tic
% Run by bins
for b = 1:n_bins
    % spike trains by bins
    
    % INIT
    % s1 + s2 push/pull init
    s1_pp_init = squeeze(bin_counts_init(b,push_idx,:));
    s2_pp_init = squeeze(bin_counts_init(b,pull_idx,:));
    % s1 + s2 right init
    s1_right_init = squeeze(bin_counts_init(b,right_idx,:));
    s2_right_init = squeeze(bin_counts_init(b,cat(2,left_idx,center_idx),:));
    % s1 + s2 left init
    s1_left_init = squeeze(bin_counts_init(b,left_idx,:));
    s2_left_init = squeeze(bin_counts_init(b,cat(2,right_idx,center_idx),:));
    % s1 + s2 center init
    s1_center_init = squeeze(bin_counts_init(b,center_idx,:));
    s2_center_init = squeeze(bin_counts_init(b,cat(2,right_idx,left_idx),:));


    % s1 + s2 push/pull reach
    s1_pp_reach = squeeze(bin_counts_reach(b,push_idx,:));
    s2_pp_reach = squeeze(bin_counts_reach(b,pull_idx,:));
    % s1 + s2 right init
    s1_right_reach = squeeze(bin_counts_reach(b,right_idx,:));
    s2_right_reach = squeeze(bin_counts_reach(b,cat(2,left_idx,center_idx),:));
    % s1 + s2 left init
    s1_left_reach = squeeze(bin_counts_reach(b,left_idx,:));
    s2_left_reach = squeeze(bin_counts_reach(b,cat(2,right_idx,center_idx),:));
    % s1 + s2 center init
    s1_center_reach = squeeze(bin_counts_reach(b,center_idx,:));
    s2_center_reach = squeeze(bin_counts_reach(b,cat(2,right_idx,left_idx),:));

    % svm accuracy
    % pp
    [bac_pp(b,1)] = svm_simple_fun_adaptedfromVK(s1_pp_init,s2_pp_init,ratio_train_val,ncv,nfold,Cvec);
    [bac_pp(b,2)] = svm_simple_fun_adaptedfromVK(s1_pp_reach,s2_pp_reach,ratio_train_val,ncv,nfold,Cvec);
    % right
    [bac_right(b,1)] = svm_simple_fun_adaptedfromVK(s1_right_init,s2_right_init,ratio_train_val,ncv,nfold,Cvec);
    [bac_right(b,2)] = svm_simple_fun_adaptedfromVK(s1_right_reach,s2_right_reach,ratio_train_val,ncv,nfold,Cvec);    
    % left
    [bac_left(b,1)] = svm_simple_fun_adaptedfromVK(s1_left_init,s2_left_init,ratio_train_val,ncv,nfold,Cvec);
    [bac_left(b,2)] = svm_simple_fun_adaptedfromVK(s1_left_reach,s2_left_reach,ratio_train_val,ncv,nfold,Cvec);
    % center
    [bac_center(b,1)] = svm_simple_fun_adaptedfromVK(s1_center_init,s2_center_init,ratio_train_val,ncv,nfold,Cvec);
    [bac_center(b,2)] = svm_simple_fun_adaptedfromVK(s1_center_reach,s2_center_reach,ratio_train_val,ncv,nfold,Cvec);

    fprintf('%s%i%s','bin ',b,' done... ')
    toc
end
fprintf('All done!! ')
toc


%% Save
svm_decoder.bin_size.bin_width_d = bin_width_d;
svm_decoder.bin_size.bin_edges_d = bin_edges_d;
svm_decoder.bin_size.n_bins = n_bins;
svm_decoder.bin_counts_init = bin_counts_init;
svm_decoder.bin_counts_reach = bin_counts_reach;
svm_decoder.svm_params.ratio_train_val = ratio_train_val;
svm_decoder.svm_params.ncv = ncv;
svm_decoder.svm_params.nfold = nfold;
svm_decoder.svm_params.Cvec = Cvec;
svm_decoder.bac_pp = bac_pp;
svm_decoder.bac_right = bac_right;
svm_decoder.bac_left = bac_left;
svm_decoder.bac_center = bac_center;

save(strcat(imec_path,filesep,'svm_decoder_PP_R_L_C.mat'),'svm_decoder');

%% Plot accuracy!

clrss = cat(1,behavior.colors.push_clr,behavior.colors.pull_clr);
clr_pp = mean(clrss,1);
axeOpt=neurons_imec.figProp.axeOpt;
time_bins = bin_edges_d(1:end-1);

figure()
subplot(121)
plot(time_bins,bac_right(:,1),'LineWidth',2,'Color',behavior.colors.right_color); hold on
plot(time_bins,bac_left(:,1),'LineWidth',2,'Color',behavior.colors.left_color); 
plot(time_bins,bac_center(:,1),'LineWidth',2,'Color',behavior.colors.center_color); 
plot(time_bins,bac_pp(:,1),'LineWidth',2,'Color',clr_pp); 

yline(0.5,'linewidth',3,'Alpha',0.7,'Color',[.9 .9 .9],'LineStyle','--');
xline(0,'linewidth',1.5,'Alpha',0.5,'Color',[.9 .9 .9],'LineStyle','-.');
axis([tm_before tm_after .48 1])
set(gca,axeOpt{:})
xlabel('time from trial initiation (s)'); ylabel('accuracy'); 

subplot(122)
plot(time_bins,bac_right(:,2),'LineWidth',2,'Color',behavior.colors.right_color); hold on
plot(time_bins,bac_left(:,2),'LineWidth',2,'Color',behavior.colors.left_color); 
plot(time_bins,bac_center(:,2),'LineWidth',2,'Color',behavior.colors.center_color); 
plot(time_bins,bac_pp(:,2),'LineWidth',2,'Color',clr_pp); 

yline(0.5,'linewidth',3,'Alpha',0.7,'Color',[.9 .9 .9],'LineStyle','--');
xline(0,'linewidth',1.5,'Alpha',0.5,'Color',[.9 .9 .9],'LineStyle','-.');
axis([tm_before tm_after .48 1])
set(gca,axeOpt{:})
xlabel('time from water collection (s)'); ylabel('accuracy'); 
legend('right','left','center','push / pull','box','off')

set(gcf,'Position', [2266 474 1313 470], 'Color','w')
saveas(gcf,strcat(out_folder,filesep,'svm_decoder_pp.png'),'png');






