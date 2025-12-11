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

for im = 0:1
    imec = num2str(im);
    fprintf('%s%i','running imec ',im)

    % Load behavior
    behavior_path = strcat(behavior_root,filesep,group_setup_behav,filesep,mouse,filesep,session);
    load(strcat(behavior_path,filesep,'behavior_session.mat'));

    % Load ehpys
    imec_path = strcat(ephys_root,filesep,ephys_sess,filesep,ephys_sess,'_imec',imec,filesep,'catGT',filesep,'kilosort4');
    neurons_imec = load(strcat(imec_path,filesep,'neurons_session.mat'));

    % Folder to save
    out_folder = strcat(imec_path,filesep,'sess_out_figures',filesep,'svm_win100_slidingWin');
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

    %% Do a sliding window of for bin counts

    tm_before = neurons_imec.neurons_params.events.tm_before + neurons_imec.neurons_params.gamma.k.paddx(1);
    tm_after = neurons_imec.neurons_params.events.tm_after + neurons_imec.neurons_params.gamma.k.paddx(2); % sec
    win_dur = tm_after-tm_before;
    bin_width = 0.1;
    bin_step = 0.01;
    n_bins = floor((win_dur/bin_step)-(bin_width/bin_step))+1;

    spk_count_init = nan(n_bins, nr_trials, nr_eg_neurons);
    spk_count_reach = nan(n_bins, nr_trials, nr_eg_neurons);
    bin_edges = nan(n_bins,2);

    for n = 1:nr_eg_neurons
        sts_init = neu_eg.st_init{n,1};
        sts_reach = neu_eg.st_reach{n,1};

        for tt = 1:nr_trials
            for i = 1:n_bins
                bin_start = tm_before + (i-1)*bin_step;
                bin_stop = tm_before + bin_width + (i-1)*bin_step;

                spk_count_init(i,tt,n) = numel(find(sts_init{tt}>bin_start & sts_init{tt}<bin_stop));
                spk_count_reach(i,tt,n) = numel(find(sts_reach{tt}>bin_start & sts_reach{tt}<bin_stop));

                bin_edges(i,1) = bin_start;
                bin_edges(i,2) = bin_stop;
            end
        end
    end

    %% %% SVM DECODER OF PUSH/PULL

    % Parameters
    ratio_train_val = 0.8;
    ncv = 100;
    nfold = 10;
    Cvec = [0.001, 0.01, 0.1, 1, 10, 100, 1000];


    bac_pp = nan(n_bins,2);
    tic
    % Run by bins
    for b = 1:n_bins
        % spike trains by bins
        s1_pp_init = squeeze(spk_count_init(b,push_idx,:));
        s2_pp_init = squeeze(spk_count_init(b,pull_idx,:));
        s1_pp_reach = squeeze(spk_count_reach(b,push_idx,:));
        s2_pp_reach = squeeze(spk_count_reach(b,pull_idx,:));

        % svm accuracy
        [bac_pp(b,1)] = svm_simple_fun_adaptedfromVK(s1_pp_init,s2_pp_init,ratio_train_val,ncv,nfold,Cvec);
        [bac_pp(b,2)] = svm_simple_fun_adaptedfromVK(s1_pp_reach,s2_pp_reach,ratio_train_val,ncv,nfold,Cvec);

        fprintf('%s%i%s','bin ',b,' done... ')
        toc
    end
    fprintf('All done!! ')
    toc

end

%% Save


svm_decoder.bin_size.wins.bin_width = bin_width;
svm_decoder.bin_size.wins.bin_step = bin_step;
svm_decoder.bin_size.wins.n_bins = n_bins;
svm_decoder.bin_size.wins.bin_edges = bin_edges;
svm_decoder.bin_size.wins.spk_count_init = spk_count_init;
svm_decoder.bin_size.wins.spk_count_reach = spk_count_reach;


svm_decoder.svm_params.ratio_train_val = ratio_train_val;
svm_decoder.svm_params.ncv = ncv;
svm_decoder.svm_params.nfold = nfold;
svm_decoder.svm_params.Cvec = Cvec;
svm_decoder.bac_pp = bac_pp;

save(strcat(imec_path,filesep,'svm_decoder_PP_slindingWin.mat'),'svm_decoder');

%% Plot accuracy!

clrss = cat(1,behavior.colors.push_clr,behavior.colors.pull_clr);
clr_pp = mean(clrss,1);
axeOpt=neurons_imec.figProp.axeOpt;
time_bins = bin_edges(:,1);

figure()
subplot(121)
% plot(time_bins,bac_right(:,1),'LineWidth',2,'Color',behavior.colors.right_color); hold on
% plot(time_bins,bac_left(:,1),'LineWidth',2,'Color',behavior.colors.left_color); 
% plot(time_bins,bac_center(:,1),'LineWidth',2,'Color',behavior.colors.center_color); 
plot(time_bins,bac_pp(:,1),'LineWidth',2,'Color',clr_pp); 

yline(0.5,'linewidth',3,'Alpha',0.7,'Color',[.9 .9 .9],'LineStyle','--');
xline(0,'linewidth',1.5,'Alpha',0.5,'Color',[.9 .9 .9],'LineStyle','-.');
axis([tm_before tm_after .48 1])
set(gca,axeOpt{:})
xlabel('time from trial initiation (s)'); ylabel('accuracy'); 

subplot(122)
% plot(time_bins,bac_right(:,2),'LineWidth',2,'Color',behavior.colors.right_color); hold on
% plot(time_bins,bac_left(:,2),'LineWidth',2,'Color',behavior.colors.left_color); 
% plot(time_bins,bac_center(:,2),'LineWidth',2,'Color',behavior.colors.center_color); 
plot(time_bins,bac_pp(:,2),'LineWidth',2,'Color',clr_pp); 

yline(0.5,'linewidth',3,'Alpha',0.7,'Color',[.9 .9 .9],'LineStyle','--');
xline(0,'linewidth',1.5,'Alpha',0.5,'Color',[.9 .9 .9],'LineStyle','-.');
axis([tm_before tm_after .48 1])
set(gca,axeOpt{:})
xlabel('time from water collection (s)'); ylabel('accuracy'); 
%legend('right','left','center','push / pull','box','off')
legend('push / pull','box','off')

set(gcf,'Position', [2266 474 1313 470], 'Color','w')
saveas(gcf,strcat(out_folder,filesep,'svm_decoder_pp_slindingWin.png'),'png');





