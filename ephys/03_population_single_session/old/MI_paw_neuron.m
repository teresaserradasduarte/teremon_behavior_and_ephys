clear; close all; clc

%% Load behavior
behavior_root ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
reaching_root = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\mat_files';
ephys_root = 'E:';
save_root = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_and_behavior\20230801_ChocolateGroup';
group_setup_behav = strcat('20230511_ChocolateGroup',filesep,'headfixed_dynamicTarget');
mouse = '4_Milka';
sess = 'R4';
imec_id = 0;
ephys_sess = '18082023_Milka_StrCer_S4_g0';
sorter_folder = 'catGT\kilosort4';

% Load behavior
behavior_path = strcat(behavior_root,filesep,group_setup_behav,filesep,mouse,filesep,sess);
load(strcat(behavior_path,filesep,'behavior_session.mat'));

% Load reaching
behavior_path = strcat(reaching_root,filesep,group_setup_behav,filesep,mouse,filesep,sess);
load(strcat(behavior_path,filesep,'session_reaching_data_paw.mat'));

% Load neurons
path_neurons  = strcat(ephys_root,filesep,ephys_sess,filesep,ephys_sess,'_imec',num2str(imec_id),filesep,sorter_folder);
load(strcat(path_neurons,filesep,'neurons_session.mat'));


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
save_path = strcat(save_root,filesep,mouse,filesep,sess,filesep,'reaches_cat1',filesep,region);
if ~exist("save_path","dir"), mkdir(save_path); end


%% Extracted reaches --> paw position and timestamps

% Window of interest around event (reaches cat1)
interval_reach = [-3 2];
pad = [-1 1];
interval_reach_wPad = interval_reach+pad;
fps_vid = session.video.frame_rate;
frames_time = interval_reach_wPad(1):1/fps_vid:interval_reach_wPad(2);
frames_before = abs(interval_reach_wPad(1))*fps_vid;
frames_after = interval_reach_wPad(2)*fps_vid;
total_nr_frames = frames_before+frames_after+1;

% Reaches only - zoomed
reaches1_all = reaches.reach_mat(:,:,reaches.cat_reach==1);
reaches1_timestamps_all = reaches.reach_timestamps_mat(reaches.cat_reach==1,:);
reaches1_peak_all = reaches1_timestamps_all(:,reaches.reach_params.reach_interval.max_reach);

% Select only the window of good recordings - only reaches within window
behav_start = behavior.behavior_duration.time_start;
behav_stop = behavior.behavior_duration.time_end + behav_start;
reaches_in_interval_flag = ...
    reaches1_peak_all + interval_reach_wPad(1) > behav_start & ...
    reaches1_peak_all + interval_reach_wPad(2) < behav_stop;
reaches1 = reaches1_all(:,:,reaches_in_interval_flag);
reaches1_timestamps = reaches1_timestamps_all(reaches_in_interval_flag,:);
reaches1_peak = reaches1_timestamps(:,reaches.reach_params.reach_interval.max_reach);
n_reaches1 = size(reaches1,3);

% PCA paw - all session
paw_dom_all = session.ppDom_paw_sess;
[zs_paw_all, paw_all_mu, paw_all_sigma] = zscore(paw_dom_all);
[coeff_paw_all, score_paw_all, latent_paw_all, tsquare_paw_all, explained_paw_all, mus_paw_all] = pca(zs_paw_all);
reaches1_full = nan(total_nr_frames,3,n_reaches1);
reaches1_pcs = nan(total_nr_frames,3,n_reaches1);

% Paw position around reach - full window of interest
siz=trials.size_trial;
frames_timestamps_in_row = [];
for i=1:length(siz)
    frames_timestamps_in_row = [frames_timestamps_in_row,trials.sync.frames_timestamps(i,1:siz(i))];
end

for i=1:n_reaches1
    pk_idx = find(frames_timestamps_in_row==reaches1_peak(i));
    reaches1_full(:,:,i) = paw_dom_all(pk_idx-frames_before:pk_idx+frames_after,:);
    reaches1_pcs(:,:,i) = score_paw_all(pk_idx-frames_before:pk_idx+frames_after,:);
end
reaches_paw_pc1 = squeeze(reaches1_pcs(:,1,:));

% Smooth individual reaches
sig_pk = 0.02;
paw_bin_width = 1/fps_vid;
k_gaus_paw = gausskernel('sig',sig_pk,'binwidth',paw_bin_width);
sm_frames_time = frames_time(1)-k_gaus_paw.paddx(1):paw_bin_width:frames_time(end)-k_gaus_paw.paddx(2)+paw_bin_width;
nr_sm_frames = length(sm_frames_time);
paw_pc1_smooth = conv2(k_gaus_paw.pdf,1,reaches_paw_pc1,"valid");



%% Plot
% Figure stuff
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10};
event_size = 5;
spk_size = 6;
rcl_clr = [44 123 182]./256;
trans = .05;

% Path to save figs
save_paw_reaches = strcat(save_path,filesep,'reaches_fig');
if ~exist(save_paw_reaches,'dir'), mkdir(save_paw_reaches); end

figure
ff = tiledlayout(3,2);
for i=1:3
    nexttile
    plot(frames_time,squeeze(reaches1_full(:,i,:)),'color',cat(2,rcl_clr,trans),'LineWidth',.5); hold on
    plot(frames_time,median(reaches1_full(:,i,:),3,'omitnan'),'color',rcl_clr,'LineWidth',2);
    hold on
    plot(reaches.tm_w,median(reaches1(:,i,:),3),'color','k','LineWidth',2);
    set(gca,axeOpt{:})
    xline(0,'--','color',[.8 .8 .8])
    xlabel('time from reack peak (s)'); ylabel(strcat('dim',num2str(i)))

    nexttile
    plot(frames_time,squeeze(reaches1_pcs(:,i,:)),'color',cat(2,rcl_clr,trans),'LineWidth',.5); hold on
    plot(frames_time,median(reaches1_pcs(:,i,:),3,'omitnan'),'color',rcl_clr,'LineWidth',2);
    set(gca,axeOpt{:})
    xline(0,'--','color',[.8 .8 .8])
    xlabel('time from reack peak (s)'); ylabel(strcat('PC',num2str(i)))
    title(sprintf('%s%.1f%s','Explained variance = ',explained_paw_all(i),'%'))
end

set(gcf,'Position',[1994 58 1743 939],'Color','w')
saveas(gcf,strcat(save_paw_reaches,filesep,'reaches_cat1_win.png'),'png');


tt=100;
figure
subplot(221)
plot(frames_time,reaches_paw_pc1(:,tt)); hold on
plot(sm_frames_time,paw_pc1_smooth(:,tt),'k'); hold off
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8])
xlabel('time from reack peak (s)'); ylabel('paw PC1');
title('example trial')

subplot(222)
plot(frames_time,mean(reaches_paw_pc1,2,'omitnan'),'LineWidth',2); hold on
plot(sm_frames_time,mean(paw_pc1_smooth,2,'omitnan'),'k','linewidth',1.5); hold on
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8])
xlabel('time from reack peak (s)'); ylabel('paw PC1');
title('mean reaches')
legend('raw','smooth','box','off')

subplot(223)
plot(frames_time,reaches_paw_pc1,'color',cat(2,rcl_clr,trans)); hold on
plot(frames_time,mean(reaches_paw_pc1,2,'omitnan'),'LineWidth',2,'color',rcl_clr); hold on
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8])
xlabel('time from reack peak (s)'); ylabel('paw PC1');
title('raw all')

subplot(224)
plot(sm_frames_time,paw_pc1_smooth,'color',cat(2,[0 0 0],trans)); hold on
plot(sm_frames_time,mean(paw_pc1_smooth,2,'omitnan'),'LineWidth',2,'color',[0 0 0]); hold off
set(gca,axeOpt{:})
xline(0,'--','color',[.8 .8 .8])
xlabel('time from reack peak (s)'); ylabel('paw PC1');
title('smooth all')

set(gcf,'Position',[1991         168        1733         762],'Color','w')
saveas(gcf,strcat(save_paw_reaches,filesep,'reaches_smooth.png'),'png');



% ---------------------------------------------
%% NEURONS
% Load neuron of interest
neu = struct2table(neurons);
eg_idx = find(neu.extraGood == 1);
n_eg = length(eg_idx);

neurons_EG_reaching = neurons(eg_idx);

% Parameters for ISI calculation
bin_width = 0.002;
bin_edges = interval_reach_wPad(1):bin_width:interval_reach_wPad(2);
n_bins = length(bin_edges);

reaches1_FR_params.interval_reach_wPad = interval_reach_wPad;
reaches1_FR_params.bin_width = bin_width;
reaches1_FR_params.bin_edges = bin_edges;

reaches1_FR = nan(n_bins,n_reaches1,n_eg);
reaches1_isi = nan(n_bins,n_reaches1,n_eg);

% Smoothing individual trials
%sig_pk = 0.02;
k_gaus = gausskernel('sig',sig_pk,'binwidth',bin_width);
sm_bin_edges = bin_edges(1)-k_gaus.paddx(1):bin_width:bin_edges(end)-k_gaus.paddx(2)+bin_width;
nr_sm_bins = length(sm_bin_edges);
reaches1_FR_smooth = nan(nr_sm_bins,n_reaches1,n_eg);

% Path to save figs
save_neu_fig_path = strcat(save_path,filesep,'neurons_fig');
if ~exist(save_neu_fig_path,'dir'), mkdir(save_neu_fig_path); end

% Loop through extra-good units
for i = 1:length(eg_idx)

    un = eg_idx(i);
    neurons_EG_reaching(i).eg_idx = i;
    neurons_EG_reaching(i).st_reaches1 = cell(n_reaches1,1);
    neurons_EG_reaching(i).spk_trials_reaches1 = cell(n_reaches1,1);

    for rr = 1:n_reaches1
        spike_reach_flags = ...
            neurons(un).st > reaches1_peak(rr) + interval_reach_wPad(1) & ...
            neurons(un).st <= reaches1_timestamps(rr,end) + interval_reach_wPad(end);
        neurons_EG_reaching(i).st_reaches1{rr} = ...
            neurons(un).st(spike_reach_flags) - reaches1_peak(rr);

        neurons_EG_reaching(i).spk_trials_reaches1{rr} = ones(size(neurons_EG_reaching(i).st_reaches1{rr}))*rr;

        reach_spikes = neurons_EG_reaching(i).st_reaches1{rr};

        if ~isempty(reach_spikes)
            % Find the spike before interval fo initial ISI
            first_spk_idx=find(spike_reach_flags,1,'first');

            if first_spk_idx == 1
                idx_in_st_all = find(neurons(un).st_all==neurons(un).st(1));
                start_instant_isi = neurons(un).st_all(idx_in_st_all) - ...
                    neurons(un).st_all(idx_in_st_all-1);
            else
                start_instant_isi = neurons(un).st(first_spk_idx) - ...
                    neurons(un).st(first_spk_idx-1);
            end
            reaches1_isi(bin_edges<=reach_spikes(1),rr,i) = start_instant_isi;
            reaches1_FR(bin_edges<=reach_spikes(1),rr,i) = 1/start_instant_isi;

            % Find spikes within the interval
            for sp=2:length(reach_spikes)
                with_isi_flag = (bin_edges > reach_spikes(sp-1)) & ...
                    bin_edges <= reach_spikes(sp);
                inst_isi = reach_spikes(sp)-reach_spikes(sp-1);

                reaches1_isi(with_isi_flag,rr,i) = inst_isi;
                reaches1_FR(with_isi_flag,rr,i) = 1/inst_isi; 
            end

            % find the spike after the interval for final ISI
             last_spk_idx=find(spike_reach_flags,1,'last');
            finish_instant_isi = neurons(un).st(last_spk_idx+1) - ...
                neurons(un).st(last_spk_idx);
            reaches1_isi(bin_edges>reach_spikes(end),rr,i) = finish_instant_isi;
            reaches1_FR(bin_edges>reach_spikes(end),rr,i) = 1/finish_instant_isi;
        else
            fprintf('%s%i%s%i\n','no spikes in reach ',rr,', neuron idxEG ',i)
        end
    end
    
    % smmothing individual trials
    reaches1_FR_smooth(:,:,i) = conv2(k_gaus.pdf,1,reaches1_FR(:,:,i),"valid");

    % store in the eg neurons struct
    neurons_EG_reaching(i).FR_reaches1 = reaches1_FR(:,:,i);
    neurons_EG_reaching(i).FR_reaches1_smooth = reaches1_FR_smooth(:,:,i);


    % Show neuron

    figure
    tt = tiledlayout(3,1);

    title(tt,sprintf('%s%s%s%s%s%s%s%s%s',...
        'Reaches (cat1) ',...
        ' | region: ',neurons(un).meta.region,...
        ' | phyID: ', num2str(neurons(un).phyID),...
        ' | idx: ', num2str(neurons(un).unitIdx),...
        ' | idxEG: ', num2str(i)),...
        'Interpreter', 'none','fontsize',12,'fontweight','normal')


    axx1 = nexttile;
    plot(sm_bin_edges,mean(reaches1_FR_smooth(:,:,i),2,'omitnan'),'color',rcl_clr,LineWidth=2)
    hold off
    xlim(interval_reach)
    %axis([tm_before tm_after 0 nr_trials])
    xlabel('time from reack peak (s)'); ylabel('firing rate sp/s')
    set(gca,axeOpt{:})

    axx2 = nexttile;
    axx2.Layout.TileSpan = [2,1];
    reaches_vec = 1:n_reaches1;
    scatter(cell2mat(neurons_EG_reaching(i).st_reaches1),cell2mat(neurons_EG_reaching(i).spk_trials_reaches1),spk_size,'k.'), hold on
    scatter(zeros(size(reaches_vec)),reaches_vec,event_size,'filled','MarkerFaceColor',rcl_clr)
    hold off
    axis tight
    %axis([tm_before tm_after 0 nr_trials])
    xlabel('time from reack peak (s)'); ylabel('reach idx across session')
    set(gca,axeOpt{:})
    xlim(interval_reach)
    set(gcf,'Position',[2157 205 539 697],'Color','w')
    shg
    saveas(gcf,strcat(save_neu_fig_path,filesep,'neuron',num2str(neurons(un).phyID),'.png'),'png');

end

%% Smooth individual trial

save_example_cell_trial = strcat(save_path,filesep,'example_trials');
if ~exist(save_example_cell_trial,"dir"), mkdir(save_example_cell_trial); end

for un = 1:n_eg
    tt=100;
    figure
    subplot(321)
    plot(bin_edges,reaches1_isi(:,tt,un));hold on
    scatter(neurons_EG_reaching(un).st_reaches1{tt},neurons_EG_reaching(un).spk_trials_reaches1{tt}./tt,10,'filled','k');
    xline(0,'Color',[.8 .8 .8]);
    xlim([-4 3])
    set(gca,axeOpt{:})
    title('isi')

    subplot(322)
    scatter(cell2mat(neurons_EG_reaching(un).st_reaches1),cell2mat(neurons_EG_reaching(un).spk_trials_reaches1),3,'filled','k');
    xline(0,'Color',[.8 .8 .8]);
    xlim([-4 3])
    ylim([0 n_reaches1])
    set(gca,axeOpt{:})
    title('spike times - all reaches')

    subplot(323)
    plot(bin_edges,reaches1_FR(:,tt,un))
    xline(0,'Color',[.8 .8 .8]);
    xlim([-4 3])
    set(gca,axeOpt{:})
    title('FR')

    subplot(324)
    plot(bin_edges,reaches1_FR(:,:,un),'color',cat(2,rcl_clr,trans),LineWidth=.5)
    xlim([-4 3])
    set(gca,axeOpt{:})
    title('FR all - reaches')

    subplot(325)
    plot(sm_bin_edges,reaches1_FR_smooth(:,tt,un))
    xline(0,'Color',[.8 .8 .8]);
    xlim([-4 3])
    set(gca,axeOpt{:})
    xlabel('time from reach peak (s)')
    title('FR smooth')

    subplot(326)
    plot(sm_bin_edges,reaches1_FR_smooth(:,:,un),'color',cat(2,rcl_clr,trans),LineWidth=.5); hold on
    plot(sm_bin_edges,mean(reaches1_FR_smooth(:,:,un),2,'omitnan'),'color','k',LineWidth=2);
    xline(0,'Color',[.8 .8 .8]);
    axis tight
    xlim([-4 3])
    set(gca,axeOpt{:})
    xlabel('time from reach peak (s)')
    title('FR smooth - all reaches')


    set(gcf,'Position',[2383          47        1174         954],'color','w')
    saveas(gcf,strcat(save_example_cell_trial,filesep,'exampl_trial',num2str(neurons(un).phyID),'.png'),'png');
end

%% Averaging trials 

save_poll_trials = strcat(save_path,filesep,'pool_trials');
if ~exist(save_poll_trials,"dir"), mkdir(save_poll_trials); end

for un = 1:n_eg
    figure
    subplot(221)
    plot(sm_bin_edges,reaches1_FR_smooth(:,:,un),'color',cat(2,rcl_clr,trans),LineWidth=.5)
    xlim([-4 3])
    xline(0,'Color',[.8 .8 .8]);
    set(gca,axeOpt{:})
    title('FR individual trials smoothed')

    subplot(222)
    plot(bin_edges,mean(reaches1_FR(:,:,un),2,'omitnan'),LineWidth=.5)
    xlim([-4 3])
    xline(0,'Color',[.8 .8 .8]);
    set(gca,axeOpt{:})
    title('mean non-smoothed FR')

    subplot(223)
    plot(bin_edges,median(reaches1_FR(:,:,un),2,'omitnan'),LineWidth=.5)
    xlim([-4 3])
    xline(0,'Color',[.8 .8 .8]);
    set(gca,axeOpt{:})
    xlabel('time from reach peak (s)')
    title('median non-smoothed FR')

    subplot(224)
    plot(sm_bin_edges,mean(reaches1_FR_smooth(:,:,un),2,'omitnan'),LineWidth=.5)
    xlim([-4 3])
    xline(0,'Color',[.8 .8 .8]);
    set(gca,axeOpt{:})
    xlabel('time from reach peak (s)')
    title('mean smoothed FR')

    set(gcf,'Position',[2383         131        1349         807],'color','w')

    saveas(gcf,strcat(save_poll_trials,filesep,'exampl_pooling',num2str(neurons(un).phyID),'.png'),'png');
end

% 
% %% Save
% save(strcat(save_path,filesep,'neurons_reaches1_',region,'.mat'),"neurons_EG_reaching","reaches1_FR_params");
% 

%% Prepare reaches and neurons for MI

interval_reach = [-2 1];
% Reaches
flag_totaltime_win_paw = sm_frames_time >= interval_reach(1) & ...
    sm_frames_time <= interval_reach(2);
time_paw = sm_frames_time(flag_totaltime_win_paw);
paw_r1 = paw_pc1_smooth(flag_totaltime_win_paw,:);

% Neurons
flag_totaltime_win_neu = sm_bin_edges >= interval_reach(1) & ...
    sm_bin_edges <= interval_reach(2);
time_neu = sm_bin_edges(flag_totaltime_win_neu);
neu_r1 = reaches1_FR_smooth(flag_totaltime_win_neu,:,:);

%%
save(strcat(save_path,filesep,'paw_and_neurons_reaches1_',region,'.mat'));


% --------------------------------------------------------------
% MUTUTAL INFORMATION
%% --------------------
% Choose
% Choose
win_size = .05;
lag_interval = [-.5 .5];
dur_interval = interval_reach(2)-interval_reach(1);
n_lags = ((lag_interval(2)-lag_interval(1))/win_size)+1;
MI_landscape_size = dur_interval/win_size;

% initialize variables
start_stop_neu = nan(MI_landscape_size,MI_landscape_size,n_lags,2);
start_stop_paw = nan(MI_landscape_size,MI_landscape_size,n_lags,2);
lag_paw_neu = nan(n_lags,1);
MI_paw_neu = nan(MI_landscape_size,MI_landscape_size,n_lags,n_eg);
% Entropy
hist_nbins = 20;

for l = 1:n_lags
    lag = lag_interval(1)+ win_size*(l-1);

    for i = 1:MI_landscape_size
        start_win = interval_reach(1) + win_size*(i-1);
        start_lagged = start_win+lag;

        if (start_lagged>=interval_reach(1) && start_lagged<interval_reach(2))

            for j = 1:MI_landscape_size
                stop_win = interval_reach(1) + win_size*j;
                stop_lagged = stop_win+lag;

                if (stop_win>start_win && stop_lagged<=interval_reach(2))

                    % Start, stop locations and lags
                    start_stop_neu(i,j,l,1) = start_win;
                    start_stop_neu(i,j,l,2) = stop_win;
                    start_stop_paw(i,j,l,1) = start_lagged;
                    start_stop_paw(i,j,l,2) = stop_lagged;
                    lag_paw_neu(l) = lag;

                    % CALCULATE MI
                    % entropy paw
                    flag_win_paw = time_paw >= start_lagged & time_paw <= stop_lagged;
                    paw_win_dist = mean(paw_r1(flag_win_paw,:),1);
                    [pPaw, pawHist_edges] = histcounts(paw_win_dist,hist_nbins,...
                        'Normalization','probability');
                    entropy_paw = -sum(pPaw.*log2(pPaw),'omitnan');

                    % Entropy neuron
                    %n = 40;
                    for n = 1:n_eg
                        flag_win_neu = time_neu >= start_win & time_neu <= stop_win;
                        neu_win_dist = mean(neu_r1(flag_win_neu,:,n),1);
                        [pNeu, neuHist_edges] = histcounts(neu_win_dist,hist_nbins,...
                            'Normalization','probability');
                        entropy_neu = -sum(pNeu.*log2(pNeu),'omitnan');

                        % Entropy joint
                        [pJointPawNeu, ~, ~] = histcounts2(paw_win_dist,neu_win_dist,hist_nbins,...
                            'Normalization','probability');
                        entropy_joint = -sum(pJointPawNeu(:).*log2(pJointPawNeu(:)),'omitnan');

                        MI_paw_neu(i,j,l,n) = entropy_paw + entropy_neu - entropy_joint;
                    end
                end
            end
        end
    end
end

%%
figure

fi = tiledlayout(4,1);

ax1 = nexttile;
plot(time_paw, paw_r1,'color',cat(2,rcl_clr,trans),'LineWidth',.5); hold on
plot(time_paw, mean(paw_r1,2),'color',rcl_clr,'LineWidth',1.5); hold off
set(gca,axeOpt{:})

ax2 = nexttile;
plot(time_neu, neu_r1(:,:,n),'color',cat(2,rcl_clr,trans),'LineWidth',.5); hold on
plot(time_neu, mean(neu_r1(:,:,n),2,'omitnan'),'color',rcl_clr,'LineWidth',1.5); hold off
set(gca,axeOpt{:})

ax3 = nexttile;
ax3.Layout.TileSpan = [2,1];
y_im = start_stop_all(:,end,1);
x_im = start_stop_all(1,:,2);
imagesc(x_im,y_im,flipud(MI_paw_neu(:,:,n)));
set(gca, 'YTickDir', 'reverse')
xline(0,'--','color',[1 1 1],'LineWidth',1.5)
ylabel('start time (s)'); xlabel('stop time (s)')
set(gca,axeOpt{:})
axis square
set(gcf,'position',[2335          47         454         954],'color','w')

%%
save_entropy = strcat(save_path,filesep,'entropy_example');
if ~exist(save_entropy,"dir"), mkdir(save_entropy); end


figure
subplot(131)
histogram(paw_win_dist,pawHist_edges,'Normalization','probability')
title(sprintf('%s%.2f','H(paw) = ',entropy_paw))
xlabel('paw pc1'); ylabel('probability')
set(gca,axeOpt{:})

subplot(132)
histogram(neu_win_dist,neuHist_edges,'Normalization','probability')
title({ ...
    sprintf('%s%.2f','I (paw;neuron) = ',MI_paw_neu), ...
    sprintf('%s%.2f','H(neuron) = ',entropy_neu)})
xlabel('neuron FR'); ylabel('probability')
set(gca,axeOpt{:})

subplot(133)
histogram2(paw_win_dist,neu_win_dist,pawHist_edges,neuHist_edges,'Normalization','probability')
title(sprintf('%s%.2f','H(paw,neuron) = ',entropy_joint))
grid off
xlabel('paw pc1'); ylabel('neuron FR'); zlabel('probability')
set(gca,axeOpt{:})
set(gcf,'position',[2094         281        1594         531],'color','w')

saveas(gcf,strcat(save_entropy,filesep,num2str(neurons_EG_reaching(n).phyID),'.png'),'png');





