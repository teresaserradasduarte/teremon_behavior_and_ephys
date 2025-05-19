%% MUTUAL INFORMARTION
% May 2025
clear; close all; clc

%% Load data
% Group and individual
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
group = '20230801_ChocolateGroup';
setup = 'headfixed_dynamicTarget';
animals = {...
    'CoteDor',...
    'Lindt',...
    'Toblerone',...
    'Milka',...
    'FerreroRocher'};
reg = 'CB';
if strcmp(reg,'BG'), load_struct = 'N_BG_all';
elseif strcmp(reg,'CB'), load_struct = 'N_CB_all';
elseif strcmp(reg,'CT'), load_struct = 'N_CT_all';
end


%% Create path to save
save_mat = fullfile(rootdir,"ephys_and_behavior","mat_files",group);
save_out = fullfile(rootdir,"ephys_and_behavior","out_files",group,'group',reg,'MI');
if ~exist(save_out,"dir"), mkdir(save_out); end
if ~exist(save_mat,"dir"), mkdir(save_mat); end
save_name = strcat('PCA_',reg,'.mat');

%% Load neurons
load(fullfile(rootdir,"ephys_and_behavior/","mat_files/",group,"all_neurons_pooled.mat"),load_struct)
load(fullfile(rootdir,"ephys_and_behavior/","mat_files/",group,"1_CoteDor/","R6/","BG/","eg_neurons.mat"),'eg_neu_FR_params');
load(fullfile(rootdir,"behavior_data/","analyzed_data/mat_files/","20230511_ChocolateGroup/","headfixed_dynamicTarget/","1_CoteDor/","R6/","behavior_session.mat"))
load(fullfile(rootdir,"behavior_data/","analyzed_data/mat_files/","20230511_ChocolateGroup/","headfixed_dynamicTarget/","1_CoteDor/","R6/","session_reaching_data_paw.mat"),'reaches')

%% Neurons table
N_tmp = struct2table(N_CT_all);
% exclude the only session with just neurons from the cerebellum
exclude_idx = strcmp(N_tmp.mouse, '1_CoteDor') & strcmp(N_tmp.sess, 'R1');
N = N_tmp(~exclude_idx, :);

%% Window selection - FR
bin_edges = eg_neu_FR_params.bin_edges;

% new window
sel_win = [-1 1];
[~,new_bin_edges_start] = min(abs(bin_edges-sel_win(1)));
[~,new_bin_edges_stop] = min(abs(bin_edges-sel_win(2)));
FR_edges_range = new_bin_edges_start:new_bin_edges_stop;

% MI landscape
dur_interval = sel_win(2)-sel_win(1);
win_size = .05;
MI_landscape_size = dur_interval/win_size;

% Reach_time_win
tm_reach = reaches.tm_w;
win_size_endpoint = win_size;
start_win_bhv = -.01;
win_stop = start_win_bhv+win_size_endpoint;
[~,r_start] = min(abs(tm_reach-start_win_bhv));
[~,r_stop] = min(abs(tm_reach-win_stop));

% Bins for distribtuion acroos reaches
hist_nbins = 20;

% Save params
MI_params.win_MI = sel_win;
MI_params.win_size = win_size;
MI_params.endpoint_win = [start_win_bhv,win_stop];
MI_params.hist_nbins = hist_nbins;


%% Figure params
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};

%% Loop of sessions - get behavior only once
[unique_sessions, ~] = unique(N(:, {'mouse', 'sess'}), 'rows');
n_unique_sessions = height(unique_sessions);


bhv_entropy = struct();

for s = 1:height(unique_sessions)
    %s=1;
    curr_mouse = unique_sessions.mouse{s};
    curr_sess = unique_sessions.sess{s};
    fprintf('Getting behavior entropy of mouse %s , session %s\n', curr_mouse, curr_sess);

    % Logical indexing to select rows matching the current mouse & session
    idx = strcmp(N.mouse, curr_mouse) & strcmp(N.sess, curr_sess);
    N_sub = N(idx, :);  % Subtable for current combo mouse/session

    % Full reach idx
    r1_idx = cell2mat(N_sub(1,:).idx_reach_cat) == 1;
    %Trajectories and endpoint
    reach_all = cell2mat(N_sub(1,:).reach_px);
    water_loc = cell2mat(N_sub(1,:).idx_reach_LCR);
    reach_r1 = reach_all(:,:,r1_idx);
    water_loc_r1 = water_loc(r1_idx);
    n_r1 = length(water_loc_r1);
    endpoint_xyz = squeeze(mean(reach_r1(r_start:r_stop,:,:),1));

    % Entropy
    endpoint_y_dist = endpoint_xyz(2,:);
    [pR1, r1Hist_edges] = histcounts(endpoint_y_dist,hist_nbins,...
        'Normalization','probability');
    entropy_r1 = -sum(pR1.*log2(pR1),'omitnan');

    % -------------------------------------------
    % Figure
    save_out_bhv = fullfile(rootdir,"ephys_and_behavior/out_files/20230801_ChocolateGroup/",char(N_sub(1,:).mouse),char(N_sub(1,:).sess),'bhv');
    save_mat_bhv = fullfile(rootdir,"ephys_and_behavior/mat_files/20230801_ChocolateGroup/",char(N_sub(1,:).mouse),char(N_sub(1,:).sess));
    if ~exist(save_out_bhv,"dir"), mkdir(save_out_bhv); end

    % Reach heatmaps with window of interest
    fig_name = 'full_reaches_MIwin';
    colors_lcr=[behavior.colors.left_color;behavior.colors.center_color;behavior.colors.right_color];
    sz = 20;

    figure
    ff = tiledlayout(3,3);
    title(ff,sprintf('%s%s%s%s','Full reaches | mouse: ',...
        char(N_sub(1,:).mouse), ' | session: ',char(N_sub(1,:).sess)),...
        'Interpreter','none');

    xyz_label = ['x';'y';'z'];
    for nxt=1:3
        ax1=nexttile(nxt);
        ax1.Layout.TileSpan = [2,1];
        imagesc(tm_reach,1:n_r1, squeeze(reach_r1(:,nxt,:))');
        xline(0,'--','Color',[1 1 1 .5],'LineWidth',2);
        xline(tm_reach(r_start),'-','Color',[.5 .5 .5 .5],'LineWidth',1);
        xline(tm_reach(r_stop),'-','Color',[.5 .5 .5 .5],'LineWidth',1);
        xlabel('time (s)'); ylabel('reach idx');
        set(gca,axeOpt{:});
        title(xyz_label(nxt))
        c1=colorbar; ylabel(c1,sprintf('%s%s',xyz_label(nxt),' ( px)'));

        nexttile(nxt+6)
        scatter(1:n_r1,endpoint_xyz(nxt,:), sz, colors_lcr(water_loc_r1, :), 'filled'); hold off
        ylabel(c1,sprintf('%s%s',xyz_label(nxt))); xlabel('reach idx');
        set(gca,axeOpt{:});
    end
    set(gcf,'Position',[1988         193        1755         732],'color','w');
    saveas(ff,strcat(save_out_bhv,filesep,fig_name,'.png'),'png');
    print(gcf,strcat(save_out_bhv,filesep,fig_name,'.pdf'), '-dpdf', '-painters');

    % Histogram
    figure
    histogram(endpoint_y_dist,hist_nbins,'Normalization','probability');
    title(sprintf('%s%.2f','Reach endpoint entropy = ',entropy_r1));
    set(gca,axeOpt{:});
    xlabel('endpoint in y (px)'); ylabel('probability');
    set(gcf,'Position',[3036         377         560         420],'color','w');
    saveas(gcf,strcat(save_out_bhv,filesep,'entropy_endpoint_y','.png'),'png');

    % Save varibles
    bhv_entropy(s).mouse = curr_mouse;
    bhv_entropy(s).sess = curr_sess;
    bhv_entropy(s).r1_idx = r1_idx;
    bhv_entropy(s).endpoint_y_dist = endpoint_y_dist;
    bhv_entropy(s).pR1 = pR1;
    bhv_entropy(s).r1Hist_edges = r1Hist_edges;
    bhv_entropy(s).entropy_r1 = entropy_r1;
    save(strcat(save_mat_bhv,filesep,'behavior_fundamentals.mat'),'bhv_entropy','-append');

end
bhv_table = struct2table(bhv_entropy);


%% Neurons loop
close all;
n_neu = height(N);
start_stop_neu = nan(MI_landscape_size,MI_landscape_size,2);
MI_pawEndpoint_neu = nan(MI_landscape_size,MI_landscape_size);
MI_landscaps_all = nan(MI_landscape_size,MI_landscape_size,n_neu);

MI_struct = struct();
%n=1;
% loop through neurons
for n = 1:n_neu
    fprintf('Running neuron %i\n', n);

    % Getting entropy and endpoint distributions from behavior
    idx_bhv = strcmp(N(n,:).mouse, bhv_table.mouse) & strcmp(N(n,:).sess, bhv_table.sess);
    entropy_paw = bhv_table(idx_bhv,:).entropy_r1;
    pawEndPoint_win_dist = cell2mat(bhv_table(idx_bhv,:).endpoint_y_dist);
    idx_reaches = cell2mat(bhv_table(idx_bhv,:).r1_idx);

    for i = 1:MI_landscape_size
        start_win = sel_win(1) + win_size*(i-1);
        if (start_win>=sel_win(1) && start_win<sel_win(2))
            for j = 1:MI_landscape_size
                stop_win = sel_win(1) + win_size*j;
                if (stop_win>start_win && stop_win<=sel_win(2))

                    % Start, stop locations
                    if n==1
                        start_stop_neu(i,j,1) = start_win;
                        start_stop_neu(i,j,2) = stop_win;
                    end

                    % Get FR of full reaches for that neuron
                    FR_reach_all = cell2mat(N(n,:).FR_reach);
                    FR_r1 = FR_reach_all(:,idx_reaches);

                    % Get the distribution within the window i,j
                    flag_win_neu = bin_edges >= start_win & bin_edges <= stop_win;
                    neu_win_dist = mean(FR_r1(flag_win_neu,:),1,"omitnan");

                    % CALCULATE MI
                    % Entropy neuron
                    [pNeu, neuHist_edges] = histcounts(neu_win_dist,hist_nbins,...
                        'Normalization','probability');
                    entropy_neu = -sum(pNeu.*log2(pNeu),'omitnan');

                    % Entropy joint
                    [pJointPawNeu, ~, ~] = histcounts2(pawEndPoint_win_dist,neu_win_dist,hist_nbins,...
                        'Normalization','probability');
                    entropy_joint = -sum(pJointPawNeu(:).*log2(pJointPawNeu(:)),'omitnan');

                    MI_pawEndpoint_neu(i,j) = entropy_paw + entropy_neu - entropy_joint;
                end
            end
        end

    end

    MI_landscaps_all(:,:,n) = MI_pawEndpoint_neu;

    MI_struct(n).pawEndPoint_win_dist = pawEndPoint_win_dist;
    MI_struct(n).neu_win_dist = neu_win_dist;
    %MI_struct(n).entropy_paw = entropy_paw;
    %MI_struct(n).entropy_neu = entropy_neu;
    MI_struct(n).entropy_joint = entropy_joint;
    MI_struct(n).MI_landscape = MI_pawEndpoint_neu;
    MI_struct(n).start_stop_neu = start_stop_neu;
    MI_struct(n).mouse = mouse;
    MI_struct(n).sess = sess;
    

    % Figures
    % Landscape figure
    figure(1)
    y_im = start_stop_neu(:,end,1);
    x_im = start_stop_neu(1,:,2);
    imagesc(x_im,y_im,MI_pawEndpoint_neu);
    axis xy;        % Invert the y-axis direction
    hold on
    %plot(stop_n,start_n,'ko','LineWidth',.5)
    %imagesc(x_im,y_im,flipud(MI_paw_neu_smooth));
    title({sprintf('%s%s','region: ',reg),...
        'Mutual information: I(paw;neuron)'})

    %set(gca, 'YTickDir', 'reverse')
    c=colorbar;
    ylabel(c,'Bits');
    xline(0,'--','color',[1 1 1],'LineWidth',1.5)
    yline(0,'--','color',[1 1 1],'LineWidth',1.5)
    ylabel('start time (s)'); xlabel('stop time (s)')
    set(gca,axeOpt{:})
    axis square
    set(gcf,'position',[2857         381         597         525],'color','w')

    saveas(gcf,strcat(save_out,filesep,'MIland_neuron_',...
        num2str(n),'.png'),'png');
end


%
save(strcat(save_mat,filesep,'MI_',reg,'.mat'),'MI_struct','MI_landscaps_all','-v7.3');
