clear; close all; clc

%% Load behavior
behavior_root ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
ephys_root = 'E:'; %group_ephys = '20230801_ChocolateGroup';
group_setup_behav = strcat('20230511_ChocolateGroup',filesep,'headfixed_dynamicTarget');
mouse = '4_Milka';
paw_pref = 'right';
session = 'R4';
imec_id = '0';
ephys_sess = '18082023_Milka_StrCer_S4_g0';
sorter_folder = 'catGT\kilosort4';

% Load behavior
behavior_path = strcat(behavior_root,filesep,group_setup_behav,filesep,mouse,filesep,session);
load(strcat(behavior_path,filesep,'behavior_session.mat'));

% Load neurons
path_neurons  = strcat(ephys_root,filesep,ephys_sess,filesep,ephys_sess,'_imec',num2str(imec_id),filesep,sorter_folder);
load(strcat(path_neurons,filesep,'neurons_session.mat'));

save_path = strcat(path_neurons,filesep,'CS_neurons');
if ~exist("save_path","dir"), mkdir(save_path); end

%% Behavior - Oganize trials
nr_trials = behavior.behavior_duration.trial_end;
trials_vec = 1:nr_trials;

push_idx = behavior.init.idx_trial_push(behavior.init.idx_trial_push<=nr_trials);
pull_idx = behavior.init.idx_trial_pull(behavior.init.idx_trial_pull<=nr_trials);
left_idx = behavior.reach.left_idx(behavior.reach.left_idx<=nr_trials);
right_idx = behavior.reach.right_idx(behavior.reach.right_idx<=nr_trials);
center_idx = behavior.reach.center_idx(behavior.reach.center_idx<=nr_trials);

% intersection
push_left_idx = intersect(push_idx,left_idx);
push_right_idx = intersect(push_idx,right_idx);
push_center_idx = intersect(push_idx,center_idx);

pull_left_idx = intersect(pull_idx,left_idx);
pull_right_idx = intersect(pull_idx,right_idx);
pull_center_idx = intersect(pull_idx,center_idx);

% push & pull bounds
push_bounds = behavior.init.push_bounds;
push_bounds(push_bounds>nr_trials)=nr_trials;
pull_bounds = behavior.init.pull_bounds;
pull_bounds(pull_bounds>nr_trials)=nr_trials;

%reach bounds
left_bounds =  behavior.reach.left_bounds;
right_bounds = behavior.reach.right_bounds;
center_bounds = behavior.reach.center_bounds;

% window of interest
xlin_w = 0.02;
tm_before = -0.5;
tm_after = 0.5;
line_x = [tm_before;tm_before]-xlin_w;

% bins for probability of complex spike
binsP_width = 0.05;
binP_edges = tm_before:binsP_width:tm_after;
bin_med = binP_edges(1:end-1)+binsP_width/2;

%% Smooting
% Create the Gaussian kernel for smoothing
sigma = 0.02;  % Standard deviation of the Gaussian kernel
kernel_size = ceil(3 * sigma / binsP_width);  % Size of the kernel (covers +/- 3 sigma)
kernel_bw = 0.01;
kernel_range = 3 * sigma;  % Range to cover ±3 sigma
num_points = ceil(kernel_range / binsP_width) * 2 + 1;  % Ensure an odd number of points

% Generate x values for the kernel
x_kernel = linspace(-kernel_range, kernel_range, num_points);

% Gaussian formula
gaussian_kernel = exp(-x_kernel.^2 / (2 * sigma^2));
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel);  % Normalize to sum to 1

%% Load neuron of interest
idx_CSs = [neurons.potentialCS] == 1;
CS_neurons = neurons(idx_CSs);
nr_cs_uns = length(CS_neurons);

for cun = 1:nr_cs_uns

    cs = CS_neurons(cun);

    %% Probability of complex spike
    binP_counts_init = nan(length(binP_edges)-1,nr_trials);
    binP_counts_reach = nan(length(binP_edges)-1,nr_trials);
    for tt = 1:nr_trials
        binP_counts_init(:,tt) = histcounts(cs.st_init{tt},binP_edges);
        binP_counts_reach(:,tt) = histcounts(cs.st_reach{tt},binP_edges);
    end

    P_CS_init_push =  sum(binP_counts_init(:,push_idx),2)/length(push_idx);
    P_CS_init_pull =  sum(binP_counts_init(:,pull_idx),2)/length(pull_idx);
    P_CS_init_left =  sum(binP_counts_init(:,left_idx),2)/length(left_idx);
    P_CS_init_center =  sum(binP_counts_init(:,center_idx),2)/length(center_idx);
    P_CS_init_right =  sum(binP_counts_init(:,right_idx),2)/length(right_idx);

    P_CS_reach_push =  sum(binP_counts_reach(:,push_idx),2)/length(push_idx);
    P_CS_reach_pull =  sum(binP_counts_reach(:,pull_idx),2)/length(pull_idx);
    P_CS_reach_left =  sum(binP_counts_reach(:,left_idx),2)/length(left_idx);
    P_CS_reach_center =  sum(binP_counts_reach(:,center_idx),2)/length(center_idx);
    P_CS_reach_right =  sum(binP_counts_reach(:,right_idx),2)/length(right_idx);

    % Intersection - init
    pCS_init_push_L =  sum(binP_counts_init(:,push_left_idx),2)/length(push_left_idx);
    pCS_init_push_R =  sum(binP_counts_init(:,push_right_idx),2)/length(push_right_idx);
    pCS_init_push_C =  sum(binP_counts_init(:,push_center_idx),2)/length(push_center_idx);

    pCS_init_pull_L =  sum(binP_counts_init(:,pull_left_idx),2)/length(pull_left_idx);
    pCS_init_pull_R =  sum(binP_counts_init(:,pull_right_idx),2)/length(pull_right_idx);
    pCS_init_pull_C =  sum(binP_counts_init(:,pull_center_idx),2)/length(pull_center_idx);

    % Intersection - reach
    pCS_reach_push_L =  sum(binP_counts_reach(:,push_left_idx),2)/length(push_left_idx);
    pCS_reach_push_R =  sum(binP_counts_reach(:,push_right_idx),2)/length(push_right_idx);
    pCS_reach_push_C =  sum(binP_counts_reach(:,push_center_idx),2)/length(push_center_idx);

    pCS_reach_pull_L =  sum(binP_counts_reach(:,pull_left_idx),2)/length(pull_left_idx);
    pCS_reach_pull_R =  sum(binP_counts_reach(:,pull_right_idx),2)/length(pull_right_idx);
    pCS_reach_pull_C =  sum(binP_counts_reach(:,pull_center_idx),2)/length(pull_center_idx);


    %% Smooth
    % Create the Gaussian kernel
    % sigma = 0.02;  % Standard deviation of the Gaussian kernel
    % kernel_size = ceil(3 * sigma / binsP_width);  % Size of the kernel (covers +/- 3 sigma)
    % kernel_bw = 0.01;
    % kernel_range = 3 * sigma;  % Range to cover ±3 sigma
    % num_points = ceil(kernel_range / binsP_width) * 2 + 1;  % Ensure an odd number of points
    %
    % % Generate x values for the kernel
    % x_kernel = linspace(-kernel_range, kernel_range, num_points);
    %
    % % Gaussian formula
    % gaussian_kernel = exp(-x_kernel.^2 / (2 * sigma^2));
    % gaussian_kernel = gaussian_kernel / sum(gaussian_kernel);  % Normalize to sum to 1
    %
    % % Convolve the data with the Gaussian kernel
    % pCS_init_push_L_smooth = conv(pCS_init_push_L, gaussian_kernel, 'same');

    % figure
    % subplot(121)
    % plot(x_kernel,gaussian_kernel);
    % subplot(122)
    % stairs(binP_edges,[pCS_init_push_L;pCS_init_push_L(end)])
    % hold on
    % bin_med = binP_edges(1:end-1)+binsP_width/2;
    % plot(bin_med,pCS_init_push_L_smooth)

    %
    % Intersection - init
    spCS_init_push_L =  conv(pCS_init_push_L, gaussian_kernel, 'same');
    spCS_init_push_R =  conv(pCS_init_push_R, gaussian_kernel, 'same');
    spCS_init_push_C = conv(pCS_init_push_C, gaussian_kernel, 'same');

    spCS_init_pull_L = conv(pCS_init_pull_L, gaussian_kernel, 'same');
    spCS_init_pull_R =  conv(pCS_init_pull_R, gaussian_kernel, 'same');
    spCS_init_pull_C =  conv(pCS_init_pull_C, gaussian_kernel, 'same');

    % Intersection - reach
    spCS_reach_push_L =  conv(pCS_reach_push_L, gaussian_kernel, 'same');
    spCS_reach_push_R =  conv(pCS_reach_push_R, gaussian_kernel, 'same');
    spCS_reach_push_C =  conv(pCS_reach_push_C, gaussian_kernel, 'same');

    spCS_reach_pull_L =  conv(pCS_reach_pull_L, gaussian_kernel, 'same');
    spCS_reach_pull_R =  conv(pCS_reach_pull_R, gaussian_kernel, 'same');
    spCS_reach_pull_C =  conv(pCS_reach_pull_C, gaussian_kernel, 'same');


    %% Save
    CS_neurons(cun).params.tm_before = tm_before;
    CS_neurons(cun).params.tm_after = tm_after;

    CS_neurons(cun).params.binsP_width = binsP_width;
    CS_neurons(cun).params.binP_edges = binP_edges;

    CS_neurons(cun).binP_counts_init = binP_counts_init;
    CS_neurons(cun).binP_counts_reach = binP_counts_reach;

    CS_neurons(cun).spCS_init_push_L = spCS_init_push_L;
    CS_neurons(cun).spCS_init_push_R = spCS_init_push_R;
    CS_neurons(cun).spCS_init_push_C = spCS_init_push_C;
    CS_neurons(cun).spCS_init_pull_L = spCS_init_pull_L;
    CS_neurons(cun).spCS_init_pull_R = spCS_init_pull_R;
    CS_neurons(cun).spCS_init_pull_C = spCS_init_pull_C;

    CS_neurons(cun).spCS_reach_push_L = spCS_reach_push_L;
    CS_neurons(cun).spCS_reach_push_R = spCS_reach_push_R;
    CS_neurons(cun).spCS_reach_push_C = spCS_reach_push_C;
    CS_neurons(cun).spCS_reach_pull_L = spCS_reach_pull_L;
    CS_neurons(cun).spCS_reach_pull_R = spCS_reach_pull_R;
    CS_neurons(cun).spCS_reach_pull_C = spCS_reach_pull_C;

    %% Clrs
    push_clr = behavior.colors.push_clr;
    pull_clr =  behavior.colors.pull_clr;
    pp_clr = [63 130 109]./256;
    right_color = behavior.colors.right_color;
    center_color = behavior.colors.center_color;
    left_color = behavior.colors.left_color;
    rcl_clr = [44 123 182]./256;
    lw=2;
    patch_wd = 8;
    spk_size = 6;
    event_size = 5;
    psth_lim = .35;

    tm_widh = (tm_after-tm_before)+xlin_w*2;
    transp_push = .8;
    transp_pull = .5;
    ex_tr = .1;

    axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
        0.05,'ticklength',[1,1]*.01,'fontsize',10};

    % Figure
    figure
    tt = tiledlayout(3,2);
    %tt = tiledlayout(4,2);

    ax1 = nexttile;

    h1 = plot(bin_med,spCS_init_push_L,'LineWidth',lw,'Color',[left_color,transp_push+ex_tr],'LineStyle','-'); hold on;
    h2 = plot(bin_med,spCS_init_push_C,'LineWidth',lw,'Color',[center_color,transp_push+ex_tr],'LineStyle','-'); hold on;
    h3 = plot(bin_med,spCS_init_push_R,'LineWidth',lw,'Color',[right_color,transp_push+ex_tr],'LineStyle','-'); hold on;

    h4 = plot(bin_med,spCS_init_pull_L,'LineWidth',lw,'Color',[left_color,transp_pull+ex_tr],'LineStyle','--');
    h5 = plot(bin_med,spCS_init_pull_C,'LineWidth',lw,'Color',[center_color,transp_pull+ex_tr],'LineStyle','--');
    h6 = plot(bin_med,spCS_init_pull_R,'LineWidth',lw,'Color',[right_color,transp_pull+ex_tr],'LineStyle','--');

    ylabel('P(CS)'); xlabel('time from trial init (s)');
    xlim([tm_before tm_after])
    ylim([0 psth_lim])
    xline(0,'-','Color',[pp_clr 0.1],'linewidth',2);
    set(gca,axeOpt{:})
    title('trial initiation (push/pull)')

    leg=legend([h1, h2, h3, h4, h5 , h6],...
        'push left', 'push center', 'push right', ...
        'pull left', 'pull center', 'pull right', ...
        'location', 'northeast', 'box', 'off');
    leg.ItemTokenSize = [20, 10];
    leg.Position = [0.42, 0.8, 0.1, 0.16];
    %leg.Layout.Tile = 'east';

    ax2 = nexttile;
    h1 = plot(bin_med,spCS_reach_push_L,'LineWidth',lw,'Color',[left_color,transp_push+ex_tr],'LineStyle','-'); hold on;
    h2 = plot(bin_med,spCS_reach_push_C,'LineWidth',lw,'Color',[center_color,transp_push+ex_tr],'LineStyle','-'); hold on;
    h3 = plot(bin_med,spCS_reach_push_R,'LineWidth',lw,'Color',[right_color,transp_push+ex_tr],'LineStyle','-'); hold on;

    h4 = plot(bin_med,spCS_reach_pull_L,'LineWidth',lw,'Color',[left_color,transp_pull+ex_tr],'LineStyle','--');
    h5 = plot(bin_med,spCS_reach_pull_C,'LineWidth',lw,'Color',[center_color,transp_pull+ex_tr],'LineStyle','--');
    h6 = plot(bin_med,spCS_reach_pull_R,'LineWidth',lw,'Color',[right_color,transp_pull+ex_tr],'LineStyle','--');


    ylabel('P(CS)'); xlabel('time from water collection (s)');
    xlim([tm_before tm_after])
    ylim([0 psth_lim])
    xline(0,'-','Color',[rcl_clr 0.1],'linewidth',2);
    set(gca,axeOpt{:})
    title('water collection (reach)')
    % leg=legend([h4, h5, h6],'pull left', 'pull center', 'pull right', ...
    %     'location', 'northeast', 'box', 'off');
    % leg.ItemTokenSize = [20, 10];

    % leg=legend([h1, h2, h3, h4, h5 , h6],...
    %     'push left', 'push center', 'push right', ...
    %     'pull left', 'pull center', 'pull right', ...
    %     'location', 'northeast', 'box', 'off');
    % leg.ItemTokenSize = [20, 10];
    % leg.Position = [0.86, 0.80, 0.1, 0.16];

    ax3 = nexttile;
    ax3.Layout.TileSpan = [2,1];
    scatter(cell2mat(cs.st_init),cell2mat(cs.spk_trials_init),spk_size,'k','Marker','|','LineWidth',1), hold on
    scatter(zeros(size(trials_vec)),trials_vec,event_size,'filled','MarkerFaceColor',pp_clr)
    plot(repmat(line_x,[1 size(right_bounds,2)]),right_bounds,'linewidth',patch_wd,'color',right_color)
    plot(repmat(line_x,[1 size(left_bounds,2)]),left_bounds,'linewidth',patch_wd,'color',left_color)
    plot(repmat(line_x,[1 size(center_bounds,2)]),center_bounds,'linewidth',patch_wd,'color',center_color)
    hold off
    axis tight;
    xlabel('time from trial init (s)'); ylabel('trials in session')
    axis([tm_before-xlin_w*2 tm_after 0 nr_trials])
    set(gca,axeOpt{:})


    for n_push =1:size(push_bounds,2)
        pos_push = [tm_before-xlin_w*2 push_bounds(1,n_push) tm_widh ...
            push_bounds(2,n_push)-push_bounds(1,n_push)];
        rectangle('Position',pos_push, 'EdgeColor',[0.2 0.2 0.2 transp_push],'LineWidth',2,...
            'LineStyle','-');

        text(tm_after+.08, push_bounds(1,n_push)+2, 'push', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', 'FontSize', 11, 'Color',[0.2 0.2 0.2 transp_push], ...
            'Rotation', 90);
    end
    for n_pull =1:size(pull_bounds,2)
        pos_pull = [tm_before-xlin_w*2 pull_bounds(1,n_pull) tm_widh ...
            pull_bounds(2,n_pull)-pull_bounds(1,n_pull)];
        rectangle('Position',pos_pull, 'EdgeColor',[0.5 0.5 0.5 transp_pull],'LineWidth',2,...
            'LineStyle','--');

        text(tm_after+.08, pull_bounds(1,n_pull)+2, 'pull', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', 'FontSize', 11, 'Color',[0.5 0.5 0.5 transp_pull], ...
            'Rotation', 90);
    end


    %
    ax4 = nexttile;
    ax4.Layout.TileSpan = [2,1];
    scatter(cell2mat(cs.st_reach),cell2mat(cs.spk_trials_reach),spk_size,'k','Marker','|','LineWidth',1), hold on
    scatter(zeros(size(trials_vec)),trials_vec,event_size,'filled','MarkerFaceColor',rcl_clr)
    plot(repmat(line_x,[1 size(right_bounds,2)]),right_bounds,'linewidth',patch_wd,'color',right_color)
    plot(repmat(line_x,[1 size(left_bounds,2)]),left_bounds,'linewidth',patch_wd,'color',left_color)
    plot(repmat(line_x,[1 size(center_bounds,2)]),center_bounds,'linewidth',patch_wd,'color',center_color)
    hold off
    axis tight;
    axis([tm_before-xlin_w*2 tm_after 0 nr_trials])
    xlabel('time from water collection (s)'); ylabel('trials in session')
    set(gca,axeOpt{:})
    %axis square
    %set(gcf,'position', [ 2761         128         839         806],'color','w');
    set(gcf,'position', [2892 331 839 543],'color','w');

    for n_push =1:size(push_bounds,2)
        pos_push = [tm_before-xlin_w*2 push_bounds(1,n_push) tm_widh ...
            push_bounds(2,n_push)-push_bounds(1,n_push)];
        rectangle('Position',pos_push, 'EdgeColor',[0.2 0.2 0.2 transp_push],'LineWidth',2,...
            'LineStyle','-');

        text(tm_after+.08, push_bounds(1,n_push)+2, 'push', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', 'FontSize', 11, 'Color',[0.2 0.2 0.2 transp_push], ...
            'Rotation', 90);
    end
    for n_pull =1:size(pull_bounds,2)
        pos_pull = [tm_before-xlin_w*2 pull_bounds(1,n_pull) tm_widh ...
            pull_bounds(2,n_pull)-pull_bounds(1,n_pull)];
        rectangle('Position',pos_pull, 'EdgeColor',[0.5 0.5 0.5 transp_pull],'LineWidth',2,...
            'LineStyle','--');

        text(tm_after+.08, pull_bounds(1,n_pull)+2, 'pull', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', 'FontSize', 11, 'Color',[0.5 0.5 0.5 transp_pull], ...
            'Rotation', 90);
    end

    saveas(gcf,strcat(save_path,filesep,'neuron',num2str(cs.phyID),'_P_CS.png'),'png');
    saveas(gcf,strcat(save_path,filesep,'neuron',num2str(cs.phyID),'_P_CS'),'epsc');

end
%% Save
save(strcat(save_path,filesep,'cs_neurons.mat'),'CS_neurons');

