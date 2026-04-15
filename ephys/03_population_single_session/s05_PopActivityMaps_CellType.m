%% POPULATION ANALYSIS - PER CELL TYPE
% heatmaps of all cells of each region separated by cell type
clear; close all; clc

%% Load Data
rootdir ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
group = '20230801_ChocolateGroup';
datapath = fullfile(rootdir,"ephys_and_behavior","mat_files",group);
reg = 'DCN';
load_file = sprintf('%s_%s_%s','all',reg,'neurons');

% Load the data
S = load(fullfile(datapath,'eg_neurons.mat'),load_file,'eg_neu_FR_params');

% Base Output Path
base_save_out = fullfile(rootdir,"ephys_and_behavior","out_files",group,'group','activityMaps',reg);
if ~exist(base_save_out,"dir"), mkdir(base_save_out); end

%% Global Table and Params
N_full = struct2table(S.(load_file));

axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};
bin_edges_orig = S.eg_neu_FR_params.bin_edges;
nr_bins_original = S.eg_neu_FR_params.nr_bins;

% Define a blue-white-red colormap for difference plots if you don't have one loaded
bwr_map = @(n) interp1([1 50 100], [0 0 1; 1 1 1; 1 0 0], linspace(1, 100, n)); 
bwr_cmap = bwr_map(64); % Use this variable for difference plots

%% MASTER LOOP: Iterate over Cell Types
% Get list of unique cell types, removing NaNs
cell_types_list = unique(string(N_full.cell_type));
cell_types_list(ismissing(cell_types_list) | cell_types_list == "NaN") = [];

fprintf('Found %d cell types: %s\n', length(cell_types_list), join(cell_types_list, ', '));

for t = 1:length(cell_types_list)
    curr_type = char(cell_types_list(t));
    fprintf('Processing Cell Type: %s ...\n', curr_type);
    
    %% 1. Filter Data for this Cell Type
    N = N_full(strcmp(N_full.cell_type, curr_type), :);
    n_neu = size(N,1);
    
    % Skip if too few neurons for PCA
    if n_neu < 3
        fprintf('  -> Skipping %s: Not enough neurons (%d).\n', curr_type, n_neu);
        continue; 
    end
    
    % Filter for Init analysis (Removing specific mouse)
    N_pp = N(~strcmp(N.mouse, '3_Toblerone'), :);
    n_neu_pp = size(N_pp,1);

    %% 2. Create Directory Structure
    % Structure: CB > CellType > init / reach
    type_save_out = fullfile(base_save_out, curr_type);
    
    save_out_init = fullfile(type_save_out, 'init');
    if ~exist(save_out_init,"dir"), mkdir(save_out_init); end
    
    save_out_reach = fullfile(type_save_out, 'reach');
    if ~exist(save_out_reach,"dir"), mkdir(save_out_reach); end

    %% ====================================================================
    %% INIT PHASE ANALYSIS
    %% ====================================================================
    if n_neu_pp > 2 % Safety check for Init phase
        
        % --- PCA ---
        FR_init_mean_tmp = cellfun(@(x) mean(x, 2), N_pp.FR_init, 'UniformOutput', false);
        FR_init_mean = cat(2, FR_init_mean_tmp{:});
        [zs_FR_init, FR_init_mu, FR_init_sigma] = zscore(FR_init_mean);
        [coeff_init, score_init, latent_init, tsquare_init, explained_init, mus_init] = pca(zs_FR_init);
        
        % Order of cells
        r_coeff_pp = sqrt(coeff_init(:,1).^2 + coeff_init(:,2).^2);
        theta_coef_pp = atan2(coeff_init(:,2), coeff_init(:,1));
        [theta_coeffPP_sorted, neu_order_pp] = sort(theta_coef_pp);
        cmap = parula(length(neu_order_pp));
        
        % Find the starting point
        [~,break_rad] = max(diff(theta_coeffPP_sorted));
        n_order_pp = [neu_order_pp(break_rad+1:end);neu_order_pp(1:break_rad)];
        
        % --- PCA Plots ---
        fig = figure('Visible', 'off'); % Set visible off to prevent window spam
        subplot(121)
        plot(cumsum(explained_init),'-ko');
        xlabel('PCs'); ylabel('explained'); title(curr_type)
        subplot(122)
        plot(cumsum(score_init(:,1:2)));
        xlabel('time'); ylabel('score')
        legend('PC1','PC2')
        set(gcf,'Position',[ 2762 134 1027 423])
        saveas(gcf,fullfile(save_out_init,'PCA_pushPull.png'),'png');
        close(fig);
        
        fig = figure('Visible', 'off');
        subplot(121)
        scatter(coeff_init(:,1),coeff_init(:,2),50,cmap,'filled');
        xlabel('PC1 coefficients'); ylabel('PC2 coefficients')
        set(gca,axeOpt{:})
        axis square
        title([curr_type ' - Polar Coords'])
        hold on
        for i=1:n_neu_pp
            subplot(122)
            polarplot(theta_coeffPP_sorted(i), r_coeff_pp(neu_order_pp(i)),...
                'o','MarkerSize', 8, 'MarkerFaceColor',cmap(i,:),'markerEdgeColor','w'); hold on
        end
        title('Angular position of PC1 vs PC2 coeffs')
        set(gcf,'Position',[2762 134 1027 734])
        saveas(gcf,fullfile(save_out_init,'PC1coeffs_vs_PC2coeffs_pushPull.png'),'png');
        close(fig);
        
        % --- SELECTION: PUSH / PULL ---
        FR_push_pull = cellfun(@(fr, idx) fr(:, idx==0 | idx==1), ...
            N_pp.FR_init, N_pp.idx_init_vPP_invPP, 'UniformOutput', false);
        zFR_push_pull = concat_zscore(FR_push_pull);
        
        FR_push_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx==0), 2, "omitnan"), ...
            zFR_push_pull, N_pp.idx_init_vPP_invPP, 'UniformOutput', false)');
        FR_pull_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx==1), 2, "omitnan"), ...
            zFR_push_pull, N_pp.idx_init_vPP_invPP, 'UniformOutput', false)');
        FR_diff_pp = FR_pull_mean_zs-FR_push_mean_zs;
        
        % --- PUSH PULL FIGURE ---
        fig_name = 'tile_push_pull';
        lim_clr = [-.7 1];
        lim_clr_diff = [-1 1];
        flip_ud = true;
        
        fig = figure('Visible', 'off');
        ax1 = subplot(131);
        imagesc(bin_edges_orig,1:n_neu_pp,FR_push_mean_zs(:,n_order_pp)',lim_clr);
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); ylabel('cells'); title({[curr_type ' Init'];'Push'})
        colormap(ax1, 'hot')
        
        ax2 = subplot(132);
        imagesc(bin_edges_orig,1:n_neu_pp,FR_pull_mean_zs(:,n_order_pp)',lim_clr)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Init'];'Pull'})
        colormap(ax2, 'hot')
        
        ax3 = subplot(133);
        imagesc(bin_edges_orig,1:n_neu_pp,FR_diff_pp(:,n_order_pp)',lim_clr_diff)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Init'];'Diff'})
        colormap(ax3, bwr_cmap)
        
        set(gcf,'position',[1972 212 1870 644],'color','w')
        saveas(gcf,fullfile(save_out_init,[fig_name,'.png']),'png');
        close(fig);
        
        % --- SELECTION: VALID / INVALID ---
        zFR_init = concat_zscore(N_pp.FR_init);
        FR_val_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx==0 | idx == 1), 2, "omitnan"), ...
            zFR_init, N_pp.idx_init_vPP_invPP, 'UniformOutput', false)');
        FR_inval_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx == 2 | idx == 3), 2, "omitnan"), ...
            zFR_init, N_pp.idx_init_vPP_invPP, 'UniformOutput', false)');
        FR_diff_val = FR_inval_mean_zs-FR_val_mean_zs;
        
        % --- VALID INVALID FIGURE ---
        fig_name = 'tile_valid_invalid';
        fig = figure('Visible', 'off');
        ax1 = subplot(131);
        imagesc(bin_edges_orig,1:n_neu_pp,FR_val_mean_zs(:,n_order_pp)',lim_clr);
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); ylabel('cells'); title({[curr_type ' Init'];'Valid'})
        colormap(ax1, 'hot')
        
        ax2 = subplot(132);
        imagesc(bin_edges_orig,1:n_neu_pp,FR_inval_mean_zs(:,n_order_pp)',lim_clr)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Init'];'Invalid'})
        colormap(ax2, 'hot')
        
        ax3 = subplot(133);
        imagesc(bin_edges_orig,1:n_neu_pp,FR_diff_val(:,n_order_pp)',lim_clr_diff)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Init'];'Diff'})
        colormap(ax3, bwr_cmap)
        
        set(gcf,'position',[1972 212 1870 644],'color','w')
        saveas(gcf,fullfile(save_out_init,[fig_name,'.png']),'png');
        close(fig);
    end
    
    %% ====================================================================
    %% REACH PHASE ANALYSIS
    %% ====================================================================
    if n_neu > 2 % Safety check for Reach phase
        
        % --- PCA ---
        FR_reach_mean_tmp = cellfun(@(x) mean(x, 2), N.FR_reach, 'UniformOutput', false);
        FR_reach_mean = cat(2, FR_reach_mean_tmp{:});
        [zs_FR_reach, FR_reach_mu, FR_reach_sigma] = zscore(FR_reach_mean);
        [coeff_reach, score_reach, latent_reach, tsquare_reach, explained_reach, mus_reach] = pca(zs_FR_reach);
        
        % Order of cells
        r_coeff_hDCnD = sqrt(coeff_reach(:,1).^2 + coeff_reach(:,2).^2);
        theta_coef_hDCnD = atan2(coeff_reach(:,2), coeff_reach(:,1));
        [theta_coeff_hDCnD_sorted, neu_order_hDCnD] = sort(theta_coef_hDCnD);
        cmap_hDCnD = parula(length(neu_order_hDCnD));
        
        % Find the starting point
        [~,break_rad_reach] = max(diff(theta_coeff_hDCnD_sorted));
        n_order_hDCnD = [neu_order_hDCnD(break_rad_reach+1:end);neu_order_hDCnD(1:break_rad_reach)];
        
        % --- PCA Plots ---
        fig = figure('Visible', 'off');
        subplot(121)
        plot(cumsum(explained_reach),'-ko');
        xlabel('PCs'); ylabel('explained'); title(curr_type)
        subplot(122)
        plot(cumsum(score_reach(:,1:2)));
        xlabel('time'); ylabel('score')
        legend('PC1','PC2')
        set(gcf,'Position',[ 2762 134 1027 423])
        saveas(gcf,fullfile(save_out_reach,'PCA_reach.png'),'png');
        close(fig);
        
        fig = figure('Visible', 'off');
        subplot(121)
        scatter(coeff_reach(:,1),coeff_reach(:,2),50,cmap_hDCnD,'filled');
        xlabel('PC1 coefficients'); ylabel('PC2 coefficients')
        set(gca,axeOpt{:})
        axis square
        title([curr_type ' - Polar Coords'])
        hold on
        for i=1:n_neu
            subplot(122)
            polarplot(theta_coef_hDCnD(i), r_coeff_hDCnD(n_order_hDCnD(i)),...
                'o','MarkerSize', 8, 'MarkerFaceColor',cmap_hDCnD(i,:),'markerEdgeColor','w'); hold on
        end
        title('Angular position of PC1 vs PC2 coeffs')
        set(gcf,'Position',[2762 134 1027 734])
        saveas(gcf,fullfile(save_out_reach,'PC1coeffs_vs_PC2coeffs_reach.png'),'png');
        close(fig);
        
        % --- DATA PREP ---
        zFR_reach = concat_zscore(N.FR_reach);
        
        % --- SELECTION: Dom/C/non-Dom ---
        FR_D_mean_zs = cell2mat(cellfun(@(fr, idx1, idx2) mean(fr(:, idx1 == 1 & idx2==1), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_hDCnD, N.idx_reach_cat, 'UniformOutput', false)');
        FR_C_mean_zs = cell2mat(cellfun(@(fr, idx1, idx2) mean(fr(:, idx1 == 2 & idx2==1), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_hDCnD, N.idx_reach_cat, 'UniformOutput', false)');
        FR_nD_mean_zs = cell2mat(cellfun(@(fr, idx1, idx2) mean(fr(:, idx1 == 3 & idx2==1), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_hDCnD, N.idx_reach_cat, 'UniformOutput', false)');
        FR_diff_DnD = FR_nD_mean_zs-FR_D_mean_zs;
        
        % --- Dom/C/nonDom FIGURE ---
        flip_ud = true;
        lim_clr = [-.8 1.1];
        fig_name = 'tile_reach_DCnD';
        
        fig = figure('Visible', 'off');
        ax1 = subplot(131);
        imagesc(bin_edges_orig,1:n_neu,FR_D_mean_zs(:,n_order_hDCnD)',lim_clr);
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); ylabel('cells'); title({[curr_type ' Reach'];'Dom'})
        colormap(ax1, 'hot')
        
        ax2 = subplot(132);
        imagesc(bin_edges_orig,1:n_neu,FR_C_mean_zs(:,n_order_hDCnD)',lim_clr)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Center'})
        colormap(ax2, 'hot')
        
        ax3 = subplot(133);
        imagesc(bin_edges_orig,1:n_neu,FR_nD_mean_zs(:,n_order_hDCnD)',lim_clr)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Non-Dom'})
        colormap(ax3,  'hot')
        
        set(gcf,'position',[1972 212 1870 644],'color','w')
        saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
        close(fig);
        
        % --- Dom/nonDom/diff FIGURE ---
        fig_name = 'tile_reach_DnD_wDiff';
        fig = figure('Visible', 'off');
        ax1 = subplot(131);
        imagesc(bin_edges_orig,1:n_neu,FR_D_mean_zs(:,n_order_hDCnD)',lim_clr);
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); ylabel('cells'); title({[curr_type ' Reach'];'Dom'})
        colormap(ax1, 'hot')
        
        ax2 = subplot(132);
        imagesc(bin_edges_orig,1:n_neu,FR_C_mean_zs(:,n_order_hDCnD)',lim_clr)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Center'})
        colormap(ax2, 'hot')
        
        ax3 = subplot(133);
        imagesc(bin_edges_orig,1:n_neu,FR_diff_DnD(:,n_order_hDCnD)',lim_clr_diff)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Diff'})
        colormap(ax3,  bwr_cmap)
        
        set(gcf,'position',[1972 212 1870 644],'color','w')
        saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
        close(fig);
        
        % --- SELECTION: L/C/R ---
        FR_L_mean_zs = cell2mat(cellfun(@(fr, idx1, idx2) mean(fr(:, idx1 == 1 & idx2==1), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_LCR, N.idx_reach_cat, 'UniformOutput', false)');
        FR_R_mean_zs = cell2mat(cellfun(@(fr, idx1, idx2) mean(fr(:, idx1 == 3 & idx2==1), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_LCR, N.idx_reach_cat, 'UniformOutput', false)');
        
        % --- L/C/R FIGURE ---
        fig_name = 'tile_reach_LCR';
        fig = figure('Visible', 'off');
        ax1 = subplot(131);
        imagesc(bin_edges_orig,1:n_neu,FR_L_mean_zs(:,n_order_hDCnD)',lim_clr);
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); ylabel('cells'); title({[curr_type ' Reach'];'Left'})
        colormap(ax1, 'hot')
        
        ax2 = subplot(132);
        imagesc(bin_edges_orig,1:n_neu,FR_C_mean_zs(:,n_order_hDCnD)',lim_clr)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Center'})
        colormap(ax2, 'hot')
        
        ax3 = subplot(133);
        imagesc(bin_edges_orig,1:n_neu,FR_R_mean_zs(:,n_order_hDCnD)',lim_clr)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Right'})
        colormap(ax3,  'hot')
        
        set(gcf,'position',[1972 212 1870 644],'color','w')
        saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
        close(fig);
        
        % --- SELECTION: Hit/miss ---
        FR_hit_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx == 1), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_hit, 'UniformOutput', false)');
        FR_miss_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx == 0), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_hit, 'UniformOutput', false)');
        FR_diff_hitMiss = FR_miss_mean_zs-FR_hit_mean_zs;
        
        % --- FIGURE HIT MISS ---
        fig_name = 'tile_hit_miss';
        lim_clr_hit = [-.6 1.3];
        fig = figure('Visible', 'off');
        ax1 = subplot(131);
        imagesc(bin_edges_orig,1:n_neu,FR_hit_mean_zs(:,n_order_hDCnD)',lim_clr_hit);
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); ylabel('cells'); title({[curr_type ' Reach'];'Hit'})
        colormap(ax1, 'hot')
        
        ax2 = subplot(132);
        imagesc(bin_edges_orig,1:n_neu,FR_miss_mean_zs(:,n_order_hDCnD)',lim_clr_hit)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Miss'})
        colormap(ax2, 'hot')
        
        ax3 = subplot(133);
        imagesc(bin_edges_orig,1:n_neu,FR_diff_hitMiss(:,n_order_hDCnD)',lim_clr_diff)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Diff'})
        colormap(ax3,  bwr_cmap)
        
        set(gcf,'position',[1972 212 1870 644],'color','w')
        saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
        close(fig);
        
        % --- SELECTION: success/fail ---
        FR_suc_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx == 1), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_succ, 'UniformOutput', false)');
        FR_fail_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx == 0), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_succ, 'UniformOutput', false)');
        FR_diff_SucFail = FR_fail_mean_zs-FR_suc_mean_zs;
        
        % --- FIGURE SUCCESS/FAIL ---
        fig_name = 'tile_succ_fail';
        fig = figure('Visible', 'off');
        ax1 = subplot(131);
        imagesc(bin_edges_orig,1:n_neu,FR_suc_mean_zs(:,n_order_hDCnD)',lim_clr);
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); ylabel('cells'); title({[curr_type ' Reach'];'Success'})
        colormap(ax1, 'hot')
        
        ax2 = subplot(132);
        imagesc(bin_edges_orig,1:n_neu,FR_fail_mean_zs(:,n_order_hDCnD)',lim_clr)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Fail'})
        colormap(ax2, 'hot')
        
        ax3 = subplot(133);
        imagesc(bin_edges_orig,1:n_neu,FR_diff_SucFail(:,n_order_hDCnD)',lim_clr_diff)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Diff'})
        colormap(ax3,  bwr_cmap)
        
        set(gcf,'position',[1972 212 1870 644],'color','w')
        saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
        close(fig);
        
        % --- SELECTION: Categories ---
        FR_r1_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx == 1), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_cat, 'UniformOutput', false)');
        FR_r2_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx == 2), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_cat, 'UniformOutput', false)');
        FR_r3_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx == 3), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_cat, 'UniformOutput', false)');
        FR_r4_mean_zs = cell2mat(cellfun(@(fr, idx) mean(fr(:, idx == 4), 2, "omitnan"), ...
            zFR_reach, N.idx_reach_cat, 'UniformOutput', false)');
        
        % --- FIGURE REACH TYPES ---
        fig_name = 'tile_reach_types';
        lim_clr_types = [-.7 1.1];
        fig = figure('Visible', 'off');
        subplot(141)
        imagesc(bin_edges_orig,1:n_neu,FR_r1_mean_zs(:,n_order_hDCnD)',lim_clr_types);
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); ylabel('cells'); title({[curr_type ' Reach'];'Full'})
        
        subplot(142)
        imagesc(bin_edges_orig,1:n_neu,FR_r2_mean_zs(:,n_order_hDCnD)',lim_clr_types)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Lift'})
        
        subplot(143)
        imagesc(bin_edges_orig,1:n_neu,FR_r3_mean_zs(:,n_order_hDCnD)',lim_clr_types)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Drink'})
        
        subplot(144)
        imagesc(bin_edges_orig,1:n_neu,FR_r4_mean_zs(:,n_order_hDCnD)',lim_clr_types)
        if flip_ud==true, axis xy; end
        xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
        set(gca,axeOpt{:},'TickDir','out')
        xlabel('time (s)'); title({[curr_type ' Reach'];'Groom'})
        
        colormap("hot")
        set(gcf,'position',[1972 212 1870 644],'color','w')
        saveas(gcf,fullfile(save_out_reach,[fig_name,'_w2.png']),'png');
        close(fig);
    end
end

fprintf('All cell types processed.\n');

%% Helper Functions

function zFR = concat_zscore(FR_all)
    zFR = cell(size(FR_all)); % Initialize cell array
    for i = 1:length(FR_all) 
        data = FR_all{i}; % Use {} to extract numeric data from the cell
        data_flat = data(:);
        zs_data_flat = zscore(data_flat);
        % Use {i} to store the numeric matrix back into the cell array
        zFR{i} = reshape(zs_data_flat, size(data));
    end
end