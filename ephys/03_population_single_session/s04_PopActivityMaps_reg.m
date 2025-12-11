%% POPULATION ANALYSIS
% heatmaps of all cells of each region
clear; close all; clc

%% Load
rootdir ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
group = '20230801_ChocolateGroup';
datapath = fullfile(rootdir,"ephys_and_behavior","mat_files",group);
reg = 'CT';

load_file = sprintf('%s_%s_%s','all',reg,'neurons');
S=load(fullfile(datapath,'eg_neurons.mat'),load_file,'eg_neu_FR_params');

% Save apth
save_out = fullfile(rootdir,"ephys_and_behavior","out_files",group,'group','activityMaps',reg);
if ~exist(save_out,"dir"), mkdir(save_out); end

%% Neurons table and properties
N = struct2table(S.(load_file));
n_neu = size(N,1);

N_pp = N(~strcmp(N.mouse, '3_Toblerone'), :);
n_neu_pp = size(N_pp,1);

axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};

%bin_width = eg_neu_FR_params.bin_width;
bin_edges_orig = S.eg_neu_FR_params.bin_edges;
nr_bins_original = S.eg_neu_FR_params.nr_bins;


%% INIT
% save_out
save_out_init = fullfile(save_out,'init');
if ~exist(save_out_init,"dir"), mkdir(save_out_init); end

%% PCA
FR_init_mean_tmp = cellfun(@(x) mean(x, 2), N_pp.FR_init, 'UniformOutput', false);
FR_init_mean = cat(2, FR_init_mean_tmp{:});   % concatenate into one matrix

[zs_FR_init, FR_init_mu, FR_init_sigma] = zscore(FR_init_mean);
[coeff_init, score_init, latent_init, tsquare_init, explained_init, mus_init] = pca(zs_FR_init);

% Order of cells
r_coeff_pp = sqrt(coeff_init(:,1).^2 + coeff_init(:,2).^2);
theta_coef_pp = atan2(coeff_init(:,2), coeff_init(:,1));
[theta_coeffPP_sorted, neu_order_pp] = sort(theta_coef_pp);
cmap = parula(length(neu_order_pp));

% Find the starting point 
[~,break_rad] = max(diff(theta_coeffPP_sorted));
figure, plot(diff(theta_coeffPP_sorted),'*')
%break_rad = 68;
n_order_pp = [neu_order_pp(break_rad+1:end);neu_order_pp(1:break_rad)];
%n_order_pp = neu_order_pp;

figure
subplot(121)
plot(cumsum(explained_init),'-ko');
xlabel('PCs'); ylabel('explained')
subplot(122)
plot(cumsum(score_init(:,1:2))); 
xlabel('time'); ylabel('score')
legend('PC1','PC2')
set(gcf,'Position',[ 2762         134        1027         423])
saveas(gcf,fullfile(save_out_init,'PCA_pushPull.png'),'png');

figure
subplot(121)
scatter(coeff_init(:,1),coeff_init(:,2),50,cmap,'filled');
xlabel('PC1 coefficients'); ylabel('PC2 coefficients')
set(gca,axeOpt{:})
title('2497         264        1109         443')
axis square
hold on
for i=1:n_neu_pp
    subplot(122)
    polarplot(theta_coeffPP_sorted(i), r_coeff_pp(neu_order_pp(i)),...
        'o','MarkerSize', 8, 'MarkerFaceColor',cmap(i,:),'markerEdgeColor','w'); hold on
end
colorbar
title('Points in Polar Coordinates');
title('Angular position of PC1 vs PC2 coeffs')
set(gcf,'Position',[2762         134        1027         734])
saveas(gcf,fullfile(save_out_init,'PC1coeffs_vs_PC2coeffs_pushPull.png'),'png');


%% NUERONS TILING!!!

% SELECTION: PUSH / PULL
% Push 
FR_push_all = cellfun(@(fr, idx) fr(:, idx == 0), ...
    N_pp.FR_init, N_pp.idx_init_vPP_invPP, ...
    'UniformOutput', false);
FR_push_mean_zs = concat_zscore_and_mean(FR_push_all);

% Pull
FR_pull_all = cellfun(@(fr, idx) fr(:, idx == 1), ...
    N_pp.FR_init, N_pp.idx_init_vPP_invPP, ...
    'UniformOutput', false);
FR_pull_mean_zs = concat_zscore_and_mean(FR_pull_all);

FR_diff_pp = FR_pull_mean_zs-FR_push_mean_zs;

%% PUSH PULL FIGURE
fig_name = 'tile_push_pull';
lim_clr = [-.7 1.2];
lim_clr_diff = [-1 1];
figure
ax1 = subplot(131);
imagesc(bin_edges_orig,1:n_neu_pp,FR_push_mean_zs(:,n_order_pp)',lim_clr);
%if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Push'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax1, 'hot')

ax2 = subplot(132);
imagesc(bin_edges_orig,1:n_neu_pp,FR_pull_mean_zs(:,n_order_pp)',lim_clr)
%if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Pull'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax2, 'hot')

ax3 = subplot(133);
imagesc(bin_edges_orig,1:n_neu_pp,FR_diff_pp(:,n_order_pp)',lim_clr_diff)
%if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Diff'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')
colormap(ax3, bwr_map)

set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_init,[fig_name,'.png']),'png');
%print(gcf, fullfile(save_out_init,[fig_name,'.pdf']), '-dpdf', '-painters');

%% SELECTION: VALID / INVALD
% Valid 
FR_val_all = cellfun(@(fr, idx) fr(:, idx == 0 | idx == 1), ...
    N_pp.FR_init, N_pp.idx_init_vPP_invPP, ...
    'UniformOutput', false);
FR_val_mean_zs = concat_zscore_and_mean(FR_val_all);

% Invalid
FR_inval_all = cellfun(@(fr, idx) fr(:, idx == 2 | idx == 3), ...
    N_pp.FR_init, N_pp.idx_init_vPP_invPP, ...
    'UniformOutput', false);
FR_inval_mean_zs = concat_zscore_and_mean(FR_inval_all);

FR_diff_val = FR_inval_mean_zs-FR_val_mean_zs;


%% VALID INVALID FIGURE
fig_name = 'tile_valid_invalid';
%lim_clr = [-.7 1.1];
lim_clr_diff = [-1 1];
figure
ax1 = subplot(131);
imagesc(bin_edges_orig,1:n_neu_pp,FR_val_mean_zs(:,n_order_pp)',lim_clr);
%if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Val'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax1, 'hot')

ax2 = subplot(132);
imagesc(bin_edges_orig,1:n_neu_pp,FR_inval_mean_zs(:,n_order_pp)',lim_clr)
if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Inval'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax2, 'hot')

ax3 = subplot(133);
imagesc(bin_edges_orig,1:n_neu_pp,FR_diff_val(:,n_order_pp)',lim_clr_diff)
if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Diff'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')
colormap(ax3, bwr_map)

set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_init,[fig_name,'.png']),'png');
%print(gcf, fullfile(save_out_init,[fig_name,'.pdf']), '-dpdf', '-painters');

%% SELECTION: Dominant, center, non-dominant paw
% Dominant paw side 
FR_domi_all = cellfun(@(fr, idx) fr(:, idx == 1), ...
    N_pp.FR_init, N_pp.DCnD_init, ...
    'UniformOutput', false);
FR_domi_mean_zs = concat_zscore_and_mean(FR_domi_all);

% Center
FR_Ci_all = cellfun(@(fr, idx) fr(:, idx == 2), ...
    N_pp.FR_init, N_pp.DCnD_init, ...
    'UniformOutput', false);
FR_Ci_mean_zs = concat_zscore_and_mean(FR_Ci_all);

% Non-dominant paw side
FR_ndomi_all = cellfun(@(fr, idx) fr(:, idx == 3), ...
    N_pp.FR_init, N_pp.DCnD_init, ...
    'UniformOutput', false);
FR_ndomi_mean_zs = concat_zscore_and_mean(FR_ndomi_all);

FR_diff_DnD_init = FR_ndomi_mean_zs-FR_domi_mean_zs;


%% WATER SIDE FIGURE
fig_name = 'tile_dom_C_nonDom_init';
%lim_clr = [-0.8 1.5];
figure
ax1 = subplot(131);
imagesc(bin_edges_orig,1:n_neu_pp,FR_domi_mean_zs(:,n_order_pp)',lim_clr);
if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Dom side'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax1, 'hot')

ax2 = subplot(132);
imagesc(bin_edges_orig,1:n_neu_pp,FR_Ci_mean_zs(:,n_order_pp)',lim_clr)
if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Center'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax2, 'hot')

ax3 = subplot(133);
imagesc(bin_edges_orig,1:n_neu_pp,FR_ndomi_mean_zs(:,n_order_pp)',lim_clr)
if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Non-Dom side'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')
colormap(ax3, 'hot')


set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_init,[fig_name,'.png']),'png');
%print(gcf, fullfile(save_out_init,[fig_name,'.pdf']), '-dpdf', '-painters');


%% WATER SIDE FIGURE
fig_name = 'tile_dom_nonDom_diff_init';
lim_clr = [-0.7 1.2];
figure
ax1 = subplot(131);
imagesc(bin_edges_orig,1:n_neu_pp,FR_domi_mean_zs(:,n_order_pp)',lim_clr);
if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Dom side'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax1, 'hot')

ax2 = subplot(132);
imagesc(bin_edges_orig,1:n_neu_pp,FR_ndomi_mean_zs(:,n_order_pp)',lim_clr)
if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Non-Dom side'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax2, 'hot')

ax3 = subplot(133);
imagesc(bin_edges_orig,1:n_neu_pp,FR_diff_DnD_init(:,n_order_pp)',lim_clr_diff)
if ~strcmp(reg,'BG'), axis xy; end
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from trial init (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Diff'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')
colormap(ax3, bwr_map)

set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_init,[fig_name,'.png']),'png');
%print(gcf, fullfile(save_out_init,[fig_name,'.pdf']), '-dpdf', '-painters');


%% -------------------------------------------------------------------------
%% REACH
% save_out
save_out_reach = fullfile(save_out,'reach');
if ~exist(save_out_reach,"dir"), mkdir(save_out_reach); end

% PCA
FR_reach_mean_tmp = cellfun(@(x) mean(x, 2), N.FR_reach, 'UniformOutput', false);
FR_reach_mean = cat(2, FR_reach_mean_tmp{:});   % concatenate into one matrix

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
%n_order_hDCnD = neu_order_hDCnD;

figure
subplot(121)
plot(cumsum(explained_reach),'-ko');
xlabel('PCs'); ylabel('explained')
subplot(122)
plot(cumsum(score_reach(:,1:2))); 
xlabel('time'); ylabel('score')
legend('PC1','PC2')
set(gcf,'Position',[ 2762         134        1027         423])
saveas(gcf,fullfile(save_out_reach,'PCA_reach.png'),'png');

figure
subplot(121)
scatter(coeff_reach(:,1),coeff_reach(:,2),50,cmap_hDCnD,'filled');
xlabel('PC1 coefficients'); ylabel('PC2 coefficients')
set(gca,axeOpt{:})
title('2497         264        1109         443')
axis square
hold on
for i=1:n_neu
    subplot(122)
    polarplot(theta_coef_hDCnD(i), r_coeff_hDCnD(n_order_hDCnD(i)),...
        'o','MarkerSize', 8, 'MarkerFaceColor',cmap_hDCnD(i,:),'markerEdgeColor','w'); hold on
end
colorbar
title('Points in Polar Coordinates');
title('Angular position of PC1 vs PC2 coeffs')
set(gcf,'Position',[2762         134        1027         734])
saveas(gcf,fullfile(save_out_init,'PC1coeffs_vs_PC2coeffs_reach.png'),'png');

%% -------------------------------------------------------------------
% Neurons tiling!!

% SELECTION: Dom/C/non-Dom
% Dom 
FR_D_all = cellfun(@(fr, idx1, idx2) fr(:, idx1 == 1 & (idx2==1)), ...
    N.FR_reach, N.idx_reach_hDCnD, N.idx_reach_cat, ...
    'UniformOutput', false);
FR_D_mean_zs = concat_zscore_and_mean(FR_D_all);

% Center
FR_C_all = cellfun(@(fr, idx1, idx2) fr(:, idx1 == 2 & (idx2==1)), ...
    N.FR_reach, N.idx_reach_hDCnD, N.idx_reach_cat, ...
    'UniformOutput', false);
FR_C_mean_zs = concat_zscore_and_mean(FR_C_all);

% non-Dom
FR_nD_all = cellfun(@(fr, idx1, idx2) fr(:, idx1 == 3 & (idx2==1)), ...
    N.FR_reach, N.idx_reach_hDCnD, N.idx_reach_cat, ...
    'UniformOutput', false);
FR_nD_mean_zs = concat_zscore_and_mean(FR_nD_all);


FR_diff_DnD = FR_nD_mean_zs-FR_D_mean_zs;

%% Dom/C/nonDom FIGURE
n_neu_all = n_neu;
fig_name = 'tile_reach_DCnD';
lim_clr = [-.8 1.5];
%bins_sel=1001:1751;
%bins_sel=1251:1751;
%bins_sel=1:length(bin_edges);
figure
ax1 = subplot(131);
imagesc(bin_edges_orig,1:n_neu_all,FR_D_mean_zs(:,n_order_hDCnD)',lim_clr);
axis xy; 
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Dom side'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax1, 'hot')

ax2 = subplot(132);
imagesc(bin_edges_orig,1:n_neu_all,FR_C_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Center'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax2, 'hot')

ax3 = subplot(133);
imagesc(bin_edges_orig,1:n_neu_all,FR_nD_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Non-Dom side'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')
colormap(ax3,  'hot')

set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
%print(gcf, fullfile(save_out_reach,[fig_name,'.pdf']), '-dpdf', '-painters');


%% Dom/nonDom/diff FIGURE 
fig_name = 'tile_reach_DnD_wDiff';
%lim_clr = [-.9 1.5];
%bins_sel=1001:1751;
%bins_sel=1251:1751;
%bins_sel=1:length(bin_edges);
figure
ax1 = subplot(131);
imagesc(bin_edges_orig,1:n_neu_all,FR_D_mean_zs(:,n_order_hDCnD)',lim_clr);
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Dom side'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax1, 'hot')

ax2 = subplot(132);
imagesc(bin_edges_orig,1:n_neu_all,FR_C_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Non-Dom side'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax2, 'hot')

ax3 = subplot(133);
imagesc(bin_edges_orig,1:n_neu_all,FR_diff_DnD(:,n_order_hDCnD)',lim_clr_diff)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Diff'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')
colormap(ax3,  bwr_map)

set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
%print(gcf, fullfile(save_out_reach,[fig_name,'.pdf']), '-dpdf', '-painters');

%% SELECTION: L/C/R
% Dom 
FR_L_all = cellfun(@(fr, idx1, idx2) fr(:, idx1 == 1 & idx2==1), ...
    N.FR_reach, N.idx_reach_LCR, N.idx_reach_cat, ...
    'UniformOutput', false);
FR_L_mean_zs = concat_zscore_and_mean(FR_L_all);

% non-Dom
FR_R_all = cellfun(@(fr, idx1, idx2) fr(:, idx1 == 3 & idx2==1), ...
    N.FR_reach, N.idx_reach_LCR, N.idx_reach_cat, ...
    'UniformOutput', false);
FR_R_mean_zs = concat_zscore_and_mean(FR_R_all);



%% L/C/R FIGURE
fig_name = 'tile_reach_LCR';
%lim_clr = [-.5 .8];
%bins_sel=1001:1751;
%bins_sel=1251:1751;
%bins_sel=1:length(bin_edges);
figure
ax1 = subplot(131);
imagesc(bin_edges_orig,1:n_neu_all,FR_L_mean_zs(:,n_order_hDCnD)',lim_clr);
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Left'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax1, 'hot')

ax2 = subplot(132);
imagesc(bin_edges_orig,1:n_neu_all,FR_C_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Center'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax2, 'hot')

ax3 = subplot(133);
imagesc(bin_edges_orig,1:n_neu_all,FR_R_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Right'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')
colormap(ax3,  'hot')

set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
%print(gcf, fullfile(save_out_reach,[fig_name,'.pdf']), '-dpdf', '-painters');

%% SELECTION: Hit/miss
% Hit 
FR_hit_all = cellfun(@(fr, idx) fr(:, idx == 1), ...
    N.FR_reach, N.idx_reach_hit, ...
    'UniformOutput', false);
FR_hit_mean_zs = concat_zscore_and_mean(FR_hit_all);

% Miss
FR_miss_all = cellfun(@(fr, idx) fr(:, idx == 0), ...
    N.FR_reach, N.idx_reach_hit, ...
    'UniformOutput', false);
FR_miss_mean_zs = concat_zscore_and_mean(FR_miss_all);

FR_diff_hitMiss = FR_miss_mean_zs-FR_hit_mean_zs;

%% FIGURE HIT MISS
fig_name = 'tile_hit_miss';
%lim_clr = [-.8 1.5];
%bins_sel=1001:1751;
%bins_sel=1251:1751;
%bins_sel=1:length(bin_edges);
figure
ax1 = subplot(131);
imagesc(bin_edges_orig,1:n_neu_all,FR_hit_mean_zs(:,n_order_hDCnD)',lim_clr);
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Hit'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax1, 'hot')

ax2 = subplot(132);
imagesc(bin_edges_orig,1:n_neu_all,FR_miss_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Miss'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax2, 'hot')

ax3 = subplot(133);
imagesc(bin_edges_orig,1:n_neu_all,FR_diff_hitMiss(:,n_order_hDCnD)',lim_clr_diff)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Diff'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')
colormap(ax3,  bwr_map)

set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
%print(gcf, fullfile(save_out_reach,[fig_name,'.pdf']), '-dpdf', '-painters');


%% SELECTION: success/fail
% Success 
FR_suc_all = cellfun(@(fr, idx) fr(:, idx == 1), ...
    N.FR_reach, N.idx_reach_succ, ...
    'UniformOutput', false);
FR_suc_mean_zs = concat_zscore_and_mean(FR_suc_all);

% Fail
FR_fail_all = cellfun(@(fr, idx) fr(:, idx == 0), ...
    N.FR_reach, N.idx_reach_succ, ...
    'UniformOutput', false);
FR_fail_mean_zs = concat_zscore_and_mean(FR_fail_all);


FR_diff_SucFail = FR_fail_mean_zs-FR_suc_mean_zs;


%% FIGURE SUCCESS/FAIL
fig_name = 'tile_succ_fail';
%lim_clr = [-.8 1.1];
%bins_sel=1001:1751;
%bins_sel=1251:1751;
%bins_sel=1:length(bin_edges);
figure
ax1 = subplot(131);
imagesc(bin_edges_orig,1:n_neu_all,FR_suc_mean_zs(:,n_order_hDCnD)',lim_clr);
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Success'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax1, 'hot')

ax2 = subplot(132);
imagesc(bin_edges_orig,1:n_neu_all,FR_fail_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Fail'})
c=colorbar;
ylabel(c,'z-scored FR')
colormap(ax2, 'hot')

ax3 = subplot(133);
imagesc(bin_edges_orig,1:n_neu_all,FR_diff_SucFail(:,n_order_hDCnD)',lim_clr_diff)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Diff'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')
colormap(ax3,  bwr_map)

set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_reach,[fig_name,'.png']),'png');
%print(gcf, fullfile(save_out_reach,[fig_name,'.pdf']), '-dpdf', '-painters');

%% SELECTION: Categories
% full reach 
FR_r1_all = cellfun(@(fr, idx) fr(:, idx == 1), ...
    N.FR_reach, N.idx_reach_cat, ...
    'UniformOutput', false);
FR_r1_mean_zs = concat_zscore_and_mean(FR_r1_all);

FR_r2_all = cellfun(@(fr, idx) fr(:, idx == 2), ...
    N.FR_reach, N.idx_reach_cat, ...
    'UniformOutput', false);
FR_r2_mean_zs = concat_zscore_and_mean(FR_r2_all);

FR_r3_all = cellfun(@(fr, idx) fr(:, idx == 3), ...
    N.FR_reach, N.idx_reach_cat, ...
    'UniformOutput', false);
FR_r3_mean_zs = concat_zscore_and_mean(FR_r3_all);

FR_r4_all = cellfun(@(fr, idx) fr(:, idx == 4), ...
    N.FR_reach, N.idx_reach_cat, ...
    'UniformOutput', false);
FR_r4_mean_zs = concat_zscore_and_mean(FR_r4_all);

%% FIGURE REACH TYPES
fig_name = 'tile_reach_types';
lim_clr = [-.7 1.1];
%bins_sel=1001:1751;
%bins_sel=1251:1751;
%bins_sel=1:length(bin_edges);
figure
subplot(141)
imagesc(bin_edges_orig,1:n_neu_all,FR_r1_mean_zs(:,n_order_hDCnD)',lim_clr);
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Full reaches'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

subplot(142)
imagesc(bin_edges_orig,1:n_neu_all,FR_r2_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Lifted paw reaches'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

subplot(143)
imagesc(bin_edges_orig,1:n_neu_all,FR_r3_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Drink / Lick'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

subplot(144)
imagesc(bin_edges_orig,1:n_neu_all,FR_r4_mean_zs(:,n_order_hDCnD)',lim_clr)
axis xy;
xline(0,'--','linewidth',1.5,'color',[1 1 1 .5])
set(gca,axeOpt{:},'TickDir','out')
xlabel('time from reach endpoint (s)'); ylabel(sprintf('%s%s%s','cells in ',reg, ' (PC1-2 sorted)'));
title({sprintf('%s%s','Mean FR in ',reg);'Grooming'})
c=colorbar;
ylabel(c,'z-scored FR (sp/s)')

colormap("hot")

set(gcf,'position',[1972 212 1870 644],'color','w')
saveas(gcf,fullfile(save_out_reach,[fig_name,'_w2.png']),'png');
%print(gcf, fullfile(save_out_reach,[fig_name,'_w2.pdf']), '-dpdf', '-painters');



