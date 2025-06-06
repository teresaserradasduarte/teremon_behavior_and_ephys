%%
mouse = '5_FerreroRocher';
sess = 'R1';
load(fullfile(rootdir,"behavior_data/","analyzed_data/mat_files/","20230511_ChocolateGroup/","headfixed_dynamicTarget/",mouse,sess,"session_reaching_data_paw.mat"),'reaches')
load(fullfile(rootdir,"behavior_data/","analyzed_data/mat_files/","20230511_ChocolateGroup/","headfixed_dynamicTarget/",mouse,sess,"behavior_session.mat"))
load(fullfile(rootdir,"ephys_and_behavior/","mat_files/","20230801_ChocolateGroup/",mouse,sess,"CB/","eg_neurons.mat"));

%%
r1_idx = neu_strct(1).idx_reach_cat==1;
reach_r1 = neu_strct(1).reach_px(:,:,r1_idx);
tm_reach = reaches.tm_w;
win_size_endpoint = .05;
start_win_bhv = -.01;
win_stop = start_win_bhv+win_size_endpoint;
[~,r_start] = min(abs(tm_reach-start_win_bhv));
[~,r_stop] = min(abs(tm_reach-win_stop));

%%
water_loc = neu_strct(1).idx_reach_LCR;
water_loc_r1 = water_loc(r1_idx);
n_r1 = length(water_loc_r1);
bin_edges = eg_neu_FR_params.bin_edges;

sel_win = [-1 1];
[~,new_bin_edges_start] = min(abs(bin_edges-sel_win(1)));
[~,new_bin_edges_stop] = min(abs(bin_edges-sel_win(2)));
FR_edges_range = new_bin_edges_start:new_bin_edges_stop;

colors_lcr=[behavior.colors.left_color;behavior.colors.center_color;behavior.colors.right_color];
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};

%% Data to MI
n = 16;

win_start_n_t = -.2;
win_stop_n_t = win_start_n_t + .1;
[~,r_start_n] = min(abs(bin_edges-win_start_n_t));
[~,r_stop_n] = min(abs(bin_edges-win_stop_n_t));

endpoint_xyz = squeeze(mean(reach_r1(r_start:r_stop,:,:),1));
neu_FR_win = squeeze(mean(neu_strct(n).FR_reach(r_start_n:r_stop_n,r1_idx),1));

%% MI landscape
N_tmp = struct2table(N_CB_all);
% exclude the only session with just neurons from the cerebellum
exclude_idx = strcmp(N_tmp.mouse, '1_CoteDor') & strcmp(N_tmp.sess, 'R1');
N = N_tmp(~exclude_idx, :);

idx_sess = strcmp(N.mouse, mouse) & strcmp(N.sess, sess);
MI_sub = MI_struct(idx_sess);

MI_land_n = MI_sub(n).MI_landscape;

%%
sz=15;
%fr_lim = [0 30];
figure
ff = tiledlayout(3,3);
%     title(ff,sprintf('%s%s%s%s','Full reaches | mouse: ',...
%         mouse, ' | session: ',sess),...
%         'Interpreter','none');

%xyz_label = ['x';'y';'z'];
ax1=nexttile;
ax1.Layout.TileSpan = [2,1];
imagesc(tm_reach,1:n_r1, squeeze(reach_r1(:,nxt,:))');
xline(0,'--','Color',[1 1 1 .5],'LineWidth',2);
xline(tm_reach(r_start),'-','Color',[.5 .5 .5 .5],'LineWidth',1.5);
xline(tm_reach(r_stop),'-','Color',[.5 .5 .5 .5],'LineWidth',1.5);
xlabel('time (s)'); ylabel('reach idx');
set(gca,axeOpt{:});
%title(xyz_label(nxt))
colormap(ax1, 'bone')
c1=colorbar; ylabel(c1,'ML (px)'); hold on
scatter(ones(size(1:n_r1)).*.125,1:n_r1, 40 ,colors_lcr(water_loc_r1, :), 'filled','Marker','square'); hold off


nexttile(7)
scatter(1:n_r1,endpoint_xyz(2,:), sz, colors_lcr(water_loc_r1, :), 'filled'); hold off
ylabel('mean position in ML (px)');
xlabel('reach idx');
xlim([1 n_r1])
set(gca,axeOpt{:});

ax2=nexttile(2);
ax2.Layout.TileSpan = [2,1];
imagesc(tm_reach,1:n_r1, squeeze(reach_r1(:,nxt,:))');
imagesc(bin_edges(FR_edges_range),1:n_r1, neu_strct(n).FR_reach(:,r1_idx)',[0 25]);
xline(0,'--','Color',[1 1 1 .5],'LineWidth',2);
xline(bin_edges(r_start_n),'-','Color',[.8 .8 .8 .5],'LineWidth',1.5);
xline(bin_edges(r_stop_n),'-','Color',[.8 .8 .8 .5],'LineWidth',1.5);
xlabel('time (s)'); ylabel('reach idx');
colormap(ax2, 'bone')
hold on
scatter(ones(size(1:n_r1)).*.99,1:n_r1, 40 ,colors_lcr(water_loc_r1, :), 'filled','Marker','square'); hold off

set(gca,axeOpt{:});
%title(xyz_label(nxt))
c2=colorbar; ylabel(c2,'instantaneous FR (sp/s)');

nexttile(8)
scatter(1:n_r1,neu_FR_win, sz ,colors_lcr(water_loc_r1, :), 'filled'); hold off
%scatter(1:n_r1,neu_FR_win,sz,'filled','MarkerFaceColor',[.3 .3 .3]); hold on
%scatter(1:n_r1,ones(size(1:n_r1)).*32, sz ,colors_lcr(water_loc_r1, :), 'filled'); hold off
ylabel('mean FR (sp(s)');
xlabel('reach idx');
xlim([1 n_r1])
ylim([-2 32]);
set(gca,axeOpt{:});


ax3=nexttile(3);
ax3.Layout.TileSpan = [3,1];
imagesc(x_im,y_im,MI_land_n);
axis xy;        % Invert the y-axis direction
hold on
plot(win_stop_n_t,win_start_n_t,'ko','LineWidth',1,'Marker','square')

alpha_data = ~isnan(MI_land_n);
set(gca, 'ALim', [0 1]);
hImg = findobj(gca, 'Type', 'image');  % get image object only
set(hImg, 'AlphaData', alpha_data);
caxis([min(MI_land_n(:), [], 'omitnan') max(MI_land_n(:), [], 'omitnan')]);
colormap(ax3, 'parula')
%set(gca, 'YTickDir', 'reverse')
c=colorbar;
ylabel(c,'Bits');
xline(0,'--','color',[1 1 1],'LineWidth',1.5)
yline(0,'--','color',[1 1 1],'LineWidth',1.5)
ylabel('start time (s)'); xlabel('stop time (s)')
set(gca,axeOpt{:})
axis square


%end
set(gcf,'Position',[2102         231        1525         505],'color','w');
%saveas(gcf,strcat(save_out_bhv,filesep,fig_name,'.png'),'png');
%print(gcf,strcat(save_out_bhv,filesep,fig_name,'.pdf'), '-dpdf', '-painters');

%
fig_name = 'bhv_and_neuron_v2';

save_out_bhv = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_and_behavior\out_files\20230801_ChocolateGroup\5_FerreroRocher\R6\bhv';
saveas(ff,strcat(save_out_bhv,filesep,fig_name,'.png'),'png');
print(gcf,strcat(save_out_bhv,filesep,fig_name,'.pdf'), '-dpdf', '-painters');
