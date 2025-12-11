% SVM all regions and sessions figure
clear; close all; clc

%% Figure defs
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};
clr_center = [72 60 50]./256;
clr_dom = [48 131 220]./256;
clr_nondom = [158 216 219]./256;
clr_init = [63 130 109]./256;
clrs = [clr_init;clr_dom;clr_center;clr_nondom];

%% Paths
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_and_behavior\';
save_out_dir = fullfile(rootdir,"out_files/20230801_ChocolateGroup/group/","SVM");
if ~exist(save_out_dir,"dir"), mkdir(save_out_dir); end
% Load CB
reg = 'BG';
load(fullfile(rootdir,"mat_files/20230801_ChocolateGroup/",['SVM_sess_',reg,'.mat']));
svm_sess = SVM_sess;
clear SVM_sess


%% Select sessions CB
% For later - check ik s1 trials and s2 trials are comparable
% now, simply select which to include manually
svm_tmp = struct2table(svm_sess);
% exclude the only session with just neurons from the cerebellum
exclude_idx = strcmp(svm_tmp.mouse, '1_CoteDor') & strcmp(svm_tmp.sess, 'R1') | ...
    strcmp(svm_tmp.mouse, '3_Toblerone') | ...
    strcmp(svm_tmp.mouse, '5_FerreroRocher') & strcmp(svm_tmp.sess, 'R6');
svm_tbl = svm_tmp(~exclude_idx, :);

%% Get mean accuracy
bac_init_tmp = cellfun(@(x) reshape(x, 1, size(x,1), size(x,2)), svm_tbl.bac_init, 'UniformOutput', false);
bac_init = cat(1, bac_init_tmp{:});  % Result: [7 x 197 x 4]

bac_reach_tmp = cellfun(@(x) reshape(x, 1, size(x,1), size(x,2)), svm_tbl.bac_reach, 'UniformOutput', false);
bac_reach = cat(1, bac_reach_tmp{:});  % Result: [7 x 197 x 4]


%% Check one mouse
bin_edges = svm_sess(1).bin_edges;
bin_mean=bin_edges(:,1)+diff(bin_edges(1:2));
figure, plot(bin_mean,squeeze(bac_init(7,:,:)));

%% Mean across mice
bac_init_mean = squeeze(mean(bac_init,1));
bac_reach_mean = squeeze(mean(bac_reach,1));
win_int = [bin_edges(1,1) bin_edges(end,2)];

%% Plot
figure()
tt = tiledlayout(1,2);

nexttile
for plt=1:4
    plot(bin_mean,bac_init_mean(:,plt),'color',clrs(plt,:),'linewidth',2); hold on
end
xlabel('time from trial init (s)'); ylabel('accuracy');
set(gca,axeOpt{:})
yline(.5,'--','color',[.7 .7 .7 .7],'LineWidth',2)
xline(0,'-.','color',[.8 .8 .8 .8])
xlim(win_int); ylim([.45 1])

nexttile
for plt=1:4
    plot(bin_mean,bac_reach_mean(:,plt),'color',clrs(plt,:),'linewidth',2); hold on
end
xlabel('time from reach endpoint (s)'); ylabel('accuracy');
set(gca,axeOpt{:})
yline(.5,'--','color',[.7 .7 .7 .7],'LineWidth',2)
xline(0,'-.','color',[.8 .8 .8 .8])
xlim(win_int); ylim([.45 1])

legend('push / pull','dominat paw side','center', 'non-dominant paw','box','off')

set(gcf,'Position',[680         585        1226         393],'Color','w');
saveas(gcf,strcat(save_out_dir,filesep,'svm_accuracy_',reg,'.png'),'png');
print(gcf,strcat(save_out_dir,filesep,'svm_accuracy_',reg,'.pdf'), '-dpdf', '-painters');




