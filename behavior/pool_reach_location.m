%% BEHAVIOR GROPU POOL - dom/center/non-dom
clear; close all; clc

%% Load
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
group = '20230801_ChocolateGroup';
group_behav = '20230511_ChocolateGroup';
setup = 'headfixed_dynamicTarget';
animals = {...
    'CoteDor',...
    'Lindt',...
    'Toblerone',...
    'Milka',...
    'FerreroRocher'};

all_sess = {['R1';'R6'],...
    ['R1';'R2';'R4'],...
    ['R3';'R6'],...
    ['R4'],...
    ['R1';'R4';'R6']};

% Initialize
load(fullfile(rootdir,"behavior_data/","analyzed_data/mat_files/","20230511_ChocolateGroup/","headfixed_dynamicTarget/","1_CoteDor/","R6/","behavior_session.mat"))


%% Create path to save
%save_mat = fullfile(rootdir,"behavior_data","analyzed_data","mat_files",group_behav,setup,'group');
save_out = fullfile(rootdir,"behavior_data","analyzed_data","output_files",group_behav,setup,'behav_group','reach_location');
if ~exist(save_out,"dir"), mkdir(save_out); end
%if ~exist(save_mat,"dir"), mkdir(save_mat); end

mat_dir_r = fullfile(rootdir,"behavior_data","analyzed_data","mat_files",group_behav,setup);
mat_dir_b = fullfile(rootdir,"ephys_and_behavior","mat_files",group);

%% Initialize variables 
n_max = 500;
n_mice = length(animals);
n_sess = 3;
reaches1_px_loc = nan(201,3,500,3,n_sess,n_mice);
n_reaches_loc = nan(3,n_sess,n_mice);
n_reaches_loc_LCR = nan(3,n_sess,n_mice);
n_reaches_cat12_loc = nan(3,2,n_sess,n_mice);
n_suc_loc = nan(3,n_sess,n_mice);
n_hit_loc = nan(3,n_sess,n_mice);
duration_f_loc = nan(n_max,3,n_sess,n_mice);

%% SVM - Lood through mice and session
for m=1:length(animals)
    animal_idx = m;
    mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});
    fprintf('%s%s%s\n','Looping through ',mouse,':') 

    for s=1:size(all_sess{m},1)
        sess = all_sess{m}(s,:);
        fprintf('%s%s%s%s%s\n','Running mouse ',mouse,', session ',sess,':')

            % Load mat files
            load(fullfile(mat_dir_r,mouse,sess,"session_reaching_data_paw.mat"),'reaches');
            load(fullfile(mat_dir_b,mouse,sess,"behavior_fundamentals.mat"),'bhv');   

            r1_D_idx = find(bhv.DCnD_inVec_idx==1 & bhv.cat_reach_inVec==1);
            r1_C_idx = find(bhv.DCnD_inVec_idx==2 & bhv.cat_reach_inVec==1);
            r1_nD_idx = find(bhv.DCnD_inVec_idx==3 & bhv.cat_reach_inVec==1);

            reaches1_px_loc(:,:,1:length(r1_D_idx),s,m) = bhv.reaches_inVec_px(:,:,r1_D_idx);
            reaches1_px_loc(:,:,1:length(r1_C_idx),s,m) = bhv.reaches_inVec_px(:,:,r1_C_idx);
            reaches1_px_loc(:,:,1:length(r1_nD_idx),s,m) = bhv.reaches_inVec_px(:,:,r1_nD_idx);

            n_reaches_loc(1,s,m) = length(find(bhv.DCnD_inVec_idx==1 ...
                & bhv.cat_reach_inVec==1 | bhv.cat_reach_inVec==2));
            n_reaches_loc(2,s,m) = length(find(bhv.DCnD_inVec_idx==2 ...
                & bhv.cat_reach_inVec==1 | bhv.cat_reach_inVec==2));
            n_reaches_loc(3,s,m) = length(find(bhv.DCnD_inVec_idx==3 ...
                & bhv.cat_reach_inVec==1 | bhv.cat_reach_inVec==2));

            n_reaches_loc_LCR(1,s,m) = length(find(bhv.LCR_inVec_idx==1));
            n_reaches_loc_LCR(2,s,m) = length(find(bhv.LCR_inVec_idx==2));
            n_reaches_loc_LCR(3,s,m) = length(find(bhv.LCR_inVec_idx==3));

            n_reaches_cat12_loc(1,1,s,m) = length(r1_D_idx);
            n_reaches_cat12_loc(2,1,s,m) = length(r1_C_idx);
            n_reaches_cat12_loc(3,1,s,m) = length(r1_nD_idx);

            n_reaches_cat12_loc(1,2,s,m) = length(find(bhv.DCnD_inVec_idx==1 & bhv.cat_reach_inVec==2));
            n_reaches_cat12_loc(2,2,s,m) = length(find(bhv.DCnD_inVec_idx==2 & bhv.cat_reach_inVec==2));
            n_reaches_cat12_loc(3,2,s,m) = length(find(bhv.DCnD_inVec_idx==3 & bhv.cat_reach_inVec==2));

            n_suc_loc(1,s,m) = length(find(bhv.DCnD_inVec_idx==1 & bhv.suc_inVec==1));
            n_suc_loc(2,s,m) = length(find(bhv.DCnD_inVec_idx==2 & bhv.suc_inVec==1));
            n_suc_loc(3,s,m) = length(find(bhv.DCnD_inVec_idx==3 & bhv.suc_inVec==1));

            n_hit_loc(1,s,m) = length(find(bhv.DCnD_inVec_idx==1 & bhv.hit_inVec==1));
            n_hit_loc(2,s,m) = length(find(bhv.DCnD_inVec_idx==2 & bhv.hit_inVec==1));
            n_hit_loc(3,s,m) = length(find(bhv.DCnD_inVec_idx==3 & bhv.hit_inVec==1));

            duration_f_tmp_inVec = reaches.start_forw.dur_forw_mat(bhv.reach_inVec_idx);
            duration_f_loc(1:length(r1_D_idx),1,s,m) = duration_f_tmp_inVec(r1_D_idx);
            duration_f_loc(1:length(r1_C_idx),2,s,m) = duration_f_tmp_inVec(r1_C_idx);
            duration_f_loc(1:length(r1_nD_idx),3,s,m) = duration_f_tmp_inVec(r1_nD_idx);

    end
end

%% Figure
% 
sz=40;
load(fullfile(mat_dir_r,mouse,sess,"behavior_session.mat"),'behavior');
clr_dom = behavior.colors.clr_dom;
clr_center = behavior.colors.center_color;
clr_nonDom = behavior.colors.clr_nondom;
clr_left = behavior.colors.left_color;
clr_right = behavior.colors.right_color;

idx_R = [1,4,5];
idx_L = [2,3];
   
x_jitterR = [...
    1 + 0.01 * (rand(length(idx_R), 1) - 0.5);  
    2 + 0.01 * (rand(length(idx_R), 1) - 0.5);  
    3 + 0.01 * (rand(length(idx_R), 1) - 0.5)]; 
x_jitterL = [...
    1 + 0.01 * (rand(length(idx_L), 1) - 0.5);  
    2 + 0.01 * (rand(length(idx_L), 1) - 0.5);  
    3 + 0.01 * (rand(length(idx_L), 1) - 0.5)]; 


clrsR=[repmat(clr_dom,length(idx_R),1);...
    repmat(clr_center,length(idx_R),1);...
    repmat(clr_nonDom,length(idx_R),1)];

clrsL=[repmat(clr_dom,length(idx_L),1);...
    repmat(clr_center,length(idx_L),1);...
    repmat(clr_nonDom,length(idx_L),1)];

clrs=[repmat(clr_dom,n_mice,1);...
    repmat(clr_center,n_mice,1);...
    repmat(clr_nonDom,n_mice,1)];

 axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
                0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};
%%

n_mean_sess = squeeze(round(mean(n_reaches_loc(:,:,:),2,'omitnan')))';
n_mean_sessR = squeeze(round(mean(n_reaches_loc(:,:,idx_R),2,'omitnan')))';
n_mean_sessL = squeeze(round(mean(n_reaches_loc(:,:,idx_L),2,'omitnan')))';

n_mean_hit = squeeze(round(mean(n_hit_loc,2,'omitnan')))';
n_mean_hitR = squeeze(round(mean(n_hit_loc(:,:,idx_R),2,'omitnan')))';
n_mean_hitL = squeeze(round(mean(n_hit_loc(:,:,idx_L),2,'omitnan')))';

n_mean_suc = squeeze(round(mean(n_suc_loc,2,'omitnan')))';
n_mean_sucR = squeeze(round(mean(n_suc_loc(:,:,idx_R),2,'omitnan')))';
n_mean_sucL = squeeze(round(mean(n_suc_loc(:,:,idx_L),2,'omitnan')))';

mean_dur = squeeze(mean(mean(duration_f_loc(:,:,:,:),1,'omitnan'),3,'omitnan'))';
mean_durR = squeeze(mean(mean(duration_f_loc(:,:,:,idx_R),1,'omitnan'),3,'omitnan'))';
mean_durL = squeeze(mean(mean(duration_f_loc(:,:,:,idx_L),1,'omitnan'),3,'omitnan'))';

x_jitter = [...
    1 + 0.2 * (rand(n_mice, 1) - 0.5);  
    2 + 0.2 * (rand(n_mice, 1) - 0.5);  
    3 + 0.2 * (rand(n_mice, 1) - 0.5)]; 

figure
ff=tiledlayout(1,4);
nexttile
boxplot(n_mean_sess,'Colors',[.9 .9 .9],'Widths',.4)
hold on
plot(n_mean_sess','color',[.5 .5 .5 .2], 'LineWidth',1.5);
scatter(x_jitterR,n_mean_sessR(:),sz,clrsR,'filled','MarkerFaceAlpha',.9); hold on
scatter(x_jitterL,n_mean_sessL(:),sz,clrsL,'filled','MarkerFaceAlpha',.9,'MarkerEdgeColor','k'); hold on
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['dominant paw side    ';'center               ';'non-dominant paw side'])
ylabel('number of reaches')


nexttile
boxplot(n_mean_hit,'Colors',[.9 .9 .9],'Widths',.4)
hold on
plot(n_mean_hit','color',[.5 .5 .5 .2], 'LineWidth',1.5);
scatter(x_jitterR,n_mean_hitR(:),sz,clrsR,'filled','MarkerFaceAlpha',.9); hold on
scatter(x_jitterL,n_mean_hitL(:),sz,clrsL,'filled','MarkerFaceAlpha',.9,'MarkerEdgeColor','k'); hold on
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['dominant paw side    ';'center               ';'non-dominant paw side'])
ylabel('number of hit reaches')


nexttile
boxplot(n_mean_suc,'Colors',[.9 .9 .9],'Widths',.4)
hold on
plot(n_mean_suc','color',[.5 .5 .5 .2], 'LineWidth',1.5);
scatter(x_jitterR,n_mean_sucR(:),sz,clrsR,'filled','MarkerFaceAlpha',.9); hold on
scatter(x_jitterL,n_mean_sucL(:),sz,clrsL,'filled','MarkerFaceAlpha',.9,'MarkerEdgeColor','k'); hold on
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['dominant paw side    ';'center               ';'non-dominant paw side'])
ylabel('number of successful reaches')

nexttile
boxplot(mean_dur,'Colors',[.9 .9 .9],'Widths',.4)
hold on
plot(mean_dur','color',[.5 .5 .5 .2], 'LineWidth',1.5);
scatter(x_jitterR,mean_durR(:),sz,clrsR,'filled','MarkerFaceAlpha',.9); hold on
scatter(x_jitterL,mean_durL(:),sz,clrsL,'filled','MarkerFaceAlpha',.9,'MarkerEdgeColor','k'); hold on
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['dominant paw side    ';'center               ';'non-dominant paw side'])
ylabel('number of successful reaches')



set(gcf,'Position',[2235         319        1075         474],'color','w');
saveas(gcf,strcat(save_out,filesep,'n_reaches_hit_suc_dur.png'),'png')
print(gcf,strcat(save_out,filesep,'n_reaches_hit_suc_dur.pdf'), '-dpdf', '-painters');


%%
tm_reach = reaches.tm_w;
reaches1_x = squeeze(reaches1_px_loc(:,1,:,:,1,:));
reaches1_y = squeeze(reaches1_px_loc(:,2,:,:,1,:));
reaches1_z = squeeze(reaches1_px_loc(:,3,:,:,1,:));
szz=40;

figure()
ff = tiledlayout(2,n_mice);

%m=1;
for m=1:n_mice

r1_D = reaches1_x(:,1:n_reaches_cat12_loc(1,1,1,m),1,m);
r1_C = reaches1_x(:,1:n_reaches_cat12_loc(2,1,1,m),2,m);
r1_nD = reaches1_x(:,1:n_reaches_cat12_loc(3,1,1,m),3,m);
r1_x=[r1_D,r1_C,r1_nD];
r1_D = reaches1_y(:,1:n_reaches_cat12_loc(1,1,1,m),1,m);
r1_C = reaches1_y(:,1:n_reaches_cat12_loc(2,1,1,m),2,m);
r1_nD = reaches1_y(:,1:n_reaches_cat12_loc(3,1,1,m),3,m);
r1_y=[r1_D,r1_C,r1_nD];
r1_D = reaches1_z(:,1:n_reaches_cat12_loc(1,1,1,m),1,m);
r1_C = reaches1_z(:,1:n_reaches_cat12_loc(2,1,1,m),2,m);
r1_nD = reaches1_z(:,1:n_reaches_cat12_loc(3,1,1,m),3,m);
r1_z=[r1_D,r1_C,r1_nD];
n_r = size(r1_x,2);

clrs_DCnD = [repmat(clr_dom,size(r1_D,2),1);...
    repmat(clr_center,size(r1_C,2),1);...
    repmat(clr_nonDom,size(r1_nD,2),1)];
if n_reaches_loc(1,1,m)==n_reaches_loc_LCR(3,1,m)
    clrs_LCR = [repmat(clr_right,size(r1_D,2),1);...
    repmat(clr_center,size(r1_C,2),1);...
    repmat(clr_left,size(r1_nD,2),1)];
else
    clrs_LCR = [repmat(clr_left,size(r1_D,2),1);...
    repmat(clr_center,size(r1_C,2),1);...
    repmat(clr_right,size(r1_nD,2),1)];
end


nexttile(m)
imagesc(tm_reach,1:n_r, r1_y'); hold on
xline(0,'--','Color',[1 1 1 .5],'LineWidth',2);
%xlabel('time (s)');
ylabel('reach idx');
set(gca,axeOpt{:});
c1=colorbar; ylabel(c1,'y (px) - ML'); hold on
scatter(ones(size(1:n_r)).*-.374,1:n_r, szz ,clrs_LCR, 'filled','Marker','square');  
scatter(ones(size(1:n_r)).*-.36,1:n_r, szz ,clrs_DCnD, 'filled','Marker','square');  

nexttile(m+n_mice)
imagesc(tm_reach,1:n_r, r1_x'); hold on
xline(0,'--','Color',[1 1 1 .5],'LineWidth',2);
xlabel('time (s)'); ylabel('reach idx');
set(gca,axeOpt{:});
c1=colorbar; ylabel(c1,'x (px) - AP'); hold on
scatter(ones(size(1:n_r)).*-.374,1:n_r, szz ,clrs_LCR, 'filled','Marker','square');  
scatter(ones(size(1:n_r)).*-.36,1:n_r, szz ,clrs_DCnD, 'filled','Marker','square');  

% nexttile(m+n_mice*2)
% imagesc(tm_reach,1:n_r, r1_z'); hold on
% xline(0,'--','Color',[1 1 1 .5],'LineWidth',2);
% xlabel('time (s)'); ylabel('reach idx');
% set(gca,axeOpt{:});
% c1=colorbar; ylabel(c1,'x (px) - AP'); hold on
% scatter(ones(size(1:n_r)).*-.374,1:n_r, szz ,clrs_LCR, 'filled','Marker','square');  
% scatter(ones(size(1:n_r)).*-.36,1:n_r, szz ,clrs_DCnD, 'filled','Marker','square');  
colormap(bone)
end
set(gcf,'Position',[1921          46        1920         963],'color','w');

saveas(gcf,strcat(save_out,filesep,'reaches_y_y.png'),'png')


%%
figure
sess1_D_dur=squeeze(duration_f_loc(:,1,1,:));
sess1_C_dur=squeeze(duration_f_loc(:,2,1,:));
sess1_nD_dur=squeeze(duration_f_loc(:,3,1,:));

n_bins = 10;
figure

ff = tiledlayout(1,n_mice);
for m=1:n_mice

    durD = sess1_D_dur(:,m);
    durC = sess1_C_dur(:,m);
    durnD = sess1_nD_dur(:,m);

    nexttile
    histogram(durD(~isnan(durD)),n_bins,'normalization','pdf','FaceColor',clr_dom,'EdgeColor','none'); hold on
    histogram(durC(~isnan(durC)),n_bins,'normalization','pdf','FaceColor',clr_center,'EdgeColor','none');
    histogram(durnD(~isnan(durnD)),n_bins,'normalization','pdf','FaceColor',clr_nonDom,'EdgeColor','none'); hold off

end
shg

%%
n_bins = 45;
dur_all_D = duration_f_loc(:,1,:,:);
dur_all_C = duration_f_loc(:,2,:,:);
dur_all_nD = duration_f_loc(:,3,:,:);


durD_flat = dur_all_D(:);
durC_flat = dur_all_C(:);
durnD_flat = dur_all_nD(:);

figure
    histogram(durD_flat(~isnan(durD_flat)),n_bins,'normalization','pdf','FaceColor',clr_dom,'EdgeColor','none'); hold on
    histogram(durC_flat(~isnan(durC_flat)),n_bins,'normalization','pdf','FaceColor',clr_center,'EdgeColor','none');
    histogram(durnD_flat(~isnan(durnD_flat)),n_bins,'normalization','pdf','FaceColor',clr_nonDom,'EdgeColor','none'); hold off


    %%
idx_R = [1,4,5];
idx_L = [2,3];
    mean_dur = squeeze(mean(mean(duration_f_loc(:,:,:,:),1,'omitnan'),3,'omitnan'))';
    mean_durR = squeeze(mean(mean(duration_f_loc(:,:,:,idx_R),1,'omitnan'),3,'omitnan'))';
    mean_durL = squeeze(mean(mean(duration_f_loc(:,:,:,idx_L),1,'omitnan'),3,'omitnan'))';
x_jitterR = [...
    1 + 0.01 * (rand(length(idx_R), 1) - 0.5);  
    2 + 0.01 * (rand(length(idx_R), 1) - 0.5);  
    3 + 0.01 * (rand(length(idx_R), 1) - 0.5)]; 
x_jitterL = [...
    1 + 0.01 * (rand(length(idx_L), 1) - 0.5);  
    2 + 0.01 * (rand(length(idx_L), 1) - 0.5);  
    3 + 0.01 * (rand(length(idx_L), 1) - 0.5)]; 

figure

clrsR=[repmat(clr_dom,length(idx_R),1);...
    repmat(clr_center,length(idx_R),1);...
    repmat(clr_nonDom,length(idx_R),1)];

clrsL=[repmat(clr_dom,length(idx_L),1);...
    repmat(clr_center,length(idx_L),1);...
    repmat(clr_nonDom,length(idx_L),1)];

boxplot(mean_dur,'Colors',[.9 .9 .9],'Widths',.4)
hold on
plot(mean_dur','color',[.5 .5 .5 .2], 'LineWidth',1.5);
scatter(x_jitterR,mean_durR(:),sz,clrsR,'filled','MarkerFaceAlpha',.9); hold on
scatter(x_jitterL,mean_durL(:),sz,clrsL,'filled','MarkerFaceAlpha',.9,'MarkerEdgeColor','k'); hold on
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['dominant paw side    ';'center               ';'non-dominant paw side'])
ylabel('number of successful reaches')




