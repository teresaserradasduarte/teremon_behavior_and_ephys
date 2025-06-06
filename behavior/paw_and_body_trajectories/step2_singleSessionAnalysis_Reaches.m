%% Session Analysis STEP2
% After extracting an categorizing reaches, display the data
% v3 -> v4 : add/replace origin from water to centre of resting bar
clear; close all; clc

%% Load data
raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\raw_data';
mat_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\mat_files';
out_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\output_files';

% group = '20201006_AtaxicMice_G2';
% setup = 'headfixed_setup1';
% mouse = 'B1_MRC-06855';

group = '20230511_ChocolateGroup';
setup = 'headfixed_dynamicTarget';
animals = {...
    'CoteDor',...
    'Lindt',...
    'Toblerone',...
    'Milka',...
    'FerreroRocher'};
animal_idx = 3;
mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});
sess = 'R3';
save_flag = true;

% Normalize reaches position
% 1 - origin is position of water
% 2 - origin is center of resting bar
% 3 - origin is video 0,0,0 in mm
% 4 - origin is video 0,0,0 in px
norm_type = 3;
if norm_type==1
    save_name='origin_water';
elseif norm_type==2
    save_name='origin_restingBar';
elseif norm_type==3
    save_name='origin_zero_mm';
else
    save_name='origin_zero_px';
end
save_mat_name = strcat('session_analysis_',save_name,'.mat');

%% Load data of single session
load(strcat(mat_folder,filesep,group,filesep,setup,filesep,mouse,filesep,sess,filesep,'session_reaching_data_paw.mat')) 
sess_name = convertCharsToStrings(sess);

% Plots colrs
if strcmpi(mouse_info.phenotype,'ctr')
    suc_clr = [74,98,116]./256;
elseif strcmpi(mouse_info.phenotype,'DTX')
    suc_clr = [234,162,33]./256;
elseif strcmpi(mouse_info.phenotype,'a2a')
    suc_clr = [216,27,96]./256;
else
    disp('phenotype was not recognized')
    suc_clr = [.5,.5,.5];
end
miss_clr = [140,86,86]./256;

% Single Session
close all
disp(strcat('Running session',{' '},sess_name,'...'))
new_folder_mat = strcat(mat_folder,filesep,group,filesep,setup,filesep,mouse,filesep,sess_name);
new_folder_out = strcat(out_folder,filesep,group,filesep,setup,filesep,mouse,filesep,char(sess_name),filesep,'reaches_singleLabel_analysis',save_name);
mkdir(new_folder_out)


%% Display the types of categories

perc1 = numel(find(reaches.cat_reach==1));
perc2 = numel(find(reaches.cat_reach==2));
perc3 = numel(find(reaches.cat_reach==3));
perc4 = numel(find(reaches.cat_reach==4));

figure(1)
perc_cat = [perc1 perc2 perc3 perc4];
explode=[0 1 1 1];
labels={'paw resting: ';'paw lifted: ';'grooming: ';'drink/reach: '};
pie_cat = pie(perc_cat,explode);
hPatch = findobj(pie_cat, 'Type', 'Patch');
pText = findobj(pie_cat,'Type','text');
percentValues = get(pText,'String');
combinedtxt = strcat(labels,percentValues);
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);

% Create legend
legend(labels,'Location','southeastoutside','Box','off');
set(gcf, 'Position', [2341 357 755 567])
if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'cat_pie_chart.png'),'png'); end
% Save mat
save(strcat(new_folder_mat,filesep,save_mat_name),'perc1','perc2','perc3','perc4');


%%  Divede purpose, hit and success dived by cat
%pp_no = numel(find(reaches.purpose_reach==0));
nopp_cat1 = numel(intersect(find(reaches.purpose_reach==0),find(reaches.cat_reach==0)));
nopp_cat2 = numel(intersect(find(reaches.purpose_reach==0),find(reaches.cat_reach==2)));
nopp_cat3 = numel(intersect(find(reaches.purpose_reach==0),find(reaches.cat_reach==3)));
nopp_cat4 = numel(intersect(find(reaches.purpose_reach==0),find(reaches.cat_reach==4)));
pp_cat1 = numel(intersect(find(reaches.purpose_reach==1),find(reaches.cat_reach==1)));
pp_cat2 = numel(intersect(find(reaches.purpose_reach==1),find(reaches.cat_reach==2)));
pp_cat3 = numel(intersect(find(reaches.purpose_reach==1),find(reaches.cat_reach==3)));
pp_cat4 = numel(intersect(find(reaches.purpose_reach==1),find(reaches.cat_reach==4)));

%hit_no = numel(find(reaches.hit_reach==0));
nohit_cat1 = numel(intersect(find(reaches.hit_reach==0),find(reaches.cat_reach==1)));
nohit_cat2 = numel(intersect(find(reaches.hit_reach==0),find(reaches.cat_reach==2)));
nohit_cat3 = numel(intersect(find(reaches.hit_reach==0),find(reaches.cat_reach==3)));
nohit_cat4 = numel(intersect(find(reaches.hit_reach==0),find(reaches.cat_reach==4)));
hit_cat1 = numel(intersect(find(reaches.hit_reach==1),find(reaches.cat_reach==1)));
hit_cat2 = numel(intersect(find(reaches.hit_reach==1),find(reaches.cat_reach==2)));
hit_cat3 = numel(intersect(find(reaches.hit_reach==1),find(reaches.cat_reach==3)));
hit_cat4 = numel(intersect(find(reaches.hit_reach==1),find(reaches.cat_reach==4)));

%suc_no = numel(find(reaches.success_reach==0));
nosuc_cat1 = numel(intersect(find(reaches.success_reach==0),find(reaches.cat_reach==1)));
nosuc_cat2 = numel(intersect(find(reaches.success_reach==0),find(reaches.cat_reach==2)));
nosuc_cat3 = numel(intersect(find(reaches.success_reach==10),find(reaches.cat_reach==3)));
nosuc_cat4 = numel(intersect(find(reaches.success_reach==10),find(reaches.cat_reach==4)));
suc_cat1 = numel(intersect(find(reaches.success_reach==1),find(reaches.cat_reach==1)));
suc_cat2 = numel(intersect(find(reaches.success_reach==1),find(reaches.cat_reach==2)));
suc_cat3 = numel(intersect(find(reaches.success_reach==1),find(reaches.cat_reach==3)));
suc_cat4 = numel(intersect(find(reaches.success_reach==1),find(reaches.cat_reach==4)));

fHand = figure(2);
aHand = axes('parent', fHand);
bar_class_cat = [...
    %    pp_no 0 0; ...
    nopp_cat1 nopp_cat2 nopp_cat3 nopp_cat4; ...
    pp_cat1 pp_cat2 pp_cat3 pp_cat4; 0 0 0 0; ...
    %    hit_no 0 0; ...
    nohit_cat1 nohit_cat2 nohit_cat3 nohit_cat4; ...
    hit_cat1 hit_cat2 hit_cat3 hit_cat4; 0 0 0 0; ...
    %    suc_no 0 0; ...
    nosuc_cat1 nosuc_cat2 nosuc_cat3 nosuc_cat4; ...
    suc_cat1 suc_cat2 suc_cat3 suc_cat4];
b1 = bar(bar_class_cat, 'stacked','FaceColor','flat');
cats_clr = zeros(4,3);
for k = 1:size(bar_class_cat,2)
    b1(k).CData = k;
    %cats_clr(k, :) = b1(k).FaceColor;
end
hold(aHand, 'on')
legend(labels,'Location','northeast');
set(gca, 'XTick', [1,2,4,5,7,8], 'XTickLabel', {'(no)','purpose','(no)','hit','(no)','success'})
shg
hold(aHand, 'off')
if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'purpose_hit_sucess.png'),'png'); end
% Save mat
cats_clr = parula(4);
cat1_clr = cats_clr(1,:);
cat2_clr = cats_clr(2,:);
cat3_clr = cats_clr(3,:);
cat4_clr = cats_clr(4,:);
save(strcat(new_folder_mat,filesep,save_mat_name),'bar_class_cat',...
    'cat1_clr','cat2_clr','cat3_clr','-append');

%% Show reaches organized by category - xyz over time
p2m = [session.video.calib.px2mm.px2mm_xz, session.video.calib.px2mm.px2mm_y, ...
    session.video.calib.px2mm.px2mm_xz];
width = session.video.width;
height = session.video.height;
tm =  reaches.tm_w;
tm_vel = tm(1:end-1);
% parts of reach
start=80;
start_mod2=136;
stop=180;

reach_choice = reaches.reach_mat;
reach_choice_otherPaw = reaches.reach_otherPaw_mat;
%reach_choice = reaches.s_reach_mat;

if strcmp(mouse_info.paw_pref,'L')
    reaches_y_flipped = abs(reach_choice(:,2,:)-width-1);
    reach_vec = reach_choice;
    reach_vec(:,2,:) = reaches_y_flipped; % dominant paw
    reaches_nonDom_y_flip = abs(reach_choice_otherPaw(:,2,:)-width-1);
    reach_nonDom_vec = reach_choice_otherPaw;
    reach_nonDom_vec(:,2,:) = reaches_nonDom_y_flip; % non-dominant paw
else
    reach_vec = reach_choice;
    reach_nonDom_vec = reach_choice_otherPaw;
end

% Normalize reaches position
% 1 - origin is position of water
% 2 - origin is center of resting bar
% 3 - origin is video 0,0,0 in mm
% 4 - origin is video 0,0,0 in px
%norm_type = 2;

if norm_type ==1  % Origin is water position
    % Convert to mm with (0,0,0)=water
    median_water_vec=nanmedian(reaches.water_reach_mat,1);
    % Replace nans by medians
    medofmed_water_vec = nanmedian(median_water_vec,3);
    median_water_vec(1,1,isnan(median_water_vec(1,1,:))) = repmat(medofmed_water_vec(1),...
        [1 1 numel(median_water_vec(1,1,isnan(median_water_vec(1,1,:))))]);
    median_water_vec(1,2,isnan(median_water_vec(1,2,:))) = repmat(medofmed_water_vec(2),...
        [1 1 numel(median_water_vec(1,2,isnan(median_water_vec(1,2,:))))]);
    median_water_vec(1,3,isnan(median_water_vec(1,3,:))) = repmat(medofmed_water_vec(3),...
        [1 1 numel(median_water_vec(1,3,isnan(median_water_vec(1,3,:))))]);

    reaches_norm=(reach_vec-repmat(median_water_vec,[size(reach_vec,1) 1 1])).*p2m;
    reaches_nonDom_norm=(reach_nonDom_vec-repmat(median_water_vec,[size(reach_vec,1) 1 1])).*p2m;
    origin_pos = (nanmean(nanmedian(reaches.water_reach_mat,1),3)).*p2m;
    origin=origin_pos-origin_pos;

elseif norm_type ==2 % Origin is center of resting bar
    % Origin
    camA_origin = [81,312];
    camB_origin = [246,381];
    camC_origin = [458,384];
    if strcmp(mouse_info.paw_pref,'L')
        x_origin = abs(camC_origin(1)-width-1);
    else
        x_origin = camA_origin(1);
    end
    y_origin = camB_origin(1);
    z_origin = abs(camB_origin(2)-height-1);
    origin_px=[x_origin,y_origin,z_origin];
    origin_rep=repmat(origin_px,[size(reach_vec,1) 1 size(reach_vec,3)]);
    origin = origin_px-origin_px;

    reaches_norm = cat(2,(reach_vec(:,1,:)-origin_rep(:,1,:)).*p2m.px2mm_x,...
        (reach_vec(:,2:3,:)-origin_rep(:,2:3,:)).*p2m.px2mm_yz);
    reaches_nonDom_norm = cat(2,(reach_nonDom_vec(:,1,:)-origin_rep(:,1,:)).*p2m.px2mm_x,...
        (reach_nonDom_vec(:,2:3,:)-origin_rep(:,2:3,:)).*p2m.px2mm_yz);

elseif norm_type ==3 % Origin is zero of video, in mm
    reaches_norm = cat(2,reach_vec(:,1,:).*p2m(1),reach_vec(:,2,:).*p2m(2),reach_vec(:,3,:).*p2m(3));
    reaches_nonDom_norm = cat(2,reach_nonDom_vec(:,1,:).*p2m(1),reach_nonDom_vec(:,2,:).*p2m(2),reach_nonDom_vec(:,3,:).*p2m(3));
    origin = [0,0,0];
else                % Origin is zero of video, in pixels
    reaches_norm = reach_vec;
    reaches_nonDom_norm = reach_nonDom_vec;
    origin = [0,0,0];
end

% Category 1
cat1_all = reaches_norm(:,:,reaches.cat_reach==1);
cat1_all_mean = nanmean(cat1_all,3);
% Category 2
cat2_all = reaches_norm(:,:,reaches.cat_reach==2);
cat2_all_mean = nanmean(cat2_all,3);
% Category 3
cat3_all = reaches_norm(:,:,reaches.cat_reach==3);
cat3_all_mean = nanmean(cat3_all,3);
% Category 3
cat4_all = reaches_norm(:,:,reaches.cat_reach==4);
cat4_all_mean = nanmean(cat4_all,3);

% -----------------
% Distance and velocity
distance_mm = sqrt(squeeze(reaches_norm(:,1,:)).^2 + ...
    squeeze(reaches_norm(:,2,:)).^2 + ...
    squeeze(reaches_norm(:,3,:)).^2);
velocity_mm = diff(distance_mm);

% Category 1
cat1_all_dist = distance_mm(:,reaches.cat_reach==1);
cat1_all_dist_mean = nanmean(cat1_all_dist,2);
cat1_all_vel = velocity_mm(:,reaches.cat_reach==1);
cat1_all_vel_mean = nanmean(cat1_all_vel,2);

% Category 2
cat2_all_dist = distance_mm(:,reaches.cat_reach==2);
cat2_all_dist_mean = nanmean(cat2_all_dist,2);
cat2_all_vel = velocity_mm(:,reaches.cat_reach==2);
cat2_all_vel_mean = nanmean(cat2_all_vel,2);

% Category 3
cat3_all_dist = distance_mm(:,reaches.cat_reach==3);
cat3_all_dist_mean = nanmean(cat3_all_dist,2);
cat3_all_vel = velocity_mm(:,reaches.cat_reach==3);
cat3_all_vel_mean = nanmean(cat3_all_vel,2);

% Category 4
cat4_all_dist = distance_mm(:,reaches.cat_reach==4);
cat4_all_dist_mean = nanmean(cat4_all_dist,2);
cat4_all_vel = velocity_mm(:,reaches.cat_reach==4);
cat4_all_vel_mean = nanmean(cat4_all_vel,2);

%Save mat
save(strcat(new_folder_mat,filesep,save_mat_name),'cat1_all','cat1_all_mean',...
    'cat2_all','cat2_all_mean','cat3_all','cat3_all_mean',...
    'cat4_all','cat4_all_mean','tm','tm_vel',...
    'cat1_clr','cat2_clr','cat3_clr','cat4_clr','cat3_clr',...
    'cat1_all_dist','cat1_all_dist_mean','cat1_all_vel','cat1_all_vel_mean',...
    'cat2_all_dist','cat2_all_dist_mean','cat2_all_vel','cat2_all_vel_mean',...
    'cat3_all_dist','cat3_all_dist_mean','cat3_all_vel','cat3_all_vel_mean',...
    'cat4_all_dist','cat4_all_dist_mean','cat4_all_vel','cat4_all_vel_mean',...
    'reaches_norm','reaches_nonDom_norm','-append');


% Plotting
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0.05,'ticklength',[1,1]*.01};


transp = 0.015;
transp_mean_proj = 1;
mean_p_lw=1.5;
mean_r_lw=2.5;
trial_lw = 1;
if norm_type==1
    x_min=-18; x_max=4; y_min=-12; y_max=2; z_min=-10; z_max=8;
elseif norm_type==2
    x_min=-5; x_max=18; y_min=-20; y_max=12; z_min=-5; z_max=19;
elseif norm_type==3
    x_min = 0;
    x_max = 40;
    y_max = 20;
    y_min = 0;
    z_max = 30;
    z_min = 0;
end

oneMat = ones(size(cat1_all,1),1);
axis_flag = true;

figure(3)
subplot(2,2,[1,3])
plot3(cat4_all_mean(:,1),cat4_all_mean(:,2),cat4_all_mean(:,3),'LineWidth',mean_r_lw,'Color',cat4_clr); hold on
plot3(cat3_all_mean(:,1),cat3_all_mean(:,2),cat3_all_mean(:,3),'LineWidth',mean_r_lw,'Color',cat3_clr); hold on
plot3(cat2_all_mean(:,1),cat2_all_mean(:,2),cat2_all_mean(:,3),'LineWidth',mean_r_lw,'Color',cat2_clr);
plot3(cat1_all_mean(:,1),cat1_all_mean(:,2),cat1_all_mean(:,3),'LineWidth',mean_r_lw,'Color',cat1_clr);
plot3(origin(1),origin(2),origin(3),'o')

%grid on
view(-45,18)
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1.5,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top');
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
axis image;
if axis_flag, axis ([x_min x_max y_min y_max z_min z_max]); end
yL = get(gca,'YLim');
zL = get(gca,'ZLim');
plot3(squeeze(cat4_all(:,1,:)), oneMat .* yL(2), squeeze(cat4_all(:,3,:)),'-','Color',cat(2,cat4_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat4_all(:,1,:)), squeeze(cat4_all(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,cat4_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat3_all(:,1,:)), oneMat .* yL(2), squeeze(cat3_all(:,3,:)),'-','Color',cat(2,cat3_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat3_all(:,1,:)), squeeze(cat3_all(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,cat3_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat2_all(:,1,:)), oneMat .* yL(2), squeeze(cat2_all(:,3,:)),'-','Color',cat(2,cat2_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat2_all(:,1,:)), squeeze(cat2_all(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,cat2_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat1_all(:,1,:)), oneMat .* yL(2), squeeze(cat1_all(:,3,:)),'-','Color',cat(2,cat1_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat1_all(:,1,:)), squeeze(cat1_all(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,cat1_clr,transp),'LineWidth',trial_lw);

plot3(cat4_all_mean(:,1), oneMat .* yL(2), cat4_all_mean(:,3),'-','Color',cat(2,cat4_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat4_all_mean(:,1), cat4_all_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,cat4_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat3_all_mean(:,1), oneMat .* yL(2), cat3_all_mean(:,3),'-','Color',cat(2,cat3_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat3_all_mean(:,1), cat3_all_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,cat3_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat2_all_mean(:,1), oneMat .* yL(2), cat2_all_mean(:,3),'-','Color',cat(2,cat2_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat2_all_mean(:,1), cat2_all_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,cat2_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat1_all_mean(:,1), oneMat .* yL(2), cat1_all_mean(:,3),'-','Color',cat(2,cat1_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat1_all_mean(:,1), cat1_all_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,cat1_clr,transp_mean_proj),'LineWidth',mean_p_lw);
hold off

subplot(2,2,2)
plot(tm,cat4_all_dist,'Color',cat(2,cat4_clr,transp)); hold on
plot(tm,cat4_all_dist_mean,'Color',cat4_clr,'Linewidth',2);
plot(tm,cat3_all_dist,'Color',cat(2,cat3_clr,transp)); hold on
plot(tm,cat3_all_dist_mean,'Color',cat3_clr,'Linewidth',2);
plot(tm,cat2_all_dist,'Color',cat(2,cat2_clr,transp));
plot(tm,cat2_all_dist_mean,'Color',cat2_clr,'Linewidth',2);
plot(tm,cat1_all_dist,'Color',cat(2,cat1_clr,transp));
plot(tm,cat1_all_dist_mean,'Color',cat1_clr,'Linewidth',2);
xline(tm(start),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(stop),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(start_mod2),'--','linewidth',2,'Color',[0.9 0.9 0.9 0.1])
hold off
xlabel('time (sec)'); ylabel('distance (mm)')
title('distance to water');
xlim([tm(1),tm(end)]); ylim([10 45])
axis square
set(gca,axeOpt{:});

transp=0.01;
subplot(2,2,4)
plot(tm_vel,cat4_all_vel,'Color',cat(2,cat4_clr,transp)); hold on
plot(tm_vel,cat3_all_vel,'Color',cat(2,cat3_clr,transp));
plot(tm_vel,cat2_all_vel,'Color',cat(2,cat2_clr,transp));
plot(tm_vel,cat1_all_vel,'Color',cat(2,cat1_clr,transp)); 
plot(tm_vel,cat4_all_vel_mean,'Color',cat4_clr,'Linewidth',2);
plot(tm_vel,cat3_all_vel_mean,'Color',cat3_clr,'Linewidth',2);
plot(tm_vel,cat2_all_vel_mean,'Color',cat2_clr,'Linewidth',2);
plot(tm_vel,cat1_all_vel_mean,'Color',cat1_clr,'Linewidth',2);
xline(tm(start),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(stop),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(start_mod2),'--','linewidth',2,'Color',[0.9 0.9 0.9 0.1])
hold off
hold off
xlabel('time (sec)'); ylabel('velocuty (mm/s)')
title('velocity');
axis square;
xlim([tm(1),tm(end)]); ylim([-.8 .8])
set(gca,axeOpt{:});

set(gcf,'color','w','Position',[2129         138        1067         820])
if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'3d_cat1234_dist_vel.png'),'png'); end
% saveas(gcf,strcat(new_folder_out,filesep,'3d_cat123'),'epsc');

%%
%transp = 0;
% 3D PLOT, CATEGORIES SEPERATED
figure(4)

subplot(221)
plot3(cat1_all_mean(:,1),cat1_all_mean(:,2),cat1_all_mean(:,3),'LineWidth',mean_r_lw,'Color',suc_clr); hold on
plot3(origin(1),origin(2),origin(3),'o')
view(-45,18)
grid on
title('SP: resting paw')
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top','Gridalpha',0.1);
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
axis ([x_min x_max y_min y_max z_min z_max])
yL = get(gca,'YLim');
zL = get(gca,'ZLim');
plot3(squeeze(cat1_all(:,1,:)), oneMat .* yL(2), squeeze(cat1_all(:,3,:)),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat1_all(:,1,:)), squeeze(cat1_all(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
plot3(cat1_all_mean(:,1), oneMat .* yL(2), cat1_all_mean(:,3),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat1_all_mean(:,1), cat1_all_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
hold off

subplot(222)
plot3(cat2_all_mean(:,1),cat2_all_mean(:,2),cat2_all_mean(:,3),'LineWidth',mean_r_lw,'Color',suc_clr); hold on
plot3(origin(1),origin(2),origin(3),'o')

view(-45,18)
title('SP: lifted paw')
grid on
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top','Gridalpha',0.1);
axis ([x_min x_max y_min y_max z_min z_max])
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
yL = get(gca,'YLim');
zL = get(gca,'ZLim');
plot3(squeeze(cat2_all(:,1,:)), oneMat .* yL(2), squeeze(cat2_all(:,3,:)),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat2_all(:,1,:)), squeeze(cat2_all(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
plot3(cat2_all_mean(:,1), oneMat .* yL(2), cat2_all_mean(:,3),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat2_all_mean(:,1), cat2_all_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(origin(1),origin(2),origin(3),'o')

hold off

subplot(223)
plot3(cat3_all_mean(:,1),cat3_all_mean(:,2),cat3_all_mean(:,3),'LineWidth',mean_r_lw,'Color',suc_clr); hold on
plot3(origin(1),origin(2),origin(3),'o')
view(-45,18)
title('grooming')
grid on
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top','Gridalpha',0.1);
axis ([x_min x_max y_min y_max z_min z_max])
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
yL = get(gca,'YLim');
zL = get(gca,'ZLim');
plot3(squeeze(cat3_all(:,1,:)), oneMat .* yL(2), squeeze(cat3_all(:,3,:)),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat3_all(:,1,:)), squeeze(cat3_all(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
plot3(cat3_all_mean(:,1), oneMat .* yL(2), cat3_all_mean(:,3),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat3_all_mean(:,1), cat3_all_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
hold off

subplot(224)
plot3(cat4_all_mean(:,1),cat4_all_mean(:,2),cat4_all_mean(:,3),'LineWidth',mean_r_lw,'Color',suc_clr); hold on
plot3(origin(1),origin(2),origin(3),'o')
view(-45,18)
title('drink / reach')
grid on
set(gca,...
    'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1,...
    'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top','Gridalpha',0.1);
axis ([x_min x_max y_min y_max z_min z_max])
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
yL = get(gca,'YLim');
zL = get(gca,'ZLim');
plot3(squeeze(cat4_all(:,1,:)), oneMat .* yL(2), squeeze(cat4_all(:,3,:)),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
plot3(squeeze(cat4_all(:,1,:)), squeeze(cat4_all(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
plot3(cat4_all_mean(:,1), oneMat .* yL(2), cat4_all_mean(:,3),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
plot3(cat4_all_mean(:,1), cat4_all_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
hold off


set(gcf,'Position',[2195         120         983         777])
if save_flag
    saveas(gcf,strcat(new_folder_out,filesep,'3d_cat1234.png'),'png');
    saveas(gcf,strcat(new_folder_out,filesep,'3d_cat1234.eps'),'epsc');
end

%% Non-dominapt paw
distance_mm_nonDom = sqrt(squeeze(reaches_nonDom_norm(:,1,:)).^2 + ...
    squeeze(reaches_nonDom_norm(:,2,:)).^2 + ...
    squeeze(reaches_nonDom_norm(:,3,:)).^2);

% Category 1 (success)
cat1_all_dist_nonDom = distance_mm_nonDom(:,reaches.cat_reach==1);
cat1_all_dist_nonDom_mean = nanmean(cat1_all_dist_nonDom,2);
% Category 2 (success)
cat2_all_dist_nonDom = distance_mm_nonDom(:,reaches.cat_reach==2);
cat2_all_dist_nonDom_mean = nanmean(cat2_all_dist_nonDom,2);
% Category 3 (all)
cat3_all_dist_nonDom = distance_mm_nonDom(:,reaches.cat_reach==3);
cat3_all_dist_nonDom_mean = nanmean(cat3_all_dist_nonDom,2);
% Category 4 (all)
cat4_all_dist_nonDom = distance_mm_nonDom(:,reaches.cat_reach==4);
cat4_all_dist_nonDom_mean = nanmean(cat4_all_dist_nonDom,2);

% Colors
dom_paw_clr = [0 0.4470 0.7410];
noDom_paw_clr = [0.8500 0.3250 0.0980];

% Save mat
save(strcat(new_folder_mat,filesep,save_mat_name),'cat1_all_dist_nonDom','cat1_all_dist_nonDom_mean',...
    'cat2_all_dist_nonDom','cat2_all_dist_nonDom_mean','cat3_all_dist_nonDom','cat3_all_dist_nonDom_mean',...
    'cat4_all_dist_nonDom','cat4_all_dist_nonDom_mean',...
    'dom_paw_clr','noDom_paw_clr','-append');

figure(5)
transp=0.1;
subplot(221)
plot(tm,cat1_all_dist_nonDom,'Color',cat(2,noDom_paw_clr,transp)); hold on
plot(tm,cat1_all_dist_nonDom_mean,'Color',noDom_paw_clr,'Linewidth',2);
plot(tm,cat1_all_dist,'Color',cat(2,dom_paw_clr,transp)); hold on
plot(tm,cat1_all_dist_mean,'Color',dom_paw_clr,'Linewidth',2);
hold off
xlabel('time (sec)'); ylabel('distance (mm)')
title('SP: paw resting'); 
xlim([tm(1),tm(end)]); 
axis square
set(gca,axeOpt{:});

subplot(222)
plot(tm,cat2_all_dist_nonDom,'Color',cat(2,noDom_paw_clr,transp)); hold on
plot(tm,cat2_all_dist_nonDom_mean,'Color',noDom_paw_clr,'Linewidth',2);
plot(tm,cat2_all_dist,'Color',cat(2,dom_paw_clr,transp)); hold on
plot(tm,cat2_all_dist_mean,'Color',dom_paw_clr,'Linewidth',2);
hold off
xlabel('time (sec)'); ylabel('distance (mm)')
title('SP: paw lifted'); 
xlim([tm(1),tm(end)]); 
axis square
set(gca,axeOpt{:});

subplot(223)
plot(tm,cat3_all_dist_nonDom,'Color',cat(2,noDom_paw_clr,transp)); hold on
plot(tm,cat3_all_dist_nonDom_mean,'Color',noDom_paw_clr,'Linewidth',2);
plot(tm,cat3_all_dist,'Color',cat(2,dom_paw_clr,transp));
plot(tm,cat3_all_dist_mean,'Color',dom_paw_clr,'Linewidth',2);
hold off
xlabel('time (sec)'); ylabel('distance (mm)')
title('grooming'); 
xlim([tm(1),tm(end)]); 
axis square
set(gca,axeOpt{:});

subplot(224)
plot(tm,cat4_all_dist_nonDom,'Color',cat(2,noDom_paw_clr,transp)); hold on
plot(tm,cat4_all_dist_nonDom_mean,'Color',noDom_paw_clr,'Linewidth',2);
plot(tm,cat4_all_dist,'Color',cat(2,dom_paw_clr,transp));
plot(tm,cat4_all_dist_mean,'Color',dom_paw_clr,'Linewidth',2);
hold off
xlabel('time (sec)'); ylabel('distance (mm)')
title('drink / reach'); 
xlim([tm(1),tm(end)]); 
axis square
set(gca,axeOpt{:});

set(gcf,'Position',[2335          90         883         771]);
if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'twoPaws_dist'),'png'); end
% saveas(gcf,strcat(new_folder_out,filesep,'twoPaws_dist'),'epsc');


%% Distance to mean
nbins = 40;

% SPACE-TIME -------------------------
% Category 1
cat1_all_distMean_spTm= sqrt((squeeze(cat1_all(:,1,:))-squeeze(cat1_all_mean(:,1))).^2 + ...
    (squeeze(cat1_all(:,2,:))-squeeze(cat1_all_mean(:,2))).^2 + ...
    (squeeze(cat1_all(:,3,:))-squeeze(cat1_all_mean(:,3))).^2);
cat1_all_distMean_spTm_int = sum(cat1_all_distMean_spTm(start:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat1_all_spTm, binsMean_cat1_all_spTm] = histcounts(cat1_all_distMean_spTm_int,nbins);
cdfMean_cat1_all_spTm = cumsum(countsMean_cat1_all_spTm); cdfMean_cat1_all_spTm(end+1)=cdfMean_cat1_all_spTm(end);
cdfMean_cat1_all_spTm_norm=cdfMean_cat1_all_spTm./cdfMean_cat1_all_spTm(end);
% grasp range
cat1_all_distMean_spTm_int_mod2 = sum(cat1_all_distMean_spTm(start_mod2:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat1_all_spTm_mod2, binsMean_cat1_all_spTm_mod2] = histcounts(cat1_all_distMean_spTm_int_mod2,nbins);
cdfMean_cat1_all_spTm_mod2 = cumsum(countsMean_cat1_all_spTm_mod2); cdfMean_cat1_all_spTm_mod2(end+1)=cdfMean_cat1_all_spTm_mod2(end);
cdfMean_cat1_all_spTm_norm_mod2=cdfMean_cat1_all_spTm_mod2./cdfMean_cat1_all_spTm_mod2(end);


% Category 2
cat2_all_distMean_spTm= sqrt((squeeze(cat2_all(:,1,:))-squeeze(cat2_all_mean(:,1))).^2 + ...
    (squeeze(cat2_all(:,2,:))-squeeze(cat2_all_mean(:,2))).^2 + ...
    (squeeze(cat2_all(:,3,:))-squeeze(cat2_all_mean(:,3))).^2);
cat2_all_distMean_spTm_int = sum(cat2_all_distMean_spTm(start:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat2_all_spTm, binsMean_cat2_all_spTm] = histcounts(cat2_all_distMean_spTm_int,nbins);
cdfMean_cat2_all_spTm = cumsum(countsMean_cat2_all_spTm); cdfMean_cat2_all_spTm(end+1)=cdfMean_cat2_all_spTm(end);
cdfMean_cat2_all_spTm_norm=cdfMean_cat2_all_spTm./cdfMean_cat2_all_spTm(end);
% grasp range
cat2_all_distMean_spTm_int_mod2 = sum(cat2_all_distMean_spTm(start_mod2:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat2_all_spTm_mod2, binsMean_cat2_all_spTm_mod2] = histcounts(cat2_all_distMean_spTm_int_mod2,nbins);
cdfMean_cat2_all_spTm_mod2 = cumsum(countsMean_cat2_all_spTm_mod2); cdfMean_cat2_all_spTm_mod2(end+1)=cdfMean_cat2_all_spTm_mod2(end);
cdfMean_cat2_all_spTm_norm_mod2=cdfMean_cat2_all_spTm_mod2./cdfMean_cat2_all_spTm_mod2(end);

% Category 3
cat3_all_distMean_spTm= sqrt((squeeze(cat3_all(:,1,:))-squeeze(cat3_all_mean(:,1))).^2 + ...
    (squeeze(cat3_all(:,2,:))-squeeze(cat3_all_mean(:,2))).^2 + ...
    (squeeze(cat3_all(:,3,:))-squeeze(cat3_all_mean(:,3))).^2);
cat3_all_distMean_spTm_int = sum(cat3_all_distMean_spTm(start:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat3_all_spTm, binsMean_cat3_all_spTm] = histcounts(cat3_all_distMean_spTm_int,nbins);
cdfMean_cat3_all_spTm = cumsum(countsMean_cat3_all_spTm); cdfMean_cat3_all_spTm(end+1)=cdfMean_cat3_all_spTm(end);
cdfMean_cat3_all_spTm_norm=cdfMean_cat3_all_spTm./cdfMean_cat3_all_spTm(end);
% grasp range
cat3_all_distMean_spTm_int_mod2 = sum(cat3_all_distMean_spTm(start_mod2:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat3_all_spTm_mod2, binsMean_cat3_all_spTm_mod2] = histcounts(cat3_all_distMean_spTm_int_mod2,nbins);
cdfMean_cat3_all_spTm_mod2 = cumsum(countsMean_cat3_all_spTm_mod2); cdfMean_cat3_all_spTm_mod2(end+1)=cdfMean_cat3_all_spTm_mod2(end);
cdfMean_cat3_all_spTm_norm_mod2=cdfMean_cat3_all_spTm_mod2./cdfMean_cat3_all_spTm_mod2(end);



% SPACE-ONLY -------------------------
n_timepoints = size(cat1_all,1);
% Category 1
n_trials_cat1_all = size(cat1_all,3);
cat1_all_distMean_sp = zeros(size(cat1_all_distMean_spTm));
for tt=1:n_trials_cat1_all
    for pt=1:n_timepoints
        % TO  THE MEAN
        tmp_dist_allMean=sqrt((cat1_all(pt,1,tt)-cat1_all_mean(:,1)).^2+...
            (cat1_all(pt,2,tt)-cat1_all_mean(:,2)).^2+...
            (cat1_all(pt,3,tt)-cat1_all_mean(:,3)).^2);
        [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
        cat1_all_distMean_sp(pt,tt)=min_dist_mean;
    end
end
cat1_all_distMean_sp_int = sum(cat1_all_distMean_sp(start:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat1_all_sp, binsMean_cat1_all_sp] = histcounts(cat1_all_distMean_sp_int,nbins);
cdfMean_cat1_all_sp = cumsum(countsMean_cat1_all_sp); cdfMean_cat1_all_sp(end+1)=cdfMean_cat1_all_sp(end);
cdfMean_cat1_all_sp_norm=cdfMean_cat1_all_sp./cdfMean_cat1_all_sp(end);
% range grasp
cat1_all_distMean_sp_int_mod2 = sum(cat1_all_distMean_sp(start_mod2:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat1_all_sp_mod2, binsMean_cat1_all_sp_mod2] = histcounts(cat1_all_distMean_sp_int_mod2,nbins);
cdfMean_cat1_all_sp_mod2 = cumsum(countsMean_cat1_all_sp_mod2); cdfMean_cat1_all_sp_mod2(end+1)=cdfMean_cat1_all_sp_mod2(end);
cdfMean_cat1_all_sp_norm_mod2=cdfMean_cat1_all_sp_mod2./cdfMean_cat1_all_sp_mod2(end);

% Category 2
n_trials_cat2_all = size(cat2_all,3);
cat2_all_distMean_sp = zeros(size(cat2_all_distMean_spTm));
for tt=1:n_trials_cat2_all
    for pt=1:n_timepoints
        % TO  THE MEAN
        tmp_dist_allMean=sqrt((cat2_all(pt,1,tt)-cat2_all_mean(:,1)).^2+...
            (cat2_all(pt,2,tt)-cat2_all_mean(:,2)).^2+...
            (cat2_all(pt,3,tt)-cat2_all_mean(:,3)).^2);
        [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
        cat2_all_distMean_sp(pt,tt)=min_dist_mean;
    end
end
cat2_all_distMean_sp_int = sum(cat2_all_distMean_sp(start:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat2_all_sp, binsMean_cat2_all_sp] = histcounts(cat2_all_distMean_sp_int,nbins);
cdfMean_cat2_all_sp = cumsum(countsMean_cat2_all_sp); cdfMean_cat2_all_sp(end+1)=cdfMean_cat2_all_sp(end);
cdfMean_cat2_all_sp_norm=cdfMean_cat2_all_sp./cdfMean_cat2_all_sp(end);
% range grasp
cat2_all_distMean_sp_int_mod2 = sum(cat2_all_distMean_sp(start_mod2:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat2_all_sp_mod2, binsMean_cat2_all_sp_mod2] = histcounts(cat2_all_distMean_sp_int_mod2,nbins);
cdfMean_cat2_all_sp_mod2 = cumsum(countsMean_cat2_all_sp_mod2); cdfMean_cat2_all_sp_mod2(end+1)=cdfMean_cat2_all_sp_mod2(end);
cdfMean_cat2_all_sp_norm_mod2 = cdfMean_cat2_all_sp_mod2./cdfMean_cat2_all_sp_mod2(end);

% Category 3
n_trials_cat3_all = size(cat3_all,3);
cat3_all_distMean_sp = zeros(size(cat3_all_distMean_spTm));
for tt=1:n_trials_cat3_all
    for pt=1:n_timepoints
        % TO  THE MEAN
        tmp_dist_allMean=sqrt((cat3_all(pt,1,tt)-cat3_all_mean(:,1)).^2+...
            (cat3_all(pt,2,tt)-cat3_all_mean(:,2)).^2+...
            (cat3_all(pt,3,tt)-cat3_all_mean(:,3)).^2);
        [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
        cat3_all_distMean_sp(pt,tt)=min_dist_mean;
    end
end
cat3_all_distMean_sp_int = sum(cat3_all_distMean_sp(start:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat3_all_sp, binsMean_cat3_all_sp] = histcounts(cat3_all_distMean_sp_int,nbins);
cdfMean_cat3_all_sp = cumsum(countsMean_cat3_all_sp); cdfMean_cat3_all_sp(end+1)=cdfMean_cat3_all_sp(end);
cdfMean_cat3_all_sp_norm=cdfMean_cat3_all_sp./cdfMean_cat3_all_sp(end);
% Range grasp
cat3_all_distMean_sp_int_mod2 = sum(cat3_all_distMean_sp(start_mod2:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat3_all_sp_mod2, binsMean_cat3_all_sp_mod2] = histcounts(cat3_all_distMean_sp_int_mod2,nbins);
cdfMean_cat3_all_sp_mod2 = cumsum(countsMean_cat3_all_sp_mod2); cdfMean_cat3_all_sp_mod2(end+1)=cdfMean_cat3_all_sp_mod2(end);
cdfMean_cat3_all_sp_norm_mod2=cdfMean_cat3_all_sp_mod2./cdfMean_cat3_all_sp_mod2(end);


% Save mat
save(strcat(new_folder_mat,filesep,save_mat_name),...
    'cat1_all_distMean_spTm','cat2_all_distMean_spTm','cat3_all_distMean_spTm',...
    'cat1_all_distMean_sp','cat2_all_distMean_sp','cat3_all_distMean_sp',...
    'cat1_all_distMean_spTm_int','cat2_all_distMean_spTm_int','cat3_all_distMean_spTm_int',...
    'cat1_all_distMean_sp_int','cat2_all_distMean_sp_int','cat3_all_distMean_sp_int',...
    'cat1_all_distMean_spTm_int_mod2','cat2_all_distMean_spTm_int_mod2','cat3_all_distMean_spTm_int_mod2',...
    'cat1_all_distMean_sp_int_mod2','cat2_all_distMean_sp_int_mod2','cat3_all_distMean_sp_int_mod2',...
    'binsMean_cat1_all_spTm','cdfMean_cat1_all_spTm_norm','binsMean_cat2_all_spTm','cdfMean_cat2_all_spTm_norm',...
    'binsMean_cat3_all_spTm','cdfMean_cat3_all_spTm_norm',...
    'binsMean_cat1_all_sp','cdfMean_cat1_all_sp_norm','binsMean_cat2_all_sp','cdfMean_cat2_all_sp_norm',...
    'binsMean_cat3_all_sp','cdfMean_cat3_all_sp_norm','nbins','start','stop',...
    'binsMean_cat1_all_spTm_mod2','cdfMean_cat1_all_spTm_norm_mod2','binsMean_cat2_all_spTm_mod2','cdfMean_cat2_all_spTm_norm_mod2',...
    'binsMean_cat3_all_spTm_mod2','cdfMean_cat3_all_spTm_norm_mod2',...
    'binsMean_cat1_all_sp_mod2','cdfMean_cat1_all_sp_norm_mod2','binsMean_cat2_all_sp_mod2','cdfMean_cat2_all_sp_norm_mod2',...
    'binsMean_cat3_all_sp_mod2','cdfMean_cat3_all_sp_norm_mod2','start_mod2',...
    '-append');

% Plotting
ymax=8; ymin=0; cdf_min=0; cdf_max=600; cdf_max_mod2=200;

cat1_pos1_mean=cat1_all_mean(:,1)+abs(min(cat1_all_mean(:,1)));
ref_plot_mean_cat1=(cat1_pos1_mean(:,1)./max(cat1_pos1_mean(:,1))).*ymax;
cat2_pos1_mean=cat2_all_mean(:,1)+abs(min(cat2_all_mean(:,1)));
ref_plot_mean_cat2=(cat2_pos1_mean(:,1)./max(cat2_pos1_mean(:,1))).*ymax;
cat3_pos1_mean=cat3_all_mean(:,1)+abs(min(cat3_all_mean(:,1)));
ref_plot_mean_cat3=(cat3_pos1_mean(:,1)./max(cat3_pos1_mean(:,1))).*ymax;

figure(6)
transp=0.05;

subplot(231)
plot(tm,cat1_all_distMean_spTm,'Color',cat(2,cat1_clr,transp),'LineWidth',0.5); hold on
plot(tm,cat2_all_distMean_spTm,'Color',cat(2,cat2_clr,transp),'LineWidth',0.5);
plot(tm,cat3_all_distMean_spTm,'Color',cat(2,cat3_clr,transp),'LineWidth',0.5);
plot(tm,nanmean(cat1_all_distMean_spTm,2),'Color',cat1_clr,'LineWidth',2);
plot(tm,nanmean(cat2_all_distMean_spTm,2),'Color',cat2_clr,'LineWidth',2);
plot(tm,nanmean(cat3_all_distMean_spTm,2),'Color',cat3_clr,'LineWidth',2);
% refs / axis
xline(tm(start),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(stop),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(start_mod2),'--','linewidth',2,'Color',[0.9 0.9 0.9 0.1])
plot(tm,ref_plot_mean_cat1,'Color',cat(2,cat1_clr,0.1),'LineWidth',3)
plot(tm,ref_plot_mean_cat2,'Color',cat(2,cat2_clr,0.1),'LineWidth',3)
plot(tm,ref_plot_mean_cat3,'Color',cat(2,cat3_clr,0.1),'LineWidth',3)
hold off
xlabel('time (sec)'); ylabel('distance to mean (mm)'); axis square
title('space-time distance to mean');
axis([tm(1),tm(end),ymin,ymax])
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);

subplot(232)
stairs(binsMean_cat1_all_spTm,cdfMean_cat1_all_spTm_norm,'Color',cat1_clr,'LineWidth',2); hold on
stairs(binsMean_cat2_all_spTm,cdfMean_cat2_all_spTm_norm,'Color',cat2_clr,'LineWidth',2);
stairs(binsMean_cat3_all_spTm,cdfMean_cat3_all_spTm_norm,'Color',cat3_clr,'LineWidth',2);
xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
set(gcf,'Position',[270 262 882 535])
axis([cdf_min cdf_max 0 1])
title('reach-grasp range')
hold off

subplot(233)
stairs(binsMean_cat1_all_spTm_mod2,cdfMean_cat1_all_spTm_norm_mod2,'Color',cat1_clr,'LineWidth',2); hold on
stairs(binsMean_cat2_all_spTm_mod2,cdfMean_cat2_all_spTm_norm_mod2,'Color',cat2_clr,'LineWidth',2);
stairs(binsMean_cat3_all_spTm_mod2,cdfMean_cat3_all_spTm_norm_mod2,'Color',cat3_clr,'LineWidth',2);
xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
set(gcf,'Position',[270 262 882 535])
axis([cdf_min cdf_max_mod2 0 1])
title('grasp range')
hold off

subplot(234)
plot(tm,cat1_all_distMean_sp,'Color',cat(2,cat1_clr,transp),'LineWidth',0.5); hold on
plot(tm,cat2_all_distMean_sp,'Color',cat(2,cat2_clr,transp),'LineWidth',0.5);
plot(tm,cat3_all_distMean_sp,'Color',cat(2,cat3_clr,transp),'LineWidth',0.5);
plot(tm,nanmean(cat1_all_distMean_sp,2),'Color',cat1_clr,'LineWidth',2);
plot(tm,nanmean(cat2_all_distMean_sp,2),'Color',cat2_clr,'LineWidth',2);
plot(tm,nanmean(cat3_all_distMean_sp,2),'Color',cat3_clr,'LineWidth',2);
% refs / axis
xline(tm(start),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(stop),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(start_mod2),'--','linewidth',2,'Color',[0.9 0.9 0.9 0.1])
plot(tm,ref_plot_mean_cat1,'Color',cat(2,cat1_clr,0.1),'LineWidth',3)
plot(tm,ref_plot_mean_cat2,'Color',cat(2,cat2_clr,0.1),'LineWidth',3)
plot(tm,ref_plot_mean_cat3,'Color',cat(2,cat3_clr,0.1),'LineWidth',3)
hold off
xlabel('time (sec)'); ylabel('distance to mean (mm)'); axis square
title('space-only distance to mean');
axis([tm(1) tm(end) ymin ymax])
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);

subplot(235)
stairs(binsMean_cat1_all_sp,cdfMean_cat1_all_sp_norm,'Color',cat1_clr,'LineWidth',2); hold on
stairs(binsMean_cat2_all_sp,cdfMean_cat2_all_sp_norm,'Color',cat2_clr,'LineWidth',2);
stairs(binsMean_cat3_all_sp,cdfMean_cat3_all_sp_norm,'Color',cat3_clr,'LineWidth',2);
xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
set(gcf,'Position',[459 262 637 535])
axis([cdf_min cdf_max 0 1])
hold off
title('reach-grasp range')

subplot(236)
stairs(binsMean_cat1_all_sp_mod2,cdfMean_cat1_all_sp_norm_mod2,'Color',cat1_clr,'LineWidth',2); hold on
stairs(binsMean_cat2_all_sp_mod2,cdfMean_cat2_all_sp_norm_mod2,'Color',cat2_clr,'LineWidth',2);
stairs(binsMean_cat3_all_sp_mod2,cdfMean_cat3_all_sp_norm_mod2,'Color',cat3_clr,'LineWidth',2);
xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
set(gcf,'Position',[459 262 637 535])
axis([cdf_min cdf_max_mod2 0 1])
hold off
title('grasp range')

set(gcf,'Position',[436 222 931 575])
if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'distMean_spaceTime_spaceOnly.png'),'png'); end
% saveas(gcf,strcat(new_folder_out,filesep,'distMean_spaceTime_spaceOnly'),'epsc');



%% Distance (space-only) to CAT2 and CAT3
% Category 1 to 3
n_timepoints=size(cat1_all,1);
n_trials_cat1_all = size(cat1_all,3);
cat1_to3_all_distMean_sp = zeros(size(size(squeeze(cat3_all(:,1,:)))));
for tt=1:n_trials_cat1_all
    for pt=1:n_timepoints
        % TO  THE MEAN
        tmp_dist_allMean=sqrt((cat1_all(pt,1,tt)-cat3_all_mean(:,1)).^2+...
            (cat1_all(pt,2,tt)-cat3_all_mean(:,2)).^2+...
            (cat1_all(pt,3,tt)-cat3_all_mean(:,3)).^2);
        [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
        cat1_to3_all_distMean_sp(pt,tt)=min_dist_mean;
    end
end
cat1_to3_all_distMean_sp_int = sum(cat1_to3_all_distMean_sp(start:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat1_to3_all_sp, binsMean_cat1_to3_all_sp] = histcounts(cat1_to3_all_distMean_sp_int,nbins);
cdfMean_cat1_to3_all_sp = cumsum(countsMean_cat1_to3_all_sp); cdfMean_cat1_to3_all_sp(end+1)=cdfMean_cat1_to3_all_sp(end);
cdfMean_cat1_to3_all_sp_norm=cdfMean_cat1_to3_all_sp./cdfMean_cat1_to3_all_sp(end);
% Range grasp
cat1_to3_all_distMean_sp_int_mod2 = sum(cat1_to3_all_distMean_sp(start_mod2:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat1_to3_all_sp_mod2, binsMean_cat1_to3_all_sp_mod2] = histcounts(cat1_to3_all_distMean_sp_int_mod2,nbins);
cdfMean_cat1_to3_all_sp_mod2 = cumsum(countsMean_cat1_to3_all_sp_mod2); cdfMean_cat1_to3_all_sp_mod2(end+1)=cdfMean_cat1_to3_all_sp_mod2(end);
cdfMean_cat1_to3_all_sp_norm_mod2=cdfMean_cat1_to3_all_sp_mod2./cdfMean_cat1_to3_all_sp_mod2(end);

% Category 2 to 3
n_trials_cat2_all = size(cat2_all,3);
cat2_to3_all_distMean_sp = zeros(size(size(squeeze(cat3_all(:,1,:)))));
for tt=1:n_trials_cat2_all
    for pt=1:n_timepoints
        % TO  THE MEAN
        tmp_dist_allMean=sqrt((cat2_all(pt,1,tt)-cat3_all_mean(:,1)).^2+...
            (cat2_all(pt,2,tt)-cat3_all_mean(:,2)).^2+...
            (cat2_all(pt,3,tt)-cat3_all_mean(:,3)).^2);
        [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
        cat2_to3_all_distMean_sp(pt,tt)=min_dist_mean;
    end
end
cat2_to3_all_distMean_sp_int = sum(cat2_to3_all_distMean_sp(start:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat2_to3_all_sp, binsMean_cat2_to3_all_sp] = histcounts(cat2_to3_all_distMean_sp_int,nbins);
cdfMean_cat2_to3_all_sp = cumsum(countsMean_cat2_to3_all_sp); cdfMean_cat2_to3_all_sp(end+1)=cdfMean_cat2_to3_all_sp(end);
cdfMean_cat2_to3_all_sp_norm=cdfMean_cat2_to3_all_sp./cdfMean_cat2_to3_all_sp(end);
% Range grasp
cat2_to3_all_distMean_sp_int_mod2 = sum(cat2_to3_all_distMean_sp(start_mod2:stop,:),1); % integral
% cumulative distribution function
[countsMean_cat2_to3_all_sp_mod2, binsMean_cat2_to3_all_sp_mod2] = histcounts(cat2_to3_all_distMean_sp_int_mod2,nbins);
cdfMean_cat2_to3_all_sp_mod2 = cumsum(countsMean_cat2_to3_all_sp_mod2); cdfMean_cat2_to3_all_sp_mod2(end+1)=cdfMean_cat2_to3_all_sp_mod2(end);
cdfMean_cat2_to3_all_sp_norm_mod2=cdfMean_cat2_to3_all_sp_mod2./cdfMean_cat2_to3_all_sp_mod2(end);

% Category 1 to 2
cat1_to2_all_distMean_sp = zeros(size(size(squeeze(cat2_all(:,1,:)))));
for tt=1:n_trials_cat1_all
    for pt=1:n_timepoints
        % TO  THE MEAN
        tmp_dist_allMean=sqrt((cat1_all(pt,1,tt)-cat2_all_mean(:,1)).^2+...
            (cat1_all(pt,2,tt)-cat2_all_mean(:,2)).^2+...
            (cat1_all(pt,3,tt)-cat2_all_mean(:,3)).^2);
        [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
        cat1_to2_all_distMean_sp(pt,tt)=min_dist_mean;
    end
end
cat1_to2_all_distMean_sp_int = sum(cat1_to2_all_distMean_sp(start:stop,:),1); % integral
% cumulative distribution function : all reach
[countsMean_cat1_to2_all_sp, binsMean_cat1_to2_all_sp] = histcounts(cat1_to2_all_distMean_sp_int,nbins);
cdfMean_cat1_to2_all_sp = cumsum(countsMean_cat1_to2_all_sp); cdfMean_cat1_to2_all_sp(end+1)=cdfMean_cat1_to2_all_sp(end);
cdfMean_cat1_to2_all_sp_norm=cdfMean_cat1_to2_all_sp./cdfMean_cat1_to2_all_sp(end);

cat1_to2_all_distMean_sp_int_mod2 = sum(cat1_to2_all_distMean_sp(start_mod2:stop,:),1); % integral
% cumulative distribution function: module 2 reach
[countsMean_cat1_to2_all_sp_mod2, binsMean_cat1_to2_all_sp_mod2] = histcounts(cat1_to2_all_distMean_sp_int_mod2,nbins);
cdfMean_cat1_to2_all_sp_mod2 = cumsum(countsMean_cat1_to2_all_sp_mod2); cdfMean_cat1_to2_all_sp_mod2(end+1)=cdfMean_cat1_to2_all_sp_mod2(end);
cdfMean_cat1_to2_all_sp_norm_mod2=cdfMean_cat1_to2_all_sp_mod2./cdfMean_cat1_to2_all_sp_mod2(end);


% Save mat
save(strcat(new_folder_mat,filesep,save_mat_name),'start','start_mod2','stop','nbins',...
    'cat1_to3_all_distMean_sp','cat1_to3_all_distMean_sp_int','binsMean_cat1_to3_all_sp','countsMean_cat1_to3_all_sp','cdfMean_cat1_to3_all_sp_norm',...
    'cat2_to3_all_distMean_sp','cat2_to3_all_distMean_sp_int','binsMean_cat2_to3_all_sp','countsMean_cat2_to3_all_sp','cdfMean_cat2_to3_all_sp_norm',...
    'cat1_to2_all_distMean_sp','cat1_to2_all_distMean_sp_int','binsMean_cat1_to2_all_sp','countsMean_cat1_to2_all_sp','cdfMean_cat1_to2_all_sp_norm',...
    'cat1_to3_all_distMean_sp_int_mod2','binsMean_cat1_to3_all_sp_mod2','countsMean_cat1_to3_all_sp_mod2','cdfMean_cat1_to3_all_sp_norm_mod2',...
    'cat2_to3_all_distMean_sp_int_mod2','binsMean_cat2_to3_all_sp_mod2','countsMean_cat2_to3_all_sp_mod2','cdfMean_cat2_to3_all_sp_norm_mod2',...
    'cat1_to2_all_distMean_sp_int_mod2','binsMean_cat1_to2_all_sp_mod2','countsMean_cat1_to2_all_sp_mod2','cdfMean_cat1_to2_all_sp_norm_mod2',...
    '-append');

% Plotting
ymax=10; ymin=0; cdf_min=0; cdf_max=800;  cdf_max_mod2=300;
transp = 0.02;
cat2_pos1_mean=cat2_all_mean(:,1)+abs(min(cat2_all_mean(:,1)));
ref_plot_mean_cat2=(cat2_pos1_mean(:,1)./max(cat2_pos1_mean(:,1))).*ymax;
cat3_pos1_mean=cat3_all_mean(:,1)+abs(min(cat3_all_mean(:,1)));
ref_plot_mean_cat3=(cat3_pos1_mean(:,1)./max(cat3_pos1_mean(:,1))).*ymax;

figure(7)
%     subplot(331)
%     plot(tm,cat3_all_distMean_sp,'Color',cat(2,cat3_clr,transp),'LineWidth',0.5); hold on
%     plot(tm,nanmean(cat3_all_distMean_sp,2),'Color',cat3_clr,'LineWidth',2);
%     plot(tm,cat1_to3_all_distMean_sp,'Color',cat(2,suc_clr,transp),'LineWidth',0.5); hold on
%     plot(tm,nanmean(cat1_to3_all_distMean_sp,2),'Color',suc_clr,'LineWidth',2);
%     % refs / axis
%     line([tm(start), tm(start)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
%     line([tm(stop), tm(stop)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
%     line([tm(start_mod2), tm(start_mod2)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2,'LineStyle','--');
%     plot(tm,ref_plot_mean_cat3,'Color',cat(2,cat3_clr,0.1),'LineWidth',3)
%     plot(tm,ref_plot_mean_cat2,'Color',cat(2,suc_clr,0.05),'LineWidth',3)
%     xlabel('time (sec)'); ylabel('distance to mean grooming (mm)'); axis square
%     title('distance SP: resting paw to grooming');
%     axis([0 0.5 ymin ymax])
%     set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
%     hold off
%     subplot(332)
%     stairs(binsMean_cat3_all_sp,cdfMean_cat3_all_sp_norm,'Color',cat(2,cat3_clr,transp),'LineWidth',2);  hold on
%     stairs(binsMean_cat1_to3_all_sp,cdfMean_cat1_to3_all_sp_norm,'Color',cat(2,suc_clr,transp),'LineWidth',2);
%     xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
%     set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
%     axis([cdf_min cdf_max 0 1])
%     hold off
%     title('reach-grasp range')
%     subplot(333)
%     stairs(binsMean_cat3_all_sp_mod2,cdfMean_cat3_all_sp_norm_mod2,'Color',cat(2,cat3_clr,transp),'LineWidth',2);  hold on
%     stairs(binsMean_cat1_to3_all_sp_mod2,cdfMean_cat1_to3_all_sp_norm_mod2,'Color',cat(2,suc_clr,transp),'LineWidth',2);
%     xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
%     set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
%     axis([cdf_min cdf_max_mod2 0 1])
%     hold off
%     title('grasp range')

subplot(231)
plot(tm,cat2_all_distMean_sp,'Color',cat(2,cat2_clr,transp),'LineWidth',0.5); hold on
plot(tm,nanmean(cat2_all_distMean_sp,2),'Color',cat2_clr,'LineWidth',2);
plot(tm,cat1_to2_all_distMean_sp,'Color',cat(2,suc_clr,transp),'LineWidth',0.5); hold on
plot(tm,nanmean(cat1_to2_all_distMean_sp,2),'Color',suc_clr,'LineWidth',2);
% refs / axis
xline(tm(start),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(stop),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(start_mod2),'--','linewidth',2,'Color',[0.9 0.9 0.9 0.1])
plot(tm,ref_plot_mean_cat2,'Color',cat(2,cat2_clr,0.1),'LineWidth',3)
plot(tm,ref_plot_mean_cat1,'Color',cat(2,suc_clr,0.05),'LineWidth',3)
xlabel('time (sec)'); ylabel('distance to mean SP: lifted (mm)'); axis square
title('distance space-only');
axis([tm(1) tm(end) ymin ymax])
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
hold off
subplot(232)
stairs(binsMean_cat2_all_sp,cdfMean_cat2_all_sp_norm,'Color',cat2_clr,'LineWidth',2);  hold on
stairs(binsMean_cat1_to2_all_sp,cdfMean_cat1_to2_all_sp_norm,'Color',suc_clr,'LineWidth',2);
xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
axis([cdf_min cdf_max 0 1])
set(gcf,'Position',[459 262 637 535])
hold off
title({'SP: resting to lifted paw';'reach-grasp range'})
subplot(233)
stairs(binsMean_cat2_all_sp_mod2,cdfMean_cat2_all_sp_norm_mod2,'Color',cat2_clr,'LineWidth',2);  hold on
stairs(binsMean_cat1_to2_all_sp_mod2,cdfMean_cat1_to2_all_sp_norm_mod2,'Color',suc_clr,'LineWidth',2);
xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
axis([cdf_min cdf_max_mod2 0 1])
set(gcf,'Position',[459 262 637 535])
hold off
title('grasp range')

subplot(234)
plot(tm,cat3_all_distMean_sp,'Color',cat(2,cat3_clr,transp),'LineWidth',0.5); hold on
plot(tm,nanmean(cat3_all_distMean_sp,2),'Color',cat3_clr,'LineWidth',2);
plot(tm,cat2_to3_all_distMean_sp,'Color',cat(2,suc_clr,transp),'LineWidth',0.5); hold on
plot(tm,nanmean(cat2_to3_all_distMean_sp,2),'Color',suc_clr,'LineWidth',2);
% refs / axis
xline(tm(start),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(stop),'linewidth',2,'Color',[0.9 0.9 0.9 0.1])
xline(tm(start_mod2),'--','linewidth',2,'Color',[0.9 0.9 0.9 0.1])
plot(tm,ref_plot_mean_cat3,'Color',cat(2,cat3_clr,0.1),'LineWidth',3)
plot(tm,ref_plot_mean_cat2,'Color',cat(2,suc_clr,0.05),'LineWidth',3)
xlabel('time (sec)'); ylabel('distance to mean grooming (mm)'); axis square
%title('distance SP: lifted paw to grooming');
axis([tm(1) tm(end) ymin ymax])
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
hold off
subplot(235)
stairs(binsMean_cat3_all_sp,cdfMean_cat3_all_sp_norm,'Color',cat3_clr,'LineWidth',2);  hold on
stairs(binsMean_cat2_to3_all_sp,cdfMean_cat2_to3_all_sp_norm,'Color',suc_clr,'LineWidth',2);
xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
axis([cdf_min cdf_max 0 1])
hold off
title('SP: lifted paw to grooming');
%title('reach-grasp range')
subplot(236)
stairs(binsMean_cat3_all_sp_mod2,cdfMean_cat3_all_sp_norm_mod2,'Color',cat3_clr,'LineWidth',2);  hold on
stairs(binsMean_cat2_to3_all_sp_mod2,cdfMean_cat2_to3_all_sp_norm_mod2,'Color',suc_clr,'LineWidth',2);
xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
axis([cdf_min cdf_max_mod2 0 1])
hold off
% title('grasp range')

set(gcf,'Position',[436 222 931 575])
if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'distMean_cat12_to23.png'),'png'); end
% saveas(gcf,strcat(new_folder_out,filesep,'distMean_cat12_to23'),'epsc');

%%  Covariance and fft
% COVARIANCE
lags=reaches.cov.lags_r_mat(:,1,1)./session.video.frame_rate;
if size(reaches.cov.cov_r_mat,2)==size(lags,1)*3
    cov_r_m = reshape(reaches.cov.cov_r_mat,[size(lags,1) 3 reaches.r_num]);
else
    cov_r_m = reaches.cov.cov_r_mat;
end

% Category 1
cat1_all_cov = cov_r_m(:,:,reaches.cat_reach==1);
cat1_all_cov_mean = nanmean(cat1_all_cov,3);
% Category 2
cat2_all_cov = cov_r_m(:,:,reaches.cat_reach==2);
cat2_all_cov_mean = nanmean(cat2_all_cov,3);
% Category 3
cat3_all_cov = cov_r_m(:,:,reaches.cat_reach==3);
cat3_all_cov_mean = nanmean(cat3_all_cov,3);
% Category 4
cat4_all_cov = cov_r_m(:,:,reaches.cat_reach==4);
cat4_all_cov_mean = nanmean(cat4_all_cov,3);

% FFT
L = size(reaches_norm,1);
freq = session.video.frame_rate*(0:(L/2))/L;
% Category 1
cat1_all_fft = reaches.fft.P1_mat(:,:,reaches.cat_reach==1);
cat1_all_fft_mean = nanmean(cat1_all_fft,3);
% Category 2
cat2_all_fft = reaches.fft.P1_mat(:,:,reaches.cat_reach==2);
cat2_all_fft_mean = nanmean(cat2_all_fft,3);
% Category 3
cat3_all_fft = reaches.fft.P1_mat(:,:,reaches.cat_reach==3);
cat3_all_fft_mean = nanmean(cat3_all_fft,3);
% Category 4
cat4_all_fft = reaches.fft.P1_mat(:,:,reaches.cat_reach==4);
cat4_all_fft_mean = nanmean(cat4_all_fft,3);

save(strcat(new_folder_mat,filesep,save_mat_name),'lags','freq',...
    'cat1_all_cov','cat1_all_cov_mean','cat2_all_cov','cat2_all_cov_mean',...
    'cat3_all_cov','cat3_all_cov_mean','cat4_all_cov','cat4_all_cov_mean',...
    'cat1_all_fft','cat1_all_fft_mean','cat2_all_fft','cat2_all_fft_mean',...
    'cat3_all_fft','cat3_all_fft_mean','cat4_all_fft','cat4_all_fft_mean',...
    '-append');

% Plotting
covX_min = -4E5;
covX_max = 4E5;
covY_min = -2E5;
covY_max = 2E5;
covZ_min = -2E5;
covZ_max = 2E5;
transp = 0.05;
fmax = 50; pmax = 30;

figure(8)
subplot(231)
plot(lags,squeeze(cat4_all_cov(:,1,:)),'Color',cat(2,cat4_clr,transp),'Linewidth',1); hold on
plot(lags,cat4_all_cov_mean(:,1),'Color',cat4_clr,'Linewidth',2);
plot(lags,squeeze(cat1_all_cov(:,1,:)),'Color',cat(2,cat1_clr,transp),'Linewidth',1); hold on
plot(lags,cat1_all_cov_mean(:,1),'Color',cat1_clr,'Linewidth',2);
plot(lags,squeeze(cat2_all_cov(:,1,:)),'Color',cat(2,cat2_clr,transp),'Linewidth',1); hold on
plot(lags,cat2_all_cov_mean(:,1),'Color',cat2_clr,'Linewidth',2);
plot(lags,squeeze(cat3_all_cov(:,1,:)),'Color',cat(2,cat3_clr,transp),'Linewidth',1); hold on
plot(lags,cat3_all_cov_mean(:,1),'Color',cat3_clr,'Linewidth',2);
hold off
xlabel('lags (sec)'); ylabel('covariance in x');
grid on;
set(gca,'linewidth',2,'box','off','GridAlpha',0.05); ylim([covX_min covX_max])
subplot(232)
plot(lags,squeeze(cat4_all_cov(:,2,:)),'Color',cat(2,cat4_clr,transp),'Linewidth',1); hold on
plot(lags,cat4_all_cov_mean(:,2),'Color',cat4_clr,'Linewidth',2);
plot(lags,squeeze(cat1_all_cov(:,2,:)),'Color',cat(2,cat1_clr,transp),'Linewidth',1); hold on
plot(lags,cat1_all_cov_mean(:,2),'Color',cat1_clr,'Linewidth',2);
plot(lags,squeeze(cat2_all_cov(:,2,:)),'Color',cat(2,cat2_clr,transp),'Linewidth',1); hold on
plot(lags,cat2_all_cov_mean(:,2),'Color',cat2_clr,'Linewidth',2);
plot(lags,squeeze(cat3_all_cov(:,2,:)),'Color',cat(2,cat3_clr,transp),'Linewidth',1); hold on
plot(lags,cat3_all_cov_mean(:,2),'Color',cat3_clr,'Linewidth',2);
hold off
xlabel('lags (sec)'); ylabel('covariance in y');
grid on;
set(gca,'linewidth',2,'box','off','GridAlpha',0.05); ylim([covY_min covY_max])
subplot(233)
plot(lags,squeeze(cat4_all_cov(:,3,:)),'Color',cat(2,cat4_clr,transp),'Linewidth',1); hold on
plot(lags,cat4_all_cov_mean(:,3),'Color',cat4_clr,'Linewidth',2);
plot(lags,squeeze(cat1_all_cov(:,3,:)),'Color',cat(2,cat1_clr,transp),'Linewidth',1); hold on
plot(lags,cat1_all_cov_mean(:,3),'Color',cat1_clr,'Linewidth',2);
plot(lags,squeeze(cat2_all_cov(:,3,:)),'Color',cat(2,cat2_clr,transp),'Linewidth',1); hold on
plot(lags,cat2_all_cov_mean(:,3),'Color',cat2_clr,'Linewidth',2);
plot(lags,squeeze(cat3_all_cov(:,3,:)),'Color',cat(2,cat3_clr,transp),'Linewidth',1); hold on
plot(lags,cat3_all_cov_mean(:,3),'Color',cat3_clr,'Linewidth',2);
hold off
xlabel('lags (sec)'); ylabel('covariance in z');
grid on;
set(gca,'linewidth',2,'box','off','GridAlpha',0.05); ylim([covZ_min covZ_max])

subplot(234)
plot(freq,squeeze(cat4_all_fft(:,1,:)),'Color',cat(2,cat4_clr,transp),'Linewidth',1); hold on
plot(freq,cat4_all_fft_mean(:,1),'Color',cat4_clr,'Linewidth',2);
plot(freq,squeeze(cat1_all_fft(:,1,:)),'Color',cat(2,cat1_clr,transp),'Linewidth',1); hold on
plot(freq,cat1_all_fft_mean(:,1),'Color',cat1_clr,'Linewidth',2);
plot(freq,squeeze(cat2_all_fft(:,1,:)),'Color',cat(2,cat2_clr,transp),'Linewidth',1); hold on
plot(freq,cat2_all_fft_mean(:,1),'Color',cat2_clr,'Linewidth',2);
plot(freq,squeeze(cat3_all_fft(:,1,:)),'Color',cat(2,cat3_clr,transp),'Linewidth',1); hold on
plot(freq,cat3_all_fft_mean(:,1),'Color',cat3_clr,'Linewidth',2);
hold off
xlabel('frequency (Hz)'); ylabel('|P1(f)| in x');
grid on; axis([0 fmax 0 pmax])
set(gca,'linewidth',2,'box','off','GridAlpha',0.05);
subplot(235)
plot(freq,squeeze(cat4_all_fft(:,2,:)),'Color',cat(2,cat4_clr,transp),'Linewidth',1); hold on
plot(freq,cat4_all_fft_mean(:,2),'Color',cat4_clr,'Linewidth',2);
plot(freq,squeeze(cat1_all_fft(:,2,:)),'Color',cat(2,cat1_clr,transp),'Linewidth',1); hold on
plot(freq,cat1_all_fft_mean(:,2),'Color',cat1_clr,'Linewidth',2);
plot(freq,squeeze(cat2_all_fft(:,2,:)),'Color',cat(2,cat2_clr,transp),'Linewidth',1); hold on
plot(freq,cat2_all_fft_mean(:,2),'Color',cat2_clr,'Linewidth',2);
plot(freq,squeeze(cat3_all_fft(:,2,:)),'Color',cat(2,cat3_clr,transp),'Linewidth',1); hold on
plot(freq,cat3_all_fft_mean(:,2),'Color',cat3_clr,'Linewidth',2);
hold off
xlabel('frequency (Hz)'); ylabel('|P1(f)| in y');
grid on; axis([0 fmax 0 pmax])
set(gca,'linewidth',2,'box','off','GridAlpha',0.05);
subplot(236)
plot(freq,squeeze(cat4_all_fft(:,3,:)),'Color',cat(2,cat4_clr,transp),'Linewidth',1); hold on
plot(freq,cat4_all_fft_mean(:,3),'Color',cat4_clr,'Linewidth',2);
plot(freq,squeeze(cat1_all_fft(:,3,:)),'Color',cat(2,cat1_clr,transp),'Linewidth',1); hold on
plot(freq,cat1_all_fft_mean(:,3),'Color',cat1_clr,'Linewidth',2);
plot(freq,squeeze(cat2_all_fft(:,3,:)),'Color',cat(2,cat2_clr,transp),'Linewidth',1); hold on
plot(freq,cat2_all_fft_mean(:,3),'Color',cat2_clr,'Linewidth',2);
plot(freq,squeeze(cat3_all_fft(:,3,:)),'Color',cat(2,cat3_clr,transp),'Linewidth',1); hold on
plot(freq,cat3_all_fft_mean(:,1),'Color',cat3_clr,'Linewidth',2);
hold off
xlabel('frequency (Hz)'); ylabel('|P1(f)| in z');
grid on; axis([0 fmax 0 pmax])
set(gca,'linewidth',2,'box','off','GridAlpha',0.05);
% 
set(gcf,'Position',[102, 255 1246 542]);

if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'cov_fft_all.png'),'png'); end
%saveas(gcf,strcat(new_folder_out,filesep,'cov_fft_all'),'epsc');



%% Show success vs miss
    new_folder_sucMiss = strcat(new_folder_out,filesep,'succ_vs_miss');
    mkdir(new_folder_sucMiss)

    % Category 1
    % success
    cat1_suc_loc = intersect(find(reaches.success_reach==1),find(reaches.cat_reach==1));
    cat1_suc = reaches_norm(:,:,cat1_suc_loc);
    cat1_suc_mean = nanmean(cat1_suc,3);
    % fail & miss
    cat1_fail_loc = intersect(find(reaches.success_reach==0),find(reaches.cat_reach==1));
    cat1_miss_loc = intersect(find(reaches.hit_reach==0),find(reaches.cat_reach==1));
    cat1_fail_miss_loc = cat(1,cat1_fail_loc,cat1_miss_loc);clc
    cat1_fail_miss= reaches_norm(:,:,cat1_fail_miss_loc);
    cat1_fail_miss_mean = nanmean(cat1_fail_miss,3);

    % Category 2
    cat2_suc_loc = intersect(find(reaches.success_reach==1),find(reaches.cat_reach==2));
    cat2_suc = reaches_norm(:,:,cat2_suc_loc);
    cat2_suc_mean = nanmean(cat2_suc,3);
    % fail & miss
    cat2_fail_loc = intersect(find(reaches.success_reach==0),find(reaches.cat_reach==2));
    cat2_miss_loc = intersect(find(reaches.hit_reach==0),find(reaches.cat_reach==2));
    cat2_fail_miss_loc = cat(1,cat2_fail_loc,cat2_miss_loc);
    cat2_fail_miss= reaches_norm(:,:,cat2_fail_miss_loc);
    cat2_fail_miss_mean = nanmean(cat2_fail_miss,3);

%     % Save mat
%     save(strcat(new_folder_mat,filesep,save_mat_name),'cat1_suc_loc','cat1_suc','cat1_suc_mean',...
%         'cat1_fail_miss_loc','cat1_fail_miss','cat1_fail_miss_mean','cat2_suc_loc','cat2_suc','cat2_suc_mean',...
%         'cat2_fail_miss_loc','cat2_fail_miss','cat2_fail_miss_mean',...
%         'suc_clr','miss_clr','-append');

    % Plotting
    transp = 0.02;
    transp_mean_proj = 1;
    mean_p_lw=1.5;
    mean_r_lw=2.5;
    trial_lw = 1;
    %x_min=-16; x_max=2; y_min=-10; y_max=4; z_min=-8; z_max=7;
    oneMat = ones(size(cat1_suc,1),1);

    figure(9)
    subplot(121)
    plot3(cat1_fail_miss_mean(:,1),cat1_fail_miss_mean(:,2),cat1_fail_miss_mean(:,3),'LineWidth',mean_r_lw,'Color',miss_clr); hold on
    plot3(cat1_suc_mean(:,1),cat1_suc_mean(:,2),cat1_suc_mean(:,3),'LineWidth',mean_r_lw,'Color',suc_clr);
    grid on
    view(-45,18)
    title('SP: resting paw')
    set(gca,...
        'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1,...
        'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top');
    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
    axis ([x_min x_max y_min y_max z_min z_max])
    yL = get(gca,'YLim');
    zL = get(gca,'ZLim');
    plot3(squeeze(cat1_fail_miss(:,1,:)), oneMat .* yL(2), squeeze(cat1_fail_miss(:,3,:)),'-','Color',cat(2,miss_clr,transp),'LineWidth',trial_lw);
    plot3(squeeze(cat1_fail_miss(:,1,:)), squeeze(cat1_fail_miss(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,miss_clr,transp),'LineWidth',trial_lw);
    plot3(squeeze(cat1_suc(:,1,:)), oneMat .* yL(2), squeeze(cat1_suc(:,3,:)),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
    plot3(squeeze(cat1_suc(:,1,:)), squeeze(cat1_suc(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);

    plot3(cat1_fail_miss_mean(:,1), oneMat .* yL(2), cat1_fail_miss_mean(:,3),'-','Color',cat(2,miss_clr,transp_mean_proj),'LineWidth',mean_p_lw);
    plot3(cat1_fail_miss_mean(:,1), cat1_fail_miss_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,miss_clr,transp_mean_proj),'LineWidth',mean_p_lw);
    plot3(cat1_suc_mean(:,1), oneMat .* yL(2), cat1_suc_mean(:,3),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
    plot3(cat1_suc_mean(:,1), cat1_suc_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
    hold off

    subplot(122)
    plot3(cat2_fail_miss_mean(:,1),cat2_fail_miss_mean(:,2),cat2_fail_miss_mean(:,3),'LineWidth',mean_r_lw,'Color',miss_clr); hold on
    plot3(cat2_suc_mean(:,1),cat2_suc_mean(:,2),cat2_suc_mean(:,3),'LineWidth',mean_r_lw,'Color',suc_clr);
    grid on
    view(-45,18)
    title('SP: lifted paw')
    set(gca,...
        'plotboxaspectratio',[1,1,1],'ticklength',[1,1]*.025,'linewidth',1,...
        'fontsize',12, 'nextplot','add','tickdir','out','box','off','layer','top');
    axis ([x_min x_max y_min y_max z_min z_max])
    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
    yL = get(gca,'YLim');
    zL = get(gca,'ZLim');
    plot3(squeeze(cat2_suc(:,1,:)), oneMat .* yL(2), squeeze(cat2_suc(:,3,:)),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
    plot3(squeeze(cat2_suc(:,1,:)), squeeze(cat2_suc(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp),'LineWidth',trial_lw);
    plot3(squeeze(cat2_fail_miss(:,1,:)), oneMat .* yL(2), squeeze(cat2_fail_miss(:,3,:)),'-','Color',cat(2,miss_clr,transp),'LineWidth',trial_lw);
    plot3(squeeze(cat2_fail_miss(:,1,:)), squeeze(cat2_fail_miss(:,2,:)), oneMat .* zL(1),'-','Color',cat(2,miss_clr,transp),'LineWidth',trial_lw);

    plot3(cat2_suc_mean(:,1), oneMat .* yL(2), cat2_suc_mean(:,3),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
    plot3(cat2_suc_mean(:,1), cat2_suc_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,suc_clr,transp_mean_proj),'LineWidth',mean_p_lw);
    plot3(cat2_fail_miss_mean(:,1), oneMat .* yL(2), cat2_fail_miss_mean(:,3),'-','Color',cat(2,miss_clr,transp_mean_proj),'LineWidth',mean_p_lw);
    plot3(cat2_fail_miss_mean(:,1), cat2_fail_miss_mean(:,2), oneMat .* zL(1),'-','Color',cat(2,miss_clr,transp_mean_proj),'LineWidth',mean_p_lw);
    hold off

      set(gcf,'Position',[2255         325         713         409])
    if save_flag, saveas(gcf,strcat(new_folder_sucMiss,filesep,'suc_vs_miss_3d_cat12.png'),'png'); end
    %saveas(gcf,strcat(new_folder_out,filesep,'xyz_suc_miss_3d_cat123_v2'),'epsc');


    %% Distance and velocity
    % Category 1
    % Success
    cat1_suc_dist = distance_mm(:,cat1_suc_loc);
    cat1_suc_dist_mean = nanmean(cat1_suc_dist,2);
    cat1_suc_vel = velocity_mm(:,cat1_suc_loc);
    cat1_suc_vel_mean = nanmean(cat1_suc_vel,2);
    % Fail and miss
    cat1_fail_miss_dist = distance_mm(:,cat1_fail_miss_loc);
    cat1_fail_miss_dist_mean = nanmean(cat1_fail_miss_dist,2);
    cat1_fail_miss_vel = velocity_mm(:,cat1_fail_miss_loc);
    cat1_fail_miss_vel_mean = nanmean(cat1_fail_miss_vel,2);

    % Category 2 (success)
    cat2_suc_dist = distance_mm(:,cat2_suc_loc);
    cat2_suc_dist_mean = nanmean(cat2_suc_dist,2);
    cat2_suc_vel = velocity_mm(:,cat2_suc_loc);
    cat2_suc_vel_mean = nanmean(cat2_suc_vel,2);
    % Fail and miss
    cat2_fail_miss_dist = distance_mm(:,cat2_fail_miss_loc);
    cat2_fail_miss_dist_mean = nanmean(cat2_fail_miss_dist,2);
    cat2_fail_miss_vel = velocity_mm(:,cat2_fail_miss_loc);
    cat2_fail_miss_vel_mean = nanmean(cat2_fail_miss_vel,2);

% 
% 
%     % Save mat
%     save(strcat(new_folder_mat,filesep,save_mat_name),'cat1_suc_dist','cat1_suc_dist_mean','cat1_suc_vel','cat1_suc_vel_mean',...
%         'cat1_fail_miss_dist','cat1_fail_miss_dist_mean','cat1_fail_miss_vel','cat1_fail_miss_vel_mean',...
%         'cat2_suc_dist','cat2_suc_dist_mean','cat2_suc_vel','cat2_suc_vel_mean',...
%         'cat2_fail_miss_dist','cat2_fail_miss_dist_mean','cat2_fail_miss_vel','cat2_fail_miss_vel_mean',...
%         '-append');

    % Plot
    figure(10)
    transp=0.08;
    subplot(221)
    plot(tm,cat1_fail_miss_dist,'Color',cat(2,miss_clr,transp)); hold on
    plot(tm,cat1_fail_miss_dist_mean,'Color',miss_clr,'Linewidth',2);
    plot(tm,cat1_suc_dist,'Color',cat(2,suc_clr,transp)); hold on
    plot(tm,cat1_suc_dist_mean,'Color',suc_clr,'Linewidth',2);
    hold off
    grid on; xlabel('time (sec)'); ylabel('distance (mm)')
    title('SP: paw resting'); xlim([tm(1),tm(end)]); axis square
    set(gca,'linewidth',2,'box','off','GridAlpha',0.05);

    subplot(222)
    plot(tm,cat2_fail_miss_dist,'Color',cat(2,miss_clr,transp)); hold on
    plot(tm,cat2_fail_miss_dist_mean,'Color',miss_clr,'Linewidth',2);
    plot(tm,cat2_suc_dist,'Color',cat(2,suc_clr,transp)); hold on
    plot(tm,cat2_suc_dist_mean,'Color',suc_clr,'Linewidth',2);
    hold off
    grid on; xlabel('time (sec)'); ylabel('distance (mm)')
    title('SP: paw lifted'); xlim([tm(1),tm(end)]); axis square
    set(gca,'linewidth',2,'box','off','GridAlpha',0.05);

    transp=0.02;
    subplot(223)
    plot(tm_vel,cat1_fail_miss_vel,'Color',cat(2,miss_clr,transp)); hold on
    plot(tm_vel,cat1_fail_miss_vel_mean,'Color',miss_clr,'Linewidth',2);
    plot(tm_vel,cat1_suc_vel,'Color',cat(2,suc_clr,transp)); hold on
    plot(tm_vel,cat1_suc_vel_mean,'Color',suc_clr,'Linewidth',2);
    hold off
    grid on; xlabel('time (sec)'); ylabel('velocity (mm/s)')
    xlim([tm(1),tm(end)]); axis square; ylim([-.5 .5])
    set(gca,'linewidth',2,'box','off','GridAlpha',0.05);

    subplot(224)
    plot(tm_vel,cat2_fail_miss_vel,'Color',cat(2,miss_clr,transp)); hold on
    plot(tm_vel,cat2_fail_miss_vel_mean,'Color',miss_clr,'Linewidth',2);
    plot(tm_vel,cat2_suc_vel,'Color',cat(2,suc_clr,transp)); hold on
    plot(tm_vel,cat2_suc_vel_mean,'Color',suc_clr,'Linewidth',2);
    hold off
    grid on; xlabel('time (sec)'); ylabel('velocity (mm/s)')
    xlim([tm(1),tm(end)]); axis square; ylim([-.5 .5])
    set(gca,'linewidth',2,'box','off','GridAlpha',0.05);

    set(gcf, 'Position', [241 232 923 565])
    if save_flag, saveas(gcf,strcat(new_folder_sucMiss,filesep,'dist_vel'),'png'); end
    % saveas(gcf,strcat(new_folder_out,filesep,'dist_vel'),'epsc');


 
    %% Distance to mean: space & time
    start = 80;
    stop = 180;
    nbins = 40;

    % Category 1
    % Success
    cat1_suc_distMean_spTm= sqrt((squeeze(cat1_suc(:,1,:))-squeeze(cat1_suc_mean(:,1))).^2 + ...
        (squeeze(cat1_suc(:,2,:))-squeeze(cat1_suc_mean(:,2))).^2 + ...
        (squeeze(cat1_suc(:,3,:))-squeeze(cat1_suc_mean(:,3))).^2);
    cat1_suc_distMean_spTm_int = sum(cat1_suc_distMean_spTm(start:stop,:),1); % integral
    % cumulative distribution function
    [countsMean_cat1_suc_spTm, binsMean_cat1_suc_spTm] = histcounts(cat1_suc_distMean_spTm_int,nbins);
    cdfMean_cat1_suc_spTm = cumsum(countsMean_cat1_suc_spTm); cdfMean_cat1_suc_spTm(end+1)=cdfMean_cat1_suc_spTm(end);
    cdfMean_cat1_suc_spTm_norm=cdfMean_cat1_suc_spTm./cdfMean_cat1_suc_spTm(end);

    % Fail/Miss
    cat1_fail_miss_distMean_spTm= sqrt((squeeze(cat1_fail_miss(:,1,:))-squeeze(cat1_fail_miss_mean(:,1))).^2 + ...
        (squeeze(cat1_fail_miss(:,2,:))-squeeze(cat1_fail_miss_mean(:,2))).^2 + ...
        (squeeze(cat1_fail_miss(:,3,:))-squeeze(cat1_fail_miss_mean(:,3))).^2);
    cat1_fail_miss_distMean_spTm_int = sum(cat1_fail_miss_distMean_spTm(start:stop,:),1); % integral
    % cumulative distribution function
    [countsMean_cat1_fail_miss_spTm, binsMean_cat1_fail_miss_spTm] = histcounts(cat1_fail_miss_distMean_spTm_int,nbins);
    cdfMean_cat1_fail_miss_spTm = cumsum(countsMean_cat1_fail_miss_spTm); cdfMean_cat1_fail_miss_spTm(end+1)=cdfMean_cat1_fail_miss_spTm(end);
    cdfMean_cat1_fail_miss_spTm_norm=cdfMean_cat1_fail_miss_spTm./cdfMean_cat1_fail_miss_spTm(end);

    % Significance
    [cat1_pMean_spTm,cat1_hMean_spTm] = ranksum(cat1_suc_distMean_spTm_int,cat1_fail_miss_distMean_spTm_int);

    % Category 2
    % Success
    cat2_suc_distMean_spTm= sqrt((squeeze(cat2_suc(:,1,:))-squeeze(cat2_suc_mean(:,1))).^2 + ...
        (squeeze(cat2_suc(:,2,:))-squeeze(cat2_suc_mean(:,2))).^2 + ...
        (squeeze(cat2_suc(:,3,:))-squeeze(cat2_suc_mean(:,3))).^2);
    cat2_suc_distMean_spTm_int = sum(cat2_suc_distMean_spTm(start:stop,:),1); % integral
    % cumulative distribution function
    [countsMean_cat2_suc_spTm, binsMean_cat2_suc_spTm] = histcounts(cat2_suc_distMean_spTm_int,nbins);
    cdfMean_cat2_suc_spTm = cumsum(countsMean_cat2_suc_spTm); cdfMean_cat2_suc_spTm(end+1)=cdfMean_cat2_suc_spTm(end);
    cdfMean_cat2_suc_spTm_norm=cdfMean_cat2_suc_spTm./cdfMean_cat2_suc_spTm(end);

    % Fail/Miss
    cat2_fail_miss_distMean_spTm= sqrt((squeeze(cat2_fail_miss(:,1,:))-squeeze(cat2_fail_miss_mean(:,1))).^2 + ...
        (squeeze(cat2_fail_miss(:,2,:))-squeeze(cat2_fail_miss_mean(:,2))).^2 + ...
        (squeeze(cat2_fail_miss(:,3,:))-squeeze(cat2_fail_miss_mean(:,3))).^2);
    cat2_fail_miss_distMean_spTm_int = sum(cat2_fail_miss_distMean_spTm(start:stop,:),1); % integral
    % cumulative distribution function
    [countsMean_cat2_fail_miss_spTm, binsMean_cat2_fail_miss_spTm] = histcounts(cat2_fail_miss_distMean_spTm_int,nbins);
    cdfMean_cat2_fail_miss_spTm = cumsum(countsMean_cat2_fail_miss_spTm); cdfMean_cat2_fail_miss_spTm(end+1)=cdfMean_cat2_fail_miss_spTm(end);
    cdfMean_cat2_fail_miss_spTm_norm=cdfMean_cat2_fail_miss_spTm./cdfMean_cat2_fail_miss_spTm(end);

    % Significance
    [cat2_pMean_spTm,cat2_hMean_spTm] = ranksum(cat2_suc_distMean_spTm_int,cat2_fail_miss_distMean_spTm_int);

    % Category 3
    % all
    cat3_all_distMean_spTm= sqrt((squeeze(cat3_all(:,1,:))-squeeze(cat3_all_mean(:,1))).^2 + ...
        (squeeze(cat3_all(:,2,:))-squeeze(cat3_all_mean(:,2))).^2 + ...
        (squeeze(cat3_all(:,3,:))-squeeze(cat3_all_mean(:,3))).^2);
    cat3_all_distMean_spTm_int = sum(cat3_all_distMean_spTm(start:stop,:),1); % integral
    % cumulative distribution function
    [countsMean_cat3_all_spTm, binsMean_cat3_all_spTm] = histcounts(cat3_all_distMean_spTm_int,nbins);
    cdfMean_cat3_all_spTm = cumsum(countsMean_cat3_all_spTm); cdfMean_cat3_all_spTm(end+1)=cdfMean_cat3_all_spTm(end);
    cdfMean_cat3_all_spTm_norm=cdfMean_cat3_all_spTm./cdfMean_cat3_all_spTm(end);

    % Save mat
%     save(strcat(new_folder_mat,filesep,save_mat_name),'start','stop','nbins',...
%         'cat1_suc_distMean_spTm','cat1_suc_distMean_spTm_int','binsMean_cat1_suc_spTm','countsMean_cat1_suc_spTm','cdfMean_cat1_suc_spTm_norm',...
%         'cat1_fail_miss_distMean_spTm','cat1_fail_miss_distMean_spTm_int','binsMean_cat1_fail_miss_spTm','countsMean_cat1_fail_miss_spTm','cdfMean_cat1_fail_miss_spTm_norm',...
%         'cat2_suc_distMean_spTm','cat2_suc_distMean_spTm_int','binsMean_cat2_suc_spTm','countsMean_cat2_suc_spTm','cdfMean_cat2_suc_spTm_norm',...
%         'cat2_fail_miss_distMean_spTm','cat2_fail_miss_distMean_spTm_int','binsMean_cat2_fail_miss_spTm','countsMean_cat2_fail_miss_spTm','cdfMean_cat2_fail_miss_spTm_norm',...
%         'cat3_all_distMean_spTm','cat3_all_distMean_spTm_int','countsMean_cat3_all_spTm','binsMean_cat3_all_spTm','cdfMean_cat3_all_spTm_norm',...
%         'cat1_pMean_spTm','cat1_hMean_spTm','cat2_pMean_spTm','cat2_hMean_spTm','-append');

    % Plot
    ymax=10; ymin=0; cdf_min=0; cdf_max=600;
    cat1_pos1_mean=cat1_suc_mean(:,1)+abs(min(cat1_suc_mean(:,1)));
    ref_plot_mean_cat1=(cat1_pos1_mean(:,1)./max(cat1_pos1_mean(:,1))).*ymax;
    cat2_pos1_mean=cat2_suc_mean(:,1)+abs(min(cat2_suc_mean(:,1)));
    ref_plot_mean_cat2=(cat2_pos1_mean(:,1)./max(cat2_pos1_mean(:,1))).*ymax;
    cat3_pos1_mean=cat3_all_mean(:,1)+abs(min(cat3_all_mean(:,1)));
    ref_plot_mean_cat3=(cat3_pos1_mean(:,1)./max(cat3_pos1_mean(:,1))).*ymax;
    transp = 0.05;

    figure(13)
    subplot(231)
    plot(tm,cat1_fail_miss_distMean_spTm,'Color',cat(2,miss_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat1_fail_miss_distMean_spTm,2),'Color',miss_clr,'LineWidth',2);
    plot(tm,cat1_suc_distMean_spTm,'Color',cat(2,suc_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat1_suc_distMean_spTm,2),'Color',suc_clr,'LineWidth',2);
    % refs / axis
    line([tm(start), tm(start)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    line([tm(stop), tm(stop)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    plot(tm,ref_plot_mean_cat1,'Color',[0 0 0 0.1],'LineWidth',3)
    xlabel('time (sec)'); ylabel('distance to mean (mm)'); axis square
    title('SP: resting paw');
    axis([tm(1) tm(end) ymin ymax])
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);

    subplot(232)
    plot(tm,cat2_fail_miss_distMean_spTm,'Color',cat(2,miss_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat2_fail_miss_distMean_spTm,2),'Color',miss_clr,'LineWidth',2);
    plot(tm,cat2_suc_distMean_spTm,'Color',cat(2,suc_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat2_suc_distMean_spTm,2),'Color',suc_clr,'LineWidth',2);
    % refs / axis
    line([tm(start), tm(start)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    line([tm(stop), tm(stop)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    plot(tm,ref_plot_mean_cat2,'Color',[0 0 0 0.1],'LineWidth',3)
    xlabel('time (sec)'); ylabel('distance to mean (mm)'); axis square
    title({'Space-time distance';'SP: lifted paw'});
    axis([tm(1) tm(end) ymin ymax])
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);

    subplot(233)
    plot(tm,cat3_all_distMean_spTm,'Color',cat(2,suc_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat3_all_distMean_spTm,2),'Color',suc_clr,'LineWidth',2);
    % refs / axis
    line([tm(start), tm(start)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    line([tm(stop), tm(stop)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    plot(tm,ref_plot_mean_cat3,'Color',[0 0 0 0.1],'LineWidth',3)
    xlabel('time (sec)'); ylabel('distance to mean (mm)'); axis square
    title('grooming');
    axis([tm(1) tm(end) ymin ymax])
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);

    subplot(234)
    stairs(binsMean_cat1_fail_miss_spTm,cdfMean_cat1_fail_miss_spTm_norm,'Color',miss_clr,'LineWidth',2); hold on
    stairs(binsMean_cat1_suc_spTm,cdfMean_cat1_suc_spTm_norm,'Color',suc_clr,'LineWidth',2);
    xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
    title(sprintf('%s%.4f','p-value=',cat1_pMean_spTm))
    axis([cdf_min cdf_max 0 1])

    subplot(235)
    stairs(binsMean_cat2_fail_miss_spTm,cdfMean_cat2_fail_miss_spTm_norm,'Color',miss_clr,'LineWidth',2); hold on
    stairs(binsMean_cat2_suc_spTm,cdfMean_cat2_suc_spTm_norm,'Color',suc_clr,'LineWidth',2);
    xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
    %title(sprintf('%s%.4f','p-value=',cat2_pMean_spTm))
    axis([cdf_min cdf_max 0 1])

    subplot(236)
    stairs(binsMean_cat3_all_spTm,cdfMean_cat3_all_spTm_norm,'Color',suc_clr,'LineWidth',2);
    xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
    set(gcf,'Position',[118 262 882 535])
    %title(sprintf('%s%.4f','p-value=',cat2_pMean_spTm))
    axis([cdf_min cdf_max 0 1])

    if save_flag, saveas(gcf,strcat(new_folder_sucMiss,filesep,'distMean_spTm.png'),'png'); end
    % saveas(gcf,strcat(new_folder_out,filesep,'distMean_spTm'),'epsc');

    %% Distance to mean: space-only
    n_timepoints = size(cat1_suc,1);

    % Category 1
    % Success
    n_trials_cat1_suc = size(cat1_suc,3);
    cat1_suc_distMean_sp = zeros(size(cat1_suc_distMean_spTm));
    for tt=1:n_trials_cat1_suc
        for pt=1:n_timepoints
            % TO  THE MEAN
            tmp_dist_allMean=sqrt((cat1_suc(pt,1,tt)-cat1_suc_mean(:,1)).^2+...
                (cat1_suc(pt,2,tt)-cat1_suc_mean(:,2)).^2+...
                (cat1_suc(pt,3,tt)-cat1_suc_mean(:,3)).^2);
            [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
            cat1_suc_distMean_sp(pt,tt)=min_dist_mean;
        end
    end
    cat1_suc_distMean_sp_int = sum(cat1_suc_distMean_sp(start:stop,:),1); % integral
    % cumulative distribution function
    [countsMean_cat1_suc_sp, binsMean_cat1_suc_sp] = histcounts(cat1_suc_distMean_sp_int,nbins);
    cdfMean_cat1_suc_sp = cumsum(countsMean_cat1_suc_sp); cdfMean_cat1_suc_sp(end+1)=cdfMean_cat1_suc_sp(end);
    cdfMean_cat1_suc_sp_norm=cdfMean_cat1_suc_sp./cdfMean_cat1_suc_sp(end);

    % Fail or miss
    n_trials_cat1_fail_miss = size(cat1_fail_miss,3);
    cat1_fail_miss_distMean_sp = zeros(size(cat1_fail_miss_distMean_spTm));
    for tt=1:n_trials_cat1_fail_miss
        for pt=1:n_timepoints
            % TO  THE MEAN
            tmp_dist_allMean=sqrt((cat1_fail_miss(pt,1,tt)-cat1_fail_miss_mean(:,1)).^2+...
                (cat1_fail_miss(pt,2,tt)-cat1_fail_miss_mean(:,2)).^2+...
                (cat1_fail_miss(pt,3,tt)-cat1_fail_miss_mean(:,3)).^2);
            [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
            cat1_fail_miss_distMean_sp(pt,tt)=min_dist_mean;
        end
    end
    cat1_fail_miss_distMean_sp_int = sum(cat1_fail_miss_distMean_sp(start:stop,:),1); % integral
    % cumulative distribution function
    
    [countsMean_cat1_fail_miss_sp, binsMean_cat1_fail_miss_sp] = histcounts(cat1_fail_miss_distMean_sp_int,nbins);
    cdfMean_cat1_fail_miss_sp = cumsum(countsMean_cat1_fail_miss_sp); cdfMean_cat1_fail_miss_sp(end+1)=cdfMean_cat1_fail_miss_sp(end);
    cdfMean_cat1_fail_miss_sp_norm=cdfMean_cat1_fail_miss_sp./cdfMean_cat1_fail_miss_sp(end);

    % Significance
    [cat1_pMean_sp,cat1_hMean_sp] = ranksum(cat1_suc_distMean_sp_int,cat1_fail_miss_distMean_sp_int);


    % Category 2
    % Success
    n_trials_cat2_suc = size(cat2_suc,3);
    cat2_suc_distMean_sp = zeros(size(cat2_suc_distMean_spTm));
    for tt=1:n_trials_cat2_suc
        for pt=1:n_timepoints
            % TO  THE MEAN
            tmp_dist_allMean=sqrt((cat2_suc(pt,1,tt)-cat2_suc_mean(:,1)).^2+...
                (cat2_suc(pt,2,tt)-cat2_suc_mean(:,2)).^2+...
                (cat2_suc(pt,3,tt)-cat2_suc_mean(:,3)).^2);
            [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
            cat2_suc_distMean_sp(pt,tt)=min_dist_mean;
        end
    end
    cat2_suc_distMean_sp_int = sum(cat2_suc_distMean_sp(start:stop,:),1); % integral
    % cumulative distribution function
    [countsMean_cat2_suc_sp, binsMean_cat2_suc_sp] = histcounts(cat2_suc_distMean_sp_int,nbins);
    cdfMean_cat2_suc_sp = cumsum(countsMean_cat2_suc_sp); cdfMean_cat2_suc_sp(end+1)=cdfMean_cat2_suc_sp(end);
    cdfMean_cat2_suc_sp_norm=cdfMean_cat2_suc_sp./cdfMean_cat2_suc_sp(end);

    % Fail or miss
    n_trials_cat2_fail_miss = size(cat2_fail_miss,3);
    cat2_fail_miss_distMean_sp = zeros(size(cat2_fail_miss_distMean_spTm));
    for tt=1:n_trials_cat2_fail_miss
        for pt=1:n_timepoints
            % TO  THE MEAN
            tmp_dist_allMean=sqrt((cat2_fail_miss(pt,1,tt)-cat2_fail_miss_mean(:,1)).^2+...
                (cat2_fail_miss(pt,2,tt)-cat2_fail_miss_mean(:,2)).^2+...
                (cat2_fail_miss(pt,3,tt)-cat2_fail_miss_mean(:,3)).^2);
            [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
            cat2_fail_miss_distMean_sp(pt,tt)=min_dist_mean;
        end
    end
    cat2_fail_miss_distMean_sp_int = sum(cat2_fail_miss_distMean_sp(start:stop,:),1); % integral
    % cumulative distribution function
    [countsMean_cat2_fail_miss_sp, binsMean_cat2_fail_miss_sp] = histcounts(cat2_fail_miss_distMean_sp_int,nbins);
    cdfMean_cat2_fail_miss_sp = cumsum(countsMean_cat2_fail_miss_sp); cdfMean_cat2_fail_miss_sp(end+1)=cdfMean_cat2_fail_miss_sp(end);
    cdfMean_cat2_fail_miss_sp_norm=cdfMean_cat2_fail_miss_sp./cdfMean_cat2_fail_miss_sp(end);

    % Significance
    [cat2_pMean_sp,cat2_hMean_sp] = ranksum(cat2_suc_distMean_sp_int,cat2_fail_miss_distMean_sp_int);


    % Category 3
    n_trials_cat3_all = size(cat3_all,3);
    cat3_all_distMean_sp = zeros(size(cat3_all_distMean_spTm));
    for tt=1:n_trials_cat3_all
        for pt=1:n_timepoints
            % TO  THE MEAN
            tmp_dist_allMean=sqrt((cat3_all(pt,1,tt)-cat3_all_mean(:,1)).^2+...
                (cat3_all(pt,2,tt)-cat3_all_mean(:,2)).^2+...
                (cat3_all(pt,3,tt)-cat3_all_mean(:,3)).^2);
            [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
            cat3_all_distMean_sp(pt,tt)=min_dist_mean;
        end
    end
    cat3_all_distMean_sp_int = sum(cat3_all_distMean_sp(start:stop,:),1); % integral
    % cumulative distribution function
    [countsMean_cat3_all_sp, binsMean_cat3_all_sp] = histcounts(cat3_all_distMean_sp_int,nbins);
    cdfMean_cat3_all_sp = cumsum(countsMean_cat3_all_sp); cdfMean_cat3_all_sp(end+1)=cdfMean_cat3_all_sp(end);
    cdfMean_cat3_all_sp_norm=cdfMean_cat3_all_sp./cdfMean_cat3_all_sp(end);

%     % Save mat
%     save(strcat(new_folder_mat,filesep,save_mat_name),...
%         'cat1_suc_distMean_sp','cat1_suc_distMean_sp_int','binsMean_cat1_suc_sp','countsMean_cat1_suc_sp','cdfMean_cat1_suc_sp_norm',...
%         'cat1_fail_miss_distMean_sp','cat1_fail_miss_distMean_sp_int','binsMean_cat1_fail_miss_sp','countsMean_cat1_fail_miss_sp','cdfMean_cat1_fail_miss_sp_norm',...
%         'cat2_suc_distMean_sp','cat2_suc_distMean_sp_int','binsMean_cat2_suc_sp','countsMean_cat2_suc_sp','cdfMean_cat2_suc_sp_norm',...
%         'cat2_fail_miss_distMean_sp','cat2_fail_miss_distMean_sp_int','binsMean_cat2_fail_miss_sp','countsMean_cat2_fail_miss_sp','cdfMean_cat2_fail_miss_sp_norm',...
%         'cat3_all_distMean_sp','cat3_all_distMean_sp_int','countsMean_cat3_all_sp','binsMean_cat3_all_sp','cdfMean_cat3_all_sp_norm',...
%         'cat1_pMean_sp','cat1_hMean_sp','cat2_pMean_sp','cat2_hMean_sp','-append');


    % Plot
    figure(14)
    subplot(231)
    plot(tm,cat1_fail_miss_distMean_sp,'Color',cat(2,miss_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat1_fail_miss_distMean_sp,2),'Color',miss_clr,'LineWidth',2);
    plot(tm,cat1_suc_distMean_sp,'Color',cat(2,suc_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat1_suc_distMean_sp,2),'Color',suc_clr,'LineWidth',2);
    % refs / axis
    line([tm(start), tm(start)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    line([tm(stop), tm(stop)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    plot(tm,ref_plot_mean_cat1,'Color',[0 0 0 0.1],'LineWidth',3)
    xlabel('time (sec)'); ylabel('distance to mean (mm)'); axis square
    title('SP: resting paw');
    axis([tm(1) tm(end) ymin ymax])
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);

    subplot(232)
    plot(tm,cat2_fail_miss_distMean_sp,'Color',cat(2,miss_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat2_fail_miss_distMean_sp,2),'Color',miss_clr,'LineWidth',2);
    plot(tm,cat2_suc_distMean_sp,'Color',cat(2,suc_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat2_suc_distMean_sp,2),'Color',suc_clr,'LineWidth',2);
    % refs / axis
    line([tm(start), tm(start)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    line([tm(stop), tm(stop)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    plot(tm,ref_plot_mean_cat2,'Color',[0 0 0 0.1],'LineWidth',3)
    xlabel('time (sec)'); ylabel('distance to mean (mm)'); axis square
    title({'Space-only distance';'SP: lifted paw'});
    axis([tm(1) tm(end) ymin ymax])
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);

    subplot(233)
    plot(tm,cat3_all_distMean_sp,'Color',cat(2,suc_clr,transp),'LineWidth',0.5); hold on
    plot(tm,nanmean(cat3_all_distMean_sp,2),'Color',suc_clr,'LineWidth',2);
    % refs / axis
    line([tm(start), tm(start)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    line([tm(stop), tm(stop)], [0,15], 'Color', [0 0 0 0.3],'LineWidth',2);
    plot(tm,ref_plot_mean_cat3,'Color',[0 0 0 0.1],'LineWidth',3)
    xlabel('time (sec)'); ylabel('distance to mean (mm)'); axis square
    title('grooming');
    axis([tm(1) tm(end) ymin ymax])
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);

    subplot(234)
    stairs(binsMean_cat1_fail_miss_sp,cdfMean_cat1_fail_miss_sp_norm,'Color',miss_clr,'LineWidth',2); hold on
    stairs(binsMean_cat1_suc_sp,cdfMean_cat1_suc_sp_norm,'Color',suc_clr,'LineWidth',2);
    xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
    title(sprintf('%s%.4f','p-value=',cat1_pMean_sp))
    axis([cdf_min cdf_max 0 1])

    subplot(235)
    stairs(binsMean_cat2_fail_miss_sp,cdfMean_cat2_fail_miss_sp_norm,'Color',miss_clr,'LineWidth',2); hold on
    stairs(binsMean_cat2_suc_sp,cdfMean_cat2_suc_sp_norm,'Color',suc_clr,'LineWidth',2);
    xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
    title(sprintf('%s%.4f','p-value=',cat2_pMean_sp))
    axis([cdf_min cdf_max 0 1])

    subplot(236)
    stairs(binsMean_cat3_all_sp,cdfMean_cat3_all_sp_norm,'Color',suc_clr,'LineWidth',2);
    xlabel('distance integral'); ylabel ('cdf'); grid on; axis square
    set(gca,'linewidth',1.5,'nextplot','add','box','off','layer','top','GridAlpha',0.05);
    set(gcf,'Position',[270 262 882 535])
    axis([cdf_min cdf_max 0 1])

    if save_flag, saveas(gcf,strcat(new_folder_sucMiss,filesep,'distMean_sp.png'),'png'); end
    % saveas(gcf,strcat(new_folder_out,filesep,'distMean_sp'),'epsc');

