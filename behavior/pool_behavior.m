%% pool behavior from logs
clear; close all; clc;


%% paths
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
group = '20230511_ChocolateGroup';
setup = 'headfixed_dynamicTarget';
mat_dir = fullfile(rootdir,"behavior_data","analyzed_data","mat_files",group,setup);

% Save output
save_name = 'behav_group';
save_out = fullfile(rootdir,"behavior_data","analyzed_data","output_files",group,setup,save_name);
if ~exist(save_out,"dir"), mkdir(save_out); end

%% Load

% Find mouse folders
folders_mice = dir(mat_dir);
folders_mice=folders_mice(3:end,:);
num_animals=size(folders_mice,1);
% Allocate space for mouse cells
mice = cell(num_animals,1);
mice_path = cell(num_animals,1);

% Allocate pace for variables
sess = ['R1';'R2';'R3';'R4';'R5';'R6';'R7'];
n_sess = size(sess,1);

blocks12 = nan(5,2,n_sess,num_animals);
blocks12_switch = nan(2,n_sess,num_animals);
block1_push = nan(n_sess,num_animals);

med_time_to_pushPull = nan(2,n_sess,num_animals);
med_time_to_LCR = nan(3,n_sess,num_animals);
med_time_to_DomCnonDom = nan(3,n_sess,num_animals);
med_reachTime_pushPull = nan(2,n_sess,num_animals);

for m = 1: num_animals
    mice(m,1)=cellstr(folders_mice(m,1).name);
    mice_path(m,1) = cellstr(strcat(mat_dir,filesep, char(mice(m,1))));

    for s=1:n_sess

        mat_file1 = char(strcat(mice_path(m,1),filesep,sess(s,:),filesep,'behavior_session.mat'));

        if exist(mat_file1,"file")
            fprintf('%s%s%s%s%s','loading mouse ',char(mice(m,1)),', session ',sess(s,:))
            fprintf('\n')

            load(mat_file1,'behavior')

            blocks12(:,1,s,m) = behavior.overview.blocks12_startStop_val;
            blocks12(:,2,s,m) = behavior.overview.blocks12_startStop_inval;
            blocks12_switch(:,1,s,m) = behavior.overview.block12_switch_val;
            blocks12_switch(:,2,s,m) = behavior.overview.block12_switch_inval;

            push_start = behavior.init.push_bounds(1);
            block1_push(s,m) = push_start==1;

            med_time_to_pushPull(1,s,m) = median(behavior.overview.time_to_push,"omitnan");
            med_time_to_pushPull(2,s,m) = median(behavior.overview.time_to_pull,"omitnan");
            
            med_time_to_LCR(1,s,m) = median(behavior.overview.time_to_left,"omitnan");
            med_time_to_LCR(2,s,m) = median(behavior.overview.time_to_center,"omitnan");
            med_time_to_LCR(3,s,m) = median(behavior.overview.time_to_right,"omitnan");

            med_time_to_DomCnonDom(1,s,m) = median(behavior.overview.time_dom_paw,"omitnan");
            med_time_to_DomCnonDom(2,s,m) = median(behavior.overview.time_to_center,"omitnan");
            med_time_to_DomCnonDom(3,s,m) = median(behavior.overview.time_nonDom,"omitnan");

            thr_time_reach = 5;
            r_time_to_push_tmp = behavior.reach.timeto.time_to_reach(...
                intersect(find(behavior.reach.flag_reach_wrogly_detected==0),behavior.init.idx_trial_push));
            r_time_to_pull_tmp = behavior.reach.timeto.time_to_reach(...
                intersect(find(behavior.reach.flag_reach_wrogly_detected==0),behavior.init.idx_trial_pull));
            r_time_to_push = r_time_to_push_tmp(r_time_to_push_tmp<thr_time_reach);
            r_time_to_pull = r_time_to_pull_tmp(r_time_to_pull_tmp<thr_time_reach);

            med_reachTime_pushPull(1,s,m) = median(r_time_to_push,"omitnan");
            med_reachTime_pushPull(2,s,m) = median(r_time_to_pull,"omitnan");

        end
    end
end
fprintf(['%s' ...
    '\n'],'loaded!')

%% Plot
m=5
s=2
clr_valid = behavior.colors.clr_valid;
clr_invalid = behavior.colors.clr_invalid;
% Figure
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};

figure()
a=plot(squeeze(blocks12(:,1,s,m)),'-o','LineWidth',2,'Color',clr_valid);
hold on
b=plot(squeeze(blocks12(:,2,s,m)),'-o','LineWidth',2,'Color',clr_invalid);
plot([2,4], squeeze(blocks12_switch(:,1,s,m)),'--','LineWidth',2,'Color',[.7 .7 .7]);
plot([2,4], squeeze(blocks12_switch(:,2,s,m)),'--','LineWidth',2,'Color',[.7 .7 .7]);
xline(3,'-','switch','Color',[.2 .2 .2],...
    'linewidth',3,'fontsize',10,'LabelHorizontalAlignment','center');
hold off
%legend([a,b],'valid','invalid','box','off','location','northeast')
xlim([.5 5.5])
set(gca,axeOpt{:})
xticks([1,2,4,5]);
xticklabels(['block1 start';'block1 end  ';'block2 start';'block2 end  '])
ylabel('counts');
set(gcf,'Position',[2409         389         367         367],'Color','w')
saveas(gcf,fullfile(save_out,['blocks_',char(mice(m)),'_',sess(s,:)]),'png');
print(gcf, fullfile(save_out, ['blocks_',char(mice(m)),'_',sess(s,:)]), '-dpdf', '-painters');


%%

valid_ratio = squeeze(blocks12(:,1,:,:)./sum(blocks12,2));
valid_ratio_switch = valid_ratio(2:2:4,:,:);
meanS_valid = squeeze(mean(valid_ratio,2,"omitnan"));
meanS_valid_switch =  meanS_valid(2:2:4,:);
transpa = .2;

figure()
ff = tiledlayout(2,3);
s=1:7;
for m = 1:num_animals
    nexttile
plot(squeeze(valid_ratio(:,s,m)),'-','LineWidth',2,'Color',cat(2,clr_valid,transpa)); hold on
plot([2,4], squeeze(valid_ratio_switch(:,s,m)),'--','LineWidth',2,'Color',[.7 .7 .7 transpa]);
xline(3,'-','switch','Color',[.2 .2 .2],...
    'linewidth',3,'fontsize',10,'LabelHorizontalAlignment','center');
yline(.5,'-.','LineWidth',2,'Color',[.8 .8 .8])
plot(squeeze(meanS_valid(:,m)),'o-','LineWidth',2,'Color',clr_valid,'MarkerFaceColor',clr_valid); hold on
plot([2,4], squeeze(meanS_valid_switch(:,m)),'--','LineWidth',2,'Color',[.7 .7 .7]);
title(char(mice(m)),'Interpreter','none')
xlim([.5 5.5])
set(gca,axeOpt{:})
xticks([1,2,4,5]);
xticklabels(['block1 start';'block1 end  ';'block2 start';'block2 end  '])
ylabel('valid ratio')
end
set(gcf,'Position',[2042         165        1267         791],'color','w');
saveas(gcf,fullfile(save_out,'valid_ratio_indivMouse_all_sess'),'png');
print(gcf, fullfile(save_out, 'valid_ratio_indivMouse_all_sess.pdf'), '-dpdf', '-painters');


%%
meanM_valid = squeeze(mean(meanS_valid,2,"omitnan"));
m_learn = [1,2,4,5];

figure()
plot(squeeze(meanS_valid(:,m_learn)),'-','LineWidth',2,'color',cat(2,clr_valid,transpa)); hold on
plot(squeeze(meanS_valid(:,m_learn)),'o','LineWidth',1,'MarkerEdgeColor',clr_valid); hold on
plot([2,4], squeeze(meanS_valid_switch(:,m_learn)),'--','LineWidth',2,'Color',[.7 .7 .7 transpa]);
plot(squeeze(meanM_valid),'o-','LineWidth',3,'Color',clr_valid,'MarkerFaceColor',clr_valid); hold on
plot(squeeze(meanM_valid),'o','LineWidth',5,'Color',clr_valid,'MarkerFaceColor',clr_valid); hold on

xline(3,'-','switch','Color',[.2 .2 .2],...
    'linewidth',3,'fontsize',10,'LabelHorizontalAlignment','center');
yline(.5,'-.','LineWidth',2,'Color',[.8 .8 .8])
xlim([.5 5.5])
set(gca,axeOpt{:})
xticks([1,2,4,5]);
xticklabels(['block1 start';'block1 end  ';'block2 start';'block2 end  '])
xtickangle(45)
ylabel('valid ratio')
set(gcf,'Position',[2359         307         603         555],'color','w');
saveas(gcf,fullfile(save_out,'valid_ratio_allMice'),'png');
print(gcf, fullfile(save_out, 'valid_ratio_allMice.pdf'), '-dpdf', '-painters');


%% Mean push/pull
clr_push = behavior.colors.push_clr;
clr_pull = behavior.colors.pull_clr;
figure()
ff = tiledlayout(2,3);
s=1:7;
for m = 1:num_animals
    nexttile
x_jitter = 1 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_time_to_pushPull(1,:,m),20,'filled','MarkerFaceColor',clr_push,'MarkerFaceAlpha',0.3); hold on
scatter(1,mean(med_time_to_pushPull(1,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_push); 

x_jitter = 2 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_time_to_pushPull(2,:,m),20,'filled','MarkerFaceColor',clr_pull,'MarkerFaceAlpha',0.3); hold on
scatter(2,mean(med_time_to_pushPull(2,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_pull); 
axis([.5 2.5 0 20])
set(gca,axeOpt{:})
xticks([1,2]);
xticklabels(['push';'pull'])
ylabel('time to init trial (s)')
title(char(mice(m)),'Interpreter','none')

end
set(gcf,'Position',[2239         148        1200         771],'color','w');
saveas(gcf,fullfile(save_out,'time_to_push_pull'),'png');
%print(gcf, fullfile(save_out, 'valid_ratio_indivMouse_all_sess.pdf'), '-dpdf', '-painters');

%% Mean push/pull - REACH -- INTERFERENCE

figure()
ff = tiledlayout(2,3);
s=1:7;
for m = 1:num_animals
    nexttile
x_jitter = 1 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_reachTime_pushPull(1,:,m),20,'filled','MarkerFaceColor',clr_push,'MarkerFaceAlpha',0.3); hold on
scatter(1,mean(med_reachTime_pushPull(1,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_push); 

x_jitter = 2 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_reachTime_pushPull(2,:,m),20,'filled','MarkerFaceColor',clr_pull,'MarkerFaceAlpha',0.3); hold on
scatter(2,mean(med_reachTime_pushPull(2,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_pull); 
axis([.5 2.5 0 4])
set(gca,axeOpt{:})
xticks([1,2]);
xticklabels(['push';'pull'])
ylabel('time to collect water (s)')
title(char(mice(m)),'Interpreter','none')

end
set(gcf,'Position',[2239         148        1200         771],'color','w');
saveas(gcf,fullfile(save_out,'time_to_REACH_push_pull'),'png');
%print(gcf, fullfile(save_out, 'valid_ratio_indivMouse_all_sess.pdf'), '-dpdf', '-painters');

%% Mean push/pull


mean_init_pp = squeeze(mean(med_time_to_pushPull(:,:,:),2,'omitnan'));
mean_reach_pp = squeeze(mean(med_reachTime_pushPull(:,:,:),2,'omitnan'));

figure()
subplot(121)
boxplot(mean_init_pp','Colors',[.9 .9 .9],'Widths',.3); hold on
x_jitter = 1 + 0.1 * (rand(num_animals, 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,mean_init_pp(1,:),20,'filled','MarkerFaceColor',clr_push,'MarkerFaceAlpha',0.3); hold on
scatter(1,mean(mean_init_pp(1,:),2,'omitnan'),60,'filled','MarkerFaceColor',clr_push); 

x_jitter = 2 + 0.1 * (rand(num_animals, 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,mean_init_pp(2,:),20,'filled','MarkerFaceColor',clr_pull,'MarkerFaceAlpha',0.3); hold on
scatter(2,mean(mean_init_pp(2,:),2,'omitnan'),60,'filled','MarkerFaceColor',clr_pull); 
axis([.5 2.5 3 11])
set(gca,axeOpt{:})
xticks([1,2]);
xticklabels(['push';'pull'])
ylabel('time to init trial (s)')

subplot(122)
boxplot(mean_reach_pp','Colors',[.9 .9 .9],'Widths',.3); hold on
x_jitter = 1 + 0.1 * (rand(num_animals, 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,mean_reach_pp(1,:),20,'filled','MarkerFaceColor',clr_push,'MarkerFaceAlpha',0.3); hold on
scatter(1,mean(mean_reach_pp(1,:),2,'omitnan'),60,'filled','MarkerFaceColor',clr_push); 

x_jitter = 2 + 0.1 * (rand(num_animals, 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,mean_reach_pp(2,:),20,'filled','MarkerFaceColor',clr_pull,'MarkerFaceAlpha',0.3); hold on
scatter(2,mean(mean_reach_pp(2,:),2,'omitnan'),60,'filled','MarkerFaceColor',clr_pull); 
axis([.5 2.5 0.7 2.5])
set(gca,axeOpt{:})
xticks([1,2]);
xticklabels(['push';'pull'])
ylabel('time to collect water (s)')

saveas(gcf,fullfile(save_out,'time_push_pull_POOLED'),'png');
print(gcf, fullfile(save_out, 'time_push_pull_POOLED.pdf'), '-dpdf', '-painters');


%% Mean LRC
clr_left = behavior.colors.left_color;
clr_center = behavior.colors.center_color;
clr_right = behavior.colors.right_color;

figure()
ff = tiledlayout(2,3);
s=1:7;
for m = 1:num_animals
    nexttile
x_jitter = 1 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_time_to_LCR(1,:,m),20,'filled','MarkerFaceColor',clr_left,'MarkerFaceAlpha',0.3); hold on
scatter(1,mean(med_time_to_LCR(1,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_left); 

x_jitter = 2 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_time_to_LCR(2,:,m),20,'filled','MarkerFaceColor',clr_center,'MarkerFaceAlpha',0.3); hold on
scatter(2,mean(med_time_to_LCR(2,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_center); 

x_jitter = 3 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_time_to_LCR(3,:,m),20,'filled','MarkerFaceColor',clr_right,'MarkerFaceAlpha',0.3); hold on
scatter(3,mean(med_time_to_LCR(3,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_right); 

axis([.5 3.5 0 5])
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['left  ';'center';'right '])
ylabel('time to collect water (s)')
title(char(mice(m)),'Interpreter','none')

end
set(gcf,'Position',[2239         148        1200         771],'color','w');
saveas(gcf,fullfile(save_out,'time_to_reach_LRC'),'png');
%print(gcf, fullfile(save_out, 'valid_ratio_indivMouse_all_sess.pdf'), '-dpdf', '-painters');


%% Mean Dominant side
clr_dom = behavior.colors.clr_dom;
clr_center = behavior.colors.center_color;
clr_nondon = behavior.colors.clr_nondom;

figure()
ff = tiledlayout(2,3);
s=1:7;
for m = 1:num_animals
    nexttile
x_jitter = 1 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_time_to_DomCnonDom(1,:,m),20,'filled','MarkerFaceColor',clr_dom,'MarkerFaceAlpha',0.3); hold on
scatter(1,mean(med_time_to_DomCnonDom(1,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_dom); 

x_jitter = 2 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_time_to_DomCnonDom(2,:,m),20,'filled','MarkerFaceColor',clr_center,'MarkerFaceAlpha',0.3); hold on
scatter(2,mean(med_time_to_DomCnonDom(2,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_center); 

x_jitter = 3 + 0.2 * (rand(length(s), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,med_time_to_DomCnonDom(3,:,m),20,'filled','MarkerFaceColor',clr_nondon,'MarkerFaceAlpha',0.3); hold on
scatter(3,mean(med_time_to_DomCnonDom(3,:,m),2,'omitnan'),60,'filled','MarkerFaceColor',clr_nondon); 

axis([.5 3.5 0 5])
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['dominant paw side    ';'center               ';'non-dominant paw side'])
ylabel('time to collect water (s)')
title(char(mice(m)),'Interpreter','none')

end
set(gcf,'Position',[2239         148        1200         771],'color','w');
saveas(gcf,fullfile(save_out,'time_to_reach_DomNonDom'),'png');
%print(gcf, fullfile(save_out, 'valid_ratio_indivMouse_all_sess.pdf'), '-dpdf', '-painters');

%% 
mean_sess_timeDom = squeeze(mean(med_time_to_DomCnonDom,2,"omitnan"));
[h_tDom, sig_tDom]=ttest2(mean_sess_timeDom(1,:),mean_sess_timeDom(3,:));

figure()
boxplot(mean_sess_timeDom(:,:)','Colors',[.9 .9 .9],'Widths',.3)
hold on
x_jitter = 1 + 0.2 * (rand(length(num_animals), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,mean_sess_timeDom(1,:),20,'filled','MarkerFaceColor',clr_dom,'MarkerFaceAlpha',0.5); hold on
scatter(1,mean(mean_sess_timeDom(1,:),2,'omitnan'),60,'filled','MarkerFaceColor',clr_dom); 

x_jitter = 2 + 0.2 * (rand(length(num_animals), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,mean_sess_timeDom(2,:),20,'filled','MarkerFaceColor',clr_center,'MarkerFaceAlpha',0.5); hold on
scatter(2,mean(mean_sess_timeDom(2,:),2,'omitnan'),60,'filled','MarkerFaceColor',clr_center); 

x_jitter = 3 + 0.2 * (rand(length(num_animals), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,mean_sess_timeDom(3,:),20,'filled','MarkerFaceColor',clr_nondon,'MarkerFaceAlpha',0.3); hold on
scatter(3,mean(mean_sess_timeDom(3,:),2,'omitnan'),60,'filled','MarkerFaceColor',clr_nondon); 

axis([.5 3.5 .5 2.5])
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['dominant paw side    ';'center               ';'non-dominant paw side'])
ylabel('time to collect water (s)')
set(gcf,'Position',[2272         230         398         607],'color','w');
ttest(mean_sess_timeDom(:,:)')
saveas(gcf,fullfile(save_out,'time_to_reach_DomNonDom_MicePooled'),'png');
print(gcf, fullfile(save_out, 'time_to_reach_DomNonDom_MicePooled.pdf'), '-dpdf', '-painters');


%%
thr_time_reach = 5;
 r_time_to_push_tmp = behavior.reach.timeto.time_to_reach(...
    intersect(find(behavior.reach.flag_reach_wrogly_detected==0),behavior.init.idx_trial_push));
r_time_to_pull_tmp = behavior.reach.timeto.time_to_reach(...
    intersect(find(behavior.reach.flag_reach_wrogly_detected==0),behavior.init.idx_trial_pull));
r_time_to_push = r_time_to_push_tmp(r_time_to_push_tmp<thr_time_reach);
r_time_to_pull = r_time_to_pull_tmp(r_time_to_pull_tmp<thr_time_reach);
