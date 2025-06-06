% Get extra-good neurons, calculate instantaneous FR, order neurons
clear; close all; clc

%% Load data & params
% Group and individual
group = '20230511_ChocolateGroup';
setup = 'headfixed_dynamicTarget';
animals = {...
    'CoteDor',...
    'Lindt',...
    'Toblerone',...
    'Milka',...
    'FerreroRocher'};
animal_idx = 4;
mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});
sess = 'R3';

% All reaches? or only the "good bhv" window considered when spike sorting?
ephys_win_only_flag = 1;
flag_dlc = 1;

if ismember(animal_idx,[1,4,5])
    paw_pref = 'R';
elseif ismember(animal_idx,[2,3])
    paw_pref = 'L';
else
    print('unknown paw pref');
end
behavior.paw_pref = paw_pref;

%% paths
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
behavior_dir = fullfile(rootdir,"behavior_data","raw_data",group,setup,mouse,sess);

% Load mat files
load(fullfile(behavior_dir,"behavior_session.mat"));

% Save output
save_name = 'behav_log_overview';
save_mat = fullfile(rootdir,"behavior_data","analyzed_data","mat_files",group,setup,mouse,sess);
save_out = fullfile(rootdir,"behavior_data","analyzed_data","output_files",group,setup,mouse,sess,save_name);
if ~exist(save_mat,"dir"), mkdir(save_mat); end
if ~exist(save_out,"dir"), mkdir(save_out); end

% Figure
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};

% Colors
clr_left = behavior.colors.left_color;
clr_right = behavior.colors.right_color;
clr_center = behavior.colors.center_color;

clr_push = behavior.colors.push_clr;
clr_pull = behavior.colors.pull_clr;

clr_valid = [63 130 109]./256;
clr_invalid = [111 26 7]./256;

clr_dom = [48 131 220]./256;
clr_nondom = [158 216 219]./256;

% Save var
behavior.colors.clr_valid = clr_valid;
behavior.colors.clr_invalid = clr_invalid;
behavior.colors.clr_dom = clr_dom;
behavior.colors.clr_nondom = clr_nondom;


%% Count valid/invalid early and late in block
behaving_idx = behavior.counts.units_slide_min<behavior.behavior_duration.time_end/60;
couts_push = behavior.counts.count_valPush_per_unit(behaving_idx);
couts_pull = behavior.counts.count_valPull_per_unit(behaving_idx);
couts_ipush = behavior.counts.count_invalPush_per_unit(behaving_idx);
couts_ipull = behavior.counts.count_invalPull_per_unit(behaving_idx);
sliding_win = 150; % in sec
sliding_step = 20;
win_delay_idx = ceil(sliding_win/sliding_step);
if behavior.init.transition_2push(1)==1
    % Fist block
    block1_val_start = couts_push(1);
    block1_inval_start = couts_ipull(1);
    last_idx_block1 = find(isnan(couts_ipull),1,'first')-1;
    block1_val_end = couts_push(last_idx_block1);
    block1_inval_end = couts_ipull(last_idx_block1);
    % Second bock
    first_idx_block2 = last_idx_block1+1+win_delay_idx;
    block2_val_start = couts_pull(first_idx_block2);
    block2_inval_start = couts_ipush(first_idx_block2);
    last_idx_block2 = find(isnan(couts_ipush(first_idx_block2:end)),1,'first')+...
        first_idx_block2-2;
    if isempty(last_idx_block2)
        block2_val_end = couts_pull(end);
        block2_inval_end = couts_ipush(end);
    else
        block2_val_end = couts_pull(last_idx_block2);
        block2_inval_end = couts_ipush(last_idx_block2);
    end

elseif behavior.init.transition_2pull(1)==1
    % Fist block
    block1_val_start =couts_pull(1);
    block1_inval_start = couts_ipush(1);
    last_idx_block1 = find(isnan(couts_ipush),1,'first')-1;
    block1_val_end = couts_pull(last_idx_block1);
    block1_inval_end = couts_ipush(last_idx_block1);
    % Second bock
    first_idx_block2 = last_idx_block1+1+win_delay_idx;
    block2_val_start = couts_push(first_idx_block2);
    block2_inval_start =couts_ipull(first_idx_block2);
    last_idx_block2 = find(isnan(couts_ipull(first_idx_block2:end)),1,'first')+...
        first_idx_block2-2;
    if isempty(last_idx_block2)
        block2_val_end = couts_push(end);
        block2_inval_end = couts_ipull(end);
    else
        block2_val_end = couts_push(last_idx_block2);
        block2_inval_end = couts_ipull(last_idx_block2);
    end

else
    disp('unknown first block!')
end

blocks12_startStop_val = [block1_val_start,block1_val_end,nan,block2_val_start,block2_val_end];
blocks12_startStop_inval = [block1_inval_start,block1_inval_end,nan,block2_inval_start,block2_inval_end];
block12_switch_val = [block1_val_end,block2_inval_start];
block12_switch_inval = [block1_inval_end,block2_val_start];

% Save var
behavior.overview.blocks12_startStop_val = blocks12_startStop_val;
behavior.overview.blocks12_startStop_inval = blocks12_startStop_inval;
behavior.overview.block12_switch_val = block12_switch_val;
behavior.overview.block12_switch_inval = block12_switch_inval;
% Figure
figure()
a=plot(blocks12_startStop_val,'-o','LineWidth',2,'Color',clr_valid);
hold on
b=plot(blocks12_startStop_inval,'-o','LineWidth',2,'Color',clr_invalid);
plot([2,4],block12_switch_val,'--','LineWidth',2,'Color',[.7 .7 .7]);
plot([2,4],block12_switch_inval,'--','LineWidth',2,'Color',[.7 .7 .7]);
xline(3,'-','switch','Color',[.2 .2 .2],...
    'linewidth',3,'fontsize',10,'LabelHorizontalAlignment','center');
hold off
legend([a,b],'valid','invalid','box','off','location','northeast')
xlim([.5 5.5])
set(gca,axeOpt{:})
xticks([1,2,4,5]);
xticklabels(['block1 start';'block1 end  ';'block2 start';'block2 end  '])
ylabel('counts')
 saveas(gcf,fullfile(save_out,'val_inval_block_start_end'),'png');

%% Time to push/pull
thr_time = 30;
time_to_push = behavior.init.timeto.time_to_init_push(behavior.init.timeto.time_to_init_push<thr_time);
time_to_pull = behavior.init.timeto.time_to_init_pull(behavior.init.timeto.time_to_init_pull<thr_time);

% Save var
behavior.overview.time_to_push = time_to_push;
behavior.overview.time_to_pull = time_to_pull;

figure()
x_jitter = 1 + 0.2 * (rand(length(time_to_push), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,time_to_push,20,'filled','MarkerFaceColor',clr_push,'MarkerFaceAlpha',0.3); hold on
scatter(1,median(time_to_push),60,'filled','MarkerFaceColor',clr_push); 

x_jitter = 2 + 0.2 * (rand(length(time_to_pull), 1) - 0.5);  % Add jitter around x = 2
scatter(x_jitter,time_to_pull,20,'filled','MarkerFaceColor',clr_pull,'MarkerFaceAlpha',0.3); hold on
scatter(2,median(time_to_pull),60,'filled','MarkerFaceColor',clr_pull); hold off
axis([.5 2.5 0 20])
set(gca,axeOpt{:})
xticks([1,2]);
xticklabels(['push';'pull'])
ylabel('time to init trial (s)')

saveas(gcf,fullfile(save_out,'time_to_init'),'png');



%% Time to reach left/center/right
%if flag_dlc == true
thr_time_reach = 5;
time_to_left_tmp = behavior.reach.timeto.time_to_reach(...
    intersect(find(behavior.reach.flag_reach_wrogly_detected==0),behavior.reach.left_idx));
time_to_center_tmp = behavior.reach.timeto.time_to_reach(...
    intersect(find(behavior.reach.flag_reach_wrogly_detected==0),behavior.reach.center_idx));
time_to_right_tmp = behavior.reach.timeto.time_to_reach(...
    intersect(find(behavior.reach.flag_reach_wrogly_detected==0),behavior.reach.right_idx));
time_to_left = time_to_left_tmp(time_to_left_tmp<thr_time_reach);
time_to_center = time_to_center_tmp(time_to_center_tmp<thr_time_reach);
time_to_right = time_to_right_tmp(time_to_right_tmp<thr_time_reach);

if strcmp(paw_pref,'R')
    time_dom_paw = time_to_right;
    time_nonDom = time_to_left;
elseif strcmp(paw_pref,'L')
    time_dom_paw = time_to_left;
    time_nonDom = time_to_right;
end

% Save var
behavior.overview.time_to_left = time_to_left;
behavior.overview.time_to_center = time_to_center;
behavior.overview.time_to_right = time_to_right;
behavior.overview.time_dom_paw = time_dom_paw;
behavior.overview.time_nonDom = time_nonDom;


figure()
x_jitter = 1 + 0.2 * (rand(length(time_to_left), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,time_to_left,20,'filled','MarkerFaceColor',clr_left,'MarkerFaceAlpha',0.3); hold on
scatter(1,median(time_to_left),60,'filled','MarkerFaceColor',clr_left); 

x_jitter = 2 + 0.2 * (rand(length(time_to_center), 1) - 0.5);  % Add jitter around x = 2
scatter(x_jitter,time_to_center,20,'filled','MarkerFaceColor',clr_center,'MarkerFaceAlpha',0.3); 
scatter(2,median(time_to_center),60,'filled','MarkerFaceColor',clr_center); 

x_jitter = 3 + 0.2 * (rand(length(time_to_right), 1) - 0.5);  % Add jitter around x = 3
scatter(x_jitter,time_to_right,20,'filled','MarkerFaceColor',clr_right,'MarkerFaceAlpha',0.3); 
scatter(3,median(time_to_right),60,'filled','MarkerFaceColor',clr_right); 
xlim([.5 3.5])
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['left  ';'center';'right '])
ylabel('time to collect water (s)')
saveas(gcf,fullfile(save_out,'time_to_reach_LCR'),'png');


figure()
x_jitter = 1 + 0.2 * (rand(length(time_dom_paw), 1) - 0.5);  % Add jitter around x = 1
scatter(x_jitter,time_dom_paw,20,'filled','MarkerFaceColor',clr_dom,'MarkerFaceAlpha',0.3); hold on
scatter(1,median(time_dom_paw),60,'filled','MarkerFaceColor',clr_dom); 

x_jitter = 2 + 0.2 * (rand(length(time_to_center), 1) - 0.5);  % Add jitter around x = 2
scatter(x_jitter,time_to_center,20,'filled','MarkerFaceColor',clr_center,'MarkerFaceAlpha',0.3); 
scatter(2,median(time_to_center),60,'filled','MarkerFaceColor',clr_center); 

x_jitter = 3 + 0.2 * (rand(length(time_nonDom), 1) - 0.5);  % Add jitter around x = 3
scatter(x_jitter,time_nonDom,20,'filled','MarkerFaceColor',clr_nondom,'MarkerFaceAlpha',0.3); 
scatter(3,median(time_nonDom),60,'filled','MarkerFaceColor',clr_nondom); 
xlim([.5 3.5])
set(gca,axeOpt{:})
xticks([1,2,3]);
xticklabels(['dominant paw side    ';'center               ';'non-dominant paw side'])
ylabel('time to collect water (s)')
saveas(gcf,fullfile(save_out,'time_to_reach_domNonDon'),'png');


%% Save
% Save bhv mat file 
save(fullfile(save_mat,"behavior_session.mat"),'behavior');



