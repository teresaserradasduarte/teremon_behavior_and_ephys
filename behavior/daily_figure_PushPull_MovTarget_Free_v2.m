%% Daily figure - from log files
clear; close all; clc

%% Load data

%person = 'teresa';
person = 'simon';
if strcmp(person,'teresa')
    raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
elseif strcmp(person,'simon')
    raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
end
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
% R:-1 |C:0 | L:1

        
% SESSION
session = 'R4';
flag_no_sync = 0;

% DISPLAY LENGTH
% Sherten session to display (discart unengaged trials)
shorten_display = 1;
%last_trial_disp_time = 38; % minutes
last_trial_disp = 116; % trial idx

%% Path
%session = char(PP_sess(s));
rootdir = strcat(raw_folder,filesep,group,filesep,setup,filesep,mouse,filesep,session);
%rootdir = pwd;

% Search string
searchstr_log = '*GlobalLogInt';
searchstr_load_cells = '*paws_loadCells';
searchstr_cam = '*Camera1Timestamp';
searchstr_water_pos = '*WaterPosition_steps';
searchstr_syncfix = '*syncfixNO'; % remove the NO if what to use


% Files
files_load = wildcardsearch(rootdir, searchstr_load_cells,true);
files_log = wildcardsearch(rootdir, searchstr_log,true);
files_cam = wildcardsearch(rootdir, searchstr_cam,true);
files_water = wildcardsearch(rootdir, searchstr_water_pos,true);
files_syncfix = wildcardsearch(rootdir, searchstr_syncfix,true);

if ~isempty(files_syncfix)
    load(char(files_syncfix))
    flag_realigned_rwd = 1;
else
    flag_realigned_rwd = 0;
end


%% READ LOG
% read log file and write on it
read_log = [];
read_log_all = [];
read_paws_all = [];
read_cam = [];
read_water = [];
ind_last_trial_available_tmp = [];
ind_first_trial_available_tmp = [];
n_logs = size(files_load,1);

for i=1:n_logs
    % read water position -------------------
    read_water_temp = csvread(char(files_water(i,1)));   
    read_water = cat(1,read_water,read_water_temp);

    % read cam ------------------------------
    read_cam_temp = dlmread(char(files_cam(i,1)),' ');
    read_cam = cat(1,read_cam,read_cam_temp);

    % read load cells ----------------------
    temp_load = csvread(char(files_load(i,1)));
    read_paws_all = cat(1,read_paws_all,temp_load);
    paws_play = read_paws_all(:,2);
    paws_time = read_paws_all(:,1);

    % read log ---------------------------------
    read_log_temp = csvread(char(files_log(i,1)));
    read_log_all = cat(1,read_log_all,read_log_temp);

    % find trials
    push_trial_ind = find(read_log_temp(:,2)==11);
    pull_trial_ind = find(read_log_temp(:,2)==10);
    trial_available_tmp = sort(cat(1,push_trial_ind,pull_trial_ind));
    
    % find trials start and end
    first_trial_available = trial_available_tmp(1);
    last_trial_available = trial_available_tmp(end);

    ind_first_trial_available_tmp = cat(1,ind_first_trial_available_tmp,...
        first_trial_available+(i-1));
    ind_last_trial_available_tmp = cat(1,ind_last_trial_available_tmp,...
        last_trial_available+(i-1));
        % read log only full trials
    read_log = cat(1,read_log,...
        read_log_temp(1:last_trial_available,:),...
        zeros(1,2)); 
end
ind_last_trial_available = cumsum(ind_last_trial_available_tmp);
if size(files_load,1)>1
    ind_first_trial_available = ind_first_trial_available_tmp+...
        cat(1,0,ind_last_trial_available(1:end-1));
else
    ind_first_trial_available = ind_first_trial_available_tmp;
end

% logs
logs = read_log(:,2);
sess_start_ind = find(logs==1001); sess_start_ind = sess_start_ind(1);
timelog = read_log(:,1)-read_log(sess_start_ind,1);
session_duration = read_log(end-1,1)-read_log(sess_start_ind,1); % in sec



%% Check sync'ing
syncCam = read_cam(:,1);

% exclude jumps (bug)
to_exclude_down = find(diff(syncCam)<-2.148E9);
to_exclude_up = find(diff(syncCam)>2.148E9);
syncCam(to_exclude_down)=syncCam(to_exclude_down-1);
syncCam(to_exclude_up)=syncCam(to_exclude_up+1);

ind_newTrial_event = find((diff([0;syncCam])~=0));
ind_newTrial_event_up = find((diff([0;syncCam])>1));
ind_newTrial_event_down = find((diff([0;syncCam])<-1));
time_newTrial_cam = read_cam(ind_newTrial_event,2);

new_trial_ind=[1;find(read_log_all(:,2)==2)];
new_reach_ind=[1;find(read_log_all(:,2)==31)];
time_newTrial_log = read_log_all(new_trial_ind,1);
time_newReach_log = read_log_all(new_reach_ind,1);
if flag_realigned_rwd
    ind_newTrial_event_inferred = find((diff([0;sync_state_hat])~=0));
    idx_reach_syncfix = 1:length(ind_newTrial_event_inferred);
end


figure
plot(syncCam); hold on
plot(ind_newTrial_event,ones(size(ind_newTrial_event)).*max(syncCam),'*');
if flag_realigned_rwd
    plot(sync_state_hat.*10e8 + 10e8)
end
hold off

figure
subplot(211), plot(diff(time_newTrial_log(2:end)));
subplot(212), plot(diff(ind_newTrial_event(2:end)));
if flag_realigned_rwd
    hold on, plot(diff(ind_newTrial_event_inferred),'.');
end

if length(time_newTrial_cam)>=length(time_newTrial_log)+2
    flag_crashed = 1;
    time_newTrial_cam_original = time_newTrial_cam;
    time_newTrial_cam = time_newTrial_cam(1:length(time_newTrial_log));
    ind_newTrial_event_original = ind_newTrial_event;
    ind_newTrial_event = ind_newTrial_event(1:length(time_newTrial_log));
    figure(2)
    subplot(211), plot(diff(time_newTrial_log(2:end)));
    subplot(212), plot(diff(ind_newTrial_event(2:end)));
    if flag_realigned_rwd
        hold on, plot(diff(ind_newTrial_event_inferred),'.');
    end
end


% Time of harp vs time of cam
figure
plot(time_newTrial_log(2:end),time_newTrial_cam(2:end)-time_newTrial_cam(2),'-'); hold on
plot(time_newTrial_log(2:end),time_newTrial_cam(2:end)-time_newTrial_cam(2),'.'); 
if flag_realigned_rwd
plot(time_newTrial_log(1+idx_reach_syncfix),ind_newTrial_event_inferred-ind_newTrial_event_inferred(1),'o'); 
end

hold off
xlabel('trial time on harp'); ylabel('trial time on  cam');




%% LOG INFO 
% some of this was computed above, repeated here to make sure everything
% works smoothly even if the session crushed midway :)

% Find relevant moments of trial
% push/pull & trial available
push_trial_ind = find(read_log_all(:,2)==11);
pull_trial_ind = find(read_log_all(:,2)==10);
trial_available_all_ind = sort(cat(1,push_trial_ind,pull_trial_ind));
% remove incomplete trials
trial_available_to_start_ind = setdiff(trial_available_all_ind,ind_last_trial_available);
end_prev_trial_ind = setdiff(trial_available_all_ind,ind_first_trial_available);

% trial init, reach, iti
trial_init_ind=find(logs==2);
water_delivered_ind=find(logs==3);
reached_really_ind=find(logs==31);
ITIs_start_ind=find(logs==4);
manual_water_ind=find(logs==300);
% invalid/jumped trials
water_cleaned_ind=find(logs==30);
reached_ind = sort(cat(1,reached_really_ind,water_cleaned_ind));
inval_push_ind = find(read_log_all(1:length(timelog)-1,2)==-11);
inval_pull_ind = find(read_log_all(1:length(timelog)-1,2)==-10);

% Nr of trials
nr_trials = length(reached_ind);

trials_vec = 1:nr_trials;
if flag_realigned_rwd==1 && idx_reach_syncfix(end)>trials_vec(end)
    idx_reach_syncfix = trials_vec;
end

nr_ITIs = length(ITIs_start_ind);

% Nr of invalid attempts
nr_inval_push = length(inval_push_ind);
nr_inval_pull = length(inval_pull_ind);

% Count of manual before init
count_manual_before_init = length(find(logs(trial_init_ind-1)==300))+...
    length(find(logs(trial_init_ind-2)==300));


%% Push/Pull identity of trials
% find index of trial of PUSH/PULL
is_push = nan(size(trial_available_to_start_ind));
for ii = 1:length(trial_available_to_start_ind)
    is_push(ii) = sum(push_trial_ind==trial_available_to_start_ind(ii));
end
idx_comp_trial_push = find(is_push==1);
idx_comp_trial_pull = find(is_push==0);

transition_push_pull=cat(1,1,find(diff(is_push)~=0));
transition_2push = transition_push_pull(is_push(transition_push_pull+1)==1);
transition_2pull = transition_push_pull(is_push(transition_push_pull+1)==0);

% Bounds
push_bounds_tmp = sort([idx_comp_trial_push(1); ...
    idx_comp_trial_push(diff(idx_comp_trial_push)>1);...
    idx_comp_trial_push(find(diff(idx_comp_trial_push)>1)+1);...
    idx_comp_trial_push(end)]);
push_bounds =  reshape(push_bounds_tmp,[2 length(push_bounds_tmp)/2]);
pull_bounds_tmp = sort([idx_comp_trial_pull(1); ...
    idx_comp_trial_pull(diff(idx_comp_trial_pull)>1);...
    idx_comp_trial_pull(find(diff(idx_comp_trial_pull)>1)+1);...
    idx_comp_trial_pull(end)]);
pull_bounds =  reshape(pull_bounds_tmp,[2 length(pull_bounds_tmp)/2]);


%% Moving water 
% Find when water was moved and where to
water_zero_pos = find(logs==0); water_zero_pos=water_zero_pos(1);
if water_zero_pos>sess_start_ind+50
    water_zero_pos = sess_start_ind;
end
water_moved_Right_all = find(logs==610); 
water_moved_Left_all = find(logs==611); 
water_moved_all_tmp = sort(cat(1,water_moved_Right_all,water_moved_Left_all));
water_moved_all_ind = water_moved_all_tmp(water_moved_all_tmp > water_zero_pos);
% if water moved many times before start the trial, exclude all but the last
if sum(water_moved_all_ind<ind_first_trial_available(1))>1 
    before_first_trial_moves = water_moved_all_ind(water_moved_all_ind<ind_first_trial_available);
    water_moved_ind = cat(1,before_first_trial_moves(end),water_moved_all_ind(water_moved_all_ind>ind_first_trial_available));
else
    water_moved_ind = water_moved_all_ind;
end
% if the first was at the center, maybe wasn't moved - make sure teh first
% position was recognized
if water_moved_ind(1)>water_zero_pos+50
    water_moved_ind = cat(1,water_zero_pos,water_moved_ind);
end

% Find the position at the transition time (where to)
water_pos_trans_ind = zeros(size(water_moved_ind));
water_pos_trans = zeros(size(water_moved_ind));
water_transition = zeros(size(water_moved_ind));
water_moved_trial_idx = zeros(size(water_moved_ind));

for i = 1: length(water_moved_ind)
    [~,water_pos_trans_ind(i)] = min(abs(read_water(:,1)-read_log(water_moved_ind(i),1)));
    water_pos_trans(i) = read_water(water_pos_trans_ind(i),2);
    if i==1
        water_moved_trial_idx(i)=1;
    else
        close_trials_tmp = water_moved_ind(i)-trial_init_ind;
        [~,water_moved_trial_idx(i)] = min(close_trials_tmp(close_trials_tmp>=0));
        water_moved_trial_idx(i) = water_moved_trial_idx(i)+1;
    end
end
water_transition(water_pos_trans<-30) = -1;
water_transition(water_pos_trans>30) = 1;  

% Find the trials of each location
% Find the trials of each location
if water_moved_trial_idx(end) == nr_trials
    water_transition_final = water_transition(1:end-1);
else
    water_moved_trial_idx = cat(1,water_moved_trial_idx,nr_trials+1);
    water_transition_final = water_transition;
end
%water_transition_final =  water_transition(1:length(water_moved_trial_idx)-2); % to make sure we ignore if the last location have no triall

Right_start_stop=sort(cat(1,water_moved_trial_idx(water_transition_final==-1),water_moved_trial_idx(find(water_transition_final==-1)+1)-1));
Left_start_stop=sort(cat(1,water_moved_trial_idx(water_transition_final==1),water_moved_trial_idx(find(water_transition_final==1)+1)-1));
Center_start_stop=sort(cat(1,water_moved_trial_idx(water_transition_final==0),water_moved_trial_idx(find(water_transition_final==0)+1)-1));

right_idx = [];
for i=1:2:numel(Right_start_stop)
    right_idx_tmp = Right_start_stop(i):Right_start_stop(i+1);
    right_idx = [right_idx, right_idx_tmp];
end
left_idx = [];
for i=1:2:numel(Left_start_stop)
    left_idx_tmp = Left_start_stop(i):Left_start_stop(i+1);
    left_idx = [left_idx, left_idx_tmp];
end
center_idx = [];
for i=1:2:numel(Center_start_stop)
    center_idx_tmp = Center_start_stop(i):Center_start_stop(i+1);
    center_idx = [center_idx, center_idx_tmp];
end

% transition trials
moved_left = Left_start_stop(1:2:end);
moved_right = Right_start_stop(1:2:end);
moved_center = Center_start_stop(1:2:end);

% bounds
if ~isempty(water_transition)
    left_bounds_tmp = sort([left_idx(1); left_idx(diff(left_idx)>1)'; ...
        left_idx(find(diff(left_idx)>1)+1)'; left_idx(end)]);
    left_bounds =  reshape(left_bounds_tmp,[2 length(left_bounds_tmp)/2]);
    right_bounds_tmp =sort([right_idx(1); right_idx(diff(right_idx)>1)';...
        right_idx(find(diff(right_idx)>1)+1)'; right_idx(end)]);
    right_bounds =  reshape(right_bounds_tmp,[2 length(right_bounds_tmp)/2]);
    center_bounds_tmp =sort([center_idx(1); center_idx(diff(center_idx)>1)';...
        center_idx(find(diff(center_idx)>1)+1)'; center_idx(end)]);
    center_bounds =  reshape(center_bounds_tmp,[2 length(center_bounds_tmp)/2]);
end

%% Time of...
% Time of reach/init
trial_available_time = timelog(trial_available_to_start_ind);
init_time = timelog(trial_init_ind);
water_delivery_time = timelog(water_delivered_ind);
reach_time = timelog(reached_ind);

% Time of reach offset corrected
offset = 0.180;
reach_time_corrected = reach_time-offset;
if flag_realigned_rwd
    reach_times_inferred_approxlog_corrected = trial_reach_times_inferred_approxlog - offset;
end


% Time of valid/invalid push pull
idx_comp_trial_push = idx_comp_trial_push<length(init_time);
valPush_time = init_time(idx_comp_trial_push);
valPull_time = init_time(idx_comp_trial_pull);
invalPush_time = timelog(inval_push_ind);
invalPull_time =  timelog(inval_pull_ind);

% time of water moving
if ~isempty(water_transition)
    moved_left_time = timelog(water_moved_ind(water_transition==1));
    moved_right_time = timelog(water_moved_ind(water_transition==-1));
    moved_center_time = timelog(water_moved_ind(water_transition==0));
end

%% Time (duration) to...
% Time to reach/init/iti (durations)
trial_available_time = trial_available_time(1:length(init_time));
time_to_init = init_time-trial_available_time;
time_to_reach = reach_time-water_delivery_time;
time_to_reach_corrected = reach_time_corrected - water_delivery_time;

if flag_realigned_rwd
    time_to_reach_inferred = trial_reach_times_inferred-trial_rwd_times(idx_reach_syncfix);
    time_to_reach_inferred_approxlog = trial_reach_times_inferred_approxlog-trial_rwd_times(idx_reach_syncfix);
    time_to_reach_inferred_approxlog_corrected = reach_times_inferred_approxlog_corrected-trial_rwd_times(idx_reach_syncfix);
end
ITI_duration = timelog(end_prev_trial_ind)-timelog(ITIs_start_ind);
ITI_duration(ITI_duration<-1) = nan;
flag_reach_wrogly_detected = time_to_reach<0.3;

% figure
% subplot(121)
% plot(trials_vec,time_to_reach,'.'); hold on
% plot(trials_vec(idx_reach_syncfix),time_to_reach_inferred,'.');
% plot(trials_vec,flag_reach_wrogly_detected,'color',[.8 .8 .8])
% ylim([0 15])
% subplot(122)
% plot(time_to_reach(idx_reach_syncfix)-time_to_reach_inferred,'.'); hold on
% plot(trials_vec,flag_reach_wrogly_detected,'color',[.8 .8 .8])
% ylim([-6 2])

% Time to init conditioned on identity
time_to_init_push=time_to_init(idx_comp_trial_push);
time_to_init_pull=time_to_init(idx_comp_trial_pull);

% Time to reqch conditioned on location
time_to_reach_right=time_to_reach(right_idx);
time_to_reach_left=time_to_reach(left_idx);
time_to_reach_center=time_to_reach(center_idx);

%% Counts of init trials and push/pulls
% --------------------------------------
% Sliding window of nr of trials per 2.5 min
% Counts per unit time (sliding window through session)
%        sliding_win = 200;
sliding_win = 150; % in sec
sliding_step = 20;
max_ind = floor((session_duration-sliding_win)/sliding_step);

count_reaches_per_unit = [];
count_init_per_unit = [];
units_slide = [];
push_slide_vec = [];
count_Push_per_unit = [];
count_Pull_per_unit = [];

Push_time = sort(cat(1,valPush_time,invalPush_time));
Pull_time = sort(cat(1,valPull_time,invalPull_time));

for i = 1:max_ind
    if i==1
        sliding_vec = 1:sliding_win;
    else
        sliding_vec = (i-1)*sliding_step+1 :(i-1)*sliding_step+sliding_win;
    end
    reaches_current_unit = numel(find(reach_time>sliding_vec(1) & reach_time<sliding_vec(end)));
    inits_current_unit = numel(find(init_time>sliding_vec(1) & init_time<sliding_vec(end)));
    Push_current_unit = numel(find(Push_time>sliding_vec(1) & Push_time<sliding_vec(end)));
    Pull_current_unit = numel(find(Pull_time>sliding_vec(1) & Pull_time<sliding_vec(end)));
    
    count_reaches_per_unit=[count_reaches_per_unit;reaches_current_unit];
    count_init_per_unit=[count_init_per_unit;inits_current_unit];
    units_slide = [units_slide;sliding_vec(end)];

    count_Push_per_unit=[count_Push_per_unit;Push_current_unit];
    count_Pull_per_unit=[count_Pull_per_unit;Pull_current_unit];

end

units_slide_min = units_slide./60;

ispush_slide_vec = zeros(size(units_slide_min));
if (size(push_bounds,2)==1 && push_bounds(1)==push_bounds(2))
    push_bounds = [push_bounds(1)-1; push_bounds(2)];
end
for j=1:size(push_bounds,2)
    ispush_slide_vec(units_slide>init_time(push_bounds(1,j)) & units_slide<init_time(push_bounds(2,j))) = 1;
end

count_valPush_per_unit = count_Push_per_unit;
count_invalPush_per_unit = count_Push_per_unit;
count_valPush_per_unit(ispush_slide_vec==0)=nan;
count_invalPush_per_unit(ispush_slide_vec==1)=nan;

count_valPull_per_unit = count_Pull_per_unit;
count_invalPull_per_unit = count_Pull_per_unit;
count_valPull_per_unit(ispush_slide_vec==1)=nan;
count_invalPull_per_unit(ispush_slide_vec==0)=nan;

% figure
% plot(units_slide_min,count_init_per_unit,'linewidth',6,'Color',[0 0 0 0.2]); hold on
% plot(units_slide_min,count_valPull_per_unit,'-','linewidth',lw,'color',pull_clr);
% plot(units_slide_min,count_valPush_per_unit,'-','linewidth',lw,'color',push_clr); hold on
% plot(units_slide_min,count_invalPush_per_unit,'--','linewidth',lw,'color',[push_clr .8]);
% plot(units_slide_min,count_invalPull_per_unit,'--','linewidth',lw,'color',[pull_clr .8]); hold off

% ---------------------------------------------------------------
%% Load cells forces around trial init
n_points_before = 1000;
n_points_after = 500;
x_around_sec = (-n_points_before:n_points_after)./1000;

% find intervals in paw forces
init_start_stop_ind = zeros(2,nr_trials);
for i = trials_vec
    [~,init_start_stop_ind(1,i)]=min(abs(paws_time-read_log(trial_available_to_start_ind(i),1))); % closest timepoint of start
    [~,init_start_stop_ind(2,i)]=min(abs(paws_time-read_log(trial_init_ind(i),1))); % closest timepoint of stop
end
nr_ind_max = max(init_start_stop_ind(2,:)-init_start_stop_ind(1,:));

% time vector for trials
trial_init = nr_ind_max+1;
x_time = ((1:trial_init+n_points_after)-trial_init)./1000;

% Paw force to initiate  trial
paw_init_play = nan(nr_ind_max+1+n_points_after,nr_trials);
for i = trials_vec
    size_trial = init_start_stop_ind(2,i)-init_start_stop_ind(1,i);
    start_trial_align_end = nr_ind_max-size_trial;
    paw_init_play(start_trial_align_end+1:nr_ind_max+1+n_points_after,i) = paws_play(init_start_stop_ind(1,i):init_start_stop_ind(2,i)+n_points_after);
end

% forces close to hit boundary to initiate trial
paw_init_play_short = paw_init_play(end-(n_points_before+n_points_after):end,:);


% ------------------------------------------
% Load cells forces around invalid PUSH trials

invalPush_loadCell_ind = zeros(nr_inval_push,1);
force_invPush = nan(n_points_before+n_points_after+1,nr_inval_push);
for i = 1:nr_inval_push
    [~,invalPush_loadCell_ind(i)]=min(abs(paws_time-read_log(inval_push_ind(i),1)));
    force_invPush(:,i) = paws_play(invalPush_loadCell_ind(i)-n_points_before:invalPush_loadCell_ind(i)+n_points_after);
end
%toc

% ------------------------------------------
% Load cells forces around invalid PULL trials
invalPull_loadCell_ind = zeros(nr_inval_pull,1);
force_invPull = nan(n_points_before+n_points_after+1,nr_inval_pull);
for i = 1:nr_inval_pull
    [~,invalPull_loadCell_ind(i)]=min(abs(paws_time-read_log(inval_pull_ind(i),1)));
    force_invPull(:,i) = paws_play(invalPull_loadCell_ind(i)-n_points_before:invalPull_loadCell_ind(i)+n_points_after);
end
%toc

% ----------------------------------------------
% Load cells forces around reach
% find index in load cell
%tic
waterReached_loadCell_ind = zeros(nr_trials,1);
force_aroundReach = nan(n_points_before+n_points_after+1,nr_trials);
for i = trials_vec
    [~,waterReached_loadCell_ind(i)]=min(abs(paws_time-read_log(reached_ind(i),1)));
    force_aroundReach(:,i) = paws_play(waterReached_loadCell_ind(i)-n_points_before:waterReached_loadCell_ind(i)+n_points_after);
end
%toc

% -------------------------------------------
% Load cell forces during ITI vs during trial init
% find intervals in paw forces
ITI_start_stop_ind = zeros(2,nr_ITIs);
for i = 1:nr_ITIs
    [~,ITI_start_stop_ind(1,i)]=min(abs(paws_time-read_log(ITIs_start_ind(i),1))); % closest timepoint of start
    [~,ITI_start_stop_ind(2,i)]=min(abs(paws_time-read_log(end_prev_trial_ind(i),1))); % closest timepoint of stop
end
nr_ind_max_iti = max(ITI_start_stop_ind(2,:)-ITI_start_stop_ind(1,:));

% Paw force during ITI
paw_ITI_play = nan(nr_ind_max_iti,nr_ITIs);
for i = 1:nr_ITIs
    size_trial_iti = ITI_start_stop_ind(2,i)-ITI_start_stop_ind(1,i);
    paw_ITI_play(1:size_trial_iti+1,i) = paws_play(ITI_start_stop_ind(1,i):ITI_start_stop_ind(2,i));
end

% General force in load cells ("movement") during ITI and during trial init
sum_gen_force_ITI_later=nansum(abs(paw_ITI_play(1000:4000,:)));
sum_gen_force_ITI_later(sum_gen_force_ITI_later==0)=nan;
sum_gen_force_init_early=nansum(abs(paw_ITI_play(1:3000,:)));



%% DAILY FIGURE - THIS WILL BE REMOVED I THINK
close all
% Figures handles
figOpt = {'color','w','position',[1921          46        1920         963]};
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0.05,'ticklength',[1,1]*.01,'fontsize',12, 'TickDir','out'};

% Colors, lines and sizes: fig properties
push_clr = [132 160 124]./256;
pull_clr = [26 94 99]./256;
%pull2_clr = [21 77 81]./256;
pp_clr = [63 130 109]./256;
right_color = [252 141 89]./256;
center_color = [72 60 50]./256;
left_color = [237 216 146]./256;
rcl_clr = [44 123 182]./256;
cy=[250,218,94]./256;
neutral_init_m = [.2 .2 .2];
neutral_init_t = [.6 .6 .6 0.1];
clrs_blue = [0 0.4470 0.7410];
transp = .2;
ylim_loadcell = [-8000 8000];
spk_size = 2;
event_size = 5;
patch_wd = 10;
%idx_comp_trial_pull = idx_comp_trial_pull(2:end);

% Legth of display
if (shorten_display == 1 && exist('last_trial_disp','var'))
    last_trial_disp_time = init_time(last_trial_disp);
elseif (shorten_display == 1 && exist('last_trial_disp_time','var'))
    last_trial_disp_time = last_trial_disp_time*60;
    tmp_min = init_time-(last_trial_disp_time);
    [~,last_trial_disp] =max(tmp_min(tmp_min<0));
     
else
    last_trial_disp = nr_trials;
    last_trial_disp_time = session_duration;
end

% FIGURE 1: PERFORMANCE ---------------------------------------
fig_performance=figure;

subplot(4,3,[1,4])
% PATCHES
patch_y = [0 30 30 0];
for l=1:size(push_bounds,2)
    patch([push_bounds(1,l) push_bounds(1,l) push_bounds(2,l) push_bounds(2,l)],...
        patch_y,push_clr,'FaceAlpha',0.2,'EdgeColor','none'); hold on
end
for l=1:size(pull_bounds,2)
    patch([pull_bounds(1,l) pull_bounds(1,l) pull_bounds(2,l) pull_bounds(2,l)],...
        patch_y,pull_clr,'FaceAlpha',0.2,'EdgeColor','none'); hold on
end

plot(trials_vec,time_to_init,'color',[.8 .8 .8 transp])
hold on
plot(trials_vec,time_to_init,'k.');
if ~isempty(idx_comp_trial_pull), xline(transition_2pull,'-','pull','Color',pull_clr,'linewidth',1.5,'fontsize',12); end
if ~isempty(idx_comp_trial_push), xline(transition_2push','-','push','Color',push_clr,'linewidth',1.5,'fontsize',12); end
%patch(patch_pull,patch_y,push_clr,'FaceColor',pull1_clr,'EdgeColor','none')
set(gca,axeOpt{:})
ylim([0 30])
xlim([0 last_trial_disp])
%xlabel('trials'); 
ylabel('time to initiate trial (s)')
title(sprintf('%s%s', 'Subject: ',mouse), 'Interpreter', 'none');

subplot(4,3,[7,10])
if ~isempty(water_transition)
    for l=1:size(left_bounds,2)
        patch([left_bounds(1,l) left_bounds(1,l) left_bounds(2,l) left_bounds(2,l)],patch_y,left_color,'FaceAlpha',0.2,'EdgeColor','none'); hold on
    end
    for l=1:size(right_bounds,2)
        patch([right_bounds(1,l) right_bounds(1,l) right_bounds(2,l) right_bounds(2,l)],patch_y,right_color,'FaceAlpha',0.2,'EdgeColor','none'); hold on
    end
    for l=1:size(center_bounds,2)
        patch([center_bounds(1,l) center_bounds(1,l) center_bounds(2,l) center_bounds(2,l)],patch_y,center_color,'FaceAlpha',0.2,'EdgeColor','none'); hold on
    end
end
plot(trials_vec,time_to_reach,'color',[.8 .8 .8 transp])
hold on
if flag_realigned_rwd
    %plot(trials_vec,time_to_reach,'b.');
    %plot(trials_vec,time_to_reach_corrected,'b.','color',[0.8 0.8 0.8]);
    %plot(trials_vec(idx_reach_syncfix),time_to_reach_inferred,'.','color',[0.5 0.5 0.5]);
    %plot(trials_vec(idx_reach_syncfix),time_to_reach_inferred_approxlog,'r.');
    plot(trials_vec(idx_reach_syncfix),time_to_reach_inferred_approxlog_corrected,'k.');
else
    plot(trials_vec,time_to_reach_corrected,'k.');
end

if exist('moved_left','var'), if ~isempty(moved_left), xline(moved_left,'-','left','Color',left_color,'linewidth',1.5,'fontsize',12); end; end
if exist('moved_right','var'), if ~isempty(moved_right), xline(moved_right,'-','right','Color',right_color,'linewidth',1.5,'fontsize',12); end; end
if exist('moved_center','var'), if ~isempty(moved_center), xline(moved_center,'-','center','Color',center_color,'linewidth',1.5,'fontsize',12); end; end
set(gca,axeOpt{:})
ylim([0 5])
xlim([0 last_trial_disp])
xlabel('trials across session'); ylabel('time to collect water (s)')


% play load cell aligned to init - valid PUSH
subplot(4,3,8)
if ~isempty(idx_comp_trial_push)
    plot(x_around_sec,paw_init_play_short(:,idx_comp_trial_push),'color',[push_clr 0.04],'linewidth',1); hold on
    plot(x_around_sec,nanmedian(paw_init_play_short(:,idx_comp_trial_push),2),'-','color',push_clr,'linewidth',2)
end
xline(0,'--','Color',cy,'linewidth',2);
ylim(ylim_loadcell); xlim([min(x_around_sec),max(x_around_sec)])
title('Valid Push')
xlabel('time from trial initiation (s)'); ylabel('force (mV/V)'); 
%xlabel('time from trial init (s)')
%ylabel('force (mV/V)'); title('Play load cell aligned to trial init');
set(gca,axeOpt{:})

% play load cell aligned to init - invalid PULL
subplot(4,3,9)
if ~isempty(inval_pull_ind)
    plot(x_around_sec,force_invPull(:,:),'color',[pull_clr 0.02],'linewidth',1); hold on
    plot(x_around_sec,nanmedian(force_invPull,2),'--','color',[pull_clr 0.8],'linewidth',2)
end
xline(0,'--','Color',cy,'linewidth',2);
ylim(ylim_loadcell); xlim([min(x_around_sec),max(x_around_sec)])
xlabel('time from trial initiation (s)'); ylabel('force (mV/V)'); 
title('Invalid Pull')
set(gca,axeOpt{:})

% play load cell aligned to init - valid PULL
subplot(4,3,11)
if ~isempty(idx_comp_trial_pull)
    plot(x_around_sec,paw_init_play_short(:,idx_comp_trial_pull),'color',[pull_clr 0.04],'linewidth',1); hold on
    plot(x_around_sec,nanmedian(paw_init_play_short(:,idx_comp_trial_pull),2),'color',pull_clr,'linewidth',2)
end
xline(0,'--','Color',cy,'linewidth',2);
ylim(ylim_loadcell); xlim([min(x_around_sec),max(x_around_sec)])
%xlabel('time from trial init (s)')
title('Valid Pull')
xlabel('time from trial initiation (s)'); ylabel('force (mV/V)'); 
%ylabel('force (mV/V)'); title('Play load cell aligned to trial init');
set(gca,axeOpt{:})

% play load cell aligned to init - invalid PUSH
subplot(4,3,12)
if ~isempty(inval_push_ind)
    plot(x_around_sec,force_invPush(:,:),'color',[push_clr 0.02],'linewidth',1); hold on
    plot(x_around_sec,nanmedian(force_invPush,2),'--','color',[push_clr 0.8],'linewidth',2)
end
xline(0,'--','Color',cy,'linewidth',2);
ylim(ylim_loadcell); xlim([min(x_around_sec),max(x_around_sec)])
xlabel('time from trial initiation (s)'); ylabel('force (mV/V)'); 
title('Invalid Push')
%ylabel('force (mV/V)'); title('Barrier load cell aligned to trial init');
set(gca,axeOpt{:})
%axis square

lw=2;
subplot(4,3,[3,5,2,6])
plot(units_slide_min,count_init_per_unit,'linewidth',6,'Color',[0 0 0 0.2]); hold on
plot(units_slide_min,count_valPull_per_unit,'-','linewidth',lw,'color',pull_clr);
plot(units_slide_min,count_valPush_per_unit,'-','linewidth',lw,'color',push_clr); hold on
plot(units_slide_min,count_invalPush_per_unit,'--','linewidth',lw,'color',[push_clr .8]);
plot(units_slide_min,count_invalPull_per_unit,'--','linewidth',lw,'color',[pull_clr .8]); hold off
legend('valid trials','valid pull','valid push','invalid push','invalid pull','box','off','AutoUpdate','off')
%legend('valid PUSH','valid PULL','invalid PUSH','valid PULL')
set(gca,axeOpt{:})
if ~isempty(transition_2push), xline(init_time(transition_2push)./60,'-','push','Color',[.5 .5 .5 .5],'linewidth',3,'fontsize',10,'LabelHorizontalAlignment','center'); end
if ~isempty(transition_2pull), xline(init_time(transition_2pull)./60,'-','pull','Color',[.5 .5 .5 .5],'linewidth',3,'fontsize',10,'LabelHorizontalAlignment','center'); end

if ~isempty(transition_2push), xline(init_time(moved_left)./60,'-.','left','Color',[.8 .8 .8],'linewidth',1,'fontsize',10,'LabelHorizontalAlignment','center'); end
if ~isempty(transition_2pull), xline(init_time(moved_right)./60,'-.','right','Color',[.8 .8 .8],'linewidth',1,'fontsize',10,'LabelHorizontalAlignment','center'); end
if ~isempty(transition_2pull), xline(init_time(moved_center)./60,'-.','center','Color',[.8 .8 .8],'linewidth',1,'fontsize',10,'LabelHorizontalAlignment','center'); end

%xlim([3, ceil(max(units_slide_min))])
xlim([2, last_trial_disp_time/60])

%title({'';sprintf('%s%.2f%s','Counts of valid/invalid push/pulls in a win of ',sliding_win./60,' minutes');})
ylabel(strcat('number of valid/invalid trials (',num2str(sliding_win/60),' min^{-1})')); 
xlabel('time across session (min)');

set(fig_performance,'position',[2454         272        1203         702],'color','w'); shg
title(sprintf('%s%s%s%s', 'Session: ',session, '   |    Nr of trials = ',num2str(nr_trials)));

shg

% SAVE 
file_name = 'daily_fig_performance';
if last_trial_disp ~= nr_trials
    file_name = strcat(file_name,'_until_trial',num2str(last_trial_disp));
end
save_path = rootdir;
eps_file = fullfile(save_path,[file_name,'.eps']);
png_file = fullfile(save_path,[file_name,'.png']);
saveas(fig_performance,png_file,'png');
%print(fig_performance,eps_file,'-depsc','-painters','-loose');
print(fig_performance, fullfile(save_path, [file_name,'.pdf']), '-dpdf', '-painters');


%% Techincal(ish) figure

fig_other=figure;

% sync
subplot(2,4,[1,2,5,6])
plot(time_newTrial_log,time_newTrial_cam,'-k'); hold on
plot(time_newTrial_log,time_newTrial_cam,'.');  hold off
xlabel('trial time on harp'); ylabel('trial time on  cam');
axis square
set(gca,axeOpt{:})
 title('Sync');

% ------------------
% ITI, REACH
% distribution of ITIs
subplot(2,4,4)
nbins = 30;
histogram(ITI_duration,nbins,'normalization','probability','facecolor',[.5 .5 .5]);
set(gca,axeOpt{:})
xlabel('ITI duration (sec)');
title('ITIs distribution')
%title({sprintf('%s%s%s',mouse,'  |  ',session);'ITI distribution'},'Interpreter', 'none')

% play load cell aligned to ITI
subplot(2,4,3)
% len_ti = numel(x_reach_sec);
% time_iti = linspace(0,1.5,len_ti);
len_ti=3001;
time_iti = linspace(0,3,len_ti);
plot(time_iti,paw_ITI_play(1:len_ti,:),'color',neutral_init_t); hold on
plot(time_iti,nanmedian(paw_ITI_play(1:len_ti,:),2),'color',neutral_init_m,'linewidth',1.5)
xline(0,'--','Color',cy,'linewidth',2);
ylim(ylim_loadcell); xlim([min(time_iti),max(time_iti)])
set(gca,axeOpt{:})
hold off
xlabel('time from ITI start (s)');
ylabel('force (mV/V)');
title('Load cell forces')

% load cell aligned to reach time
subplot(2,4,7)
plot(x_around_sec,force_aroundReach(:,:)','color',neutral_init_t); hold on
plot(x_around_sec,nanmedian(force_aroundReach,2),'color',neutral_init_m,'linewidth',1.5)
xline(0,'--','Color',cy,'linewidth',2);
ylim(ylim_loadcell)
set(gca,axeOpt{:})
hold off
xlabel('time from water reached (s)')
ylabel('force (mV/V)');

% sum of force per trial during ITI vs init
subplot(2,4,8)
histogram(sum_gen_force_init_early,'BinWidth',3e5,'normalization','probability','facecolor',clrs_blue);
hold on
histogram(sum_gen_force_ITI_later,'BinWidth',3e5,'normalization','probability','facecolor',[.5 .5 .5]);
set(gca,axeOpt{:})
xlim([-1e5 6e6])
legend('mov @ trial init','mov @ ITI','linewidth',.5,'box','off')
xlabel('sum of force per trial')
shg
set(fig_other,'position',[2454 442 1243 495],'color','w'); shg

% SAVE 
file_name2 = 'daily_fig_others';
%eps_file2 = fullfile(save_path,[file_name2,'.eps']);
png_file2 = fullfile(save_path,[file_name2,'.png']);
saveas(fig_other,png_file2,'png');
%print(fig_performance,eps_file2,'-depsc','-painters','-loose');

%% Save: behavior structure!
% Path and identity
behavior.mouse = mouse;
behavior.session = session;
behavior.group = group;
behavior.setup = setup;
behavior.path = rootdir;

% session general
behavior.nr_trials = nr_trials;
behavior.session_duration = session_duration;
behavior.behavior_duration.shorten_sess = shorten_display;
behavior.behavior_duration.trial_end = last_trial_disp;
behavior.behavior_duration.time_end = last_trial_disp_time; % timelog
behavior.behavior_duration.time_start = read_log(sess_start_ind,1); %harptime
behavior.flag_no_sync = flag_no_sync;

% inputs
behavior.inputs.real_log_all = read_log_all;
behavior.inputs.read_log = read_log;
behavior.inputs.timelog = timelog;
behavior.inputs.read_cam = read_cam;
behavior.inputs.paws_all = read_paws_all;
behavior.inputs.read_water = read_water;

% logs
behavior.logs.trial_available_all_ind = trial_available_all_ind;
behavior.logs.init_start_stop_ind = init_start_stop_ind;

behavior.logs.push_trial_ind = push_trial_ind;
behavior.logs.pull_trial_ind = pull_trial_ind;
behavior.logs.trial_init_ind = trial_init_ind;
behavior.logs.inval_push_ind = inval_push_ind;
behavior.logs.inval_pull_ind = inval_pull_ind;
behavior.logs.reached_ind = reached_ind;
behavior.logs.reached_really_ind = reached_really_ind;
behavior.logs.ITIs_start_ind = ITIs_start_ind;
behavior.logs.manual_water_ind = manual_water_ind;
behavior.logs.sess_start_ind = sess_start_ind;
behavior.logs.water_delivered_ind = water_delivered_ind;

% sync
behavior.sync.time_newTrial_log = time_newTrial_log;
behavior.sync.time_newTrial_cam = time_newTrial_cam;

% init trial
behavior.init.trial_init_ind = trial_init_ind;
behavior.init.is_push = is_push;
behavior.init.idx_trial_push = idx_comp_trial_push;
behavior.init.idx_trial_pull = idx_comp_trial_pull;
behavior.init.transition_2push = transition_2push;
behavior.init.transition_2pull = transition_2pull;
behavior.init.push_bounds = push_bounds;
behavior.init.pull_bounds = pull_bounds;
% init time of...
behavior.init.timeof.init_time = init_time;
behavior.init.timeof.valPush_time = valPush_time;
behavior.init.timeof.valPull_time = valPull_time;
behavior.init.timeof.invalPush_time = invalPush_time;
behavior.init.timeof.invalPull_time = invalPull_time;
% init times to (durations)...
behavior.init.timeto.time_to_init_push = time_to_init_push;
behavior.init.timeto.time_to_init_pull = time_to_init_pull;
% load cells around init
behavior.init.loadcells.x_around_sec = x_around_sec;
behavior.init.loadcells.paw_init_play_short = paw_init_play_short;
behavior.init.loadcells.force_invPull = force_invPull;
behavior.init.loadcells.force_invPush = force_invPush;

% reach trial
behavior.reach.left_idx = left_idx;                 
behavior.reach.right_idx = right_idx;                 
behavior.reach.center_idx = center_idx;  
behavior.reach.water_transition = water_transition;
behavior.reach.moved_left = moved_left;
behavior.reach.moved_right = moved_right;
behavior.reach.moved_center = moved_center;
behavior.reach.left_bounds = left_bounds;
behavior.reach.right_bounds = right_bounds;
behavior.reach.center_bounds = center_bounds;
behavior.reach.flag_reach_wrogly_detected = flag_reach_wrogly_detected;
behavior.reach.flag_realigned_rwd = flag_realigned_rwd;
behavior.reach.offset_water_detection = offset;
% reach time of...
behavior.reach.timeof.reach_time = reach_time;
behavior.reach.timeof.reach_time_corrected = reach_time_corrected;
if flag_realigned_rwd
    behavior.reach.timeof.reach_time_inferred = trial_reach_times_inferred;
    behavior.reach.timeof.reach_time_inferred_approxlog = trial_reach_times_inferred_approxlog;
    behavior.reach.timeof.reach_time_inferred_approxlog_corrected = reach_times_inferred_approxlog_corrected;
end
behavior.reach.timeof.moved_left_time = moved_left_time;
behavior.reach.timeof.moved_right_time = moved_right_time;
behavior.reach.timeof.moved_center_time = moved_center_time;
% reach times to (durations)...
behavior.reach.timeto.time_to_reach = time_to_reach;
behavior.reach.timeto.time_to_reach_corrected = time_to_reach_corrected;
if flag_realigned_rwd
    behavior.reach.timeto.time_to_reach_inferred = time_to_reach_inferred;
    behavior.reach.timeto.time_to_reach_inferred_approxlog  = time_to_reach_inferred_approxlog;
    behavior.reach.timeto.time_to_reach_inferred_approxlog_corrected  = time_to_reach_inferred_approxlog_corrected;
end
behavior.reach.timeto.time_to_reach_right = time_to_reach_right;
behavior.reach.timeto.time_to_reach_left = time_to_reach_left;
behavior.reach.timeto.time_to_reach_center = time_to_reach_center;
% load cells around reach
behavior.reach.loadcells.x_around_sec = x_around_sec;
behavior.reach.loadcells.force_aroundReach = force_aroundReach;

% ITI
behavior.iti.ITI_duration = ITI_duration;
behavior.iti.loadcells.time_iti = time_iti;
behavior.iti.loadcells.paw_ITI_play = paw_ITI_play;
behavior.iti.loadcells.sum_gen_force_init_early = sum_gen_force_init_early;
behavior.iti.loadcells.sum_gen_force_ITI_later = sum_gen_force_ITI_later;

% Counts of trials
behavior.counts.units_slide_min = units_slide_min;
behavior.counts.count_init_per_unit = count_init_per_unit;
behavior.counts.count_valPull_per_unit = count_valPull_per_unit;
behavior.counts.count_valPush_per_unit = count_valPush_per_unit;
behavior.counts.count_invalPush_per_unit = count_invalPush_per_unit;
behavior.counts.count_invalPull_per_unit = count_invalPull_per_unit;

% colrs to plot
behavior.colors.push_clr = push_clr;
behavior.colors.pull_clr = pull_clr;
behavior.colors.right_color = right_color;
behavior.colors.center_color = center_color;
behavior.colors.left_color = left_color;

save(strcat(rootdir,filesep,'behavior_session.mat'),'behavior');


