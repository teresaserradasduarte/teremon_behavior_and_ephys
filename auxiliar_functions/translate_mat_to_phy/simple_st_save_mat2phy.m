%% Simplify dataset - just spike times and behavior logs & converto to python
clear; close all; clc

%% Load
% ehpys
ephys_root = 'E:\'; %group_ephys = '20230801_ChocolateGroup';
ephys_sess = '26082023_Lindt_StrCer_S4_g0';

% Behavior
behavior_root ='D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\data\TD\behavior_data\raw_data';
group_setup_behav = strcat('20230511_ChocolateGroup',filesep,'headfixed_dynamicTarget');
mouse = '2_Lindt';
session = 'R4';
imec = '0';

% Load behavior
behavior_path = strcat(behavior_root,filesep,group_setup_behav,filesep,mouse,filesep,session);
load(strcat(behavior_path,filesep,'behavior_session.mat'));

% Load ehpys
imec_path = strcat(ephys_root,filesep,ephys_sess,filesep,ephys_sess,'_imec',imec,filesep,'catGT',filesep,'kilosort4');
neurons_imec = load(strcat(imec_path,filesep,'neurons_session.mat'));

% Create folder
simplified_data = strcat(imec_path,filesep,'simplified_data');
if ~exist('simplified_data','dir'), mkdir(simplified_data); end

%% Get behavior logs and windows
behav_info.was_sess_shortened = boolean(behavior.behavior_duration.shorten_sess);
behav_info.last_trial_considered = behavior.behavior_duration.trial_end;
behav_info.time_win =  [behavior.behavior_duration.time_start behavior.behavior_duration.time_end];
behav_info.logs = behavior.inputs.real_log_all;


%% Get the spike_times
neu = struct2table(neurons_imec.neurons);
good_idx = find(neu.quality == 2);
neurons_spk_times = neu(good_idx,["st_raw_all","st","phyID","extraGood"]);
struct_spikes = table2struct(neurons_spk_times);

%% Neurons and behavior 
save(strcat(simplified_data,filesep,'behav_info.mat'),'behav_info');
save(strcat(simplified_data,filesep,'struct_spikes.mat'), 'struct_spikes');

%writetable(neurons_spk_times, strcat(simplified_data,filesep,'neurons_spk_times.csv'));



