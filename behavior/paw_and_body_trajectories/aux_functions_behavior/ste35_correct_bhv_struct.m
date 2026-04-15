%% CORRECT THE BHV STRCUT BUG
% recalculate the is_hit, is_success etc
clear; close all; clc

%% %% Load data
% Group and individual
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\';
group = '20230801_ChocolateGroup';
group_bhv = '20230511_ChocolateGroup';
setup = 'headfixed_dynamicTarget';
animals = {...
    'CoteDor',...
    'Lindt',...
    'Toblerone',...
    'Milka',...
    'FerreroRocher'};
n_mice = length(animals);
n_sess = 7;

%%

for m=1:length(animals)
    %m=1;
    animal_idx = m;
    mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});
    fprintf('%s%s%s\n','Running ',mouse,'...')

    mouse_path = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse);

    for s=1:n_sess
        %s=1;
        sess = sprintf('%s%i','R',s);
        sess_path = fullfile(mouse_path,sess);

        % Run through sessions
        if exist(sess_path,'dir')~=0
            fprintf('%s%s:\n','Running session ',sess)

            reaching_dir = fullfile(rootdir,"behavior_data","analyzed_data","mat_files",group_bhv,setup,mouse,sess);
            load(fullfile(sess_path,'behavior_fundamentals.mat'))
            load(fullfile(reaching_dir,"session_reaching_data_paw.mat"),'reaches');

            reach_inVec_idx = bhv.reach_inVec_idx;
            bhv.hit_inVec = reaches.hit_reach(reach_inVec_idx);
            bhv.suc_inVec =  reaches.success_reach(reach_inVec_idx);

            save_mat = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse,sess);
            save(fullfile(save_mat,'behavior_fundamentals.mat'),'bhv');

        end
    end
end

disp('all correct bhv structs saved!')
