% Get data from individual session struct, join in a group big struct
clear; close all; clc

%% Load data
% Group and individual
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD';
group = '20230801_ChocolateGroup';
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


%% Create structure with all neurons
N_BG_all = [];
N_CT_all = [];
N_CB_all = [];

for m=1:length(animals)
    animal_idx = m;
    mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});
    fprintf('%s%s%s\n','Running ',mouse,'...')

    for s=1:size(all_sess{m},1)

        sess = all_sess{m}(s,:);
        fprintf('%s%s%s%s%s\n','Running mouse ',mouse,', session ',sess,':')

        % BG
        reg = 'BG';
        ephys_struct_dir = fullfile(rootdir,"ephys_and_behavior/","mat_files/",group,mouse,sess,reg);
        if ~exist(ephys_struct_dir,'dir')
            fprintf('%s%s%s\n','region ',reg,' has no data. Skipping...')
        else
            % Load mat files
            load(fullfile(ephys_struct_dir,['neu_',reg,'_struct.mat'])');
            N_BG_all = [N_BG_all; neu_strct];
            fprintf('%s%s%s\n','region ',reg,': data loaded and added!')
        end
        clear neu_strct

        % CT
        reg = 'CT';
        ephys_struct_dir = fullfile(rootdir,"ephys_and_behavior/","mat_files/",group,mouse,sess,reg);
        if ~exist(ephys_struct_dir,'dir')
            fprintf('%s%s%s\n','region ',reg,' has no data. Skipping...')
        else
            % Load mat files
            load(fullfile(ephys_struct_dir,['neu_',reg,'_struct.mat'])');
            N_CT_all = [N_CT_all; neu_strct];
            fprintf('%s%s%s\n','region ',reg,': data loaded and added!')
        end
        clear neu_strct

        reg = 'CB';
        ephys_struct_dir = fullfile(rootdir,"ephys_and_behavior/","mat_files/",group,mouse,sess,reg);
        if ~exist(ephys_struct_dir,'dir')
            fprintf('%s%s%s\n','region ',reg,' has no data. Skipping...')
        else
            % Load mat files
            load(fullfile(ephys_struct_dir,['neu_',reg,'_struct.mat'])');
            N_CB_all = [N_CB_all; neu_strct];
            fprintf('%s%s%s\n','region ',reg,': data loaded and added!')
        end
        clear neu_strct

    end
end
fprintf('Done, all neurons structs for each region are created!!')

%% Save mat
save_mat = fullfile(rootdir,"ephys_and_behavior/","mat_files/",group);
%save(fullfile(save_mat,'all_neurons_pooled.mat'),'N_BG_all','N_CT_all','N_CB_all','-v7.3');
save(fullfile(save_mat,'all_neurons_pooled.mat'),'N_CB_all','-append','-v7.3');

%%

% 
% neu_CT = neu(matches(neu.reg,"CTX"),:);
% 



