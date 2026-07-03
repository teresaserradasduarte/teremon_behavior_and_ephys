%% Create struct with st for every region instead of FR
clear; close all; clc

%% Load data
% Group and individual
rootdir = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\';
group = '20230801_ChocolateGroup';
ephys_path = fullfile(rootdir, "ephys_data/",group);
animals = {...
    'CoteDor',...
    'Lindt',...
    'Toblerone',...
    'Milka',...
    'FerreroRocher'};
n_mice = length(animals);
n_sess = 7;

%% Select parameters: fileds to keep
keep_list = {'phyID','st_init', 'st_reach','st','brain_region'};


%% Load
% load list of paths
load(fullfile(ephys_path,'teremon_paths.mat'))

% Initialize struct for cell types
all_CTX_sts = [];
all_CP_sts = [];

% Select mouse and session
for animal_idx = 1:n_mice
    %animal_idx = 5;

    for s = 1:n_sess
        %s = 1;
        mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});
        mouse_path = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse);
        sess = sprintf('%s%i','R',s);
        sess_path = fullfile(mouse_path,sess);

        fprintf('%s%s%s%i\n','Running ',mouse,' session R',s)


        if exist(fullfile(sess_path,'eg_neurons_BG.mat'),'file') ~=0

            load(fullfile(sess_path,'eg_neurons_BG.mat'))
            load(fullfile(sess_path,'behavior_fundamentals.mat'))

            %Convert to table
            neu_eg = struct2table(eg_neurons);
            neu = neu_eg(:, keep_list);

            % Add behavior info and mouse / session id
            neu.mouse = repmat(mouse,[size(neu,1),1]);
            neu.sess = repmat(sess,[size(neu,1),1]);

            % Neu struct for each region, less info but with bhv
            CTX_idx = neu.brain_region == "CTX";
            CP_idx = neu.brain_region == "CP";

            neu_CTX = neu(CTX_idx,:);
            neu_CP = neu(CP_idx,:);

            sts_CTX = table2struct(neu_CTX);
            sts_CP = table2struct(neu_CP);


            % Mat with all neurons together
            all_CTX_sts = [all_CTX_sts;sts_CTX];
            all_CP_sts = [all_CP_sts;sts_CP];


            clear eg_neurons neu_eg neu


        end % if exist

      
    end % session
end % mouse

%% Save
load(fullfile(rootdir,"ephys_and_behavior","mat_files",group,'4_Milka/R1/eg_neurons_BG.mat'),'eg_neu_FR_params');
fprintf('%s','Saving struct extra-good neurons spike times (instrad of FR)...')
save_path = fullfile(rootdir,"ephys_and_behavior","mat_files",group);
save(fullfile(save_path,'st_eg_neurons_GIL'),...
    'all_CTX_sts','all_CP_sts',...
    'eg_neu_FR_params','-v7.3')
fprintf('%s','All done!!')
