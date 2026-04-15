%% Create a new region structre with selected fields: CB
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
keep_list = {'phyID','st_init', 'st_reach','brain_region','cell_type'};


%% Load
% load list of paths
load(fullfile(ephys_path,'teremon_paths.mat'))

% Initialize struct for cell types
all_CB_sts = [];
all_DCN_sts = [];

% Select mouse and session
for animal_idx = 1:n_mice
    %animal_idx = 5;

    for s = 1:n_sess
        %s = 1;
        mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});
        mouse_path = fullfile(rootdir,"ephys_and_behavior","mat_files",group,mouse);
        sess = sprintf('%s%i','R',s);
        sess_path = fullfile(mouse_path,sess);

        if exist(fullfile(sess_path,'eg_neurons_CB.mat'),'file') ~=0

            fprintf('%s%s%s%i\n','Running ',mouse,' session R',s)
            load(fullfile(sess_path,'eg_neurons_CB.mat'))
            load(fullfile(sess_path,'behavior_fundamentals.mat'))

            %Convert to table
            neu_eg = struct2table(eg_neurons);
            neu = neu_eg(:, keep_list);

            % Add behavior info and mouse / session id
            neu.idx_init_vPP_invPP = repmat(bhv.valPushPull_invalPushPull_idx,[1 size(neu,1)])';
            neu.DCnD_init = repmat(bhv.DCnD_init,[1 size(neu,1)])';
            neu.idx_reach_cat = repmat(bhv.cat_reach_inVec,[1 size(neu,1)])';
            neu.idx_reach_hit = repmat(bhv.hit_inVec,[1 size(neu,1)])';
            neu.idx_reach_succ = repmat(bhv.suc_inVec,[1 size(neu,1)])';
            neu.idx_reach_hDCnD =  repmat(bhv.DCnD_inVec_idx,[1 size(neu,1)])';
            neu.idx_reach_LCR =  repmat(bhv.LCR_inVec_idx,[1 size(neu,1)])';
            neu.idx_reach_PP =  repmat(bhv.pp_inVec_idx,[1 size(neu,1)])';
            neu.reach_px = repmat({bhv.reaches_inVec_px},size(neu,1),1);
            neu.mouse = repmat(mouse,[size(neu,1),1]);
            neu.sess = repmat(sess,[size(neu,1),1]);

            % Neu struct for each region, less info but with bhv
            CB_idx = neu.brain_region == "CB";
            DCN_idx = neu.brain_region == "DCN";

            neu_CB = neu(CB_idx,:);
            neu_DCN = neu(DCN_idx,:);

            sts_CB = table2struct(neu_CB);
            sts_DCN = table2struct(neu_DCN);


            % Mat with all neurons together
            all_CB_sts = [all_CB_sts;sts_CB];
            all_DCN_sts = [all_DCN_sts;sts_DCN];

            clear eg_neurons neu_eg neu


        end % if exist
    end % session
end % mouse

%% Save
load(fullfile(rootdir,"ephys_and_behavior","mat_files",group,'4_Milka/R1/eg_neurons_BG.mat'),'eg_neu_FR_params');
fprintf('%s','Saving struct with CB and DCN extra-good neurons spike times...')
save_path = fullfile(rootdir,"ephys_and_behavior","mat_files",group);
save(fullfile(save_path,'st_CB_DCN'),...
    'all_CB_sts','all_DCN_sts',...
    'eg_neu_FR_params','-v7.3')
fprintf('%s','All done!!')
