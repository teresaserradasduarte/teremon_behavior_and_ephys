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
keep_list = {'phyID','st_init', 'st_reach','brain_region','cell_type'};


%% Load
% load list of paths
load(fullfile(ephys_path,'teremon_paths.mat'))

% Initialize struct for cell types
all_CTX_sts = [];
all_CP_sts = [];

all_CB_sts = [];
all_DCN_sts = [];
all_BS_sts = [];
all_MB_sts = [];
all_P_sts = [];


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
            CTX_idx = neu.brain_region == "CTX";
            CP_idx = neu.brain_region == "CP";

            neu_CTX = neu(CTX_idx,:);
            neu_CP = neu(CP_idx,:);

            sts_CTX = table2struct(neu_CTX);
            sts_CP = table2struct(neu_XP);


            % Mat with all neurons together
            all_CTX_sts = [all_CTX_sts;sts_CTX];
            all_CP_sts = [all_CP_sts;sts_CP];


            clear eg_neurons neu_eg neu


        end % if exist

        if exist(fullfile(sess_path,'eg_neurons_CB.mat'),'file') ~=0

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
            BS_idx = neu.brain_region == "BS";
            MB_idx = neu.brain_region == "MB";
            P_idx = neu.brain_region == "P";

            neu_CB = neu(CB_idx,:);
            neu_DCN = neu(DCN_idx,:);
            neu_BS = neu(BS_idx,:);
            neu_MB = neu(MB_idx,:);
            neu_P = neu(P_idx,:);

            sts_CB = table2struct(neu_CB);
            sts_DCN = table2struct(neu_DCN);
            sts_BS = table2struct(neu_BS);
            sts_MB = table2struct(neu_MB);
            sts_P = table2struct(neu_P);


            % Mat with all neurons together
            all_CB_sts = [all_CB_sts;sts_CB];
            all_DCN_sts = [all_DCN_sts;sts_DCN];
            all_BS_sts = [all_BS_sts;sts_BS];
            all_MB_sts = [all_MB_sts;sts_MB];
            all_P_sts = [all_P_sts;sts_P];

            clear eg_neurons neu_eg neu


        end % if exist
    end % session
end % mouse

%% Save
load(fullfile(rootdir,"ephys_and_behavior","mat_files",group,'4_Milka/R1/eg_neurons_BG.mat'),'eg_neu_FR_params');
fprintf('%s','Saving struct extra-good neurons spike times (instrad of FR)...')
save_path = fullfile(rootdir,"ephys_and_behavior","mat_files",group);
save(fullfile(save_path,'st_eg_neurons'),...
    'all_CB_sts','all_DCN_sts','all_P_sts','all_BS_sts','all_MB_sts',...
    'all_CTX_sts','all_CP_sts',...
    'eg_neu_FR_params','-v7.3')
fprintf('%s','All done!!')
