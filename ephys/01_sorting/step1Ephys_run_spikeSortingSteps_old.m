%% EPHYS PREPROCESSING --> ALIGN AND SPIKE SORTING
% 1 - run CatGT: correct time shifts, optionally CAR for artifact removal
% 2 - run kilisort to get units' clusters
clear; close all; clc

%% Select dataset & probe
% copy dataset fom local disk/cloud to SSD (ephys_root)
%probe = '1'; %imec1
probe = '0'; %imec0
run_catGT = 1;
run_KS = 0;
saveKS_to_matlab = 0;
create_dest_folder = 1;
double_spike_removal_flag = 1;
dest_folder_name = strcat('catGT');
tic
%% Code and data paths - CHANGE ACCORDINGLY
kilosort_root = 'C:\Users\Teresa\Documents\MATLAB\Kilosort-2.5'; % path to kilosort folder
npy_root = 'C:\Users\Teresa\Documents\MATLAB\npy-matlab'; % path to npy folder
catgt_root = 'C:\Users\Teresa\Documents\CatGTWinApp\CatGT-win'; % path to catGT folder
ephys_root = 'E:'; % source data
pathToYourConfigFile = strcat(kilosort_root,filesep,'configFiles'); % path to config files

% Add to path
addpath(genpath(kilosort_root))
addpath(npy_root)

% Find dataset in path
ephys_folder = '20082023_Ferrero_StrCer_S6_g0';%'26082023_Lindt_StrCer_S4_g0';
%ephys_to_run = dir(ephys_root);
%ephys_folder = ephys_to_run(4).name; %ephys_folder = ephys_to_run(3:end).name;
run_name = ephys_folder(1:end-3); % date_mouse_recordingSite_session

% Destination folder (pre-processed)
rawdata_folder = strcat(ephys_root,filesep,ephys_folder,filesep,ephys_folder,'_imec',probe);

%rawdata_folder = strcat(ephys_root,filesep,ephys_folder,filesep,ephys_folder,'_imec',probe);
if create_dest_folder == 1
    destination_folder = strcat(rawdata_folder,filesep,dest_folder_name);

    if ~exist(destination_folder,'dir'), mkdir(destination_folder); end
else
    destination_folder = rawdata_folder;
end

% -------------------------------------------------------------------------
%% run CatGT  -------------------------------------------------------------
if run_catGT==1
    %% Params for running CatGT
    % Parameters to choose
    %artifact_correction = ' -gblcar';
    %artifact_correction = ' -gbldmx';
    artifact_correction = '';
    %ap_filter = ' -apfilter=butter,6,300,10000';
    ap_filter = '';
    destination_folder_cmd = sprintf('%s%s%s',' -dest=',destination_folder,...
        ' -no_catgt_fld');
    %destination_folder_cmd = '';

    % fixed
    catgt_cmd = strcat(catgt_root,filesep,'CatGT');
    g='0';
    t ='0';

    %% Send messages to command line
    cmd_catGT = sprintf('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s',catgt_cmd,' ',...
        '-dir=',ephys_root,' -run=',run_name,...
        ' -g=',g,' -t=',t,' -prb_fld -prb=',probe,' -ap',...
        ap_filter,artifact_correction,destination_folder_cmd);
    status=system(cmd_catGT);
end


% -------------------------------------------------------------------------
%% run KILOSORT -----------------------------------------------------------
% Paths
rootZ = destination_folder; % source data to sort
rootH = rootZ; % path to temporary binary file

if run_KS==1
    %% Configurations for kilosort
    chanMapFile = 'neuropixPhase3B1_kilosortChanMap.mat'; % channel map NP1.0
    configFile_name = 'StandardConfig_TSD.m';

    ops.trange    = [0 Inf]; % time range to sort
    ops.NchanTOT  = 385; % total number of channels in your recording

    run(fullfile(pathToYourConfigFile, configFile_name))
    ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
    ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);

    %% Run kilosort
    % find the binary file
    fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
    ops.fbinary = fullfile(rootZ, fs(1).name);

    % Preprocessing
    % preprocess data to create temp_wh.dat
    rez = preprocessDataSub(ops);
    % do data registration
    rez = datashift2(rez, 1); % last input is for shifting data
    % order of batches random, controlled by random number generator
    iseed = 1;

    % Kilosort
    % main tracking and template matching algorithm
    rez = learnAndSolve8b(rez, iseed);

    % remove double-counted spikes
    overlap_s = 0.166e-3;
    channel_separation_um = 60;
    if double_spike_removal_flag == 1
        rez = remove_ks2_duplicate_spikes(rez,'overlap_s',overlap_s,'channel_separation_um',char(60));
    end

    % final merges
    rez = find_merges(rez, 1);

    % final splits by SVD
    rez = splitAllClusters(rez, 1);

    % decide on cutoff
    rez = set_cutoff(rez);

    % eliminate widely spread waveforms (likely noise)
    rez.good = get_good_units(rez);
    fprintf('found %d good units \n', sum(rez.good>0))

    % Saving data to Phy
    % write to Phy
    fprintf('Saving results to Phy  \n')
    rezToPhy(rez, rootZ);

    %% if you want to save the results to a Matlab file...
    if saveKS_to_matlab == 1
        % discard features in final rez file (too slow to save)
        rez.cProj = [];
        rez.cProjPC = [];

        % final time sorting of spikes, for apps that use st3 directly
        [~, isort]   = sortrows(rez.st3);
        rez.st3      = rez.st3(isort, :);

        % Ensure all GPU arrays are transferred to CPU side before saving to .mat
        rez_fields = fieldnames(rez);
        for i = 1:numel(rez_fields)
            field_name = rez_fields{i};
            if(isa(rez.(field_name), 'gpuArray'))
                rez.(field_name) = gather(rez.(field_name));
            end
        end

        % save final results as rez2
        fprintf('Saving final results in rez2  \n')
        fname = fullfile(rootZ, 'rez2.mat');
        save(fname, 'rez', '-v7.3');
    end
end
toc
