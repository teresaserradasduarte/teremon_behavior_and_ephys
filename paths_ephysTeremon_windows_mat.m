% Teremon_Path_Builder
clear; close all; clc

%% 1. Define Root Path
root_path = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\ephys_data\20230801_ChocolateGroup';

%% 2. Define Sub-paths (CB)
% Note: Using cell array {...} for lists of strings
subpaths_CB = {
    % Côte d'Or
    '1_CoteDor\02082023_CoteDorS1CB_g0\02082023_CoteDorS1CB_g0_imec0\catGT\kilosort4';
    '1_CoteDor\04082023_CoteDor_StrCer_S4_g0\04082023_CoteDor_StrCer_S4_g0_imec0\ibl_sorter_results';
    '1_CoteDor\04082023_CoteDor_Cer_S5_g0\04082023_CoteDor_Cer_S5_g0_imec0\ibl_sorter_results_driftAdapt40';
    '1_CoteDor\08082023_CoteDor_StrCer_S6_g0\08082023_CoteDor_StrCer_S6_g0_imec0\catGT\kilosort4';
    '1_CoteDor\08082023_CoteDor_StrCer_S7_g0\08082023_CoteDor_StrCer_S7_g0_imec0\ibl_sorter_results';
    % Lindt
    '2_Lindt\23082023_Lindt_StrCer_S1_true_g0\23082023_Lindt_StrCer_S1_true_g0_imec1\catGT\kilosort4';
    '2_Lindt\23082023_Lindt_StrCer_S2_g0\23082023_Lindt_StrCer_S2_g0_imec1\catGT\kilosort4';
    '2_Lindt\25082023_Lindt_StrCer_S3_g0\25082023_Lindt_StrCer_S3_g0_imec1\ibl_sorter_results';
    '2_Lindt\26082023_Lindt_StrCer_S4_g0\26082023_Lindt_StrCer_S4_g0_imec1\catGT\kilosort4';
    '2_Lindt\27082023_Lindt_StrCer_S5_g0\27082023_Lindt_StrCer_S5_g0_imec1\ibl_sorter_results';
    '2_Lindt\28082023_Lindt_StrCer_S6_g0\28082023_Lindt_StrCer_S6_g0_imec1\ibl_sorter_results';
    '2_Lindt\29082023_Lindt_StrCer_S7_g0\29082023_Lindt_StrCer_S7_g0_imec1\ibl_sorter_results_driftAdapt20';
    % Toblerone
    '3_Toblerone\23082023_Toblerone_StrCer_S1_g0\23082023_Toblerone_StrCer_S1_g0_imec1\ibl_sorter_results';
    '3_Toblerone\24082023_Toblerone_StrCer_S2_g0\24082023_Toblerone_StrCer_S2_g0_imec1\ibl_sorter_results';
    '3_Toblerone\25082023_Toblerone_StrCer_S3_g0\25082023_Toblerone_StrCer_S3_g0_imec1\catGT\kilosort4';
    '3_Toblerone\26082023_Toblerone_StrCer_S4_g0\26082023_Toblerone_StrCer_S4_g0_imec1\ibl_sorter_results';
    '3_Toblerone\27082023_Toblerone_StrCer_S5_g0\27082023_Toblerone_StrCer_S5_g0_imec1\ibl_sorter_results';
    '3_Toblerone\28082023_Toblerone_StrCer_S6_g0\28082023_Toblerone_StrCer_S6_g0_imec1\catGT\kilosort4';
    % Milka
    '4_Milka\15082023_Milka_StrCer_S1_g0\15082023_Milka_StrCer_S1_g0_imec0\ibl_sorter_results';
    '4_Milka\18082023_Milka_StrCer_S4_g0\18082023_Milka_StrCer_S4_g0_imec0\catGT\kilosort4';
    '4_Milka\19082023_Milka_StrCer_S5_g0\19082023_Milka_StrCer_S5_g0_imec0\ibl_sorter_results';
    '4_Milka\20082023_Milka_StrCer_S6_g0\20082023_Milka_StrCer_S6_g0_imec0\ibl_sorter_results';
    '4_Milka\21082023_Milka_StrCer_S7_g0\21082023_Milka_StrCer_S7_g0_imec0\ibl_sorter_results';
    % Ferrero Rocher
    '5_FerreroRocher\15082023_Ferrero_StrCer_S1_g0\15082023_Ferrero_StrCer_S1_g0_imec0\catGT\kilosort4';
    '5_FerreroRocher\18082023_Ferrero_StrCer_S4_g0\18082023_Ferrero_StrCer_S4_g0_imec0\catGT\kilosort4';
    '5_FerreroRocher\19082023_Ferrero_StrCer_S5_g0\19082023_Ferrero_StrCer_S5_g0_imec0\ibl_sorter_results';
    '5_FerreroRocher\20082023_Ferrero_StrCer_S6_g0\20082023_Ferrero_StrCer_S6_g0_imec0\catGT\kilosort4';
};

%% 3. Define Sub-paths (BG)
subpaths_BG = {
    % Côte d'Or
    '1_CoteDor\02082023_CoteDorS2Str_g0\02082023_CoteDorS2Str_g0_imec0\ibl_sorter_results';
    '1_CoteDor\04082023_CoteDor_Str_S3_g0\04082023_CoteDor_Str_S3_g0_imec0\ibl_sorter_results';
    '1_CoteDor\04082023_CoteDor_StrCer_S4_g0\04082023_CoteDor_StrCer_S4_g0_imec1\ibl_sorter_results_driftAdapt2';
    '1_CoteDor\08082023_CoteDor_StrCer_S6_g0\08082023_CoteDor_StrCer_S6_g0_imec1\catGT\kilosort4';
    '1_CoteDor\08082023_CoteDor_StrCer_S7_g0\08082023_CoteDor_StrCer_S7_g0_imec1\ibl_sorter_results';
    % Lindt
    '2_Lindt\23082023_Lindt_StrCer_S1_true_g0\23082023_Lindt_StrCer_S1_true_g0_imec0\catGT\kilosort4';
    '2_Lindt\23082023_Lindt_StrCer_S2_g0\23082023_Lindt_StrCer_S2_g0_imec0\catGT\kilosort4';
    '2_Lindt\25082023_Lindt_StrCer_S3_g0\25082023_Lindt_StrCer_S3_g0_imec0\ibl_sorter_results';
    '2_Lindt\26082023_Lindt_StrCer_S4_g0\26082023_Lindt_StrCer_S4_g0_imec0\catGT\kilosort4';
    '2_Lindt\27082023_Lindt_StrCer_S5_g0\27082023_Lindt_StrCer_S5_g0_imec0\ibl_sorter_results';
    '2_Lindt\28082023_Lindt_StrCer_S6_g0\28082023_Lindt_StrCer_S6_g0_imec0\ibl_sorter_results';
    '2_Lindt\29082023_Lindt_StrCer_S7_g0\29082023_Lindt_StrCer_S7_g0_imec0\ibl_sorter_results';
    % Toblerone
    '3_Toblerone\23082023_Toblerone_StrCer_S1_g0\23082023_Toblerone_StrCer_S1_g0_imec0\ibl_sorter_results';
    '3_Toblerone\24082023_Toblerone_StrCer_S2_g0\24082023_Toblerone_StrCer_S2_g0_imec0\ibl_sorter_results';
    '3_Toblerone\25082023_Toblerone_StrCer_S3_g0\25082023_Toblerone_StrCer_S3_g0_imec0\catGT\kilosort4';
    '3_Toblerone\26082023_Toblerone_StrCer_S4_g0\26082023_Toblerone_StrCer_S4_g0_imec0\catGT\kilosort4';
    '3_Toblerone\27082023_Toblerone_StrCer_S5_g0\27082023_Toblerone_StrCer_S5_g0_imec0\ibl_sorter_results';
    '3_Toblerone\28082023_Toblerone_StrCer_S6_g0\28082023_Toblerone_StrCer_S6_g0_imec0\catGT\kilosort4';
    % Milka
    '4_Milka\15082023_Milka_StrCer_S1_g0\15082023_Milka_StrCer_S1_g0_imec1\ibl_sorter_results';
    '4_Milka\18082023_Milka_StrCer_S4_g0\18082023_Milka_StrCer_S4_g0_imec1\catGT\kilosort4';
    '4_Milka\19082023_Milka_StrCer_S5_g0\19082023_Milka_StrCer_S5_g0_imec1\ibl_sorter_results';
    '4_Milka\20082023_Milka_StrCer_S6_g0\20082023_Milka_StrCer_S6_g0_imec1\ibl_sorter_results';
    '4_Milka\21082023_Milka_StrCer_S7_g0\21082023_Milka_StrCer_S7_g0_imec1\ibl_sorter_results';
    % Ferrero Rocher
    '5_FerreroRocher\15082023_Ferrero_StrCer_S1_g0\15082023_Ferrero_StrCer_S1_g0_imec1\catGT\kilosort4';
    '5_FerreroRocher\18082023_Ferrero_StrCer_S4_g0\18082023_Ferrero_StrCer_S4_g0_imec1\catGT\kilosort4';
    '5_FerreroRocher\19082023_Ferrero_StrCer_S5_g0\19082023_Ferrero_StrCer_S5_g0_imec1\ibl_sorter_results';
    '5_FerreroRocher\20082023_Ferrero_StrCer_S6_g0\20082023_Ferrero_StrCer_S6_g0_imec1\catGT\kilosort4';
};

%% 4. Combine into Full Paths
% This creates a new cell array where every entry is the full absolute path
files_curation_CB = fullfile(root_path, subpaths_CB);
files_curation_BG = fullfile(root_path, subpaths_BG);

%% 5. Save to MAT file
% Saves only the final full path variables
save('teremon_paths.mat', 'files_curation_CB', 'files_curation_BG');

disp('Paths saved to teremon_paths.mat');