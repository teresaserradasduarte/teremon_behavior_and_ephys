% Get extra-good neurons, calculate instantaneous FR, order neurons
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
    ['R3'],...
    ['R4'],...
    ['R1';'R4';'R6']};

reg = 'BG';

% Initialize
SVM_sess = struct();
cc = 0; % counter of sessions (all mice)

% Save mat dir
save_mat = fullfile(rootdir,"ephys_and_behavior","mat_files",group);

%% Window
win_int = [-2.5 2.5];
tm_before = win_int(1);
tm_after = win_int(2); % sec
win_dur = tm_after-tm_before;
bin_width = 0.1;
bin_step = 0.025;
n_bins = floor((win_dur/bin_step)-(bin_width/bin_step))+1;

%% Figure defs
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',...
    0.05,'ticklength',[1,1]*.01,'fontsize',10, 'TickDir', 'out'};
clr_center = [72 60 50]./256;
clr_dom = [48 131 220]./256;
clr_nondom = [158 216 219]./256;
clr_init = [63 130 109]./256;
clrs = [clr_init;clr_dom;clr_center;clr_nondom];

%% SVM - Lood through mice and session
for m=1:length(animals)
    animal_idx = m;
    mouse = sprintf('%i_%s',animal_idx,animals{animal_idx});
    fprintf('%s%s%s\n','Looping through ',mouse,':')

    for s=1:size(all_sess{m},1)
        sess = all_sess{m}(s,:);
        fprintf('%s%s%s%s%s\n','Running mouse ',mouse,', session ',sess,':')

        ephys_struct_dir = fullfile(rootdir,"ephys_and_behavior/","mat_files/",group,mouse,sess,reg);
        if ~exist(ephys_struct_dir,'dir')
            fprintf('%s%s%s\n','region ',reg,' has no data. Skipping...')
        else

            % Load mat files
            load(fullfile(ephys_struct_dir,['neu_',reg,'_struct.mat'])');
            cc = cc+1;

            % SVM struct
            SVM_sess(cc).mouse = neu_strct.mouse;
            SVM_sess(cc).sess = neu_strct.sess;

            % Number of neurons
            n_neu = length(neu_strct);
            SVM_sess(cc).n_neu = n_neu;


            %%
            i_pp_trials_all = neu_strct.idx_init_vPP_invPP;
            i_DCnD_trials_all = neu_strct.DCnD_init;
            r_hit_trials_all = neu_strct.idx_reach_hit;
            r_DCnD_trials_all = neu_strct.idx_reach_hDCnD;
            r_pp_trials_all = neu_strct.idx_reach_PP;

            n_trials_pp = length(i_pp_trials_all);
            n_trials_rr = length(r_hit_trials_all);

            spk_count_init = nan(n_bins, n_trials_pp,n_neu);
            spk_count_reach = nan(n_bins, n_trials_rr,n_neu);


            for n=1:n_neu

                bin_edges = nan(n_bins,2);
                sts_init = neu_strct(n).st_init;
                sts_reach = neu_strct(n).st_reach;

                for i = 1:n_bins
                    bin_start = tm_before + (i-1)*bin_step;
                    bin_stop = tm_before + bin_width + (i-1)*bin_step;
                    for tt = 1:n_trials_pp
                        spk_count_init(i,tt,n) = numel(find(sts_init{tt}>bin_start & sts_init{tt}<bin_stop));
                    end

                    for tt = 1:n_trials_rr
                        spk_count_reach(i,tt,n) = numel(find(sts_reach{tt}>bin_start & sts_reach{tt}<bin_stop));
                    end
                    bin_edges(i,1) = bin_start;
                    bin_edges(i,2) = bin_stop;
                end
            end
            SVM_sess(cc).bin_edges = bin_edges;

            %% Init time
            s1_i_pp = spk_count_init(:,i_pp_trials_all == 0,:);
            s2_i_pp = spk_count_init(:,i_pp_trials_all == 1,:);
            s1_i_D = spk_count_init(:,i_DCnD_trials_all == 1,:);
            s2_i_D = spk_count_init(:,i_DCnD_trials_all == 2 | i_DCnD_trials_all == 3,:);
            s1_i_C = spk_count_init(:,i_DCnD_trials_all == 2,:);
            s2_i_C = spk_count_init(:,i_DCnD_trials_all == 1 | i_DCnD_trials_all == 3,:);
            s1_i_nD = spk_count_init(:,i_DCnD_trials_all == 3,:);
            s2_i_nD = spk_count_init(:,i_DCnD_trials_all == 1 | i_DCnD_trials_all == 2,:);

            % Reach time
            s1_r_pp = spk_count_reach(:,r_pp_trials_all == 0 & r_hit_trials_all == 1,:);
            s2_r_pp = spk_count_reach(:,r_pp_trials_all == 1 & r_hit_trials_all == 1,:);
            s1_r_D = spk_count_reach(:,r_DCnD_trials_all == 1 & r_hit_trials_all == 1,:);
            s2_r_D = spk_count_reach(:,(r_DCnD_trials_all == 2 | r_DCnD_trials_all == 3) ...
                & r_hit_trials_all == 1,:);
            s1_r_C = spk_count_reach(:,r_DCnD_trials_all == 2 & r_hit_trials_all == 1,:);
            s2_r_C = spk_count_reach(:,(r_DCnD_trials_all == 1 | r_DCnD_trials_all == 3) ...
                & r_hit_trials_all == 1,:);
            s1_r_nD = spk_count_reach(:,r_DCnD_trials_all == 3 & r_hit_trials_all == 1,:);
            s2_r_nD = spk_count_reach(:,(r_DCnD_trials_all == 1 | r_DCnD_trials_all == 2) ...
                & r_hit_trials_all == 1,:);

            SVM_sess(cc).init_s1s2.s1_i_pp = s1_i_pp;
            SVM_sess(cc).init_s1s2.s2_i_pp = s2_i_pp;
            SVM_sess(cc).init_s1s2.s1_i_D = s1_i_D;
            SVM_sess(cc).init_s1s2.s2_i_D = s2_i_D;
            SVM_sess(cc).init_s1s2.s1_i_C = s1_i_C;
            SVM_sess(cc).init_s1s2.s2_i_C = s2_i_C;
            SVM_sess(cc).init_s1s2.s1_i_nD = s1_i_nD;
            SVM_sess(s).init_s1s2.s2_i_nD = s2_i_nD;

            SVM_sess(cc).reach_s1s2.s1_r_pp = s1_r_pp;
            SVM_sess(cc).reach_s1s2.s2_r_pp = s2_r_pp;
            SVM_sess(cc).reach_s1s2.s1_r_D = s1_r_D;
            SVM_sess(cc).reach_s1s2.s2_r_D = s2_r_D;
            SVM_sess(cc).reach_s1s2.s1_r_C = s1_r_C;
            SVM_sess(cc).reach_s1s2.s2_r_C = s2_r_C;
            SVM_sess(cc).reach_s1s2.s1_r_nD = s1_r_nD;
            SVM_sess(cc).reach_s1s2.s2_r_nD = s2_r_nD;

            %% %% RUN SVM DECODER

            % Parameters
            ratio_train_val = 0.8;
            ncv = 100;
            nfold = 10;
            Cvec = [0.001, 0.01, 0.1, 1, 10, 100, 1000];

            % store svm classifier accuracy
            bac_init = nan(n_bins,4);
            bac_reach = nan(n_bins,4);

            tic
            fprintf('SVM - running 8 classifiers...\n');

            % Run by bins
            for b = 1:n_bins
                % INIT
                [bac_init(b,1)] = svm_simple_fun_adaptedfromVK(...
                    squeeze(s1_i_pp(b,:,:)),...
                    squeeze(s2_i_pp(b,:,:)),...
                    ratio_train_val,ncv,nfold,Cvec);
                [bac_init(b,2)] = svm_simple_fun_adaptedfromVK(...
                    squeeze(s1_i_D(b,:,:)),...
                    squeeze(s2_i_D(b,:,:)),...
                    ratio_train_val,ncv,nfold,Cvec);
                [bac_init(b,3)] = svm_simple_fun_adaptedfromVK(...
                    squeeze(s1_i_C(b,:,:)),...
                    squeeze(s2_i_C(b,:,:)),...
                    ratio_train_val,ncv,nfold,Cvec);
                [bac_init(b,4)] = svm_simple_fun_adaptedfromVK(...
                    squeeze(s1_i_nD(b,:,:)),...
                    squeeze(s2_i_nD(b,:,:)),...
                    ratio_train_val,ncv,nfold,Cvec);

                % REACH
                [bac_reach(b,1)] = svm_simple_fun_adaptedfromVK(...
                    squeeze(s1_r_pp(b,:,:)),...
                    squeeze(s2_r_pp(b,:,:)),...
                    ratio_train_val,ncv,nfold,Cvec);
                [bac_reach(b,2)] = svm_simple_fun_adaptedfromVK(...
                    squeeze(s1_r_D(b,:,:)),...
                    squeeze(s2_r_D(b,:,:)),...
                    ratio_train_val,ncv,nfold,Cvec);
                [bac_reach(b,3)] = svm_simple_fun_adaptedfromVK(...
                    squeeze(s1_r_C(b,:,:)),...
                    squeeze(s2_r_C(b,:,:)),...
                    ratio_train_val,ncv,nfold,Cvec);
                [bac_reach(b,4)] = svm_simple_fun_adaptedfromVK(...
                    squeeze(s1_r_nD(b,:,:)),...
                    squeeze(s2_r_nD(b,:,:)),...
                    ratio_train_val,ncv,nfold,Cvec);

                fprintf('\rRunning bin %d of %d (%.1f%%)', b, n_bins, 100 * b / n_bins);
                %toc
            end
            fprintf('\ndone!!\n')
            toc

            % Save
            SVM_sess(cc).bac_init = bac_init;
            SVM_sess(cc).bac_reach = bac_reach;

            %% FIGURE

            %Save out dir
            save_out_dir = fullfile(rootdir,"ephys_and_behavior/","out_files/",group,mouse,sess,reg,'SVM');
            if ~exist(save_out_dir,"dir"), mkdir(save_out_dir); end

            %%
            bin_mean=bin_edges(:,1)+diff(bin_edges(1:2));

            figure()
            subplot(121)
            for plt=1:4
                plot(bin_mean,bac_init(:,plt),'color',clrs(plt,:),'linewidth',2); hold on
            end
            xlabel('time from trial init (s)'); ylabel('accuracy');
            set(gca,axeOpt{:})
            yline(.5,'--','color',[.7 .7 .7 .7],'LineWidth',1)
            xlim(win_int); ylim([.45 1])
            title(mouse,'Interpreter','none');

            subplot(122)
            for plt=1:4
                plot(bin_mean,bac_reach(:,plt),'color',clrs(plt,:),'linewidth',2); hold on
            end
            xlabel('time from reach endpoint (s)'); ylabel('accuracy');
            set(gca,axeOpt{:})
            yline(.5,'--','color',[.7 .7 .7 .7],'LineWidth',1)
            xlim(win_int); ylim([.45 1])
            title(sess)

            set(gcf,'Position',[680         585        1226         393],'Color','w');
            saveas(gcf,strcat(save_out_dir,filesep,'svm_accuracy.png'),'png');

            % Params
            SVM_sess(cc).params.win_int = win_int;
            SVM_sess(cc).params.ratio_train_val = ratio_train_val;
            SVM_sess(cc).params.ncv = ncv;
            SVM_sess(cc).params.nfold = nfold;
            SVM_sess(cc).params.Cvec = Cvec;
            SVM_sess(cc).params.bin_width = bin_width;
            SVM_sess(cc).params.bin_step = bin_step;
            % Save mat

        end
    end
end
fprintf('Saving...\n')

save(strcat(save_mat,filesep,'SVM_sess_',reg,'.mat'),"SVM_sess");

fprintf('Add done!!\n')




