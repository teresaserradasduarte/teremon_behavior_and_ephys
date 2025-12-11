%% MI - part2
%% Do Mean and PCA
mean_MI = mean(MI_landscaps_all,3,'omitnan');

% pca
n_components = 3;
[MI_pc_score,explained_MI] = pca2d_nanMat(MI_landscaps_all,n_components);

% Save append
save(strcat(save_mat,filesep,'MI_',reg,'.mat'),'x_im','y_im',...
    'mean_MI','MI_pc_score','explained_MI','-append');

%% Figures
% y_im = start_stop_neu(:,end,1);
% x_im = start_stop_neu(1,:,2);

figure
imagesc(x_im,y_im,mean_MI); hold on
axis xy;
title('Mean MI');
c=colorbar;
ylabel(c,'mean MI');
xline(0,'--','color',[1 1 1],'LineWidth',1.5)
yline(0,'--','color',[1 1 1],'LineWidth',1.5)
ylabel('start time (s)'); xlabel('stop time (s)')
axis square
title('Mean MI across neurons')

figure
ff = tiledlayout(1,3);
for pc=1:n_components
    nexttile
    imagesc(x_im,y_im,MI_pc_score(:,:,pc)); hold on
    axis xy;
    title(sprintf('%s%i%s%.2f%s','PC',pc,' (',explained_MI(pc),'%)'));
    c=colorbar;
    ylabel(c,'pc score');
    xline(0,'--','color',[1 1 1],'LineWidth',1.5)
    yline(0,'--','color',[1 1 1],'LineWidth',1.5)
    ylabel('start time (s)'); xlabel('stop time (s)')
end


%% LEAD LAG
% whithout zero
MI_landscape_size = size(MI_landscaps_all,1);
lead_win_start = 1:MI_landscape_size/2;
lead_win_stop = 1:MI_landscape_size/2;
lag_win_start = (MI_landscape_size/2)+1:MI_landscape_size;
lag_win_stop = (MI_landscape_size/2)+1:MI_landscape_size;
large_win_start = 1:MI_landscape_size/2;
large_win_stop =(MI_landscape_size/2)+1:MI_landscape_size;
MI_lead = permute(squeeze(MI_landscaps_all(lead_win_start,lead_win_stop,:)),[3 1 2]);
MI_lag = permute(squeeze(MI_landscaps_all(lag_win_start,lag_win_stop,:)),[3 1 2]);
diff_lead_lag = mean(MI_lead(:,:),2,'omitnan') -  mean(MI_lag(:,:),2,'omitnan');

% with zero
n_neu = size(MI_landscaps_all,3);
[X, Y] = meshgrid(x_im, y_im); 
idx_lead = X + Y < 0;
idx_diag = -X==Y;
tol = 1e-10;

mask_lead = zeros(size(mean_MI));
mask_lag = zeros(size(mean_MI));
mask_lead(idx_lead) = true;
mask_lag(~idx_lead) = true;
mask_lag(abs(X + Y) < tol) = false;
[row_idx_lead, col_idx_lead] = find(mask_lead);
[row_idx_lag, col_idx_lag] = find(mask_lag);

MI_lead_L = nan(n_neu,length(x_im),length(y_im));
MI_lag_L = nan(n_neu,length(x_im),length(y_im));
for nn=1:n_neu
    MI_tmp_lead = MI_landscaps_all(:,:,nn);
    MI_tmp_lag = MI_landscaps_all(:,:,nn);
    MI_tmp_lead(mask_lead==0) = nan;
    MI_tmp_lag(mask_lag==0) = nan;
    MI_lead_L(nn,:,:) = MI_tmp_lead;
    MI_lag_L(nn,:,:) = MI_tmp_lag;
end
diff_lead_lag_L = mean(MI_lead_L(:,:),2,'omitnan') -  mean(MI_lag_L(:,:),2,'omitnan');

% Entropy of landscape
   n_bins_ent = 40;
MI_per = permute(MI_landscaps_all,[3 1 2]);
entropy_MIland = nan(n_neu,1);
for n = 1:n_neu
    pLarge = histcounts(MI_per(n,:),n_bins_ent,'Normalization','probability');
    entropy_MIland(n) = -sum(pLarge.*log2(pLarge),'omitnan');
end

% Save append
save(strcat(save_mat,filesep,'MI_',reg,'.mat'),...
    'diff_lead_lag','diff_lead_lag_L','entropy_MIland','-append');

%% Figure
bin_size = 0.03; 
figure
subplot(311)
histogram(diff_lead_lag,'BinWidth',bin_size)

subplot(312)
histogram(diff_lead_lag_L,'BinWidth',bin_size)

subplot(313)
histogram(entropy_MIland,'BinWidth',bin_size)

%% slope lead lag
n_bins_ll = 20;
pLead_neu = nan(n_bins_ll,n_neu);
pLag_neu = nan(n_bins_ll,n_neu);

for n = 1:n_neu
    leadd = MI_lead_L(n,:);
    lagg = MI_lag_L(n,:);
    pLead_neu(:,n) = histcounts(leadd(~isnan(leadd)),n_bins_ll,'Normalization','probability');
    pLag_neu(:,n) = histcounts(lagg(~isnan(lagg)),n_bins_ll,'Normalization','probability');
end



%%
figure
histogram(leadd(~isnan(leadd)),n_bins_ll); hold on
histogram(lagg(~isnan(lagg)),n_bins_ll); hold off
xlim([0 2])

%%
n=3;
figure
scatter(pLead_neu(:,n),pLag_neu(:,n),'filled'); hold on
f=fit(pLead_neu(:,n),pLag_neu(:,n),'poly1');

x_vals = xlim;
y_vals = f.p1 * x_vals + f.p2;



plot(xlim, xlim, '--k')
xlabel('lead'); ylabel('lag');
plot(xlim,y_vals,'-r')

