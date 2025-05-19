%% Check on individual windows the distribution of FR used for MI
%save_fig_name = strcat('entropy_dist',filesep,'entropy_neu');
tt = tiledlayout(1,3);
    title(ff,sprintf('%s%s%s%s%s%s%s%i','Entropy in ',reg,...
        ' | mouse: ',char(N(n,:).mouse),...
        ' | session: ',char(N(n,:).sess),...
        ' | neuron phyID: ',N(n,:).phyID),...
        'Interpreter','none');
nexttile
histogram(pawEndPoint_win_dist,hist_nbins,'Normalization','probability');
    title(sprintf('%s%.2f','Reach endpoint entropy = ',entropy_paw));
    set(gca,axeOpt{:});
    xlabel('endpoint in y (px)'); ylabel('probability');

nexttile
histogram(neu_win_dist,hist_nbins,'Normalization','probability');
    title(sprintf('%s%.2f','Neuron entropy = ',entropy_neu));
    set(gca,axeOpt{:});
    xlabel('endpoint in y (px)'); ylabel('probability');

nexttile
histogram2(pawEndPoint_win_dist,neu_win_dist,hist_nbins,'Normalization','probability');
    title(sprintf('%s%.2f','Reach endpoint entropy = ',entropy_neu));
    set(gca,axeOpt{:});
    xlabel('endpoint in y (px)'); ylabel('probability');

    set(gcf,'Position',[3036         377         560         420],'color','w');
    %saveas(gcf,strcat(save_out,filesep,'save_fig_name','.png'),'png');