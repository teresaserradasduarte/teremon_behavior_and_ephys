function plot_reach_raster_and_psth(eg_neurons,i,bin_edges,win_interest,figProp,trialsIdx)

      tt = tiledlayout(3,1);
    title(tt,sprintf('%s%s%s%s%s%s%s%s%s',...
        figProp.name,...
        ' | region: ',eg_neurons(i).reg,...
        ' | phyID: ', num2str(eg_neurons(i).phyID),...
        ' | idx: ', num2str(eg_neurons(i).unitIdx),...
        ' | idxEG: ',num2str(eg_neurons(i).egIdx)),...
        'Interpreter', 'none','fontsize',12,'fontweight','normal')


    % Get spike times and corresponding trial indices
    spike_times = eg_neurons(i).st_reach(trialsIdx)';
    trials_idx_scat = cell2mat(arrayfun(@(i) i * ones(1, numel(spike_times{i})),...
    1:length(spike_times), 'UniformOutput', false));

    nexttile;
    plot(bin_edges,mean(eg_neurons(i).FR_reach(:,trialsIdx),2,'omitnan'),'Color',figProp.rcl_clr,LineWidth=2)
    hold off
    xlim(win_interest)
    xlabel('time from reack peak (s)'); 
    ylabel('firing rate sp/s')
    set(gca,figProp.axeOpt{:})

    axx2 = nexttile;
    axx2.Layout.TileSpan = [2,1];
    scatter(cell2mat(spike_times),trials_idx_scat,figProp.spk_size,'k.'), hold on
    %scatter(zeros(size(reaches_vec)),reaches_vec,event_size,'filled','MarkerFaceColor',rcl_clr)
    hold off
    axis tight
    xlim(win_interest)
    xlabel('time from reack peak (s)');
    set(gca,figProp.axeOpt{:})
    xlim(win_interest)
    set(gcf,'position',[2621 207 539 697],'Color','w')

    yticks_to_label = round(linspace(1, length(trialsIdx), 10));
    ytick_labels = figProp.ticks_y(yticks_to_label);
    ylim([1 length(trialsIdx)])
    yticks(yticks_to_label);
    yticklabels(string(ytick_labels));
    ylabel(figProp.y_name)
    shg