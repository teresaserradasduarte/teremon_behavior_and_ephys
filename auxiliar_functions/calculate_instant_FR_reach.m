    
function [eg_neurons,FR_reach_raw] = calculate_instant_FR_reach(eg_neurons,i,event_tm,win_padded,bin_width)

% Dimensions of FR matrix
bin_edges_padded = win_padded(1):bin_width:win_padded(2);
n_bins_padded = length(bin_edges_padded);
n_reach_inVec = length(event_tm);
% initialize FR
FR_reach_raw = zeros(n_bins_padded,n_reach_inVec);

for j = 1:n_reach_inVec
    spike_reach_flags = ...
        eg_neurons(i).st > event_tm(j) + win_padded(1) & ...
        eg_neurons(i).st <= event_tm(j) + win_padded(2);

    eg_neurons(i).st_reach{j} = ...
        eg_neurons(i).st(spike_reach_flags) - event_tm(j);
    eg_neurons(i).st_reachIdx{j} = ones(size(eg_neurons(i).st_reach{j}))*j;
    reach_spikes = eg_neurons(i).st_reach{j};

    if ~isempty(reach_spikes)
        % Find the spike before interval fo initial ISI
        first_spk_idx=find(spike_reach_flags,1,'first');
        if first_spk_idx == 1
            FR_reach_raw(bin_edges_padded<=reach_spikes(1),j) = 0;
        else
            start_instant_isi = eg_neurons(i).st(first_spk_idx) - ...
                eg_neurons(i).st(first_spk_idx-1);
            FR_reach_raw(bin_edges_padded<=reach_spikes(1),j) = 1/start_instant_isi;
        end
        % Find spikes within the interval
        for sp=2:length(reach_spikes)
            with_isi_flag = (bin_edges_padded > reach_spikes(sp-1)) & ...
                bin_edges_padded <= reach_spikes(sp);
            inst_isi = reach_spikes(sp)-reach_spikes(sp-1);
            FR_reach_raw(with_isi_flag,j) = 1/inst_isi;
        end

        % find the spike after the interval for final ISI
        last_spk_idx=find(spike_reach_flags,1,'last');
        finish_instant_isi = eg_neurons(i).st(last_spk_idx+1) - ...
            eg_neurons(i).st(last_spk_idx);
        FR_reach_raw(bin_edges_padded>reach_spikes(end),j) = 1/finish_instant_isi;
    end
end
end