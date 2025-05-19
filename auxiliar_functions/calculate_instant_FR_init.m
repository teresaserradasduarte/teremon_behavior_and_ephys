    
function [eg_neurons,FR_init_raw] = calculate_instant_FR_init(eg_neurons,i,event_tm,win_padded,bin_width)

% Dimensions of FR matrix
bin_edges_padded = win_padded(1):bin_width:win_padded(2);
n_bins_padded = length(bin_edges_padded);
n_events = length(event_tm);
% initialize FR
FR_init_raw = zeros(n_bins_padded,n_events);

for j = 1:n_events
    spike_init_flags = ...
        eg_neurons(i).st > event_tm(j) + win_padded(1) & ...
        eg_neurons(i).st <= event_tm(j) + win_padded(2);

    eg_neurons(i).st_init{j} = ...
        eg_neurons(i).st(spike_init_flags) - event_tm(j);
    eg_neurons(i).st_initIdx{j} = ones(size(eg_neurons(i).st_init{j}))*j;
    init_spikes = eg_neurons(i).st_init{j};

    if ~isempty(init_spikes)
        % Find the spike before interval fo initial ISI
        first_spk_idx=find(spike_init_flags,1,'first');
        if first_spk_idx == 1
            FR_init_raw(bin_edges_padded<=init_spikes(1),j) = 0;
        else
            start_instant_isi = eg_neurons(i).st(first_spk_idx) - ...
                eg_neurons(i).st(first_spk_idx-1);
            FR_init_raw(bin_edges_padded<=init_spikes(1),j) = 1/start_instant_isi;
        end
        % Find spikes within the interval
        for sp=2:length(init_spikes)
            with_isi_flag = (bin_edges_padded > init_spikes(sp-1)) & ...
                bin_edges_padded <= init_spikes(sp);
            inst_isi = init_spikes(sp)-init_spikes(sp-1);
            FR_init_raw(with_isi_flag,j) = 1/inst_isi;
        end

        % find the spike after the interval for final ISI
        last_spk_idx=find(spike_init_flags,1,'last');
        if last_spk_idx==length(eg_neurons(i).st)
            FR_init_raw(bin_edges_padded>init_spikes(end),j) = 0;
        else
            finish_instant_isi = eg_neurons(i).st(last_spk_idx+1) - ...
                eg_neurons(i).st(last_spk_idx);
            FR_init_raw(bin_edges_padded>init_spikes(end),j) = 1/finish_instant_isi;
        end
    end
end
end