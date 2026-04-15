function plot_cell_metrics(neu_table, color_var, out_path, region_name)
% PLOT_CELL_METRICS Plots waveforms and CCGs color-coded by a metric and saves the figure.
%   Args:
%       neu_table: Table containing cell metrics and waveforms
%       color_var: String, e.g., 'pred_probability', 'confidence_ratio', 'model_votes'
%       out_path:  String, folder path where the image will be saved
%       region_name: String, name of the region (e.g., 'CB') for the filename

    %% 1. Data Preparation & Color Variable Setup
    
    % Handle 'model_votes' conversion
    if contains(color_var, 'votes')
        votes_str = neu_table.model_votes;
        votes_num = zeros(height(neu_table), 1);
        for i = 1:height(neu_table)
            if isempty(votes_str{i})
                votes_num(i) = NaN; 
            else
                parts = split(votes_str{i}, '/');
                votes_num(i) = str2double(parts{1}) / str2double(parts{2});
            end
        end
        neu_table.votes_numeric = votes_num; 
        target_col = 'votes_numeric'; 
    else
        target_col = color_var;
    end

    % --- CUSTOM LIMIT LOGIC ---
    if strcmp(color_var, 'pred_probability')
        clims = [0.2 1.0];
        
    elseif contains(color_var, 'votes')
        clims = [0.8 1.0];
        
    else
        % Default/Confidence Ratio: Dynamic (0 to 95th percentile)
        all_vals = neu_table.(target_col);
        valid_vals = all_vals(~isnan(all_vals) & ~isinf(all_vals));
        
        if isempty(valid_vals)
            clims = [0 1];
        else
            d_max = quantile(valid_vals, 0.95); 
            if d_max == 0; d_max = 1; end 
            clims = [0 d_max];
        end
    end
    % --------------------------

    % Filter data (Remove empty types and NaNs in target metric)
    is_valid = ~cellfun(@isempty, neu_table.predicted_cell_type) & ...
               ~isnan(neu_table.(target_col));
    neu_plot = neu_table(is_valid, :);

    cell_types = unique(neu_plot.predicted_cell_type);
    num_types = length(cell_types);

    % Setup Colormap
    cmap = turbo(256); 
    
    % Mapping Function
    map_val_to_idx = @(v) max(1, min(256, round( (v - clims(1)) / (clims(2) - clims(1)) * 256 )));
    get_color = @(val) cmap(map_val_to_idx(val), :);

    fs = 30000; 
    ccg_lims_default = [-50 50];
    ccg_lims_pkc_cs  = [-150 150];

    %% 2. Setup Figure
    fig_height = max(500, 200 * num_types);
    f = figure('Name', ['Overview: ' color_var], 'Color', 'w', 'Position', [2785, 50, 1020, fig_height]);
    
    tl = tiledlayout(num_types, 2, 'TileSpacing', 'tight', 'Padding', 'compact');
    title(tl, ['Waveforms and CCGs (Color = ' strrep(color_var, '_', ' ') ')'], 'Interpreter', 'none');

    %% 3. Loop through Cell Types
    for t = 1:num_types
        current_type = cell_types{t};
        
        rows = strcmp(neu_plot.predicted_cell_type, current_type);
        subset_neu = neu_plot(rows, :);
        
        % --- LEFT COLUMN: Waveforms ---
        ax1 = nexttile(tl);
        hold(ax1, 'on');
        
        for i = 1:height(subset_neu)
            wf_data = subset_neu.templateWaveforms{i};
            peak_chs = subset_neu.templatePeakCh{i};
            val = subset_neu.(target_col)(i);
            
            clr = get_color(val);
            time_wf = (0 : size(wf_data, 1) - 1) / fs * 1000;
            
            primary_peak_ch = peak_chs(1); 
            if ndims(wf_data) == 3
                 trace = wf_data(:, primary_peak_ch, 1); 
            else
                 trace = wf_data(:, primary_peak_ch);
            end
            
            plot(ax1, time_wf, trace, 'Color', [clr, 0.5], 'LineWidth', 1.2);
        end
        
        ylabel(ax1, {current_type; 'Amplitude (\muV)'}, 'FontWeight', 'bold', 'Interpreter', 'none');
        if t == 1; title(ax1, 'Waveforms'); end
        if t == num_types; xlabel(ax1, 'Time (ms)'); end
        axis(ax1, 'tight'); box(ax1, 'off');
        
        % --- RIGHT COLUMN: CCGs ---
        ax2 = nexttile(tl);
        hold(ax2, 'on');
        
        for i = 1:height(subset_neu)
            ccg_data = subset_neu.CCG{i};
            ccg_bins = subset_neu.CCG_bins{i};
            val = subset_neu.(target_col)(i);
            clr = get_color(val);
            
            time_ccg = ccg_bins * 1000; 
            plot(ax2, time_ccg, ccg_data, 'Color', [clr, 0.6], 'LineWidth', 1.5);
        end
        
        if t == 1; title(ax2, 'Auto-Correlograms'); end
        if t == num_types; xlabel(ax2, 'Lag (ms)'); end
        ylabel(ax2, 'Rate');
        
        % --- CUSTOM LIMITS FOR PkC_cs ---
        if contains(current_type, 'PkC_cs', 'IgnoreCase', true)
            xlim(ax2, ccg_lims_pkc_cs);
        else
            xlim(ax2, ccg_lims_default);
        end
        
        box(ax2, 'off');
    end

    %% 4. Final Polish (Colorbar & Save)
    colormap(f, cmap);
    
    cb = colorbar(ax2);       
    cb.Layout.Tile = 'east';  
    
    cb.Label.String = strrep(color_var, '_', ' ');
    cb.Label.Interpreter = 'none';
    cb.Label.FontWeight = 'bold';
    clim(clims);

    % --- SAVE FIGURE ---
    if ~exist(out_path, 'dir')
       mkdir(out_path);
    end
    
    % Construct filename: e.g., "CB_pred_probability.png"
    filename = sprintf('%s_%s.png', region_name, color_var);
    full_path = fullfile(out_path, filename);
    
    saveas(f, full_path);
    fprintf('Figure saved to: %s\n', full_path);
end