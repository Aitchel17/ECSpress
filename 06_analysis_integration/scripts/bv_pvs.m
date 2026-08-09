clc, clear
addpath(genpath('g:\03_program\01_ecspress\06_analysis_integration\01_tablemanager'));
masterDirTable_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl\00_igkl_dirtable.xlsx';
mtable_FWHMsleep = tableManager.load_recon(masterDirTable_path, "mtable_FWHMsleep.mat");
%% State summary analysis

mtable_FWHMsleep.analysis_table = mtable_FWHMsleep.subTables.state_summary; % State summary to be target table
mtable_FWHMsleep.filtLogics = [];
mtable_FWHMsleep.filtLogics.vType = contains(mtable_FWHMsleep.analysis_table.VesselID, "PA", 'IgnoreCase', true); % Artery filter
mtable_FWHMsleep.filtLogics.depth = mtable_FWHMsleep.analysis_table.NumericDepth <70; % L1 filter
mtable_FWHMsleep.filtLogics.bout_dur = mtable_FWHMsleep.analysis_table.bout_duration > 30; % 30 seconds sleep filter
mtable_FWHMsleep.filtLogics.state = ismember(mtable_FWHMsleep.analysis_table.state_name,["awake","drowsy","nrem","rem"] );
mtable_FWHMsleep.apply_filter

%%

bv_pvsAnalyzer = tableAnalyzer(mtable_FWHMsleep.filtered_table, mtable_FWHMsleep.action_log);
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
bv_pvsAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames,"division");
%%
tmp.id_logic = {"MouseID","Date","VesselID","bout_idx"};
tmp.correlation = {"thickness_bv","thickness_totalpvs"};
tmp.state = "rem";

%%
%% Collect raw_data per (MouseID, Date, VesselID, state, bout_idx)
T = bv_pvsAnalyzer.filtered_table;

T_bv  = T(strcmp(T.DataType, 'thickness_bv'), :);
T_pvs = T(strcmp(T.DataType, 'thickness_totalpvs'), :);

results = table();
for i = 1:height(T_bv)
    row_bv = T_bv(i,:);
    
    % Find matching pvs row with same identity and bout_idx
    mask_pvs = strcmp(T_pvs.MouseID,     row_bv.MouseID) & ...
               strcmp(T_pvs.Date,        row_bv.Date) & ...
               strcmp(T_pvs.VesselID,    row_bv.VesselID) & ...
               strcmp(T_pvs.state_name,  row_bv.state_name) & ...
               T_pvs.bout_idx ==         row_bv.bout_idx;

    if ~any(mask_pvs), continue; end

    % Extract data as 1x1 cells
    bv_data = row_bv.raw_data;
    if height(T_pvs(mask_pvs,:)) > 1
        % Take the first match in case of duplicates
        all_matches = T_pvs.raw_data(mask_pvs);
        pvs_data = all_matches(1); 
    else
        pvs_data = T_pvs.raw_data(mask_pvs);
    end

    row = table(row_bv.MouseID, row_bv.Date, row_bv.VesselID, ...
                row_bv.state_name, row_bv.bout_idx, ...
                bv_data, pvs_data, ...
                'VariableNames', {'MouseID','Date','VesselID','state_name','bout_idx','raw_bv','raw_totalpvs'});
    results = [results; row];
end

%% Plot all bouts colored by state (filtered by vessel and date)
vessel_to_plot = 'PA01';
date_to_plot = '251016';

figure('Position', [2048, 431, 286, 385])
hold on;
plot_states = unique(results.state_name);
colors = lines(numel(plot_states));

for s = 1:numel(plot_states)
    st = plot_states(s);
    
    % Assign color based on state
    switch lower(char(st))
        case 'awake'
            % Gray
            c = clee.lch(55, 0, 0);
        case 'rem'
            % Scarlet
            c = clee.lch(55, 80, 10);
        case 'nrem'
            % Skyblue
            c = clee.lch(55, 80, 240);
        case 'drowsy'
            c = clee.lch(75, 80, 135);
        otherwise
            c = clee.lch(55, 0, 0); % fallback
    end
    
    % Get all bouts for this state AND target vessel AND target date
    t_state = results(strcmp(results.state_name, st) & ...
                      strcmp(results.VesselID, vessel_to_plot) & ...
                      strcmp(results.Date, date_to_plot), :);
    
    if isempty(t_state) || height(t_state) == 0
        continue;
    end
    
    % Flatten and combine all bouts into single vectors
    all_bv  = cell2mat(cellfun(@(x) x(:), t_state.raw_bv, 'UniformOutput', false));
    all_pvs = cell2mat(cellfun(@(x) x(:), t_state.raw_totalpvs, 'UniformOutput', false));
    
    % Reduce marker alpha to make scatter less busy
    scatter(all_bv, all_pvs, 5, c, 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', char(st))
    
    % Add robust fitted line
    valid_idx = ~isnan(all_bv) & ~isnan(all_pvs) & ~isinf(all_bv) & ~isinf(all_pvs);
    if sum(valid_idx) > 2
        b = robustfit(all_bv(valid_idx), all_pvs(valid_idx));
        x_fit = [min(all_bv(valid_idx)), max(all_bv(valid_idx))];
        y_fit = b(1) + b(2)*x_fit;
        plot(x_fit, y_fit, 'Color', c, 'LineWidth', 2, 'HandleVisibility', 'off');
    end
end

xlabel('BV Thickness')
ylabel('Total PVS Thickness')
%title(sprintf('All States: BV vs Total PVS (Vessel: %s | Date: %s)', vessel_to_plot, date_to_plot))
grid off;
hold off;







