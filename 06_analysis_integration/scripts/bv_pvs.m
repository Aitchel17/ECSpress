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
%% Collect raw_data per (MouseID, Date, VesselID, state) for thickness_bv and thickness_totalpvs
T = bv_pvsAnalyzer.filtered_table;

% Separate by DataType
T_bv  = T(strcmp(T.DataType, 'thickness_bv'), :);
T_pvs = T(strcmp(T.DataType, 'thickness_totalpvs'), :);

states = unique(T.state_name);
id_vars = {'MouseID', 'Date', 'VesselID'};

% Get unique vessel identities
[unique_ids, ~, ic_bv] = unique(T_bv(:, id_vars), 'rows', 'stable');

results = table();
for i = 1:height(unique_ids)
    for s = 1:numel(states)
        state = states(s);

        % Match rows in T_bv
        mask_bv = ic_bv == i & strcmp(T_bv.state_name, state);
        % Match rows in T_pvs for same identity
        mask_pvs = strcmp(T_pvs.MouseID, unique_ids.MouseID(i)) & ...
                   strcmp(T_pvs.Date, unique_ids.Date(i)) & ...
                   strcmp(T_pvs.VesselID, unique_ids.VesselID(i)) & ...
                   strcmp(T_pvs.state_name, state);

        if ~any(mask_bv) || ~any(mask_pvs), continue; end

        % Concatenate all raw_data cells into a single vector
    raw_bv  = T_bv.raw_data(mask_bv);    % cell array, one entry per bout
    raw_pvs = T_pvs.raw_data(mask_pvs);  % cell array, one entry per bout

        row = [unique_ids(i,:), table(state, {raw_bv}, {raw_pvs}, ...
               'VariableNames', {'state_name','raw_bv','raw_totalpvs'})];
        results = [results; row];
    end
end

