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
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
mtable_FWHMsleep.apply_resolution("NumericResolution",data_colnames,numeric_colnames);
mtable_FWHMsleep.meanFrom2("raw_data","Q2Q3_mean",0.25,0.75)
mtable_FWHMsleep.addPrctilecol("raw_data","prctile_95", 95);
mtable_FWHMsleep.addPrctilecol("raw_data","prctile_5",5);
mtable_FWHMsleep.numeric_tables.filtered_table = mtable_FWHMsleep.filtered_table;
mtable_FWHMsleep.get_numericsummary("bout_idx","filtered_table") % intra session
mtable_FWHMsleep.get_numericsummary("Date","bout_idx_ave") % inter session
mtable_FWHMsleep.get_numericsummary("VesselID","Date_ave") % Vessel
mtable_FWHMsleep.get_numericsummary("MouseID","VesselID_ave") % mouse

result.state_summary.abs_numeric = mtable_FWHMsleep.numeric_tables;


%%
tmp.awakenumtable = result.state_summary.abs_numeric.bout_idx_ave(result.state_summary.abs_numeric.bout_idx_ave.state_name=="awake",:);

% 1. Extract the denominator baseline
tmp.key_cols = ["MouseID","Date","VesselID","DataType"];
tmp.awake_baseline = tmp.awakenumtable(:, [tmp.key_cols,"Q2Q3_mean"]);

summarystat_awakenorm = outerjoin(result.state_summary.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);







%% transition analysis
mtable_FWHMsleep.analysis_table = mtable_FWHMsleep.subTables.transition; % State summary to transition table 
mtable_FWHMsleep.filtLogics = [];
mtable_FWHMsleep.filtLogics.vType = contains(mtable_FWHMsleep.analysis_table.VesselID, "PA", 'IgnoreCase', true); % Artery filter
mtable_FWHMsleep.filtLogics.depth = mtable_FWHMsleep.analysis_table.NumericDepth <70; % L1 filter
mtable_FWHMsleep.apply_filter
data_colnames = {"data"};
numeric_colnames = {'pre_mean','pre_median','pre_q1','pre_q3', 'pre_var',...
                'post_mean','post_median','post_q1','post_q3', 'post_var'};
mtable_FWHMsleep.apply_resolution("NumericResolution",data_colnames,numeric_colnames);
mtable_FWHMsleep.numeric_tables.filtered_table = mtable_FWHMsleep.filtered_table;
mtable_FWHMsleep.get_numericsummary("bout_idx","filtered_table") % intra session bout average
mtable_FWHMsleep.get_numericsummary("Date","bout_idx_ave") % inter session average
mtable_FWHMsleep.get_numericsummary("VesselID","Date_ave") % within mouse vessel average
mtable_FWHMsleep.get_numericsummary("MouseID","VesselID_ave") % mouse average
mtable_FWHMsleep.get_datasummary("bout_idx","filtered_table")
mtable_FWHMsleep.get_datasummary("Date","bout_idx_ave")
mtable_FWHMsleep.get_datasummary("VesselID","Date_ave")
mtable_FWHMsleep.get_datasummary("MouseID","VesselID_ave")
result.transition.abs_numeric = mtable_FWHMsleep.numeric_tables;
result.transition.abs_data = mtable_FWHMsleep.data_tables;
%%

%% Normalize transition data by awake Q2Q3_mean


% 2. Normalize the Transition Table 
% We outerjoin keeping the filtered table keys.
trans_norm = outerjoin(trans.abs_numeric.filtered_table, awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', key_cols, ...
    'RightKeys', key_cols);

% Calculate relative changes for NUMERIC columns
numeric_colnames = {'pre_mean','pre_median','pre_q1','pre_q3', 'pre_var', ...
                    'post_mean','post_median','post_q1','post_q3', 'post_var'};
for num_idx = 1:numel(numeric_colnames)
    col_name = numeric_colnames{num_idx};
    trans_norm.("rel_" + col_name) = trans_norm.(col_name) ./ trans_norm.awake_baseline_Q2Q3;
end

% Calculate relative changes for DATA array columns
% The numeric 'filtered_table' contains the 'data' cell column as well!
data_colnames = {"data"};
for dat_idx = 1:numel(data_colnames)
    col_name = data_colnames{dat_idx};
    rel_data_cell = cell(height(trans_norm), 1);
    
    for row = 1:height(trans_norm)
        baseline_val = trans_norm.awake_baseline_Q2Q3(row);
        arr_val = trans_norm.(col_name){row};
        
        if ~isnan(baseline_val) && baseline_val ~= 0 && ~isempty(arr_val)
            rel_data_cell{row} = arr_val ./ baseline_val;
        else
            rel_data_cell{row} = arr_val; % Keep unnormalized or you could assign NaN
        end
    end
    trans_norm.("rel_" + col_name) = rel_data_cell;
end
