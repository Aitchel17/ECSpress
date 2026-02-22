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
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
myAnalyzer = TableAnalyzer(mtable_FWHMsleep.filtered_table, mtable_FWHMsleep.action_log);
myAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames);
myAnalyzer.meanFrom2("raw_data","Q2Q3_mean",0.25,0.75)
myAnalyzer.addPrctilecol("raw_data","prctile_95", 95);
myAnalyzer.addPrctilecol("raw_data","prctile_5",5);
myAnalyzer.numeric_tables.filtered_table = myAnalyzer.filtered_table;
myAnalyzer.get_numericsummary("bout_idx","filtered_table") % intra session
myAnalyzer.get_numericsummary("Date","bout_idx_ave") % inter session
myAnalyzer.get_numericsummary("VesselID","Date_ave") % Vessel
myAnalyzer.get_numericsummary("MouseID","VesselID_ave") % mouse

result.state_summary.abs_numeric = myAnalyzer.numeric_tables;


%%
tmp.awakenumtable = result.state_summary.abs_numeric.bout_idx_ave(result.state_summary.abs_numeric.bout_idx_ave.state_name=="awake",:);

% 1. Extract the denominator baseline
tmp.key_cols = ["MouseID","Date","VesselID","DataType"];
tmp.awake_baseline = tmp.awakenumtable(:, [tmp.key_cols,"Q2Q3_mean"]);
tmp.awake_baseline.Properties.VariableNames("Q2Q3_mean") = "awake_boutq2q3mean";
tmp.awake_baseline.awake_scale = 1./tmp.awake_baseline.awake_boutq2q3mean;

summarystat_awakenorm = outerjoin(result.state_summary.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);

myAnalyzer2 = TableAnalyzer(summarystat_awakenorm, mtable_FWHMsleep.action_log);
myAnalyzer2.scale_table("awake_scale",data_colnames,numeric_colnames);
myAnalyzer2.meanFrom2("raw_data","Q2Q3_mean",0.25,0.75)
myAnalyzer2.addPrctilecol("raw_data","prctile_95", 95);
myAnalyzer2.addPrctilecol("raw_data","prctile_5",5);
myAnalyzer2.numeric_tables.filtered_table = myAnalyzer2.filtered_table;
myAnalyzer2.get_numericsummary("bout_idx","filtered_table") % intra session
myAnalyzer2.get_numericsummary("Date","bout_idx_ave") % inter session
myAnalyzer2.get_numericsummary("VesselID","Date_ave") % Vessel
myAnalyzer2.get_numericsummary("MouseID","VesselID_ave") % mouse




%% transition analysis
mtable_FWHMsleep.analysis_table = mtable_FWHMsleep.subTables.transition; % State summary to transition table 
mtable_FWHMsleep.filtLogics = [];
mtable_FWHMsleep.filtLogics.vType = contains(mtable_FWHMsleep.analysis_table.VesselID, "PA", 'IgnoreCase', true); % Artery filter
mtable_FWHMsleep.filtLogics.depth = mtable_FWHMsleep.analysis_table.NumericDepth <70; % L1 filter
mtable_FWHMsleep.apply_filter
data_colnames = {"data"};
numeric_colnames = {'pre_mean','pre_median','pre_q1','pre_q3', 'pre_var',...
                'post_mean','post_median','post_q1','post_q3', 'post_var'};
myTransAnalyzer = TableAnalyzer(mtable_FWHMsleep.filtered_table, mtable_FWHMsleep.action_log);
myTransAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames);
myTransAnalyzer.numeric_tables.filtered_table = myTransAnalyzer.filtered_table;
myTransAnalyzer.get_numericsummary("bout_idx","filtered_table") % intra session bout average
myTransAnalyzer.get_numericsummary("Date","bout_idx_ave") % inter session average
myTransAnalyzer.get_numericsummary("VesselID","Date_ave") % within mouse vessel average
myTransAnalyzer.get_numericsummary("MouseID","VesselID_ave") % mouse average
myTransAnalyzer.get_datasummary("bout_idx","filtered_table")
myTransAnalyzer.get_datasummary("Date","bout_idx_ave")
myTransAnalyzer.get_datasummary("VesselID","Date_ave")
myTransAnalyzer.get_datasummary("MouseID","VesselID_ave")
result.transition.abs_numeric = myTransAnalyzer.numeric_tables;
result.transition.abs_data = myTransAnalyzer.data_tables;
%%

%% Normalize transition data by awake Q2Q3_mean

transition_awakenorm = outerjoin(result.transition.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);


%% Apply normalization using awake_scale
myTransAnalyzer2 = TableAnalyzer(transition_awakenorm, mtable_FWHMsleep.action_log);
myTransAnalyzer2.scale_table("awake_scale", data_colnames, numeric_colnames);

% Re-calculate numeric and data summaries
myTransAnalyzer2.numeric_tables.filtered_table = myTransAnalyzer2.filtered_table;
myTransAnalyzer2.get_numericsummary("bout_idx", "filtered_table") % intra session bout average
myTransAnalyzer2.get_numericsummary("Date", "bout_idx_ave") % inter session average
myTransAnalyzer2.get_numericsummary("VesselID", "Date_ave") % within mouse vessel average
myTransAnalyzer2.get_numericsummary("MouseID", "VesselID_ave") % mouse average

myTransAnalyzer2.get_datasummary("bout_idx", "filtered_table")
myTransAnalyzer2.get_datasummary("Date", "bout_idx_ave")
myTransAnalyzer2.get_datasummary("VesselID", "Date_ave")
myTransAnalyzer2.get_datasummary("MouseID", "VesselID_ave")

result.transition_norm.abs_numeric = myTransAnalyzer2.numeric_tables;
result.transition_norm.abs_data = myTransAnalyzer2.data_tables;
