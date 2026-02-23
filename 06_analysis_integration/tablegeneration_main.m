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
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
summarystatAnalyzer = tableAnalyzer(mtable_FWHMsleep.filtered_table, mtable_FWHMsleep.action_log);
summarystatAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames);
summarystatAnalyzer.meanFrom2("raw_data","Q2Q3_mean",0.25,0.75)
summarystatAnalyzer.addPrctilecol("raw_data","prctile_95", 95);
summarystatAnalyzer.addPrctilecol("raw_data","prctile_5",5);
summarystatAnalyzer.numeric_tables.filtered_table = summarystatAnalyzer.filtered_table;
summarystatAnalyzer.get_numericsummary("bout_idx","filtered_table") % intra session
summarystatAnalyzer.get_numericsummary("Date","bout_idx_ave") % inter session
summarystatAnalyzer.get_numericsummary("VesselID","Date_ave") % Vessel
summarystatAnalyzer.get_numericsummary("MouseID","VesselID_ave") % mouse
result.state_summary.abs_numeric = summarystatAnalyzer.numeric_tables;


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

awakescaled_summarystatAnalyser = tableAnalyzer(summarystat_awakenorm, mtable_FWHMsleep.action_log);
awakescaled_summarystatAnalyser.scale_table("awake_scale",data_colnames,numeric_colnames);
awakescaled_summarystatAnalyser.meanFrom2("raw_data","Q2Q3_mean",0.25,0.75)
awakescaled_summarystatAnalyser.addPrctilecol("raw_data","prctile_95", 95);
awakescaled_summarystatAnalyser.addPrctilecol("raw_data","prctile_5",5);
awakescaled_summarystatAnalyser.numeric_tables.filtered_table = awakescaled_summarystatAnalyser.filtered_table;
awakescaled_summarystatAnalyser.get_numericsummary("bout_idx","filtered_table") % intra session
awakescaled_summarystatAnalyser.get_numericsummary("Date","bout_idx_ave") % inter session
awakescaled_summarystatAnalyser.get_numericsummary("VesselID","Date_ave") % Vessel
awakescaled_summarystatAnalyser.get_numericsummary("MouseID","VesselID_ave") % mouse

result.state_summary.awkaenorm_numeric = awakescaled_summarystatAnalyser.numeric_tables;


%% transition analysis
mtable_FWHMsleep.analysis_table = mtable_FWHMsleep.subTables.transition; % State summary to transition table 
mtable_FWHMsleep.filtLogics = [];
mtable_FWHMsleep.filtLogics.vType = contains(mtable_FWHMsleep.analysis_table.VesselID, "PA", 'IgnoreCase', true); % Artery filter
mtable_FWHMsleep.filtLogics.depth = mtable_FWHMsleep.analysis_table.NumericDepth <70; % L1 filter
mtable_FWHMsleep.apply_filter
data_colnames = {"data"};
numeric_colnames = {'pre_mean','pre_median','pre_q1','pre_q3', 'pre_var',...
                'post_mean','post_median','post_q1','post_q3', 'post_var'};
transitionAnalyzer = tableAnalyzer(mtable_FWHMsleep.filtered_table, mtable_FWHMsleep.action_log);
transitionAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames);
transitionAnalyzer.numeric_tables.filtered_table = transitionAnalyzer.filtered_table;
transitionAnalyzer.get_numericsummary("bout_idx","filtered_table") % intra session bout average
transitionAnalyzer.get_numericsummary("Date","bout_idx_ave") % inter session average
transitionAnalyzer.get_numericsummary("VesselID","Date_ave") % within mouse vessel average
transitionAnalyzer.get_numericsummary("MouseID","VesselID_ave") % mouse average
transitionAnalyzer.get_datasummary("bout_idx","filtered_table")
transitionAnalyzer.get_datasummary("Date","bout_idx_ave")
transitionAnalyzer.get_datasummary("VesselID","Date_ave")
transitionAnalyzer.get_datasummary("MouseID","VesselID_ave")
result.transition.abs_numeric = transitionAnalyzer.numeric_tables;
result.transition.abs_data = transitionAnalyzer.data_tables;
%%

%% Normalize transition data by awake Q2Q3_mean

transition_awakenorm = outerjoin(result.transition.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);


% Apply normalization using awake_scale
awakescaled_transitionAnalyzer = tableAnalyzer(transition_awakenorm, mtable_FWHMsleep.action_log);
awakescaled_transitionAnalyzer.scale_table("awake_scale", data_colnames, numeric_colnames);

% Re-calculate numeric and data summaries
awakescaled_transitionAnalyzer.numeric_tables.filtered_table = awakescaled_transitionAnalyzer.filtered_table;
awakescaled_transitionAnalyzer.get_numericsummary("bout_idx", "filtered_table") % intra session bout average
awakescaled_transitionAnalyzer.get_numericsummary("Date", "bout_idx_ave") % inter session average
awakescaled_transitionAnalyzer.get_numericsummary("VesselID", "Date_ave") % within mouse vessel average
awakescaled_transitionAnalyzer.get_numericsummary("MouseID", "VesselID_ave") % mouse average

awakescaled_transitionAnalyzer.get_datasummary("bout_idx", "filtered_table")
awakescaled_transitionAnalyzer.get_datasummary("Date", "bout_idx_ave")
awakescaled_transitionAnalyzer.get_datasummary("VesselID", "Date_ave")
awakescaled_transitionAnalyzer.get_datasummary("MouseID", "VesselID_ave")

result.transition.awakenorm_numeric = awakescaled_transitionAnalyzer.numeric_tables;
result.transition.awakenorm_data = awakescaled_transitionAnalyzer.data_tables;

%% saving mechanism

tmp.savefieldnames = fieldnames(result);
tmp.savefilenames = strcat(string(tmp.savefieldnames),'.mat');
tmp.savepaths =  fullfile(mtable_FWHMsleep.save_dir , string(tmp.savefilenames));

for save_idx = 1:numel(tmp.savefieldnames)
    save_content = result.(tmp.savefieldnames{save_idx});
    save(tmp.savepaths{save_idx}, "save_content")
end















