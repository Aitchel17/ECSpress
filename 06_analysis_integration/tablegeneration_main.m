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
%% Absolute analysis
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
summarystatAnalyzer = tableAnalyzer(mtable_FWHMsleep.filtered_table, mtable_FWHMsleep.action_log);
summarystatAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames,"division");
summarystatAnalyzer = common_summarystat_analysis(summarystatAnalyzer);
result.state_summary.abs_numeric = summarystatAnalyzer.numeric_tables;


% Awake normalization/subtraction
tmp.awakenumtable = result.state_summary.abs_numeric.bout_idx_ave(result.state_summary.abs_numeric.bout_idx_ave.state_name=="awake",:);

% 0. Extract the denominator baseline
tmp.key_cols = ["MouseID","Date","VesselID","DataType"];
tmp.awake_baseline = tmp.awakenumtable(:, [tmp.key_cols,"Q2Q3_mean"]);
tmp.awake_baseline.Properties.VariableNames("Q2Q3_mean") = "awake_boutq2q3mean";
tmp.awake_baseline.awake_scale = 1./tmp.awake_baseline.awake_boutq2q3mean;

summarystat_awakenorm = outerjoin(result.state_summary.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
% 1. get awake normalized table
awakescaled_summarystatAnalyzer = tableAnalyzer(summarystat_awakenorm, mtable_FWHMsleep.action_log);
awakescaled_summarystatAnalyzer.scale_table("awake_scale",data_colnames,numeric_colnames,"division");
awakescaled_summarystatAnalyzer = common_summarystat_analysis(awakescaled_summarystatAnalyzer);
result.state_summary.awakenorm_numeric = awakescaled_summarystatAnalyzer.numeric_tables;

% 2. get awake subtracted table
awakesubtract_summarystatAnalyzer = tableAnalyzer(summarystat_awakenorm, mtable_FWHMsleep.action_log);
awakesubtract_summarystatAnalyzer.scale_table("awake_boutq2q3mean",data_colnames,numeric_colnames,"subtraction");
awakesubtract_summarystatAnalyzer = common_summarystat_analysis(awakesubtract_summarystatAnalyzer);
result.state_summary.awakesubtract_numeric = awakesubtract_summarystatAnalyzer.numeric_tables;



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
transitionAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames,"division");
transitionAnalyzer = common_transition_analysis(transitionAnalyzer);
result.transition.abs_numeric = transitionAnalyzer.numeric_tables;
result.transition.abs_data = transitionAnalyzer.data_tables;

% Normalize transition data by awake Q2Q3_mean
transition_awakenorm = outerjoin(result.transition.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
% Apply normalization using awake_scale
awakescaled_transitionAnalyzer = tableAnalyzer(transition_awakenorm, mtable_FWHMsleep.action_log);
awakescaled_transitionAnalyzer.scale_table("awake_scale", data_colnames, numeric_colnames,"division");
% Re-calculate numeric and data summaries
awakescaled_transitionAnalyzer = common_transition_analysis(awakescaled_transitionAnalyzer);
result.transition.awakenorm_numeric = awakescaled_transitionAnalyzer.numeric_tables;
result.transition.awakenorm_data = awakescaled_transitionAnalyzer.data_tables;
%%
clc
awakesubtracted_transitionAnalyzer = tableAnalyzer(transition_awakenorm, mtable_FWHMsleep.action_log);
awakesubtracted_transitionAnalyzer.scale_table("awake_boutq2q3mean", data_colnames, numeric_colnames,"subtraction");
awakesubtracted_transitionAnalyzer = common_transition_analysis(awakesubtracted_transitionAnalyzer);
result.transition.awakesubtract_numeric = awakesubtracted_transitionAnalyzer.numeric_tables;
result.transition.awakesubtract_data = awakesubtracted_transitionAnalyzer.data_tables;







%% saving mechanism

tmp.savefieldnames = fieldnames(result);
tmp.savefilenames = strcat(string(tmp.savefieldnames),'.mat');
tmp.savepaths =  fullfile(mtable_FWHMsleep.save_dir , string(tmp.savefilenames));

for save_idx = 1:numel(tmp.savefieldnames)
    save_content = result.(tmp.savefieldnames{save_idx});
    save(tmp.savepaths{save_idx}, "save_content")
end


function tableAnalyzer = common_summarystat_analysis(tableAnalyzer)
    tableAnalyzer.meanFrom2("raw_data","Q2Q3_mean",0.25,0.75)
    tableAnalyzer.addPrctilecol("raw_data","prctile_95", 95);
    tableAnalyzer.addPrctilecol("raw_data","prctile_5",5);
    tableAnalyzer.numeric_tables.filtered_table = tableAnalyzer.filtered_table;
    tableAnalyzer.get_numericsummary("bout_idx","filtered_table") % intra session
    tableAnalyzer.get_numericsummary("Date","bout_idx_ave") % inter session
    tableAnalyzer.get_numericsummary("VesselID","Date_ave") % Vessel
    tableAnalyzer.get_numericsummary("MouseID","VesselID_ave") % mouse
end



function tableAnalyzer = common_transition_analysis(tableAnalyzer)
    tableAnalyzer.numeric_tables.filtered_table = tableAnalyzer.filtered_table;
    tableAnalyzer.get_numericsummary("bout_idx","filtered_table") % intra session bout average
    tableAnalyzer.get_numericsummary("Date","bout_idx_ave") % inter session average
    tableAnalyzer.get_numericsummary("VesselID","Date_ave") % within mouse vessel average
    tableAnalyzer.get_numericsummary("MouseID","VesselID_ave") % mouse average
    tableAnalyzer.get_datasummary("bout_idx","filtered_table")
    tableAnalyzer.get_datasummary("Date","bout_idx_ave")
    tableAnalyzer.get_datasummary("VesselID","Date_ave")
    tableAnalyzer.get_datasummary("MouseID","VesselID_ave")
end









