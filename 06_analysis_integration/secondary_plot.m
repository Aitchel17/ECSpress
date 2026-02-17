clc, clear
masterDirTable_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl\00_igkl_dirtable.xlsx';
mtable_FWHMsleep = tableManager.load_recon(masterDirTable_path, "mtable_FWHMsleep.mat");
%%
transition_table = mtable_FWHMsleep.subTables.transition;
%% logic generation
% translogic = transition_table.state_name =="nr_trans";
vessellogic = contains(transition_table.VesselID, "PA", 'IgnoreCase', true);
depthlogic = transition_table.Depth <70;
combinedlogic = all([vessellogic, depthlogic],2);
%%
L1pa_transtable =   transition_table(combinedlogic,:);
data_col = "data";
numeric_cols = {'pre_mean','pre_median','pre_q1','pre_q3', 'post_mean', 'post_median', 'post_q1', 'post_q3'};
key_cols = ['MouseID', 'VesselID', "state_name"];
sessAve_transitionTable = average_table(L1pa_transtable,key_cols, numeric_cols,data_col);
%%
summary_table = mtable_FWHMsleep.subTables.state_summary;
vessellogic = contains(summary_table.VesselID, "PA", 'IgnoreCase', true);
depthlogic = summary_table.Depth <70;
lengthlogic = summary_table.bout_duration > 30;
combinedlogic = all([vessellogic, depthlogic,lengthlogic],2);
L1pa_summarytable = summary_table(combinedlogic,:);
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
% Find numeric columns (excluding 'data' which needs special handling)
sessAve_summaryTable = groupsummary(L1pa_summarytable, ["VesselID", "MouseID", "state_name"], "mean", numeric_colnames);
%% Awake Drowsy NREM REM summary stat

awake_summary =  sessAve_summaryTable(sessAve_summaryTable.state_name == "awake",:);
drowsy_summary = sessAve_summaryTable(sessAve_summaryTable.state_name == "drowsy",:);
nrem_summary = sessAve_summaryTable(sessAve_summaryTable.state_name == "nrem",:);
rem_summary = sessAve_summaryTable(sessAve_summaryTable.state_name == "rem",:);

%%

% Normalize data by 'awake' baseline
key_vars = {'MouseID', 'VesselID'};
target_var = 'mean_raw_median';

normalize_table = @(tbl, baseline) ...
    innerjoin(tbl, baseline, 'Keys', key_vars, 'RightVariables', {target_var}, 'LeftVariables', {target_var});

% Join tables to align rows
drowsy_joined = normalize_table(drowsy_summary, awake_summary);
nrem_joined = normalize_table(nrem_summary, awake_summary);
rem_joined = normalize_table(rem_summary, awake_summary);

%% Calculate normalized values
norm_drowsy = drowsy_joined.mean_raw_median_tbl ./ drowsy_joined.mean_raw_median_baseline;
norm_nrem = nrem_joined.mean_raw_median_tbl ./ nrem_joined.mean_raw_median_baseline;
norm_rem = rem_joined.mean_raw_median_tbl ./ rem_joined.mean_raw_median_baseline;

% % Combine into a matrix (Note: rows might not align if vessels are missing in some states)
% % Better to work with tables or ensure equal rows
%%
[mean(norm_drowsy) mean(norm_nrem) mean(norm_rem)]

%%
figure()

plot(mean_dnr,"Color",'r')
hold on

plot(norm_dnr(:,:)',"Color",[0.8 0.8 0.8])
