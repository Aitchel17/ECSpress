%% Load
clc, clear
masterDirTable_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl\00_igkl_dirtable.xlsx';
mtable_FWHMsleep = tableManager.load_recon(masterDirTable_path, "mtable_FWHMsleep.mat");
summary_table = mtable_FWHMsleep.subTables.state_summary; % load state summary
%% Filter table
% logic generation
tmp.vessellogic = contains(summary_table.VesselID, "PA", 'IgnoreCase', true); % Artery filter
tmp.depthlogic = summary_table.Depth <70; % L1 filter
tmp.lengthlogic = summary_table.bout_duration > 30; % 30 seconds sleep filter
tmp.dtypelogic = summary_table.DataType == "thickness_bv";
combinedlogic = all([tmp.vessellogic, tmp.depthlogic,tmp.lengthlogic,tmp.dtypelogic],2);  % make combined logic
L1pa_summarytable = summary_table(combinedlogic,:); % apply logic

%% Need to calculate each raw_data last 50% and add column
raw_last50_mean = zeros(height(L1pa_summarytable), 1);
for i = 1:height(L1pa_summarytable)
    trace = L1pa_summarytable.raw_data{i};
    num_points = length(trace);
    start_idx = floor(num_points * 0.25) + 1;
    end_idx = floor(num_points * 0.75);
    raw_last50_mean(i) = mean(trace(start_idx:end_idx), 'omitnan');
end
L1pa_summarytable.raw_last50_mean = raw_last50_mean;
%%


L1pa_summarytable(L1pa_summarytable.VesselID == "PA02",:)


%% Find numeric columns (excluding 'data' which needs special handling)

numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var', 'raw_last50_mean'};
sessAve_summaryTable = groupsummary(L1pa_summarytable, ["VesselID", "MouseID", "state_name"], "mean", numeric_colnames);
%%














%% Awake Drowsy NREM REM summary stat

awake_summary =  sessAve_summaryTable(sessAve_summaryTable.state_name == "awake",:);
drowsy_summary = sessAve_summaryTable(sessAve_summaryTable.state_name == "drowsy",:);
nrem_summary = sessAve_summaryTable(sessAve_summaryTable.state_name == "nrem",:);
rem_summary = sessAve_summaryTable(sessAve_summaryTable.state_name == "rem",:);

%%

% Normalize data by 'awake' baseline
key_vars = {'MouseID', 'VesselID'};
target_var = 'mean_raw_last50_mean';

% Helper function to normalize a table against awake_summary using Left Outer Join
% This preserves all rows from the target table, filling missing baseline data with NaN
normalize_table = @(tbl, baseline) ...
    outerjoin(tbl, baseline, 'Keys', key_vars, 'Type', 'left', 'MergeKeys', true, ...
    'RightVariables', {target_var}, 'LeftVariables', {target_var});

% Join tables to align rows (Left Join ensures we keep all state rows)
drowsy_joined = normalize_table(drowsy_summary, awake_summary);
nrem_joined = normalize_table(nrem_summary, awake_summary);
rem_joined = normalize_table(rem_summary, awake_summary);

%% Calculate normalized values (will range from 0 to Inf, or NaN if baseline missing)


% Calculate ratios
norm_drowsy = drowsy_joined.mean_raw_last50_mean_left ./ drowsy_joined.mean_raw_last50_mean_right;
norm_nrem = nrem_joined.mean_raw_last50_mean_left ./ nrem_joined.mean_raw_last50_mean_right;
norm_rem = rem_joined.mean_raw_last50_mean_left ./ rem_joined.mean_raw_last50_mean_right;

%% Combine into a single matrix for plotting (aligning by joining all to a master list)
% To get a [N x 3] matrix where rows align, we need a master list of all unique vessels across all states
all_keys = unique([drowsy_summary(:, key_vars); nrem_summary(:, key_vars); rem_summary(:, key_vars); awake_summary(:, key_vars)]);

% Join each state to the master list
drowsy_all = outerjoin(all_keys, drowsy_joined, 'Keys', key_vars, 'Type', 'left', 'MergeKeys', true);
nrem_all   = outerjoin(all_keys, nrem_joined, 'Keys', key_vars, 'Type', 'left', 'MergeKeys', true);
rem_all    = outerjoin(all_keys, rem_joined, 'Keys', key_vars, 'Type', 'left', 'MergeKeys', true);

% Re-calculate ratios ensuring we use the aligned tables
nd = drowsy_all.mean_raw_last50_mean_left ./ drowsy_all.mean_raw_last50_mean_right;
nn = nrem_all.mean_raw_last50_mean_left ./ nrem_all.mean_raw_last50_mean_right;
nr = rem_all.mean_raw_last50_mean_left ./ rem_all.mean_raw_last50_mean_right;

norm_dnr = [nd, nn, nr];
mean_dnr = mean(norm_dnr, 1, 'omitnan');

disp('Normalized Means (Drowsy, NREM, REM):');
disp(mean_dnr);
%%
figure()

plot(mean_dnr,"Color",'r', 'LineWidth', 2)
hold on
plot(norm_dnr',"Color",[0.8 0.8 0.8])
plot(mean_dnr,"Color",'r', 'LineWidth', 2) % Replot mean on top
xticks([1 2 3])
xticklabels({'Drowsy', 'NREM', 'REM'})
ylabel('Normalized Median Thickness')
title('Vessel Thickness Normalized to Awake')
hold off



