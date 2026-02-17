%% Load
clc, clear
masterDirTable_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl\00_igkl_dirtable.xlsx';
mtable_FWHMsleep = tableManager.load_recon(masterDirTable_path, "mtable_FWHMsleep.mat");
summary_table = mtable_FWHMsleep.subTables.state_summary; % load state summary
% Filter table
% logic generation
tmp.vessellogic = contains(summary_table.VesselID, "PA", 'IgnoreCase', true); % Artery filter
tmp.depthlogic = summary_table.NumericDepth <70; % L1 filter
tmp.lengthlogic = summary_table.bout_duration > 30; % 30 seconds sleep filter
tmp.dtypelogic = summary_table.DataType == "thickness_bv";
combinedlogic = all([tmp.vessellogic, tmp.depthlogic,tmp.lengthlogic,tmp.dtypelogic],2);  % make combined logic
L1pa_summarytable = summary_table(combinedlogic,:); % apply logic

% Need to calculate each raw_data last 50% and add column
raw_last50_mean = zeros(height(L1pa_summarytable), 1);
for i = 1:height(L1pa_summarytable)
    trace = L1pa_summarytable.raw_data{i};
    num_points = length(trace);
    start_idx = floor(num_points * 0.25) + 1;
    end_idx = floor(num_points * 0.75);
    raw_last50_mean(i) = mean(trace(start_idx:end_idx), 'omitnan') * L1pa_summarytable.NumericResolution(i);
end
L1pa_summarytable.raw_last50_mean = raw_last50_mean;

% Match the unit
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var', 'raw_last50_mean'};
for i = 1:length(numeric_colnames)
L1pa_summarytable.(numeric_colnames{i}) = L1pa_summarytable.(numeric_colnames{i}) .* L1pa_summarytable.NumericResolution;
end
%%




L1pa_summarytable(L1pa_summarytable.VesselID == "PA02",:)


%% Find numeric columns (excluding 'data' which needs special handling)

sessAve_summaryTable = groupsummary(L1pa_summarytable, ["VesselID", "MouseID", "state_name"], "mean", numeric_colnames);

% Reshape to wide format with states as columns
% Unstack the table to have states as columns for 'mean_raw_last50_mean'
wide_summary = unstack(sessAve_summaryTable(:, {'VesselID', 'MouseID', 'state_name', 'mean_raw_last50_mean'}), 'mean_raw_last50_mean', 'state_name');

%% Calculate normalized values (relative to 'awake')
% Normalize by awake baseline
if ismember('awake', wide_summary.Properties.VariableNames)
    norm_drowsy = wide_summary.drowsy ./ wide_summary.awake;
    % Drowsy normalization might produce Inf if awake is 0, handle if necessary (unlikely for thickness)
    
    norm_nrem = wide_summary.nrem ./ wide_summary.awake;
    norm_rem = wide_summary.rem ./ wide_summary.awake;
    
    norm_dnr = [norm_drowsy, norm_nrem, norm_rem];
else
    warning('Awake state not found in data, cannot normalize.');
    norm_dnr = [];
end

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



