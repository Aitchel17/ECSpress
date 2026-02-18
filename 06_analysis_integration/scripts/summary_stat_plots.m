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
% tmp.dtypelogic = summary_table.DataType == "thickness_bv";
combinedlogic = all(table2array(struct2table(tmp)), 2); % Dynamically Combine all logic fields in tmp struct 
L1pa_summarytable = summary_table(combinedlogic,:); % apply logic

% Need to calculate each raw_data last 50% and add column
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

%% Match the unit
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var', 'raw_last50_mean'};
for i = 1:length(numeric_colnames)
L1pa_summarytable.(numeric_colnames{i}) = L1pa_summarytable.(numeric_colnames{i}) .* L1pa_summarytable.NumericResolution;
end
% Find numeric columns (excluding 'data' which needs special handling)

sessAve_summaryTable = groupsummary(L1pa_summarytable, ["VesselID", "MouseID", "state_name"], "mean", numeric_colnames);

% Reshape to wide format with states as columns
% Unstack the table to have states as columns for 'mean_raw_last50_mean'
wide_summary = unstack(sessAve_summaryTable(:, {'VesselID', 'MouseID', 'state_name', 'mean_raw_last50_mean'}), 'mean_raw_last50_mean', 'state_name');

% Calculate normalized values (relative to 'awake')
% Normalize by awake baseline
if ismember('awake', wide_summary.Properties.VariableNames)
    norm_drowsy = wide_summary.drowsy   %./ wide_summary.awake;
    % Drowsy normalization might produce Inf if awake is 0, handle if necessary (unlikely for thickness)
    
    norm_nrem = wide_summary.nrem  %./ wide_summary.awake;
    norm_rem = wide_summary.rem %./ wide_summary.awake;
    
    norm_dnr = [norm_drowsy, norm_nrem, norm_rem];
else
    warning('Awake state not found in data, cannot normalize.');
    norm_dnr = [];
end

% Calculate stats
n = sum(~isnan(norm_dnr));
sem_dnr = std(norm_dnr, 0, 1, 'omitnan') ./ sqrt(n);
ts = tinv(0.975, n-1); % 95% CI
ci95 = sem_dnr .* ts;

mean_dnr = mean(norm_dnr, 1, 'omitnan');

disp('Means (Drowsy, NREM, REM):');
disp(mean_dnr);
%%
figure('Name', 'State Summary Plot', 'Color', 'w');

x = 1:3;
y = mean_dnr;
yu = y + ci95;
yl = y - ci95;

% Plot CI patch
patch([x fliplr(x)], [yu fliplr(yl)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on

% Plot mean line
plot(x, y, "Color", 'b', 'LineWidth', 2, 'Marker', 'o', 'MarkerFaceColor', 'b');

% Plot individual lines (optional, commented out as per request for cleaner plot with CI, 
% or keep if "similar to line 67" implies keeping individual traces but adding CI. 
% User said "I want to have confidence interval as patch", implying addition or replacement.
% Usually patch replaces the spaghetti plot or sits behind it. I'll put it behind.
% If user wants "similar to line 67" which had `plot(norm_dnr'...)`, maybe they want both?
% "The final figure I want to have is similar to line 67~ but I want to have confidence interval..."
% I'll keep the individual traces but make them very faint if they want "information". 
% Actually, "have information" might refer to the title. 
% I will exclude individual traces to make the CI clear, as is standard with patches.
% If they want individual traces, they can uncomment.

xticks([1 2 3])
xticklabels({'Drowsy', 'NREM', 'REM'})
ylabel('Median Thickness (\mum)')

% Dynamic Title
if ~isempty(L1pa_summarytable)
    dtype = string(L1pa_summarytable.DataType(1));
else
    dtype = "Unknown";
end
title(replace(dtype, '_', ' ') + " (Mean \pm 95% CI)")

hold off
grid on

%% Statistical Tests

% Bonferroni Correction
% Total number of tests: 3 (vs 1) + 3 (between states) = 6
n_tests = 6;
alpha = 0.05;
alpha_corr = alpha / n_tests;

disp('--------------------------------------------------');
fprintf('Statistical Analysis (Bonferroni Corrected alpha = %.4f for %d tests)\n', alpha_corr, n_tests);
disp('--------------------------------------------------');

disp('1. One-sample t-test vs 1 [Awake Baseline]');
states = {'Drowsy', 'NREM', 'REM'};
for i = 1:3
    [h, p, ci, stats] = ttest(norm_dnr(:, i), 1);
    p_corr = min(p * n_tests, 1); % Bonferroni corrected p-value
    sig_str = "";
    if p < alpha_corr
        sig_str = "*";
    end
    fprintf('%s vs 1: p_raw = %.4f, p_corr = %.4f, t(%d) = %.4f %s\n', states{i}, p, p_corr, stats.df, stats.tstat, sig_str);
end

disp('--------------------------------------------------');
disp('2. Paired t-tests between states');
pairs = [1 2; 2 3; 1 3]; % Drowsy-NREM, NREM-REM, Drowsy-REM
for i = 1:size(pairs, 1)
    s1 = pairs(i, 1);
    s2 = pairs(i, 2);
    [h, p, ci, stats] = ttest(norm_dnr(:, s1), norm_dnr(:, s2));
    p_corr = min(p * n_tests, 1); % Bonferroni corrected p-value
    sig_str = "";
    if p < alpha_corr
        sig_str = "*";
    end
    fprintf('%s vs %s: p_raw = %.4f, p_corr = %.4f, t(%d) = %.4f %s\n', states{s1}, states{s2}, p, p_corr, stats.df, stats.tstat, sig_str);
end
disp('--------------------------------------------------');



