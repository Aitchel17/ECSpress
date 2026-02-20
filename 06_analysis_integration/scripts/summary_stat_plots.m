%% Load
clc, clear
masterDirTable_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl\00_igkl_dirtable.xlsx';
mtable_FWHMsleep = tableManager.load_recon(masterDirTable_path, "mtable_FWHMsleep.mat");
mtable_FWHMsleep.analysis_table = mtable_FWHMsleep.subTables.state_summary; % State summary to be target table
%%
mtable_FWHMsleep.filtLogics = [];
mtable_FWHMsleep.filtLogics.vType = contains(mtable_FWHMsleep.analysis_table.VesselID, "PA", 'IgnoreCase', true); % Artery filter
mtable_FWHMsleep.filtLogics.depth = mtable_FWHMsleep.analysis_table.NumericDepth <70; % L1 filter
mtable_FWHMsleep.filtLogics.bout_dur = mtable_FWHMsleep.analysis_table.bout_duration > 30; % 30 seconds sleep filter
mtable_FWHMsleep.filtLogics.state = ismember(mtable_FWHMsleep.analysis_table.state_name,["awake","drowsy","nrem","rem"] );
%mtable_FWHMsleep.filtLogics.dtype = mtable_FWHMsleep.analysis_table.DataType == "thickness_bv"; % target analysis: thickness_bv
mtable_FWHMsleep.apply_filter
%%
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
mtable_FWHMsleep.apply_resolution("NumericResolution",data_colnames,numeric_colnames);
%%
mtable_FWHMsleep.meanFrom2("raw_data","Q2Q3_mean",0.25,0.75)
mtable_FWHMsleep.addPrctilecol("raw_data","prctile_95", 95);
mtable_FWHMsleep.addPrctilecol("raw_data","prctile_5",5);
%%
mtable_FWHMsleep.get_numericsummary("Date","filtered_table")
%%
mtable_FWHMsleep.get_numericsummary("VesselID","Date_ave")
mtable_FWHMsleep.get_numericsummary("MouseID","VesselID_ave")
%%

mtable = mtable_FWHMsleep.numeric_tables.MouseID_ave;
%%


bv_table = ;
%%
figure("Name","")
%%
cla
plot(mtable(mtable.DataType == "thickness_bv",:),"Q2Q3_mean",color = 'r',Marker='x')
hold on
plot(mtable(mtable.DataType == "thickness_totalpvs",:),"Q2Q3_mean",color = 'g',Marker='x')
ylabel("Absolute thickness (\mum)")
%%
figure("Name","")

plot(mtable(mtable.DataType == "thickness_eps",:),"Q2Q3_mean",color = 'k',Marker='x')






%%
ttable = mtable_FWHMsleep.numeric_tables.mean(mtable_FWHMsleep.numeric_tables.mean.DataType =="thickness_bv",:);


%%


%% Find numeric columns (excluding 'data' which needs special handling)
sessAve_summaryTable = groupsummary(L1pa_summarytable, ["VesselID", "MouseID", "state_name","DataType"], "mean", numeric_colnames);

%% Reshape to wide format with states as columns
% Unstack the table to have states as columns for 'mean_raw_last50_mean'
summary_struct = struct();
dtype_names =  unique(sessAve_summaryTable.DataType);

for dtype_idx = 1:numel(dtype_names)
    target_table = sessAve_summaryTable(sessAve_summaryTable.DataType == dtype_names(dtype_idx),:);
    summary_struct.(dtype_names(dtype_idx)) = unstack(target_table(:, {'VesselID', 'MouseID', 'state_name', 'mean_raw_last50_mean'}), 'mean_raw_last50_mean', 'state_name');
end

%%
wide_summary = unstack(sessAve_summaryTable(:, {'VesselID', 'MouseID', 'state_name', 'mean_raw_last50_mean'}), 'mean_raw_last50_mean', 'state_name');

%% Calculate normalized values (relative to 'awake')
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



