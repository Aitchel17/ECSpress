%MERGE_DIRTABLE  Two cohort dirtables -> one merged sheet.
%   map_directory writes a dirtable per cohort; this joins them into the one sheet
%   centralize_primary walks. That is all it does now -- it used to go on and build
%   mtable_FWHMsleep.mat by opening every session's paxfwhm_state.mat, and that job
%   moved to centralize_primary and centralize_state, which are incremental.
%
%   why  a script and not a function : it is read as a record of which cohorts went
%        into a dataset
%   caution  this is the ONLY place the merged dataset comes from. It lived in a
%        scratch folder until 2026-08-08, which meant the 32-session result could
%        not be rebuilt from the repo

clc, clear
addpath(genpath('g:\03_program\01_ecspress'));

% integrate_analysisresult exports these; they survive the clear above. Run this
% file on its own and the defaults apply.
dirs.secondary_root = getenv('ECSPRESS_ROOT');
if isempty(dirs.secondary_root)
    dirs.secondary_root = ['E:\OneDrive - The Pennsylvania State University\' ...
        '2023ecspress\02_secondary_analysis'];
end
param.dataset = string(getenv('ECSPRESS_DATASET'));   % output folder and sheet name
if param.dataset == ""
    param.dataset = "merged_igkl_igkltdt";
end
param.sources      = ["00_igkl", "01_igkltdt"];   % cohort folders, each holding <name>_dirtable.xlsx
dirs.working_root = "G:\tmp";                    % the live data tree

%% 1. Join the cohort sheets
% err  E:\ is a BACKUP drive, not a second data tree. write_dirtable merges an old
%      sheet with a fresh scan, so an igkltdt sheet that was once written while
%      pointed at E:\ keeps those rows, and they duplicate sessions that also sit
%      under G:\tmp. Every row not on dirs.working_root is dropped
part = cell(1, numel(param.sources));
for source_idx = 1:numel(param.sources)
    sheet_path = fullfile(dirs.secondary_root, param.sources(source_idx), ...
        param.sources(source_idx) + "_dirtable.xlsx");
    import_opts = detectImportOptions(sheet_path, 'Sheet', 'reference');
    cohort_table = readtable(sheet_path, import_opts);
    on_working = startsWith(string(cohort_table.Directory), dirs.working_root);
    % which line a row came from, kept for per-cohort counts downstream
    cohort_table.Cohort = repmat(param.sources(source_idx), height(cohort_table), 1);
    fprintf('%-12s %4d rows | on %s %4d | dropped %d\n', param.sources(source_idx), ...
        height(cohort_table), dirs.working_root, nnz(on_working), nnz(~on_working));
    part{source_idx} = cohort_table(on_working, :);
end
merged_table = vertcat(part{:});

% err  the key is Directory, NOT MouseID|Date|SessionID : sibling folders share a
%      session number (HQL090_sleep251012_005_notanalyzable_branch and the like),
%      so a three-column key silently collapses two different recordings
[~, first_of_directory] = unique(string(merged_table.Directory), 'stable');
if numel(first_of_directory) ~= height(merged_table)
    fprintf('dropping %d duplicate Directory rows\n', ...
        height(merged_table) - numel(first_of_directory));
    merged_table = merged_table(first_of_directory, :);
end
fprintf('\nmerged %d rows | %d distinct Directory | mice %s\n', ...
    height(merged_table), numel(unique(string(merged_table.Directory))), ...
    strjoin(unique(string(merged_table.MouseID))', ', '));

% What centralize_primary will find. The State_PaxFWHM column is still written by
% mapdirstruct and still counted here, but nothing reads it any more -- the state
% analysis is computed into centralized_paxfwhm_state.mat, not written per session.
has_fwhm  = string(merged_table.Primary_paxFWHM) == "paxfwhm.mat";
has_score = string(merged_table.Peripheral_SleepScore) == "sleep_score.mat";
fprintf('  paxfwhm %d | sleep_score %d | both %d\n', ...
    nnz(has_fwhm), nnz(has_score), nnz(has_fwhm & has_score));

%% 2. Write the merged sheet
% caution  writetable does not clear the sheet first, so a shorter table would
%          leave the old tail behind. The file has to go -- back it up, do not
%          just delete, the way fwhm_rederive does
dataset_dir = fullfile(dirs.secondary_root, param.dataset);
if ~isfolder(dataset_dir); mkdir(dataset_dir); end
sheet_out = fullfile(dataset_dir, param.dataset + "_dirtable.xlsx");
if isfile(sheet_out)
    backup_path = fullfile(dataset_dir, param.dataset + "_dirtable_backup_" + ...
        string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) + ".xlsx");
    copyfile(sheet_out, backup_path);
    delete(sheet_out);
    fprintf('backed up previous sheet to %s\n', backup_path);
end
writetable(merged_table, sheet_out, 'Sheet', 'reference');
fprintf('wrote %s\n', sheet_out);

%% 3. What came out
fprintf('\nmerged sheet %d rows | mice %s\n', height(merged_table), ...
    strjoin(unique(string(merged_table.MouseID))', ', '));
session_type = string(merged_table.SessionType);
type_name = unique(session_type);
for k = 1:numel(type_name)
    fprintf('   %-12s %d\n', type_name(k), sum(session_type == type_name(k)));
end
fprintf('centralize_primary reads this sheet next\n');
