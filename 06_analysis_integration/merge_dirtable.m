%MERGE_DIRTABLE  Two cohort dirtables -> one merged sheet -> one mtable.
%   Stage 3 of the FWHM table chain. map_directory writes a dirtable per cohort
%   and batch_state_analysis writes paxfwhm_state.mat per session; this joins the
%   cohorts and runs the secondary_main chain on the join. tablegeneration_main
%   then reads what this leaves behind.
%
%   why  merging at the DIRTABLE level and not at the mtable : explode_nest
%        stamps SessionIndex as the row number of the sheet it walked
%        (explode_nest line 39). One merged sheet therefore gives a consistent
%        SessionIndex for free and no class has to change. Merging two mtables
%        would leave two independent 1..N indices pointing at different rows
%   why  a script and not a function : it is the same shape as secondary_main and
%        is read as a record of which cohorts went into a dataset
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
param.nest_names   = ["state_summary", "transition", "band_decomposition", ...
                "powerdensity", "peak_trough"];

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

has_fwhm  = string(merged_table.Primary_paxFWHM) == "paxfwhm.mat";
has_state = string(merged_table.State_PaxFWHM)   == "paxfwhm_state.mat";
fprintf('  paxfwhm %d | paxfwhm_state %d | both %d\n', ...
    nnz(has_fwhm), nnz(has_state), nnz(has_fwhm & has_state));

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

%% 3. The secondary_main chain, on the merged sheet
% note  identical to secondary_main.m except for the path it is pointed at. If
%       that file's steps change, change them here too
mtable_FWHMsleep = tableManager(char(sheet_out));
mtable_FWHMsleep.filter_refTable("Primary_paxFWHM", "paxfwhm");
mtable_FWHMsleep.filter_refTable("State_PaxFWHM",   "paxfwhm");
mtable_FWHMsleep.parseParams();
mtable_FWHMsleep.aggregateData("State_PaxFWHM");
mtable_FWHMsleep.addnest2subtable(param.nest_names);
mtable_FWHMsleep.save2disk("mtable_FWHMsleep.mat");
% save2disk empties subTables before writing -- they are derived from
% aggregatedTable and load2disk rebuilds them. tableManager is a handle class, so
% the object left in the workspace is emptied too and everything after this line
% would read []. see CLAUDE_LOG.md
mtable_FWHMsleep.addnest2subtable(param.nest_names);

%% 4. What came out
fprintf('\nrefTable %d rows | mice %s\n', height(mtable_FWHMsleep.refTable), ...
    strjoin(unique(string(mtable_FWHMsleep.refTable.MouseID))', ', '));
fprintf('subTables : %s\n', strjoin(fieldnames(mtable_FWHMsleep.subTables)', ', '));
transition_table = mtable_FWHMsleep.subTables.transition;
fprintf('transition %d rows | SessionIndex 1..%d | %d DataTypes\n', ...
    height(transition_table), max(transition_table.SessionIndex), ...
    numel(unique(string(transition_table.DataType))));
% caution  more than one trace length here means a session ran at a different fps
%          or under older code. tableAnalyzer.match_tracelength downsamples to the
%          shortest, never up, but a surprise here is worth chasing to its source
fprintf('trace lengths present : %s\n', ...
    mat2str(unique(cellfun(@numel, transition_table.data))'));
