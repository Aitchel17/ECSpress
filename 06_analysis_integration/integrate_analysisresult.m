%INTEGRATE_ANALYSISRESULT  Directory table -> merged dataset -> figures, block by block.
%   The whole secondary path in one file. Every cell runs on its own, in order or
%   not, and says what it needs before it needs it. Nothing here computes: each
%   cell calls the stage that owns the work.
%
%   THE CHAIN, and what each stage leaves on disk
%     1  map_directory        <cohort>_dirtable.xlsx        per cohort, in dirs.secondary_root
%        batch_state_analysis <session>/state_analysis/paxfwhm_state.mat
%     2  secondary_main       mtable_FWHMsleep.mat          per cohort
%     3  merge_dirtable       <param.dataset>_dirtable.xlsx + mtable_FWHMsleep.mat
%     4  tablegeneration_main state_summary.mat, transition.mat
%     5  figures              svg, beside each session or in the dataset folder
%
%   Stages 1-4 WRITE. 5 only reads. Cell 0 is the only place a path is decided;
%   everything below reads it, so pointing this at another dataset is one edit.
%
%   Re-running a stage overwrites its output. Stages 1-3 back up first; stage 4
%   does not -- see the note in its cell.
%
%   The figure scripts are standalone (they clc/clear and load their own table).
%   Cell 0 exports the dataset name to the environment, which survives their
%   clear, so they follow this file without being edited. Run one directly and it
%   falls back to its own default.

%% 0. Configuration -- the only cell that decides anything
clc, clear

% Where the analysis products live. Not the data.
dirs.secondary_root = ['E:\OneDrive - The Pennsylvania State University\' ...
    '2023ecspress\02_secondary_analysis'];

% The live data tree. E:\ is the BACKUP drive -- a cohort path under it produces
% rows that duplicate the same sessions under G:\tmp. see CLAUDE_LOG.md
dirs.working_root = 'G:\tmp';

% Cohorts, as <folder under dirs.working_root> -> <name of its dirtable folder>
param.cohorts = ["00_igkl", "01_igkltdt"];

% The joined dataset stages 3-5 work on
param.dataset = "merged_igkl_igkltdt";

% Nested tables explode_nest pulls out of the aggregated rows
param.nest_names = ["state_summary", "transition", "band_decomposition", ...
              "powerdensity", "peak_trough"];

% What map_directory looks for in each session folder
param.primary_map = { ...
    'RadonResult',  'radon_result.mat'; ...
    'RoiList',      'roilist.mat'; ...
    'PolarCluster', 'polarcluster.mat'; ...
    'PaxFWHM',      'paxfwhm.mat'};
param.peripheral_map = { ...
    'AnalogAnalysis',   'analysis_analog.mat'; ...
    'BehaviorAnalysis', 'analysis_camera.mat'; ...
    'SleepScore',       'sleep_score.mat'};
param.stateanalysis_map = {'paxfwhm_state', 'paxfwhm_state.mat'};

addpath(genpath('g:\03_program\01_ecspress'));
setenv('ECSPRESS_ROOT',    dirs.secondary_root);
setenv('ECSPRESS_DATASET', param.dataset);

dataset_dir = fullfile(dirs.secondary_root, param.dataset);
fprintf('root    %s\n', dirs.secondary_root);
fprintf('data    %s\n', dirs.working_root);
fprintf('dataset %s\n', param.dataset);

%% 0.1 What already exists -- run this before anything, it writes nothing
% A stage whose input is missing fails deep inside a loop. This says so first.
fprintf('\n%-34s %-10s %-12s %s\n', 'file', 'exists', 'size', 'modified');
check_paths = strings(0);
for c = 1:numel(param.cohorts)
    check_paths(end+1) = fullfile(dirs.secondary_root, param.cohorts(c), param.cohorts(c) + "_dirtable.xlsx"); %#ok<SAGROW>
    check_paths(end+1) = fullfile(dirs.secondary_root, param.cohorts(c), "mtable_FWHMsleep.mat");        %#ok<SAGROW>
end
check_paths(end+1) = fullfile(dataset_dir, param.dataset + "_dirtable.xlsx");
check_paths(end+1) = fullfile(dataset_dir, "mtable_FWHMsleep.mat");
check_paths(end+1) = fullfile(dataset_dir, "state_summary.mat");
check_paths(end+1) = fullfile(dataset_dir, "transition.mat");
for k = 1:numel(check_paths)
    [~, nm, ex] = fileparts(check_paths(k));
    if isfile(check_paths(k))
        d = dir(check_paths(k));
        fprintf('%-34s %-10s %8.1f MB  %s\n', nm + ex, 'yes', d.bytes/1e6, ...
            string(d.date));
    else
        fprintf('%-34s %-10s\n', nm + ex, 'MISSING');
    end
end

%% 1. Directory table, per cohort   [WRITES to dirs.secondary_root and to the data tree]
% Walks one cohort, records which analysis files each session has, and writes the
% sheet every later stage keys on. batch_state_analysis then fills in whatever
% paxfwhm_state.mat is missing, which is the slow part.
%   Set param.cohort_idx and run this cell once per cohort.
param.cohort_idx = 1;

cohort       = param.cohorts(param.cohort_idx);
cohort_dir   = fullfile(dirs.working_root, cohort);
cohort_out   = fullfile(dirs.secondary_root, cohort);
dirtable_dir = fullfile(cohort_out, cohort + "_dirtable.xlsx");
dircache_dir = fullfile(cohort_out, cohort + "_dircache.mat");
if ~isfolder(cohort_dir)
    error('integrate:noCohort', 'cohort folder not found: %s', cohort_dir);
end
if ~isfolder(cohort_out)
    mkdir(cohort_out)
end

% Legacy eye.avi / whisker.avi -> prefix_eye.avi. Renames files in the DATA tree
rename_avi_files(char(cohort_dir));

dirstruct_table = mapdirstruct(char(cohort_dir), param.primary_map, param.peripheral_map, ...
    param.stateanalysis_map, char(dircache_dir));
write_dirtable(dirstruct_table, char(dirtable_dir));
fprintf('%s : %d sessions -> %s\n', cohort, height(dirstruct_table), dirtable_dir);

%% 1.1 Batch state analysis   [WRITES paxfwhm_state.mat into every session folder]
% The slow one. Skips sessions that already have their state file.
batch_state_analysis(dirstruct_table, char(cohort_dir), param.primary_map, ...
    param.peripheral_map, param.stateanalysis_map, char(dirtable_dir), true, char(dircache_dir));

%% 1.2 Field drift in stored objects -- run when stage 3 fails on a merge
% One session saved before getdiameter gained eps/epschanges is enough to make
% struct2table wrap every row in a cell, and explode_nest then fails.
% see CLAUDE_LOG.md. dryrun lists what would change and writes nothing.
fwhm_rederive(char(fullfile(dirs.working_root, param.cohorts(param.cohort_idx))), dryrun = true);

%% 2. One cohort through the table chain   [WRITES mtable_FWHMsleep.mat]
% Only needed to look at a cohort on its own. Stage 3 does this again on the join,
% so a merged run does not need it.
param.cohort_idx = 1;

cohort_sheet = fullfile(dirs.secondary_root, param.cohorts(param.cohort_idx), ...
    param.cohorts(param.cohort_idx) + "_dirtable.xlsx");
mtable_cohort = tableManager(char(cohort_sheet));
mtable_cohort.filter_refTable("Primary_paxFWHM", "paxfwhm");
mtable_cohort.filter_refTable("State_PaxFWHM",   "paxfwhm");
mtable_cohort.parseParams();
mtable_cohort.aggregateData("State_PaxFWHM");
mtable_cohort.addnest2subtable(param.nest_names);
mtable_cohort.save2disk("mtable_FWHMsleep.mat");

%% 3. Join the cohorts and build the merged mtable   [WRITES, backs up first]
% Joins at the SHEET level, not at the mtable: explode_nest stamps SessionIndex as
% the row number of the sheet it walked, so one merged sheet gives a consistent
% index for free. SOURCES is still set inside merge_dirtable -- change it there.
run(fullfile('g:\03_program\01_ecspress\06_analysis_integration', 'merge_dirtable.m'));

%% 4. Summary and transition tables   [WRITES state_summary.mat, transition.mat]
% The filters live in tablegeneration_main: PA only, depth < 70, bout > 30 s,
% four states, sleep sessions only. They are the analysis decision, so they stay
% in that file rather than moving up here.
%   No backup is written. state_summary/transition are pure functions of the
%   mtable and of those filters, so the committed script IS the record -- as long
%   as it is committed. see CLAUDE_LOG.md
run(fullfile('g:\03_program\01_ecspress\06_analysis_integration', 'tablegeneration_main.m'));

%% 4.1 What came out
summary_path    = fullfile(dataset_dir, "state_summary.mat");
transition_path = fullfile(dataset_dir, "transition.mat");
if ~isfile(transition_path)
    error('integrate:noTables', 'run cell 4 first; %s is missing', transition_path);
end
loaded_transition = load(transition_path);
transition_result = loaded_transition.save_content;
n_row    = height(transition_result.abs_numeric.filtered_table);
n_mouse  = numel(unique(string(transition_result.abs_numeric.MouseID_ave.MouseID)));
n_vessel = height(transition_result.abs_numeric.VesselID_ave);
fprintf('transition : %d rows | %d vessels | %d mice\n', n_row, n_vessel, n_mouse);
fprintf('trace lengths : %s\n', ...
    mat2str(unique(cellfun(@numel, transition_result.abs_numeric.filtered_table.data))'));

%% 5. Figures   [READ ONLY]
% Each is standalone and reads ECSPRESS_DATASET, which cell 0 set. Run a cell on
% its own or the block; they do not depend on each other.

%% 5.1 Transition traces, four transitions x three vessels
run(fullfile('g:\03_program\01_ecspress\07_figuregeneration', 'transition_fig.m'));

%% 5.2 Pre/post scatter for the same transitions
run(fullfile('g:\03_program\01_ecspress\07_figuregeneration', 'transition_prepost.m'));

%% 5.3 Tonic state summary
run(fullfile('g:\03_program\01_ecspress\07_figuregeneration', 'tonic_state.m'));

%% 5.4 Vessel-averaged polar contours
% Needs polar_ave.mat, which polarcluster_integration writes and which is NOT
% part of the table chain -- it walks the sessions itself. Run that first if the
% file is missing.
run(fullfile('g:\03_program\01_ecspress\07_figuregeneration', 'polar_ave_fig.m'));

%% 5.0 Polar contour extraction, only when polar_ave.mat is missing   [WRITES]
% Minutes, not seconds. polar_ave_fig does the alignment and averaging, so this
% only has to run when the contours themselves change.
run(fullfile('g:\03_program\01_ecspress\06_analysis_integration', ...
    'polarcluster_integration.m'));
