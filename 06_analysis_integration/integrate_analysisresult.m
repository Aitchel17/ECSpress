%INTEGRATE_ANALYSISRESULT  Directory table -> merged dataset -> figures, block by block.
%   The whole secondary path in one file. Every cell runs on its own, in order or
%   not, and says what it needs before it needs it. Nothing here computes: each
%   cell calls the stage that owns the work.
%
%   THE CHAIN, and what each stage leaves on disk
%     1  map_directory        <cohort>_dirtable.xlsx        per cohort, in dirs.secondary_root
%     2  merge_dirtable       <param.dataset>_dirtable.xlsx
%     3  centralize_primary   centralized/centralized_<product>.mat, five of them
%     4  centralize_state     centralized/centralized_paxfwhm_state.mat
%     5  tablegeneration_main state_summary.mat, transition.mat
%     6  figures              svg, beside each session or in the dataset folder
%
%   Stages 1-5 WRITE, and ONLY under dirs.secondary_root. Nothing below stage 1
%   opens a session folder: 3 is the last stage that reads the data tree at all,
%   and 4 onward work from what it gathered. Cell 0 is the only place a path is
%   decided, so pointing this at another dataset is one edit.
%
%   Stages 3 and 4 are INCREMENTAL. Each stamps every row with the size and
%   modified time of the file it was built from and recomputes only the rows whose
%   source has moved since. A re-run that changed nothing costs seconds.
%     err  the state analysis used to be written into the session folders by
%          batch_state_analysis, which skipped a session whenever the file already
%          existed -- so a re-scored recording kept its old answer and no one could
%          see it. Two sessions were found that way. see CLAUDE_LOG.md
%
%   Re-running a stage overwrites its output. Stages 1-2 back up first; 3-5 do not
%   -- 3 and 4 are a cache of files that still exist, and 5 is a pure function of 4
%   and of the filters written into it.
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

% The joined dataset stages 2-6 work on
param.dataset = "merged_igkl_igkltdt";

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
central_dir = fullfile(dirs.secondary_root, 'centralized');
central_products = ["paxfwhm", "polarcluster", "roilist", "sleep_score", ...
    "analysis_analog", "paxfwhm_state"];
check_paths = [ ...
    fullfile(dirs.secondary_root, param.cohorts, param.cohorts + "_dirtable.xlsx"), ...
    fullfile(dataset_dir, param.dataset + "_dirtable.xlsx"), ...
    fullfile(central_dir, "centralized_" + central_products + ".mat"), ...
    fullfile(dataset_dir, ["state_summary.mat", "transition.mat"])];
fprintf('\n');
check_exists = util_checkdirstat(check_paths);
fprintf('%d of %d missing -- the stage that writes each is in the header above\n', ...
    sum(~check_exists), numel(check_exists));

%% 1. Directory table, per cohort   [WRITES to dirs.secondary_root and to the data tree]
% Walks one cohort, records which analysis files each session has, and writes the
% sheet every later stage keys on.
%   Set param.cohort_idx and run this cell once per cohort.
param.cohort_idx = 2;

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

%% 1.1 Field drift in stored objects -- run when stage 5 fails on a flatten
% One session saved before getdiameter gained eps/epschanges is enough to make
% struct2table wrap every row in a cell, and explode_nest then fails.
% see CLAUDE_LOG.md. dryrun lists what would change and writes nothing.
fwhm_rederive(char(fullfile(dirs.working_root, param.cohorts(param.cohort_idx))), dryrun = true);

%% 2. Join the cohort sheets into one   [WRITES, backs up first]
% One sheet for the whole dataset, which is what centralize_primary walks.
% SOURCES is still set inside merge_dirtable -- change it there.
run(fullfile('g:\03_program\01_ecspress\06_analysis_integration', 'merge_dirtable.m'));

%% 3. Gather the per-session products   [WRITES under dirs.secondary_root\centralized]
% Reads paxfwhm.mat, polarcluster.mat, roilist.mat, sleep_score.mat and
% analysis_analog.mat out of every session the merged sheet lists, and writes one
% table per product. The LAST stage that opens a session folder. Incremental: a row
% is rebuilt only when its source file's size or modified time has moved.
%   First run is minutes and writes several GB. A re-run that finds nothing changed
%   is seconds and does not rewrite the files.
run(fullfile('g:\03_program\01_ecspress\06_analysis_integration', 'centralize_primary.m'));

%% 4. State analysis for every session   [WRITES centralized_paxfwhm_state.mat]
% state_integration -> state_linefwhm over the nine series, computed from the
% centralized inputs. Nothing is read from or written to the data tree.
%   This replaced batch_state_analysis, which wrote paxfwhm_state.mat into each
%   session folder and skipped any session that already had one -- so a re-scored
%   recording silently kept its old answer. see CLAUDE_LOG.md
run(fullfile('g:\03_program\01_ecspress\06_analysis_integration', 'centralize_state.m'));

%% 5. Summary and transition tables   [WRITES state_summary.mat, transition.mat]
% Flattens the two nests out of the centralized state analysis and applies the
% filters: PA only, depth < 70, bout > 30 s, four states, sleep sessions only for
% the transitions. Those are the analysis decision, so they stay in that file
% rather than moving up here.
%   Nothing caches what this produces. state_summary/transition are pure functions
%   of the centralized state table and of those filters, so the committed script IS
%   the record -- as long as it is committed. see CLAUDE_LOG.md
run(fullfile('g:\03_program\01_ecspress\06_analysis_integration', 'tablegeneration_main.m'));

%% 5.1 What came out
summary_path    = fullfile(dataset_dir, "state_summary.mat");
transition_path = fullfile(dataset_dir, "transition.mat");
if ~isfile(transition_path)
    error('integrate:noTables', 'run cell 5 first; %s is missing', transition_path);
end
loaded_transition = load(transition_path);
transition_result = loaded_transition.save_content;
n_row    = height(transition_result.abs_numeric.filtered_table);
n_mouse  = numel(unique(string(transition_result.abs_numeric.MouseID_ave.MouseID)));
n_vessel = height(transition_result.abs_numeric.VesselID_ave);
fprintf('transition : %d rows | %d vessels | %d mice\n', n_row, n_vessel, n_mouse);
fprintf('trace lengths : %s\n', ...
    mat2str(unique(cellfun(@numel, transition_result.abs_numeric.filtered_table.data))'));

%% 6. Figures   [READ ONLY]
% Each is standalone and reads ECSPRESS_DATASET, which cell 0 set. Run a cell on
% its own or the block; they do not depend on each other.
%   Only polar_ave_fig writes files. The other three draw on screen.

%% 6.1 Transition traces, four transitions x three vessels
run(fullfile('g:\03_program\01_ecspress\07_figuregeneration', 'transition_fig.m'));

%% 6.2 Pre/post scatter for the same transitions
run(fullfile('g:\03_program\01_ecspress\07_figuregeneration', 'transition_prepost.m'));

%% 6.3 Tonic state summary
run(fullfile('g:\03_program\01_ecspress\07_figuregeneration', 'tonic_state.m'));

%% 6.4 Vessel-averaged polar contours   [WRITES svg]
% Needs polar_ave.mat, which polarcluster_integration writes and which is NOT
% part of the table chain -- it walks the sessions itself. Run that first if the
% file is missing.
run(fullfile('g:\03_program\01_ecspress\07_figuregeneration', 'polar_ave_fig.m'));

%% 6.0 Polar contour extraction, only when polar_ave.mat is missing   [WRITES]
% Minutes, not seconds. polar_ave_fig does the alignment and averaging, so this
% only has to run when the contours themselves change.
run(fullfile('g:\03_program\01_ecspress\06_analysis_integration', ...
    'polarcluster_integration.m'));
