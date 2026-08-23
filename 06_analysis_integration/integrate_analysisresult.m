%INTEGRATE_ANALYSISRESULT  Directory table -> merged dataset -> figures, block by block.
%   The whole secondary path in one file. Every cell runs on its own, in order or
%   not, and says what it needs before it needs it. Nothing here computes: each
%   cell calls the stage that owns the work.
%
%   THE CHAIN, and what each stage leaves on disk
%     1  rename_avi_files     nothing -- it renames video files in the DATA tree
%     2  centralize_primary   <cohort>_dirtable.xlsx, <dataset>_dirtable.xlsx and
%                             centralized/centralized_<product>.mat, five of them
%     3  centralize_state     centralized/centralized_paxfwhm_state.mat
%     4  tablegeneration_main state_summary.mat, transition.mat
%     5  figures              svg in the dataset folder, and figure windows
%
%   CELL 1 IS THE ONLY ONE THAT WRITES INTO THE DATA TREE, and it only renames
%   videos. Everything else writes under dirs.secondary_root alone. Stage 2 is the
%   last stage that so much as opens a session folder; 3 onward work from what it
%   gathered. Cell 0 is the only place a path is decided, so pointing this at
%   another dataset is one edit.
%
%   Stage 2 rebuilds the directory sheets before it reads them, by running
%   merge_dirtable, so a recording added since the last run cannot go unnoticed.
%   The rescan does not overwrite metadata corrected by hand: write_dirtable
%   updates only the file columns of a row that already exists.
%
%   Stages 2 and 3 are INCREMENTAL. Each stamps every row with the size and
%   modified time of the file it was built from and recomputes only the rows whose
%   source has moved since. A re-run that changed nothing costs seconds.
%     err  the state analysis used to be written into the session folders by
%          batch_state_analysis, which skipped a session whenever the file already
%          existed -- so a re-scored recording kept its old answer and no one could
%          see it. Two sessions were found that way. see CLAUDE_LOG.md
%
%   Re-running a stage overwrites its output. The merged sheet is backed up first;
%   the centralized files are not -- they are a cache of files that still exist,
%   and stage 4 is a pure function of stage 3 and of the filters written into it.
%
%   The figure scripts are standalone (they clc/clear and load their own table).
%   Cell 0 exports the dataset name to the environment, which survives their
%   clear, so they follow this file without being edited. Run one directly and it
%   falls back to its own default.

%% 0. Configuration -- the only cell that decides anything
clc, clear

% Where the analysis products live. Not the data.
dirs = util_centraldirs();

% The live data tree comes from util_centraldirs and is "" off the rig.
% E:\ is the BACKUP drive -- a cohort path under it produces rows that
% duplicate the same sessions under G:\tmp. see CLAUDE_LOG.md

% Cohorts, as <folder under dirs.working_root>. Only cell 1 reads this, and only
% to pick which tree to rename videos in. Which cohorts go into the DATASET is
% centralize_primary's param.cohorts, which it hands to merge_dirtable.
param.cohorts = ["00_igkl", "01_igkltdt"];

% The joined dataset stages 2-5 work on
param.dataset = "merged_igkl_igkltdt";

addpath('g:\03_program\01_ecspress\functions');   % where util_ecspath lives
util_ecspath;                                     % three roots, minus zz_notinuse
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

%% 1. Legacy AVI names   [WRITES INTO THE DATA TREE]
% eye.avi / whisker.avi -> prefix_eye.avi, so the scan can tell which session a
% video belongs to. The only cell in this file that touches G:, which is why it is
% on its own and never runs as part of anything else.
%   Set param.cohort_idx and run this cell once per cohort. It only has to run
%   after an acquisition that produced the legacy names.
param.cohort_idx = 2;
cohort_dir = fullfile(dirs.working_root, param.cohorts(param.cohort_idx));
if ~isfolder(cohort_dir)
    error('integrate:noCohort', 'cohort folder not found: %s', cohort_dir);
end
rename_avi_files(char(cohort_dir));

%% 1.1 Field drift in stored objects -- run when stage 4 fails on a flatten
% One session saved before getdiameter gained eps/epschanges is enough to make
% struct2table wrap every row in a cell, and explode_nest then fails.
% see CLAUDE_LOG.md. dryrun lists what would change and writes nothing.
fwhm_rederive(char(fullfile(dirs.working_root, param.cohorts(param.cohort_idx))), dryrun = true);

%% 2. Gather the per-session products   [WRITES under dirs.secondary_root]
% Scans both cohorts and rebuilds the sheets first (merge_dirtable, about a
% second), then reads paxfwhm.mat, polarcluster.mat, roilist.mat, sleep_score.mat
% and analysis_analog.mat out of every session that sheet lists and writes one
% table per product. The LAST stage that opens a session folder. Incremental: a
% row is rebuilt only when its source file's size or modified time has moved.
%   First run is minutes and writes several GB. A re-run that finds nothing
%   changed is seconds and does not rewrite the files.
%   note  the scan does not overwrite a metadata value corrected by hand in the
%        reference sheet -- write_dirtable updates only the file columns
run(which('centralize_primary.m'));

%% 3. State analysis for every session   [WRITES centralized_paxfwhm_state.mat]
% state_integration -> state_linefwhm over the nine series, computed from the
% centralized inputs. Nothing is read from or written to the data tree.
%   This replaced batch_state_analysis, which wrote paxfwhm_state.mat into each
%   session folder and skipped any session that already had one -- so a re-scored
%   recording silently kept its old answer. see CLAUDE_LOG.md
run(which('centralize_state.m'));

%% 4. Summary and transition tables   [WRITES state_summary.mat, transition.mat]
% Flattens the two nests out of the centralized state analysis and applies the
% filters: PA only, depth < 70, bout > 30 s, four states, sleep sessions only for
% the transitions. Those are the analysis decision, so they stay in that file
% rather than moving up here.
%   Nothing caches what this produces. state_summary/transition are pure functions
%   of the centralized state table and of those filters, so the committed script IS
%   the record -- as long as it is committed. see CLAUDE_LOG.md
run(which('tablegeneration_main.m'));

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
%   Only polar_ave_fig writes files. The other three draw on screen.

%% 5.1 Transition traces, four transitions x three vessels
run(which('transition_fig.m'));

%% 5.2 Pre/post scatter for the same transitions
run(which('transition_prepost.m'));

%% 5.3 Tonic state summary
run(which('tonic_state.m'));

%% 5.4 Vessel-averaged polar contours   [WRITES svg]
% Needs polar_ave.mat, which polarcluster_integration writes and which is NOT
% part of the table chain -- it walks the sessions itself. Run that first if the
% file is missing.
run(which('polar_ave_fig.m'));

%% 5.0 Polar contour extraction, only when polar_ave.mat is missing   [WRITES]
% Minutes, not seconds. polar_ave_fig does the alignment and averaging, so this
% only has to run when the contours themselves change.
run(which('polarcluster_integration.m'));
