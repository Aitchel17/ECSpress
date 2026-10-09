%INTEGRATE_ANALYSISRESULT  Runs the dataset chain one cell per stage: directory sheets, centralized products, state analysis, tables, figures.
%   Each cell runs the stage script that owns the work; nothing is computed here and
%   nothing writes into the data tree. Cell 0 exports ECSPRESS_ROOT and ECSPRESS_DATASET,
%   which survive the stages' clear, so every stage and figure follows it.

%% 0. Configuration -- the only cell that decides anything
clc, clear
dirs = dirs_central();
param.dataset = "merged";   % str  the dataset; its cohorts are dirtable_main's param.dataset_cohort

addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
setenv('ECSPRESS_ROOT',    dirs.secondary_root);
setenv('ECSPRESS_DATASET', param.dataset);

dataset_dir = fullfile(dirs.secondary_root, param.dataset);
fprintf('root    %s\n', dirs.secondary_root);
fprintf('data    %s\n', strjoin(dirs.working_root, ' + '));
fprintf('dataset %s\n', param.dataset);

%% 0.1 What already exists -- run this before anything, it writes nothing
central_dir = fullfile(dirs.secondary_root, 'centralized');
central_products = ["paxfwhm", "polarcluster", "roilist", "sleep_score", ...
    "analysis_analog", "bout", "drift", "paxfwhm_state", "drift_state"];
check_paths = [ ...
    fullfile(dataset_dir, param.dataset + "_dirtable.xlsx"), ...
    fullfile(central_dir, "centralized_" + central_products + ".mat"), ...
    fullfile(dataset_dir, ["state_summary.mat", "transition.mat", ...
    "drift_summary.mat", "drift_transition.mat"])];
fprintf('\n');
check_exists = util_checkdirstat(check_paths);
fprintf('%d of %d missing -- each cell title below says what it writes\n', ...
    sum(~check_exists), numel(check_exists));

%% 1. Directory tables   [WRITES sheets under dirs.secondary_root]
% every stage below reads the dataset sheet this writes, so a new recording waits for it
run(which('dirtable_main.m'));

%% 1.1 Field drift in stored objects -- run when stage 4 fails on a flatten
% one paxfwhm.mat saved before getdiameter gained eps is enough to break explode_nest;
% the dry run lists what would change and writes nothing
dirs = dirs_central();   % again: cell 1 ran a script that clears the workspace
param.rederive_cohort = "01_igkltdt";   % str  the cohort folder the dry run walks, one per run
cohort_dir = fullfile(dirs.working_root, param.rederive_cohort);   % one per root
for root_idx = find(isfolder(cohort_dir))
    fwhm_rederive(char(cohort_dir(root_idx)), dryrun = true);
end

%% 2. Gather the per-session products   [WRITES under dirs.secondary_root]
% the five .mat products out of the session folders; a row is re-read only when its file moved.
% First run takes minutes and several GB, a re-run with nothing moved takes seconds
run(which('centralize_primary.m'));

%% 2.1 Bouts   [WRITES centralized_bout.mat]
% the state and transition windows of every session, from the scoring cell 2 gathered;
% both state stages below read them, so a re-scored session needs this before cell 3
run(which('centralize_bout.m'));

%% 2.2 Global motion   [WRITES centralized_drift.mat]
% each session's *_motion.txt, the frame shift the extraction applied
run(which('centralize_drift.m'));

%% 3. State analysis for every session   [WRITES centralized_paxfwhm_state.mat]
% from the centralized inputs only; a row is recomputed when its source or its bouts moved
run(which('centralize_paxfwhmstate.m'));

%% 3.1 State analysis of the global motion   [WRITES centralized_drift_state.mat]
run(which('centralize_driftstate.m'));

%% 4. Summary and transition tables   [WRITES state_summary.mat, transition.mat]
% the analysis filters (artery, depth, bout length, states) are set in tablegeneration_main
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

%% 4.2 Global motion tables   [WRITES drift_summary.mat, drift_transition.mat]
% cell 4's filters, on recordings the FWHM pipeline judged; reads what cells 3 and 3.1 wrote
run(which('tablegeneration_drift.m'));

%% 5. Figures   [READ ONLY]
% each reads ECSPRESS_DATASET from cell 0 and none depends on another; only polar_ave_fig writes files

%% 5.1 Transition traces, four transitions x three vessels
run(which('transition_fig.m'));

%% 5.2 Pre/post scatter for the same transitions
run(which('transition_prepost.m'));

%% 5.3 Tonic state summary
run(which('tonic_state.m'));

%% 5.4 Vessel-averaged polar contours   [WRITES svg]
% needs polar_ave.mat from cell 5.0, which walks the session folders itself
run(which('polar_ave_fig.m'));

%% 5.0 Polar contour extraction, only when polar_ave.mat is missing   [WRITES]
% minutes; rerun only when the contours themselves change
run(which('polarcluster_integration.m'));
