%REBUILD_STATE_ANALYSIS  Re-run batch_state_analysis for a SUBSET, with backups.
%   map_directory runs the whole tree. This runs the sessions that a specific
%   thing went wrong for, which is what every repair during the 2026-08-06 rebuild
%   actually needed. Pick a reason, look at what it selected, then let it write.
%
%   why      re-running the whole tree to fix a handful of sessions costs most of
%            an hour and rewrites every file that was already right
%            see FINDINGS.md
%   caution  this WRITES paxfwhm_state.mat under the data root. Every selected
%            file is copied to <name>_backup_<backup_stamp>.mat first -- do not remove
%            that, it is the only undo
%   caution  param.dryrun is true. Nothing is written until it is set false, and the
%            selection is printed either way
%   note     paxfwhm.mat itself is NOT touched here. Backfilling thickness.eps is
%            fwhm_rederive's job and has to run BEFORE this

clc, clear
addpath(genpath('g:\03_program\01_ecspress'));

param.dryrun  = true;                        % false = actually re-run and overwrite
backup_stamp   = string(datetime('now', 'Format', 'yyyyMMdd'));
param.dataset = "merged_igkl_igkltdt";
dirs.secondary_root = ['E:\OneDrive - The Pennsylvania State University\' ...
    '2023ecspress\02_secondary_analysis'];
dirs.data_root = 'G:\tmp';

% Which sessions. Each reason is a rule, not a hand-typed list, so re-running it
% later finds whatever is broken THEN rather than what was broken once.
%   'tracelength'  stored transition traces are not 2*round(tw*fs/2)+1 long. fs and
%                  transition_window agree across the tree, so any other length is
%                  the pre-fix code path -- those traces span 50 s where the
%                  current ones span 25, and resampling would put two different
%                  WINDOWS on one axis
%   'missing'      has paxfwhm.mat but no paxfwhm_state.mat
%   'session'      a name substring, for a session rebuilt by hand
param.reason       = "tracelength";
param.session_name = "";                     % only read when param.reason == "session"

% the four maps batch_state_analysis needs, same as map_directory
primary_map = {
    'RadonResult',  'radon_result.mat';
    'RoiList',      'roilist.mat';
    'PolarCluster', 'polarcluster.mat';
    'PaxFWHM',      'paxfwhm.mat'};
peripheral_map = {
    'AnalogAnalysis',   'analysis_analog.mat';
    'BehaviorAnalysis', 'analysis_camera.mat';
    'SleepScore',       'sleep_score.mat'};
stateanalysis_map = {'paxfwhm_state', 'paxfwhm_state.mat'};

%% 1. Read the merged sheet
sheet_path = fullfile(dirs.secondary_root, param.dataset, param.dataset + "_dirtable.xlsx");
import_opts = detectImportOptions(sheet_path, 'Sheet', 'reference');
dir_table = readtable(sheet_path, import_opts);
has_state = string(dir_table.State_PaxFWHM) == "paxfwhm_state.mat";
has_fwhm  = string(dir_table.Primary_paxFWHM) == "paxfwhm.mat";
fprintf('%d rows | paxfwhm %d | paxfwhm_state %d\n', ...
    height(dir_table), nnz(has_fwhm), nnz(has_state));

%% 2. Select
select = false(height(dir_table), 1);
switch param.reason
    case "tracelength"
        candidate = find(has_state);
        fprintf('\n%-52s %6s %5s %6s %7s\n', 'session', 'fs', 'tw', 'npts', 'expect');
        for k = candidate'
            state_path = fullfile(char(dir_table.Directory(k)), 'state_analysis', ...
                'paxfwhm_state.mat');
            loaded = load(state_path, 'state_linefwhm');
            obj = loaded.state_linefwhm;
            trans_fields = fieldnames(obj.transition);
            n_point = numel(obj.transition.(trans_fields{1}).data{1});
            window_s = obj.sleep_obj.param.transition_window;
            fs = obj.param.fs;
            n_expect = 2*round(window_s*fs/2) + 1;
            select(k) = n_point ~= n_expect;
            if select(k)
                fprintf('%-52s %6.3f %5g %6d %7d  <- rerun\n', ...
                    extractAfter(string(dir_table.Directory(k)), "tmp\"), ...
                    fs, window_s, n_point, n_expect);
            end
        end

    case "missing"
        select = has_fwhm & ~has_state;

    case "session"
        if param.session_name == ""
            error('rebuild_state_analysis:noName', ...
                'param.reason is "session" but param.session_name is empty.');
        end
        select = contains(string(dir_table.Directory), param.session_name) & has_fwhm;
end

fprintf('\nreason "%s" selected %d of %d sessions\n', param.reason, nnz(select), height(dir_table));
for k = find(select)'
    fprintf('  %s\n', extractAfter(string(dir_table.Directory(k)), "tmp\"));
end
if ~any(select); fprintf('nothing to do\n'); return; end

%% 3. Back up, then run
% note  checkcode marks the else branch UNRCH. That warning is CORRECT and is the
%       point : with param.dryrun at its default nothing below can run. Do not restructure
%       to silence it -- the day it stops appearing is the day the guard is off
if param.dryrun
    fprintf('\nDRY RUN -- nothing written. Set param.dryrun = false to run.\n');
else

n_backed = 0;
for k = find(select)'
    state_path = fullfile(char(dir_table.Directory(k)), 'state_analysis', ...
        'paxfwhm_state.mat');
    if isfile(state_path)
        copyfile(state_path, ...
            char(extractBefore(string(state_path), strlength(state_path)-3) + ...
                 "_backup_" + backup_stamp + ".mat"));
        n_backed = n_backed + 1;
    end
end
fprintf('\nbacked up %d existing state files (stamp %s)\n', n_backed, backup_stamp);

batch_state_analysis(dir_table(select, :), dirs.data_root, primary_map, peripheral_map, ...
    stateanalysis_map, '', true, '');

% 4. Rebuild the mtable so the sheet and the tables agree again
% caution  skipping this leaves mtable_FWHMsleep.mat holding the OLD traces for
%          the sessions just re-run, and every figure downstream reads the mtable
% note     a %% section break cannot go here : this block is inside the param.dryrun
%          else, and running it as a cell on its own would skip the guard
fprintf('\n=== rebuild mtable ===\n');
mtable_FWHMsleep = tableManager(char(sheet_path));
mtable_FWHMsleep.filter_refTable("Primary_paxFWHM", "paxfwhm");
mtable_FWHMsleep.filter_refTable("State_PaxFWHM",   "paxfwhm");
mtable_FWHMsleep.parseParams();
mtable_FWHMsleep.aggregateData("State_PaxFWHM");
mtable_FWHMsleep.addnest2subtable(["state_summary", "transition", ...
    "band_decomposition", "powerdensity", "peak_trough"]);
mtable_FWHMsleep.save2disk("mtable_FWHMsleep.mat");

transition_table = mtable_FWHMsleep.subTables.transition;
fprintf('transition %d rows | trace lengths %s\n', height(transition_table), ...
    mat2str(unique(cellfun(@numel, transition_table.data))'));

end   % param.dryrun
