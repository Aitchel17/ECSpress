function info = bout_drift(central_dir, rebuild)
%BOUT_DRIFT  combine centralized_drift and centralized_bout using state1d
%   Caller: d_centralize
%
% IN   central_dir  1 x 1 str    where the centralized_*.mat files live
%      rebuild      1 x 1 bool   true recomputes every session
% OUT  info         1 x 1 struct  n_rows, n_computed, n_reused, out_path

% Bump when what is taken off state_1d changes, or when either method changes what it
% produces. The inputs are untouched by such a change, so no stamp can see it.
%   2  the bout tables are centralized_bout's and are no longer copied in
%   3  state windows are half-open, and only the transition windows lose the end margin;
%      param carries state_1d's settings (butter, psd) beside fs; input_stamp is split into
%      input_source_modified and input_bout_computed; a transition window cut by the record's first or
%      last frame is dropped instead of clipped (bout_trim is gone)
schema_version = 3;

% 1. check the two centralized tables exist
in_path = fullfile(central_dir, "centralized_" + ["drift", "bout"] + ".mat");
[in_exists, ~, ~] = util_checkdirstat(in_path);
if ~in_exists(1)
    error('bout_drift:noDrift', 'run centralize_drift first, centralized_drift.mat is missing');
end
if ~in_exists(2)
    error('bout_drift:noBout', 'run centralize_bout first, centralized_bout.mat is missing');
end
% 2. load them
drift_table = load(in_path(1)).central;
bout_table = load(in_path(2)).central;
fprintf('drift %d rows, bout %d\n', height(drift_table), height(bout_table));
% 3. session key of each table
[key_drift, key_columns] = util_sessionkey(drift_table);
key_bout = util_sessionkey(bout_table);

% 4. load the previous output, when its schema matches
out_path = fullfile(central_dir, "centralized_drift_state.mat");
loaded_table = table();
if isfile(out_path) && ~rebuild
    stored = load(out_path, 'central', 'schema_version');
    if isfield(stored, 'schema_version') && stored.schema_version == schema_version
        loaded_table = stored.central;
    else
        fprintf('schema stamp does not match %d, recomputing every session\n', schema_version);
    end
    clear stored
end
% 4.1 session key of the loaded rows
key_loaded = strings(0, 1);
if ~isempty(loaded_table)
    key_loaded = util_sessionkey(loaded_table);
end

% 5. one row per session, split again or reused
scan = struct('drift_i', {}, 'data', {}, 'input_source_modified', {}, 'input_bout_computed', {}, 'reused', {});
for k = 1:height(drift_table)
    % 5.1 matching bout row
    hit_boutidx = find(key_bout == key_drift(k), 1);
    if isempty(hit_boutidx)
        fprintf('   no bout tables %s\n', key_drift(k));
        continue
    end
    % 5.2 its bout tables
    bout_data = bout_table.data{hit_boutidx};
    % 5.3 the two input stamps
    source_modified = drift_table.source_modified(k);
    bout_computed = bout_table.computed(hit_boutidx);
    % 5.4 this session's row in the loaded output, if any
    hit = find(key_loaded == key_drift(k), 1);
    stay_row = loaded_table(hit, :);
    % 5.5 is the loaded row still valid
    is_reusable = false;   % reuse only when neither the drift source moved nor the bout row was recomputed
    if ~isempty(stay_row)
        is_samesource = stay_row.input_source_modified == source_modified;
        is_samebout = stay_row.input_bout_computed == bout_computed;
        is_reusable = is_samesource && is_samebout;
    end
    % 5.6 the row: the loaded result, a new split, or on failure the loaded row and its stamps
    row = struct('drift_i', k, 'data', [], 'input_source_modified', source_modified, ...
        'input_bout_computed', bout_computed, 'reused', is_reusable);

    if is_reusable
        row.data = stay_row.data{1}; % reuse
    else
        drift = drift_table.data{k};
        try
            row.data = get_driftstate(drift, bout_data.time_table, bout_data.transition_window);
        catch err
            fprintf('   FAILED %s -- %s\n', key_drift(k), err.message);
            if isempty(stay_row)
                continue
            end
            % a session that will not split keeps the loaded row, and its stamp
            row.data = stay_row.data{1};
            row.input_source_modified = stay_row.input_source_modified;
            row.input_bout_computed = stay_row.input_bout_computed;
            row.reused = true;
        end
    end
    % 5.7 add to the output struct array
    scan(end+1) = row; %#ok<AGROW>
    if mod(k, 10) == 0
        fprintf('   %d/%d scanned, %d held\n', k, height(drift_table), numel(scan));
    end
end

% 6. make output
% 6.1 rows kept
kept = [scan.drift_i];
% 6.2 key columns, the columns the filters read, Depth and Resolution as numbers (stored as '70um', '0.19um'), and the frame rate
central = drift_table(kept, [key_columns, "VesselID", "Depth", "Resolution", "Cohort", "Directory"]);
central.NumericDepth = util_unit2double(drift_table.Depth(kept));
central.NumericResolution = util_unit2double(drift_table.Resolution(kept));
central.drift_fps = drift_table.drift_fps(kept);

% 6.3 data column and the input stamps
central.data = {scan.data}';
central.input_source_modified = [scan.input_source_modified]';   % the stamps the next run reuses this row against
central.input_bout_computed = [scan.input_bout_computed]';
save(out_path, 'central', 'schema_version', '-v7.3');

% 6.4 summary
listing = dir(out_path);
reused = [scan.reused];
fprintf('\n%d rows: %d computed, %d reused\n', ...
    height(central), nnz(~reused), nnz(reused));
fprintf('%.0f MB -> %s\n', listing.bytes/1e6, out_path);
info = struct('n_rows', height(central), 'n_computed', nnz(~reused), 'n_reused', nnz(reused), 'out_path', out_path);
end

%% ---------------------------------------------------------------- helpers
function data = get_driftstate(drift, time_table, transition_window)
%GET_DRIFTSTATE  One session's two drift components split by its bouts through state_1d: the results and the axes.
%   Caller: bout_drift
%   IN   drift              1x1 struct     its centralized_drift data: disp_row, disp_col (N x 1 px), drift_fps
%        time_table         1x1 struct     its bout tables, centralized_bout's data.time_table
%        transition_window  1x1 float s    centralized_bout's data.transition_window, 0 for an analog row
%   OUT  data               1x1 struct     state_summary, transition, t_axis, state_idx, param
    % 1. hold the bouts on the drift's time axis and map them onto its frames
    t_drift = (0:numel(drift.disp_row) - 1) / drift.drift_fps;
    drift_state = state_1d(time_table, transition_window, t_drift);
    drift_state.get_state_indices();
    % 2. flatten each component by one straight line over the session, 1 x N like the FWHM series
    row_flat = detrend(drift.disp_row)';
    col_flat = detrend(drift.disp_col)';
    % 3. run each through summary and transition
    drift_state.get_summary("drift_row", row_flat);
    drift_state.get_transitionsummary("drift_row", row_flat);
    drift_state.get_summary("drift_col", col_flat);
    drift_state.get_transitionsummary("drift_col", col_flat);
    % 4. keep the results and the axes
    data = struct();
    for f = ["state_summary", "transition", "t_axis", "state_idx", "param"]
        data.(f) = drift_state.(f);
    end
end
