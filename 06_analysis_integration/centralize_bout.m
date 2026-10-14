function info = centralize_bout(central_dir, rebuild)
%CENTRALIZE_BOUT  Read centralized sleep_score and analysis_analog to make fusion bouttable
% whiskerb,p,v session require analysis_analog table, sleep session require sleep_score

% IN   central_dir  1 x 1 str    where the centralized_*.mat files live; centralized_sleep_score
%                                and centralized_analysis_analog are read from it
%      rebuild      1 x 1 bool   true recomputes every session
% OUT  info         1 x 1 struct  n_rows, n_computed, n_reused, out_path

% Bump when bout_sleep or bout_puffmove changes what a time_table holds or how a bout is cut.
% The inputs are untouched by such a change, so no timestamp can see it; every row is recomputed,
% computed moves, and bout_paxfwhm and bout_drift follow it.
%   2  puff_quietend_sec : the pre-puff quiet window ends 0.1 s before the onset
%   3  speed rule out, the puff itself is the go window, long chunks out, puffs off the air channel
%   4  data.transition_window on every row, 0 for an analog row; param.dur<X> without rounded_duration;
%      loco_rest_puff / _nopuff are loco_restpuff / loco_restnopuff
schema_version = 4;


% 1. load centralized table of sleep_score and analysis_anlog
in_path = fullfile(central_dir, "centralized_" + ["sleep_score", "analysis_analog"] + ".mat");
[in_exists, ~, ~] = util_checkdirstat(in_path);
if ~all(in_exists)
    error('centralize_bout:noScoring', ...
        'run centralize_mat first, centralized_sleep_score.mat or centralized_analysis_analog.mat is missing');
end
sleepscore_table = load(in_path(1)).central;
analog_table = load(in_path(2)).central;
fprintf('sleep_score %d rows, analysis_analog %d\n', height(sleepscore_table), height(analog_table));

% 2. Extract session key of sleepscore table and analog table
[key_sleepscore, key_columns] = util_sessionkey(sleepscore_table); % get key name together
key_analog = util_sessionkey(analog_table);

% 3. join sleepscore and analog table to make fusion table
session_table = [sleepscore_table(:, key_columns); analog_table(:, key_columns)]; % sleep session must come first for analog purge purpose
% 3.2 remove analysis analog row of sleep session from session table
key_session = util_sessionkey(session_table);
[key_session, first_idx] = unique(key_session, 'stable'); % if sleepscore.mat exist, analysis analog row automatically dropped
session_table = session_table(first_idx, :);

% 4. load centralized_bout.mat
out_path = fullfile(central_dir, "centralized_bout.mat");
loaded_bouttable = table();
if isfile(out_path) && ~rebuild
    stored = load(out_path, 'central', 'schema_version');
    if isfield(stored, 'schema_version') && stored.schema_version == schema_version
        loaded_bouttable = stored.central;
    else
        fprintf('schema timestamp does not match %d, recomputing every session\n', schema_version);
    end
    clear stored
end

% 5. Get session key from loaded bout table
key_previous = strings(0, 1);
if ~isempty(loaded_bouttable)
    key_previous = util_sessionkey(loaded_bouttable);
end

% 6. Fill struct array using loaded_bouttable or from sleepscore.mat or newly compute from centralized data
scan_strarray = struct('session_idx', {}, 'bout_category', {}, 'data', {}, 'input_stamp', {}, ...
    'computed', {}, 'reused', {});
for session_idx = 1:height(session_table)
    % 6.1 make timestamp and read sleepscore or analog primary
    sleepscore = [];
    primary_analog = [];
    timestamp = []; % time stamp of centralized time
    hit_sleepscore = find(key_sleepscore == key_session(session_idx), 1); 
    hit_analog = find(key_analog == key_session(session_idx), 1);
    if ~isempty(hit_sleepscore)
        sleepscore = sleepscore_table.data{hit_sleepscore}; 
        timestamp = sleepscore_table.source_modified(hit_sleepscore); 
    elseif ~isempty(hit_analog) && isfield(analog_table.data{hit_analog}, 'air')
        primary_analog = analog_table.data{hit_analog};
        timestamp = analog_table.source_modified(hit_analog);
    end
    if isempty(timestamp)
        fprintf('   no bout_category %s\n', key_session(session_idx));
        continue
    end
    % 6.2 check whether row already inserted or not
    hit = find(key_previous == key_session(session_idx), 1);
    stay_row = loaded_bouttable(hit, :);
    % 6.3 check centralized data changed or not
    is_reusable = ~isempty(stay_row) && stay_row.input_stamp == timestamp;

    % 6.4 empty row
    row = struct('session_idx', session_idx, 'bout_category', "", 'data', [], 'input_stamp', timestamp, ...
        'computed', convertTo(datetime('now'), 'datenum'), 'reused', is_reusable);   % the clock source_modified is on
    % 6.5 fill the row
    if is_reusable % reuse stay_row if matched row exist, and row computed from matched timestamp of file
        row.bout_category = string(stay_row.bout_category);
        row.data = stay_row.data{1}; % reuse
        row.computed = stay_row.computed;
    else
        try
            [row.data, row.bout_category] = get_bouttable(sleepscore, primary_analog); % calculate 
        catch err
            % if row preexist leave with empty row, if exist fill with existing
            fprintf('   FAILED %s -- %s\n', key_session(session_idx), err.message);
            if isempty(stay_row)
                continue
            end
            row.bout_category = string(stay_row.bout_category);
            row.data = stay_row.data{1};
            row.input_stamp = stay_row.input_stamp;
            row.computed = stay_row.computed;
            row.reused = true;
        end
    end
    scan_strarray(end+1) = row; %#ok<AGROW>
    if mod(session_idx, 25) == 0
        fprintf('   %d/%d scanned, %d held\n', session_idx, height(session_table), numel(scan_strarray));
    end
end

%% out
existingdata_idx = [scan_strarray.session_idx];
reused = [scan_strarray.reused];
central = session_table(existingdata_idx, :);
central.bout_category = [scan_strarray.bout_category]';
central.data = {scan_strarray.data}';
central.input_stamp = [scan_strarray.input_stamp]';
central.computed = [scan_strarray.computed]';
save(out_path, 'central', 'schema_version', '-v7.3');

listing = dir(out_path);
fprintf('\n%d rows: %d computed, %d reused\n', height(central), nnz(~reused), nnz(reused));
fprintf('%.0f MB -> %s\n', listing.bytes/1e6, out_path);
info = struct('n_rows', height(central), 'n_computed', nnz(~reused), 'n_reused', nnz(reused), 'out_path', out_path);
end

%% ---------------------------------------------------------------- helpers
function [data, bout_category] = get_bouttable(sleepscore, primary_analog)
%GET_BOUTTABLE  One session's bout tables, from its sleep sleepscore or else its air channel.
%   Caller: centralize_bout
%   IN   sleepscore           1x1 struct | []   the sleep score, if it has one
%        primary_analog  1x1 struct | []   the whisker bout_category instead
%   OUT  data            1x1 struct        time_table, param, transition_window (s, 0 for an analog row)
%        bout_category         1x1 str           "sleep" or "analog"
    if ~isempty(sleepscore)
        bout_category = "sleep";
        [time_table, param] = bout_sleep(sleepscore.behavState, sleepscore.statecodes, sleepscore.binwidth_sec);
        transition_window = param.transition_window;
    else
        bout_category = "analog";
        transition_window = 0;   % s  a puff row has no *_trans table, so no window
        [time_table, param] = bout_puffmove(primary_analog.air.on, primary_analog.air.fs, ...
            primary_analog.ball.resampled_ball, primary_analog.ball.rs_taxis, primary_analog.ball.ds_fps);
    end
    data = struct('time_table', time_table, 'param', param, 'transition_window', transition_window);
end
