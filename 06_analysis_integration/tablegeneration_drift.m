%TABLEGENERATION_DRIFT  drift_statetable : one row per (session x state x bout) of tissue drift.
%   Reads centralized_drift.mat (the motion-correction series, sign flipped to tissue
%   motion by io_readmotion) and centralized_sleep_score.mat, puts each session's
%   scored bouts on the drift sample axis, and writes drift_statetable.mat beside the
%   dataset. Every quantity is a difference; no absolute position leaves this file.
%
%   ROW = (session, state, bout)
%     lvl_row / lvl_col   bout median minus the baseline : the quiet samples of the
%                         flanks, param.flank_s either side, same-state bouts excluded
%     net_row / net_col   median of the last edge_frac of the bout minus the first
%     speed_med           median per-sample step, px/s
%     net_mag / net_angle hypot and atan2d of net_row, net_col. Average the components,
%                         never the angle
%     flank_ok            both flanks lie inside the recording
%
%   AXES  +col = column index up = image right = lateral; +row = row index up = image
%         down = caudal. rostral is -row and is applied by the figure layer, not here.
%     caution  the OUTPUT lands in the dataset folder
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

%% settings
param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
param.state_names = ["awake", "nrem", "rem", "drowsy", "uarousal"];   % time_table fields read (str)
param.flank_s     = 180;      % Flank either side of a bout, the baseline is read there (s)
param.quiet_pct   = 15;       % Flank samples kept for the baseline, by step size, both axes (pct)
param.min_quiet   = 10;       % Quiet samples a baseline needs, else NaN (count)
param.edge_frac   = 0.10;     % Share of the bout at each end that net_* compares (frac)
param.control_state = "nrem"; % The state a figure reads REM against (str)

filt.session_type = "sleep";  % Session type kept (str), "" = every type
filt.drop_backup  = true;     % Drop a row whose folder says backup (logical)

dirs = dirs_central();
dirs.out = fullfile(dirs.secondary_root, param.dataset);

%% Read the two centralized products
in_path = fullfile(dirs.central, "centralized_" + ["drift", "sleep_score"] + ".mat");
for k = 1:numel(in_path)
    if ~isfile(in_path(k))
        error('tablegeneration_drift:noCentral', 'run centralize_primary / centralize_drift first; %s is missing', in_path(k));
    end
end
drift_table = load(in_path(1)).central;
score_table = load(in_path(2)).central;
fprintf('centralized drift %d rows, sleep_score %d rows\n', height(drift_table), height(score_table));

%% Which sessions the table is about
cond.right_type  = filt.session_type == "" | string(drift_table.SessionType) == filt.session_type;
cond.not_backup  = ~filt.drop_backup | ~contains(string(drift_table.Directory), "backup", 'IgnoreCase', true);
key_drift = util_sessionkey(drift_table);
key_score = util_sessionkey(score_table);
cond.has_score   = ismember(key_drift, key_score);
keep_session = cond.right_type & cond.not_backup & cond.has_score;
fprintf('drift_statetable : %d of %d sessions kept (type %d | not backup %d | scored %d)\n', ...
    sum(keep_session), numel(keep_session), sum(cond.right_type), sum(cond.not_backup), sum(cond.has_score));
filtered = drift_table(keep_session, :);
key_filtered = key_drift(keep_session);

%% One row per bout, every session
gathered = cell(height(filtered), 1);
for k = 1:height(filtered)
    series = filtered.data{k};
    score = score_table.data{key_score == key_filtered(k)};
    % the drift sample axis : seconds from the recording's first frame, the clock
    % the score's time tables are on
    t_drift = (0:series.drift_n - 1) / series.drift_fps;
    state_integrate = state_integration(score, []);
    evalc('state_integrate.trim_to_duration(t_drift(end));');   % it prints three lines
    img_state = state_image(state_integrate);
    evalc('img_state.get_state_indices(t_drift, series.drift_fps);');   % add_taxis prints every field
    gathered{k} = drift_bouts(series.disp_row, series.disp_col, series.regerror, ...
        series.drift_fps, img_state.state_idx, param);
    if mod(k, 25) == 0
        fprintf('   %d/%d sessions\n', k, height(filtered));
    end
end
n_row_session = cellfun(@numel, gathered);

%% The table
key_names = ["MouseID", "Date", "SessionType", "SessionID", "VesselID", "Depth", "Cohort", "Directory"];
key_names = intersect(key_names, string(filtered.Properties.VariableNames), 'stable');
session_of = repelem((1:height(filtered))', n_row_session);
keys = filtered(session_of, key_names);
keys.vessel_key = string(keys.MouseID) + "_" + string(keys.VesselID);
bout_rows = struct2table([gathered{:}]');
drift_statetable = [keys, bout_rows];
fprintf('drift_statetable : %d rows over %d sessions, %d vessels\n', height(drift_statetable), ...
    height(filtered), numel(unique(drift_statetable.vessel_key)));
for state_name = param.state_names
    on_state = drift_statetable.state_name == state_name;
    fprintf('   %-9s %4d bouts | flank_ok %4d | baseline %4d\n', state_name, sum(on_state), ...
        sum(on_state & drift_statetable.flank_ok), sum(on_state & isfinite(drift_statetable.lvl_row)));
end

%% Save
%   caution  this overwrites, and the copy on disk is the only record of which param built it
if ~isfolder(dirs.out)
    mkdir(dirs.out)
end
save_content = struct('drift_statetable', drift_statetable, 'param', param, 'filt', filt);
out_path = fullfile(dirs.out, 'drift_statetable.mat');
save(out_path, "save_content")
listing = dir(out_path);
fprintf('%.1f MB -> %s\n', listing.bytes / 1e6, out_path);

%% ---------------------------------------------------------------- helpers
function rows = drift_bouts(disp_row, disp_col, regerror, fps, state_idx, param)
%DRIFT_BOUTS  One session's rows : every bout of every state in param.state_names.
%   IN   disp_row, disp_col  N x 1 double  tissue displacement, px
%        regerror            N x 1 double  registration residual
%        fps                 1 x 1 double  drift samples per second
%        state_idx           struct        <state> = B x 2 sample indices, from state_image
%        param               struct        flank_s, quiet_pct, min_quiet, edge_frac, state_names
%   OUT  rows                1 x R struct  the columns of drift_statetable after the keys
    n_sample = numel(disp_row);
    flank_n = round(param.flank_s * fps);
    step_row = [NaN; abs(diff(disp_row(:)))];
    step_col = [NaN; abs(diff(disp_col(:)))];
    rows = struct('state_name', {}, 'bout_idx', {}, 'n_bout', {}, 'bout_start_s', {}, ...
        'bout_end_s', {}, 'bout_duration_s', {}, 'n_sample', {}, 'n_quiet', {}, 'flank_ok', {}, ...
        'lvl_row', {}, 'lvl_col', {}, 'net_row', {}, 'net_col', {}, 'net_mag', {}, ...
        'net_angle', {}, 'speed_med', {}, 'regerror_med', {}, 'regerror_p95', {}, 'drift_fps', {});
    for state_name = param.state_names
        if ~isfield(state_idx, state_name)
            continue
        end
        bouts = state_idx.(state_name);
        n_bout = size(bouts, 1);
        % every sample inside any bout of THIS state, so a flank never reads one
        in_state = false(n_sample, 1);
        for b = 1:n_bout
            in_state(bouts(b, 1):bouts(b, 2)) = true;
        end
        for b = 1:n_bout
            first = bouts(b, 1);
            last = bouts(b, 2);
            inside = first:last;
            % 1. the flanks, minus the same state
            before = max(1, first - flank_n) : first - 1;
            after = last + 1 : min(n_sample, last + flank_n);
            flank = [before, after];
            flank = flank(~in_state(flank));
            % 2. the quiet samples of the flanks : small steps on both axes
            quiet = false(size(flank));
            if ~isempty(flank)
                cut_row = prctile(step_row(flank), param.quiet_pct);
                cut_col = prctile(step_col(flank), param.quiet_pct);
                quiet = step_row(flank) <= cut_row & step_col(flank) <= cut_col;
            end
            baseline_row = NaN;
            baseline_col = NaN;
            if nnz(quiet) >= param.min_quiet
                baseline_row = median(disp_row(flank(quiet)));
                baseline_col = median(disp_col(flank(quiet)));
            end
            % 3. the two edges of the bout
            n_edge = max(1, round(param.edge_frac * numel(inside)));
            head = inside(1:n_edge);
            tail = inside(end - n_edge + 1:end);
            net_row = median(disp_row(tail)) - median(disp_row(head));
            net_col = median(disp_col(tail)) - median(disp_col(head));
            step_inside = hypot(diff(disp_row(inside)), diff(disp_col(inside)));
            rows(end+1) = struct( ...
                'state_name',      state_name, ...
                'bout_idx',        b, ...
                'n_bout',          n_bout, ...
                'bout_start_s',    (first - 1) / fps, ...
                'bout_end_s',      (last - 1) / fps, ...
                'bout_duration_s', (last - first) / fps, ...
                'n_sample',        numel(inside), ...
                'n_quiet',         nnz(quiet), ...
                'flank_ok',        first - flank_n >= 1 && last + flank_n <= n_sample, ...
                'lvl_row',         median(disp_row(inside)) - baseline_row, ...
                'lvl_col',         median(disp_col(inside)) - baseline_col, ...
                'net_row',         net_row, ...
                'net_col',         net_col, ...
                'net_mag',         hypot(net_row, net_col), ...
                'net_angle',       atan2d(net_row, net_col), ...
                'speed_med',       median(step_inside) * fps, ...
                'regerror_med',    median(regerror(inside)), ...
                'regerror_p95',    prctile(regerror(inside), 95), ...
                'drift_fps',       fps); %#ok<AGROW>
        end
    end
end
