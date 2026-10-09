%CENTRALIZE_DRIFTSTATE  State analysis over the two drift components, every session,
%   written as one table. The drift counterpart of centralize_paxfwhmstate: the
%   same class, the same two methods, the same output shape, so tablegeneration_main
%   and the figure scripts read drift_row / drift_col the way they read a thickness.
%
%   Run centralize_drift first. This reads what that wrote.

%% settings
param.key_columns  = ["MouseID", "Date", "SessionType", "SessionID"];
param.rebuild      = false;         % bool  true recomputes every session

% Bump when what is taken off state_1d changes, or when either method changes what
% it produces. The inputs are untouched by such a change, so no stamp can see it.
schema_version = 1;

addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
dirs = dirs_central();

%% the inputs
in_path = fullfile(dirs.central, "centralized_" + ["drift", "sleep_score", "analysis_analog"] + ".mat");
[in_exists, ~, ~] = util_checkdirstat(in_path);
if ~in_exists(1)
    error('centralize_driftstate:noDrift', 'run centralize_drift first, centralized_drift.mat is missing');
end
drift_table = load(in_path(1)).central;
score_table = table();
analog_table = table();
if in_exists(2)
    score_table = load(in_path(2)).central;
end
if in_exists(3)
    analog_table = load(in_path(3)).central;
end
fprintf('drift %d rows, sleep_score %d, analysis_analog %d\n', ...
    height(drift_table), height(score_table), height(analog_table));

key_drift = util_sessionkey(drift_table);
key_score = util_sessionkey(score_table);
key_analog = util_sessionkey(analog_table);

%% what is already computed
out_path = fullfile(dirs.central, "centralized_drift_state.mat");
previous = table();
if isfile(out_path) && ~param.rebuild
    stored = load(out_path, 'central', 'schema_version');
    if isfield(stored, 'schema_version') && stored.schema_version == schema_version
        previous = stored.central;
    else
        fprintf('schema stamp does not match %d, recomputing every session\n', schema_version);
    end
    clear stored
end
key_previous = strings(0, 1);
if ~isempty(previous)
    key_previous = util_sessionkey(previous);
end

%% one row per session that has a scoring to split by
scan = struct('drift_i', {}, 'data', {}, 'input_stamp', {}, 'reused', {});
for k = 1:height(drift_table)
    [score, primary_analog, stamp] = centralize_scoring(score_table, key_score, ...
        analog_table, key_analog, key_drift(k), drift_table.source_modified(k));
    if isempty(stamp)
        fprintf('   no scoring %s\n', key_drift(k));
        continue
    end
    hit = find(key_previous == key_drift(k), 1);
    old = previous(hit, :);
    series = drift_table.data{k};
    [row, failure] = centralize_driftsession(k, series.disp_row, series.disp_col, ...
        series.drift_fps, score, primary_analog, stamp, old);
    if failure ~= ""
        fprintf('   FAILED %s -- %s\n', key_drift(k), failure);
    end
    if ~isempty(row)
        scan(end+1) = row; %#ok<SAGROW>
    end
    if mod(k, 10) == 0
        fprintf('   %d/%d scanned, %d held\n', k, height(drift_table), numel(scan));
    end
end

%% out
kept = [scan.drift_i];
reused = [scan.reused];
central = drift_table(kept, param.key_columns);
extra = intersect(["VesselID", "Depth", "Resolution", "Cohort", "Directory"], ...
    string(drift_table.Properties.VariableNames));
central = [central, drift_table(kept, extra)];

% Depth reads '70um' and Resolution '0.19um'. Both are read as numbers downstream,
% Depth by the layer filter and Resolution by the pixel-to-micrometre scaling.
central.NumericDepth = util_unit2double(drift_table.Depth(kept));
central.NumericResolution = util_unit2double(drift_table.Resolution(kept));

central.data = {scan.data}';
central.input_stamp = [scan.input_stamp]';
central.drift_fps = drift_table.drift_fps(kept);
save(out_path, 'central', 'schema_version', '-v7.3');

listing = dir(out_path);
fprintf('\n%d rows: %d computed, %d reused\n', ...
    height(central), nnz(~reused), nnz(reused));
fprintf('%.0f MB -> %s\n', listing.bytes/1e6, out_path);
