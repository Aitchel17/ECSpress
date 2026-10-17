function info = centralize_drift(central_dir, dirtable_path, rebuild)
%CENTRALIZE_DRIFT  Each session's *_motion.txt, as one centralized_drift.mat under central_dir.
%   Caller: d_centralize
%   A session is read again only when its motion file moved since the last run; the table is saved every run.
%
% IN   central_dir    1 x 1 str  where the centralized_*.mat files live
%      dirtable_path  1 x 1 str  the dataset's dirtable workbook, as c_dirtable writes it; its first sheet lists the sessions
%      rebuild        1 x 1 bool true re-reads every session
% OUT  info           1 x 1 struct  n_rows, n_read, n_reused, out_path

% Bump when what is taken off the file changes; the stamp cannot see a code change
schema_version = 1;

% 1. read directory table and verify existance
dir_table = readtable(dirtable_path, 'VariableNamingRule', 'preserve');
alive = false(height(dir_table), 1);
for k = 1:height(dir_table)
    alive(k) = isfolder(string(dir_table.Directory(k)));
end
confirmed_dirtable = dir_table(alive, :);
key_value = util_sessionkey(confirmed_dirtable);
fprintf('dirtable %s\n   %d rows, %d point at a folder\n', ...
    dirtable_path, height(dir_table), sum(alive));

% 2. load centralized_drift table if exist
out_path = fullfile(central_dir, "centralized_drift.mat");
loaded_centraltable = table();
if isfile(out_path) && ~rebuild
    loaded = load(out_path, 'central', 'schema_version');
    if isfield(loaded, 'schema_version') && loaded.schema_version == schema_version
        loaded_centraltable = loaded.central;
    else
        fprintf('schema stamp does not match %d, reading every session again\n', schema_version);
    end
    clear loaded
end
% 3. get session key from loaded centralized table
key_previous = strings(0, 1);
if ~isempty(loaded_centraltable)
    key_previous = util_sessionkey(loaded_centraltable);
end

% 4. fill structarray
scan_strarray =struct('confirmed_dirtableidx', {}, 'data', {}, 'source_bytes', {}, 'source_modified', {}, 'reused', {});
for dirtable_idx = 1:height(confirmed_dirtable)
    session_dir = string(confirmed_dirtable.Directory(dirtable_idx));
    hit = find(key_previous == key_value(dirtable_idx), 1);
    stay_row = loaded_centraltable(hit, :);
    % 4.1 find the motion file; globbed, since the stem is the MDF file's name, not the session's
    found = dir(fullfile(session_dir, '*_motion.txt'));
    if isempty(found)
        continue
    end
    if numel(found) > 1
        fprintf('   %d motion files in %s, took %s\n', numel(found), session_dir, found(1).name);
    end
    source_info = found(1);
    source_path = fullfile(source_info.folder, source_info.name);
    % 4.2 compare its size and modified time with the loaded row
    is_reusable = ~isempty(stay_row) && stay_row.source_bytes == source_info.bytes && stay_row.source_modified == source_info.datenum;
    row = struct('confirmed_dirtableidx', dirtable_idx, 'data', [], 'source_bytes', source_info.bytes, ...
        'source_modified', source_info.datenum, 'reused', is_reusable);
    % 4.3 reuse the loaded data, or read the file again
    if is_reusable
        row.data = stay_row.data{1};
    else
        try
            row.data = get_driftsummary(source_path);
        catch err
            fprintf('   FAILED %s -- %s\n', session_dir, err.message);
            if isempty(stay_row)
                continue
            end
            % a file that will not parse keeps the loaded row, and its stamp
            row.data = stay_row.data{1};
            row.source_bytes = stay_row.source_bytes;
            row.source_modified = stay_row.source_modified;
            row.reused = true;
        end
    end
    scan_strarray(end+1) = row; %#ok<AGROW>
    if mod(dirtable_idx, 25) == 0
        fprintf('   %d/%d scanned, %d held\n', dirtable_idx, height(confirmed_dirtable), numel(scan_strarray));
    end
end

%% out
existingdata_idx = [scan_strarray.confirmed_dirtableidx];
central = confirmed_dirtable(existingdata_idx, :);
central.data = {scan_strarray.data}';
central.source_bytes = [scan_strarray.source_bytes]';       % the stamps bout_drift rebuilds its rows against
central.source_modified = [scan_strarray.source_modified]';

% the scalars as flat columns, so a filt names one without opening data
summary = [scan_strarray.data];
central.drift_fps = [summary.drift_fps]';
central.drift_n = [summary.drift_n]';
for name = ["disp_row_p95", "disp_row_max", "disp_row_stepmax", ...
            "disp_col_p95", "disp_col_max", "disp_col_stepmax", ...
            "regerror_p50", "regerror_p95"]
    central.(name) = [summary.(name)]';
end

key_out = key_value(existingdata_idx);
if numel(unique(key_out)) < numel(key_out)
    [counted, name] = groupcounts(key_out);
    for k = find(counted > 1)'
        fprintf('   COLLISION %s\n', name(k));
    end
end

save(out_path, 'central', 'schema_version', '-v7.3');
source_info = dir(out_path);
reused = [scan_strarray.reused];

fprintf('\n%d rows: %d read, %d reused\n', height(central), nnz(~reused), nnz(reused));
fprintf('%.1f MB -> %s\n', source_info.bytes/1e6, out_path);
info = struct('n_rows', height(central), 'n_read', nnz(~reused), 'n_reused', nnz(reused), 'out_path', out_path);
end

%% ---------------------------------------------------------------- helpers
function summary = get_driftsummary(source_path)
%GET_DRIFTSUMMARY  One *_motion.txt read : the series, and the scalars a filt reads without opening them.
%   Caller: centralize_drift
%   Scalars are taken on magnitudes, so a sign flip in the file moves none of them.
%   OUT  summary  1x1 struct   series as columns, then the scalars
    [disp_row, disp_col, regerror, drift_fps] = io_readmotion(source_path); % sign inverted in here (correction -> actual movement)
    summary = struct( ...
        'disp_row',         disp_row, ...      % N x 1 double  px
        'disp_col',         disp_col, ...      % N x 1 double  px
        'regerror',         regerror, ...      % N x 1 double  dimensionless
        'drift_fps',        drift_fps, ...     % 1 x 1 double  Hz
        'drift_n',          numel(disp_row), ...
        'disp_row_p95',     prctile(abs(disp_row), 95), ...
        'disp_row_max',     max(abs(disp_row)), ...
        'disp_row_stepmax', max(abs(diff(disp_row))), ...
        'disp_col_p95',     prctile(abs(disp_col), 95), ...
        'disp_col_max',     max(abs(disp_col)), ...
        'disp_col_stepmax', max(abs(diff(disp_col))), ...
        'regerror_p50',     median(regerror), ...
        'regerror_p95',     prctile(regerror, 95));
end
