%CENTRALIZE_DRIFT  Every session's *_motion.txt as one table: dataset sheet -> motion files -> centralized_drift.mat.
%   Reads the sheet dirtable_main writes; run that first when a recording is new.
%   Why this is not a centralize_primary product: see CLAUDE_LOG.md

%% settings
param.dataset  = "merged";      % str   the dataset whose sheet lists the sessions
param.rebuild  = false;         % bool  true re-reads every session
schema_version = 1;             % double  bump when what is taken off the file changes

addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
dirs = dirs_central();

%% the sessions
sheet_path = fullfile(dirs.secondary_root, param.dataset, param.dataset + "_dirtable.xlsx");
dir_table = readtable(sheet_path, 'VariableNamingRule', 'preserve');

alive = false(height(dir_table), 1);
for k = 1:height(dir_table)
    alive(k) = isfolder(string(dir_table.Directory(k)));
end
live_table = dir_table(alive, :);
key_value = util_sessionkey(live_table);
fprintf('dirtable %s\n   %d rows, %d point at a folder\n', ...
    sheet_path, height(dir_table), sum(alive));

%% what is already read
out_path = fullfile(dirs.central, "centralized_drift.mat");
previous = table();
if isfile(out_path) && ~param.rebuild
    stored = load(out_path, 'central', 'schema_version');
    if isfield(stored, 'schema_version') && stored.schema_version == schema_version
        previous = stored.central;
    else
        fprintf('schema stamp does not match %d, reading every session again\n', schema_version);
    end
    clear stored
end
key_previous = strings(0, 1);
if ~isempty(previous)
    key_previous = util_sessionkey(previous);
end

%% one row per live session that carries a motion file
scan = struct('live_i', {}, 'data', {}, 'source_bytes', {}, ...
    'source_modified', {}, 'reused', {});
for i = 1:height(live_table)
    hit = find(key_previous == key_value(i), 1);
    old = previous(hit, :);
    [row, failure] = scan_session(i, string(live_table.Directory(i)), old);
    if failure ~= ""
        fprintf('   FAILED %s -- %s\n', string(live_table.Directory(i)), failure);
    end
    if ~isempty(row)
        scan(end+1) = row; %#ok<SAGROW>
    end
    if mod(i, 25) == 0
        fprintf('   %d/%d scanned, %d held\n', i, height(live_table), numel(scan));
    end
end

%% out
kept = [scan.live_i];
reused = [scan.reused];
central = live_table(kept, :);
central.data = {scan.data}';
central.source_bytes = [scan.source_bytes]';
central.source_modified = [scan.source_modified]';

% the scalars as flat columns, so a filt names one without opening data
summary = [scan.data];
central.drift_fps = [summary.drift_fps]';
central.drift_n = [summary.drift_n]';
for name = ["disp_row_p95", "disp_row_max", "disp_row_stepmax", ...
            "disp_col_p95", "disp_col_max", "disp_col_stepmax", ...
            "regerror_p50", "regerror_p95"]
    central.(name) = [summary.(name)]';
end

key_out = key_value(kept);
if numel(unique(key_out)) < numel(key_out)
    [counted, name] = groupcounts(key_out);
    for k = find(counted > 1)'
        fprintf('   COLLISION %s\n', name(k));
    end
end

save(out_path, 'central', 'schema_version', '-v7.3');
listing = dir(out_path);
fprintf('\n%d rows: %d read, %d reused\n', height(central), nnz(~reused), nnz(reused));
fprintf('%.1f MB -> %s\n', listing.bytes/1e6, out_path);

%% ---------------------------------------------------------------- helpers
function [row, failure] = scan_session(i, session_dir, old)
%SCAN_SESSION  One session's row, or 0x0 if it carries no motion file.
%   IN   i            1 x 1 double     which row of live_table this is
%        session_dir  1 x 1 str        the session folder
%        old          0|1 row table    what the last run wrote for this session
%   OUT  row          0x0 | 1x1 struct
%        failure      1 x 1 str        "" unless the read threw
    row = struct('live_i', {}, 'data', {}, 'source_bytes', {}, ...
        'source_modified', {}, 'reused', {});
    failure = "";

    % globbed: the stem is the MDF file's name, not the session's; see CLAUDE_LOG.md
    found = dir(fullfile(session_dir, '*_motion.txt'));
    if isempty(found)
        return
    end
    if numel(found) > 1
        failure = sprintf('%d motion files, took %s', numel(found), found(1).name);
    end
    listing = found(1);
    source_path = fullfile(listing.folder, listing.name);

    unmoved = ~isempty(old) && old.source_bytes == listing.bytes ...
        && old.source_modified == listing.datenum;
    if ~unmoved
        try
            [disp_row, disp_col, regerror, drift_fps] = io_readmotion(source_path);
            row(1).live_i = i;
            row(1).data = drift_summary(disp_row, disp_col, regerror, drift_fps);
            row(1).source_bytes = listing.bytes;
            row(1).source_modified = listing.datenum;
            row(1).reused = false;
            return
        catch err
            % a file that will not parse keeps its old row; the failure is named
            failure = string(err.message);
        end
    end
    if isempty(old)
        return
    end
    row(1).live_i = i;
    row(1).data = old.data{1};
    row(1).source_bytes = old.source_bytes;
    row(1).source_modified = old.source_modified;
    row(1).reused = true;
end

function summary = drift_summary(disp_row, disp_col, regerror, drift_fps)
%DRIFT_SUMMARY  The series, and the scalars a filt reads without opening them.
%   Scalars are taken on magnitudes, so a sign flip in the file moves none of them.
%   OUT  summary  1x1 struct   series as columns, then the scalars
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
