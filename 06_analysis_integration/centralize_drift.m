%CENTRALIZE_DRIFT  Every session's *_motion.txt as one table, series and all.
%   The sixth centralized product. Reads the acquisition's drift estimate off the
%   session folders and writes centralized_drift.mat beside the other five.
%
%   IT CANNOT JOIN centralize_primary'S PRODUCT LOOP, for three reasons that are all
%   in that file: :46 welds ".mat" onto the stem, :68 builds a fixed stem inside a
%   fixed subfolder, and centralize_scansession :22 reads with load() and :29 throws
%   on anything that is not an object. The motion file is text, it sits at the
%   session root, and its stem is per session.
%
%   THE FILENAME IS NOT THE SESSION. mdf_xymovie.m:216 builds the path from the MDF
%   FILE's name rather than the folder's, so the stem is found by globbing and the
%   row's identity comes from the dirtable. One folder under hql073 holds a file
%   named for HQL063; keying off the filename would file that session under the
%   wrong mouse.
%
%   Run centralize_primary first only if the dirtable needs refreshing. This reads
%   the sheet, not the centralized products.

%% settings
param.dataset  = "merged_igkl_igkltdt";
param.rebuild  = false;         % bool  true re-reads every session

% which sessions the product is about
cohorts = ["00_igkl", "01_igkltdt"];   % 1 x C str  membership, as centralize_primary

% Bump when what is taken off the file changes. The sources are untouched by such a
% change, so no stamp can see it.
schema_version = 1;

addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
dirs = dirs_central();

%% the sessions
sheet_info = dirs_mergetable(dirs.working_root, dirs.secondary_root, ...
    cohorts, param.dataset);
dir_table = readtable(sheet_info.sheet_path, 'VariableNamingRule', 'preserve');

alive = false(height(dir_table), 1);
for k = 1:height(dir_table)
    alive(k) = isfolder(string(dir_table.Directory(k)));
end
live_table = dir_table(alive, :);
key_value = util_sessionkey(live_table);
fprintf('dirtable %s\n   %d rows, %d point at a folder\n', ...
    sheet_info.sheet_path, height(dir_table), sum(alive));

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

% the scalars every consumer reads without opening a series. Flat columns rather
% than another level inside data, so a filt can name one directly
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
            % a file that will not parse is no reason to drop what was already read
            % from it, so the old row stands and the failure is named
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
%   Every scalar is taken on a MAGNITUDE, so none of them moves if the stored sign
%   convention is ever reversed.
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
