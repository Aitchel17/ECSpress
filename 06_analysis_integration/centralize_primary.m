% Centralize files in each session file
%% Directory setup
clc; clear;
addpath('g:\03_program\01_ecspress\09_dirstruct'); dirs_ecspath; % add all function in path                                     
dirs = dirs_central; % centralize path setup
if ~isfolder(dirs.central)
    mkdir(dirs.central);
    fprintf('created %s\n', dirs.central);
end
%% refresh the directory table
% cohorts is the only place the dataset's membership is decided -- dirs_mergetable
% takes it as an argument rather than keeping a second copy
dataset     = "merged_igkl_igkltdt";
cohorts     = ["00_igkl", "01_igkltdt"];
sheet_info = dirs_mergetable(dirs.working_root, dirs.secondary_root, ...
    cohorts, dataset);
%% read directory table and extract keys
dir_table = readtable(sheet_info.sheet_path, 'VariableNamingRule', 'preserve');

alive = false(height(dir_table), 1);
for k = 1:height(dir_table)
    alive(k) = isfolder(string(dir_table.Directory(k)));
end
live_table = dir_table(alive, :);
key_value = util_sessionkey(live_table);

fprintf('dirtable %s\n   %d rows\n', sheet_info.sheet_path, height(dir_table));
fprintf('   %d rows point at a folder, %d do not\n', sum(alive), sum(~alive));


%% Centralize each analysis output

products    = [ ...            % Px2 str  <file stem>, <folder under the session>
    "paxfwhm",      "primary_analysis"; ...
    "polarcluster", "primary_analysis"; ...
    "roilist",      "primary_analysis"; ...
    "sleep_score",  "peripheral"; ...
    "analysis_analog", "peripheral"];
param.rebuild = false;         % bool     true reads every source again
schema_version = 1;            % double   bumping it forces a full rebuild

%% one product at a time, so peak memory is one table and not three
for p = 1:size(products, 1)
    product = products(p, 1);
    subfolder = products(p, 2);
    file_name = product + ".mat";
    out_path = fullfile(dirs.central, "centralized_" + file_name);
    fprintf('\n=== %s ===\n', product);
    %% Load session key if previous struct exist
    previous = table();
    if isfile(out_path) && ~param.rebuild % if rebuild skip loading
        stored = load(out_path, 'central', 'schema_version');
        if isfield(stored, 'schema_version') && stored.schema_version == schema_version
            previous = stored.central;
        else
            fprintf('   schema stamp does not match %d, rebuilding in full\n', schema_version);
        end
        clear stored
    end
    key_previous = strings(0, 1);
    if ~isempty(previous)
        key_previous = util_sessionkey(previous);
    end
    %% one row per live session that carries this product
    scan = struct('live_i', {}, 'data', {}, 'source_bytes', {}, ...
        'source_modified', {}, 'reused', {});
    for i = 1:height(live_table)
        source_path = fullfile(string(live_table.Directory(i)), subfolder, file_name);
        hit = find(key_previous == key_value(i), 1);
        old = previous(hit, :);
        [row, failure] = centralize_scansession(i, source_path, product, old);
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

    % the descriptive columns are taken fresh from the sheet every run; only data
    % and its stamp are ever carried over
    kept = [scan.live_i];
    reused = [scan.reused];
    central = live_table(kept, :);
    central.data = {scan.data}';
    central.source_bytes = [scan.source_bytes]';
    central.source_modified = [scan.source_modified]';

    % the key has to separate the rows that are actually in this table
    key_out = key_value(kept);
    if numel(unique(key_out)) < numel(key_out)
        [counted, name] = groupcounts(key_out);
        for k = find(counted > 1)'
            fprintf('   COLLISION %s\n', name(k));
            same_key = kept(key_out == name(k));
            fprintf('      %s\n', string(live_table.Directory(same_key)));
        end
    end

    % Writing 2.4 GB back unchanged is the whole cost of a run where nothing moved,
    % so the save is skipped when the new table matches the old one everywhere the
    % data is not: same rows, same order, same descriptive columns, nothing
    % recomputed and nothing dropped. Comparing without the data column keeps that
    % test cheap.
    dropped = key_previous(~ismember(key_previous, key_out));
    unchanged = all(reused) && isempty(dropped) && ~isempty(previous) ...
        && height(previous) == height(central) ...
        && isequal(previous(:, ~strcmp(previous.Properties.VariableNames, 'data')), ...
                   central(:, ~strcmp(central.Properties.VariableNames, 'data')));
    if unchanged
        fprintf('   %d rows, all reused and identical -- left as it is\n', height(central));
    else
        save(out_path, 'central', 'schema_version', '-v7.3');
    end

    listing = dir(out_path);
    fprintf('   %d rows: %d reused, %d read from source, %d dropped\n', ...
        height(central), nnz(reused), nnz(~reused), numel(dropped));
    fprintf('   %.0f MB -> %s\n', listing.bytes/1e6, out_path);
    fprintf('   DROPPED %s -- no live session carries this key any more\n', dropped);
    clear central scan previous
end

