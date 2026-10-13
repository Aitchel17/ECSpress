function info = centralize_mat(central_dir, dirtable_path, rebuild)
%CENTRALIZE_MAT  Each session's product .mat files, as one centralized_<product>.mat per product under central_dir.
% Same as directory table with data, byte (key for nextrun), modified (key for nextrun)
%   Caller: d_centralize
%
% IN   central_dir  1 x 1 str    centralized_*.mat out path
%      dirtable_path   1 x 1 str    dirtable path
%      rebuild      1 x 1 bool   true reads every source again
% OUT  info         1 x 1 struct  product -> n_rows, n_read, n_reused, n_dropped, out_path

% Bump when what centralize_strip keeps out of a product changes. Every product is rebuilt in full.
%   2  the air channel (air.on, air.fs) in place of 02_othersignal's airtable
schema_version = 2;



% 1. Define list [file_name folder_name] of file to collect session/folder_name/file_name. mat
products = [ ...               % P x 2 str  <file stem>, <folder under the session>
    "paxfwhm",      "primary_analysis"; ...
    "polarcluster", "primary_analysis"; ...
    "roilist",      "primary_analysis"; ...
    "sleep_score",  "peripheral"; ...
    "analysis_analog", "peripheral"];

% 2. Read directory table and check session directory exist
dir_table = readtable(dirtable_path, 'VariableNamingRule', 'preserve');
alive = false(height(dir_table), 1);
for k = 1:height(dir_table)
    alive(k) = isfolder(string(dir_table.Directory(k)));
end
confirmed_dirtable = dir_table(alive, :);

% 3. Make session key
key_value = util_sessionkey(confirmed_dirtable);
fprintf('dirtable %s\n   %d rows, %d point at a folder\n', dirtable_path, height(dir_table), sum(alive));

% 4. walk through one product at each loop
info = struct();
for p = 1:size(products, 1)
    % 4.1 select filename and subfolder
    product = products(p, 1);
    subfolder = products(p, 2);
    file_name = product + ".mat";
    out_path = fullfile(central_dir, "centralized_" + file_name);
    fprintf('\n=== %s ===\n', product);
    % 4.2 load last centralized product if its schema still holds
    loaded_centraltable = table();
    if isfile(out_path) && ~rebuild
        loaded = load(out_path, 'central', 'schema_version');
        if isfield(loaded, 'schema_version') && loaded.schema_version == schema_version
            loaded_centraltable = loaded.central;
        else
            fprintf('   schema stamp does not match %d, rebuilding in full\n', schema_version);
        end
        clear loaded
    end
    % 4.3 get session key of product table
    key_previous = strings(0, 1);
    if ~isempty(loaded_centraltable)
        key_previous = util_sessionkey(loaded_centraltable);
    end
    % 4.4 struct array generation for centralized struct array
    scan_strarray = struct('confirmed_dirtableidx', {}, 'data', {}, 'source_bytes', {}, ...
        'source_modified', {}, 'reused', {});
    for dirtable_idx = 1:height(confirmed_dirtable)
        session_dir = string(confirmed_dirtable.Directory(dirtable_idx));
        source_path = fullfile(session_dir, subfolder, file_name);
        hit = find(key_previous == key_value(dirtable_idx), 1);
        stay_row = loaded_centraltable(hit, :);
        if ~isfile(source_path)
            continue
        end
        % 4.4.1 Read source file size and modified datetime and make change flag
        source_info = dir(source_path);
        is_reusable = ~isempty(stay_row) && stay_row.source_bytes == source_info.bytes && stay_row.source_modified == source_info.datenum;
        % 4.4.2 Make row with actual source_byte and modified datte
        row = struct('confirmed_dirtableidx', dirtable_idx, 'data', [], 'source_bytes', source_info.bytes, ...
            'source_modified', source_info.datenum, 'reused', is_reusable);
        % 4.4.3 Reuse data from loaded table or read data to reduce read time
        if is_reusable
            row.data = stay_row.data{1}; % reuse
        else
            try
                row.data = read_product(source_path, product); % try to load
            catch err
                fprintf('   FAILED %s -- %s\n', session_dir, err.message);
                if isempty(stay_row)
                    continue
                end
                % if load fail make it same as original
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

    % 4.5 Centralized table generation based on directory table
    existingdata_idx = [scan_strarray.confirmed_dirtableidx];   % sheet rows that got a row
    central = confirmed_dirtable(existingdata_idx, :);   % copy dirtable with extract data existing rows
    central.data = {scan_strarray.data}'; % join data column to dirtable
    central.source_bytes = [scan_strarray.source_bytes]'; % join byte column to dirtable
    central.source_modified = [scan_strarray.source_modified]'; % join modification date column to table

    % 4.6 list duplicated session keys
    key_out = key_value(existingdata_idx);
    if numel(unique(key_out)) < numel(key_out)
        [counted, name] = groupcounts(key_out);
        for k = find(counted > 1)'
            fprintf('   COLLISION %s\n', name(k));
            same_key = existingdata_idx(key_out == name(k));
            fprintf('      %s\n', string(confirmed_dirtable.Directory(same_key)));
        end
    end

    % 4.7 save table
    save(out_path, 'central', 'schema_version', '-v7.3');

    % 4.8 count reused, read and dropped rows and report them
    reused = [scan_strarray.reused];
    dropped = key_previous(~ismember(key_previous, key_out));   % in the last file, no live row today
    source_info = dir(out_path);
    fprintf('   %d rows: %d reused, %d read from source, %d dropped\n', ...
        height(central), nnz(reused), nnz(~reused), numel(dropped));
    fprintf('   %.0f MB -> %s\n', source_info.bytes/1e6, out_path);
    fprintf('   DROPPED %s -- no live session carries this key any more\n', dropped);
    info.(product) = struct('n_rows', height(central), 'n_read', nnz(~reused), ...
        'n_reused', nnz(reused), 'n_dropped', numel(dropped), 'out_path', out_path);
    % 4.9 free the tables before the next product
    clear central scan_strarray loaded_centraltable
end
end

%% ---------------------------------------------------------------- helpers
function data = read_product(source_path, product)
%READ_PRODUCT  One product .mat : its one variable (or the whole file when it holds several), stripped.
%   Caller: centralize_mat
%   An object whose class is not on the path loads as an empty or numeric value; that throws
%   here instead of being loaded.
    % 1. load the session's .mat
    loaded = load(source_path);
    % 2. take its one variable, or the whole struct when it holds several
    field_list = fieldnames(loaded);
    if isscalar(field_list)
        content = loaded.(field_list{1});
    else
        content = loaded;
    end
    % 3. stop on an object whose class is not on the path
    if isnumeric(content) || isempty(content)
        error('centralize_mat:classMissing', ...
            'loaded as %s %s, not an object -- its class is not on the path: %s', ...
            mat2str(size(content)), class(content), source_path);
    end
    % 4. keep what centralize_strip keeps for this product
    data = centralize_strip(product, content);
end
