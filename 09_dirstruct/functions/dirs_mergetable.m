function info = dirs_mergetable(working_root, secondary_root, cohorts, dataset)
%DIRS_MERGETABLE  Scan each cohort, write its sheet, join them into one.

%
%   IN   working_root    1x1 str    the live data tree, one folder per cohort
%        secondary_root  1x1 str    where the sheets go. NOT the data
%        cohorts         1xN str    cohort folder names, in the order they join
%        dataset         1x1 str    names the output folder and the merged sheet
%   OUT  info  struct
%          .sheet_path   1x1 str    the merged sheet this wrote
%          .n_scanned    1xN double sessions the walk found, per cohort
%          .n_offtree    1xN double rows dropped as not under working_root
%          .n_duplicate  1x1 double rows dropped as a repeated Directory
%          .n_stale      1x1 double rows dropped as a dead twin of a live session
%          .collision    1xC str    keys that name more than one LIVE folder, so
%                                   they were left alone. Empty when there are none
%          .n_row        1x1 double rows in the merged sheet

    arguments
        working_root   (1,1) string
        secondary_root (1,1) string
        cohorts        (1,:) string
        dataset        (1,1) string
    end


    primary_map = { ...
        'RadonResult',  'radon_result.mat'; ...
        'RoiList',      'roilist.mat'; ...
        'PolarCluster', 'polarcluster.mat'; ...
        'PaxFWHM',      'paxfwhm.mat'};
    peripheral_map = { ...
        'AnalogAnalysis',   'analysis_analog.mat'; ...
        'BehaviorAnalysis', 'analysis_camera.mat'; ...
        'SleepScore',       'sleep_score.mat'};
    stateanalysis_map = {'paxfwhm_state', 'paxfwhm_state.mat'};

    n_cohort = numel(cohorts);
    n_scanned = zeros(1, n_cohort);
    n_offtree = zeros(1, n_cohort);

    %% 1. Scan each cohort and write its sheet
    % READ ONLY under working_root. The only writes are the per-cohort workbook and
    % its MDF metadata cache, both under secondary_root.
    cohort_sheet = strings(1, n_cohort);
    for cohort_idx = 1:n_cohort
        cohort = cohorts(cohort_idx);
        cohort_dir = fullfile(working_root, cohort);
        cohort_out = fullfile(secondary_root, cohort);
        if ~isfolder(cohort_dir)
            error('dirs_mergetable:noCohort', 'cohort folder not found: %s', cohort_dir);
        end
        if ~isfolder(cohort_out)
            mkdir(cohort_out)
        end
        cohort_sheet(cohort_idx) = fullfile(cohort_out, cohort + "_dirtable.xlsx");
        cache_path = fullfile(cohort_out, cohort + "_dircache.mat");

        scan_start = tic;
        scanned_table = dirs_mapstruct(char(cohort_dir), primary_map, ...
            peripheral_map, stateanalysis_map, char(cache_path));
        dirs_writetable(scanned_table, char(cohort_sheet(cohort_idx)));
        n_scanned(cohort_idx) = height(scanned_table);
        fprintf('%-12s scanned %d sessions in %.1f s\n', cohort, ...
            n_scanned(cohort_idx), toc(scan_start));
    end

    %% 2. Join the cohort sheets

    part = cell(1, n_cohort);
    for cohort_idx = 1:n_cohort
        import_opts = detectImportOptions(cohort_sheet(cohort_idx), 'Sheet', 'reference');
        cohort_table = readtable(cohort_sheet(cohort_idx), import_opts);
        on_working = startsWith(string(cohort_table.Directory), working_root);
        n_offtree(cohort_idx) = nnz(~on_working);
        % which line a row came from, kept for per-cohort counts downstream
        cohort_table.Cohort = repmat(cohorts(cohort_idx), height(cohort_table), 1);
        fprintf('%-12s %4d rows | on %s %4d | dropped %d\n', cohorts(cohort_idx), ...
            height(cohort_table), working_root, nnz(on_working), n_offtree(cohort_idx));
        part{cohort_idx} = cohort_table(on_working, :);
    end
    merged_table = vertcat(part{:});

    % Exact repeats first: the same Directory listed twice
    [~, first_of_directory] = unique(string(merged_table.Directory), 'stable');
    n_duplicate = height(merged_table) - numel(first_of_directory);
    if n_duplicate > 0
        fprintf('dropping %d repeated Directory rows\n', n_duplicate);
        merged_table = merged_table(first_of_directory, :);
    end

    %% 3. One row per session key
    % A folder that gets a suffix (_piv, _notanalyzable, _lowcontrast, ...) is
    % rescanned under its new name while dirs_writetable keeps the row under the old
    % one, so the sheet ends up with two rows for one recording and both carry the
    % same key. Downstream every stage looks a session up by that key with
    % find(...,1) and takes whichever came first, which is a coin toss.
    %   caution  a key that names more than one LIVE folder is NOT a duplicate. It
    %        is two recordings the key cannot tell apart -- HQL072_250325_005_PA01_1
    %        and HQL072_250325_007_PA03_1 both parse to SessionID 1 because the
    %        parser takes the trailing _1. Those are kept and named, because
    %        dropping one loses a recording. see CLAUDE_LOG.md
    row_key = util_sessionkey(merged_table);
    row_alive = false(height(merged_table), 1);
    for k = 1:height(merged_table)
        row_alive(k) = isfolder(string(merged_table.Directory(k)));
    end

    keep = true(height(merged_table), 1);
    collision = strings(1, 0);
    [unique_key, ~, group_of] = unique(row_key, 'stable');
    for group_idx = 1:numel(unique_key)
        member = find(group_of == group_idx);
        if numel(member) < 2
            continue
        end
        live_member = member(row_alive(member));
        if isscalar(live_member)
            keep(member) = false;
            keep(live_member) = true;
        elseif isempty(live_member)
            keep(member) = false;
            keep(member(1)) = true;
        else
            collision(end+1) = unique_key(group_idx); %#ok<AGROW>
        end
    end
    n_stale = sum(~keep);
    if n_stale > 0
        fprintf('dropping %d rows whose folder is gone and whose session key is live elsewhere\n', ...
            n_stale);
    end
    merged_table = merged_table(keep, :);

    for k = 1:numel(collision)
        member = find(util_sessionkey(merged_table) == collision(k));
        fprintf('KEY COLLISION %s -- %d live folders, kept both:\n', ...
            collision(k), numel(member));
        for m = 1:numel(member)
            fprintf('   %s\n', string(merged_table.Directory(member(m))));
        end
    end

    % What centralize_primary will find
    has_fwhm  = string(merged_table.Primary_paxFWHM) == "paxfwhm.mat";
    has_score = string(merged_table.Peripheral_SleepScore) == "sleep_score.mat";
    fprintf('\nmerged %d rows | paxfwhm %d | sleep_score %d | both %d\n', ...
        height(merged_table), nnz(has_fwhm), nnz(has_score), nnz(has_fwhm & has_score));

    %% 4. Write the merged sheet
    % caution  writetable does not clear the sheet first, so a shorter table would
    %          leave the old tail behind. The file has to go -- back it up, do not
    %          just delete, the way fwhm_rederive does
    dataset_dir = fullfile(secondary_root, dataset);
    if ~isfolder(dataset_dir)
        mkdir(dataset_dir)
    end
    sheet_path = fullfile(dataset_dir, dataset + "_dirtable.xlsx");
    if isfile(sheet_path)
        stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        backup_path = fullfile(dataset_dir, dataset + "_dirtable_backup_" + stamp + ".xlsx");
        copyfile(sheet_path, backup_path);
        delete(sheet_path);
        fprintf('backed up previous sheet to %s\n', backup_path);
    end
    writetable(merged_table, sheet_path, 'Sheet', 'reference');
    fprintf('wrote %s\n', sheet_path);

    session_type = string(merged_table.SessionType);
    type_name = unique(session_type);
    for k = 1:numel(type_name)
        fprintf('   %-12s %d\n', type_name(k), sum(session_type == type_name(k)));
    end

    info = struct();
    info.sheet_path = sheet_path;
    info.n_scanned = n_scanned;
    info.n_offtree = n_offtree;
    info.n_duplicate = n_duplicate;
    info.n_stale = n_stale;
    info.collision = collision;
    info.n_row = height(merged_table);
end
