function info = merge_dirtable(working_root, secondary_root, cohorts, dataset)
%MERGE_DIRTABLE  Scan each cohort, write its sheet, join them into one.
%   Everything that produces a directory table is here: the walk over each cohort,
%   the per-cohort sheet, and the merged sheet centralize_primary reads. It used to
%   go on and build mtable_FWHMsleep.mat by opening every session's
%   paxfwhm_state.mat, and that job moved to centralize_primary and
%   centralize_state, which are incremental.
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
%
%   THE SCAN IS SAFE TO RE-RUN, and that is why centralize_primary calls it before
%   reading. write_dirtable updates only the FILE columns of a row that already
%   exists; VesselID, Depth, Resolution, Notes and the rest of the metadata are
%   left alone, so a value corrected by hand in the reference sheet survives every
%   later scan. It also leaves the workbook untouched when nothing changed. With
%   the MDF metadata cache warm the whole walk is about a second.
%
%   why  a function and not a script : the caller says which cohorts go into the
%        dataset, so that decision lives in ONE place instead of being repeated
%        here and in whatever ran this
%   caution  this is the ONLY place the merged dataset comes from. It lived in a
%        scratch folder until 2026-08-08, which meant the 32-session result could
%        not be rebuilt from the repo
%   note  rename_avi_files is deliberately NOT here. It renames files in the DATA
%        tree, and nothing that runs automatically should do that
    arguments
        working_root   (1,1) string
        secondary_root (1,1) string
        cohorts        (1,:) string
        dataset        (1,1) string
    end


    % What the scan looks for in each session folder, as <column name>, <file name>.
    % Not a setting: it is the list of products this pipeline makes, so it lives
    % with the code that looks for them.
    %   note  paxfwhm_state.mat is still looked for and still written into the
    %         sheet, and nothing reads that column -- the state analysis is
    %         computed into centralized_paxfwhm_state.mat now, not written per
    %         session. The column stays so a session that still has the old file
    %         is recorded rather than silently unlisted
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
            error('merge_dirtable:noCohort', 'cohort folder not found: %s', cohort_dir);
        end
        if ~isfolder(cohort_out)
            mkdir(cohort_out)
        end
        cohort_sheet(cohort_idx) = fullfile(cohort_out, cohort + "_dirtable.xlsx");
        cache_path = fullfile(cohort_out, cohort + "_dircache.mat");

        scan_start = tic;
        scanned_table = mapdirstruct(char(cohort_dir), primary_map, ...
            peripheral_map, stateanalysis_map, char(cache_path));
        write_dirtable(scanned_table, char(cohort_sheet(cohort_idx)));
        n_scanned(cohort_idx) = height(scanned_table);
        fprintf('%-12s scanned %d sessions in %.1f s\n', cohort, ...
            n_scanned(cohort_idx), toc(scan_start));
    end

    %% 2. Join the cohort sheets
    % err  E:\ is a BACKUP drive, not a second data tree. write_dirtable merges an
    %      old sheet with a fresh scan, so an igkltdt sheet that was once written
    %      while pointed at E:\ keeps those rows, and they duplicate sessions that
    %      also sit under G:\tmp. Every row not on working_root is dropped
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
    % rescanned under its new name while write_dirtable keeps the row under the old
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
