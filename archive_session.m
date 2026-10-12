%ARCHIVE_SESSION  Mirror, hash, then dormant recordings off the work disks onto the archive disks, cell by cell.
% Role:
    % work: NVME, active sessions
    % archive: RAID0, intermediate speed, store mdfextract that unchanged for filt.min_age_month
    % mirror: RAID5, low speed, mirror archieve + work  
% Aim: Empty work drives and ensure integrity based on directory table
% 1. Auto track drives by beacon.txt
% 2. Copy session files (and original .mdf files) with filt criteria      
% 3. Integrity test: Record and compare hash (SHA-256) (merged dirtable)
% 4. Delete if integrity test green
%   How to schedule (in powershell): 
%       schtasks /Create 
%       /TN archive_session 
%       /SC DAILY /MO 3 /ST 03:00

%       /TR "\"C:\Program Files\MATLAB\R2024b\bin\matlab.exe\" -sd g:\03_program\01_ecspress -logfile E:\ecspress_datadriveindex\archive_session_batch.log -batch archive_session"

%% 0. Configuration -- the only cell that decides anything
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;

tmp.dirs = dirs_central();
dirs.work_root    = tmp.dirs.work_root;         % 1 x K str  every disk declaring work, by its beacon
dirs.archive_root = tmp.dirs.archive_root;       % 1 x K str  every disk declaring archive, by its beacon
dirs.mirror_root  = tmp.dirs.mirror_root;        % 1 x M str  every disk declaring mirror
dirs.robocopy_log = fullfile(tempdir, 'archive_session.log');   % char

filt.min_age_month   = 2;           % archive a recording this many months after it was made
filt.product_file    = ["paxfwhm.mat", "polarcluster.mat", "radon_result.mat", "radon_analysis.mat", ...
    "pivtransition.mat", "pivpolar.mat"];   % 1 x n str  a recording holding any of these is kept

param.dataset        = "merged";                % str    the dirtable whose backup sheet says what the mirror verified
param.n_scrub        = NaN;                     % count  sessions read back from the mirror; NaN = one in thirty
dirs.dirtable = fullfile(tmp.dirs.secondary_root, param.dataset, param.dataset + "_dirtable.xlsx");   % derived

% Robocopy:
% eg: robocopy "target dir" "copy dir" /E /R:2 /W:5 /J /ETA /TEE /LOG+:logpath
% /E (entire subdirectory), /R (# retry), /W (# wait between retry sec)
% /J No cache buffer, % /TEE show progress at console /LOG+ add log
param.robocopy_flags = {'/E', '/R:2', '/W:5', '/J', '/ETA', '/TEE'};   % system command "robocopy" arguments

fprintf('work    %s\narchive %s\nmirror  %s\n', strjoin(dirs.work_root, ' | '), ...
    strjoin(dirs.archive_root, ' | '), strjoin(dirs.mirror_root, ' | '));

%% 1. Mirror (robocopy all files under /tmp to mirror flagged drive independently)
% main: dirs_mirror
live_dirs = [dirs.work_root, dirs.archive_root];
dirs_mirror(live_dirs, dirs.mirror_root); % Live directoris copied to mirror root

%% 2. Hash (Look up directory table and add hash for backup sheet, little by little scrub)
% main: dirs_hashstate()
live_tmp = [dirs.work_root, dirs.archive_root];
hash_info = dirs_hashstate(dirs.dirtable, live_tmp, fullfile(dirs.mirror_root, "tmp"),... % 
    param.n_scrub);
fprintf('hashed %d | scrubbed %d | mismatch %d | skipped on error %d | %.0f s\n', ...
    hash_info.n_hashed, hash_info.n_scrubbed, hash_info.n_mismatch, hash_info.n_failed, hash_info.seconds);
if batchStartupOptionUsed
    fprintf('unattended run: mirror and hash done, nothing copied to the archive\n');
    return
end

%% 3. Make latest session file table using dirs and read backup_sheet from directory table
if ~any(isfolder(dirs.archive_root))
    error('archive_session:noArchive', 'no archive root is mounted');
end
record_part = cell(1, numel(dirs.work_root));
for root_idx = 1:numel(dirs.work_root)
    record_part{root_idx} = recordingtable(dirs.work_root(root_idx), filt.product_file);
end
record_table = vertcat(record_part{:});
verified_rel = strings(0, 1);
verified_hash = strings(0, 1);
if isfile(dirs.dirtable) && ismember("backup", sheetnames(dirs.dirtable))
    backup_sheet = readtable(dirs.dirtable, 'Sheet', 'backup', 'TextType', 'string');
    is_match = string(backup_sheet.HashMatch) == "true";
    verified_rel = backup_sheet.RelativePath(is_match);
    verified_hash = string(backup_sheet.SessionHash(is_match));
end
fprintf('%d recordings, %.1f GB | mirror verified %d\n', height(record_table), ...
    sum(record_table.GB), numel(verified_rel));

%% 4. Compare session file table and directory table 
is_old       = record_table.AgeMonth >= filt.min_age_month;
has_product  = record_table.HasProduct;
has_mdf      = record_table.MdfGB > 0;
is_verified  = ismember(lower(record_table.Rel), lower(verified_rel));
keep_archive = is_old & ~has_product & ~has_mdf & is_verified;

fprintf('older than %d months %d | holds a product %d | carries .mdf %d | mirror verified %d\n', ...
    filt.min_age_month, sum(is_old), sum(has_product), sum(has_mdf), sum(is_verified));
fprintf('archive : %d of %d recordings, %.1f GB\n', sum(keep_archive), ...
    height(record_table), sum(record_table.GB(keep_archive)));

move_table = record_table(keep_archive, :);
% the hash the mirror copy was scrubbed against; cell 6 holds the work copy to it
[~, sheet_idx] = ismember(lower(move_table.Rel), lower(verified_rel));
move_table.SheetHash = verified_hash(sheet_idx);
move_table.Target = strings(height(move_table), 1);
for move_idx = 1:height(move_table)
    target_root = dirs_pickroot(dirs.archive_root, move_table.Rel(move_idx));
    if strlength(target_root) == 0
        error('archive_session:archiveGone', 'the archive pool is no longer mounted');
    end
    move_table.Target(move_idx) = fullfile(target_root, move_table.Rel(move_idx));
end
move_table = sortrows(move_table, 'GB', 'descend');
disp(move_table(:, {'GB', 'AgeMonth', 'Rel'}))

%% 5. Copy to the archive   [WRITES under dirs.archive_root]
% Nothing is removed here. Re-running skips what robocopy finds already identical.
% robocopy status 0-7 is success, 8 and up a failure.
move_table.Copied = false(height(move_table), 1);
copy_start = tic;
done_gb = 0;
for move_idx = 1:height(move_table)
    fprintf('[%3d/%d] %6.1f GB  %s   done %.1f of %.1f GB, %.1f h elapsed\n', ...
        move_idx, height(move_table), move_table.GB(move_idx), move_table.Rel(move_idx), ...
        done_gb, sum(move_table.GB), toc(copy_start) / 3600);
    copy_status = copyrecording(move_table.Source(move_idx), move_table.Target(move_idx), ...
        param.robocopy_flags, dirs.robocopy_log);
    move_table.Copied(move_idx) = copy_status >= 0 && copy_status < 8;
    if ~move_table.Copied(move_idx)
        fprintf('COPY FAILED status %d  %s\n', copy_status, move_table.Rel(move_idx));
    end
    done_gb = done_gb + move_table.GB(move_idx);
end
fprintf('copy finished in %.1f h, %d of %d copied\n', toc(copy_start) / 3600, ...
    sum(move_table.Copied), height(move_table));

%% 6. Verify: size and modified time, then one hash of the work and archive copies   [READ ONLY but a re-copy]
% A copy whose sizes or times do not match the source is copied once more before it is hashed.
verify_start = tic;
move_table.Verified = false(height(move_table), 1);
for move_idx = find(move_table.Copied)'
    source_path = move_table.Source(move_idx);
    target_path = move_table.Target(move_idx);
    if ~samefiles(source_path, target_path)
        copy_status = copyrecording(source_path, target_path, param.robocopy_flags, dirs.robocopy_log);
        if copy_status < 0 || copy_status >= 8 || ~samefiles(source_path, target_path)
            fprintf('NOT VERIFIED  sizes or times still differ   %s\n', move_table.Rel(move_idx));
            continue
        end
    end
    source_hash = dirs_queryhash(source_path);
    archive_hash = dirs_queryhash(target_path);
    if archive_hash ~= source_hash
        fprintf('NOT VERIFIED  archive copy content differs  %s\n', move_table.Rel(move_idx));
        continue
    end
    if source_hash ~= move_table.SheetHash(move_idx)
        fprintf('NOT VERIFIED  changed since the mirror scrub %s\n', move_table.Rel(move_idx));
        continue
    end
    move_table.Verified(move_idx) = true;
    fprintf('verified      %6.1f GB  %s\n', move_table.GB(move_idx), move_table.Rel(move_idx));
end
fprintf('verified : %d of %d, %.1f GB, in %.1f h\n', sum(move_table.Verified), height(move_table), ...
    sum(move_table.GB(move_table.Verified)), toc(verify_start) / 3600);

%% 7. Remove the source   [DELETES IN THE DATA TREE]
% Rows cell 6 verified, and only while all three copies still agree on size and modified time.
% A row that moved since cell 6 is skipped; the next run copies and hashes it again.
remove_table = move_table(move_table.Verified, :);
fprintf('about to remove %d recording folders, %.1f GB, under %s\n', ...
    height(remove_table), sum(remove_table.GB), strjoin(dirs.work_root, ' | '));
if ~strcmp(input('type REMOVE to proceed: ', 's'), 'REMOVE')
    fprintf('nothing removed\n');
    return
end
for remove_idx = 1:height(remove_table)
    source_path = remove_table.Source(remove_idx);
    mirror_rel  = fullfile("tmp", remove_table.Rel(remove_idx));
    mirror_disk = dirs_pickroot(dirs.mirror_root, mirror_rel);
    mirror_path = fullfile(mirror_disk, mirror_rel);
    if strlength(mirror_disk) == 0 || ~samefiles(source_path, mirror_path)
        fprintf('SKIP mirror copy missing or out of date   %s\n', remove_table.Rel(remove_idx));
        continue
    end
    if ~samefiles(source_path, remove_table.Target(remove_idx))
        fprintf('SKIP changed since cell 6                 %s\n', remove_table.Rel(remove_idx));
        continue
    end
    rmdir(source_path, 's');
    fprintf('removed %s\n', source_path);
end

%% ------------------------------------------------------------------------
function record_table = recordingtable(work_root, product_file)
    %RECORDINGTABLE  one row per recording folder under <cohort>\<mouse>\<date>\
    %   Caller: archive_session
    %
    % IN   work_root     1 x n char    the tree to walk
    %      product_file  1 x m str     file names that mark a recording as worked on
    % OUT  record_table  n x 6 table   Source Rel GB MdfGB HasProduct AgeMonth; HasProduct is
    %                                  true when any product_file sits anywhere inside;
    %                                  AgeMonth is NaN when the date folder name carries no yyMMdd
    run_date = datetime('today');
    row = struct('Source', {}, 'Rel', {}, 'GB', {}, 'MdfGB', {}, 'HasProduct', {}, 'AgeMonth', {});
    date_folder = dir(fullfile(work_root, '*', '*', '*'));
    date_folder = date_folder([date_folder.isdir] & ~startsWith({date_folder.name}, '.'));
    for date_idx = 1:numel(date_folder)
        age_month = foldermonths(date_folder(date_idx).name, run_date);
        parent_path = fullfile(date_folder(date_idx).folder, date_folder(date_idx).name);
        record_folder = dir(parent_path);
        record_folder = record_folder([record_folder.isdir] & ~startsWith({record_folder.name}, '.'));
        for record_idx = 1:numel(record_folder)
            record_path = fullfile(parent_path, record_folder(record_idx).name);
            file_list = dir(fullfile(record_path, '**', '*'));
            file_list = file_list(~[file_list.isdir]);
            if isempty(file_list)
                continue
            end
            file_name = lower(string({file_list.name}));
            is_mdf = endsWith(file_name, ".mdf");
            row(end+1).Source   = string(record_path); %#ok<AGROW>
            row(end).Rel        = string(extractAfter(record_path, strlength(work_root) + 1));
            row(end).GB         = sum([file_list.bytes]) / 2^30;
            row(end).MdfGB      = sum([file_list(is_mdf).bytes]) / 2^30;
            row(end).HasProduct = any(ismember(file_name, lower(product_file)));
            row(end).AgeMonth   = age_month;
        end
    end
    record_table = struct2table(row);
end

function age_month = foldermonths(folder_name, run_date)
    %FOLDERMONTHS  whole months from the yyMMdd a date folder opens with, to run_date
    %   Caller: recordingtable
    %
    % IN   folder_name  1 x n char        e.g. 250923_hql088_sleep
    %      run_date     1 x 1 datetime
    % OUT  age_month    1 x 1 double      NaN when the name carries no yyMMdd
    token = regexp(folder_name, '^(\d{6})', 'tokens', 'once');
    age_month = NaN;
    if isempty(token)
        return
    end
    made = datetime(token{1}, 'InputFormat', 'yyMMdd');
    age_month = calmonths(between(made, run_date, 'months'));
end

function copy_status = copyrecording(source_path, target_path, flags, log_path)
    %COPYRECORDING  one recording folder by robocopy; status 0-7 is success, 8 and up a failure
    %   Caller: archive_session
    %
    % IN   source_path  1 x 1 string
    %      target_path  1 x 1 string
    %      flags        1 x n cell     robocopy switches
    %      log_path     1 x n char     appended to
    % OUT  copy_status  1 x 1 double   robocopy's exit code; negative when it did not run
    flag_text = strjoin([flags, {['/LOG+:' log_path]}], ' ');
    command = sprintf('robocopy "%s" "%s" %s', source_path, target_path, flag_text);
    copy_status = system(command);
end

function same = samefiles(source_path, target_path)
    %SAMEFILES  every file of source_path present under target_path at the same size and modified time
    %   Caller: archive_session
    %
    % IN   source_path  1 x 1 string
    %      target_path  1 x 1 string
    % OUT  same         1 x 1 logical   false when the target lacks a file or a size or time differs
    source_map = foldermap(source_path);
    target_map = foldermap(target_path);
    if isempty(source_map)
        same = false;
        return
    end
    same = true;
    source_key = keys(source_map);
    for key_idx = 1:numel(source_key)
        if ~isKey(target_map, source_key{key_idx})
            same = false;
            return
        end
        if ~isequal(target_map(source_key{key_idx}), source_map(source_key{key_idx}))
            same = false;
            return
        end
    end
end

function file_map = foldermap(folder_path)
    %FOLDERMAP  every file under folder_path as relative path -> [bytes, modified datenum]
    %   Caller: samefiles
    %
    % IN   folder_path  1 x 1 string
    % OUT  file_map     1 x 1 containers.Map   char key, 1 x 2 double value; empty when absent
    file_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    if ~isfolder(folder_path)
        return
    end
    file_list = dir(fullfile(folder_path, '**', '*'));
    file_list = file_list(~[file_list.isdir]);
    for file_idx = 1:numel(file_list)
        full_path = fullfile(file_list(file_idx).folder, file_list(file_idx).name);
        relative = extractAfter(full_path, strlength(folder_path) + 1);
        file_map(char(relative)) = [file_list(file_idx).bytes, file_list(file_idx).datenum];
    end
end
