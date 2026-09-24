%ARCHIVE_SESSION  Mirror, hash, then dormant recordings off the work disks onto the archive disks, cell by cell.
%   Scheduled: schtasks /Create /TN archive_session /SC DAILY /ST 03:00
%              /TR "\"C:\Program Files\MATLAB\R2024b\bin\matlab.exe\" -batch archive_session"
%   Unattended it mirrors, hashes, copies and verifies; the removal cell wants a typed
%   word and returns without it, so a daily run leaves the data tree as it was.
%   The work volumes are the disks every stage reads from. What clogs them is not the
%   finished work but the recordings that were extracted and then never taken up:
%   a recording centralize_primary has not gathered months after it was made is not
%   being worked on, and its TIFFs are what makes it large.
%
%   A recording folder is archived when all four hold:
%     it was made at least filt.min_age_month months ago
%     centralized_<filt.central_product>.mat does not list it
%     no .mdf sits inside it -- the .mdf go to the queue drive, not here
%     the backup sheet of the dirtable says the mirror holds it byte-identical, so a
%     verified copy exists somewhere before the work copy goes
%
%   Cell 1 writes to the mirror disks, cell 2 the backup sheet of the dirtable, cell 5
%   under the archive disks. Cells 0, 3, 4 and 6 touch nothing. CELL 7 IS THE ONLY ONE
%   THAT DELETES, it deletes in the DATA TREE, and it deletes only what cell 6 verified
%   in the same MATLAB session.
%
%   The date comes from the date folder above the recording, so a recording under
%   <mouse>\zstack\ has none and is never archived. dirs_mergetable has to read the
%   archive root as well, or a moved recording drops off the sheet at the next
%   centralize run.

%% 0. Configuration -- the only cell that decides anything
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;

resolved = dirs_central();
dirs.central      = resolved.central;            % char   where centralized_*.mat live
dirs.work_root    = resolved.work_root;          % 1 x K str  every disk declaring work, by its beacon
dirs.archive_root = resolved.archive_root;       % 1 x K str  every disk declaring archive, by its beacon
dirs.lto_root     = resolved.lto_root;           % 1 x K str  the lto folders of those
dirs.mirror_root  = resolved.mirror_root;        % 1 x M str  every disk declaring mirror
dirs.robocopy_log = fullfile(tempdir, 'archive_session.log');   % char

filt.min_age_month   = 2;           % archive a recording this many months after it was made
filt.central_product = "paxfwhm";   % a recording listed in centralized_<this>.mat is kept

param.dataset        = "merged_igkl_igkltdt";   % str    the dirtable whose backup sheet says what the mirror verified
param.n_scrub        = NaN;                     % count  sessions read back from the mirror; NaN = one in thirty
param.robocopy_flags = {'/E', '/R:2', '/W:5', '/ETA', '/TEE'};   % 1 x n cell
dirs.dirtable = fullfile(resolved.secondary_root, param.dataset, param.dataset + "_dirtable.xlsx");   % derived

if ~any(isfolder(dirs.archive_root))
    error('archive_session:noArchive', 'no archive root is mounted');
end
fprintf('work    %s\narchive %s\nmirror  %s\n', strjoin(dirs.work_root, ' | '), ...
    strjoin(dirs.archive_root, ' | '), strjoin(dirs.mirror_root, ' | '));

%% 1. Mirror: every live item onto the mirror disks   [WRITES under dirs.mirror_root]
live_root = [dirs.work_root, dirs.archive_root, dirs.lto_root];
mirror_info = dirs_mirror(live_root, dirs.mirror_root);

%% 2. Hash: a content hash per session and the mirror's state, into the backup sheet   [WRITES dirs.dirtable]
live_tmp = [dirs.work_root, dirs.archive_root];
hash_info = dirs_hashstate(dirs.dirtable, live_tmp, fullfile(dirs.mirror_root, "tmp"), param.n_scrub);
fprintf('mirror %d items | hashed %d | scrubbed %d | mismatch %d\n', numel(mirror_info.pass), ...
    hash_info.n_hashed, hash_info.n_scrubbed, hash_info.n_mismatch);

%% 3. Every recording on the work disks, what centralize gathered, what the mirror holds verified   [READ ONLY]
record_part = cell(1, numel(dirs.work_root));
for root_idx = 1:numel(dirs.work_root)
    record_part{root_idx} = recordingtable(dirs.work_root(root_idx));
end
record_table = vertcat(record_part{:});
central_path = fullfile(dirs.central, "centralized_" + filt.central_product + ".mat");
stored = load(central_path, 'central');
gathered_dir = string(stored.central.Directory);
verified_rel = strings(0, 1);
if isfile(dirs.dirtable) && ismember("backup", sheetnames(dirs.dirtable))
    backup_sheet = readtable(dirs.dirtable, 'Sheet', 'backup', 'TextType', 'string');
    is_match = string(backup_sheet.HashMatch) == "true";
    verified_rel = backup_sheet.RelativePath(is_match);
end
fprintf('%d recordings, %.1f GB | %s lists %d | mirror verified %d\n', height(record_table), ...
    sum(record_table.GB), filt.central_product, numel(gathered_dir), numel(verified_rel));

%% 4. Which recordings archive   [READ ONLY]
is_old       = record_table.AgeMonth >= filt.min_age_month;
is_gathered  = ismember(lower(record_table.Source), lower(gathered_dir));
has_mdf      = record_table.MdfGB > 0;
is_verified  = ismember(lower(record_table.Rel), lower(verified_rel));
keep_archive = is_old & ~is_gathered & ~has_mdf & is_verified;

fprintf('older than %d months %d | gathered %d | carries .mdf %d | mirror verified %d\n', ...
    filt.min_age_month, sum(is_old), sum(is_gathered), sum(has_mdf), sum(is_verified));
fprintf('archive : %d of %d recordings, %.1f GB\n', sum(keep_archive), ...
    height(record_table), sum(record_table.GB(keep_archive)));

move_table = record_table(keep_archive, :);
move_table.Target = strings(height(move_table), 1);
for move_idx = 1:height(move_table)
    target_root = dirs_pickroot(dirs.archive_root, move_table.Rel(move_idx));
    move_table.Target(move_idx) = fullfile(target_root, move_table.Rel(move_idx));
end
move_table = sortrows(move_table, 'GB', 'descend');
disp(move_table(:, {'GB', 'AgeMonth', 'Rel'}))

%% 5. Copy to the archive   [WRITES under dirs.archive_root]
% Nothing is removed here. Re-running skips what robocopy finds already identical.
copy_start = tic;
done_gb = 0;
for move_idx = 1:height(move_table)
    fprintf('[%3d/%d] %6.1f GB  %s   done %.1f of %.1f GB, %.1f h elapsed\n', ...
        move_idx, height(move_table), move_table.GB(move_idx), move_table.Rel(move_idx), ...
        done_gb, sum(move_table.GB), toc(copy_start) / 3600);
    flags = strjoin([param.robocopy_flags, {['/LOG+:' dirs.robocopy_log]}], ' ');
    command = sprintf('robocopy "%s" "%s" %s', move_table.Source(move_idx), ...
        move_table.Target(move_idx), flags);
    system(command);
    done_gb = done_gb + move_table.GB(move_idx);
end
fprintf('copy finished in %.1f h\n', toc(copy_start) / 3600);

%% 6. Verify every file, name and size   [READ ONLY]
move_table.Verified = false(height(move_table), 1);
for move_idx = 1:height(move_table)
    move_table.Verified(move_idx) = samefiles(move_table.Source(move_idx), ...
        move_table.Target(move_idx));
    fprintf('%-62s %6.1f GB  %s\n', move_table.Rel(move_idx), move_table.GB(move_idx), ...
        string(move_table.Verified(move_idx)));
end
fprintf('verified : %d of %d, %.1f GB\n', sum(move_table.Verified), height(move_table), ...
    sum(move_table.GB(move_table.Verified)));

%% 7. Remove the source   [DELETES IN THE DATA TREE]
% Only rows cell 6 verified, and each is checked again at the moment it is removed.
remove_table = move_table(move_table.Verified, :);
fprintf('about to remove %d recording folders, %.1f GB, under %s\n', ...
    height(remove_table), sum(remove_table.GB), strjoin(dirs.work_root, ' | '));
if batchStartupOptionUsed
    fprintf('unattended run: nothing removed\n');
    return
end
if ~strcmp(input('type REMOVE to proceed: ', 's'), 'REMOVE')
    fprintf('nothing removed\n');
    return
end
for remove_idx = 1:height(remove_table)
    if ~samefiles(remove_table.Source(remove_idx), remove_table.Target(remove_idx))
        fprintf('SKIP no longer verifies  %s\n', remove_table.Rel(remove_idx));
        continue
    end
    rmdir(remove_table.Source(remove_idx), 's');
    fprintf('removed %s\n', remove_table.Source(remove_idx));
end

%% ------------------------------------------------------------------------
function record_table = recordingtable(work_root)
    %RECORDINGTABLE  one row per recording folder under <cohort>\<mouse>\<date>\
    %   Caller: archive_session
    %
    % IN   work_root     1 x n char    the tree to walk
    % OUT  record_table  n x 5 table   Source Rel GB MdfGB AgeMonth; AgeMonth is NaN
    %                                  when the date folder name carries no yyMMdd
    run_date = datetime('today');
    row = struct('Source', {}, 'Rel', {}, 'GB', {}, 'MdfGB', {}, 'AgeMonth', {});
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
            is_mdf = endsWith(string({file_list.name}), ".mdf", 'IgnoreCase', true);
            row(end+1).Source = string(record_path); %#ok<AGROW>
            row(end).Rel      = string(extractAfter(record_path, strlength(work_root) + 1));
            row(end).GB       = sum([file_list.bytes]) / 2^30;
            row(end).MdfGB    = sum([file_list(is_mdf).bytes]) / 2^30;
            row(end).AgeMonth = age_month;
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

function same = samefiles(source_path, target_path)
    %SAMEFILES  every file of source_path present under target_path at the same size
    %   Caller: archive_session
    %
    % IN   source_path  1 x 1 string
    %      target_path  1 x 1 string
    % OUT  same         1 x 1 logical   false when the target lacks a file or a size differs
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
        if target_map(source_key{key_idx}) ~= source_map(source_key{key_idx})
            same = false;
            return
        end
    end
end

function file_map = foldermap(folder_path)
    %FOLDERMAP  every file under folder_path as relative path -> bytes
    %   Caller: samefiles
    %
    % IN   folder_path  1 x 1 string
    % OUT  file_map     1 x 1 containers.Map   char key, double value; empty when absent
    file_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    if ~isfolder(folder_path)
        return
    end
    file_list = dir(fullfile(folder_path, '**', '*'));
    file_list = file_list(~[file_list.isdir]);
    for file_idx = 1:numel(file_list)
        full_path = fullfile(file_list(file_idx).folder, file_list(file_idx).name);
        relative = extractAfter(full_path, strlength(folder_path) + 1);
        file_map(char(relative)) = file_list(file_idx).bytes;
    end
end
