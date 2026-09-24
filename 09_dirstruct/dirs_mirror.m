function info = dirs_mirror(live_root, mirror_root, opt)
%DIRS_MIRROR  Robocopy every live item onto the mirror pool, one item at a time. Adds and updates, never deletes.
%   Caller: archive_session
%   An item is a session folder four levels under a tmp root, or a top-level entry (a
%   date folder or a .7z) under an lto root. Each goes to the mirror root that holds it
%   already, else to the emptiest, so the pool is one union tree, an item that moves
%   between live roots keeps its mirror, and a full mirror disk is extended by adding
%   one. /E only: nothing on a mirror is removed here. What a mirror holds and no live
%   root does is listed, not deleted.
%
% IN   live_root    1 x K str      the tmp and lto folders of the work and archive disks;
%                                  one that does not exist today is skipped
%      mirror_root  1 x M str      the roots of the mirror disks
%      opt.dryrun   1 x 1 logical  robocopy /L, lists only (default false)
%      opt.log_dir  1 x 1 str      robocopy logs and the orphan list (default <first mirror>\mirror_log)
% OUT  info.pass     N x 1 struct  source, target, status, seconds, one per item; robocopy
%                                  status 0-7 is success, 8 and up failure
%      info.log      1 x 1 str     this run's robocopy log
%      info.orphan   N x 1 str     session folders on a mirror and on no live root
arguments
    live_root (1, :) string
    mirror_root (1, :) string
    opt.dryrun (1, 1) logical = false
    opt.log_dir (1, 1) string = ""
end
mirror_root = mirror_root(isfolder(mirror_root));
if isempty(mirror_root)
    error('dirs_mirror:noMirror', 'no mirror root is mounted');
end
if strlength(opt.log_dir) == 0
    opt.log_dir = fullfile(mirror_root(1), "mirror_log");
end
if ~isfolder(opt.log_dir)
    mkdir(opt.log_dir);
end
stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmm'));
info.log = fullfile(opt.log_dir, "mirror_" + stamp + ".log");

%% 1. Every live item onto the mirror root that takes it
% one stream, no /MT: parallel streams onto the HDD mirror measured several times slower, see CLAUDE_LOG.md
flags = '/E /R:2 /W:5 /NP /NFL /NDL /XD $RECYCLE.BIN "System Volume Information"';
if opt.dryrun
    flags = [flags ' /L'];
end
info.pass = struct('source', {}, 'target', {}, 'status', {}, 'seconds', {});
for k = 1:numel(live_root)
    source_root = live_root(k);
    if ~isfolder(source_root)
        fprintf('skip    %s  (not mounted or absent)\n', source_root);
        continue
    end
    [~, leaf] = fileparts(source_root);
    item = live_items(source_root, string(leaf));
    gb_total = sum([item.bytes]) / 2^30;
    fprintf('mirror  %s : %d items, %.1f GB\n', source_root, numel(item), gb_total);
    started = tic;
    gb_done = 0;
    n_failed = 0;
    for j = 1:numel(item)
        relative = fullfile(leaf, item(j).relative);
        target_root = dirs_pickroot(mirror_root, relative);
        target = fullfile(target_root, relative);
        pass = robocopy_item(item(j).source, target, item(j).is_file, flags, info.log);
        info.pass(end + 1) = pass;
        gb_done = gb_done + item(j).bytes / 2^30;
        fprintf('  [%3d/%d] %6.1f GB  %-60s %5.0f s   %.1f of %.1f GB, %.2f h\n', j, numel(item), ...
            item(j).bytes / 2^30, item(j).relative, pass.seconds, gb_done, gb_total, toc(started) / 3600);
        if pass.status >= 8
            n_failed = n_failed + 1;
            fprintf('        FAILED status %d  %s -> %s\n', pass.status, pass.source, pass.target);
        end
    end
    fprintf('        %d failed, %.0f s\n', n_failed, toc(started));
end

%% 2. Orphans: session folders on a mirror and on no live tmp root today
live_tmp = live_root(endsWith(live_root, "tmp") & isfolder(live_root));
info.orphan = strings(0, 1);
for m = 1:numel(mirror_root)
    mirror_tmp = fullfile(mirror_root(m), "tmp");
    if ~isfolder(mirror_tmp)
        continue
    end
    listing = dir(fullfile(mirror_tmp, '*', '*', '*', '*'));
    listing = listing([listing.isdir] & ~startsWith({listing.name}, '.'));
    for k = 1:numel(listing)
        session = fullfile(listing(k).folder, listing(k).name);
        relative = extractAfter(session, strlength(mirror_tmp) + 1);
        on_live = arrayfun(@(root) isfolder(fullfile(root, relative)), live_tmp);
        if ~any(on_live)
            info.orphan(end + 1, 1) = string(session);
        end
    end
end
orphan_path = fullfile(opt.log_dir, "orphan_" + stamp + ".txt");
writelines(info.orphan, orphan_path);
fprintf('orphans %d session folders on a mirror and on no live root -> %s\n', numel(info.orphan), orphan_path);
end

%% ------------------------------------------------------------------------
function item = live_items(source_root, leaf)
% the items of one live root: session folders four levels under tmp, top-level entries under lto
if leaf == "tmp"
    listing = dir(fullfile(source_root, '*', '*', '*', '*'));
    listing = listing([listing.isdir]);
else
    listing = dir(fullfile(source_root, '*'));
end
listing = listing(~startsWith({listing.name}, '.'));
item = struct('source', {}, 'relative', {}, 'is_file', {}, 'bytes', {});
for k = 1:numel(listing)
    source = string(fullfile(listing(k).folder, listing(k).name));
    relative = extractAfter(source, strlength(source_root) + 1);
    bytes = item_bytes(source, listing(k));
    item(end + 1) = struct('source', source, 'relative', relative, 'is_file', ~listing(k).isdir, 'bytes', bytes); %#ok<AGROW>
end
end

function bytes = item_bytes(source, entry)
% the size of an item: the file's own, or every file under the folder
if ~entry.isdir
    bytes = entry.bytes;
    return
end
files = dir(fullfile(source, '**', '*'));
files = files(~[files.isdir]);
bytes = sum([files.bytes]);
end

function pass = robocopy_item(source, target, is_file, flags, log_path)
% one robocopy call; a file goes as <its folder> <target folder> <name>, without /E
started = tic;
if is_file
    [source_dir, name, ext] = fileparts(source);
    target_dir = fileparts(target);
    file_flags = replace(flags, "/E ", "");
    command = sprintf('robocopy "%s" "%s" "%s" %s /LOG+:"%s"', source_dir, target_dir, name + ext, file_flags, log_path);
else
    command = sprintf('robocopy "%s" "%s" %s /LOG+:"%s"', source, target, flags, log_path);
end
[status, ~] = system(command);
pass = struct('source', source, 'target', target, 'status', status, 'seconds', toc(started));
end
