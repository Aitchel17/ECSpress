function info = dirs_hashstate(dirtable_path, live_root, mirror_root, n_scrub)
%DIRS_HASHSTATE  A content hash per session and the mirror's state, as the backup sheet of the dirtable.
%   Caller: archive_session, after dirs_mirror
%   Reads the reference sheet, writes the backup sheet, one row per session, keyed by its path
%   below \tmp\ (the session key beside it is for lookup; it names two folders in places).
%   A session is read again only when its file count, byte count or newest file changed; the
%   rest keep their hash. A session is found where its VolumeSerial says, else under whichever
%   live root holds it today, never by the letter the sheet was written with. Nothing is
%   copied, moved or deleted.
%
% IN   dirtable_path  1 x n char   the merged dirtable workbook, reference sheet in, backup sheet out
%      live_root      1 x K str    the tmp folders mounted today, every work volume and archive
%      mirror_root    1 x M str    the tmp folders of the mirror disks, together the union of the live roots
%      n_scrub        1 x 1 double sessions read back from the mirror this run; NaN = one in
%                                  thirty, 0 = none
% OUT  info.n_session   1 x 1 double  rows in the reference sheet
%      info.n_hashed    1 x 1 double  sessions read this run
%      info.gb_hashed   1 x 1 double  bytes read for that (GB)
%      info.n_mirrored  1 x 1 double  sessions the mirror holds complete, by name and size
%      info.n_scrubbed  1 x 1 double  sessions read back from the mirror
%      info.n_mismatch  1 x 1 double  of those, sessions whose content differed
%      info.seconds     1 x 1 double  wall time

started = tic;
columns = ["SessionKey", "RelativePath", "NumFiles", "Bytes", "GB", "NewestFile", ...
    "SessionHash", "HashedAt", "MirroredAt", "VerifiedAt", "HashMatch", "Mismatch"];

%% 1. Sessions from the reference sheet, and what the backup sheet already knows
reference = readtable(dirtable_path, 'Sheet', 'reference', 'TextType', 'string');
reference.SessionKey = util_sessionkey(reference);
known = table();
if ismember("backup", sheetnames(dirtable_path))
    known = readtable(dirtable_path, 'Sheet', 'backup', 'TextType', 'string');
end
has_known = ~isempty(known) && ismember("RelativePath", string(known.Properties.VariableNames));
if has_known
    % an empty cell comes back from the sheet as <missing>, which compares as nothing and breaks fprintf
    for name = string(known.Properties.VariableNames)
        if isstring(known.(name))
            known.(name)(ismissing(known.(name))) = "";
        end
    end
end

n_session = height(reference);
backup = table('Size', [n_session, numel(columns)], ...
    'VariableTypes', ["string", "string", "double", "double", "double", "string", ...
    "string", "string", "string", "string", "string", "string"], ...
    'VariableNames', columns);
backup.SessionKey = reference.SessionKey;
backup.RelativePath = extractAfter(reference.Directory, '\tmp\');
by_serial = strings(n_session, 1);
if ismember("VolumeSerial", string(reference.Properties.VariableNames))
    by_serial = dirs_sessionpath(reference.VolumeSerial, reference.Directory);
end

%% 2. Each session: size and newest file; the hash again only when one of them moved
n_hashed = 0;
bytes_hashed = 0;
for s = 1:n_session
    old = table();
    if has_known
        old = known(known.RelativePath == backup.RelativePath(s), :);
    end
    session_dir = find_session(by_serial(s), live_root, backup.RelativePath(s));
    if isempty(session_dir)
        % not on any mounted root today: the row keeps what was last known
        if ~isempty(old)
            backup(s, columns(3:end)) = old(1, columns(3:end));
        end
        continue
    end
    listing = session_files(session_dir);
    backup.NumFiles(s) = numel(listing);
    backup.Bytes(s) = sum([listing.bytes]);
    backup.GB(s) = backup.Bytes(s) / 2^30;
    backup.NewestFile(s) = newest_text(listing);

    unchanged = ~isempty(old) ...
        && old.NumFiles(1) == backup.NumFiles(s) ...
        && old.Bytes(1) == backup.Bytes(s) ...
        && old.NewestFile(1) == backup.NewestFile(s);
    if unchanged
        backup(s, columns(7:end)) = old(1, columns(7:end));
    else
        backup.SessionHash(s) = dirs_queryhash(session_dir);
        backup.HashedAt(s) = now_text();
        n_hashed = n_hashed + 1;
        bytes_hashed = bytes_hashed + backup.Bytes(s);
    end

    mirror_dir = find_session("", mirror_root, backup.RelativePath(s));
    if ~isempty(mirror_dir) && mirror_complete(session_dir, listing, mirror_dir)
        backup.MirroredAt(s) = now_text();
    elseif ~isempty(old)
        backup.MirroredAt(s) = old.MirroredAt(1);
    end
end
n_mirrored = nnz(backup.MirroredAt ~= "");
fprintf('sessions %d | hashed %d (%.1f GB) | complete on the mirror %d\n', ...
    n_session, n_hashed, bytes_hashed / 2^30, n_mirrored);

%% 3. Scrub: read the mirror back for the sessions verified longest ago
if isnan(n_scrub)
    n_scrub = ceil(n_session / 30);
end
is_candidate = backup.MirroredAt ~= "";
candidate_idx = find(is_candidate);
[~, order] = sort(backup.VerifiedAt(candidate_idx));     % "" sorts first: never verified
scrub_idx = candidate_idx(order(1:min(n_scrub, numel(order))));
n_mismatch = 0;
for s = scrub_idx'
    session_dir = find_session(by_serial(s), live_root, backup.RelativePath(s));
    mirror_dir = find_session("", mirror_root, backup.RelativePath(s));
    if isempty(mirror_dir)
        backup.MirroredAt(s) = "";
        continue
    end
    [mirror_hash, mirror_files] = dirs_queryhash(mirror_dir);
    backup.VerifiedAt(s) = now_text();
    if mirror_hash == backup.SessionHash(s)
        backup.HashMatch(s) = "true";
        backup.Mismatch(s) = "";
    else
        backup.HashMatch(s) = "false";
        [~, source_files] = dirs_queryhash(session_dir);
        backup.Mismatch(s) = first_difference(source_files, mirror_files);
        n_mismatch = n_mismatch + 1;
        fprintf('MISMATCH %s : %s\n', backup.SessionKey(s), backup.Mismatch(s));
    end
end
fprintf('scrub %d of %d candidates | mismatch %d\n', numel(scrub_idx), numel(candidate_idx), n_mismatch);

%% 4. The backup sheet, this sheet only
writetable(backup, dirtable_path, 'Sheet', 'backup', 'WriteMode', 'overwritesheet');
info = struct('n_session', n_session, 'n_hashed', n_hashed, 'gb_hashed', bytes_hashed / 2^30, ...
    'n_mirrored', n_mirrored, 'n_scrubbed', numel(scrub_idx), 'n_mismatch', n_mismatch, ...
    'seconds', toc(started));
fprintf('backup sheet written to %s in %.0f s\n', dirtable_path, info.seconds);
end

%% ------------------------------------------------------------------------
function folder = find_session(by_serial, live_root, relative_path)
% where the sheet's serial says, else the first live root that holds it today; '' when none does
folder = '';
if strlength(by_serial) > 0 && isfolder(by_serial)
    folder = char(by_serial);
    return
end
for k = 1:numel(live_root)
    candidate = fullfile(char(live_root(k)), char(relative_path));
    if isfolder(candidate)
        folder = candidate;
        return
    end
end
end

function listing = session_files(folder)
% every file under folder, any depth
listing = dir(fullfile(folder, '**', '*'));
listing = listing(~[listing.isdir]);
end

function text = newest_text(listing)
if isempty(listing)
    text = "";
    return
end
newest = datetime(max([listing.datenum]), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
text = string(newest);
end

function text = now_text()
text = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end

function ok = mirror_complete(folder, listing, mirror_dir)
% every file present on the mirror at the same size; content is the scrub's job
ok = ~isempty(listing);
for k = 1:numel(listing)
    full = fullfile(listing(k).folder, listing(k).name);
    mirrored = fullfile(mirror_dir, extractAfter(full, strlength(folder) + 1));
    copy = dir(mirrored);
    if isempty(copy) || copy.bytes ~= listing(k).bytes
        ok = false;
        return
    end
end
end

function text = first_difference(source_files, mirror_files)
% the first file that tells the two apart, from the file tables dirs_queryhash returned
missing = setdiff(source_files.RelativePath, mirror_files.RelativePath);
if ~isempty(missing)
    text = missing(1) + " missing on the mirror";
    return
end
extra = setdiff(mirror_files.RelativePath, source_files.RelativePath);
if ~isempty(extra)
    text = extra(1) + " only on the mirror";
    return
end
[~, source_idx, mirror_idx] = intersect(source_files.RelativePath, mirror_files.RelativePath);
is_different = source_files.Hash(source_idx) ~= mirror_files.Hash(mirror_idx);
differs = source_files.RelativePath(source_idx(is_different));
if isempty(differs)
    text = "no file differs now";
    return
end
text = differs(1) + " differs";
end
