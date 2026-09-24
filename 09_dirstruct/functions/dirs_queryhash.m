function [folder_hash, file_hash] = dirs_queryhash(folder_path, relative_path)
%DIRS_QUERYHASH  SHA-256 of a folder's files, one per file and one over all of them.
%   Caller: dirs_hashstate
%   The folder hash is the SHA-256 of the sorted "relative path|file hash" lines, so the
%   order the files are listed in cannot move it. A file is opened for reading with others
%   allowed only to read, so a file being written throws instead of hashing half of it.
%
% IN   folder         1 x 1 str     the folder whose files are hashed
%      relative_path  N x 1 str     the files to hash, relative to folder; default every
%                                   file under it, any depth
% OUT  folder_hash    1 x 1 str     64 hex characters
%      file_hash      N x 3 table   RelativePath (str) | Bytes (double) | Hash (str); a file
%                                   that does not exist has Bytes NaN and Hash ""
arguments
    folder_path (1, 1) string
    relative_path (:, 1) string = strings(0, 1)
end
if nargin < 2
    relative_path = listed_files(folder_path);
end

n_file = numel(relative_path);
file_hash = table(relative_path, nan(n_file, 1), strings(n_file, 1), ...
    'VariableNames', {'RelativePath', 'Bytes', 'Hash'});
for k = 1:n_file
    path = fullfile(folder_path, relative_path(k));
    if ~isfile(path)
        continue
    end
    entry = dir(path);
    file_hash.Bytes(k) = entry.bytes;
    file_hash.Hash(k) = sha256_file(path);
end
lines = sort(file_hash.RelativePath + "|" + file_hash.Hash);
folder_hash = sha256_text(strjoin(lines, newline));
end

%% ------------------------------------------------------------------------
function relative_path = listed_files(folder)
% every file under folder, any depth, relative to it
listing = dir(fullfile(folder, '**', '*'));
listing = listing(~[listing.isdir]);
relative_path = strings(numel(listing), 1);
for k = 1:numel(listing)
    full = string(fullfile(listing(k).folder, listing(k).name));
    relative_path(k) = extractAfter(full, strlength(folder) + 1);
end
end

function hex = sha256_file(path)
% hashed inside .NET; MATLAB only opens the stream
sha = System.Security.Cryptography.SHA256.Create();
stream = System.IO.FileStream(path, System.IO.FileMode.Open, System.IO.FileAccess.Read, ...
    System.IO.FileShare.Read, int32(4 * 1048576));
digest = sha.ComputeHash(stream);
stream.Close();
hex = string(strrep(char(System.BitConverter.ToString(digest)), '-', ''));
end

function hex = sha256_text(text)
sha = System.Security.Cryptography.SHA256.Create();
bytes = System.Text.Encoding.UTF8.GetBytes(char(text));
digest = sha.ComputeHash(bytes);
hex = string(strrep(char(System.BitConverter.ToString(digest)), '-', ''));
end
