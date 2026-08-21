function [exists, bytes, modified] = util_checkdirstat(paths, opt)
%UTIL_CHECKDIRSTAT  What is on disk at each of these paths, printed and returned.
%   A stage whose input is missing fails deep inside a loop. Run this before the
%   loop and it says so first.
%
% IN   paths     1xN string or cellstr, files or folders
%      label     1xN string shown in place of the path, default the last two
%                components: one list often holds the same file name under
%                several folders, and the name alone cannot tell them apart
%      quiet     true prints nothing, default false
% OUT  exists    1xN logical
%      bytes     1xN double; total over the contents for a folder, NaN if absent
%      modified  1xN datetime; NaT if absent
%
%   The paths are not echoed back: the caller already has them, and an output
%   that repeats its own argument cannot be told apart from a measurement.

arguments
    paths       (1,:) string
    opt.label   (1,:) string = strings(1,0)
    opt.quiet   (1,1) logical = false
end

n_path = numel(paths);
exists = false(1, n_path);
bytes = nan(1, n_path);
modified = NaT(1, n_path);

shown = opt.label;
if isempty(shown)
    shown = strings(1, n_path);
    for k = 1:n_path
        [parent, stem, ext] = fileparts(paths(k));
        [~, parent_stem, parent_ext] = fileparts(parent);
        shown(k) = fullfile(parent_stem + parent_ext, stem + ext);
    end
end

for k = 1:n_path
    if isfolder(paths(k))
        exists(k) = true;
        listing = dir(fullfile(paths(k), '**', '*'));
        listing = listing(~[listing.isdir]);
        bytes(k) = sum([listing.bytes]);
        here = dir(paths(k));
        modified(k) = datetime(here(1).datenum, 'ConvertFrom', 'datenum');
    elseif isfile(paths(k))
        exists(k) = true;
        here = dir(paths(k));
        bytes(k) = here.bytes;
        modified(k) = datetime(here.datenum, 'ConvertFrom', 'datenum');
    end
end

if opt.quiet
    return
end
width_name = max(34, max(strlength(shown)) + 2);
fprintf('%-*s %-9s %12s  %s\n', width_name, 'file', 'exists', 'size', 'modified');
for k = 1:n_path
    if exists(k)
        fprintf('%-*s %-9s %9.1f MB  %s\n', width_name, shown(k), 'yes', ...
            bytes(k)/1e6, string(modified(k), 'yyyy-MM-dd HH:mm'));
    else
        fprintf('%-*s %-9s\n', width_name, shown(k), 'MISSING');
    end
end
end
