function root = dirs_pickroot(pool_root, relative_path)
%DIRS_PICKROOT  The root of a pool that takes an item: the one holding it already, else the emptiest.
%   Caller: dirs_mirror, archive_session
%   A pool is the disks of one role, each item on exactly one of them. Two roots holding
%   the same item is a cloned disk, and stops here rather than one being picked.
%
% IN   pool_root      1 x K str   the roots of one role; one that is not mounted is ignored
%      relative_path  1 x 1 str   the item's path below a root, a folder or a file
% OUT  root           1 x 1 str   the root to use; "" when none of the pool is mounted

arguments
    pool_root (1, :) string
    relative_path (1, 1) string
end
root = "";
pool_root = pool_root(isfolder(pool_root));
if isempty(pool_root)
    return
end

candidate = fullfile(pool_root, relative_path);
holds = isfolder(candidate) | isfile(candidate);
if nnz(holds) > 1
    error('dirs_pickroot:twoRoots', '%s is under more than one root: %s', ...
        relative_path, strjoin(pool_root(holds), ', '));
end
if any(holds)
    root = pool_root(holds);
    return
end

free_bytes = zeros(size(pool_root));
for k = 1:numel(pool_root)
    folder = java.io.File(char(pool_root(k)));
    free_bytes(k) = double(folder.getUsableSpace());
end
[~, emptiest] = max(free_bytes);
root = pool_root(emptiest);
end
