function bouts = stateidx2bouts(stateidx)
%STATEIDX2BOUTS Separate indices into continuous bouts
%   bouts = stateidx2bouts(stateidx)
%   stateidx: sorted vector of indices
%   bouts: cell array where each cell contains a continuous segment of indices

if isempty(stateidx)
    bouts = {};
    return;
end

% Ensure sorted column vector
stateidx = sort(stateidx(:));

% Find discontinuities (diff > 1)
% breaks are the indices in 'stateidx' where the jump occurs
% i.e., stateidx(k+1) - stateidx(k) > 1
jump_locs = find(diff(stateidx) > 1);

% Define start and end indices for each bout within the stateidx vector
bout_starts = [1; jump_locs + 1];
bout_ends = [jump_locs; length(stateidx)];

% Create cell array
bouts = cell(length(bout_starts), 1);
for i = 1:length(bout_starts)
    bouts{i} = stateidx(bout_starts(i):bout_ends(i));
end
end

