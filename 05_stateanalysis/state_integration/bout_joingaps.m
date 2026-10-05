function joined = bout_joingaps(bout_table, gap_s)
%BOUT_JOINGAPS  Bouts closer than gap_s merged into one.
%   Caller: state_integration.awake_integration
%
% IN   bout_table  N x 2 float s   [start, end], sorted by start, non-overlapping
%      gap_s       1 x 1 float s   a gap shorter than this joins its two neighbours
% OUT  joined      M x 2 float s   M <= N
    joined = zeros(0, 2);
    for k = 1:size(bout_table, 1)
        if ~isempty(joined) && bout_table(k, 1) - joined(end, 2) < gap_s
            joined(end, 2) = max(joined(end, 2), bout_table(k, 2));
        else
            joined(end + 1, :) = bout_table(k, :); %#ok<AGROW>
        end
    end
end
