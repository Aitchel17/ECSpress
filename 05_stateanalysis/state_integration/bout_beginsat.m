function hit = bout_beginsat(bout_table, time_s)
%BOUT_BEGINSAT  Per bout, whether it begins at one of the times, within a millisecond.
%   Caller: bout_puffmove.zhang_loco
%   A rest that was going when a puff came is cut at the puff onset, so its in-window piece begins
%   exactly there; a rest that began later inside the window does not.
%
% IN   bout_table  N x 2 float s   [start, end] bouts
%      time_s      K x 1 float s   candidate onsets
% OUT  hit         N x 1 bool
    hit = false(size(bout_table, 1), 1);
    for k = 1:size(bout_table, 1)
        hit(k) = any(abs(time_s(:) - bout_table(k, 1)) < 1e-3);
    end
end
