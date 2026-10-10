function [hit, miss] = bout_splitat(bout_table, window_table, min_s)
%BOUT_SPLITAT  Bouts sorted against windows by where they BEGIN: one begun in a window is the window's whole, any other is cut at the window edges.
%   Caller: bout_puffmove.cut_at_puffs
%   A run the puff starts stays the puff's until the mouse stops, so it is not cut. A run that was
%   already going when the puff came keeps what lies before the window and after it as its own,
%   and hands the window itself to the puff.
%
% IN   bout_table    N x 2 float s   [start, end] bouts
%      window_table  M x 2 float s   [start, end] windows; may overlap one another
%      min_s         1 x 1 float s   a piece shorter than this is dropped
% OUT  hit           P x 2 float s   whole bouts begun in a window, and the in-window pieces of the rest
%      miss          Q x 2 float s   the out-of-window pieces of bouts begun outside every window
    hit = zeros(0, 2);
    miss = zeros(0, 2);
    if isempty(bout_table)
        return
    end
    if isempty(window_table)
        miss = bout_table;
        return
    end
    begun_inside = false(size(bout_table, 1), 1);
    for k = 1:size(bout_table, 1)
        begun_inside(k) = any(window_table(:, 1) <= bout_table(k, 1) & bout_table(k, 1) <= window_table(:, 2));
    end
    hit = bout_table(begun_inside, :);
    [cut_in, cut_out] = bout_cutwindows(bout_table(~begun_inside, :), window_table, min_s);
    hit = sortrows([hit; cut_in], 1);
    miss = cut_out;
end
