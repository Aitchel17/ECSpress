function hit = bout_hasoverlap(bout_table, other_table, pre_s, post_s)
%BOUT_HASOVERLAP  Per row of bout_table, whether any row of other_table touches [start - pre_s, end + post_s].
%   Caller: state_integration.awake_integration
%
% IN   bout_table   N x 2 float s   [start, end] bouts to label
%      other_table  M x 2 float s   [start, end] bouts to look for
%      pre_s        1 x 1 float s   window extension before each bout
%      post_s       1 x 1 float s   window extension after each bout
% OUT  hit          N x 1 bool
    hit = false(size(bout_table, 1), 1);
    if isempty(bout_table) || isempty(other_table)
        return
    end
    window_start = bout_table(:, 1) - pre_s;
    window_end = bout_table(:, 2) + post_s;
    for k = 1:size(bout_table, 1)
        hit(k) = any(other_table(:, 1) <= window_end(k) & other_table(:, 2) >= window_start(k));
    end
end
