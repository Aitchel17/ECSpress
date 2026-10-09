function trimmed_table = bout_trim(time_table, max_time_s, margin_s)
%BOUT_TRIM  The *_trans windows clipped to a measurement's end, less one transition window.
%   Caller: centralize_statesession, centralize_driftsession, stateintegration_test
%   Only a transition window needs the margin: it keeps get_transitionsummary's centre
%   +/- half-window inside the record. Every other table is left whole, because
%   bout_frames keeps only the frames inside the record, so the tail is scored too.
%
% IN   time_table     1x1 struct    table name -> N x 2 [start end] bouts in seconds
%      max_time_s     1x1 float s   the measurement's last time point
%      margin_s       1x1 float s   taken off max_time_s first, the transition window
% OUT  trimmed_table  1x1 struct    the same names; a transition window begun past the
%                                   limit is dropped, one running past it ends there
    limit_s = max(0, max_time_s - margin_s);
    fprintf('[bout_trim] measurement ends %.1f s, margin %.1f s, transition windows clipped to %.1f s\n', ...
        max_time_s, margin_s, limit_s);
    trimmed_table = time_table;
    table_names = fieldnames(time_table);
    is_transition = endsWith(table_names, "_trans");
    for k = find(is_transition)'
        bouts = time_table.(table_names{k});
        if isempty(bouts)
            continue
        end
        keep_bout = bouts(:, 1) < limit_s;
        bouts = bouts(keep_bout, :);
        is_overrun = bouts(:, 2) > limit_s;
        bouts(is_overrun, 2) = limit_s;
        trimmed_table.(table_names{k}) = bouts;
    end
end
