function [inside, outside] = bout_cutwindows(bout_table, window_table, min_s)
%BOUT_CUTWINDOWS  Each bout split at the windows: the pieces inside any window, and the pieces outside all.
%   Caller: bout_puffmove.zhang_loco
%
% IN   bout_table    N x 2 float s   [start, end] bouts to cut
%      window_table  M x 2 float s   [start, end] windows; may overlap one another
%      min_s         1 x 1 float s   a piece shorter than this is dropped
% OUT  inside        P x 2 float s   pieces that lie in a window
%      outside       Q x 2 float s   pieces that lie in no window
    inside = zeros(0, 2);
    outside = zeros(0, 2);
    if isempty(bout_table)
        return
    end
    if isempty(window_table)
        outside = bout_table;
        return
    end
    % windows merged into disjoint spans, so each bout is walked once
    sorted_window = sortrows(window_table, 1);
    span = sorted_window(1, :);
    for k = 2:size(sorted_window, 1)
        if sorted_window(k, 1) <= span(end, 2)
            span(end, 2) = max(span(end, 2), sorted_window(k, 2));
        else
            span(end+1, :) = sorted_window(k, :); %#ok<AGROW>
        end
    end
    for k = 1:size(bout_table, 1)
        bout_start = bout_table(k, 1);
        bout_end = bout_table(k, 2);
        cursor = bout_start;
        touching = span(span(:, 2) > bout_start & span(:, 1) < bout_end, :);
        for j = 1:size(touching, 1)
            piece_start = max(touching(j, 1), bout_start);
            piece_end = min(touching(j, 2), bout_end);
            if piece_start - cursor >= min_s
                outside(end+1, :) = [cursor, piece_start]; %#ok<AGROW>
            end
            if piece_end - piece_start >= min_s
                inside(end+1, :) = [piece_start, piece_end]; %#ok<AGROW>
            end
            cursor = piece_end;
        end
        if bout_end - cursor >= min_s
            outside(end+1, :) = [cursor, bout_end]; %#ok<AGROW>
        end
    end
end
