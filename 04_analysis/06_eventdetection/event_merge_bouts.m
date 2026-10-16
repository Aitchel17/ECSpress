function pairs = event_merge_bouts(state_idx, merge, t_axis, search_s)
%EVENT_MERGE_BOUTS  Pair adjacent state bouts and build their two search windows.
%   Caller: analysis_event.arousaltransition
%   Structural step: reads the bout tables only, never the diameter trace.
%   The transition frame is search_from(2) == search_to(1).
%
% IN   state_idx  struct of [nBout x 2] frame-index bout tables, timetable2frame per table
%      merge      Nx4 cell {earlierState, laterState, label, polarity}
%      t_axis     1xT frame times (s); also sets the frame rate
%      search_s   search window each side of the transition (s), default 30
% OUT  pairs      1xM struct, chronological: state, pol, bout, start_f, end_f,
%                 start_sec, end_sec, search_from, search_to

arguments
    state_idx struct
    merge     cell
    t_axis    (1,:) double
    search_s  (1,1) double {mustBePositive} = 30
end

% 0. Setup
% 0.1 Check merge table shape
if size(merge, 2) ~= 4
    error('event_merge_bouts:mergeShape', ...
        'merge must be Nx4 {earlier,later,label,polarity}, got %d columns.', size(merge, 2));
end
% 0.2 Search half-width in frames
fps     = 1 / mean(diff(t_axis));
nsearch = max(1, round(search_s * fps));
% 0.3 Output template
pairs = struct('state',{},'pol',{},'bout',{},'start_f',{},'end_f',{}, ...
               'start_sec',{},'end_sec',{},'search_from',{},'search_to',{});

% 1. Pair bouts, one merge row at a time
for m = 1:size(merge, 1)
    % 1.1 Extract from container
    pre_state = merge{m,1};
    post_state = merge{m,2};
    transition_label = merge{m,3};
    polarity = merge{m,4};
    % 1.2 Bout tables of the earlier and the later state
    pre_bouts = state_idx.(pre_state);
    post_bouts = state_idx.(post_state);
    for i = 1:size(pre_bouts, 1)
        % 1.3 Find the later bout starting exactly where this earlier bout ends
        j = find(post_bouts(:,1) == pre_bouts(i,2), 1);
        if isempty(j); continue; end
        % 1.4 The two bouts' ends, in frames
        pre_first = pre_bouts(i,1);
        pre_last = pre_bouts(i,2);
        post_first = post_bouts(j,1);
        post_last = post_bouts(j,2);
        % 1.5 Search windows butted against the transition, nsearch frames at most, inside their own bout
        search_pre_start = max(pre_first, pre_last - nsearch + 1);
        search_post_end = min(post_last, post_first + nsearch - 1);
        % 1.6 The pair, filled in the template's field order
        pair = struct();
        pair.state = transition_label;
        pair.pol = polarity;
        pair.bout = i;
        pair.start_f = pre_first;
        pair.end_f = post_last;
        pair.start_sec = t_axis(pre_first);
        pair.end_sec = t_axis(post_last);
        pair.search_from = [search_pre_start, pre_last];
        pair.search_to = [post_first, search_post_end];
        pairs(end+1) = pair; %#ok<AGROW>
    end
end

% 2. Sort pairs chronologically
if ~isempty(pairs)
    [~, ord] = sort([pairs.start_f]);
    pairs = pairs(ord);
end
end
