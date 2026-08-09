function pairs = event_merge_bouts(state_idx, merge, t_axis, search_s)
%EVENT_MERGE_BOUTS  Pair adjacent state bouts and build their two search windows.
%   Structural step: reads the bout tables only, never the diameter trace.
%   The transition frame is search_from(2) == search_to(1).
%
% IN   state_idx  struct of [nBout x 2] frame-index bout tables (img_state.state_idx)
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
    % 1.1 Bout tables of the earlier (A) and later (B) state
    A = state_idx.(merge{m,1});
    B = state_idx.(merge{m,2});
    for i = 1:size(A, 1)
        % 1.2 Find the B bout starting exactly where this A bout ends
        j = find(B(:,1) == A(i,2), 1);
        if isempty(j); continue; end
        % 1.3 Build the pair: search windows butted against the transition,
        %     clamped to their own bout
        pairs(end+1) = struct( ...
            'state', merge{m,3}, 'pol', merge{m,4}, 'bout', i, ...
            'start_f', A(i,1), 'end_f', B(j,2), ...
            'start_sec', t_axis(A(i,1)), 'end_sec', t_axis(B(j,2)), ...
            'search_from', [max(A(i,1), A(i,2) - nsearch + 1), A(i,2)], ...
            'search_to',   [B(j,1), min(B(j,2), B(j,1) + nsearch - 1)]); %#ok<AGROW>
    end
end

% 2. Sort pairs chronologically
if ~isempty(pairs)
    [~, ord] = sort([pairs.start_f]);
    pairs = pairs(ord);
end
end
