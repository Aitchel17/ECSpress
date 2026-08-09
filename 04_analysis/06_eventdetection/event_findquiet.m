function [quietlist, info] = event_findquiet(dtrace, t_axis, state_idx, opt)
%EVENT_FINDQUIET  Stable stretches inside each state bout, as eventlist rows.
%   The control every transition event needs. Rows carry the SAME fields
%   event_pick_excursions produces, so they concatenate into one eventlist.
%   'quiet' is the name, not the claim: stable end to end, not still.
%
%   THREE RULES, each replacing a scheme that was tried first -- CLAUDE_LOG.md
%     1. Search search_frac of each bout, the middle by default
%     2. Match the two ENDS, |dsg(to) - dsg(from)| <= max_dd. NOT flatness
%     3. Pick at RANDOM among those that qualify, not the flattest. seed fixed
%   max_range caps the excursion between the ends; Inf by default, because
%   gating on it brings rule 2's problem back. Filter the reported range instead.
%
%   SEVERAL LENGTHS PER BOUT, because most of the PIV error accumulates with LAG.
%
%   A row is (state, pol) -- two axes, so no third category is needed:
%     'nrem2rem', 'dilation'      a transition event
%     'nrem',     'none'          a stable stretch inside one nrem bout
%     'rem',      'dilation'      room for a within-state event later
%
% IN   dtrace      1 x T float   raw trace, for diameter_change_raw
%      t_axis      1 x T float   frame times (s)
%      state_idx   struct        [nBout x 2] frame-index bout tables per state
%      states      cell          which state_idx fields to search. [] = all
%      search_frac 1 x 2 float   the portion of each bout to search, as fractions
%                                of its span. Default [1/3 2/3]
%      max_dd      float         the most the two ENDS may differ, in the UNIT OF
%                                dtrace. Required -- there is no defensible default
%                                without knowing what an event looks like;
%                                line_fwhm.quietdetection derives it from the
%                                transitions it already holds
%      max_range   float         optional cap on the excursion BETWEEN the ends.
%                                Inf (default) = do not gate on it; it is reported
%                                on every row either way
%      len_frac    1 x K float   window length as a fraction of the search
%                                region's longest unexcluded run. Default [1 .5 .25]
%      len_min_s   float s       a window shorter than this is not reported
%      exclude     1 x T bool    frames no window may touch. [] = none
%      dsg         1 x T float   smoothed trace; range and diameter_change are read
%                                off it so the rows match the transitions'.
%                                [] = use dtrace
%      seed        int           RandStream seed. Its own stream, so the caller's
%                                rng state is untouched
% OUT  quietlist   1 x M struct  state, pol, bout, start_f, end_f, start_sec,
%                                end_sec, search_from, search_to, from, to, rise_s,
%                                diameter_change, diameter_change_raw, range, sd
%                                three nested spans: start_f/end_f is the BOUT,
%                                search_from/search_to the band scanned inside it,
%                                from/to the window itself
%      info        struct        per state: bouts, windows, and why lengths dropped
%
%   n is BOUTS, not windows. Lengths from one bout share its tissue and its local
%   conditions even when they no longer share frames.

arguments
    dtrace          (1,:) double
    t_axis          (1,:) double
    state_idx       struct
    opt.states      cell = {}
    opt.search_frac (1,2) double {mustBeInRange(opt.search_frac, 0, 1)} = [1/3 2/3]
    opt.max_dd      (1,1) double {mustBePositive}
    opt.max_range   (1,1) double {mustBePositive} = Inf
    opt.len_frac    (1,:) double {mustBePositive} = [1 0.5 0.25]
    opt.len_min_s   (1,1) double {mustBePositive} = 2
    opt.exclude     = []
    opt.dsg         (1,:) double = []
    opt.seed        (1,1) double = 0
end

% 0. Setup
T = numel(dtrace);
if numel(t_axis) ~= T
    error('event_findquiet:lengthMismatch', ...
        'dtrace has %d samples but t_axis has %d.', T, numel(t_axis));
end
dsg = opt.dsg;
if isempty(dsg)
    dsg = dtrace;
end
fps     = 1 / mean(diff(t_axis));
len_min = max(2, round(opt.len_min_s * fps));
stream  = RandStream('twister', 'Seed', opt.seed);
% 0.1 Forbidden frames, as a mask
excluded = false(1, T);
if ~isempty(opt.exclude)
    if islogical(opt.exclude); excluded = opt.exclude(:).' | excluded;
    else; excluded(opt.exclude(opt.exclude >= 1 & opt.exclude <= T)) = true;
    end
end
% 0.2 Which states to walk
states = opt.states;
if isempty(states)
    states = fieldnames(state_idx)';
end
% 0.3 Output template, field for field as event_pick_excursions leaves them
quietlist = struct('state',{},'pol',{},'bout',{},'start_f',{},'end_f',{}, ...
    'start_sec',{},'end_sec',{},'search_from',{},'search_to',{}, ...
    'from',{},'to',{},'rise_s',{},'diameter_change',{},'diameter_change_raw',{}, ...
    'range',{},'sd',{});
counts = struct();

% 1. One state at a time
for s = 1:numel(states)
    name = states{s};
    if ~isfield(state_idx, name)
        continue
    end
    bouts = state_idx.(name);
    counts.(name) = struct('bouts', size(bouts,1), 'too_short', 0, ...
        'no_matched_pair', 0, 'windows', 0);
    for i = 1:size(bouts, 1)
        b1 = max(1, bouts(i,1));   b2 = min(T, bouts(i,2));
        if b2 <= b1
            continue
        end
        % 1.1 The middle band of this bout, minus anything forbidden
        span = b2 - b1;
        s1 = b1 + round(opt.search_frac(1) * span);
        s2 = b1 + round(opt.search_frac(2) * span);
        free = false(1, T);   free(s1:s2) = true;   free = free & ~excluded;
        run_max = longest_run(free);
        if run_max < len_min
            counts.(name).too_short = counts.(name).too_short + 1;   continue
        end
        % 1.2 One window per requested length, at a RANDOM qualifying position
        len_list = unique(max(len_min, round(opt.len_frac * run_max)), 'stable');
        for len = len_list(len_list <= run_max)
            start = qualifying_starts(dsg, free, len, opt.max_dd, opt.max_range);
            if isempty(start)
                counts.(name).no_matched_pair = counts.(name).no_matched_pair + 1;
                continue
            end
            f = start(randi(stream, numel(start)));
            t = f + len - 1;
            seg = f:t;
            quietlist(end+1) = struct( ...
                'state', name, 'pol', 'none', 'bout', i, ...
                'start_f', b1, 'end_f', b2, ...
                'start_sec', t_axis(b1), 'end_sec', t_axis(b2), ...
                'search_from', [s1 s2], 'search_to', [s1 s2], ...
                'from', f, 'to', t, ...
                'rise_s', t_axis(t) - t_axis(f), ...
                'diameter_change',     dsg(t)    - dsg(f), ...
                'diameter_change_raw', dtrace(t) - dtrace(f), ...
                'range', max(dsg(seg)) - min(dsg(seg)), ...
                'sd',    std(dsg(seg))); %#ok<AGROW>
            counts.(name).windows = counts.(name).windows + 1;
        end
    end
end

% 2. Chronological, as event_pick_excursions leaves its own
if ~isempty(quietlist)
    [~, ord] = sort([quietlist.from]);   quietlist = quietlist(ord);
end
info = struct('counts', counts, 'len_frac', opt.len_frac, ...
    'len_min_s', opt.len_min_s, 'search_frac', opt.search_frac, ...
    'max_dd', opt.max_dd, 'max_range', opt.max_range, 'seed', opt.seed, 'fps', fps, ...
    'n_excluded', nnz(excluded), ...
    'n_bout', numel(unique(arrayfun(@(r) ...
        string(r.state) + "#" + r.bout, quietlist))));
end

% ---------------------------------------------------------------------------
function start = qualifying_starts(dsg, free, len, max_dd, max_range)
%QUALIFYING_STARTS  Window starts whose two ENDS read the same caliber.
%   Every position is scored at once rather than ranked -- the ranking is what
%   used to force the nesting.
    T = numel(dsg);
    if len > T
        start = [];
        return
    end
    first = 1:T-len+1;
    % 1. The pair. This is the whole criterion: the endpoint scheme compares the
    %    frames at first against the frames at first+len-1, so those two are what
    %    has to match
    pair = abs(dsg(first + len - 1) - dsg(first)) <= max_dd;
    % 2. Optional cap on what happened in between. Off by default -- see the note
    %    in the header on why gating this brings the flatness problem back
    if isfinite(max_range)
        span = movmax(dsg, [0 len-1]) - movmin(dsg, [0 len-1]);
        pair = pair & span(first) <= max_range;
    end
    % 3. Wholly inside the free band
    csum  = cumsum([0, double(free)]);
    whole = (csum(first + len) - csum(first)) == len;
    start = first(whole & pair);
end

% ---------------------------------------------------------------------------
function n = longest_run(mask)
%LONGEST_RUN  Longest stretch of consecutive true in a logical row.
    e = diff([false, mask(:).', false]);
    n = max([0, find(e == -1) - find(e == 1)]);
end
