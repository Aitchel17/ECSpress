function [swinglist, info] = event_findswing(dtrace, t_axis, state_idx, opt)
%EVENT_FINDSWING  The largest excursions INSIDE each state, as eventlist rows.
%   The row event_findquiet's header left room for: (state, pol) where state is a
%   state name and pol is a direction, not 'none'. Rows carry the same fields
%   event_pick_excursions produces, so they concatenate into one eventlist.
%
%   It answers a question a transition cannot. A transition row has a large
%   diameter change AND a state change, so a slope drawn through transitions and
%   controls cannot say which of the two moved the tissue. These rows hold the
%   state fixed and let the amplitude vary.
%
%   IT RANKS, IT DOES NOT THRESHOLD. There is no defensible amplitude threshold
%   for "an event happened", so candidates are ordered by how far the diameter
%   moved and the caller takes as many as it wants. n_keep is a CEILING, never a
%   quota -- an ineligible candidate is skipped, not promoted, and a state that
%   cannot fill n_keep reports the shortfall rather than reaching further down.
%
%   RUN IT ON THE PURE STATE FIELDS, not the rough ones. state_integration builds
%   both: roughnrem is NREM with uArousal folded in, so a microarousal sits inside
%   a bout and its constriction would come back as the largest swing in NREM.
%   Plain nrem splits the bout at the microarousal instead, so it cannot happen.
%
% IN   dtrace      1 x T float   raw trace, for diameter_change_raw
%      t_axis      1 x T float   frame times (s)
%      state_idx   struct        [nBout x 2] frame-index bout tables per state
%      states      cell          which state_idx fields to search. [] = all
%      n_keep      int           most rows per (state, direction). Default 10
%      min_rise_s  float s       an excursion shorter than this is not an event.
%                                Default 3
%      edge_margin int frames    frames an endpoint must keep clear of the bout
%                                edge. Required, because it MUST cover the PIV
%                                halfwin -- a smaller one lets the correlation
%                                window reach outside the bout, and the row is
%                                then not inside one state. There is no
%                                defensible default without knowing that halfwin
%
%   AN EXCURSION'S TWO ENDS ARE BOTH TURNING POINTS, so a rise that is still
%   rising when the searched band ends is not a candidate. That is on purpose:
%   its peak is outside the band, so its amplitude would be a truncated
%   measurement, and the shape it describes -- a bout that begins or ends
%   mid-dilation -- is a transition tail rather than a swing. A bout whose band
%   holds fewer than two turning points contributes nothing and is counted in
%   no_extrema, so it is distinguishable from one that was never searched.
%      exclude     1 x T bool    frames no excursion may touch. [] = none
%      dsg         1 x T float   smoothed trace; the extrema and every reported
%                                amplitude are read off it so the rows match the
%                                transitions'. [] = use dtrace
% OUT  swinglist   1 x M struct  state, pol, bout, start_f, end_f, start_sec,
%                                end_sec, search_from, search_to, from, to, rise_s,
%                                diameter_change, diameter_change_raw, range, sd
%                                three nested spans, as event_findquiet: start_f/
%                                end_f is the BOUT, search_from/search_to the band
%                                inside it the margins leave, from/to the
%                                excursion itself
%      info        struct        per state and direction: how many candidates were
%                                found, how many each rule rejected, how many were
%                                kept, and whether n_keep was reached. Also
%                                <direction>_amp, EVERY eligible candidate's
%                                amplitude in descending order, kept or not --
%                                which is what says whether the kept rows are a
%                                separated few or the top of one continuous
%                                distribution. searched_s makes the counts
%                                comparable between states of different length
%
%   THE RANK IS NOT A FIELD. Within one (state, pol) the rows come back ordered by
%   amplitude, and sorting on diameter_change reproduces that at any later point.
%   Storing it as well would put a function of an existing column beside it.
%
%   see also EVENT_FINDQUIET, EVENT_PICK_EXCURSIONS, ANALYSIS_EVENT

arguments
    dtrace          (1,:) double
    t_axis          (1,:) double
    state_idx       struct
    opt.states      cell = {}
    opt.n_keep      (1,1) double {mustBeInteger, mustBePositive} = 10
    opt.min_rise_s  (1,1) double {mustBePositive} = 3
    opt.edge_margin (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    opt.exclude     = []
    opt.dsg         (1,:) double = []
end

% 0. Setup
T = numel(dtrace);
if numel(t_axis) ~= T
    error('event_findswing:lengthMismatch', ...
        'dtrace has %d samples but t_axis has %d.', T, numel(t_axis));
end
dsg = opt.dsg;
if isempty(dsg)
    dsg = dtrace;
end
fps      = 1 / mean(diff(t_axis));
min_rise = max(1, round(opt.min_rise_s * fps));
% 0.1 Forbidden frames, as a mask
excluded = false(1, T);
if ~isempty(opt.exclude)
    if islogical(opt.exclude)
        excluded = opt.exclude(:).' | excluded;
    else
        excluded(opt.exclude(opt.exclude >= 1 & opt.exclude <= T)) = true;
    end
end
% 0.2 Which states to walk
states = opt.states;
if isempty(states)
    states = fieldnames(state_idx)';
end
% 0.3 Output template, field for field as event_pick_excursions leaves them.
%     Held as its own variable because each state's candidate pool has to START
%     empty -- initialising it from swinglist would carry the previous state's
%     kept rows into the next state's ranking, and they would be taken twice
row_template = struct('state',{},'pol',{},'bout',{},'start_f',{},'end_f',{}, ...
    'start_sec',{},'end_sec',{},'search_from',{},'search_to',{}, ...
    'from',{},'to',{},'rise_s',{},'diameter_change',{},'diameter_change_raw',{}, ...
    'range',{},'sd',{});
swinglist = row_template;
counts = struct();
directions = ["dilation", "constriction"];

% 1. One state at a time. Candidates pool across the state's bouts, because the
%    ranking is over the state and not over one bout
for s = 1:numel(states)
    name = states{s};
    if ~isfield(state_idx, name)
        continue
    end
    bouts = state_idx.(name);
    tally = struct('bouts', size(bouts,1), 'candidates', 0, ...
        'too_short', 0, 'near_edge', 0, 'excluded', 0, 'searched_s', 0);
    pool = row_template;

    for i = 1:size(bouts, 1)
        b1 = max(1, bouts(i,1));
        b2 = min(T, bouts(i,2));
        % 1.1 The band the margins leave. An endpoint outside it puts the PIV
        %     window across the bout edge, which is another state
        s1 = b1 + opt.edge_margin;
        s2 = b2 - opt.edge_margin;
        if s2 - s1 < min_rise
            continue
        end
        tally.searched_s = tally.searched_s + (t_axis(s2) - t_axis(s1));
        % 1.2 Every rise and every fall between consecutive turning points
        [ext_loc, ext_is_max] = alternating_extrema(dsg, s1, s2);
        tally.candidates = tally.candidates + max(0, numel(ext_loc) - 1);
        for k = 1:numel(ext_loc) - 1
            f = ext_loc(k);
            t = ext_loc(k+1);
            if t - f < min_rise
                tally.too_short = tally.too_short + 1;
                continue
            end
            if f < s1 || t > s2
                tally.near_edge = tally.near_edge + 1;
                continue
            end
            seg = f:t;
            if any(excluded(seg))
                tally.excluded = tally.excluded + 1;
                continue
            end
            if ext_is_max(k)
                pol = 'constriction';
            else
                pol = 'dilation';
            end
            pool(end+1) = struct( ...
                'state', name, 'pol', pol, 'bout', i, ...
                'start_f', b1, 'end_f', b2, ...
                'start_sec', t_axis(b1), 'end_sec', t_axis(b2), ...
                'search_from', [s1 s2], 'search_to', [s1 s2], ...
                'from', f, 'to', t, ...
                'rise_s', t_axis(t) - t_axis(f), ...
                'diameter_change',     dsg(t)    - dsg(f), ...
                'diameter_change_raw', dtrace(t) - dtrace(f), ...
                'range', max(dsg(seg)) - min(dsg(seg)), ...
                'sd',    std(dsg(seg))); %#ok<AGROW>
        end
    end

    % 1.3 Rank within each direction and take at most n_keep. Fewer is reported,
    %     never filled from candidates that failed a rule above
    for d = directions
        pick = pool(strcmp({pool.pol}, d));
        amp  = abs([pick.diameter_change]);
        [amp, ord] = sort(amp, 'descend');
        pick = pick(ord);
        n_take = min(opt.n_keep, numel(pick));
        tally.(d + "_eligible") = numel(pick);
        tally.(d + "_kept")     = n_take;
        tally.(d + "_short")    = opt.n_keep - n_take;
        % Every eligible candidate's amplitude, kept or not. Whether the rows
        % that were taken are a separated few or the top of one continuous
        % distribution is not answerable from a count, and rerunning to find out
        % would mean re-deciding the rules that produced the count
        tally.(d + "_amp") = amp;
        swinglist = [swinglist, pick(1:n_take)]; %#ok<AGROW>
    end
    counts.(name) = tally;
end

% 2. Chronological, as event_pick_excursions leaves its own. The amplitude order
%    is not lost -- sorting on diameter_change inside a (state, pol) recovers it
if ~isempty(swinglist)
    [~, ord] = sort([swinglist.from]);
    swinglist = swinglist(ord);
end

info = struct('counts', counts, 'n_keep', opt.n_keep, ...
    'min_rise_s', opt.min_rise_s, 'edge_margin', opt.edge_margin, ...
    'fps', fps, 'n_excluded', nnz(excluded), 'states', {states});
end

% ---------------------------------------------------------------------------
function [loc, is_max] = alternating_extrema(y, i1, i2)
%ALTERNATING_EXTREMA  Turning points of y(i1:i2), in original indices.
%   Maxima and minima interleave by construction -- between two maxima there is a
%   minimum -- so consecutive entries bound one rise or one fall.
    seg = y(i1:i2);
    [~, loc_max] = findpeaks(seg);
    [~, loc_min] = findpeaks(-seg);
    loc_max = loc_max(:).' + i1 - 1;
    loc_min = loc_min(:).' + i1 - 1;
    loc     = [loc_max, loc_min];
    is_max  = [true(1, numel(loc_max)), false(1, numel(loc_min))];
    [loc, ord] = sort(loc);
    is_max = is_max(ord);
end
