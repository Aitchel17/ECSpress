function [events, info] = piv_detect_events(dtrace, t_axis, state_idx, opt)
%PIV_DETECT_EVENTS  One diameter excursion per adjacent sleep-state bout pair.
%   Extracted from piv_integration_testbed so the detection can be re-run and
%   re-tuned without touching the PIV stage.
%
%   [events, info] = piv_detect_events(dtrace, t_axis, state_idx, ...)
%
%   WHAT IT DOES
%     For every pair of ADJACENT state bouts (A then B) named in opt.merge, it emits
%     exactly one event. The transition type PRESCRIBES the polarity -- arousal
%     (-> awake) constricts, awake -> nrem and nrem -> rem dilate -- and the data only
%     supplies the size. Detecting extrema instead (findpeaks etc.) makes every event
%     the same polarity by construction and silently inverts constricting transitions.
%
%   HOW THE TWO ENDPOINTS ARE CHOSEN
%     1. Smooth hard (Savitzky-Golay, sg_win_s) -> the picking trace.
%     2. SYMMETRIC search windows of search_s seconds each, both butted against the
%        transition: the LAST search_s of bout A and the FIRST search_s of bout B,
%        each clamped to its bout. A bout shorter than search_s is therefore searched
%        WHOLE -- which is what short states need (uArousal runs 10-15 s here).
%     3. Inside each window keep every sample in the extreme percentile band
%        (<= pct_lo for a minimum, >= pct_hi for a maximum) ...
%     4. ... and take the candidate CLOSEST TO THE TRANSITION: the last qualifying
%        sample in A, the first in B. Both picks shorten the from->to span, so the
%        event is the tightest excursion consistent with both endpoints being extreme.
%     Step 3-4 exist because a plain min/max over the window picks the single global
%     extreme, which for a long later bout lands far past the transition and swallows
%     whole sub-excursions (measured: a rem>awake TO landed ~115 s into the wake bout,
%     giving a 149 s "constriction"). Selecting a BAND and then choosing by POSITION
%     anchors both endpoints to the state change.
%     SG overshoot near steps is tolerable here precisely because no single peak is
%     ever localised -- only a percentile band, then a position.
%
%   INPUTS
%     dtrace     1xT absolute vessel diameter (um). NOT a median-subtracted change
%                trace: the percentile bands are taken on this directly.
%     t_axis     1xT time axis (s), same sampling as dtrace.
%     state_idx  struct of [nBout x 2] FRAME-index bout tables, i.e.
%                img_state.state_idx after get_state_indices. Call
%                sleep_integrate.trim_to_duration(t_axis(end)) FIRST: timetable2frame
%                resolves times by min|t - t_axis|, so untrimmed tail bouts all
%                collapse onto the last frame and silently swallow late events.
%   OPTIONS
%     fps        sampling rate; [] = derive from t_axis (default []).
%     merge      Nx4 cell {earlierField, laterField, label, polarity}
%                polarity is 'dilation' or 'constriction'. Default = the six sleep
%                transitions. Fields must exist in state_idx.
%                NB the NREM side is the RAW 'nrem', not 'roughnrem': roughnrem folds
%                uArousal in, so a uArousal bout sits INSIDE it and 'uarousal'->
%                'roughnrem' has zero adjacencies (measured). Treating uArousal as its
%                own state is what makes the micro-arousal pair constructible.
%                The awake side stays 'roughawake' (Awake+Drowsy) because Drowsy sits
%                between Awake and NREM: raw 'awake'->'nrem' has zero adjacencies.
%     search_s   search-window length (s) on EACH side, default 30. Clamped to the
%                bout, so short bouts are searched whole.
%     min_rise_s minimum from->to separation (s), default 3.
%     smooth_s   movmean window (s) for the REPORTED amplitude trace, default 1.
%     pick_filter smoother for the PICKING trace: 'none' (DEFAULT -- pick straight off
%                the measured trace; the percentile band is the noise buffer),
%                'sgolay' (shape-preserving but overshoots at steps) or 'movmean'.
%                With 'none' the picking trace IS dtrace, so events.dD and .dD_raw
%                differ only by the light smooth_s movmean used for reporting.
%     sg_ord     Savitzky-Golay order for the PICKING trace, default 3.
%     sg_win_s   Savitzky-Golay window (s) for the picking trace, default 15.
%     band_mode  how the extreme band is defined, default 'amplitude':
%                'amplitude' -- a LEVEL band_frac of the way from the window median
%                    to the window extreme:  lvl = med +/- band_frac*(ext - med).
%                    Anchored to the excursion depth, so it is stable on a plateau.
%                'prctile' -- the older RANK band (pct_lo/pct_hi). On a flat stretch
%                    the rank threshold lands on the data value itself (measured: a
%                    pick cleared its p20 threshold by 0.007 um), so which frame wins
%                    is decided by noise rather than by the rule.
%     band_frac  amplitude-mode fraction toward the extreme, default 0.80.
%     pct_lo/hi  percentile band edges for band_mode='prctile', default 20 / 80.
%     adj_tol_s  bout end vs next bout start must match within this (s), default 1.
%                Bouts normally abut exactly in frame space; this only guards drift.
%     min_bout_s skip a merge if EITHER bout is shorter than this (s), default 0.
%     verbose    print the summary + checks, default true.
%
%   OUTPUTS
%     events  1xN struct, chronological:
%       .from .to        frame indices, ALWAYS from < to (time order -- the PIV pair
%                        order sets the sign of the measured displacement)
%       .pol             'dilation' | 'constriction' (prescribed, not measured)
%       .rise_s          (to-from)/fps
%       .dD              amplitude on the PICKING trace (selection-consistent)
%       .dD_raw          amplitude on the lightly-smoothed trace -- use THIS to
%                        report physical amplitude; strong SG flattens peaks
%       .state .bout     transition label and index of bout A
%       .bout_win        [A_start B_end] merged window (for plots)
%       .search_from/.search_to   the two search windows
%       .A_s .B_s        bout durations (s)
%     info    struct: dsg, smoothd, params, and the diagnostics below.
%
%   DIAGNOSTICS in info (also printed when verbose)
%     .sign_ok    every dilation has dD>0 and every constriction dD<0
%     .clipped    events whose extremum sat ON a search-window edge -> the excursion
%                 probably continues past it and dD is a lower bound
%     .n_skipped  merges rejected, by reason

arguments
    dtrace     (1,:) double
    t_axis     (1,:) double
    state_idx  struct
    opt.fps        double  = []
    opt.merge      cell    = { ...
%        earlier(A)    later(B)      label             polarity
        'nrem',       'rem',        'nrem2rem',       'dilation'     ; ...
        'uarousal',   'nrem',       'uarousal2nrem',  'dilation'     ; ...
        'roughawake', 'nrem',       'awake2nrem',     'dilation'     ; ...
        'rem',        'roughawake', 'rem2awake',      'constriction' ; ...
        'nrem',       'roughawake', 'nrem2awake',     'constriction' ; ...
        'nrem',       'uarousal',   'nrem2uarousal',  'constriction' }
    opt.search_s   (1,1) double  = 30
    opt.min_rise_s (1,1) double  = 3.0
    opt.smooth_s   (1,1) double  = 1.0
    opt.pick_filter (1,:) char {mustBeMember(opt.pick_filter, {'none','sgolay','movmean'})} = 'sgolay'
    opt.sg_ord     (1,1) double  = 3
    opt.sg_win_s   (1,1) double  = 3
    opt.band_mode  (1,:) char {mustBeMember(opt.band_mode, {'amplitude','prctile'})} = 'amplitude'
    opt.band_frac  (1,1) double  = 0.80
    opt.pct_lo     (1,1) double  = 20
    opt.pct_hi     (1,1) double  = 80
    opt.adj_tol_s  (1,1) double  = 1.0
    opt.min_bout_s (1,1) double  = 0
    opt.verbose    (1,1) logical = true
end

nT = numel(dtrace);
if numel(t_axis) < nT
    error('piv_detect_events:axisShort', ...
        't_axis has %d samples but dtrace has %d.', numel(t_axis), nT);
end
fps = opt.fps;
if isempty(fps); fps = 1 / mean(diff(t_axis(1:nT))); end

if size(opt.merge, 2) ~= 4
    error('piv_detect_events:mergeShape', ...
        'merge must be Nx4 {earlier,later,label,polarity}, got %d columns.', ...
        size(opt.merge, 2));
end
need = unique(opt.merge(:, 1:2));
gone = need(~isfield(state_idx, need));
if ~isempty(gone)
    error('piv_detect_events:missingBout', ...
        'state_idx has no field(s): %s', strjoin(gone', ', '));
end

%% traces
smoothd  = movmean(dtrace, max(3, round(opt.smooth_s * fps)), 'omitnan');
sg_frame = 2*floor(opt.sg_win_s * fps / 2) + 1;          % sgolayfilt needs an odd frame
if ~strcmp(opt.pick_filter, 'none') && sg_frame <= opt.sg_ord
    error('piv_detect_events:sgWindow', ...
        'sg_win_s=%g gives frame %d, which must exceed sg_ord=%d.', ...
        opt.sg_win_s, sg_frame, opt.sg_ord);
end
% MINIMAL smoothing by default (sg_win_s = 3 s). The width matters far more than the
% filter type. Scored by "is the picked frame actually extreme on the measured trace"
% (0 = perfect, 2 = opposite), over this session's events:
%     sgolay  3 s -> p90 0.40, |dD-dD_raw| max 0.12 um     <- default
%     sgolay  5 s -> p90 0.44, max 0.29
%     sgolay 10 s -> p90 0.61, max 0.87
%     sgolay 15 s -> p90 1.16, max 1.58   (endpoints the data does not support)
%     none        -> p90 0.73, max 0.47
% Events run 5-48 s, so a 15 s window was wider than many of them and dragged the
% pick off. Unsmoothed is not the answer either: the percentile band absorbs a spike
% into the band, but the endpoint is then chosen by POSITION, so a lone noisy frame
% can still win. A 3 s window (9 samples here) removes exactly that without touching
% excursion shape. sgolay also OVERSHOOTS at steps -- harmless this narrow, real at 15 s.
switch opt.pick_filter
    case 'none';    dsg = dtrace;
    case 'sgolay';  dsg = sgolayfilt(dtrace, opt.sg_ord, sg_frame);
    case 'movmean'; dsg = movmean(dtrace, sg_frame, 'omitnan');
end
mind    = max(1, round(opt.min_rise_s * fps));
nsearch = max(1, round(opt.search_s * fps));   % symmetric search length (frames)
tol  = max(1, round(opt.adj_tol_s * fps));
minb = round(opt.min_bout_s * fps);

%% one event per adjacent bout pair
events  = struct('from',{},'to',{},'pol',{},'rise_s',{},'dD',{},'dD_raw',{}, ...
                 'state',{},'bout',{},'bout_win',{},'search_from',{},'search_to',{}, ...
                 'A_s',{},'B_s',{});
skipped = struct('not_adjacent',0, 'short_bout',0, 'empty_band',0, 'too_brief',0);
for m = 1:size(opt.merge, 1)
    A = state_idx.(opt.merge{m,1});
    B = state_idx.(opt.merge{m,2});
    if isempty(A) || isempty(B); continue; end
    isdil = strcmp(opt.merge{m,4}, 'dilation');
    for i = 1:size(A, 1)
        j = find(abs(B(:,1) - A(i,2)) <= tol, 1);
        if isempty(j); skipped.not_adjacent = skipped.not_adjacent + 1; continue; end
        La = A(i,2) - A(i,1) + 1;   Lb = B(j,2) - B(j,1) + 1;
        if La < minb || Lb < minb;  skipped.short_bout = skipped.short_bout + 1; continue; end

        % symmetric absolute windows butted against the transition, clamped to the
        % bout -- a bout shorter than search_s is searched whole
        wa = max(A(i,1), A(i,2) - nsearch + 1) : A(i,2);   % last search_s of A
        wb = B(j,1) : min(B(j,2), B(j,1) + nsearch - 1);   % first search_s of B

        % extreme percentile band, then the candidate nearest the state change:
        % the transition is at the END of A (take the last) and the START of B (first)
        va = dsg(wa);   vb = dsg(wb);
        if isdil
            ca = wa(band_mask(va, 'lo', opt));
            cb = wb(band_mask(vb, 'hi', opt));
        else
            ca = wa(band_mask(va, 'hi', opt));
            cb = wb(band_mask(vb, 'lo', opt));
        end
        if isempty(ca) || isempty(cb); skipped.empty_band = skipped.empty_band + 1; continue; end
        f = ca(end);   t = cb(1);
        if (t - f) < mind; skipped.too_brief = skipped.too_brief + 1; continue; end

        events(end+1) = struct('from', f, 'to', t, 'pol', opt.merge{m,4}, ...
            'rise_s', (t - f)/fps, 'dD', dsg(t) - dsg(f), ...
            'dD_raw', smoothd(t) - smoothd(f), ...
            'state', opt.merge{m,3}, 'bout', i, 'bout_win', [A(i,1) B(j,2)], ...
            'search_from', [wa(1) wa(end)], 'search_to', [wb(1) wb(end)], ...
            'A_s', La/fps, 'B_s', Lb/fps); %#ok<AGROW>
    end
end
if ~isempty(events)
    [~, ord] = sort([events.from]);   events = events(ord);
end

%% diagnostics
sign_bad = false(1, numel(events));
clipped  = false(1, numel(events));
for e = 1:numel(events)
    isdil = strcmp(events(e).pol, 'dilation');
    sign_bad(e) = (isdil && events(e).dD <= 0) || (~isdil && events(e).dD >= 0);
    clipped(e)  = events(e).from == events(e).search_from(1) || ...
                  events(e).to   == events(e).search_to(2);
end
info = struct('dsg', dsg, 'smoothd', smoothd, 'fps', fps, 'sg_frame', sg_frame, ...
    'mind', mind, 'merge', {opt.merge}, 'opt', opt, ...
    'sign_ok', ~any(sign_bad), 'sign_bad', find(sign_bad), ...
    'clipped', find(clipped), 'n_skipped', skipped);

if ~opt.verbose; return; end
fprintf('\n--- piv_detect_events: %d events from %d merge rules ---\n', ...
    numel(events), size(opt.merge, 1));
for m = 1:size(opt.merge, 1)
    inM = strcmp({events.state}, opt.merge{m,3});
    if ~any(inM); fprintf('  %-11s none\n', opt.merge{m,3}); continue; end
    fprintf(['  %-11s %-12s n=%2d   dur %5.1f-%5.1f s (med %5.1f)   ' ...
             '|dD| %.2f-%.2f um (med %.2f)\n'], opt.merge{m,3}, opt.merge{m,4}, nnz(inM), ...
        min([events(inM).rise_s]), max([events(inM).rise_s]), median([events(inM).rise_s]), ...
        min(abs([events(inM).dD])), max(abs([events(inM).dD])), median(abs([events(inM).dD])));
end
if any(structfun(@(x) x > 0, skipped))
    fprintf('  skipped: not_adjacent %d, short_bout %d, empty_band %d, too_brief %d\n', ...
        skipped.not_adjacent, skipped.short_bout, skipped.empty_band, skipped.too_brief);
end
if info.sign_ok
    fprintf('  sign check : all %d events have dD consistent with their polarity\n', numel(events));
else
    warning('piv_detect_events:signMismatch', ...
        'events %s have dD opposite to the prescribed polarity', mat2str(info.sign_bad));
end
if isempty(info.clipped)
    fprintf('  edge check : no extremum sat on a search-window edge\n');
else
    fprintf('  edge check : %d event(s) hit a search-window edge (dD is a lower bound): %s\n', ...
        numel(info.clipped), mat2str(info.clipped));
end
end

% ---------------------------------------------------------------------------
function m = band_mask(v, side, opt)
%BAND_MASK  Which samples of v count as "extreme" on the requested side.
%   amplitude: a level band_frac of the way from the median to the extreme. The
%   median is a robust centre and (ext - med) is the excursion's own depth, so the
%   level scales with the event instead of with the sample distribution. A flat
%   window collapses to lvl == med and everything qualifies, which is the honest
%   answer there -- position then decides, as it should.
    md = median(v, 'omitnan');
    switch opt.band_mode
        case 'amplitude'
            switch side
                case 'lo'; m = v <= md - opt.band_frac * (md - min(v));
                case 'hi'; m = v >= md + opt.band_frac * (max(v) - md);
            end
        case 'prctile'
            switch side
                case 'lo'; m = v <= prctile(v, opt.pct_lo);
                case 'hi'; m = v >= prctile(v, opt.pct_hi);
            end
    end
end
