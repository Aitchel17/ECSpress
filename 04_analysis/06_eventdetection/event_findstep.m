function [steplist, info] = event_findstep(dtrace, t_axis, state_idx, opt)
%EVENT_FINDSTEP  Graded movements of the caliber, found by velocity, not amplitude.
%   Rows carry the same fields event_pick_excursions produces, so they concatenate
%   into one eventlist. Most frames belong to no event: the record is partitioned
%   into moving and resting and only the moving part can become a row.
%   Every stage is blind to sign, so the constrictions are the null for the
%   dilations. A transition sits in no single bout and is dropped; those belong to
%   event_pick_excursions.
%
% IN   dtrace      1 x T float   raw caliber trace, PX
%      t_axis      1 x T float   frame times (s); fps derives from this
%      state_idx   struct        [nBout x 2] frame-index bout tables per state
%      states      cell          state_idx fields, in precedence order. {} = all
%      sg_win_s    float s       Savitzky-Golay window, order 3
%      base_win_s  float s       medfilt1 baseline window, 'truncate'
%      lp_hz       float Hz      band edge; the resolution floor is 1/(2*lp_hz)
%      amp_frac    float         level a merged object must cross, as a fraction of
%                                the local baseline caliber. 0 = off
%      min_area    float frac*s  smallest swept area an object may have. 0 = off
%      gate_pct    float 0..100  percentile of abs(d) that becomes the gate.
%                                REQUIRED; it fixes what fraction of the record rests
%      edge_margin int frames    clearance from the bout edge. REQUIRED; it must
%                                cover whatever halfwin a later stage uses
%      crossing    char          'level' | 'zero', which selection rule runs
%      baseline    char          'record' | 'bout', what the residual is against
%      dsg         1 x T float   trace the reported amplitudes read off.
%                                [] = use dtrace
% OUT  steplist    1 x M struct  state, pol, bout, start_f, end_f, start_sec,
%                                end_sec, search_from, search_to, from, to, rise_s,
%                                diameter_change, diameter_change_raw, range, sd
%      info        struct        derived windows, the measured gate, the rejection
%                                tally, every object's height and area, and arrays
%                                aligned to the returned rows
%
%   rise_s at or below info.floor_s is a lower bound, flagged info.row.censored.
%   Rows flagged info.row.edge sit where medfilt1('truncate') ran a different
%   baseline estimator. Identity is the index into eventlist, so adding rows
%   renumbers every downstream piv_run id.

arguments
    dtrace          (1,:) double
    t_axis          (1,:) double
    state_idx       struct
    opt.states      cell = {}
    opt.sg_win_s    (1,1) double {mustBePositive} = 1.5
    opt.base_win_s  (1,1) double {mustBePositive} = 300
    opt.lp_hz       (1,1) double {mustBePositive} = 0.5
    opt.amp_frac    (1,1) double {mustBeNonnegative} = 0.10
    opt.min_area    (1,1) double {mustBeNonnegative} = 0
    opt.gate_pct    (1,1) double {mustBeInRange(opt.gate_pct, 0, 100)}
    opt.edge_margin (1,1) double {mustBeInteger, mustBeNonnegative}
    opt.crossing    (1,:) char {mustBeMember(opt.crossing, {'level','zero'})} = 'level'
    opt.baseline    (1,:) char {mustBeMember(opt.baseline, {'record','bout'})} = 'record'
    opt.dsg         (1,:) double = []
end
sg_ord = 3;     % Savitzky-Golay polynomial order

% 0. Setup
T = numel(dtrace);
if numel(t_axis) ~= T
    error('event_findstep:lengthMismatch', ...
        'dtrace has %d samples but t_axis has %d.', T, numel(t_axis));
end
dsg = opt.dsg;
if isempty(dsg)
    dsg = dtrace;
end
if numel(dsg) ~= T
    error('event_findstep:dsgLengthMismatch', ...
        'dsg has %d samples but dtrace has %d.', numel(dsg), T);
end
fps = 1 / mean(diff(t_axis));
sg_frame = 2*floor(opt.sg_win_s * fps / 2) + 1;
if sg_frame <= sg_ord
    error('event_findstep:sgWindowShort', ...
        'sg_win_s %g s is %d frames at %.4f Hz; order %d needs more.', ...
        opt.sg_win_s, sg_frame, fps, sg_ord);
end
nyquist_hz = fps / 2;
if opt.lp_hz >= nyquist_hz
    error('event_findstep:lpAboveNyquist', ...
        'lp_hz %g Hz is at or above the %g Hz Nyquist of this trace.', ...
        opt.lp_hz, nyquist_hz);
end
base_frame = 2*floor(opt.base_win_s * fps / 2) + 1;
edge_frames = floor(base_frame / 2);
floor_s = 1 / (2 * opt.lp_hz);
min_frames = max(2, round(floor_s * fps));
states = opt.states;
if isempty(states)
    states = fieldnames(state_idx)';
end
steplist = new_rowtemplate();
tally = struct('runs_raw', 0, 'absorbed', 0, 'dropped_short', 0, ...
    'dropped_open', 0, 'dropped_nocross', 0, 'objects', 0, ...
    'height_inverted', 0, 'no_level_cross', 0, 'no_area', 0, 'spans_bout', 0);

% 1. De-quantise, subtract a baseline, band limit, differentiate
x = sgolayfilt(dtrace, sg_ord, sg_frame);
[span_lo, span_hi, in_span] = baseline_spans(state_idx, states, T, opt.baseline);
[b_lp, a_lp] = butter(2, opt.lp_hz / nyquist_hz, 'low');
base = median(x) * ones(1, T);
resid = zeros(1, T);
min_filt = 3 * 4;
for s = 1:numel(span_lo)
    seg = span_lo(s):span_hi(s);
    if numel(seg) <= min_filt
        continue
    end
    win = min(base_frame, 2*floor(numel(seg)/2) - 1);
    base_seg = medfilt1(x(seg), win, 'truncate');
    base(seg) = base_seg;
    resid(seg) = filtfilt(b_lp, a_lp, x(seg) - base_seg);
end
d = gradient(resid) * fps;      % gradient, not diff: it states the frame mapping

% 2. Partition into moving and resting. One mask and its complement, so a tie at
%    the gate goes to moving
gate = prctile(abs(d(in_span)), opt.gate_pct);
is_still = abs(d) < gate;
is_still(~in_span) = true;

% 3. Moving runs
[m1, m2, dir_run] = signed_runs(d, is_still);
is_still = frames_outside(m1, m2, T);
tally.runs_raw = numel(m1);

% 4. Absorb sub-floor stationary holes, same direction only. A turning point is a
%    short hole flanked by opposite directions and has to survive
for k = 1:numel(m1)-1
    gap_first = m2(k) + 1;
    gap_last = m1(k+1) - 1;
    gap_len = gap_last - gap_first + 1;
    if gap_len >= min_frames
        continue
    end
    if dir_run(k) ~= dir_run(k+1)
        continue
    end
    is_still(gap_first:gap_last) = false;
    tally.absorbed = tally.absorbed + 1;
end

% 5. Drop sub-floor moving runs. One pass is a fixed point, so do not loop
[m1, m2, ~] = signed_runs(d, is_still);
is_still = frames_outside(m1, m2, T);
for k = 1:numel(m1)
    run_len = m2(k) - m1(k) + 1;
    if run_len >= min_frames
        continue
    end
    is_still(m1(k):m2(k)) = true;
    tally.dropped_short = tally.dropped_short + 1;
end
[m1, m2, dir_run] = signed_runs(d, is_still);
is_still = frames_outside(m1, m2, T);
n_run_all = numel(m1);

% 6. Stationary regions and their heights. R runs, R+1 regions, region j precedes
%    run j. The height is the median of the whole region, not of one end
q1 = [1, m2 + 1];
q2 = [m1 - 1, T];
n_region = numel(q1);
hfrac = nan(1, n_region);
for j = 1:n_region
    if q2(j) < q1(j)
        continue
    end
    seg = q1(j):q2(j);
    h_here = median(resid(seg));
    b_here = median(base(seg));
    hfrac(j) = h_here / b_here;
end

% 6.1 A run with an empty flanking region has no anchor
run_ok = true(1, n_run_all);
for k = 1:n_run_all
    if isnan(hfrac(k)) || isnan(hfrac(k+1))
        run_ok(k) = false;
        tally.dropped_open = tally.dropped_open + 1;
    end
end
% 6.2 'zero' also demands the run change the side of the baseline
if strcmp(opt.crossing, 'zero')
    for k = 1:n_run_all
        if ~run_ok(k)
            continue
        end
        side_before = sign(hfrac(k));
        side_after = sign(hfrac(k+1));
        if side_before == side_after
            run_ok(k) = false;
            tally.dropped_nocross = tally.dropped_nocross + 1;
        end
    end
end
keep_idx = find(run_ok);
R = numel(keep_idx);

% 7. Merge. Greedy and non-overlapping; each pause must be shorter than the last,
%    which is what makes the chain terminate
n_run_row = [];
area_row = [];
area_all = [];
h_start_row = [];
h_end_row = [];
edge_row = [];
censored_row = [];
net_all = [];
a = 1;
while a <= R
    ka = keep_idx(a);
    b = a;
    prev_gap = t_axis(q1(ka+1)) - t_axis(q2(ka));
    while b < R
        kb = keep_idx(b);
        kn = keep_idx(b+1);
        if dir_run(kn) ~= dir_run(ka)
            break
        end
        gap = t_axis(q1(kn)) - t_axis(q2(kb));
        if ~(gap < prev_gap)
            break
        end
        if dir_run(ka) > 0 && ~(hfrac(kn+1) > hfrac(kb+1))
            break
        end
        if dir_run(ka) < 0 && ~(hfrac(kn+1) < hfrac(kb+1))
            break
        end
        b = b + 1;
        prev_gap = gap;
    end
    kb = keep_idx(b);

    % 8. Merged object. net_area is referenced to the object's own opening level,
    %    not to the baseline, which carries the state's offset
    from = q2(ka);
    to = q1(kb+1);
    h_start = hfrac(ka);
    h_end = hfrac(kb+1);
    net = h_end - h_start;
    obj_seg = from:to;
    obj_cal = median(base(obj_seg));
    net_area = trapz(t_axis(obj_seg), resid(obj_seg) - resid(from)) / obj_cal;
    tally.objects = tally.objects + 1;
    net_all(end+1) = abs(net); %#ok<AGROW>
    area_all(end+1) = abs(net_area); %#ok<AGROW>
    a = b + 1;

    % 8.1 The opening step is tested here, not in the loop, so a failure is counted
    %     rather than silently truncating the chain
    if dir_run(ka) > 0 && ~(h_end > h_start)
        tally.height_inverted = tally.height_inverted + 1;
        continue
    end
    if dir_run(ka) < 0 && ~(h_end < h_start)
        tally.height_inverted = tally.height_inverted + 1;
        continue
    end
    if dir_run(ka) > 0
        pol = 'dilation';
    else
        pol = 'constriction';
    end

    % 9. Selection
    keep = level_crossed(h_start, h_end, opt.amp_frac, opt.crossing);
    if ~keep
        tally.no_level_cross = tally.no_level_cross + 1;
        continue
    end
    if opt.min_area > 0 && abs(net_area) < opt.min_area
        tally.no_area = tally.no_area + 1;
        continue
    end

    % 10. State and bout, first match wins. An object is one event, never two rows
    [matched, name, i_bout, b1, b2, s1, s2] = ...
        first_bout(state_idx, states, from, to, T, opt.edge_margin);
    if ~matched
        tally.spans_bout = tally.spans_bout + 1;
        continue
    end

    % 11. Fill the row. Amplitudes read off dsg, never the detection trace
    seg = from:to;
    steplist(end+1) = struct( ...
        'state', name, 'pol', pol, 'bout', i_bout, ...
        'start_f', b1, 'end_f', b2, ...
        'start_sec', t_axis(b1), 'end_sec', t_axis(b2), ...
        'search_from', [s1 s2], 'search_to', [s1 s2], ...
        'from', from, 'to', to, ...
        'rise_s', t_axis(to) - t_axis(from), ...
        'diameter_change', dsg(to) - dsg(from), ...
        'diameter_change_raw', dtrace(to) - dtrace(from), ...
        'range', max(dsg(seg)) - min(dsg(seg)), ...
        'sd', std(dsg(seg))); %#ok<AGROW>
    in_chain = keep_idx >= ka & keep_idx <= kb;
    n_run_row(end+1) = nnz(in_chain); %#ok<AGROW>
    area_row(end+1) = net_area; %#ok<AGROW>
    h_start_row(end+1) = h_start; %#ok<AGROW>
    h_end_row(end+1) = h_end; %#ok<AGROW>
    edge_row(end+1) = (from <= edge_frames) || (to > T - edge_frames); %#ok<AGROW>
    censored_row(end+1) = (t_axis(to) - t_axis(from)) <= floor_s; %#ok<AGROW>
end

% 12. Chronological on exit; the aligned arrays follow the same order
if ~isempty(steplist)
    [~, ord] = sort([steplist.from]);
    steplist = steplist(ord);
    n_run_row = n_run_row(ord);
    area_row = area_row(ord);
    h_start_row = h_start_row(ord);
    h_end_row = h_end_row(ord);
    edge_row = edge_row(ord);
    censored_row = censored_row(ord);
end
net_all = sort(net_all, 'descend');
counts = count_states(steplist, states);

info = struct( ...
    'fps', fps, ...
    'sg_frame', sg_frame, ...
    'base_frame', base_frame, ...
    'edge_frames', edge_frames, ...
    'floor_s', floor_s, ...
    'min_frames', min_frames, ...
    'gate', gate, ...
    'still_frac', nnz(is_still) / T, ...
    'tally', tally, ...
    'net_all', net_all, ...
    'area_all', sort(area_all, 'descend'), ...
    'counts', counts, ...
    'row', struct('n_run', n_run_row, 'net_area', area_row, 'h_start', h_start_row, ...
        'h_end', h_end_row, 'edge', logical(edge_row), ...
        'censored', logical(censored_row)));
end


function template = new_rowtemplate()
%NEW_ROWTEMPLATE  The 16-field eventlist contract, as a 1 x 0 struct.
template = struct('state', {}, 'pol', {}, 'bout', {}, 'start_f', {}, 'end_f', {}, ...
    'start_sec', {}, 'end_sec', {}, 'search_from', {}, 'search_to', {}, ...
    'from', {}, 'to', {}, 'rise_s', {}, 'diameter_change', {}, ...
    'diameter_change_raw', {}, 'range', {}, 'sd', {});
end


function [r1, r2] = free_runs(mask)
%FREE_RUNS  First and last index of each run of true.
edges = diff([false, mask, false]);
r1 = find(edges == 1);
r2 = find(edges == -1) - 1;
end


function still = frames_outside(m1, m2, T)
%FRAMES_OUTSIDE  The complement of the runs, which signed_runs no longer is.
still = true(1, T);
for k = 1:numel(m1)
    still(m1(k):m2(k)) = false;
end
end


function [span_lo, span_hi, in_span] = baseline_spans(state_idx, states, T, mode)
%BASELINE_SPANS  The stretches the baseline is computed inside.
%   'record' is one span, 'bout' one per bout in precedence order. A bout that
%   overlaps a claimed stretch keeps the free part rather than being dropped whole.
if strcmp(mode, 'record')
    span_lo = 1;
    span_hi = T;
    in_span = true(1, T);
    return
end
claimed = false(1, T);
span_lo = [];
span_hi = [];
for s = 1:numel(states)
    this_name = states{s};
    if ~isfield(state_idx, this_name)
        continue
    end
    bouts = state_idx.(this_name);
    for i = 1:size(bouts, 1)
        lo = max(1, bouts(i,1));
        hi = min(T, bouts(i,2));
        if hi <= lo
            continue
        end
        free = false(1, T);
        free(lo:hi) = ~claimed(lo:hi);
        [f1, f2] = free_runs(free);
        for r = 1:numel(f1)
            claimed(f1(r):f2(r)) = true;
            span_lo(end+1) = f1(r); %#ok<AGROW>
            span_hi(end+1) = f2(r); %#ok<AGROW>
        end
    end
end
[span_lo, ord] = sort(span_lo);
span_hi = span_hi(ord);
in_span = claimed;
end


function [m1, m2, dir_run] = signed_runs(d, is_still)
%SIGNED_RUNS  Moving runs, split wherever the direction reverses inside one.
%   A run is above the gate AND of one sign. At a sharp turn two frames straddle
%   the zero crossing with the gate cleared on both sides, so the split is needed.
label = sign(d);
label(is_still) = 0;
% the frame before a reversal IS the turning point; donate it to stillness so the
% two runs each get a region to anchor on
flip_at = find(label(1:end-1) ~= 0 & label(2:end) ~= 0 & ...
    label(1:end-1) ~= label(2:end));
label(flip_at) = 0;
starts = find(label ~= 0 & [true, diff(label) ~= 0]);
m1 = starts;
m2 = zeros(1, numel(starts));
dir_run = label(starts);
for k = 1:numel(starts)
    stop = starts(k);
    while stop < numel(label) && label(stop+1) == label(starts(k))
        stop = stop + 1;
    end
    m2(k) = stop;
end
end


function keep = level_crossed(h_start, h_end, amp, mode)
%LEVEL_CROSSED  Whether a merged object qualifies, under the selected rule.
if amp <= 0
    keep = true;
    return
end
if strcmp(mode, 'zero')
    keep = (h_end >= amp) || (h_end <= -amp);
    return
end
cross_up_hi = h_start < amp && h_end >= amp;
cross_dn_hi = h_start >= amp && h_end < amp;
cross_dn_lo = h_start > -amp && h_end <= -amp;
cross_up_lo = h_start <= -amp && h_end > -amp;
keep = cross_up_hi || cross_dn_hi || cross_dn_lo || cross_up_lo;
end


function [matched, name, i_bout, b1, b2, s1, s2] = ...
    first_bout(state_idx, states, from, to, T, edge_margin)
%FIRST_BOUT  The first bout holding both landmarks inside its margin band.
matched = false;
name = '';
i_bout = 0;
b1 = 0;
b2 = 0;
s1 = 0;
s2 = 0;
for s = 1:numel(states)
    this_name = states{s};
    if ~isfield(state_idx, this_name)
        continue
    end
    bouts = state_idx.(this_name);
    for i = 1:size(bouts, 1)
        lo = max(1, bouts(i,1));
        hi = min(T, bouts(i,2));
        band_lo = lo + edge_margin;
        band_hi = hi - edge_margin;
        if from < band_lo || to > band_hi
            continue
        end
        matched = true;
        name = this_name;
        i_bout = i;
        b1 = lo;
        b2 = hi;
        s1 = band_lo;
        s2 = band_hi;
        return
    end
end
end


function counts = count_states(steplist, states)
%COUNT_STATES  Rows per state and per direction.
counts = struct();
if isempty(steplist)
    row_state = {};
    row_pol = {};
else
    row_state = {steplist.state};
    row_pol = {steplist.pol};
end
for s = 1:numel(states)
    this_name = states{s};
    in_state = strcmp(row_state, this_name);
    is_dil = in_state & strcmp(row_pol, 'dilation');
    is_con = in_state & strcmp(row_pol, 'constriction');
    counts.(this_name) = struct('n_row', nnz(in_state), ...
        'n_dilation', nnz(is_dil), 'n_constriction', nnz(is_con));
end
end
