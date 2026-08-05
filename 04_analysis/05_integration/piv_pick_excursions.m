function [events, info] = piv_pick_excursions(dtrace, t_axis, pairs, sg_win_s, peaktol)
%PIV_PICK_EXCURSIONS  Locate the two endpoints of each paired-bout excursion.
%   piv_merge_bouts fixed which stretches to search and what polarity to expect;
%   this locates where on the trace each endpoint lands. Polarity is prescribed by
%   the transition, never measured.
%
% IN   dtrace    1xT absolute vessel diameter (um)
%      t_axis    1xT frame times (s); also sets the frame rate
%      pairs     output of piv_merge_bouts
%      sg_win_s  Savitzky-Golay window (s) for the picking trace; order is 3
%      peaktol   band level, from the window median toward its extreme
% OUT  events    1xN struct, chronological: every pairs field plus from, to,
%                rise_s, diameter_change (picking trace), diameter_change_raw
%      info      dsg, sg_frame, P, sign_ok, sign_bad, clipped, dropped

arguments
    dtrace   (1,:) double
    t_axis   (1,:) double
    pairs    struct
    sg_win_s (1,1) double {mustBePositive}
    peaktol  (1,1) double {mustBeInRange(peaktol, 0, 1)}
end
SG_ORD = 3;     % Savitzky-Golay polynomial order

% 0. Check inputs
% 0.1 Polarity labels
badpol = unique({pairs(~ismember({pairs.pol}, {'dilation','constriction'})).pol});
if ~isempty(badpol)
    error('piv_pick_excursions:badPolarity', ...
        'polarity must be dilation|constriction, got: %s', strjoin(badpol, ', '));
end
% 0.2 The two series must describe the same frames
if numel(dtrace) ~= numel(t_axis)
    error('piv_pick_excursions:lengthMismatch', ...
        'dtrace has %d samples but t_axis has %d.', numel(dtrace), numel(t_axis));
end

% 1. Build the picking trace
nT       = numel(dtrace);
fps      = 1 / mean(diff(t_axis));
sg_frame = 2*floor(sg_win_s * fps / 2) + 1;         % 1.1 frame length, odd
dsg      = sgolayfilt(dtrace, SG_ORD, sg_frame);    % 1.2 smooth

% 2. Locate one endpoint pair per merged bout pair
% 2.1 Output template and drop counters
efields = [fieldnames(pairs); ...
           {'from'; 'to'; 'rise_s'; 'diameter_change'; 'diameter_change_raw'}];
events  = cell2struct(cell(numel(efields), 0), efields, 1);
dropped = struct('empty_band', 0, 'out_of_range', 0);
for k = 1:numel(pairs)
    % 2.2 The two search windows
    Q  = pairs(k);
    wa = Q.search_from(1):Q.search_from(2);
    wb = Q.search_to(1):Q.search_to(2);
    % 2.3 Drop windows running off the trace
    if wa(1) < 1 || wb(end) > nT
        dropped.out_of_range = dropped.out_of_range + 1;  continue
    end
    % 2.4 Band each side: a dilation starts low and ends high, constriction reversed
    sA = 1 - 2 * strcmp(Q.pol, 'dilation');
    ca = wa(band_mask(dsg(wa),  sA, peaktol));
    cb = wb(band_mask(dsg(wb), -sA, peaktol));
    if isempty(ca) || isempty(cb)
        dropped.empty_band = dropped.empty_band + 1;  continue
    end
    % 2.5 Take the banded sample closest to the transition on each side
    f = ca(end);   t = cb(1);
    % 2.6 Fill the event record
    E = Q;
    E.from = f;
    E.to   = t;
    E.rise_s              = t_axis(t) - t_axis(f);
    E.diameter_change     = dsg(t)    - dsg(f);
    E.diameter_change_raw = dtrace(t) - dtrace(f);
    events(end+1) = E; %#ok<AGROW>
end
% 2.7 Sort events chronologically
if ~isempty(events)
    [~, ord] = sort([events.from]);   events = events(ord);
end

% 3. Diagnostics
sign_bad = false(1, numel(events));
clipped  = false(1, numel(events));
for e = 1:numel(events)
    isdil = strcmp(events(e).pol, 'dilation');
    % 3.1 Measured amplitude against the prescribed polarity
    sign_bad(e) = (isdil && events(e).diameter_change <= 0) || ...
                  (~isdil && events(e).diameter_change >= 0);
    % 3.2 Endpoint on a search-window edge, so the amplitude is a lower bound
    clipped(e)  = events(e).from == events(e).search_from(1) || ...
                  events(e).to   == events(e).search_to(2);
end
% 3.3 Pack
info = struct('dsg', dsg, 'sg_frame', sg_frame, ...
    'P', struct('sg_win_s', sg_win_s, 'peaktol', peaktol, 'fps', fps), ...
    'sign_ok', ~any(sign_bad), 'sign_bad', find(sign_bad), ...
    'clipped', find(clipped), 'dropped', dropped);
end

% ---------------------------------------------------------------------------
function m = band_mask(v, s, frac)
%BAND_MASK  Samples inside the extreme amplitude band. s = +1 upper, -1 lower.
%   level = median + frac*(extreme - median), an amplitude fraction, not a rank.
    v  = s * v;
    md = median(v, 'omitnan');
    m  = v >= md + frac * (max(v) - md);
end
