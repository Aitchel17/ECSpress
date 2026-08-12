function [idx, info] = eventdetect_findtransition(v, direction, side, frac)
%EVENTDETECT_FINDTRANSITION  Where the trace turns, on one side of a scored boundary.
%   The scored bout boundary is nominal: it says which frame the scorer called the
%   state change, not which frame the diameter started moving. This bands one
%   half-window at median + frac*(extreme - median) toward `direction` and returns
%   the banded sample on the boundary end of it, which is where the trace leaves
%   its plateau. event_pick_excursions calls it twice, once per side.
%
% IN   v          1 x N float  ONE half-window of the picking trace, already
%                              smoothed. median, min and max are all local to it
%      direction  str          "max" bands toward the window maximum, "min" toward
%                              its minimum. NOT the polarity -- a dilation wants
%                              "min" before the boundary and "max" after it, a
%                              constriction the other way round
%      side       str          "pre"  boundary at v(end), take the last banded
%                              "post" boundary at v(1),   take the first banded
%      frac       float 0..1   band level, from the window median toward its extreme
% OUT  idx        int          index INTO v; [] when the band is empty
%      info       struct       thr, n_band, band_span -- what the banding measured
%
% UNIT  whatever v was in. idx does not depend on it; info.thr and info.band_span do
%   why    median + frac*(extreme-median) is linear, so a constant scale cancels
%          out of the level it compares
%   see FINDINGS.md
%   why    the statistics are local to the half-window on purpose. One median over
%          both sides of a transition lands mid-transition, and then the lower band
%          swallows the whole starting plateau and the upper band the whole ending
%          one, so neither edge is found
%   note   the pipeline hands down analysis_event.m opt.peaktol and does not fall
%          through to the default here. Two copies of 0.80 -- do not let them drift

arguments
    v         (1,:) double
    direction (1,1) string {mustBeMember(direction, ["min" "max"])}
    side      (1,1) string {mustBeMember(side, ["pre" "post"])}
    frac      (1,1) double {mustBeInRange(frac, 0, 1)} = 0.80
end

% 1. Band the window. Flipping to always work upward lets one comparison serve
%    both directions
if direction == "min"
    sign_flip = -1;
else
    sign_flip = 1;
end
v_up    = sign_flip * v;
md      = median(v_up, 'omitnan');
peak    = max(v_up);
thr     = md + frac * (peak - md);
in_band = v_up >= thr;
banded  = find(in_band);

% 2. The banded sample on the boundary end of the window
if isempty(banded)
    idx = [];
elseif side == "pre"
    idx = banded(end);
else
    idx = banded(1);
end

% 3. What the banding measured
info = struct('thr', sign_flip * thr, 'n_band', numel(banded), 'band_span', NaN);
if ~isempty(banded)
    v_banded = v(banded);
    info.band_span = max(v_banded) - min(v_banded);
end
end
