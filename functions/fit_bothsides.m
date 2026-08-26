function [slope_low, slope_high, info] = fit_bothsides(x_value, y_value, min_column)
%FIT_BOTHSIDES  A straight line either side of zero, on a span forced symmetric.
%   x is expected to be centred already -- on a modal diameter, an event onset,
%   whatever the caller's anchor is -- so zero is a real place and the two sides
%   mean something.
%
%   WHICH OUTPUT IS WHICH IS FIXED BY THE SIGN OF X AND NEVER BY MAGNITUDE. A
%   version that sorted the two and called the smaller one anything would be
%   zz_notinuse/tertiary_analysis/secondary_afterproceesing.m:47-51 reintroduced:
%   sorting before comparing puts every point on one side of the diagonal
%   whatever the data says.
%   see STRUCTURE.md
%
%   THE SPAN IS FORCED SYMMETRIC. The two sides rarely reach equally far, and a
%   fit over each side's whole reach compares a short span against a long one --
%   which is a difference in leverage, not in slope. The longer side is truncated.
%
% IN   x_value    1 x n double  the centred abscissa, any unit
%      y_value    1 x n double  what is fitted, any unit
%      min_column 1 x 1 int     points a side needs before it is fitted
% OUT  slope_low   1 x 1 double  fitted where x < 0. NaN if that side is too thin
%      slope_high  1 x 1 double  fitted where x > 0
%      info  struct
%        .span     1 x 1 double  the half width both sides were truncated to
%        .n_low    1 x 1 int     points the low side was fitted on
%        .n_high   1 x 1 int     points the high side was fitted on
%        .bend     1 x 1 double  slope_high - slope_low, the paired quantity
%        .intercept_low  1 x 1 double  so the fitted line can be drawn again
%        .intercept_high 1 x 1 double
    arguments
        x_value    (1,:) double
        y_value    (1,:) double
        min_column (1,1) double {mustBeInteger, mustBePositive}
    end

slope_low = NaN;
slope_high = NaN;
info = struct('span', NaN, 'n_low', 0, 'n_high', 0, 'bend', NaN, ...
    'intercept_low', NaN, 'intercept_high', NaN);

usable = isfinite(x_value) & isfinite(y_value);
if ~any(usable & x_value < 0) || ~any(usable & x_value > 0)
    return
end
reach_low = abs(min(x_value(usable)));
reach_high = max(x_value(usable));
info.span = min(reach_low, reach_high);

on_low = usable & x_value < 0 & x_value >= -info.span;
on_high = usable & x_value > 0 & x_value <= info.span;
info.n_low = sum(on_low);
info.n_high = sum(on_high);
if info.n_low < min_column || info.n_high < min_column
    return
end

coef_low = polyfit(x_value(on_low), y_value(on_low), 1);
coef_high = polyfit(x_value(on_high), y_value(on_high), 1);
slope_low = coef_low(1);
slope_high = coef_high(1);
info.intercept_low = coef_low(2);
info.intercept_high = coef_high(2);
info.bend = slope_high - slope_low;
end
