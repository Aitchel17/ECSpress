function slope = heatmap_modeslope(x_base, mode_curve, min_column)
%HEATMAP_MODESLOPE  The mode curve's slope, whole and either side of zero.
%   IN   x_base      1 x N double  column centres, on the anchored axis
%        mode_curve  1 x N double  each column's most-occupied row, same axis
%        min_column  1x1 double    columns a one-sided fit needs
%   OUT  slope       1x1 struct    whole          1x1  robustfit over the whole curve
%                                  constricted    1x1  fitted where x < 0
%                                  dilated        1x1  where x > 0, on the same span
%                                  n_constricted  1x1  columns that fit was made on
%                                  n_dilated      1x1
%                                  bend           1x1  dilated minus constricted
slope = struct('whole', NaN, 'constricted', NaN, 'dilated', NaN, ...
    'n_constricted', 0, 'n_dilated', 0, 'bend', NaN);

if numel(x_base) >= 4
    whole_fit = robustfit(x_base, mode_curve);
    slope.whole = whole_fit(2);
end

[slope_low, slope_high, side_info] = fit_bothsides(x_base, mode_curve, min_column);
slope.constricted = slope_low;
slope.dilated = slope_high;
slope.n_constricted = side_info.n_low;
slope.n_dilated = side_info.n_high;
slope.bend = side_info.bend;
end
