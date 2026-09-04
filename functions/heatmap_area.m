function area = heatmap_area(x_base, bv_column, eps_column, baseline_at)
%HEATMAP_AREA  The annulus area along the mode curve, against its baseline column.
%   Caller: tablegeneration_fwhmrelation.sample_heatmap
%   IN   x_base       1 x N double  column centres, on the anchored axis
%        bv_column    1 x N double  lumen diameter at each column, um, absolute
%        eps_column   1 x N double  outer diameter at each column, um, absolute
%        baseline_at  1x1 double    the baseline column, on the anchored axis
%   OUT  area         1x1 struct    curve         1 x N  (pi/4)(eps^2 - bv^2), um^2
%                                   baseline      1x1  curve at baseline_at, um^2
%                                   bv_baseline   1x1  lumen diameter there, um
%                                   eps_baseline  1x1  outer diameter there, um
%                                   gap_constricted      1x1  mean departure, x < baseline_at
%                                   gap_dilated          1x1  x > baseline_at
%                                   gapfrac_constricted  1x1  the same over baseline
%                                   gapfrac_dilated      1x1
area.curve = (pi/4) * (eps_column.^2 - bv_column.^2);
area.baseline = interp1(x_base, area.curve, baseline_at);
area.bv_baseline = interp1(x_base, bv_column, baseline_at);
area.eps_baseline = sqrt(area.bv_baseline^2 + 4 * area.baseline / pi);

departure = area.curve - area.baseline;
on_constricted = x_base < baseline_at;
on_dilated = x_base > baseline_at;
area.gap_constricted = mean(departure(on_constricted), "omitnan");
area.gap_dilated = mean(departure(on_dilated), "omitnan");
area.gapfrac_constricted = area.gap_constricted / area.baseline;
area.gapfrac_dilated = area.gap_dilated / area.baseline;
end
