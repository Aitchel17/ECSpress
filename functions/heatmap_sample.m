function row = heatmap_sample(bv_px, eps_px, y_series, um_per_px, sample_mask)
%HEATMAP_SAMPLE  One heatmap over a set of samples, and what is read off it.
%   Caller: heatmap_paxfwhm
%   IN   bv_px       1xT double   lumen diameter, px
%        eps_px      1xT double   outer diameter, px
%        y_series    1x1 str      which series goes on y. "eps" or "totalpvs"
%        um_per_px   1x1 double   this session's scale
%        sample_mask 1xT logical  which samples to read. [] = all of them
%   OUT  row         1x1 struct   the twenty-three fields, declared on the first line
    row = struct('n_sample', 0, 'bv_range', NaN, 'cond_sd', NaN, 'sample_slope', NaN, ...
        'sample_r', NaN, 'mode_slope', NaN, 'slope_constricted', NaN, ...
        'slope_dilated', NaN, 'n_constricted', 0, 'n_dilated', 0, ...
        'beta_ratio', NaN, 'beta_gap', NaN, 'bv_sessionmedian', NaN, 'y_sessionmedian', NaN, ...
        'eps_sessionmedian', NaN, ...
        'area_baseline_sample', NaN, ...
        'area_gap_constricted', NaN, 'area_gap_dilated', NaN, ...
        'area_gapfrac_constricted', NaN, 'area_gapfrac_dilated', NaN, ...
        'bv_heatmedian', NaN, 'baseline_bv', NaN, 'baseline_y', NaN, 'heat', []);

    % COLUMNS from here down. robustfit and corr both read a column and answer a
    % matrix, not a number, on a row
    bv_px = bv_px(:);
    eps_px = eps_px(:);
    pvs_px = eps_px - bv_px;

    % what the histogram bins on y, and the lumen it dropped to get there
    if y_series == "totalpvs"
        y_px = pvs_px;
        lumen_removed = 1;
    else
        y_px = eps_px;
        lumen_removed = 0;
    end

    % Filter NaN
    both_finite = isfinite(bv_px) & isfinite(eps_px);
    good = both_finite;
    if ~isempty(sample_mask)
        good = good & sample_mask(:);
    end
    bv_good = bv_px(good); y_good = y_px(good); eps_good = eps_px(good);
    row.n_sample = sum(good);

    [bin_counts, x_edge, y_edge] = histcounts2(bv_good, y_good, ...
        'BinWidth', [1 1], 'Normalization', 'percentage');

    % the shared zero for every state of this session: computed before the mask
    bv_sessionmedian = median(bv_px(both_finite)) * um_per_px;
    y_sessionmedian = median(y_px(both_finite)) * um_per_px;
    % y is eps only when y_series says so, and bv/eps needs eps either way
    eps_sessionmedian = median(eps_px(both_finite)) * um_per_px;

    % the map, then the three measurements, each on its own input
    heat = heatmap_postprocessing(bin_counts, x_edge * um_per_px, y_edge * um_per_px, ...
        bv_sessionmedian, y_sessionmedian);

    % the mode curve in absolute um, and the eps it implies
    bv_heatmedian = median(bv_good) * um_per_px - bv_sessionmedian;
    bv_column = bv_sessionmedian + heat.x_baseceneters;
    y_column = y_sessionmedian + heat.mode_curve;
    eps_column = y_column + lumen_removed * bv_column;

    slope = heatmap_modeslope(heat.x_baseceneters, heat.mode_curve, 2);
    area = heatmap_area(heat.x_baseceneters, bv_column, eps_column, bv_heatmedian);
    fitted = sample_fit(bv_good, y_good, eps_good);

    row.bv_sessionmedian = bv_sessionmedian;
    row.y_sessionmedian = y_sessionmedian;
    row.eps_sessionmedian = eps_sessionmedian;
    row.bv_range = range(heat.x_baseceneters);
    row.cond_sd = heat.column_std;
    row.baseline_bv = heat.x_mode;
    row.baseline_y = heat.y_mode;
    row.bv_heatmedian = bv_heatmedian;

    row.mode_slope = slope.whole;
    row.slope_constricted = slope.constricted;
    row.slope_dilated = slope.dilated;
    row.n_constricted = slope.n_constricted;
    row.n_dilated = slope.n_dilated;

    row.beta_ratio = fitted.beta_ratio;
    row.sample_slope = fitted.slope;
    row.sample_r = fitted.r;

    row.area_gap_constricted = area.gap_constricted;
    row.area_gap_dilated = area.gap_dilated;
    row.area_gapfrac_constricted = area.gapfrac_constricted;
    row.area_gapfrac_dilated = area.gapfrac_dilated;
    row.beta_gap = (fitted.beta_ratio - 1) * area.bv_baseline / area.eps_baseline;
    % from the median sample geometry, not from the curve above
    row.area_baseline_sample = (pi/4) * ((median(eps_good) * um_per_px)^2 - ...
        (median(bv_good) * um_per_px)^2);

    heat.area_curve = area.curve;
    row.heat = heat;
end

function fitted = sample_fit(bv_good, y_good, eps_good)
%SAMPLE_FIT  The fits made on every sample rather than on the map's columns.
%   Scale free, so pixels need no conversion.
%   IN   bv_good     Nx1 double  lumen diameter, the samples that survived the mask
%        y_good      Nx1 double  the series on y, same samples
%        eps_good    Nx1 double  outer diameter, same samples
%   OUT  fitted      1x1 struct  beta_ratio  1x1  slope of eps^2 against bv^2
%                                slope       1x1  slope of y against bv
%                                r           1x1  their correlation
fitted = struct('beta_ratio', NaN, 'slope', NaN, 'r', NaN);
if numel(bv_good) < 4
    return
end

area_fit = robustfit(bv_good.^2, eps_good.^2);
fitted.beta_ratio = area_fit(2);

series_fit = robustfit(bv_good, y_good);
fitted.slope = series_fit(2);
fitted.r = corr(bv_good, y_good);
end
