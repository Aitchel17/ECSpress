%FIG_FWHM_EPSBV  Every vessel's eps-against-bv curve, each against its OWN no-flux curve.
%   One line per vessel on the shared pooling grid. Area conservation is
%   d(eps)/d(bv) = bv/eps, that vessel's own geometry, so each line gets its own
%   reference. Block 3 subtracts it: zero is area conserved, the panel a
%   vessel-level claim reads off. Nothing is saved unless block 4 is uncommented.

%% 1. Read
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');
dirs_ecspath;
clee = color_lee();

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
param.column_stat = "median";      % Which stored curve to draw (str), mode/centroid/median
param.y_show = "totalpvs";         % Axis to draw on (str), "totalpvs" or "eps"
param.fontsize = 11;
param.min_column = 4;              % Columns a vessel needs before its line is drawn (count)

dirs = dirs_central();
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);
relation = load(fullfile(dirs.save_dir, 'fwhmrelation.mat')).save_content;
vessel = relation.vesselbend;
grid_x = relation.grid_x;
reach = relation.vessel_pool.vessel_count > 0;      % V x N, columns that vessel touched
n_vessel = height(vessel);
fprintf('%d vessels | grid %+.2f to %+.2f um, %d columns\n', ...
    n_vessel, min(grid_x), max(grid_x), numel(grid_x));

%% 2. Every vessel, measured against its own no-flux curve
% each grey curve is the level curve of (pi/4)(eps^2 - bv^2) through that vessel's
% own rest point; a line above its grey is an annulus that gained area
[measured, reference, drawn, eps_curve] = vessel_curves(vessel, grid_x, reach, param);
fprintf('%d of %d vessels reach %d columns\n', sum(drawn), n_vessel, param.min_column);

fig_curve = figure('Color', 'w', 'Position', [80 80 620 560], 'Name', 'fig_fwhm_epsbv_curve');
ax = axes(fig_curve);
hold(ax, 'on')
for v = find(drawn)'
    plot(ax, grid_x, reference(v, :), '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.8);
end
for v = find(drawn)'
    plot(ax, grid_x, measured(v, :), '-', 'Color', clee.clist.nrem, 'LineWidth', 1);
end
xline(ax, 0, '-', 'Color', [0.6 0.6 0.6]);
xlabel(ax, 'lumen diameter, um from that vessel''s own median')
ylabel(ax, y_label(param.y_show))
set(ax, 'FontSize', param.fontsize, 'Box', 'on')
title(ax, sprintf('%d vessels, %s curve | grey = that vessel''s own no-flux', ...
    sum(drawn), param.column_stat), 'FontWeight', 'normal')

%% 3. The same curves with each vessel's own reference taken out
% zero is area conserved; above it the annulus gained area, below it gave area up.
% The same number on either axis : y_show subtracts x from curve and reference alike
departure = measured - reference;

fig_gap = figure('Color', 'w', 'Position', [120 120 620 560], 'Name', 'fig_fwhm_epsbv_gap');
ax = axes(fig_gap);
hold(ax, 'on')
for v = find(drawn)'
    plot(ax, grid_x, departure(v, :), '-', 'Color', [0.75 0.80 0.86], 'LineWidth', 0.8);
end
% the median over the vessels that reached each column, and how many that was
column_n = sum(isfinite(departure(drawn, :)), 1);
column_median = median(departure(drawn, :), 1, 'omitnan');
enough = column_n >= param.min_column;
plot(ax, grid_x(enough), column_median(enough), '-', 'Color', clee.clist.nrem, ...
    'LineWidth', 2.4);
yline(ax, 0, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
xline(ax, 0, '-', 'Color', [0.6 0.6 0.6]);
xlabel(ax, 'lumen diameter, um from that vessel''s own median')
ylabel(ax, 'outer PVS diameter, um from where no flux would put it')
set(ax, 'FontSize', param.fontsize, 'Box', 'on')
title(ax, sprintf('zero = area conserved | thick line is the median of %d vessels', ...
    sum(drawn)), 'FontWeight', 'normal')

fprintf('%8s %8s %10s %10s\n', 'x um', 'vessels', 'median', 'IQR');
for x_at = [-2 -1 0 1 2 3 4]
    [~, k] = min(abs(grid_x - x_at));
    one = departure(drawn, k);
    one = one(isfinite(one));
    if numel(one) < param.min_column
        continue
    end
    quart = quantile(one, [0.25 0.75]);
    fprintf('%+8.1f %8d %10.3f %5.2f..%.2f\n', grid_x(k), numel(one), median(one), ...
        quart(1), quart(2));
end

%% 3b. The same departure in area, where the squaring shows
% A = (pi/4)(eps^2 - bv^2) is flat at zero under area conservation, so no reference is
% subtracted. eps^2 carries 2*eps times the boundary error, so the fraction is printed beside it
bv_at = vessel.bv_vesselmedian + grid_x;
eps_at = vessel.eps_vesselmedian + eps_curve;
area_at = (pi/4) * (eps_at.^2 - bv_at.^2);
area_rest = (pi/4) * (vessel.eps_vesselmedian.^2 - vessel.bv_vesselmedian.^2);
area_gap = area_at - area_rest;

fig_area = figure('Color', 'w', 'Position', [160 160 620 560], 'Name', 'fig_fwhm_epsbv_area');
ax = axes(fig_area);
hold(ax, 'on')
for v = find(drawn)'
    plot(ax, grid_x, area_gap(v, :), '-', 'Color', [0.86 0.80 0.76], 'LineWidth', 0.8);
end
area_n = sum(isfinite(area_gap(drawn, :)), 1);
area_med = median(area_gap(drawn, :), 1, 'omitnan');
enough_area = area_n >= param.min_column;
plot(ax, grid_x(enough_area), area_med(enough_area), '-', 'Color', clee.clist.awake, ...
    'LineWidth', 2.4);
yline(ax, 0, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
xline(ax, 0, '-', 'Color', [0.6 0.6 0.6]);
xlabel(ax, 'lumen diameter, um from that vessel''s own median')
ylabel(ax, 'annulus area, um^2 from its resting value')
set(ax, 'FontSize', param.fontsize, 'Box', 'on')
title(ax, sprintf('zero = area conserved | %d vessels', sum(drawn)), 'FontWeight', 'normal')

% how much of the spread is the vessel size doing the multiplying
fprintf('\n%8s %8s %12s %14s %12s %12s\n', 'x um', 'vessels', 'median um^2', ...
    'IQR um^2', 'IQR width', 'IQR frac');
for x_at = [-2 -1 0 1 2 3 4]
    [~, k] = min(abs(grid_x - x_at));
    one = area_gap(drawn, k);
    frac = area_gap(drawn, k) ./ area_rest(drawn);
    one = one(isfinite(one));
    frac = frac(isfinite(frac));
    if numel(one) < param.min_column
        continue
    end
    q = quantile(one, [0.25 0.75]);
    qf = quantile(frac, [0.25 0.75]);
    fprintf('%+8.1f %8d %12.2f %6.1f..%-6.1f %12.2f %12.4f\n', grid_x(k), numel(one), ...
        median(one), q(1), q(2), q(2) - q(1), qf(2) - qf(1));
end

%% 4. Save, when a panel above is worth keeping
% save_figure(fig_curve, dirs.save_dir, "fig_fwhm_epsbv_curve")
% save_figure(fig_gap,   dirs.save_dir, "fig_fwhm_epsbv_gap")
% save_figure(fig_area,  dirs.save_dir, "fig_fwhm_epsbv_area")

%% ---------------------------------------------------------------- helpers
function [measured, reference, drawn, eps_curve] = vessel_curves(vessel, grid_x, reach, param)
%VESSEL_CURVES  Each vessel's stored curve and the level curve through its own rest.
%   IN   vessel     Vx1 table    vesselbend, holding curve_* and the two medians
%        grid_x     1xN double   the shared abscissa, um from each vessel's own median
%        reach      VxN logical  columns that vessel's sessions actually touched
%        param      1x1 struct   column_stat, y_show, min_column
%   OUT  measured   VxN double   NaN where that vessel did not reach
%        reference  VxN double   the same, area conserved through its own rest
%        drawn      Vx1 logical  vessels with enough columns to draw
%        eps_curve  VxN double   the same curve on the eps axis, whatever y_show says
    stored = vessel.("curve_" + param.column_stat);
    stored(~reach) = NaN;

    bv_rest = vessel.bv_vesselmedian;
    eps_rest = vessel.eps_vesselmedian;
    bv_value = bv_rest + grid_x;
    eps_noflux = sqrt(eps_rest.^2 + bv_value.^2 - bv_rest.^2);

    eps_curve = stored;
    measured = stored;
    reference = eps_noflux - eps_rest;
    if param.y_show == "totalpvs"
        measured = measured - grid_x;
        reference = reference - grid_x;
    end
    reference(~isfinite(measured)) = NaN;
    drawn = sum(isfinite(measured), 2) >= param.min_column;
end

function label = y_label(y_show)
%Y_LABEL  What the ordinate is, in the unit the block drew it in.
    if y_show == "totalpvs"
        label = "PVS thickness, um from that vessel's own median";
    else
        label = "outer PVS diameter, um from that vessel's own median";
    end
end

function save_figure(fig, save_dir, fig_name) %#ok<DEFNU>
% R2023b exportgraphics cannot write svg, so the vector copy goes through print.
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    print(fig, fullfile(save_dir, fig_name + ".svg"), '-dsvg', '-vector');
    exportgraphics(fig, fullfile(save_dir, fig_name + ".png"), 'Resolution', 200);
    fprintf('wrote %s\n', fullfile(save_dir, fig_name));
end
