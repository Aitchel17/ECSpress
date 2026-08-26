%PVS_DIAMETER_FIG  PVS extent against lumen diameter, as a 2-D heatmap.
%   Reads fwhmrelation.mat and puts it on axes. Every heatmap here is drawn by
%   plot_heatpanel; this file decides which rows go on which tile and what each
%   panel is labelled with.
%
%   Nothing is measured below. tablegeneration_fwhmrelation built the densities,
%   the per-session numbers, the pooled map and the per-vessel bend, so changing a
%   colour or a tile count costs a redraw and changing what is measured is a run
%   of that file.
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
clee = color_lee();

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
dirs = dirs_central();
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);

% How the rows are drawn. None of these changes the sample: that was settled in
% tablegeneration_fwhmrelation and travels here inside relation.param.
param.colormap = clee.gradient.inferno;   % a heatmap, so a perceptually ordered ramp
param.fontsize = 11;
param.n_tile = Inf;        % sessions the grid shows, widest diameter range first.
param.n_col = 10;          % columns in that grid

%% What is drawn
relation_path = fullfile(dirs.save_dir, 'fwhmrelation.mat');
relation = load(relation_path).save_content;
%%
session = relation.session;
pooled = relation.pooled;
vesselbend = relation.vesselbend;

drawn = session(session.drawn, :);

% The two references, and neither is a number typed in. Both come out of the same
% expression eps^2 = beta*bv^2 + c: beta = 1 IS the area-conserving hypothesis, and
% beta = the table's fitted beta_ratio is what these sessions did. They are CURVES,
% built per panel because each panel has its own abscissa and its own operating
% point, and because the area-conserving reference is a level curve whose slope
% runs 0.42 to 0.65 across the diameters one panel spans. see FINDINGS.md
ref_name = ["no flux, area conserved", "measured"];
ref_beta = [1, median(drawn.beta_ratio)];
ref_r = [NaN, median(drawn.sample_r)];

ref_noflux = find(ref_name == "no flux, area conserved", 1);
ref_measured = find(ref_name == "measured", 1);
y_label = relation.param.y_series + ", um from its mode";

fprintf('fwhmrelation %s\n   %d sessions, %d drawn, %d vessels pooled\n', ...
    relation_path, height(session), height(drawn), height(vesselbend));

%% One session, with its two marginals
[~, pick] = max(drawn.bv_range);
one = drawn(pick, :);
heat = one.heat{1};
fprintf('\nrepresentative %s\n   %d samples | %.3f um/px | diameter range %.2f um\n', ...
    string(one.Directory), one.n_sample, one.NumericResolution, one.bv_range);

fig_one = figure('Color', 'w', 'Position', [80 60 980 780], 'Name', 'pvs_diameter_heatmap');
layout_one = tiledlayout(fig_one, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

axis_main = nexttile(layout_one, 5, [3 3]);
[ref_main, slope_main] = reference_curve(heat.x_baseceneters, one.bv_anchor, ...
    one.pvs_anchor, ref_beta, relation.param.y_series);
ref_label_main = ref_name + label_reference(slope_main, ref_r);
handle_main = plot_heatpanel(axis_main, heat.x_baseceneters, heat.y_baseceneters, ...
    heat.xy_counts_clean, 'curve', heat.modepvs(:)', 'ref_curve', ref_main, ...
    'colormap', param.colormap, 'fontsize', param.fontsize);
xlabel(axis_main, 'lumen diameter, um from this session''s mode')
ylabel(axis_main, y_label)
% the fit is labelled where it is drawn, so a panel can be read without the console
mode_label = sprintf('mode per diameter bin | sample fit %+.3f, r %+.2f', ...
    one.sample_slope, one.sample_r);
legend(axis_main, [handle_main.curve, handle_main.ref], [string(mode_label), ref_label_main], ...
    'Location', 'southwest', 'TextColor', 'w', 'Color', [0 0 0], 'Box', 'off', ...
    'FontSize', param.fontsize - 2)

axis_top = nexttile(layout_one, 1, [1 3]);
bar(axis_top, heat.x_baseceneters, sum(heat.xy_counts_clean, 1, 'omitmissing'), 1, ...
    'FaceColor', [0.2 0.6 1], 'EdgeColor', 'none');
ylabel(axis_top, '% of samples')
title(axis_top, 'lumen diameter')
set(axis_top, 'XTickLabel', [], 'FontSize', param.fontsize - 1, 'Box', 'off')
xlim(axis_top, [min(heat.x_baseceneters) max(heat.x_baseceneters)])

axis_right = nexttile(layout_one, 8, [3 1]);
barh(axis_right, heat.y_baseceneters, sum(heat.xy_counts_clean, 2, 'omitmissing'), 1, ...
    'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none');
xlabel(axis_right, '% of samples')
title(axis_right, relation.param.y_series)
set(axis_right, 'YTickLabel', [], 'FontSize', param.fontsize - 1, 'Box', 'off')
ylim(axis_right, [min(heat.y_baseceneters) max(heat.y_baseceneters)])

title(layout_one, sprintf('%s %s %s   |   %d samples, %.3f um/px', string(one.MouseID), ...
    string(one.Date), string(one.VesselID), one.n_sample, one.NumericResolution))
save_figure(fig_one, dirs.save_dir, 'pvs_diameter_heatmap');

%% Many sessions, each on its own origin
tiled = sortrows(drawn, 'bv_range', 'descend');
n_show = min(param.n_tile, height(tiled));
tiled = tiled(1:n_show, :);
n_grid_row = ceil(n_show / param.n_col);
fig_grid = figure('Color', 'w', 'Name', 'pvs_diameter_grid', ...
    'Position', [20 20 220 * param.n_col, 190 * n_grid_row + 90]);
layout_grid = tiledlayout(fig_grid, n_grid_row, param.n_col, ...
    'TileSpacing', 'tight', 'Padding', 'compact');
fprintf('grid : %d sessions in %d x %d\n', n_show, n_grid_row, param.n_col);
for k = 1:n_show
    tile_heat = tiled.heat{k};
    axis_tile = nexttile(layout_grid);
    ref_tile = reference_curve(tile_heat.x_baseceneters, tiled.bv_anchor(k), ...
        tiled.pvs_anchor(k), ref_beta([ref_measured ref_noflux]), relation.param.y_series);
    plot_heatpanel(axis_tile, tile_heat.x_baseceneters, tile_heat.y_baseceneters, ...
        tile_heat.xy_counts_clean, 'curve', tile_heat.modepvs(:)', ...
        'ref_curve', ref_tile, ...
        'ref_style', ["--", ":"], 'ref_color', [1 0.85 0.3; 0.6 0.8 1], ...
        'colormap', param.colormap, 'curve_width', 1.5, 'fontsize', 6, 'ticks', false);
    % the slope and the scatter, because those are the two things the eye is judging
    title(axis_tile, sprintf('%s %s %s  %.2f r%.2f  sd%.2f', string(tiled.MouseID(k)), ...
        string(tiled.Date(k)), string(tiled.VesselID(k)), tiled.sample_slope(k), ...
        tiled.sample_r(k), tiled.cond_sd(k)), 'FontSize', 6)
    axis(axis_tile, 'tight')
end
grid_legend = sprintf(['white = mode per bin | yellow dashed = %s, beta %.3f | ' ...
    'blue dotted = %s, beta %.0f'], ref_name(ref_measured), ref_beta(ref_measured), ...
    ref_name(ref_noflux), ref_beta(ref_noflux));
title(layout_grid, [char(relation.param.y_series) ...
    ' against lumen diameter, one panel per session, each on its own mode' ...
    newline grid_legend ...
    newline 'beside each label: the sample-level fitted slope, its r, and the PVS scatter at a fixed diameter. Ticks are dropped; every panel spans its own range'])
save_figure(fig_grid, dirs.save_dir, 'pvs_diameter_grid');

%% Pooled over vessels
kept_x = pooled.grid_x(~pooled.thin_column);
fig_pool = figure('Color', 'w', 'Position', [70 50 900 820], 'Name', 'pvs_diameter_pooled');
layout_pool = tiledlayout(fig_pool, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

axis_pool = nexttile(layout_pool, 1, [3 1]);
% the pooled map is an average of maps that were each shifted onto their own
% operating point, so its anchor is the median of theirs
in_pool = session(session.in_pool, :);
[ref_pool, slope_pool] = reference_curve(pooled.grid_x, median(in_pool.bv_anchor), ...
    median(in_pool.pvs_anchor), ref_beta, relation.param.y_series);
% the label belongs to the curve it sits on: these are the POOLED anchors, and the
% main panel's were one session's. Reusing one label for both put that session's
% numbers under the pooled curves
ref_label_pool = ref_name + label_reference(slope_pool, ref_r);
handle_pool = plot_heatpanel(axis_pool, pooled.grid_x, pooled.grid_y, pooled.pooled_map, ...
    'curve', pooled.pooled_mode, 'ref_curve', ref_pool, 'ref_style', ["--", "-"], ...
    'colormap', param.colormap, 'fontsize', param.fontsize);

% the two sides of the modal diameter, each fitted over the SAME half width so the
% pair is one measurement made twice. see fit_bothsides
on_constricted = pooled.grid_x < 0 & pooled.grid_x >= -pooled.sym_span;
on_dilated = pooled.grid_x > 0 & pooled.grid_x <= pooled.sym_span;
fit_handle = plot(axis_pool, pooled.grid_x(on_constricted), ...
    polyval([pooled.sym_constricted, pooled.sym_constricted_intercept], pooled.grid_x(on_constricted)), ...
    '-', 'Color', [0.3 1 0.5], 'LineWidth', 2.5);
plot(axis_pool, pooled.grid_x(on_dilated), ...
    polyval([pooled.sym_dilated, pooled.sym_dilated_intercept], pooled.grid_x(on_dilated)), ...
    '-', 'Color', [0.3 1 0.5], 'LineWidth', 2.5);

ylabel(axis_pool, y_label)
xlabel(axis_pool, 'lumen diameter, um from each session''s own mode')
fit_label = sprintf('either side of the mode: %+.3f (r %+.2f) | %+.3f (r %+.2f)', ...
    pooled.sym_constricted, pooled.r_sym_constricted, pooled.sym_dilated, pooled.r_sym_dilated);
legend(axis_pool, [handle_pool.curve, fit_handle, handle_pool.ref], ...
    [relation.param.column_stat + " of the pooled column", string(fit_label), ref_label_pool], ...
    'Location', 'southwest', 'TextColor', 'w', 'Color', [0 0 0], 'Box', 'off', ...
    'FontSize', param.fontsize - 2)
set(axis_pool, 'XTickLabel', [])
xlim(axis_pool, [min(kept_x) max(kept_x)])

axis_count = nexttile(layout_pool, 4, [1 1]);
bar(axis_count, pooled.grid_x, pooled.vessel_reaches, 1, ...
    'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none');
hold(axis_count, 'on')
yline(axis_count, relation.param.min_vessel_n, '-', 'drawn above this', 'LineWidth', 1);
xlabel(axis_count, 'lumen diameter, um from each session''s own mode')
ylabel(axis_count, 'vessels')
set(axis_count, 'FontSize', param.fontsize - 1, 'Box', 'off')
xlim(axis_count, [min(kept_x) max(kept_x)])

title(layout_pool, [char(relation.param.y_series + " given lumen diameter, pooled over vessels") ...
    newline 'each diameter column is a distribution that sums to one, averaged within vessel then across'])
save_figure(fig_pool, dirs.save_dir, 'pvs_diameter_pooled');

%% The bend, one point per vessel
fitted = vesselbend(vesselbend.fitted, :);
vessel_constricted = fitted.("constricted_" + relation.param.column_stat);
vessel_dilated = fitted.("dilated_" + relation.param.column_stat);
vessel_bend = vessel_dilated - vessel_constricted;

fig_bend = figure('Color', 'w', 'Position', [90 70 940 440], 'Name', 'pvs_diameter_bend');
layout_bend = tiledlayout(fig_bend, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axis_pair = nexttile(layout_bend);
hold(axis_pair, 'on')
axis_limit = [min([vessel_constricted; vessel_dilated]) - 0.1, max([vessel_constricted; vessel_dilated]) + 0.1];
plot(axis_pair, axis_limit, axis_limit, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
% here the axes ARE slopes, so the reference is the level curve slope at the
% operating point rather than the curve itself
xline(axis_pair, slope_pool(ref_noflux), '--', 'no flux', 'Color', [0.4 0.4 0.4]);
yline(axis_pair, slope_pool(ref_noflux), '--', 'no flux', 'Color', [0.4 0.4 0.4]);
scatter(axis_pair, vessel_constricted, vessel_dilated, 42, 'MarkerFaceColor', clee.clist.nrem, ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.75);
xlabel(axis_pair, 'slope on the constricted side')
ylabel(axis_pair, 'slope on the dilated side')
title(axis_pair, 'one point per vessel')
axis(axis_pair, 'equal')
xlim(axis_pair, axis_limit)
ylim(axis_pair, axis_limit)
set(axis_pair, 'FontSize', param.fontsize, 'Box', 'on')
axis_pair.Toolbar.Visible = 'off';

axis_hist = nexttile(layout_bend);
histogram(axis_hist, vessel_bend, 'BinWidth', 0.1, 'FaceColor', clee.clist.nrem, ...
    'EdgeColor', 'none');
xline(axis_hist, 0, '-', 'no bend', 'LineWidth', 1.5);
xline(axis_hist, median(vessel_bend), '--', sprintf('median %+.3f', median(vessel_bend)), ...
    'Color', clee.clist.rem, 'LineWidth', 1.5);
xlabel(axis_hist, 'dilated slope minus constricted slope')
ylabel(axis_hist, 'vessels')
title(axis_hist, sprintf('%d of %d above zero', sum(vessel_bend > 0), numel(vessel_bend)))
set(axis_hist, 'FontSize', param.fontsize, 'Box', 'off')
axis_hist.Toolbar.Visible = 'off';

title(layout_bend, ['the PVS gives up less per micrometre once the vessel is already wide' ...
    newline 'each vessel fitted on its own span, equal length either side of its modal diameter'])
save_figure(fig_bend, dirs.save_dir, 'pvs_diameter_bend');

%% ---------------------------------------------------------------- helpers
function save_figure(fig, save_dir, fig_name)
    % R2023b exportgraphics cannot write svg, so the vector copy goes through
    % print, which is what make_fig.save2svg does. see CLAUDE.md
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    print(fig, fullfile(save_dir, fig_name + ".svg"), '-dsvg', '-vector');
    exportgraphics(fig, fullfile(save_dir, fig_name + ".png"), 'Resolution', 200);
    fprintf('wrote %s.svg and .png to %s\n', fig_name, save_dir);
end

function [curve, slope_at_rest] = reference_curve(x_value, bv_rest, y_rest, beta_list, y_series)
%REFERENCE_CURVE  eps^2 = beta*bv^2 + c, anchored so it passes through the origin.
%   One row per beta. beta = 1 is the area-conserving level curve; the fitted
%   beta_ratio is the measured one. Anchoring at the operating point rather than at
%   the least-squares intercept is what puts every reference through (0,0), which is
%   what plot_heatpanel draws and what the panel's zero means.
% IN   x_value   1xN double   the shifted abscissa, um from the operating point
%      bv_rest   1x1 double   lumen diameter there, um
%      y_rest    1x1 double   the plotted series there, um
%      beta_list 1xR double   one reference per entry
%      y_series  1x1 str      "eps" | "totalpvs"
% OUT  curve         RxN double  on the plotted axis, zero at x = 0
%      slope_at_rest 1xR double  d(y)/d(bv) there, for a legend or an axis line
    if y_series == "totalpvs"
        eps_rest = bv_rest + y_rest;
    else
        eps_rest = y_rest;
    end
    bv_value = bv_rest + x_value;
    curve = zeros(numel(beta_list), numel(x_value));
    slope_at_rest = zeros(1, numel(beta_list));
    for r = 1:numel(beta_list)
        eps_value = sqrt(eps_rest^2 + beta_list(r) * (bv_value.^2 - bv_rest^2));
        slope_at_rest(r) = beta_list(r) * bv_rest / eps_rest;
        if y_series == "totalpvs"
            curve(r, :) = (eps_value - bv_value) - y_rest;
            slope_at_rest(r) = slope_at_rest(r) - 1;
        else
            curve(r, :) = eps_value - eps_rest;
        end
    end
end

function label = label_reference(slope_list, r_list)
    % what is appended to each reference name: its slope always, and its correlation
    % only where one exists. A hypothesis has no r, and printing one would invite the
    % reader to compare a measured tightness against a number nobody measured.
    label = strings(1, numel(slope_list));
    for k = 1:numel(slope_list)
        label(k) = sprintf(' %+.3f', slope_list(k));
        if isfinite(r_list(k))
            label(k) = label(k) + sprintf(', r %+.2f', r_list(k));
        end
    end
end
