%PVS_DIAMETER_FIG  PVS extent against lumen diameter, as a 2-D density.

clc, clear
addpath('g:\03_program\01_ecspress\functions');   % where util_ecspath lives
util_ecspath;                                     % three roots, minus zz_notinuse
clee = color_lee();

%% Which dataset, and where it lives
% Set by integrate_analysisresult cell 0. The fallback is what makes this script
% runnable on its own; nothing below this block addresses the dataset.
param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
dirs = util_centraldirs();
dirs.central = fullfile(dirs.secondary_root, 'centralized');
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);

%% What gets in

param.depth_thr = 70;      % um       L1 / L2-3 boundary, the stratum every result uses
param.min_frame = 900;     % frames   fewer than this does not fill a 2-D histogram
param.min_bv_range = 2;    % um       diameter span a session must cover to be fitted
param.max_cond_sd = Inf;   % um       PVS scatter at a fixed diameter. Inf = no gate

param.pool_resolution = 0.285;   % um/px  a session is pooled only at this setting
param.resolution_tol = 0.01;     % um/px  how close counts as the same setting

param.min_vessel_n = 10;   % vessels a diameter column needs before it is drawn
param.min_column = 4;      % columns a side needs before that vessel is fitted

%% Parameter section
% How the surviving rows are reduced and drawn. Changing any of these changes the
% picture, never the sample.
param.pool_step = [];      % um       common grid spacing, [] = the coarsest pixel
param.column_stat = "median";
param.y_series = "eps";

param.colormap = clee.gradient.inferno;   % density, so a perceptually ordered ramp
param.stat_name = ["mode", "centroid", "median"];
param.fontsize = 11;
param.n_tile = Inf;        % sessions the grid shows, widest diameter range first.
                           % Inf draws every one, which is the point of the grid
param.n_col = 10;          % columns in that grid


param.ref_slope = [0.550, 0.434];
param.ref_name = [ "no flux, area conserved", ...
    "measured"];

param.ref_r = [NaN, 0.612];
if param.y_series == "totalpvs"
    param.ref_slope = param.ref_slope - 1;
    param.ref_r = [NaN, -0.749];
end

ref_noflux = find(param.ref_name == "no flux, area conserved", 1);
ref_measured = find(param.ref_name == "measured", 1);
param.y_label = param.y_series + ", um from its mode";

central = load(fullfile(dirs.central, 'centralized_paxfwhm.mat')).central;
fprintf('centralized_paxfwhm %d sessions\n', height(central));

%% Which rows the result is about
vessel_id = string(central.VesselID);
depth_um = lead_number(central.Depth);
resolution_um = lead_number(central.Resolution);
is_artery = contains(vessel_id, "PA", 'IgnoreCase', true);
is_layer1 = depth_um < param.depth_thr;
has_resolution = isfinite(resolution_um);
keep_session = is_artery & is_layer1 & has_resolution;
fprintf('centralized_paxfwhm : %d of %d rows kept\n', ...
    sum(keep_session), numel(keep_session));
fprintf('   PA %d | depth < %d um %d | resolution known %d\n', ...
    sum(is_artery), param.depth_thr, sum(is_layer1), sum(has_resolution));
row_list = find(keep_session);

%% One density per session, and the two per-session numbers the panels are labelled with
n_session = numel(row_list);
heat = cell(n_session, 1);
drawn = false(n_session, 1);
wide_enough = false(n_session, 1);
n_frame = zeros(n_session, 1);
bv_range = zeros(n_session, 1);
cond_sd = nan(n_session, 1);
% The frame-level fit, which is what the density is a picture OF. The mode-curve
% slope below is fitted to one point per column and answers a different question,
% so both are carried and both are labelled.
%   note  slope needs no unit conversion -- it is a ratio of two thicknesses on the
%        same pixel axis, so px/px equals um/um
frame_slope = nan(n_session, 1);
frame_r = nan(n_session, 1);
for k = 1:n_session
    pax_fwhm = central.data{row_list(k)};
    if ~isfield(pax_fwhm, 'thickness')
        continue
    end
    has_series = isfield(pax_fwhm.thickness, 'bv') && isfield(pax_fwhm.thickness, param.y_series);
    if ~has_series
        continue
    end
    bv_px = double(pax_fwhm.thickness.bv(:));
    pvs_px = double(pax_fwhm.thickness.(param.y_series)(:));
    good = isfinite(bv_px) & isfinite(pvs_px);
    n_frame(k) = sum(good);
    if n_frame(k) < param.min_frame
        continue
    end
    % Bin in PIXELS, then convert the axes to micrometres afterwards. Both series
    % are differences of integer row indices, so they take only integer values;
    % binning at the micrometre width lets a bin edge fall between two allowed
    % values and leaves empty rows and columns ruled across the density. One-pixel
    % bins catch one value each. The axes below are still micrometres.
    raw_heat = xy2heatmap(bv_px(good), pvs_px(good), 1);
    this_heat = heatmap_postprocessing(raw_heat);
    um_per_px = resolution_um(row_list(k));
    this_heat.x_baseceneters = this_heat.x_baseceneters * um_per_px;
    this_heat.y_baseceneters = this_heat.y_baseceneters * um_per_px;
    this_heat.modepvs = this_heat.modepvs * um_per_px;
    heat{k} = this_heat;
    bv_range(k) = range(this_heat.x_baseceneters);

    % the PVS scatter left once each frame's own diameter bin is accounted for
    bin_edge = min(bv_px(good)) - 0.5 : 1 : max(bv_px(good)) + 0.5;
    which_bin = discretize(bv_px(good), bin_edge);
    bin_median = accumarray(which_bin, pvs_px(good), [], @median, NaN);
    residual_px = pvs_px(good) - bin_median(which_bin);
    cond_sd(k) = std(residual_px, 'omitnan') * um_per_px;

    frame_fit = robustfit(bv_px(good), pvs_px(good));
    frame_slope(k) = frame_fit(2);
    frame_r(k) = corr(bv_px(good), pvs_px(good));

    wide_enough(k) = bv_range(k) >= param.min_bv_range;
    drawn(k) = wide_enough(k) && cond_sd(k) <= param.max_cond_sd;
end
fprintf('   %d had fewer than %d frames | %d spanned less than %g um\n', ...
    sum(n_frame < param.min_frame), param.min_frame, ...
    sum(n_frame >= param.min_frame & ~wide_enough & bv_range > 0), param.min_bv_range);
fprintf('   %d sessions drawn, %d held back for PVS scatter over %g um at a fixed diameter\n', ...
    sum(drawn), sum(wide_enough & cond_sd > param.max_cond_sd), param.max_cond_sd);

%% One slope per session off the mode curve
% Split from the loop above because it costs a different amount: that one builds a
% density per session, this one fits a line to a handful of points. Anything that
% only re-labels the panels belongs here, where it can be re-run for free.
%   caution  the mode of a diameter bin holding three frames is noise, so this
%        slope is dragged toward zero wherever the density is thin. It is a LABEL
%        on the panel, not the leg's estimate of d(totalpvs)/d(bv) -- that comes
%        off the state summary by regression. see FINDINGS.md
% Fitted for every session that spans enough diameter, INCLUDING the ones the
% scatter gate holds back, so the two populations can be printed against each other.
mode_slope_of = nan(n_session, 1);
for k = 1:n_session
    if ~wide_enough(k)
        continue
    end
    fit_coef = robustfit(heat{k}.x_baseceneters(:), heat{k}.modepvs(:));
    mode_slope_of(k) = fit_coef(2);
end

%% One session, with its two marginals
[~, piqueck] = max(bv_range .* drawn);
pick_row = row_list(pick);
this_heat = heat{pick};
fprintf('\nrepresentative %s\n   %d frames | %.3f um/px | diameter range %.2f um\n', ...
    string(central.Directory(pick_row)), n_frame(pick), ...
    resolution_um(pick_row), bv_range(pick));

fig_one = figure('Color', 'w', 'Position', [80 60 980 780], 'Name', 'pvs_diameter_density');
layout_one = tiledlayout(fig_one, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

axis_main = nexttile(layout_one, 5, [3 3]);
x_um = this_heat.x_baseceneters;
y_um = this_heat.y_baseceneters;
counts = this_heat.xy_counts_clean;
pcolor(axis_main, x_um, y_um, log10(counts + 1e-4));
shading(axis_main, 'flat');
colormap(axis_main, param.colormap);
hold(axis_main, 'on')
mode_handle = plot(axis_main, x_um, this_heat.modepvs, '-', 'Color', [1 1 1], 'LineWidth', 2);
ref_handle = gobjects(1, numel(param.ref_slope));
for s = 1:numel(param.ref_slope)
    ref_handle(s) = plot(axis_main, x_um, param.ref_slope(s) * x_um, '--', ...
        'LineWidth', 1.2, 'Color', [1 1 1] * (0.30 + 0.16 * s));
end
xlabel(axis_main, 'lumen diameter, um from this session''s mode')
ylabel(axis_main, param.y_label)
% the fit is labelled where it is drawn, so a panel can be read without the console
mode_label = sprintf('mode per diameter bin | frame fit %+.3f, r %+.2f', ...
    frame_slope(pick), frame_r(pick));
ref_label = param.ref_name + label_reference(param.ref_slope, param.ref_r);
legend(axis_main, [mode_handle, ref_handle], [string(mode_label), ref_label], ...
    'Location', 'southwest', 'TextColor', 'w', 'Color', [0 0 0], 'Box', 'off', ...
    'FontSize', param.fontsize - 2)
set(axis_main, 'Color', [0 0 0], 'FontSize', param.fontsize, 'Layer', 'top', 'Box', 'on')
axis_main.Toolbar.Visible = 'off';

axis_top = nexttile(layout_one, 1, [1 3]);
bar(axis_top, x_um, sum(counts, 1, 'omitmissing'), 1, ...
    'FaceColor', [0.2 0.6 1], 'EdgeColor', 'none');
ylabel(axis_top, '% of frames')
title(axis_top, 'lumen diameter')
set(axis_top, 'XTickLabel', [], 'FontSize', param.fontsize - 1, 'Box', 'off')
xlim(axis_top, [min(x_um) max(x_um)])

axis_right = nexttile(layout_one, 8, [3 1]);
barh(axis_right, y_um, sum(counts, 2, 'omitmissing'), 1, ...
    'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none');
xlabel(axis_right, '% of frames')
title(axis_right, param.y_series)
set(axis_right, 'YTickLabel', [], 'FontSize', param.fontsize - 1, 'Box', 'off')
ylim(axis_right, [min(y_um) max(y_um)])

session_label = string(central.MouseID(pick_row)) + " " + string(central.Date(pick_row)) + ...
    " " + string(central.VesselID(pick_row));
title(layout_one, sprintf('%s   |   %d frames, %.3f um/px', ...
    session_label, n_frame(pick), resolution_um(pick_row)))
save_figure(fig_one, dirs.save_dir, 'pvs_diameter_density');

%% Many sessions, each on its own origin
[~, order] = sort(bv_range .* drawn, 'descend');
n_show = min(param.n_tile, sum(drawn));
tile_row = order(1:n_show);
n_grid_row = ceil(n_show / param.n_col);
fig_grid = figure('Color', 'w', 'Name', 'pvs_diameter_grid', ...
    'Position', [20 20 220 * param.n_col, 190 * n_grid_row + 90]);
layout_grid = tiledlayout(fig_grid, n_grid_row, param.n_col, ...
    'TileSpacing', 'tight', 'Padding', 'compact');
fprintf('grid : %d sessions in %d x %d\n', n_show, n_grid_row, param.n_col);
for k = 1:numel(tile_row)
    idx = tile_row(k);
    axis_tile = nexttile(layout_grid);
    tile_heat = heat{idx};
    pcolor(axis_tile, tile_heat.x_baseceneters, tile_heat.y_baseceneters, ...
        log10(tile_heat.xy_counts_clean + 1e-4));
    shading(axis_tile, 'flat');
    colormap(axis_tile, param.colormap);
    hold(axis_tile, 'on')
    plot(axis_tile, tile_heat.x_baseceneters, tile_heat.modepvs, '-', ...
        'Color', [1 1 1], 'LineWidth', 1.5);
    % the same two references the other panels rule, taken from param rather than
    % written again -- a second copy is a number that goes stale on its own
    tile_measured = param.ref_slope(ref_measured) * tile_heat.x_baseceneters;
    tile_noflux = param.ref_slope(ref_noflux) * tile_heat.x_baseceneters;
    plot(axis_tile, tile_heat.x_baseceneters, tile_measured, '--', ...
        'Color', [1 0.85 0.3], 'LineWidth', 1.2);
    plot(axis_tile, tile_heat.x_baseceneters, tile_noflux, ':', ...
        'Color', [0.6 0.8 1], 'LineWidth', 1.2);
    set(axis_tile, 'Color', [0 0 0], 'FontSize', 6, 'Layer', 'top', 'Box', 'on', ...
        'XTick', [], 'YTick', [])
    axis_tile.Toolbar.Visible = 'off';   % it exports into the panel otherwise
    tile_label = string(central.MouseID(row_list(idx))) + " " + ...
        string(central.Date(row_list(idx))) + " " + string(central.VesselID(row_list(idx)));
    % the slope and the scatter, because those are the two things the eye is judging
    title(axis_tile, sprintf('%s  %.2f r%.2f  sd%.2f', tile_label, frame_slope(idx), ...
        frame_r(idx), cond_sd(idx)), 'FontSize', 6)
    axis(axis_tile, 'tight')
end
grid_legend = sprintf('white = mode per bin | yellow dashed = %s %.3f | blue dotted = %s %.3f', ...
    param.ref_name(ref_measured), param.ref_slope(ref_measured), param.ref_name(ref_noflux), param.ref_slope(ref_noflux));
title(layout_grid, [char(param.y_series) ' against lumen diameter, one panel per session, each on its own mode' ...
    newline grid_legend ...
    newline 'beside each label: the frame-level fitted slope, its r, and the PVS scatter at a fixed diameter. Ticks are dropped; every panel spans its own range'])
save_figure(fig_grid, dirs.save_dir, 'pvs_diameter_grid');

%% What the mode curve says, gated and ungated
% The gate is reported against the population it removed, because whether the
% answer moves IS the evidence that the gate is not choosing it.
fprintf('\nmode-curve slope d(%s)/d(bv)\n', param.y_series);
report_slopes('every session that spans enough diameter', ...
    mode_slope_of(wide_enough), param.ref_slope, param.ref_name);
held_back = wide_enough & ~drawn;
if any(held_back)
    report_slopes(sprintf('scatter at a fixed diameter under %g um', param.max_cond_sd), ...
        mode_slope_of(drawn), param.ref_slope, param.ref_name);
    report_slopes('the ones the scatter gate held back', ...
        mode_slope_of(held_back), param.ref_slope, param.ref_name);
end

%% Pool at the vessel level
% Sessions cannot be averaged as they stand. Each carries a different number of
% frames and reaches a different distance along the diameter axis, so a plain mean
% would let the widest-ranging sessions set the wings and the narrowest set the
% centre, and the pooled band's thickness would be a picture of which session
% reached where. Normalising each DIAMETER COLUMN to sum to one first turns every
% column into the distribution of PVS thickness AT that diameter, and the average
% of those is a quantity with a meaning: where the PVS sits at this much dilation,
% averaged over vessels.
%   why  vessel and not session : one vessel here contributes five sessions and
%        most contribute one, so a session-level mean would weight that vessel five
%        times. This is the same bout -> session -> vessel -> mouse collapse the
%        rest of the leg uses
%   note  a column is drawn only where at least param.min_vessel_n vessels reach,
%        the rule polar_ave_fig applies to its angles
in_pool = wide_enough;
if ~isempty(param.pool_resolution)
    same_setting = abs(resolution_um(row_list) - param.pool_resolution) <= param.resolution_tol;
    in_pool = in_pool & same_setting;
    fprintf('\npooling only the sessions at %.3f um/px : %d of %d\n', ...
        param.pool_resolution, sum(in_pool), sum(wide_enough));
end
if ~any(in_pool)
    error('pvs_diameter_fig:noPool', ...
        'no session at %.3f um/px; the settings present are %s', ...
        param.pool_resolution, mat2str(unique(round(resolution_um(row_list(wide_enough)), 3))'));
end

pool_step = param.pool_step;
if isempty(pool_step)
    pool_step = max(resolution_um(row_list(in_pool)));
end
x_edge = 0;
y_edge = 0;
for k = 1:n_session
    if ~in_pool(k)
        continue
    end
    x_edge = max(x_edge, max(abs(heat{k}.x_baseceneters)));
    y_edge = max(y_edge, max(abs(heat{k}.y_baseceneters)));
end
grid_x = -x_edge : pool_step : x_edge;
grid_y = -y_edge : pool_step : y_edge;
[mesh_x, mesh_y] = meshgrid(grid_x, grid_y);
fprintf('\npooling on a %.3f um grid, %d x %d\n', pool_step, numel(grid_y), numel(grid_x));
%%
vessel_key = string(central.MouseID(row_list)) + "_" + string(central.VesselID(row_list));
[unique_vessel, ~, vessel_of] = unique(vessel_key(in_pool), 'stable');
n_vessel = numel(unique_vessel);
vessel_of_all = nan(n_session, 1);
for k = 1:n_session
    if ~in_pool(k)
        continue
    end
    vessel_of_all(k) = find(unique_vessel == vessel_key(k), 1);
end
vessel_sum = zeros(numel(grid_y), numel(grid_x), n_vessel);
vessel_count = zeros(n_vessel, numel(grid_x));

for k = 1:n_session
    if ~in_pool(k)
        continue
    end
    source_x = heat{k}.x_baseceneters;
    source_y = heat{k}.y_baseceneters;
    counts = heat{k}.xy_counts_clean;
    counts(~isfinite(counts)) = 0;
    column_total = sum(counts, 1);
    filled = column_total > 0;
    counts(:, filled) = counts(:, filled) ./ column_total(filled);

    on_grid = interp2(source_x, source_y, counts, mesh_x, mesh_y, 'linear', 0);
    reaches = interp1(source_x, double(filled), grid_x, 'nearest', 0) > 0.5;
    v = vessel_of_all(k);
    vessel_sum(:, :, v) = vessel_sum(:, :, v) + on_grid;
    vessel_count(v, reaches) = vessel_count(v, reaches) + 1;
end

% one map per vessel, each column divided by the sessions that reached it
vessel_map = nan(numel(grid_y), numel(grid_x), n_vessel);
for v = 1:n_vessel
    reached = vessel_count(v, :) > 0;
    if ~any(reached)
        continue
    end
    scaled = vessel_sum(:, reached, v) ./ vessel_count(v, reached);
    vessel_map(:, reached, v) = scaled;
end
%%
vessel_reaches = sum(vessel_count > 0, 1);
pooled = mean(vessel_map, 3, 'omitnan');
thin_column = vessel_reaches < param.min_vessel_n;
pooled(:, thin_column) = NaN;
fprintf('   %d vessels from %d mice | columns kept %d of %d\n', n_vessel, ...
    numel(unique(extractBefore(unique_vessel, "_"))), sum(~thin_column), numel(grid_x));
kept_x = grid_x(~thin_column);
fprintf('   drawn from %+.1f to %+.1f um, where at least %d vessels reach\n', ...
    min(kept_x), max(kept_x), param.min_vessel_n);

% one number per column of the pooled distribution, by whichever statistic the
% per-vessel fit below uses, so the two figures are reading the same thing
pooled_mode = nan(1, numel(grid_x));
for j = 1:numel(grid_x)
    if thin_column(j)
        continue
    end
    weight = pooled(:, j);
    weight(~isfinite(weight) | weight < 0) = 0;
    total = sum(weight);
    if total <= 0
        continue
    end
    weight = weight / total;
    switch param.column_stat
        case "mode"
            [~, at_row] = max(weight);
            pooled_mode(j) = grid_y(at_row);
        case "centroid"
            pooled_mode(j) = sum(grid_y(:) .* weight);
        otherwise
            pooled_mode(j) = half_mass(grid_y, weight);
    end
end

fig_pool = figure('Color', 'w', 'Position', [70 50 900 820], 'Name', 'pvs_diameter_pooled');
layout_pool = tiledlayout(fig_pool, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

axis_pool = nexttile(layout_pool, 1, [3 1]);
pcolor(axis_pool, grid_x, grid_y, pooled);
shading(axis_pool, 'flat');
colormap(axis_pool, param.colormap);
hold(axis_pool, 'on')
pool_handle = plot(axis_pool, grid_x, pooled_mode, '-', 'Color', [1 1 1], 'LineWidth', 2);

% The bend, read off the pooled map rather than off one session. Fitting each side
% of the modal diameter separately is the claim the retired chain tried to make and
% could not, because it split the BIN LIST in half instead of the diameter axis.
%   caution  the two sides do not reach equally far. Vessels sit near the low end
%        of their own range, so the dilated side runs about twice as far as the
%        constricted one, and a fit over each side's whole reach compares a short
%        span against a long one. The SYMMETRIC fit below is the one drawn, because
%        it is the only version where the two slopes are the same measurement made
%        twice; the full-reach numbers are printed beside it
has_mode = isfinite(pooled_mode);
narrow_col = has_mode & grid_x < 0;
wide_col = has_mode & grid_x > 0;
pool_narrow = polyfit(grid_x(narrow_col), pooled_mode(narrow_col), 1);
pool_wide = polyfit(grid_x(wide_col), pooled_mode(wide_col), 1);

span = min(abs(min(grid_x(has_mode))), max(grid_x(has_mode)));
narrow_sym = narrow_col & grid_x >= -span;
wide_sym = wide_col & grid_x <= span;
sym_narrow = polyfit(grid_x(narrow_sym), pooled_mode(narrow_sym), 1);
sym_wide = polyfit(grid_x(wide_sym), pooled_mode(wide_sym), 1);

fit_handle = plot(axis_pool, grid_x(narrow_sym), polyval(sym_narrow, grid_x(narrow_sym)), ...
    '-', 'Color', [0.3 1 0.5], 'LineWidth', 2.5);
plot(axis_pool, grid_x(wide_sym), polyval(sym_wide, grid_x(wide_sym)), ...
    '-', 'Color', [0.3 1 0.5], 'LineWidth', 2.5);
%%
fprintf('   over each side''s whole reach\n');
fprintf('      narrower %.3f (%d columns) | wider %.3f (%d columns) | bend %+.3f\n', ...
    pool_narrow(1), sum(narrow_col), pool_wide(1), sum(wide_col), ...
    pool_wide(1) - pool_narrow(1));
fprintf('   over the symmetric span, +/- %.2f um\n', span);
fprintf('      narrower %.3f (%d columns) | wider %.3f (%d columns) | bend %+.3f\n', ...
    sym_narrow(1), sum(narrow_sym), sym_wide(1), sum(wide_sym), ...
    sym_wide(1) - sym_narrow(1));

pool_ref = gobjects(1, numel(param.ref_slope));
tmp.lineshape = ["--","-"];
for s = 1:numel(param.ref_slope)
    pool_ref(s) = plot(axis_pool, kept_x, param.ref_slope(s) * kept_x, tmp.lineshape(s), ...
        'LineWidth', 1.2, 'Color', [1 1 1] * (0.30 + 0.16 * s));
end
%%
ylabel(axis_pool, param.y_label)
% r of each side fit, against the column statistic the line was fitted to
r_narrow = corr(grid_x(narrow_sym)', pooled_mode(narrow_sym)');
r_wide = corr(grid_x(wide_sym)', pooled_mode(wide_sym)');
fit_label = sprintf('either side of the mode: %+.3f (r %+.2f) | %+.3f (r %+.2f)', ...
    sym_narrow(1), r_narrow, sym_wide(1), r_wide);
legend(axis_pool, [pool_handle, fit_handle, pool_ref], ...
    [param.column_stat + " of the pooled column", string(fit_label), ...
    param.ref_name + label_reference(param.ref_slope, param.ref_r)], ...
    'Location', 'southwest', 'TextColor', 'w', 'Color', [0 0 0], 'Box', 'off', ...
    'FontSize', param.fontsize - 2)
set(axis_pool, 'Color', [0 0 0], 'FontSize', param.fontsize, 'Layer', 'top', ...
    'Box', 'on', 'XTickLabel', [])
axis_pool.Toolbar.Visible = 'off';
xlim(axis_pool, [min(kept_x) max(kept_x)])
xlabel(axis_pool, 'lumen diameter, um from each session''s own mode')

%%
axis_count = nexttile(layout_pool, 4, [1 1]);
bar(axis_count, grid_x, vessel_reaches, 1, 'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none');
hold(axis_count, 'on')
yline(axis_count, param.min_vessel_n, '-', 'drawn above this', 'LineWidth', 1);
xlabel(axis_count, 'lumen diameter, um from each session''s own mode')
ylabel(axis_count, 'vessels')
set(axis_count, 'FontSize', param.fontsize - 1, 'Box', 'off')
xlim(axis_count, [min(kept_x) max(kept_x)])

title(layout_pool, [char(param.y_series + " given lumen diameter, pooled over vessels") newline ...
    'each diameter column is a distribution that sums to one, averaged within vessel then across'])
save_figure(fig_pool, dirs.save_dir, 'pvs_diameter_pooled');

%% The bend, once per vessel, so it has a distribution and not just a value
% The pooled map gives one number and no spread. Fitting each vessel's own map on
% its own symmetric span makes the two slopes a paired measurement within that
% vessel, so the difference has a null of exactly zero and the vessels can be
% counted.
%   why  each vessel gets ITS OWN span rather than a common one : a common span
%        would be set by the vessel that reached least far and would throw away
%        most of the others. Within a vessel the two sides are still the same
%        length, which is the only symmetry the comparison needs
%   note  the two axes below are the CONSTRICTED and DILATED sides, fixed by
%        geometry. The retired chain plotted the smaller slope against the larger
%        after sorting them, which puts every point on one side of the diagonal
%        whatever the data says. see STRUCTURE.md
n_stat = numel(param.stat_name);
vessel_narrow_all = nan(n_vessel, n_stat);
vessel_wide_all = nan(n_vessel, n_stat);
vessel_span = nan(n_vessel, 1);
for v = 1:n_vessel
    map_v = vessel_map(:, :, v);
    covered = find(vessel_count(v, :) > 0);
    stat_v = nan(n_stat, numel(grid_x));
    for j = covered
        weight = map_v(:, j);
        weight(~isfinite(weight) | weight < 0) = 0;
        total = sum(weight);
        if total <= 0
            continue
        end
        weight = weight / total;
        [~, at_row] = max(weight);
        stat_v(1, j) = grid_y(at_row);
        stat_v(2, j) = sum(grid_y(:) .* weight);
        stat_v(3, j) = half_mass(grid_y, weight);
    end
    has_v = all(isfinite(stat_v), 1);
    if ~any(has_v & grid_x < 0) || ~any(has_v & grid_x > 0)
        continue
    end
    span_v = min(abs(min(grid_x(has_v))), max(grid_x(has_v)));
    on_narrow = has_v & grid_x < 0 & grid_x >= -span_v;
    on_wide = has_v & grid_x > 0 & grid_x <= span_v;
    if sum(on_narrow) < param.min_column || sum(on_wide) < param.min_column
        continue
    end
    for s = 1:n_stat
        narrow_fit_v = polyfit(grid_x(on_narrow), stat_v(s, on_narrow), 1);
        wide_fit_v = polyfit(grid_x(on_wide), stat_v(s, on_wide), 1);
        vessel_narrow_all(v, s) = narrow_fit_v(1);
        vessel_wide_all(v, s) = wide_fit_v(1);
    end
    vessel_span(v) = span_v;
end

fprintf('\nbend per vessel, the same fit off three different column statistics\n');
fprintf('%-10s %8s %10s %12s %10s %10s\n', 'statistic', 'vessels', 'constricted', ...
    'dilated', 'bend', 'p');
for s = 1:n_stat
    usable = isfinite(vessel_narrow_all(:, s)) & isfinite(vessel_wide_all(:, s));
    difference = vessel_wide_all(usable, s) - vessel_narrow_all(usable, s);
    fprintf('%-10s %8d %10.3f %12.3f %+10.3f %10.4g\n', param.stat_name(s), ...
        sum(usable), median(vessel_narrow_all(usable, s)), ...
        median(vessel_wide_all(usable, s)), median(difference), signrank(difference));
end

which_stat = find(param.stat_name == param.column_stat, 1);
vessel_narrow = vessel_narrow_all(:, which_stat);
vessel_wide = vessel_wide_all(:, which_stat);

fitted = isfinite(vessel_narrow) & isfinite(vessel_wide);
vessel_bend = vessel_wide(fitted) - vessel_narrow(fitted);
fprintf('\nbend per vessel, each on its own symmetric span, %d of %d vessels\n', ...
    sum(fitted), n_vessel);
fprintf('   span used : median %.2f um | range %.2f to %.2f\n', ...
    median(vessel_span(fitted)), min(vessel_span(fitted)), max(vessel_span(fitted)));
fprintf('   constricted side : median %.3f | IQR %.3f to %.3f\n', ...
    median(vessel_narrow(fitted)), prctile(vessel_narrow(fitted), 25), ...
    prctile(vessel_narrow(fitted), 75));
fprintf('   dilated side     : median %.3f | IQR %.3f to %.3f\n', ...
    median(vessel_wide(fitted)), prctile(vessel_wide(fitted), 25), ...
    prctile(vessel_wide(fitted), 75));
fprintf('   difference, null 0 : median %+.3f | IQR %+.3f to %+.3f | %d of %d positive\n', ...
    median(vessel_bend), prctile(vessel_bend, 25), prctile(vessel_bend, 75), ...
    sum(vessel_bend > 0), numel(vessel_bend));
fprintf('   signed rank p %.4g\n', signrank(vessel_bend));
fprintf('   against the no-flux slope %.3f : constricted below it %d of %d, dilated %d of %d\n', ...
    param.ref_slope(ref_noflux), sum(vessel_narrow(fitted) < param.ref_slope(ref_noflux)), sum(fitted), ...
    sum(vessel_wide(fitted) < param.ref_slope(ref_noflux)), sum(fitted));

fig_bend = figure('Color', 'w', 'Position', [90 70 940 440], 'Name', 'pvs_diameter_bend');
layout_bend = tiledlayout(fig_bend, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axis_pair = nexttile(layout_bend);
hold(axis_pair, 'on')
axis_limit = [min([vessel_narrow(fitted); vessel_wide(fitted)]) - 0.1, ...
    max([vessel_narrow(fitted); vessel_wide(fitted)]) + 0.1];
plot(axis_pair, axis_limit, axis_limit, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
xline(axis_pair, param.ref_slope(ref_noflux), '--', 'no flux', 'Color', [0.4 0.4 0.4]);
yline(axis_pair, param.ref_slope(ref_noflux), '--', 'no flux', 'Color', [0.4 0.4 0.4]);
scatter(axis_pair, vessel_narrow(fitted), vessel_wide(fitted), 42, ...
    'MarkerFaceColor', clee.clist.nrem, 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.75);
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

%% Does the curve bend? Fit the two sides of the modal diameter separately
% Each session is centred on its own modal diameter, so x = 0 is a real anchor:
% below it the vessel is narrower than it usually is, above it wider. The retired
% chain split at floor(n/2) instead, which is the middle of the BIN LIST and moves
% with the bin count rather than sitting anywhere on the diameter axis.
% see STRUCTURE.md
slope_narrow = nan(n_session, 1);
slope_wide = nan(n_session, 1);
for k = 1:n_session
    if ~wide_enough(k)
        continue
    end
    x_um = heat{k}.x_baseceneters(:);
    y_um = heat{k}.modepvs(:);
    on_narrow = x_um < 0;
    on_wide = x_um > 0;
    if sum(on_narrow) < 5 || sum(on_wide) < 5
        continue
    end
    narrow_fit = robustfit(x_um(on_narrow), y_um(on_narrow));
    wide_fit = robustfit(x_um(on_wide), y_um(on_wide));
    slope_narrow(k) = narrow_fit(2);
    slope_wide(k) = wide_fit(2);
end
both_sides = isfinite(slope_narrow) & isfinite(slope_wide);
bend = slope_wide(both_sides) - slope_narrow(both_sides);
fprintf('\nbend across the modal diameter, %d sessions\n', sum(both_sides));
fprintf('   narrower than the mode : median %.3f\n', median(slope_narrow(both_sides)));
fprintf('   wider than the mode    : median %.3f\n', median(slope_wide(both_sides)));
fprintf('   difference, null 0     : median %+.3f | IQR %+.3f to %+.3f | %d of %d positive\n', ...
    median(bend), prctile(bend, 25), prctile(bend, 75), sum(bend > 0), numel(bend));
fprintf('   signed rank p %.3g\n', signrank(bend));
fprintf('   for scale, the PVS scatter at a fixed diameter is %.2f um at the median\n', ...
    median(cond_sd(wide_enough), 'omitnan'));

fprintf('\nPVS scatter at a fixed diameter, um : ');
fprintf('%.2f ', prctile(cond_sd(wide_enough), [5 25 50 75 95]));
fprintf('  (5 / 25 / 50 / 75 / 95)\n');

%% ---------------------------------------------------------------- helpers
function value = half_mass(centre, weight)
    % Where the column's cumulative weight crosses one half, read between bins so
    % the answer is not pinned to the grid the way the brightest bin is.
    cumulative = cumsum(weight(:));
    above = find(cumulative >= 0.5, 1);
    if isempty(above)
        value = NaN;
        return
    end
    if above == 1
        value = centre(1);
        return
    end
    below_mass = cumulative(above - 1);
    step_mass = cumulative(above) - below_mass;
    if step_mass <= 0
        value = centre(above);
        return
    end
    fraction = (0.5 - below_mass) / step_mass;
    value = centre(above - 1) + fraction * (centre(above) - centre(above - 1));
end

function report_slopes(label, slope_list, ref_slope, ref_name)
    slope_list = slope_list(isfinite(slope_list));
    if isempty(slope_list)
        return
    end
    fprintf('\n   %s : %d sessions\n', label, numel(slope_list));
    fprintf('      median %.3f | IQR %.3f to %.3f | range %.3f to %.3f\n', ...
        median(slope_list), prctile(slope_list, 25), prctile(slope_list, 75), ...
        min(slope_list), max(slope_list));
    for s = 1:numel(ref_slope)
        fprintf('      below %-24s (%6.3f) : %d of %d\n', ref_name(s), ...
            ref_slope(s), sum(slope_list < ref_slope(s)), numel(slope_list));
    end
end

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

function value = lead_number(text_column)
    % Depth reads '70um' and Resolution '0.19um', both transcribed off the
    % acquisition info.txt, so the number is at the front and the unit follows
    text_column = string(text_column);
    value = nan(numel(text_column), 1);
    for k = 1:numel(text_column)
        head = regexp(text_column(k), '^[\d.]+', 'match', 'once');
        if strlength(head) > 0
            value(k) = str2double(head);
        end
    end
end
