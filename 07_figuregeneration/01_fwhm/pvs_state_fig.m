%PVS_STATE_FIG  The PVS-diameter coupling split by arousal state, and read as area.
%   pvs_diameter_fig draws the coupling pooled over a whole recording. This asks the
%   two questions that need the recording taken apart: whether the arousal state
%   changes it, and what it looks like when the same relation is written as the
%   annulus AREA rather than as a slope.
%
%   WHY THE AREA AXIS EXISTS. On the eps axis the no-flux null is a ray from the
%   origin, so a curve shallower than the null sits BELOW it for x > 0 and ABOVE it
%   for x < 0 -- one panel read two ways. In area the null is a horizontal line and
%   both halves read the same direction. see CLAUDE_LOG.md
%
%   Reads centralized_paxfwhm.mat and centralized_sleep_score.mat, which is 2.6 MB.
%   It does NOT open centralized_paxfwhm_state.mat: that file's state_idx has
%   overlapping fields, consecutive bouts sharing a boundary frame, and a tail left
%   unlabelled by trim_to_duration, and the scoring it was built from is right here.
%
%   note  two protections an earlier design made mandatory were measured and
%        dropped. The transition band costs a third of the vessels and moves
%        nothing, and label contamination can only shrink a state difference, never
%        create one. The histogram matching was guarding against the bend alone, and
%        the bend does not track the operating-point gap. Both are still reachable
%        through param, and both should be switched back ON before any positive
%        result is believed. see FINDINGS.md
clc, clear
addpath('g:\03_program\01_ecspress\functions');   % where util_ecspath lives
util_ecspath;                                     % three roots, minus zz_notinuse
clee = color_lee();
dirs = util_centraldirs();

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);

param.depth_thr = 70;        % um     L1 / L2-3 boundary, the stratum every result uses
param.pool_resolution = 0.285;   % um/px  pool one pixel size only, nothing resampled
param.resolution_tol = 0.01;     % um/px  how close counts as the same setting
param.min_frame = 900;       % frames a cell needs before it is fitted
param.min_col = 20;          % frames a one-pixel column needs, per cell
param.min_side = 4;          % columns each side of the mode needs, for the bend
param.min_vessel = 6;        % vessels a pooled column needs before it is drawn
param.column_stat = "median";    % see column_stat; median is read between lattice points
param.fontsize = 10;

% WHICH THICKNESS THE FIGURE IS ABOUT.
%   "eps"   the outer PVS edges. Shares no boundary row with bv, so the errors are
%           independent and the pure-noise limit of the slope is 0
%   "area"  (pi/4)(eps^2 - bv^2). The same content with the null flattened, and the
%           only form in which "fluid left the plane" is a quantity
%   caution area is nonlinear in both, so it is computed per frame in ABSOLUTE
%        micrometres and referred to the vessel's resting value afterwards. It
%        cannot be computed from mode-centred values
param.y_space = "eps";

% HOW THE STATES ARE CUT.
%   "nrem_awake"       the fine codes, NREM against NotSleep. Drowsy and uArousal
%                      are in neither
%   "roughnrem_awake"  NREM + uArousal against NotSleep, which is what
%                      state_integration calls roughnrem. uArousal is a brief
%                      arousal INSIDE an NREM bout, so folding it in adds
%                      constrictions to the NREM cell and narrows its operating point
param.state_partition = "roughnrem_awake";

% Frames within this many of a label change are dropped from both cells. The 5 s
% scoring grid locates a transition only to +/- 2.5 s, which is 8 frames at 3 Hz.
%   0 keeps everything, which is what the sweep argued for while the answers are
%   null. see FINDINGS.md
param.edge_frame = 0;

% Subsample the two cells to a common bv column histogram before comparing. OFF:
% it removed a third of the vessels and the confound it guards was measured absent.
param.match_columns = false;

fprintf('y space %s | partition %s | edge band %d frames | matching %d\n\n', ...
    param.y_space, param.state_partition, param.edge_frame, param.match_columns);

central = load(fullfile(dirs.central, 'centralized_paxfwhm.mat')).central;
score_table = load(fullfile(dirs.central, 'centralized_sleep_score.mat')).central;
key_central = util_sessionkey(central);
key_score = util_sessionkey(score_table);

%% The stratum
depth_um = lead_number(central.Depth);
resolution_um = lead_number(central.Resolution);
is_artery = contains(string(central.VesselID), "PA", 'IgnoreCase', true);
same_setting = abs(resolution_um - param.pool_resolution) <= param.resolution_tol;
in_stratum = is_artery & depth_um < param.depth_thr & same_setting;
is_scored = ismember(key_central, key_score);
row_list = find(in_stratum & is_scored);
fprintf('PA, depth < %d um, %.3f um/px : %d sessions | with a scoring : %d\n', ...
    param.depth_thr, param.pool_resolution, sum(in_stratum), numel(row_list));

%% One row per session, two cells each
cell_name = ["sleep", "awake"];
n_session = numel(row_list);
grid_x = -4 : param.pool_resolution : 4;
n_grid = numel(grid_x);

session = struct( ...
    'vessel',      strings(n_session, 1), ...   % Nx1 str    mouse_vessel
    'slope',       nan(n_session, 2), ...       % Nx2 float  d(y)/d(bv), per cell
    'noflux',      nan(n_session, 2), ...       % Nx2 float  bv/eps, per cell
    'operating',   nan(n_session, 2), ...       % Nx2 float um  median bv, per cell
    'n_frame',     nan(n_session, 2), ...       % Nx2 int    frames per cell
    'offset',      nan(n_session, 1), ...       % Nx1 float  cell 1 minus cell 2 at one bv
    'bend',        nan(n_session, 1), ...       % Nx1 float  wide side minus narrow side
    'resting',     nan(n_session, 1));          % Nx1 float um2  annulus area at the mode
curve = nan(2, n_grid, n_session);              % 2 x G x N  the column statistic per cell
pooled = zeros(2, n_grid);
reached = zeros(2, n_grid);

for a = 1:n_session
    k = row_list(a);
    pax_fwhm = central.data{k};
    um_per_px = resolution_um(k);
    session.vessel(a) = string(central.MouseID(k)) + "_" + string(central.VesselID(k));

    bv_px = double(pax_fwhm.thickness.bv(:));
    bv_um = bv_px * um_per_px;
    eps_um = double(pax_fwhm.thickness.eps(:)) * um_per_px;
    finite_pair = isfinite(bv_um) & isfinite(eps_um);

    % the y series. area is absolute here on purpose -- see the param block
    if param.y_space == "area"
        y_value = (pi / 4) * (eps_um.^2 - bv_um.^2);
    else
        y_value = eps_um;
    end

    score_row = score_table.data{find(key_score == key_central(k), 1)};
    [state_code, state_info] = state_perframe(pax_fwhm.t_axis, ...
        double(score_row.behavState), score_row.binwidth_sec, param.edge_frame);
    usable = finite_pair & ~state_info.near_change;

    in_sleep = usable & state_code == score_row.statecodes.NREM;
    if param.state_partition == "roughnrem_awake" && isfield(score_row.statecodes, 'uarousal')
        in_sleep = in_sleep | (usable & state_code == score_row.statecodes.uarousal);
    end
    in_awake = usable & state_code == score_row.statecodes.NotSleep;
    cell_mask = {in_sleep, in_awake};
    session.n_frame(a, :) = [sum(in_sleep), sum(in_awake)];
    if any(session.n_frame(a, :) < param.min_frame)
        continue
    end

    if param.match_columns
        [keep_sleep, keep_awake] = match_columncounts(bv_px, in_sleep, in_awake, param.min_col);
        cell_mask = {keep_sleep, keep_awake};
    end

    % ONE origin over both cells, so the grid means the same thing for each
    origin_px = mode(bv_px(usable));

    for c = 1:2
        picked = cell_mask{c};
        fit_coef = robustfit(bv_um(picked), y_value(picked));
        session.slope(a, c) = fit_coef(2);
        session.noflux(a, c) = median(bv_um(picked)) / median(eps_um(picked));
        session.operating(a, c) = median(bv_um(picked));

        [stat_here, stat_info] = column_stat(bv_px(picked), y_value(picked), ...
            param.column_stat, param.min_col);
        centred_um = (stat_info.x_px - origin_px) * um_per_px;
        curve(c, :, a) = interp1(centred_um, stat_here, grid_x, 'linear', NaN);
    end

    % the resting value, and the two derived quantities that need both cells
    [all_stat, all_info] = column_stat(bv_px(usable), y_value(usable), ...
        param.column_stat, param.min_col);
    all_centred = (all_info.x_px - origin_px) * um_per_px;
    session.resting(a) = interp1(all_centred, all_stat, 0, 'linear', NaN);
    [~, ~, bend_info] = fit_bothsides(all_centred, all_stat, param.min_side);
    session.bend(a) = bend_info.bend;

    shared = isfinite(curve(1, :, a)) & isfinite(curve(2, :, a));
    if any(shared)
        session.offset(a) = median(curve(1, shared, a) - curve(2, shared, a));
    end

    % pool each cell's curve, referred to the resting value where that makes sense
    for c = 1:2
        this_curve = curve(c, :, a);
        if param.y_space == "area" && isfinite(session.resting(a)) && session.resting(a) > 0
            this_curve = 100 * (this_curve - session.resting(a)) / session.resting(a);
            curve(c, :, a) = this_curve;
        end
        filled = isfinite(this_curve);
        pooled(c, filled) = pooled(c, filled) + this_curve(filled);
        reached(c, filled) = reached(c, filled) + 1;
    end
end
clear central

%% Collapse to the vessel, then report
fitted = isfinite(session.slope(:, 1)) & isfinite(session.slope(:, 2));
[vessel_name, ~, vessel_of] = unique(session.vessel(fitted), 'stable');
n_vessel = numel(vessel_name);
mouse_name = extractBefore(vessel_name, "_");
vessel = struct('slope', nan(n_vessel, 2), 'noflux', nan(n_vessel, 2), ...
    'operating', nan(n_vessel, 2), 'offset', nan(n_vessel, 1), 'bend', nan(n_vessel, 1));
sub = structfun(@(v) v(fitted, :), session, 'UniformOutput', false);
for v = 1:n_vessel
    picked = vessel_of == v;
    for f = ["slope", "noflux", "operating", "offset", "bend"]
        vessel.(f)(v, :) = median(sub.(f)(picked, :), 1, 'omitnan');
    end
end
fprintf('fitted %d sessions | %d vessels | %d mice\n', sum(fitted), n_vessel, ...
    numel(unique(mouse_name)));
fprintf('frames per session : %s %.0f | %s %.0f\n\n', ...
    cell_name(1), median(sub.n_frame(:, 1)), cell_name(2), median(sub.n_frame(:, 2)));

fprintf('%-26s %10s %10s %11s\n', '', cell_name(1), cell_name(2), 'difference');
report_pair('slope d(y)/d(bv)', vessel.slope(:, 1), vessel.slope(:, 2));
report_pair('no-flux null bv/eps', vessel.noflux(:, 1), vessel.noflux(:, 2));
report_pair('operating point, um', vessel.operating(:, 1), vessel.operating(:, 2));

fprintf('\npaired across vessels, null 0\n');
report_signed('slope difference', vessel.slope(:, 1) - vessel.slope(:, 2), mouse_name);
report_signed('offset at one bv', vessel.offset, mouse_name);
report_signed('bend, both cells', vessel.bend, mouse_name);

% is a slope difference just the bend, read at two operating points
gap = vessel.operating(:, 1) - vessel.operating(:, 2);
difference = vessel.slope(:, 1) - vessel.slope(:, 2);
usable_pair = isfinite(gap) & isfinite(difference);
[rho, p_rho] = corr(gap(usable_pair), difference(usable_pair), 'type', 'Spearman');
fprintf('\ncorr(operating-point gap, slope difference) : rho %+0.3f, p %.3g, n %d\n', ...
    rho, p_rho, sum(usable_pair));
fprintf('a slope difference that is really the bend would track this\n');
fprintf('detectable slope difference at n = %d : %.3f\n', n_vessel, ...
    std(difference(usable_pair)) * 2.80 / sqrt(n_vessel));

%% Figure
mean_curve = pooled ./ max(reached, 1);
mean_curve(reached < param.min_vessel) = NaN;
if param.y_space == "area"
    y_text = 'change in annulus area, % of resting';
    null_text = 'no flux';
else
    y_text = 'eps, um from the mode';
    null_text = 'no flux, bv/eps';
end

fig_state = figure('Color', 'w', 'Name', 'pvs_state', 'Position', [40 40 1180 430]);
layout_state = tiledlayout(fig_state, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cell_color = [clee.clist.nrem; 0.45 0.45 0.45];

axis_all = nexttile(layout_state, 1);
hold(axis_all, 'on')
for a = find(fitted)'
    plot(axis_all, grid_x, curve(1, :, a), '-', 'Color', [0.78 0.78 0.80], 'LineWidth', 0.7);
end
draw_null(axis_all, param, vessel, grid_x, null_text);
plot(axis_all, grid_x, mean_curve(1, :), '-', 'Color', cell_color(1, :), 'LineWidth', 3);
finish_axis(axis_all, param, 'lumen diameter, um from the mode', y_text, ...
    sprintf('%s, one grey line per session', cell_name(1)));

axis_both = nexttile(layout_state, 2);
hold(axis_both, 'on')
draw_null(axis_both, param, vessel, grid_x, null_text);
for c = 1:2
    plot(axis_both, grid_x, mean_curve(c, :), '-', 'Color', cell_color(c, :), 'LineWidth', 2.6);
end
legend(axis_both, [null_text, cell_name], 'Location', 'southwest', 'Box', 'off', ...
    'FontSize', param.fontsize - 1)
finish_axis(axis_both, param, 'lumen diameter, um from the mode', y_text, ...
    sprintf('both cells, %d vessels', n_vessel));

axis_pair = nexttile(layout_state, 3);
hold(axis_pair, 'on')
axis_limit = [min(vessel.slope, [], 'all'), max(vessel.slope, [], 'all')] + [-0.05 0.05];
plot(axis_pair, axis_limit, axis_limit, '-', 'Color', [0.6 0.6 0.6]);
scatter(axis_pair, vessel.slope(:, 2), vessel.slope(:, 1), 46, ...
    'MarkerFaceColor', cell_color(1, :), 'MarkerEdgeColor', 'none');
axis(axis_pair, 'equal')
xlim(axis_pair, axis_limit)
ylim(axis_pair, axis_limit)
finish_axis(axis_pair, param, sprintf('slope in %s', cell_name(2)), ...
    sprintf('slope in %s', cell_name(1)), 'one point per vessel');

title(layout_state, sprintf(['%s against lumen diameter, split by arousal state' newline ...
    'partition %s | transition band %d frames | column matching %d'], ...
    param.y_space, param.state_partition, param.edge_frame, param.match_columns), ...
    'Interpreter', 'none')
save_figure(fig_state, dirs.save_dir, "pvs_state_" + param.y_space);

%% ---------------------------------------------------------------- helpers
function [keep_a, keep_b] = match_columncounts(x_px, mask_a, mask_b, min_col)
    % A GATE: it returns verdicts and the caller applies them. Each one-pixel column
    % keeps the same number of frames from both cells, so the two share an x
    % distribution, a mean and a variance -- and therefore the same attenuation.
    rng(11);                       % fixed, so the same subsample comes back tomorrow
    keep_a = false(size(x_px));
    keep_b = false(size(x_px));
    for column_value = unique(x_px(mask_a | mask_b))'
        in_a = find(mask_a & x_px == column_value);
        in_b = find(mask_b & x_px == column_value);
        take = min(numel(in_a), numel(in_b));
        if take < min_col
            continue
        end
        keep_a(in_a(randperm(numel(in_a), take))) = true;
        keep_b(in_b(randperm(numel(in_b), take))) = true;
    end
end

function draw_null(target_axis, param, vessel, grid_x, null_text)
    % in area the null is a horizontal line; on the eps axis it is each cell's own
    % ray, and the two cells do not share one -- their geometry differs
    if param.y_space == "area"
        yline(target_axis, 0, '--', null_text, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);
        return
    end
    for c = 1:2
        null_slope = median(vessel.noflux(:, c), 'omitnan');
        plot(target_axis, grid_x, null_slope * grid_x, '--', ...
            'Color', [0.55 0.55 0.55] + 0.15 * (c - 1), 'LineWidth', 1.2);
    end
end

function finish_axis(target_axis, param, x_text, y_text, title_text)
    grid(target_axis, 'on')
    set(target_axis, 'FontSize', param.fontsize, 'Box', 'on')
    target_axis.Toolbar.Visible = 'off';   % it exports into the panel otherwise
    xlabel(target_axis, x_text)
    ylabel(target_axis, y_text)
    title(target_axis, title_text, 'FontSize', param.fontsize)
end

function report_pair(label, value_a, value_b)
    median_a = median(value_a, 'omitnan');
    median_b = median(value_b, 'omitnan');
    fprintf('%-26s %+10.4f %+10.4f %+11.4f\n', label, median_a, median_b, median_a - median_b);
end

function report_signed(label, difference, mouse_name)
    keep = isfinite(difference);
    difference = difference(keep);
    if isempty(difference)
        fprintf('%-22s (none)\n', label);
        return
    end
    [unique_mouse, ~, mouse_of] = unique(mouse_name(keep));
    per_mouse = accumarray(mouse_of, difference, [], @median);
    fprintf('%-22s median %+8.4f | IQR %+7.3f .. %+7.3f | %2d/%2d pos | p %.4g | mouse p %.4g (n %d)\n', ...
        label, median(difference), prctile(difference, 25), prctile(difference, 75), ...
        sum(difference > 0), numel(difference), signrank(difference), ...
        signrank(per_mouse), numel(unique_mouse));
end

function save_figure(fig, save_dir, fig_name)
    % R2023b exportgraphics cannot write svg, so the vector copy goes through print,
    % which is what make_fig.save2svg does. see CLAUDE.md
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    print(fig, fullfile(save_dir, fig_name + ".svg"), '-dsvg', '-vector');
    exportgraphics(fig, fullfile(save_dir, fig_name + ".png"), 'Resolution', 200);
    fprintf('wrote %s.svg and .png to %s\n', fig_name, save_dir);
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
