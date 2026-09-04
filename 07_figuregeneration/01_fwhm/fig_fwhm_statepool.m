%FIG_FWHM_STATEPOOL  The pooled PVS-diameter map, one panel per arousal state.
%   Reads fwhmrelation.mat and puts state_pool on axes; nothing is measured below.
%   state_pool is normalised over the whole map (brightness is occupancy), unlike
%   vessel_pool's per-column rule in fig_fwhm_diameter. The three panels share one
%   colour axis and both limits, and the zero is the SESSION median taken before the
%   state mask, so the shift between panels is a measurement.
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

param.colormap = clee.gradient.inferno;   % a heatmap, so a perceptually ordered ramp
param.fontsize = 11;
param.state_label = ["Awake", "NREM", "REM"];   % what each panel is called (str)
param.color_decade = 2.5;   % decades of weight the ramp spans, below the pooled max

%% What is drawn
relation_path = fullfile(dirs.save_dir, 'fwhmrelation.mat');
relation = load(relation_path).save_content;
state_pool = relation.state_pool;
grid_x = relation.grid_x;
grid_y = relation.vessel_pool.grid_y;
y_label = relation.param.y_series + ", um from its mode";

fprintf('fwhmrelation %s\n   state_pool : %s\n', relation_path, ...
    strjoin(state_pool.state', ', '));
disp(state_pool(:, ["state", "n_session", "n_vessel", "n_mouse", ...
    "bv_statemedian", "y_statemedian"]));

%% One colour axis over the three panels
% the maps go in as multiples of a uniform map, one constant for every panel, so they
% clear plot_heatpanel's log floor without any ratio moving
uniform_weight = 1 / (numel(grid_x) * numel(grid_y));
n_state = height(state_pool);
drawn_map = cell(n_state, 1);
for k = 1:n_state
    drawn_map{k} = state_pool.map{k} / uniform_weight;
end
% the ramp runs a fixed number of decades below the pooled maximum; the smallest
% positive cell is an interpolation tail
pooled_value = vertcat(drawn_map{:});
pooled_value = pooled_value(isfinite(pooled_value) & pooled_value > 0);
color_top = log10(max(pooled_value) + 1e-4);
color_limit = [color_top - param.color_decade, color_top];
x_limit = [min(grid_x), max(grid_x)];
y_limit = [min(grid_y), max(grid_y)];
fprintf('colour axis log10 %.2f to %.2f, in multiples of a uniform map\n', color_limit);

%% The board
fig_state = figure('Color', 'w', 'Position', [60 60 1220 700], ...
    'Name', 'fig_fwhm_statepool');
layout_state = tiledlayout(fig_state, 2, n_state, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

axis_map = gobjects(1, n_state);
for k = 1:n_state
    axis_map(k) = nexttile(layout_state, k);
    % the null each panel is read against, anchored at this state's operating point
    ref_noflux = reference_curve(grid_x, state_pool.bv_statemedian(k), ...
        state_pool.y_statemedian(k), 1, relation.param.y_series);
    handle_map = plot_heatpanel(axis_map(k), grid_x, grid_y, drawn_map{k}, ...
        'curve', state_pool.mode_curve{k}, 'ref_curve', ref_noflux, ...
        'ref_style', "--", 'colormap', param.colormap, 'fontsize', param.fontsize);
    clim(axis_map(k), color_limit);
    xlim(axis_map(k), x_limit);
    ylim(axis_map(k), y_limit);
    title(axis_map(k), sprintf('%s   %d vessels, %d mice', param.state_label(k), ...
        state_pool.n_vessel(k), state_pool.n_mouse(k)));
    if k == 1
        ylabel(axis_map(k), y_label);
        legend(axis_map(k), [handle_map.curve, handle_map.ref], ...
            [relation.param.column_stat + " of the column", "no flux, area conserved"], ...
            'Location', 'northwest', 'Box', 'off', 'TextColor', 'w');
    end
end
color_bar = colorbar(axis_map(n_state));
color_bar.Label.String = 'log_{10} weight, uniform map = 0';

axis_count = gobjects(1, n_state);
for k = 1:n_state
    axis_count(k) = nexttile(layout_state, n_state + k);
    bar(axis_count(k), grid_x, state_pool.vessel_reaches{k}, 1, ...
        'FaceColor', clee.clist.gray, 'EdgeColor', 'none');
    yline(axis_count(k), relation.filt.min_vessel, '-', 'drawn above this', ...
        'LineWidth', 1, 'FontSize', param.fontsize - 3);
    xlim(axis_count(k), x_limit);
    ylim(axis_count(k), [0, max(state_pool.n_vessel) + 2]);
    set(axis_count(k), 'FontSize', param.fontsize, 'Box', 'off');
    xlabel(axis_count(k), 'lumen diameter, um from this session''s median');
    if k == 1
        ylabel(axis_count(k), 'vessels reaching');
    end
end

title(layout_state, [char(relation.param.y_series) ...
    ' given lumen diameter, pooled over vessels within each arousal state'], ...
    'FontWeight', 'bold', 'FontSize', param.fontsize + 2);
save_figure(fig_state, dirs.save_dir, 'fig_fwhm_statepool');

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

function curve = reference_curve(x_value, bv_baseline, y_baseline, beta, y_series)
%REFERENCE_CURVE  eps^2 = beta*bv^2 + c, anchored so it passes through the origin.
%   The same expression fig_fwhm_diameter uses, cut down to one beta because a state
%   panel is read against the hypothesis and not against a second fit.
% IN   x_value      1xN double  the shifted abscissa, um from the operating point
%      bv_baseline  1x1 double  lumen diameter there, um
%      y_baseline   1x1 double  the plotted series there, um
%      beta         1x1 double  1 = area conserved
%      y_series     1x1 str     "eps" | "totalpvs"
% OUT  curve        1xN double  on the plotted axis, zero at x = 0
    if y_series == "totalpvs"
        eps_baseline = bv_baseline + y_baseline;
    else
        eps_baseline = y_baseline;
    end
    bv_value = bv_baseline + x_value;
    eps_value = sqrt(eps_baseline^2 + beta * (bv_value.^2 - bv_baseline^2));
    if y_series == "totalpvs"
        curve = (eps_value - bv_value) - y_baseline;
    else
        curve = eps_value - eps_baseline;
    end
end
