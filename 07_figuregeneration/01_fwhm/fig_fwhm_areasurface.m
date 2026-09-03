%FIG_FWHM_AREASURFACE  The annulus area as a surface over the two radii.
%   A = pi(r_eps^2 - r_bv^2) is a hyperbolic paraboloid, and the two hypotheses
%   are two curves ON it: area conservation is a level curve, the measured one is
%   the same curve opened by the fitted beta_ratio. Drawn as a surface because the
%   distance between them is the area given up, which a slope-against-slope plot
%   cannot show.
%     caution  the trajectories used to be drawn STRAIGHT, which is the level
%              curve's tangent at the operating point stretched over the whole
%              walk and doubles the apparent area loss. Both are curves now

clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
% Style
clee = color_lee(); sty.fontsize = 11; sty.view = [-38 26];



param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
dirs = dirs_central();
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);


%% The surface
gen.r_max = 15; % um       
gen.n_grid = 220; gen.n_level = 14;
r_bv = linspace(0, gen.r_max, gen.n_grid); r_eps = linspace(0, gen.r_max, gen.n_grid);
[grid_bv, grid_eps] = meshgrid(r_bv, r_eps);
% pvs area
pvs_area = pi * (grid_eps.^2 - grid_bv.^2); pvs_area(grid_eps < grid_bv) = NaN; % non-negative

%% Where the vessels actually sit, and the two ways out of there
% Read from the table, not typed in. Rebuild fwhmrelation and this figure follows;
% a number here would have to be kept in step with it by hand.
relation = load(fullfile(dirs.save_dir, 'fwhmrelation.mat')).save_content;
drawn = relation.session(relation.session.drawn, :);
param.walk = 4;             % um       how far along bv each trajectory is drawn

% the surface is drawn in RADII and the table holds diameters
start_bv = median(drawn.bv_sessionmedian) / 2;
start_eps = median(drawn.y_sessionmedian) / 2;
beta_ratio = median(drawn.beta_ratio);
bvwalk = linspace(start_bv - param.walk, start_bv + param.walk / 2, 120);
start_area = pi * (start_eps^2 - start_bv^2);

% ONE expression for both trajectories: eps^2 = beta*bv^2 + c, anchored on the
% operating point. beta = 1 is the level curve that keeps the area, and the fitted
% beta_ratio is what the vessels did. Drawing the measured one straight, as this
% file used to, is that curve's tangent stretched over the whole walk.
noflux_epswalk = sqrt(start_eps^2 + (bvwalk.^2 - start_bv^2));
noflux_area = pi * (noflux_epswalk.^2 - bvwalk.^2);
slope_noflux = start_bv / start_eps;       % the level curve slope AT the operating point

measured_eps = sqrt(start_eps^2 + beta_ratio * (bvwalk.^2 - start_bv^2));
measured_area = pi * (measured_eps.^2 - bvwalk.^2);
slope_measured = beta_ratio * slope_noflux;


%%

fprintf('surface : %d x %d over 0..%g um, area 0..%.0f um^2\n', gen.n_grid, ...
    gen.n_grid, gen.r_max, max(pvs_area(:)));

fprintf('operating point : r_bv %.2f, r_eps %.2f um, area %.1f um^2\n', ...
    start_bv, start_eps, start_area);
fprintf('   slope of the level curve there : %.4f  (the vessels'' median %.4f)\n', ...
    start_bv / start_eps, slope_noflux);
fprintf('   after %+.1f um of DIAMETER : no flux %.1f, measured %.1f, given up %.1f um^2\n', ...
    param.walk, noflux_area(end), measured_area(end), noflux_area(end) - measured_area(end));

%% Draw
fig = figure('Color', 'w', 'Position', [60 50 1180 520], 'Name', 'fig_fwhm_areasurface');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

axis_surface = nexttile(layout);
surf(axis_surface, grid_bv, grid_eps, pvs_area, 'EdgeColor', 'none', 'FaceAlpha', 0.92);
colormap(axis_surface, clee.gradient.inferno);
hold(axis_surface, 'on')

ridge = plot3(axis_surface, r_bv, r_bv, zeros(size(r_bv)), '-', ...
    'Color', [0 0 0], 'LineWidth', 2.5);
noflux_line = plot3(axis_surface, bvwalk, noflux_epswalk, noflux_area, '-', ...
    'Color', [0.35 0.85 1], 'LineWidth', 3);
measured_line = plot3(axis_surface, bvwalk, measured_eps, measured_area, '-', ...
    'Color', [0.3 1 0.5], 'LineWidth', 3);
plot3(axis_surface, start_bv, start_eps, start_area, 'o', 'MarkerSize', 9, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

xlabel(axis_surface, 'r_{bv}, um')
ylabel(axis_surface, 'r_{\epsilon}, um')
zlabel(axis_surface, 'PVS area, um^2')
view(axis_surface, sty.view)
grid(axis_surface, 'on')
set(axis_surface, 'FontSize', sty.fontsize, 'Box', 'on')
axis_surface.Toolbar.Visible = 'off';
legend(axis_surface, [ridge, noflux_line, measured_line], ...
    ["r_\epsilon = r_{bv}, the sheath closed onto the wall", ...
     sprintf('no flux, the level curve (%.3f here)', start_bv / start_eps), ...
     sprintf('measured %.3f', slope_measured)], ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', sty.fontsize - 2)
title(axis_surface, 'A = \pi(r_\epsilon^2 - r_{bv}^2)')

%% The same thing from above, where the level curves are the point
axis_plan = nexttile(layout);
contourf(axis_plan, grid_bv, grid_eps, pvs_area, gen.n_level, 'LineColor', [1 1 1 ]*0.4);
colormap(axis_plan, clee.gradient.inferno);
hold(axis_plan, 'on')
plot(axis_plan, r_bv, r_bv, '-', 'Color', [0 0 0], 'LineWidth', 2.5);
plot(axis_plan, bvwalk, noflux_epswalk, '-', 'Color', [0.35 0.85 1], 'LineWidth', 3);
plot(axis_plan, bvwalk, measured_eps, '-', 'Color', [0.3 1 0.5], 'LineWidth', 3);
plot(axis_plan, start_bv, start_eps, 'o', 'MarkerSize', 9, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
xlabel(axis_plan, 'r_{bv}, um')
ylabel(axis_plan, 'r_{\epsilon}, um')
axis(axis_plan, 'equal')
xlim(axis_plan, [0 gen.r_max])
ylim(axis_plan, [0 gen.r_max])
set(axis_plan, 'FontSize', sty.fontsize, 'Box', 'on', 'Layer', 'top')
axis_plan.Toolbar.Visible = 'off';
colorbar(axis_plan, 'eastoutside');
title(axis_plan, 'from above: a contour is a constant PVS area')

title(layout, ['the annulus area over the two radii, and the two ways out of one vessel' ...
    newline 'the measured slope crosses contours downhill, which is the area the PVS gives up'])
save_figure(fig, dirs.save_dir, 'fig_fwhm_areasurface');

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
