%FIG_REVIEW_SESSION  One session on one board: the traces down the left, what the
%   FWHM and the cluster contours say about it down the right.
%
%   Everything here is READ and DRAWN. Nothing is measured -- the heatmap, the mode
%   curve and the annulus area were all computed by tablegeneration_fwhmrelation and
%   are read back out of this session's row; the polar profiles come from
%   analysis_clusterpolar_polarplot, which needs only the four contour ROIs.
%
%   BEFORE RUNNING, in this order:
%     02_othersignal/main_analog.m   cells 1 to the force/ball block, which is
%                                    everything before the first uiwait. That leaves
%                                    figStruct in the workspace. Without it the
%                                    ECoG, EMG and pupil cells are skipped
%     the same sessiondir in both files
%
%   Writes png and not svg. print -dsvg -vector on a board carrying the ECoG
%   spectrogram walks every pixel of it into vector geometry; see CLAUDE_LOG.md

clc
addpath('G:\03_program\01_ecspress\09_dirstruct');
dirs_ecspath;

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
param.bin_deg  = 10;     % deg    angular block the 1-deg polar bins collapse into
param.fontsize = 7;      % pt     one size over panels drawn at different figure sizes

sessiondir = 'G:\tmp\00_igkl\hql090\251016_hql090_sleep\HQL090_sleep251016_005';

%% 1. the session
session = ECSSession(sessiondir);
session = session.load_primary_results;
sleepscore = load(fullfile(sessiondir, "peripheral/sleep_score.mat"));

%% 2. FWHM panels
pax_fig = analysis_pax_makefig(session.pax_fwhm, session.pax_fwhm.t_axis, ...
    session.img_param.pixel2um, session.dir_struct.figures_fwhm);

%% 3. cluster panels
analysis_cluster_makefig(session.polarcluster, session.roilist, session.pax_fwhm, ...
    session.pax_fwhm.t_axis, session.img_param.pixel2um, ...
    session.dir_struct.figures_polarcluster);

%% 4. polar contours
% analysis_clusterpolar_polarplot gates on the four contour labels in roilist, not on
% polarcluster.manual_roi, so a session contoured in an earlier run goes straight
% through and no polygon has to be redrawn
contour_labels = ["constrictedBV_contour", "constrictedPVS_contour", ...
                  "dilatedBV_contour",     "dilatedPVS_contour"];
have_polar = false;
if all(ismember(contour_labels, session.roilist.list()))
    session.polarcluster = analysis_clusterpolar_polarplot(session.polarcluster, session.roilist);
    have_polar = isfield(session.polarcluster, 'polar_profiles');
    if have_polar
        analysis_polar_makefig(session.polarcluster, session.dir_struct.figures_polarcluster);
    end
end
fprintf('polar contours available : %d\n', have_polar);

%% 5. this session's row of fwhmrelation
dirs = dirs_central();
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);
relation = load(fullfile(dirs.save_dir, "fwhmrelation.mat")).save_content;
row = relation.session(relation.session.Directory == string(sessiondir), :);
if isempty(row)
    error('fig_review_session:noRow', '%s has no row in fwhmrelation.mat', sessiondir);
end
heat = row.heat{1};
x_label = 'lumen diameter, um from this session''s median';
y_label = [char(relation.param.y_series) ', um from this session''s median'];

%% 5a. the density, and nothing over it
fig_heat = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5 4], ...
    'Name', 'session_epsbv_heatmap');
ax_heat = axes(fig_heat);
plot_heatpanel(ax_heat, heat.x_baseceneters, heat.y_baseceneters, ...
    heat.xy_counts_clean, 'fontsize', param.fontsize);
xlabel(ax_heat, x_label)
ylabel(ax_heat, y_label)
title(ax_heat, sprintf('%s %s %s', string(row.MouseID), string(row.Date), ...
    string(row.VesselID)))

%% 5b. what area conservation demands, against what the vessel did
% Area conservation is d(eps)/d(bv) = bv/eps, so it can be evaluated at every column
% from that column's own measured bv and eps. Integrated along bv it has a closed
% form, eps = sqrt(eps0^2 + bv^2 - bv0^2), which takes ONE anchor -- the session's own
% measured median, not a fit. Black is the measured mode; no fitted beta is drawn.
bv_column   = row.bv_sessionmedian + heat.x_baseceneters;
eps_noflux  = sqrt(row.y_sessionmedian^2 + bv_column.^2 - row.bv_sessionmedian^2);
curve_noflux = eps_noflux - row.y_sessionmedian;
fprintf('mode slope %+.3f | constricted %+.3f | dilated %+.3f\n', ...
    row.mode_slope, row.slope_constricted, row.slope_dilated);

fig_ref = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5 4], ...
    'Name', 'session_noflux_curve');
ax_ref = axes(fig_ref);
hold(ax_ref, 'on')
plot(ax_ref, heat.x_baseceneters, curve_noflux, '-', 'Color', [1 0 0], 'LineWidth', 2);
plot(ax_ref, heat.x_baseceneters, heat.mode_curve(:)', '-', 'Color', [0 0 0], 'LineWidth', 2);
xlabel(ax_ref, x_label)
ylabel(ax_ref, y_label)
title(ax_ref, 'red = no flux, black = measured mode')
set(ax_ref, 'FontSize', param.fontsize, 'Box', 'off')

%% 5c. the same thing as an area
% heat.area_curve is (pi/4)(eps^2 - bv^2) along the mode curve, measured by
% heatmap_area. Area conservation is a FLAT line here, drawn at the curve's own
% value where the diameter sits at the session median
[~, at_median] = min(abs(heat.x_baseceneters));
area_baseline = heat.area_curve(at_median);

fig_area = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5 4], ...
    'Name', 'session_area_curve');
ax_area = axes(fig_area);
hold(ax_area, 'on')
plot(ax_area, heat.x_baseceneters, heat.area_curve, '-', 'Color', [0 0 0], 'LineWidth', 2);
yline(ax_area, area_baseline, 'r--', 'LineWidth', 1.2);
xlabel(ax_area, x_label)
ylabel(ax_area, 'PVS annulus area, um^2')
title(ax_area, sprintf('red = area conserved, %.1f um^2', area_baseline))
set(ax_area, 'FontSize', param.fontsize, 'Box', 'off')

%% 6. dilated radius minus constricted radius, per angle
% The quantity is polar_ave_fig's, a plain difference of radii per bin. The 1-deg
% bins are collapsed by MEDIAN into param.bin_deg blocks because a contour at r px
% has only about 2*pi*r pixels to spread over 360 bins. A block with no support
% stays NaN and polarplot breaks the line there, which is the point.
if have_polar
    n_block = 360 / param.bin_deg;
    profile = session.polarcluster.polar_profiles;
    to_block = @(v) median(reshape(v, param.bin_deg, n_block), 1, 'omitnan');
    change_bv  = to_block(profile(3).binned_r) - to_block(profile(1).binned_r);
    change_pvs = to_block(profile(4).binned_r) - to_block(profile(2).binned_r);
    fprintf('change blocks with support : BV %d of %d, PVS %d of %d\n', ...
        sum(isfinite(change_bv)), n_block, sum(isfinite(change_pvs)), n_block);

    theta_block = mean(reshape(session.polarcluster.polar_theta + ...
        deg2rad(session.polarcluster.polar_binstart_deg), param.bin_deg, n_block), 1);
    closed = @(v) [v, v(1)];
    clee = color_lee;
    fig_change = make_fig('session_polar_change', 'polar');
    fig_change.update_figsize([6 6]);
    ax_change = fig_change.ax;
    hold(ax_change, 'on');
    polarplot(ax_change, closed(theta_block), closed(1 + change_bv), '-o', ...
        'Color', clee.clist.magenta, 'LineWidth', 1.6, 'MarkerSize', 3);
    polarplot(ax_change, closed(theta_block), closed(1 + change_pvs), '-o', ...
        'Color', clee.clist.green, 'LineWidth', 1.6, 'MarkerSize', 3);
    polarplot(ax_change, deg2rad(0:360), ones(1, 361), 'k--', 'LineWidth', 0.8);
    title(ax_change, sprintf('%d deg median  (1 = no change)', param.bin_deg));
end

%% 7. the board
% 80 rows so the left panels can carry different heights and the strip can land on
% half the pupil panel: 5 + 19 + 18 + 18 + 10 + 10. The right side splits the same
% 80 into 40 + 40. No colorbars anywhere.
% the variable surviving is not enough: closing figures leaves figStruct in place
% with every axes handle deleted, and add() then fails on an invalid parent
have_analog = exist('figStruct', 'var') == 1 && ...
    all(isgraphics([figStruct.spectrogram.ax, figStruct.emg_power.ax, figStruct.pupil.ax]));
if ~have_analog
    disp('no live figStruct -- run main_analog for the ECoG, EMG and pupil panels')
end
sleepstrip = plot_sleep_strip(sleepscore);
board = fusion_figure('grid', [80 9], 'width', 18, 'height', 10.5, ...
    'fontsize', param.fontsize, 'spacing', "tight");

board.add(sleepstrip,                  'at', [ 1 1], 'span', [ 5 3], 'label', "0");
board.add(pax_fig.pvs.ax,              'at', [ 6 1], 'span', [19 3], 'label', "1");
board.add(pax_fig.totalpvs_changes.ax, 'at', [25 1], 'span', [18 3], 'label', "2");
if have_analog
    board.add(figStruct.spectrogram.ax, 'at', [43 1], 'span', [18 3], 'label', "3");
    board.add(figStruct.emg_power.ax,   'at', [61 1], 'span', [10 3], 'label', "4");
    board.add(figStruct.pupil.ax,       'at', [71 1], 'span', [10 3], 'label', "5");
end

board.add(ax_heat, 'at', [1 4], 'span', [40 2], 'label', "G");
board.add(ax_ref,  'at', [1 6], 'span', [40 2], 'label', "H");
board.add(ax_area, 'at', [1 8], 'span', [40 2], 'label', "H2");

if have_polar
    board.add(pick_axes("cluster_polar_contours", ""), 'at', [41 4], 'span', [40 2], 'label', "I1");
    board.add(ax_change,                               'at', [41 6], 'span', [40 2], 'label', "I2");
end
board.add(pick_axes("cluster_constricted_images", "Merged"), 'at', [41 8], 'span', [20 2], 'label', "J");
board.add(pick_axes("cluster_dilated_images", "Merged"),     'at', [61 8], 'span', [20 2], 'label', "K");

board.link_x(1);
% plot_sleep_strip takes its own axes past the end of the recording and link_x adopts
% that, which leaves the traces inside an empty right margin. 1800 s is 30 min.
board.panel(1).ax.XLim = [0 1800];

%% 8. write it
exportgraphics(board.fig, fullfile(dirs.save_dir, "fig_review_session.png"), ...
    'Resolution', 200);
fprintf('wrote fig_review_session.png to %s | %d panels\n', dirs.save_dir, ...
    numel(board.panel));

%% ---------------------------------------------------------------- helpers
function ax = pick_axes(figure_name, panel_title)
%PICK_AXES  One axes out of an open figure, by figure name and panel title.
%   The makefigs return nothing, so the axes they drew are reached the way
%   fusion_figure.available lists them: by walking the open figures.
%   IN   figure_name  1x1 str   the make_fig name
%        panel_title  1x1 str   "" takes the first axes on that figure
%   OUT  ax           1x1 axes | polaraxes
    fig = findall(groot, 'Type', 'figure', 'Name', figure_name);
    if isempty(fig)
        error('fig_review_session:noFigure', 'no open figure named %s', figure_name);
    end
    found = findall(fig(1), 'Type', 'axes', '-or', 'Type', 'polaraxes');
    if panel_title == ""
        ax = found(1);
        return
    end
    for k = 1:numel(found)
        if strcmpi(join(string(found(k).Title.String), " "), panel_title)
            ax = found(k);
            return
        end
    end
    error('fig_review_session:noPanel', 'no panel titled %s in %s', panel_title, figure_name);
end
