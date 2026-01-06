function analysis_clusterpolar_makefig(pax_cluster, pax_fwhm, t_axis, pixel2um, save_dir)
% ANALYSIS_CLUSTERPOLAR_MAKEFIG Generates figures for cluster analysis
%
% Inputs:
%   pax_cluster: Struct containing cluster analysis results
%   pax_fwhm: Object/Struct containing FWHM analysis results (line_fwhm)
%   t_axis: Time axis vector
%   pixel2um: Conversion factor
%   save_dir: Directory to save figures

%% Check sorting result
fig_paxsortedbv = make_fig('paxsortedBV_figure');
fig_paxsortedbv.bring_fig
fig_paxsortedbv.update_figsize([8 3])
fig_paxsortedbv.reset_axis()
fig_paxsortedbv.resolution = pixel2um;

fig_paxsortedbv.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed(:,pax_cluster.kgph_lumen_columnidx), t_axis)

% Handle property name differences if any
if isfield(pax_fwhm.idx, 'clean_upperBVboundary')
    up_bnd = pax_fwhm.idx.clean_upperBVboundary;
    low_bnd = pax_fwhm.idx.clean_lowerBVboundary;
else
    up_bnd = pax_fwhm.idx.clean_upperboundary;
    low_bnd = pax_fwhm.idx.clean_lowerboundary;
end

fig_paxsortedbv.plot_line(up_bnd(pax_cluster.kgph_lumen_columnidx),'r');
fig_paxsortedbv.plot_line(low_bnd(pax_cluster.kgph_lumen_columnidx),'r');

fig_paxsortedbv.plot_xline(pax_cluster.clusterboundary(:,2),'m')
fig_paxsortedbv.put_yaxistitle('Length (\mum)')
fig_paxsortedbv.put_xaxistitle('Total duration (sec)')

text_x_dilate = fig_paxsortedbv.loc.x(pax_cluster.clusterboundary(pax_cluster.dilated_cluster_idx, 2));
text_x_const = fig_paxsortedbv.loc.x(pax_cluster.clusterboundary(pax_cluster.constricted_cluster_idx, 2));

text_y = fig_paxsortedbv.loc.y(end)+1;

text(fig_paxsortedbv.ax,text_x_const,text_y,'Constricted')
text(fig_paxsortedbv.ax,text_x_dilate,text_y,'Dilated')
fig_paxsortedbv.save2svg(save_dir)

%% Prepare Contours
has_roi = isfield(pax_cluster, 'manual_roi') && ~isempty(pax_cluster.manual_roi);
masks = struct('c_bv', [], 'c_pvs', [], 'd_bv', [], 'd_pvs', []);
if has_roi
    % 1: CBV, 2: CPVS, 3: DBV, 4: DPVS
    if length(pax_cluster.manual_roi) >= 1, masks.c_bv = pax_cluster.manual_roi(1).EdgeMask; end
    if length(pax_cluster.manual_roi) >= 2, masks.c_pvs = pax_cluster.manual_roi(2).EdgeMask; end
    if length(pax_cluster.manual_roi) >= 3, masks.d_bv = pax_cluster.manual_roi(3).EdgeMask; end
    if length(pax_cluster.manual_roi) >= 4, masks.d_pvs = pax_cluster.manual_roi(4).EdgeMask; end
end

clee = color_lee;
%% Median Images - Constricted
constricted_bv = double(pax_cluster.cluster_bvcsf_constdil(:,:,1,1));
constricted_pvs = double(pax_cluster.cluster_bvcsf_constdil(:,:,2,1));

fig_cluster_c = make_fig('cluster_constricted_images');
fig_cluster_c.update_figsize([12 4]);




% 1. BV (Gray)
ax1 = subplot(1, 3, 1, 'Parent', fig_cluster_c.fig);
imagesc(ax1, constricted_bv);
axis(ax1, 'image', 'off');
colormap(ax1, 'gray');
title(ax1, 'BV');
hold(ax1, 'on');
plot_contour_overlay(ax1, masks.c_bv, clee.clist.red, false); % Single col

% 2. PVS (Gray)
ax2 = subplot(1, 3, 2, 'Parent', fig_cluster_c.fig);
imagesc(ax2, constricted_pvs);
axis(ax2, 'image', 'off');
colormap(ax2, 'gray');
title(ax2, 'PVS');
hold(ax2, 'on');
plot_contour_overlay(ax2, masks.c_pvs, clee.clist.darkgreen, false);

% 3. Merged
ax3 = subplot(1, 3, 3, 'Parent', fig_cluster_c.fig);
rgb_c = plot_make_rgb(constricted_bv, constricted_pvs);
imagesc(ax3, rgb_c);
axis(ax3, 'image', 'off');
title(ax3, 'Merged');
hold(ax3, 'on');
% yellow bg + color fg
plot_contour_overlay(ax3, masks.c_bv, clee.clist.red, true);
plot_contour_overlay(ax3, masks.c_pvs, clee.clist.darkgreen, true);

fig_cluster_c.save2svg(save_dir);


%% Median Images - Dilated
dilated_bv = double(pax_cluster.cluster_bvcsf_constdil(:,:,1,2));
dilated_pvs = double(pax_cluster.cluster_bvcsf_constdil(:,:,2,2));

fig_cluster_d = make_fig('cluster_dilated_images');
fig_cluster_d.update_figsize([12 4]);

% 1. BV (Gray)
ax1 = subplot(1, 3, 1, 'Parent', fig_cluster_d.fig);
imagesc(ax1, dilated_bv);
axis(ax1, 'image', 'off');
colormap(ax1, 'gray');
title(ax1, 'BV');
hold(ax1, 'on');
plot_contour_overlay(ax1, masks.d_bv, clee.clist.magenta, false);

% 2. PVS (Gray)
ax2 = subplot(1, 3, 2, 'Parent', fig_cluster_d.fig);
imagesc(ax2, dilated_pvs);
axis(ax2, 'image', 'off');
colormap(ax2, 'gray');
title(ax2, 'PVS');
hold(ax2, 'on');
plot_contour_overlay(ax2, masks.d_pvs, clee.clist.green, false);

% 3. Merged
ax3 = subplot(1, 3, 3, 'Parent', fig_cluster_d.fig);
rgb_d = plot_make_rgb(dilated_bv, dilated_pvs);
imagesc(ax3, rgb_d);
axis(ax3, 'image', 'off');
title(ax3, 'Merged');
hold(ax3, 'on');
plot_contour_overlay(ax3, masks.d_bv, clee.clist.magenta, true);
plot_contour_overlay(ax3, masks.d_pvs, clee.clist.green, true);

fig_cluster_d.save2svg(save_dir);


    function import_colors
        assignin('caller', 'clee', color_lee);
    end

    function plot_contour_overlay(ax, mask, color, use_bg)
        if isempty(mask), return; end
        [B, ~] = bwboundaries(mask);
        for k = 1:length(B)
            boundary = B{k};
            % If use_bg, plot yellow thick line first
            if use_bg
                plot(ax, boundary(:,2), boundary(:,1), 'Color', [1 1 0], 'LineWidth', 2.5); % Yellow
            end
            % Foreground line
            plot(ax, boundary(:,2), boundary(:,1), 'Color', color, 'LineWidth', 1);
        end
    end

end
