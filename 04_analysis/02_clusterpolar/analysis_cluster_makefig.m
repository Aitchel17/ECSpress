function analysis_clusterpolar_makefig(polarcluster, roilist, pax_fwhm, t_axis, pixel2um, save_dir)
% ANALYSIS_CLUSTERPOLAR_MAKEFIG Generates figures for cluster analysis
%
% Inputs:
%   polarcluster: Struct containing cluster analysis results
%   roilist: roi_handle object (for contours)
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

fig_paxsortedbv.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed(:,polarcluster.kgph_lumen_columnidx), t_axis)

% Handle property name differences if any
if isfield(pax_fwhm.idx, 'clean_upperBVboundary')
    up_bnd = pax_fwhm.idx.clean_upperBVboundary;
    low_bnd = pax_fwhm.idx.clean_lowerBVboundary;
else
    up_bnd = pax_fwhm.idx.clean_upperboundary;
    low_bnd = pax_fwhm.idx.clean_lowerboundary;
end

fig_paxsortedbv.plot_line(up_bnd(polarcluster.kgph_lumen_columnidx),'r');
fig_paxsortedbv.plot_line(low_bnd(polarcluster.kgph_lumen_columnidx),'r');

fig_paxsortedbv.plot_xline(polarcluster.clusterboundary(:,2),'m')
fig_paxsortedbv.put_yaxistitle('Length (\mum)')
fig_paxsortedbv.put_xaxistitle('Total duration (sec)')

text_x_dilate = fig_paxsortedbv.loc.x(polarcluster.clusterboundary(polarcluster.dilated_cluster_idx, 2));
text_x_const = fig_paxsortedbv.loc.x(polarcluster.clusterboundary(polarcluster.constricted_cluster_idx, 2));

text_y = fig_paxsortedbv.loc.y(end)+1;

text(fig_paxsortedbv.ax,text_x_const,text_y,'Constricted')
text(fig_paxsortedbv.ax,text_x_dilate,text_y,'Dilated')
fig_paxsortedbv.save2svg(save_dir)

%% Prepare Contours
%% Prepare Contours from roilist
masks = struct('c_bv', [], 'c_pvs', [], 'd_bv', [], 'd_pvs', []);
if ~isempty(roilist)
    masks.c_bv = get_mask_safe(roilist, 'constrictedBV_contour');
    masks.c_pvs = get_mask_safe(roilist, 'constrictedPVS_contour');
    masks.d_bv = get_mask_safe(roilist, 'dilatedBV_contour');
    masks.d_pvs = get_mask_safe(roilist, 'dilatedPVS_contour');
end

clee = color_lee;
%% Median Images - Constricted
constricted_bv = double(polarcluster.cluster_bvcsf_constdil(:,:,1,1));
constricted_pvs = double(polarcluster.cluster_bvcsf_constdil(:,:,2,1));

fig_cluster_c = make_fig('cluster_constricted_images');
fig_cluster_c.update_figsize([12 4]);




% 1. BV (Gray)
ax1 = subplot(1, 3, 1, 'Parent', fig_cluster_c.fig);
imagesc(ax1, constricted_bv);
axis(ax1, 'image', 'off');
colormap(ax1, 'gray');
title(ax1, 'BV');
hold(ax1, 'on');
plot_contour_overlay(ax1, constricted_bv, masks.c_bv, clee.clist.red, false); % Single col

% 2. PVS (Gray)
ax2 = subplot(1, 3, 2, 'Parent', fig_cluster_c.fig);
imagesc(ax2, constricted_pvs);
axis(ax2, 'image', 'off');
colormap(ax2, 'gray');
title(ax2, 'PVS');
hold(ax2, 'on');
plot_contour_overlay(ax2, constricted_pvs, masks.c_pvs, clee.clist.darkgreen, false);

% 3. Merged
ax3 = subplot(1, 3, 3, 'Parent', fig_cluster_c.fig);
rgb_c = plot_make_rgb(constricted_bv, constricted_pvs);
imagesc(ax3, rgb_c);
axis(ax3, 'image', 'off');
title(ax3, 'Merged');
hold(ax3, 'on');
% yellow bg + color fg
plot_contour_overlay(ax3, constricted_bv, masks.c_bv, clee.clist.red, true);
plot_contour_overlay(ax3, constricted_pvs, masks.c_pvs, clee.clist.darkgreen, true);

fig_cluster_c.save2svg(save_dir);


%% Median Images - Dilated
dilated_bv = double(polarcluster.cluster_bvcsf_constdil(:,:,1,2));
dilated_pvs = double(polarcluster.cluster_bvcsf_constdil(:,:,2,2));

fig_cluster_d = make_fig('cluster_dilated_images');
fig_cluster_d.update_figsize([12 4]);

% 1. BV (Gray)
ax1 = subplot(1, 3, 1, 'Parent', fig_cluster_d.fig);
imagesc(ax1, dilated_bv);
axis(ax1, 'image', 'off');
colormap(ax1, 'gray');
title(ax1, 'BV');
hold(ax1, 'on');
plot_contour_overlay(ax1, dilated_bv, masks.d_bv, clee.clist.magenta, false);

% 2. PVS (Gray)
ax2 = subplot(1, 3, 2, 'Parent', fig_cluster_d.fig);
imagesc(ax2, dilated_pvs);
axis(ax2, 'image', 'off');
colormap(ax2, 'gray');
title(ax2, 'PVS');
hold(ax2, 'on');
plot_contour_overlay(ax2, dilated_pvs, masks.d_pvs, clee.clist.green, false);

% 3. Merged
ax3 = subplot(1, 3, 3, 'Parent', fig_cluster_d.fig);
rgb_d = plot_make_rgb(dilated_bv, dilated_pvs);
imagesc(ax3, rgb_d);
axis(ax3, 'image', 'off');
title(ax3, 'Merged');
hold(ax3, 'on');
plot_contour_overlay(ax3, dilated_bv, masks.d_bv, clee.clist.magenta, true);
plot_contour_overlay(ax3, dilated_pvs, masks.d_pvs, clee.clist.green, true);

fig_cluster_d.save2svg(save_dir);




    function plot_contour_overlay(ax, img, mask, color, use_bg)
        % Draw only the auto-detected contour that falls inside the manual
        % selection polygon (mask) -- NOT the polygon selection region itself.
        % Edge detection must match analysis_clusterpolar_contour.m.
        if isempty(mask) || isempty(img), return; end
        img_smooth = imgaussfilt(double(img), 3);
        edges = edge(img_smooth, 'log', 0.005);
        sel = edges & logical(mask);          % keep edges within the selection
        [Ey, Ex] = find(sel);
        if isempty(Ex), return; end
        % If use_bg, plot a yellow halo behind for contrast on the merged view
        if use_bg
            plot(ax, Ex, Ey, '.', 'Color', [1 1 0], 'MarkerSize', 6); % Yellow
        end
        % Foreground contour, same per-panel color
        plot(ax, Ex, Ey, '.', 'Color', color, 'MarkerSize', 3);
    end

    function mask = get_mask_safe(roilist, name)
        mask = [];
        % Try standard name only
        if any(strcmp(roilist.list(), name))
            mask = roilist.getmask(name);
        end
    end

end
