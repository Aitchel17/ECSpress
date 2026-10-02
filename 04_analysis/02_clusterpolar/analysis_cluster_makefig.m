function analysis_cluster_makefig(polarcluster, roilist, pax_fwhm, t_axis, pixel2um, save_dir)
%ANALYSIS_CLUSTER_MAKEFIG  The sorted kymograph with its three picks, and each state's median images with their contours.
%   Caller: analysis_main, fig_review_session, review_session
%
% Inputs:
%   polarcluster: analysis_clusterpolar
%   roilist: roi_handle object (for contours)
%   pax_fwhm: line_fwhm
%   t_axis: Time axis vector
%   pixel2um: Conversion factor
%   save_dir: Directory to save figures

state_names  = ["constricted", "intermediate", "dilated"];
state_labels = ["Constricted", "Intermediate", "Dilated"];

%% Sorted kymograph: cluster boundaries, and the cluster each state was taken from
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

text_y = fig_paxsortedbv.loc.y(end)+1;
for state_idx = 1:numel(state_names)
    cluster_idx = polarcluster.(state_names(state_idx) + "_cluster_idx");
    text_x = fig_paxsortedbv.loc.x(polarcluster.clusterboundary(cluster_idx, 2));
    text(fig_paxsortedbv.ax, text_x, text_y, state_labels(state_idx))
end
fig_paxsortedbv.save2svg(save_dir)

%% Median images, one figure per state: BV, PVS and merged, each with its contour
% The contour colours are the pairs analysis_polar_makefig draws the profiles in
clee = color_lee;
color_bv  = {clee.clist.red,       clee.clist.pink,       clee.clist.magenta};
color_pvs = {clee.clist.darkgreen, clee.clist.lightgreen, clee.clist.green};

for state_idx = 1:numel(state_names)
    state = state_names(state_idx);
    img_bv   = double(polarcluster.stateimage("bv",  state));
    img_pvs  = double(polarcluster.stateimage("pvs", state));
    mask_bv  = get_mask_safe(roilist, char(state + "BV_contour"));
    mask_pvs = get_mask_safe(roilist, char(state + "PVS_contour"));

    fig_cluster = make_fig(char("cluster_" + state + "_images"));
    fig_cluster.update_figsize([12 4]);

    % 1. BV (Gray)
    ax_bv = subplot(1, 3, 1, 'Parent', fig_cluster.fig);
    imagesc(ax_bv, img_bv);
    axis(ax_bv, 'image', 'off');
    colormap(ax_bv, 'gray');
    title(ax_bv, 'BV');
    hold(ax_bv, 'on');
    plot_contour_overlay(ax_bv, img_bv, mask_bv, color_bv{state_idx}, false);

    % 2. PVS (Gray)
    ax_pvs = subplot(1, 3, 2, 'Parent', fig_cluster.fig);
    imagesc(ax_pvs, img_pvs);
    axis(ax_pvs, 'image', 'off');
    colormap(ax_pvs, 'gray');
    title(ax_pvs, 'PVS');
    hold(ax_pvs, 'on');
    plot_contour_overlay(ax_pvs, img_pvs, mask_pvs, color_pvs{state_idx}, false);

    % 3. Merged, yellow halo behind each contour for contrast
    ax_merged = subplot(1, 3, 3, 'Parent', fig_cluster.fig);
    rgb = plot_make_rgb(img_bv, img_pvs);
    imagesc(ax_merged, rgb);
    axis(ax_merged, 'image', 'off');
    title(ax_merged, 'Merged');
    hold(ax_merged, 'on');
    plot_contour_overlay(ax_merged, img_bv,  mask_bv,  color_bv{state_idx},  true);
    plot_contour_overlay(ax_merged, img_pvs, mask_pvs, color_pvs{state_idx}, true);

    fig_cluster.save2svg(save_dir);
end


    function plot_contour_overlay(ax, img, mask, color, use_bg)
        % Draw only the auto-detected contour that falls inside the manual
        % selection polygon (mask) -- NOT the polygon selection region itself.
        % Edge detection must match analysis_clusterpolar.drawcontour.
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
        if isempty(roilist)
            return
        end
        if any(strcmp(roilist.list(), name))
            mask = roilist.getmask(name);
        end
    end

end
