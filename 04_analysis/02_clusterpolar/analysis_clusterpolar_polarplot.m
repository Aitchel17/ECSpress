function analysis_clusterpolar_polarplot(pax_cluster, roilist, save_dir)
% ANALYSIS_CLUSTERPOLAR_POLARPLOT Generates polar plots from manual contours.
%   Uses roivtx2polar and analyze_polar to ensure alignment with PAX axis.
%
%   Inputs:
%       pax_cluster: Struct with manual_roi
%       roilist: roi_handle object (needed for PAX angle)
%       save_dir: Output directory

if ~isfield(pax_cluster, 'manual_roi') || isempty(pax_cluster.manual_roi)
    warning('pax_cluster.manual_roi is empty. Skipping polar plot.');
    return;
end

%% 1. Define Center (from Constricted BV - Index 1) and Reference Angle (PAX)
if length(pax_cluster.manual_roi) >= 1 && ~isempty(pax_cluster.manual_roi(1).SelectionMask)
    mask_cbv = pax_cluster.manual_roi(1).SelectionMask;
    img_size = size(mask_cbv);
    props = regionprops(mask_cbv, 'Centroid');
    if ~isempty(props)
        Center = props(1).Centroid; % [x, y]
    else
        Center = size(mask_cbv, [2, 1]) / 2;
    end
else
    warning('Constricted BV mask not found. Cannot determine center.');
    return;
end

% Get PAX Angle for rotation alignment
if ~isempty(roilist)
    try
        pax_vertices = roilist.getvertices('pax');
        if ~isempty(pax_vertices) && size(pax_vertices,1) >= 2
            pax_vec = pax_vertices(2,:) - pax_vertices(1,:);
            pax_angle = atan2d(pax_vec(2), pax_vec(1));
        else
            warning('PAX ROI not found or invalid in roilist. Using 0 degrees.');
            pax_angle = 0;
        end
    catch
        warning('Could not retrieve PAX vertices. Using 0 degrees.');
        pax_angle = 0;
    end
else
    pax_angle = 0;
end

%% 2. Generate Polar Maps (analyze_polar)
% We need angle_map and radius_map. analyze_polar does this.
% It needs a stack (dummy), Center, Start_angle, n_angles (dummy), mode.
dummy_stack = zeros(img_size(1), img_size(2), 1);
[~, radius_map, angle_map] = analyze_polar(dummy_stack, Center, pax_angle, 24, true);
% Note: analyze_polar returns centered_anglemap as 3rd arg.

%% 3. Process and Plot via roivtx2polar
clee = color_lee;
plot_colors = {clee.clist.red, clee.clist.darkgreen, clee.clist.magenta, clee.clist.green};
% Labels corresponding to: CBV, CPVS, DBV, DPVS

% Initialize Figure
fig = make_fig('cluster_polar_contours', 'polar');
fig.bring_fig();
fig.update_figsize([6, 6]);
hold(fig.ax, 'on');
%%
for i = 1:4
    if i > length(pax_cluster.manual_roi)
        continue;
    end

    roi = pax_cluster.manual_roi(i);
    if isempty(roi.EdgeMask)
        continue;
    end

    target_mask = roi.EdgeMask;

    [rows, cols] = find(target_mask);
    if isempty(rows)
        continue;
    end

    % Convert to vertices [x, y]
    vertices = [cols, rows];

    % Call roivtx2polar
    % Returns: [Angle(rad); Angle(deg); Radius]
    tr = roivtx2polar(vertices, angle_map, radius_map);

    angles_deg = tr(2,:);
    radii = tr(3,:);

    % Binning (0 to 359) for Distal Point Selection
    % Bin into 1-degree steps
    deg_bins = 0:360;
    binned_r = nan(1, 360);

    for k = 1:360
        bin_start = deg_bins(k);
        bin_end = deg_bins(k+1);

        % Matches in this bin
        idx = angles_deg >= bin_start & angles_deg < bin_end;

        if any(idx)
            binned_r(k) = max(radii(idx)); % Distal point
        end
    end

    % Prepare for plotting (removing NaNs usually handled by plot, but explicit is good)
    valid = ~isnan(binned_r);
    theta_centers = deg2rad(deg_bins(1:360) + 0.5); % Center of bin

    % Plotting
    % fig.plot_polar takes (theta, r)
    fig.plot_polar(theta_centers(valid), binned_r(valid), plot_colors{i}, 'x');

end
%%

title(fig.ax, 'Cluster Contours (PAX Aligned)');

if ~isempty(save_dir)
    fig.save2svg(save_dir);
end

    function import_colors
        assignin('caller', 'clee', color_lee);
    end

end
