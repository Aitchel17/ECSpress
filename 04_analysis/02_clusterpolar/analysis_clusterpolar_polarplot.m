function analysis_clusterpolar_polarplot(polarcluster, roilist, save_dir)
% ANALYSIS_CLUSTERPOLAR_POLARPLOT Generates polar plots from manual contours.
%   Uses roivtx2polar and analyze_polar to ensure alignment with PAX axis.
%
%   Inputs:
%       polarcluster: Struct with manual_roi
%       roilist: roi_handle object (needed for PAX angle)
%       save_dir: Output directory

if ~isfield(polarcluster, 'manual_roi') || isempty(polarcluster.manual_roi)
    warning('polarcluster.manual_roi is empty. Skipping polar plot.');
    return;
end

%% 1. Define Center (from Constricted BV) and Reference Angle (PAX)
target_roi_name = 'constrictedBV_contour';
if ~any(strcmp(roilist.list(), target_roi_name))
    % Check legacy or alternative naming if needed, but per instruction we enforce new naming.
    % Actually, if we just renamed it in the contour script, we might not have it yet if the user hasn't run it?
    % But we assume the pipeline will be run.
    % warning('constrictedBV_contour not found. Center calculation might fail.');
end

try
    mask_cbv = roilist.getmask(target_roi_name);
    if ~isempty(mask_cbv)
        img_size = size(mask_cbv);
        props = regionprops(mask_cbv, 'Centroid');
        if ~isempty(props)
            Center = props(1).Centroid; % [x, y]
        else
            Center = size(mask_cbv, [2, 1]) / 2;
        end
    else
        warning('Constricted BV mask found but empty. Using defaults.');
        % We don't have image size here easily if mask is empty, but analyze_polar needs Center.
        % Assuming standard size or failing?
        % Let's try to get image size from pax_cluster if possible
        if isfield(polarcluster, 'cluster_bvcsf_constdil')
            sz = size(polarcluster.cluster_bvcsf_constdil);
            Center = [sz(2), sz(1)] / 2;
            img_size = [sz(1), sz(2)];
        else
            error('Cannot determine center: Constricted BV mask empty and no cluster data size.');
        end
    end
catch
    warning('Constricted BV mask not found in roilist. Cannot determine center.');
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
% Iterate through 4 conditions to plot:
% 1: CBV (State 1, Ch 1) -> Green
% 2: CPVS (State 1, Ch 2) -> Dark Green
% 3: DBV (State 2, Ch 1) -> Magenta
% 4: DPVS (State 2, Ch 2) -> Green (Note: Original code used Green for 4th, likely intended separate color? Preserving logic)

% Map for Colors/Order matching original script:
% Index 1 (CBV) -> Color{1} (Red) ?? Wait, verify original colors.
% Original script: plot_colors = {red, darkgreen, magenta, green}
% 1: CBV -> Red
% 2: CPVS -> DarkGreen
% 3: DBV -> Magenta
% 4: DPVS -> Green

state_indices = [1, 1, 2, 2]; % 1=Constricted, 2=Dilated
ch_indices    = [1, 2, 1, 2]; % 1=BV, 2=PVS
roi_names_map = { ...
    'constrictedBV_contour', ...  % 1
    'constrictedPVS_contour', ... % 2
    'dilatedBV_contour', ...      % 3
    'dilatedPVS_contour' ...      % 4
    };

% Ensure cluster data exists
if ~isfield(polarcluster, 'cluster_bvcsf_constdil')
    warning('polarcluster.cluster_bvcsf_constdil missing.');
    return;
end

for i = 1:4
    roi_label = roi_names_map{i};
    st_idx = state_indices(i);
    ch_idx = ch_indices(i);

    % 1. Get raw image (for edge detection)
    img = polarcluster.cluster_bvcsf_constdil(:,:,ch_idx,st_idx);

    % 2. Re-run Edge Detection (Same parameters as contour script)
    try
        img_smooth = imgaussfilt(double(img), 3);
        edges = edge(img_smooth, 'log', 0.005);
    catch
        edges = false(size(img));
    end

    % 3. Get Selection Mask from ROILIST (The Polygon)
    try
        % Find ROI in roilist
        % Note: setup_rois creates "manual_" prefixed names.
        % We need to check if the user is using "manual_" prefix map or standard names.
        % The contour script generated names like "constricted_bv" (without manual_ prefix??)
        % Wait, contour script used: contour_name = sprintf('%s_%s', state_names{state_idx}, ch_names{ch_idx});
        % which results in "constricted_bv", "dilated_pvs", etc.
        % setup_rois used "manual_constricted_bv".
        % Check if "constricted_bv" exists, if not try "manual_constricted_bv" or vice versa.

        target_roi_name = roi_label;
        if ~any(strcmp(roilist.list(), target_roi_name))
            % No fallback to manual_
            % continue;
        end

        roi_vertices = roilist.getvertices(target_roi_name);
        roi_mask = roilist.getmask(target_roi_name); % Need getmask method or recreate mask
        % roilist.getmask is available based on previous read of roi_handle.m

        if isempty(roi_mask)
            % If mask is empty but vertices exist, recreate?
            % roi_handle.getmask returns obj.ROIs(i).Mask
            continue;
        end
    catch
        % ROI not found
        continue;
    end

    % 4. Filter Edges
    target_mask = edges & roi_mask;

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



end
