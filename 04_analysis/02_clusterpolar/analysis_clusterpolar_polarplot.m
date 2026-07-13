function polarcluster = analysis_clusterpolar_polarplot(polarcluster, roilist)
% ANALYSIS_CLUSTERPOLAR_POLARPLOT Cartesian->polar transform of manual contours.
%   Converts the 4 manual contours (CBV/CPVS/DBV/DPVS) into PAX-aligned polar
%   profiles and stores them on polarcluster for downstream plotting.
%   Rendering is done separately by analysis_polar_makefig.
%
%   Inputs:
%       polarcluster: Struct with manual_roi + cluster_bvcsf_constdil
%       roilist: roi_handle object (needed for PAX angle + contour masks)
%   Output:
%       polarcluster: same struct with added fields
%           .polar_theta    (1x360) bin-center angles [rad], shared by all conditions
%           .polar_profiles (1x4 struct) .roi_label, .binned_r (1x360, NaN=missing)

if ~isfield(polarcluster, 'manual_roi') || isempty(polarcluster.manual_roi)
    warning('polarcluster.manual_roi is empty. Skipping polar transform.');
    return;
end

if ~isfield(polarcluster, 'cluster_bvcsf_constdil')
    warning('polarcluster.cluster_bvcsf_constdil missing. Skipping polar transform.');
    return;
end

%% 1. Define Center (from Constricted BV) and Reference Angle (PAX)
target_roi_name = 'constrictedBV_contour';

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
        warning('Constricted BV mask found but empty. Using cluster data size for center.');
        sz = size(polarcluster.cluster_bvcsf_constdil);
        Center = [sz(2), sz(1)] / 2;
        img_size = [sz(1), sz(2)];
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

%% 2. Convert cartesian to polar (analyze_polar builds PAX-aligned maps)
% angle_map / radius_map are intermediate: consumed by roivtx2polar below,
% not stored (regenerable from Center + pax_angle + img_size).
dummy_stack = zeros(img_size(1), img_size(2), 1);
[~, radius_map, angle_map] = analyze_polar(dummy_stack, Center, pax_angle, 24, true);

%% 3. Per-condition transform (heavy loop) -> binned polar profiles
% Order matches the color order used in analysis_polar_makefig:
%   1: CBV, 2: CPVS, 3: DBV, 4: DPVS
state_indices = [1, 1, 2, 2]; % 1=Constricted, 2=Dilated
ch_indices    = [1, 2, 1, 2]; % 1=BV, 2=PVS
roi_names_map = { ...
    'constrictedBV_contour', ...  % 1
    'constrictedPVS_contour', ... % 2
    'dilatedBV_contour', ...      % 3
    'dilatedPVS_contour' ...      % 4
    };

theta_centers = deg2rad((0:359) + 0.5); % bin-center angles [rad], shared by all conditions

polar_profiles = struct('roi_label', {}, 'binned_r', {});
for i = 1:4
    roi_label = roi_names_map{i};

    % Contour mask for THIS condition (must be fetched per condition)
    try
        roi_mask = roilist.getmask(roi_label);
    catch
        roi_mask = [];
    end

    % Raw image + edge detection (same params as contour script)
    img = polarcluster.cluster_bvcsf_constdil(:, :, ch_indices(i), state_indices(i));
    try
        edges = edge(imgaussfilt(double(img), 3), 'log', 0.005);
    catch
        edges = false(size(img));
    end

    % Missing mask or no edge pixels -> all-NaN profile (kept so array length stays 4)
    if isempty(roi_mask)
        polar_profiles(i) = struct('roi_label', roi_label, 'binned_r', nan(1, 360));
        continue;
    end
    [rows, cols] = find(edges & roi_mask);
    if isempty(rows)
        polar_profiles(i) = struct('roi_label', roi_label, 'binned_r', nan(1, 360));
        continue;
    end

    % Cartesian vertices [x, y] -> polar via pre-built maps
    tr = roivtx2polar([cols, rows], angle_map, radius_map); % [rad; deg; radius]
    binned_r = bin_distal(tr(2, :), tr(3, :));

    polar_profiles(i) = struct('roi_label', roi_label, 'binned_r', binned_r);
end

%% 4. Store on polarcluster (small vectors; persisted via existing polarcluster save)
polarcluster.polar_theta    = theta_centers;
polarcluster.polar_profiles = polar_profiles;
% Offset to ADD to polar_theta to rotate back to the ORIGINAL IMAGE frame:
% analyze_polar centered the map as centered = mod(theta_map - bin_start, 360),
% so theta_map (true image-screen angle) = centered + bin_start. angle_range =
% 360/n_angles with n_angles = 24 (as passed to analyze_polar above).
polarcluster.polar_binstart_deg = mod(-pax_angle, 180) - (360/24)/2;
% PAX line angle (deg) for the makefig title + direction marker
polarcluster.pax_angle = pax_angle;

end


function binned_r = bin_distal(angles_deg, radii)
% Bin angles into 1-degree steps (0..359) and take the distal (max radius)
% point per bin. Empty bins stay NaN so the polar line breaks at missing
% angular sectors instead of connecting across them.
deg_bins = 0:360;
binned_r = nan(1, 360);
for k = 1:360
    idx = angles_deg >= deg_bins(k) & angles_deg < deg_bins(k+1);
    if any(idx)
        binned_r(k) = max(radii(idx)); % distal point
    end
end
end
