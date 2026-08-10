function [wedgemap, info] = sector_wedge(coremask, n_wedge)
%SECTOR_WEDGE  Angular partition about a mask centroid, without the rings.
%   The angular half of sector_polar, on its own.
%
%   Why it is worth separating. sector_polar fuses ring and wedge into one label
%   channel, and its ring edges are not a fixed physical scale -- 'equidistance'
%   takes its ring COUNT from the farthest in-band pixel, which sits at an image
%   corner, and 'equalarea' takes its EDGES from the sorted distance values of
%   this particular mask. Either way ring k is a different distance in a
%   different recording, so rows cannot be pooled by ring index. A caller that
%   bins radius itself on a fixed micron scale wants the wedges and the distance
%   map, and nothing else.
%
% IN   coremask   H x W bool   true inside the vessel; its boundary is r = 0
%      n_wedge    int          angular bins. Equal ANGLE, not equal area
% OUT  wedgemap   H x W double 1..n_wedge, NaN inside the core
%      info       what was measured about the partition
%        centroid          1 x 2 float   [cy cx] px, the mask's centre of mass
%        radius_map        H x W float   px from the core WALL (bwdist), 0 inside
%        theta_map         H x W float   rad about the centroid, [0 2pi)
%        wedge_edge_rad    1 x n_wedge+1 float  bin edges, derived from n_wedge
%        wedge_center_rad  1 x n_wedge float    MEASURED circular mean of the
%                          in-band pixels in each wedge. NOT the bin centre: a
%                          wedge that runs off the frame keeps only the pixels
%                          that fit, so where its measurement actually sits is
%                          not where its bin says it should
%        n_pixel           1 x n_wedge int      in-band pixels per wedge. A wedge
%                          reaching the frame edge sooner than its neighbours
%                          shows up here and nowhere else
%
%   see also SECTOR_POLAR, SECTOR_BLOCK

arguments
    coremask  logical
    n_wedge   (1,1) double {mustBePositive, mustBeInteger}
end

[H, W] = size(coremask);
if ~any(coremask, 'all')
    error('sector_wedge:noCore', 'coremask is empty.');
end

% 1. Distance from the wall, angle about the centroid. Same construction as
%    sector_polar, so the two agree pixel for pixel on the angular channel
radius_map = bwdist(coremask);
[ry, rx] = find(coremask);
cy = mean(ry);
cx = mean(rx);
[X, Y] = meshgrid(1:W, 1:H);
[theta_raw, ~] = cart2pol(X - cx, Y - cy);
theta_map = mod(theta_raw, 2*pi);

% 2. The angular bin. Equal angle, so the edges are exact and derived
wedge_width = 2*pi / n_wedge;
wedge_idx = floor(theta_map / wedge_width) + 1;
% cart2pol can return exactly 2*pi after the mod on a pixel sitting on the cut,
% which would index one past the last wedge
wedge_idx = min(wedge_idx, n_wedge);

inband = ~coremask;
wedgemap = nan(H, W);
wedgemap(inband) = wedge_idx(inband);

% 3. Where each wedge's pixels actually sit. A circular mean, because the wedge
%    that straddles theta = 0 would otherwise average to the opposite side
wedge_edge_rad = (0:n_wedge) * wedge_width;
wedge_center_rad = nan(1, n_wedge);
n_pixel = zeros(1, n_wedge);
theta_inband = theta_map(inband);
idx_inband = wedge_idx(inband);
for k = 1:n_wedge
    here = idx_inband == k;
    n_pixel(k) = nnz(here);
    if n_pixel(k) == 0
        continue
    end
    theta_here = theta_inband(here);
    mean_sin = mean(sin(theta_here));
    mean_cos = mean(cos(theta_here));
    ang = atan2(mean_sin, mean_cos);
    wedge_center_rad(k) = mod(ang, 2*pi);
end

info = struct('centroid', [cy cx], ...
    'radius_map', radius_map, ...
    'theta_map', theta_map, ...
    'wedge_edge_rad', wedge_edge_rad, ...
    'wedge_center_rad', wedge_center_rad, ...
    'n_pixel', n_pixel);
end
