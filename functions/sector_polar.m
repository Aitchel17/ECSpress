function [sectormap, info] = sector_polar(coremask, ring_width, n_sectors)
% IN   coremask    logical HxW, r = 0 is boudary of coremask
%      ring_width  ring thickness (px)
%      n_sectors   angular bins
% OUT  sectormap   HxW uint32 label = sub2ind([n_rings n_sectors], ring, wedge).
%                  0 = no sector: inside the core or beyond the last ring
%      info        grid [n_rings n_sectors], n_labels, ring_width, ring_edges,
%                  centroid [cy cx], radius_map (px from the core wall),
%                  theta_map (rad about the centroid)

arguments
    coremask   logical
    ring_width (1,1) double {mustBePositive}
    n_sectors  (1,1) double {mustBePositive}
end

% 0. Setup
[H, W] = size(coremask);
if ~any(coremask, 'all')
    error('sector_polar:noCore', 'coremask is empty.');
end

% 1.1. Measuring distance from the core wall and angle about its centroid
[X, Y]     = meshgrid(1:W, 1:H);
radius_map = bwdist(coremask);
[ry, rx]   = find(coremask);
cx = mean(rx);
cy = mean(ry);
theta_map  = mod(atan2(Y - cy, X - cx), 2*pi);
% 1.2. Counting the rings that reach the frame edge
nr = max(1, ceil(max(radius_map(~coremask)) / ring_width));
ns = n_sectors;

% 2.1. Binning each pixel into a ring and a wedge
ring_idx = floor(radius_map / ring_width) + 1;
sect_idx = floor(theta_map / (2*pi / ns)) + 1;
% 2.2. Fusing them into one label, ordered by ring then wedge
inband    = ~coremask & ring_idx >= 1 & ring_idx <= nr;
sectormap = zeros(H, W, 'uint32');
sectormap(inband) = uint32(sub2ind([nr, ns], ring_idx(inband), sect_idx(inband)));

% 3. Packing
info = struct('grid', [nr ns], 'n_labels', nr*ns, 'ring_width', ring_width, ...
    'ring_edges', (0:nr) * ring_width, 'centroid', [cy cx], ...
    'radius_map', radius_map, 'theta_map', theta_map);
end
