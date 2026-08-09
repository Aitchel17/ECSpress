function [sectormap, info] = sector_block(coremask, box_px)
% IN   coremask   logical HxW, layers grow outward from its wall
%      box_px     target box side in px; both the layer thickness and the arc a
%                 box spans are held at this
% OUT  sectormap  HxW uint32 label, ordered layer first then angle, so
%                 consecutive labels are neighbouring boxes. 0 = inside the core
%      info       n_per (boxes in each layer), n_labels, offset (label before the
%                 first of each layer), box_px, centroid, radius_map, theta_map
%
%   Layers come from sector_polar, so both partitions share one radial geometry.
%   sector_polar holds the wedge COUNT fixed, keeping angles registered across
%   rings; this one holds the BOX SIZE fixed, splitting each layer into as many
%   equal-area boxes as its arc length allows, so net flux in and out of a box is
%   comparable at every distance.

arguments
    coremask logical
    box_px   (1,1) double {mustBePositive}
end

% 0. Setup
% n_sectors = 1, so a cell IS a layer and channel 1 already equals the ring index
[layers, geo] = sector_polar(coremask, box_px, 1);
layer = layers(:,:,1);
% ring count off the map itself : sector_polar no longer echoes a grid, and Inf
% marks the region past the last ring so it has to come out of the max
ring_ch = layers(:,:,2);
n_lay   = max(ring_ch(isfinite(ring_ch)), [], 'all');
n_per = zeros(1, n_lay);
pos   = zeros(size(layer));

% 1. Splitting each layer into equal-area boxes
for L = 1:n_lay
    here = layer == L;
    npx  = nnz(here);
    if npx == 0; continue; end
    % 1.1. Choosing the count that puts each box closest to box_px squared
    n_per(L) = max(1, round(npx / box_px^2));
    % 1.2. Cutting at equal pixel counts along the angle, so the boxes come out
    %      equal-area even where the frame clips the layer
    th = geo.theta_map(here);
    [~, ord]  = sort(th);
    rank      = zeros(npx, 1);
    rank(ord) = 1:npx;
    pos(here) = min(n_per(L), floor((rank - 1) * n_per(L) / npx) + 1);
end

% 2. Fusing layer and position into one label
offset    = cumsum([0, n_per(1:end-1)]);
inband    = layer > 0;
sectormap = zeros(size(layer), 'uint32');
sectormap(inband) = uint32(offset(layer(inband))' + pos(inband));

% 3. Packing
% box_px is the caller's own argument, so it is not echoed back here
info = struct('n_per', n_per, 'n_labels', sum(n_per), 'offset', offset, ...
    'centroid', geo.centroid, ...
    'radius_map', geo.radius_map, 'theta_map', geo.theta_map);
end
