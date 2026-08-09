function [div_cells, anchor_xy] = vfield_divergence(flow, sectormap, n_labels)
% IN   flow       HxWx2 field, last dim 1 = u (horizontal), 2 = v (vertical);
%                 HxWxNx2 is summed over time first
%      sectormap  HxW label map from sector_block / sector_polar
%      n_labels   how many labels the partition DEFINES. sector_block gives an
%                 exact count; for sector_polar it is n_rings * n_sectors.
%                 [] = read it off the map, which comes up SHORT when the
%                 highest-numbered cell is empty
% OUT  div_cells  n_labels x 1 divergence du/dx+dv/dy, NaN where the sector held
%                 no fittable samples. Units per-frame (px displacement / px)
%      anchor_xy  n_labels x 2 [x y] centroid of each sector's samples
%
%   Computation only: it fits, it does not gate and it does not paint. Prune with
%   sector_validate first, spread with sector_paint after.

arguments
    flow      {mustBeNumeric, mustBeNonempty}
    sectormap {mustBeNumeric}
    n_labels  double = []
end

% 0. Setup
% 0.1. Splitting the flow into u and v
if ndims(flow) == 4
    u = sum(flow(:,:,:,1), 3);
    v = sum(flow(:,:,:,2), 3);
else
    f = squeeze(flow);
    u = f(:,:,1);
    v = f(:,:,2);
end
[H, W] = size(u);
lab    = double(sectormap);
if ~isequal(size(lab), [H, W])
    error('vfield_divergence:mapSize', 'sectormap must be %dx%d.', H, W);
end
% 0.2. Sizing the outputs
if isempty(n_labels); n_labels = max(lab(:)); end
div_cells = NaN(n_labels, 1);
anchor_xy = NaN(n_labels, 2);
[X, Y]    = meshgrid(1:W, 1:H);
sample    = ~isnan(u) & ~isnan(v);

% 1. Fitting one least-squares plane per sector
for k = 1:n_labels
    sel = sample & lab == k;
    if ~any(sel, 'all'); continue; end
    xi = X(sel);  yi = Y(sel);
    % 1.1. vfield_divergence_plane returns ok = false on a rank-deficient spread
    [dval, ok] = vfield_divergence_plane(xi, yi, u(sel), v(sel));
    if ~ok; continue; end
    div_cells(k)    = dval;
    anchor_xy(k, :) = [mean(xi), mean(yi)];
end
end
