function [sectormap, info] = sector_polar(coremask, ring_size, n_sectors, mode)
% input
%      coremask    logical HxW, r = 0 is boudary of coremask
%      ring_size   'equidistance' : ring thickness (px)
%                  'equalarea'    : pixels per ring
%      n_sectors   angular bins
%      mode        'equidistance' fixed thickness, ring count follows the frame
%                  'equalarea'    fixed pixel count, thickness shrinks outward
% output
%      sectormap   HxW x 3 double, one channel per way of naming a cell
%                    (:,:,1) label  unique per cell, 0 = no cell
%                    (:,:,2) ring   1..n_rings, NaN inside the core
%                    (:,:,3) wedge  1..n_sectors, NaN inside the core
%      info        ring_edges,
%                  centroid [cy cx],
%                  radius_map (px from the core wall),
%                  theta_map (rad about the centroid)
%
%   why  three channels and not one label : fusing (ring, wedge) into a single
%        number creates an ordering that every consumer then has to re-derive, and
%        getting it backwards transposes the result silently. Channels 2 and 3
%        answer "which ring / which wedge" without anyone knowing the ordering,
%        and channel 1 keeps the label the accumulation helpers need. The ordering
%        is now an implementation detail of THIS file
%   note  channel 1 uses 0 for "no cell" because sector_paint, sector_validate and
%        vfield_divergence all skip 0. Channels 2 and 3 use NaN instead, so
%        isfinite(sectormap(:,:,2)) is the cell test and 0 is never ambiguous
%        with ring 0
%   caution  TWO different counts, and they are not interchangeable
%          max(sectormap(:,:,1))   how far a `for k = 1:...` has to go. Correct
%                                  for looping: labels above it do not exist
%          n_rings * n_sectors     how many cells the GRID has. Required to
%                                  reshape a per-cell vector into rings x wedges,
%                                  and larger than the first whenever the last
%                                  cell holds no pixels -- see FINDINGS.md
%        n_rings is max(sectormap(:,:,2)), which is always exact because the
%        farthest pixel is always in the last ring. n_sectors is the caller's own
%        argument, so neither is returned in info
%   note  the two axes have DIFFERENT origins on purpose : radius is bwdist from
%        the core WALL, so rings are wall-parallel bands that follow a
%        non-circular vessel, while theta is about the mask CENTROID. A cell is
%        "this far out from the wall, in this direction from the centre", not an
%        annulus of a circle
%   caution  n_rings is not a free choice, it is ceil(max distance / ring_width)
%        and that maximum sits at an IMAGE CORNER. Change the crop and the ring
%        count changes on identical data, and the outermost rings can carry no
%        vectors at all -- see FINDINGS.md

arguments
    coremask  logical
    ring_size (1,1) double {mustBePositive}
    n_sectors (1,1) double {mustBePositive}
    mode      (1,:) char {mustBeMember(mode, {'equidistance','equalarea'})} = 'equidistance'
end

% 0. Setup
[H, W] = size(coremask);
if ~any(coremask, 'all')
    error('sector_polar:noCore', 'coremask is empty.');
end

% 1.1. Measuring distance from the core wall and angle about its centroid
radius_map = bwdist(coremask);
%% Theta map generation, Center of mask > grid matrix
[ry, rx]   = find(coremask);
cx = mean(rx);
cy = mean(ry);
[X, Y]     = meshgrid(1:W, 1:H);
[theta_map, ~] = cart2pol(X - cx, Y - cy);
theta_map = mod(theta_map, 2*pi);
inband = ~coremask;

% 1.2. Ring edges, and with them the ring count. The two modes differ ONLY here
switch mode
    case 'equidistance'
        % Fixed thickness. n_rings follows the frame, so the crop decides it.
        n_rings    = max(1, ceil(max(radius_map(inband)) / ring_size));
        ring_edges = (0:n_rings) * ring_size;
        % min() closes the last ring on its upper edge, as histcounts does.
        % see CLAUDE_LOG.md -- without it sub2ind stops with IndexOutOfRange
        ring_idx = min(floor(radius_map / ring_size) + 1, n_rings);

    case 'equalarea'
        % Fixed pixel count per ring, so thickness shrinks outward and every ring
        % carries comparable n. Ported from analyze_polar, which binned this way
        % per wedge; here the edges are global so a ring is still one radius band
        % across all wedges and rings stay comparable to each other.
        radius_sorted = sort(radius_map(inband), 'ascend');
        n_rings       = max(1, floor(numel(radius_sorted) / ring_size));
        edge_at       = ring_size * (1:n_rings);
        ring_edges    = [0, double(radius_sorted(edge_at)).'];
        % the remainder above the last cut joins the last ring rather than
        % falling out as NaN, which is what min() does in the other mode
        ring_edges(end) = double(max(radius_sorted));
        ring_idx        = discretize(radius_map, ring_edges);
        ring_idx(isnan(ring_idx)) = n_rings;
end

% 2.1. The angular bin. Independent of the mode
sect_idx = floor(theta_map / (2*pi / n_sectors)) + 1;
% note  ring_idx >= 1 needs no test : bwdist is non-negative, so floor(r/w)+1 is
%       at least 1 everywhere in the frame

% 2.2. Channel 1, the fused label. sub2ind fixes the ordering HERE and nowhere
%      else -- nothing outside this file needs to know it is ring-major
ring_inband  = ring_idx(inband);
wedge_inband = sect_idx(inband);
label_ch = zeros(H, W);
label_ch(inband) = sub2ind([n_rings, n_sectors], ring_inband, wedge_inband);

% 2.3. Channels 2 and 3, the two indices on their own
% note  there is no "past the last ring" region. n_rings is derived FROM the
%       farthest non-core pixel, so the rings cover the whole frame by
%       construction and NaN (inside the core) is the only sentinel that occurs.
%       An Inf branch here would be guarding a case that cannot happen; it only
%       becomes reachable if n_rings ever turns into a caller's argument
ring_ch  = nan(H, W);
wedge_ch = nan(H, W);
ring_ch(inband)  = ring_idx(inband);
wedge_ch(inband) = sect_idx(inband);

sectormap = cat(3, label_ch, ring_ch, wedge_ch);

% 3. Packing. Only what was MEASURED -- ring_width and n_sectors are the caller's
%    own arguments and echoing them here would hide where they were decided
info = struct('ring_edges', ring_edges, ...
    'centroid', [cy cx], 'radius_map', radius_map, 'theta_map', theta_map);
end
