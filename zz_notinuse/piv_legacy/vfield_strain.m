function [strain_map, strain_cells, info] = vfield_strain(uv, opt)
%VFIELD_STRAIN  One u,v map in, a strain or divergence map out, by any of three routes.
%   The three ways of differentiating a PIV field disagree, and the disagreement
%   is the measurement -- so they take the SAME input and give the SAME three
%   outputs, and the only thing that changes between them is method=.
%
%     for m = ["fd" "tri" "radial_fit"]
%         [map, cells, info] = vfield_strain(uv, method = m, coremask = core);
%     end
%
%   Router, not an engine. It normalises the input, picks the partition when the
%   method needs one and the caller has not given one, and merges the engine's
%   own info into the returned struct so nothing is buried a level down. Every
%   number comes from vfield_divergence_fd / _tri / _polar unchanged.
%
%   WHICH ONE TO BELIEVE
%     fd          differentiates point by point, so it AMPLIFIES the vector noise.
%                 Cheapest, no geometry, and the one to distrust on a sparse grid
%     tri         area-weighted over triangles = the boundary flux of their union
%                 (the divergence theorem: every interior edge cancels). Integrates
%                 rather than differentiates, so it averages the noise DOWN
%     radial_fit  not a divergence at all -- eps_rr = d<v_r>/dr, one component of
%                 the strain tensor. The only one that isolates PVS compression
%   see FINDINGS.md
%
%   'radial_tri' NEEDS TWO THINGS, AND WITHOUT EITHER IT LOOKS LIKE NOISE
%
%   1. COARSE WEDGES. A per-triangle gradient is the noisiest estimate there is --
%      three neighbouring vectors differenced. What redeems it is the area-weighted
%      sum over the cell, and a cell with a handful of triangles has nothing to
%      average. 'tri' gets away with it because box_partition defaults to 4x the
%      vector spacing; a 10 px ring x 30 deg wedge is far smaller
%      see FINDINGS.md
%      note      info.n_tri_cells is the triangle count. Check it before going finer
%
%   2. DROP THE WALL-ADJACENT RING. Ring 1 is unstable because it is thin, not
%      because its windows are clipped
%      err       an earlier version of this comment blamed the coremask : the mask
%                is excluded from the correlation, so a window centred just outside
%                the wall would have half its area blanked. That was checked and is
%                NOT what happens -- only 10 of 1423 vectors have >=25% core
%                overlap, and 78 of ring 1's 108 vectors never touch the PVS at all
%      see FINDINGS.md
%      err       including it is what made eps_rr look like a population trend.
%                Its sign is opposite to rings 2+, so a median across 0-30 px
%                cancels the signal against the noise
%
%   WITH BOTH, on 32 events and 43 matched-pair controls, scored by how far the
%   dilations and constrictions separate in units of the control scatter:
%     estimator                band        effect   sign ok   binom p
%     radial_fit  ring 10x12   1-3 (old)     0.19     16/32     0.57
%     radial_fit  ring 10x12   2-6           1.96     29/32     1.3e-06
%     radial_tri  ring 10x12   2-4           0.71     26/32     2.7e-04
%     radial_tri  ring 20x 6   1-3           2.33     29/32     1.3e-06
%     radial_tri  ring 20x 6   2-4           2.64     31/32     7.7e-09
%   see FINDINGS.md
%   caution   many bands were tried. Bonferroni over ~20 gives 1.5e-7, and the
%             coarse route holds 1.64-2.64 across EVERY band that excludes ring 1,
%             so this is not one lucky combination
%   note      for scale, wall motion v_r scores 4.00 on the same rows. eps_rr is
%             within a factor of 2 of it, not the order of magnitude an earlier
%             pass of this comment claimed
%   err  comparing fd against radial without saying so compares a trace against a
%        single component, and they are not the same quantity
%
% IN   uv         H x W x 2 float  (:,:,1) = u, (:,:,2) = v, NaN off the vector
%                 grid. H x W x N x 2 is summed over N first
%      method     which quantity, and by which route:
%                   'fd'          du/dx+dv/dy, finite difference on the grid
%                   'tri'         du/dx+dv/dy per triangle, area-weighted per cell
%                   'radial_tri'  eps_rr per TRIANGLE, area-weighted onto the polar
%                                 wedges. The strain tensor route. READ THE WEDGE
%                                 SIZE NOTE BELOW BEFORE USING IT             [core]
%                   'plane'       du/dx+dv/dy per polar wedge, one LSQ plane   [core]
%                   'radial'      eps_rr, centred difference between rings      [core]
%                   'radial_fit'  eps_rr, LSQ over a radial window              [core]
%                   'hoop'        eps_tt = v_r/r + (1/r) dv_th/dth              [core]
%                   'divergence'  eps_rr + eps_tt, areal strain                 [core]
%      coremask   H x W bool  r = 0 wall. REQUIRED by every method marked [core]
%      exclmask   H x W bool  true = never use this pixel
%      sectormap  H x W int   'tri' only; [] builds a square box partition
%      box_px     float px    'tri' box side when sectormap is []. NaN = 4 x the
%                             vector spacing, which puts 8-32 triangles in a cell
%      n_sectors / ring_width / min_valid    polar partition
%      max_edge / min_angle                  'tri' geometry filters
% OUT  strain_map    H x W float    NaN where nothing was measured
%      strain_cells  per region; [] for 'fd', which has no regions
%      info          struct  method, unit, sectormap, n_labels, plus every field
%                    the engine reported, merged in flat
%
% caution  units are px displacement per px, PER EVENT for an endpoint field --
%          analysis_pivensemble has already multiplied by result.scale
% caution  'fd' tiles its answer out by nearest neighbour, so its map looks dense
%          and is not: the resolution is still the vector spacing

arguments
    uv             {mustBeNumeric, mustBeNonempty}
    opt.method     (1,:) char {mustBeMember(opt.method, ...
        {'fd','tri','plane','radial','radial_fit','hoop','divergence', ...
         'radial_tri'})} = 'fd'
    opt.coremask   logical = []
    opt.exclmask   logical = []
    opt.sectormap  {mustBeNumeric} = []
    opt.box_px     (1,1) double = NaN
    opt.n_sectors  (1,1) double = 12
    opt.ring_width (1,1) double = 10
    opt.min_valid  (1,1) double = 5
    opt.max_edge   double = []
    opt.min_angle  (1,1) double = 15
end

% 0. Setup
% 0.1 One H x W x 2 field, whatever came in
if ndims(uv) == 4
    flow = cat(3, sum(uv(:,:,:,1), 3), sum(uv(:,:,:,2), 3));
else
    flow = squeeze(uv);
end
[H, W, ~] = size(flow);
% vfield_gradient_tri works on the interrogation grid now, so the two tri
% methods hand it the vectors with their own coordinates. This function's own
% input is still the dense map, and rebuilding the grid from it recovers the
% coordinates only to the pixel -- exactly what the dense form threw away. The
% numbers therefore do not move, which is the point while this file is a
% candidate for archiving; a caller that has the grid should use vfield_polar
% see PIV_PLAN.md 5.2b
[fy, fx] = find(~isnan(flow(:,:,1)) & ~isnan(flow(:,:,2)));
grid_idx = sub2ind([H W], fy, fx);
u_flat = flow(:,:,1);
v_flat = flow(:,:,2);
xyuv_from_dense = cat(3, fx(:), fy(:), u_flat(grid_idx), v_flat(grid_idx));
is_polar = ismember(opt.method, ...
    {'plane','radial','radial_fit','hoop','divergence','radial_tri'});
if is_polar && isempty(opt.coremask)
    error('vfield_strain:noCoremask', ...
        'method ''%s'' measures distance from a wall; pass coremask.', opt.method);
end

% 1. Route
switch opt.method
    case 'fd'
        % 1.1 No partition, no geometry. divergence() on the sparse grid
        [strain_map, div_mean, div_sum] = vfield_divergence_fd(flow, ...
            exclmask = opt.exclmask);
        strain_cells = [];
        info = struct('mean', div_mean, 'sum', div_sum, ...
            'sectormap', [], 'n_labels', 0);

    case 'tri'
        % 1.2 Square boxes when the caller has none. Big enough that a cell holds
        %     several triangles: one triangle is a bare gradient, a cell of them
        %     is the flux through that cell's boundary, which is the whole point
        sectormap = double(opt.sectormap);
        if isempty(sectormap)
            box = opt.box_px;
            if isnan(box)
                box = 4 * vfield_grid_step(flow);
            end
            sectormap = box_partition(H, W, box);
        end
        n_labels = max(sectormap(:));
        % vfield_gradient_tri returns the triangles and does not aggregate, so
        % the label per triangle and the area weighting are done here
        [tri, tri_info] = vfield_gradient_tri(xyuv_from_dense, ...
            max_edge = opt.max_edge, min_angle = opt.min_angle, ...
            exclmask = opt.exclmask);
        label_tri = tri_labels(sectormap, tri.cxy);
        strain_cells = accumulate_tri(label_tri, tri.area, tri.divergence, n_labels);
        strain_map = sector_paint(sectormap, strain_cells);
        info = tri_info;
        info.tri       = tri.conn;
        info.div_tri   = tri.divergence;
        info.label_tri = label_tri;
        info.sectormap = sectormap;
        info.n_labels  = n_labels;

    case 'radial_tri'
        % 1.3 eps_rr from the strain TENSOR, per triangle, onto the polar wedges.
        %     why  on a triangle the linear interpolant's gradient is CONSTANT, so
        %          all four partials are exact there and eps_rr = r'*E*r needs no
        %          stencil and no radial window. 'radial_fit' has to choose one
        %     why  area-weighted onto a wedge is the same integral 'tri' takes for
        %          the trace, so the noise averages down instead of amplifying
        %     err  eps_rr is a COMPONENT, not the trace. Using div_tri for it would
        %          silently answer eps_rr + eps_tt
        % channel 1 is the fused label, which is what tri_labels reads and
        % sector_paint accumulates on. Channels 2-3 are not needed here
        [part, geo] = sector_polar(opt.coremask, opt.ring_width, opt.n_sectors);
        wedges  = part(:,:,1);
        ring_ch = part(:,:,2);
        % ring count off the map, since sector_polar no longer echoes a grid.
        % n_sectors is the caller's own argument and stays theirs
        n_rings = max(ring_ch(isfinite(ring_ch)), [], 'all');
        % the GRID cell count, which is what the reshapes below need. NOT
        % max(wedges(:)) -- that is only how far a for-loop has to go, and it
        % comes up short whenever the last cell holds no pixels
        n_gridcell = n_rings * opt.n_sectors;
        if ~isempty(opt.exclmask)
            wedges(opt.exclmask) = 0;
        end
        [tri, tri_info] = vfield_gradient_tri(xyuv_from_dense, ...
            max_edge = opt.max_edge, min_angle = opt.min_angle, ...
            exclmask = opt.exclmask);
        label_tri = tri_labels(wedges, tri.cxy);
        % 1.3.1 Radial unit vector at each triangle, from the core's centroid
        rx = tri.cxy(:,1) - geo.centroid(2);
        ry = tri.cxy(:,2) - geo.centroid(1);
        rn = max(hypot(rx, ry), eps);
        rx = rx ./ rn;   ry = ry ./ rn;
        % 1.3.2 eps_rr = r' * E * r, E the symmetric part of the gradient
        g       = tri.grad;                             % [du/dx du/dy dv/dx dv/dy]
        shear   = (g(:,2) + g(:,3)) / 2;
        eps_rr  = g(:,1).*rx.^2 + 2*shear.*rx.*ry + g(:,4).*ry.^2;
        eps_tt  = g(:,1).*ry.^2 - 2*shear.*rx.*ry + g(:,4).*rx.^2;
        % 1.3.3 Area-weighted onto the wedges
        strain_cells = accumulate_tri(label_tri, tri.area, eps_rr, n_gridcell);
        strain_map   = sector_paint(wedges, strain_cells);
        info = tri_info;
        info.tri        = tri.conn;
        info.div_tri    = tri.divergence;
        info.label_tri  = label_tri;
        info.eps_rr_tri = eps_rr;
        info.eps_tt_tri = eps_tt;
        info.eps_tt     = accumulate_tri(label_tri, tri.area, eps_tt, n_gridcell);
        info.centroid   = geo.centroid;
        info.n_rings    = n_rings;
        info.sectormap  = part;
        info.n_gridcell = n_gridcell;
        strain_cells    = reshape(strain_cells, n_rings, opt.n_sectors);
        info.eps_tt     = reshape(info.eps_tt,  n_rings, opt.n_sectors);

    otherwise
        % 1.4 Rings x wedges from the wall outward. vfield_divergence_polar owns
        %     which of its five measures comes back in div_cells
        [strain_map, strain_cells, polar_info] = vfield_divergence_polar(flow, ...
            opt.coremask, method = opt.method, n_sectors = opt.n_sectors, ...
            ring_width = opt.ring_width, min_valid = opt.min_valid, ...
            exclmask = opt.exclmask);
        info = polar_info;
        info.sectormap = [];
        info.n_labels  = numel(strain_cells);
end

% 2. Stamp what was asked for onto the answer, so a saved result says what it is
info.method = opt.method;
switch opt.method
    case {'fd','tri','plane'};  info.unit = 'du/dx + dv/dy   (px per px, per event)';
    case {'radial','radial_fit','radial_tri'}
                                info.unit = 'eps_rr   (px per px, per event)';
    case 'hoop';                info.unit = 'eps_tt   (px per px, per event)';
    case 'divergence';          info.unit = 'eps_rr + eps_tt   (px per px, per event)';
end
end

% ---------------------------------------------------------------------------
function label = tri_labels(labelmap, cxy)
%TRI_LABELS  Which partition cell each triangle falls in, read at its centroid.
%   The linear interpolant's gradient is constant over a triangle, so the
%   centroid is where that gradient applies and it is the point the label is
%   read at. A triangle straddling two cells therefore lands wholly in one of
%   them, which is the assumption the area weighting in accumulate_tri rests on.
%
%   vfield_divergence_tri used to do this internally. vfield_gradient_tri
%   returns the triangles and stops there, so the labelling moved out to the
%   caller that owns the partition. See CLAUDE_LOG.md
    cx = round(cxy(:,1));
    cy = round(cxy(:,2));
    ci = sub2ind(size(labelmap), cy, cx);
    label = labelmap(ci);
end

% ---------------------------------------------------------------------------
function cells = accumulate_tri(label, area, value, n_labels)
%ACCUMULATE_TRI  Area-weighted mean of a per-triangle quantity, per label.
%   The weighting is what makes this an integral over the cell rather than an
%   average of gradients: a sliver contributes what its area is worth.
    cells = NaN(n_labels, 1);
    ok = label > 0 & label <= n_labels & isfinite(value);
    if ~any(ok)
        return
    end
    w = accumarray(label(ok), area(ok),              [n_labels 1]);
    s = accumarray(label(ok), area(ok).*value(ok),   [n_labels 1]);
    has = w > 0;
    cells(has) = s(has) ./ w(has);
end

% ---------------------------------------------------------------------------
function lab = box_partition(H, W, box_px)
%BOX_PARTITION  Square boxes over the whole frame, labelled row-major.
%   The partition of last resort: no wall, no centre, nothing to choose but the
%   size. Rows first, so consecutive labels are horizontal neighbours.
    box = max(1, round(box_px));
    [X, Y] = meshgrid(1:W, 1:H);
    col = floor((X - 1) / box);
    row = floor((Y - 1) / box);
    lab = row * (floor((W - 1) / box) + 1) + col + 1;
end
