function [div_map, div_cells, info] = vfield_divergence_polar(flow, coremask, opt)
%VFIELD_DIVERGENCE_POLAR  Perivascular divergence in polar (ring x sector) wedges.
%   Composition: sector_polar partitions, sector_validate prunes, vfield_divergence
%   fits, sector_paint spreads.
%   Adds what only makes sense in polar coordinates -- the mean radial and
%   tangential velocity per wedge, and the 'radial' method d<v_r>/dr along the
%   radius. Units per-frame (px displacement / px), same as vfield_divergence.
%   A wedge clipped by the frame edge still gives the right v_r and v_theta:
%   each pixel is projected before averaging, so partial angular coverage does
%   not bias them (measured on a synthetic pure-radial field, 17% coverage,
%   0.0% error). Only a wedge's MEAN DISPLACEMENT VECTOR is coverage-dependent,
%   which is why none is reported.
%
% IN   flow        H x W x 2 (or H x W x 1 x 2) displacement field, NaN off-grid
%      coremask    logical H x W, the blob rings start at. This pipeline passes
%                  manual_dilated_pvs, so r = 0 is the PVS outer wall
%      n_sectors   angular bins, default 4
%      ring_width  ring thickness (px), default 20
%      min_valid   vectors required per cell, default 4
%      exclmask    logical H x W extra exclusion, default []
%      method      which measure div_map / div_cells carry. NOTE that only
%                  'plane' and 'divergence' are actually divergences; the file
%                  name is historical and the others are single strain-tensor
%                  components:
%                    'plane'       du/dx+dv/dy per wedge, one least-squares plane.
%                                  A divergence, and assumes no axisymmetry
%                    'radial'      eps_rr = d<v_r>/dr, centred difference between
%                                  ring means. RADIAL STRAIN, not a divergence
%                    'radial_fit'  the same quantity, least squares over the
%                                  vectors in a radial window instead of two ring
%                                  means. Use this one: a ring is 10 px wide and
%                                  the scatter inside it runs larger than the step
%                                  between rings
%                    'hoop'        eps_tt = v_r/r + (1/r) dv_theta/dtheta.
%                                  CIRCUMFERENTIAL STRAIN
%                    'divergence'  eps_rr + eps_tt, the trace. Areal strain
%                  info carries all of them whatever is chosen
% OUT  div_map     H x W dartboard, every pixel painted with its cell value
%      div_cells   n_rings x n_sectors, NaN where invalid
%      info        centroid [cy cx], ring_width, n_rings, n_sectors, counts,
%                  ring_edges, cell_xy, vr_cells (+ = outward), vt_cells
%                  (+ = counter-clockwise), rcen, vr_map, vt_map

arguments
    flow          {mustBeNumeric, mustBeNonempty}
    coremask      logical
    opt.n_sectors (1,1) double = 4
    opt.ring_width (1,1) double = 20
    opt.min_valid (1,1) double = 4
    opt.exclmask  logical = []
    opt.method (1,:) char {mustBeMember(opt.method, ...
        {'plane','radial','radial_fit','hoop','divergence'})} = 'plane'
    opt.fit_halfwidth (1,1) double = 15   % 'radial_fit' radial half-window (px)
    opt.fit_min       (1,1) double = 6    % vectors a 'radial_fit' window needs
end

% 0. Setup
% 0.1 Split the flow into u and v
if ndims(flow) == 4
    u = sum(flow(:,:,:,1), 3);
    v = sum(flow(:,:,:,2), 3);
else
    f = squeeze(flow);
    u = f(:,:,1);
    v = f(:,:,2);
end
[H, W] = size(u);
if ~isequal(size(coremask), [H, W])
    error('vfield_divergence_polar:maskSize', 'coremask must be %dx%d.', H, W);
end
% 1.1. Partitioning into wedges -- geometry only, identical for every event
% sector_polar returns H x W x 3 : label, ring, wedge. Everything below works on
% the LABEL channel, because sector_validate / vfield_divergence / sector_paint
% all accumulate by label. Channels 2-3 are for callers that select by (ring,wedge)
[part, geo] = sector_polar(coremask, opt.ring_width, opt.n_sectors);
wedges  = part(:,:,1);
ring_ch = part(:,:,2);
nr = max(ring_ch(isfinite(ring_ch)), [], 'all');
ns = opt.n_sectors;
% the GRID cell count, which is what every reshape below needs. Not the same as
% max(wedges(:)) : that is only how far a for-loop has to go, and it comes up
% short whenever the last cell holds no pixels
n_gridcell = nr * ns;
% 1.2. Dropping excluded pixels from the partition
if ~isempty(opt.exclmask)
    wedges(opt.exclmask) = 0;
end
% 1.3. Clearing the wedges that cannot carry the reported measure. 'plane' solves
%      a + bx + cy per component, so it needs the samples to span two rows and two
%      columns; 'radial' only averages a scalar, and forcing that rank condition on
%      it drops a third of the wedges for nothing. Measured on ev2: the drop is
%      uneven across rings, which by itself moved ring 2 from 1.135 to 1.374 px
min_span = 2;
if strcmp(opt.method, 'radial')
    min_span = 1;
end
sample = ~isnan(u) & ~isnan(v);
[wedges, counts] = sector_validate(wedges, sample, n_gridcell, opt.min_valid, min_span);

% 2. Fitting one plane per wedge
[plane_cells, anchor_xy] = vfield_divergence(flow, wedges, n_gridcell);

% 3. Per-wedge radial quantities
% 3.1. Projecting onto the radial and tangential directions
vrad_field = u .* cos(geo.theta_map) + v .* sin(geo.theta_map);
vtan_field = -u .* sin(geo.theta_map) + v .* cos(geo.theta_map);
% 3.2. Averaging both, and distance-from-wall, per wedge
lab      = double(wedges);
vr_cells = NaN(n_gridcell, 1);
vt_cells = NaN(n_gridcell, 1);
rcen     = NaN(n_gridcell, 1);
for k = 1:n_gridcell
    sel = sample & lab == k;
    if ~any(sel, 'all')
        continue
    end
    vr_cells(k) = mean(vrad_field(sel));
    vt_cells(k) = mean(vtan_field(sel));
    rcen(k)     = mean(geo.radius_map(sel));
end
vr_cells = reshape(vr_cells, nr, ns);
vt_cells = reshape(vt_cells, nr, ns);
rcen     = reshape(rcen, nr, ns);

% 3.3. The strain components, so nothing downstream has to assemble them.
%      eps_rr = d<v_r>/dr is ONE COMPONENT of the strain tensor, not a
%      divergence; the divergence is its trace, eps_rr + eps_tt. Both terms of
%      eps_tt are kept: v_r/r, which is positive wherever the field points
%      outward and so is geometrically forced, and (1/r) dv_theta/dtheta, which
%      is not. r is measured from the CENTROID for these, since that is the
%      radius the polar form divides by -- radius_map is distance from the wall,
%      so r0 (the core's equivalent-circle radius) is added back. That circle is
%      an approximation: the mask is not round, and its error is worst near the
%      wall, exactly where eps_tt is largest
r0      = sqrt(nnz(coremask)/pi);
r_cen   = rcen + r0;
eps_rr  = radial_gradient(vr_cells, rcen);
hoop_r  = vr_cells ./ r_cen;
hoop_th = tangential_gradient(vt_cells, r_cen);
eps_tt  = hoop_r + hoop_th;

% 4. Pick the reported measure, then spread it over the frame
switch opt.method
    case 'plane';  div_cells = reshape(plane_cells, nr, ns);
    case 'radial'; div_cells = eps_rr;                            % d<v_r>/dr
    case 'hoop';   div_cells = eps_tt;                            % v_r/r + ...
    case 'divergence'
        % The trace. Positive = the annulus gains area, which is what the
        % gliovascular model calls tangential stretch winning over radial
        % compression
        div_cells = eps_rr + eps_tt;
    case 'radial_fit'
        % Least squares over the individual vectors in a radial window, instead
        % of differencing two ring means. A ring is 10 px wide and the scatter
        % inside it runs larger than the step between rings, so the two-point
        % difference is dominated by that scatter; a window of 3 rings fitted at
        % once is the same quantity, conditioned.
        sect_idx  = floor(geo.theta_map / (2*pi/ns)) + 1;
        div_cells = radial_fit(vrad_field, geo.radius_map, sect_idx, sample, ...
            rcen, nr, ns, opt.fit_halfwidth, opt.fit_min);
end
div_map = sector_paint(wedges, div_cells(:));

% 5. Pack
% note  ring_width is not echoed back : it came in as opt.ring_width and the
%       caller still has it. ring_edges below IS derived, so it stays
info = struct('centroid', geo.centroid, ...
    'n_rings', nr, 'n_sectors', ns, 'counts', reshape(counts, nr, ns), ...
    'ring_edges', geo.ring_edges, 'cell_xy', reshape(anchor_xy, nr, ns, 2), ...
    'vr_cells', vr_cells, 'vt_cells', vt_cells, 'rcen', rcen, 'r_cen', r_cen, ...
    'eps_rr', eps_rr, 'eps_tt', eps_tt, 'hoop_r', hoop_r, 'hoop_th', hoop_th, ...
    'divergence', eps_rr + eps_tt, ...
    'vr_map', sector_paint(wedges, vr_cells(:)), ...
    'vt_map', sector_paint(wedges, vt_cells(:)));
end

% ---------------------------------------------------------------------------
function g = tangential_gradient(vt_all, r_all)
%TANGENTIAL_GRADIENT  (1/r) d<v_theta>/dtheta, the second term of eps_tt.
%   Centred in theta and wrapping, since the sectors close on themselves. Small
%   on average -- measured on ev2 its median was -0.00018 against +0.0059 for
%   v_r/r -- but per cell it reaches a third of that term, so dropping it is not
%   free. NaN where either neighbour is missing.
    [nr, ns] = size(vt_all);
    dth = 2*pi/ns;
    g   = NaN(nr, ns);
    for ir = 1:nr
        for is = 1:ns
            im = mod(is-2, ns) + 1;
            ip = mod(is,   ns) + 1;
            if any(isnan([vt_all(ir,im), vt_all(ir,ip), r_all(ir,is)]))
                continue
            end
            g(ir, is) = (vt_all(ir,ip) - vt_all(ir,im)) / (2*dth) / r_all(ir,is);
        end
    end
end

% ---------------------------------------------------------------------------
function g = radial_fit(vrad, rmap, sect, sample, rcen, nr, ns, halfw, minn)
%RADIAL_FIT  d<v_r>/dr as the slope of v_r against distance, per sector.
%   One straight-line fit per (ring, sector) over every vector of that SECTOR
%   lying within halfw px of the ring, so the fit crosses ring boundaries and
%   sees a real spread in r. NaN where the window is too thin to fit.
    g = NaN(nr, ns);
    for is = 1:ns
        here = sample & sect == is;
        if ~any(here, 'all')
            continue
        end
        d = rmap(here);   y = vrad(here);
        for ir = 1:nr
            r0 = rcen(ir, is);
            if isnan(r0)
                continue
            end
            in = abs(d - r0) <= halfw;
            if nnz(in) < minn || range(d(in)) < halfw
                continue
            end
            p = [ones(nnz(in), 1), d(in)] \ y(in);
            g(ir, is) = p(2);
        end
    end
end

% ---------------------------------------------------------------------------
function g = radial_gradient(vp_all, rp_all)
%RADIAL_GRADIENT  d<v_r>/dr per sector: centred difference, else forward, else backward.
    [nr, ns] = size(vp_all);
    g = NaN(nr, ns);
    for is = 1:ns
        vp = vp_all(:, is);   rp = rp_all(:, is);
        for ir = 1:nr
            lo = ir - 1;  hi = ir + 1;
            if hi <= nr && lo >= 1 && ~any(isnan([vp(lo) vp(hi) rp(lo) rp(hi)])) && rp(hi) ~= rp(lo)
                g(ir, is) = (vp(hi) - vp(lo)) / (rp(hi) - rp(lo));
            elseif hi <= nr && ~any(isnan([vp(ir) vp(hi) rp(ir) rp(hi)])) && rp(hi) ~= rp(ir)
                g(ir, is) = (vp(hi) - vp(ir)) / (rp(hi) - rp(ir));
            elseif lo >= 1 && ~any(isnan([vp(ir) vp(lo) rp(ir) rp(lo)])) && rp(ir) ~= rp(lo)
                g(ir, is) = (vp(ir) - vp(lo)) / (rp(ir) - rp(lo));
            end
        end
    end
end
