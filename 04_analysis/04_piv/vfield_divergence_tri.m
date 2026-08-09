function [div_cells, info] = vfield_divergence_tri(flow, sectormap, n_labels, opt)
%VFIELD_DIVERGENCE_TRI  Divergence per sector, by triangulating the vectors.
%   Delaunay-triangulates the scattered vectors, takes the gradient on each
%   triangle, then aggregates to the sector by area.
%
%   Why this and not a finite difference. On a triangle the linear interpolant of
%   (u,v) has a CONSTANT gradient, so du/dx + dv/dy there is exact for that
%   interpolant -- no stencil, no grid spacing, nothing to choose. And the
%   area-weighted sum over a set of triangles equals the boundary flux of their
%   union, because every interior edge is traversed twice with opposite normals
%   and cancels (the divergence theorem). So aggregating small triangles IS the
%   integral form, and it comes out of the same code as the local values.
%
%   That is what makes it worth the trouble: differentiating scattered data
%   amplifies noise, integrating it averages the noise down, and this gives the
%   integral for free while still reporting per-triangle values.
%
%   Delaunay will happily connect vectors across a hole, and those slivers span
%   ground where nothing was measured -- across the masked vessel, or across a
%   patch the gates emptied. Their gradients are ill-conditioned and their area
%   is large, so they would dominate an area-weighted sum. Hence the two filters
%   below; both are geometry, neither looks at the values.
%
% IN   flow       H x W x 2 displacement field, NaN off the vector grid
%      sectormap  H x W label map from sector_polar / sector_block
%      n_labels   how many labels the partition DEFINES. sector_block gives an
%                 exact count; for sector_polar it is n_rings * n_sectors.
%                 [] = read it off the map, which comes up SHORT when the
%                 highest-numbered cell is empty
%      max_edge   longest triangle edge kept (px). [] = 3x the median edge, which
%                 on a regular vector grid keeps the neighbours and drops the
%                 bridges over holes
%      min_angle  smallest interior angle kept (deg), default 15. A sliver has a
%                 nearly singular gradient
%      exclmask   H x W, true = never interpolate here; a triangle whose centroid
%                 lands inside is dropped, default []
% OUT  div_cells  n_labels x 1 area-weighted divergence, NaN where the sector
%                 held no usable triangle
%      info       tri (kept connectivity), xy, div_tri, area_tri, label_tri,
%                 area_cells, n_tri_cells, dropped (counts per filter), plus
%                 grad_tri (nTri x 4 [du/dx du/dy dv/dx dv/dy]) and cxy_tri
%                 (nTri x 2 centroids). div_tri is grad_tri(:,1)+grad_tri(:,4);
%                 the other two are there so a caller can build the strain TENSOR
%                 and take a directional component instead of the trace

arguments
    flow      {mustBeNumeric, mustBeNonempty}
    sectormap {mustBeNumeric}
    n_labels  double = []
    opt.max_edge  double = []
    opt.min_angle (1,1) double = 15
    opt.exclmask  = []
end

% 0. Setup
u = flow(:,:,1);   v = flow(:,:,2);
lab = double(sectormap);
if isempty(n_labels); n_labels = max(lab(:)); end

% 0.1 The vectors, as scattered points
[yi, xi] = find(~isnan(u) & ~isnan(v));
idx = sub2ind(size(u), yi, xi);
xy  = [xi, yi];
uu  = u(idx);   vv = v(idx);
if numel(uu) < 3
    div_cells = NaN(n_labels, 1);
    info = struct('tri', [], 'xy', xy, 'div_tri', [], 'area_tri', [], ...
        'label_tri', [], 'area_cells', zeros(n_labels,1), ...
        'n_tri_cells', zeros(n_labels,1), 'dropped', struct());
    return
end

% 1. Triangulate
DT  = delaunayTriangulation(xy);
tri = DT.ConnectivityList;
P1  = xy(tri(:,1), :);   P2 = xy(tri(:,2), :);   P3 = xy(tri(:,3), :);

% 2. Geometry filters
% 2.1 Edge lengths; the default cut is relative, so it follows the grid step
e1 = vecnorm(P2-P1, 2, 2);   e2 = vecnorm(P3-P2, 2, 2);   e3 = vecnorm(P1-P3, 2, 2);
maxe = max([e1 e2 e3], [], 2);
lim  = opt.max_edge;
if isempty(lim); lim = 3 * median(min([e1 e2 e3], [], 2)); end
keep_edge = maxe <= lim;
% 2.2 Smallest interior angle, by the law of cosines on the shortest side
ang = @(a,b,c) acosd(min(1, max(-1, (b.^2 + c.^2 - a.^2) ./ (2*b.*c))));
amin = min([ang(e1,e3,e2), ang(e2,e1,e3), ang(e3,e2,e1)], [], 2);
keep_ang = amin >= opt.min_angle;
% 2.3 Centroid, used to assign a sector
cx = (P1(:,1) + P2(:,1) + P3(:,1)) / 3;
cy = (P1(:,2) + P2(:,2) + P3(:,2)) / 3;
ci = sub2ind(size(u), round(cy), round(cx));
% 2.4 Mask test on all three VERTICES, not the centroid. A triangle can have its
%     centroid clear of the mask while one corner sits on the boundary, and that
%     corner's window is shared with the masked region, so its vector carries the
%     boundary's bias straight into the fit
keep_mask = true(size(cx));
if ~isempty(opt.exclmask)
    vx = [P1(:,1) P2(:,1) P3(:,1)];
    vy = [P1(:,2) P2(:,2) P3(:,2)];
    vi = sub2ind(size(u), round(vy), round(vx));
    keep_mask = ~any(opt.exclmask(vi), 2);
end

keep = keep_edge & keep_ang & keep_mask;
dropped = struct('long_edge', nnz(~keep_edge), 'sliver', nnz(~keep_ang), ...
                 'masked', nnz(~keep_mask), 'kept', nnz(keep), 'total', numel(keep));

tri = tri(keep, :);
P1 = P1(keep,:);  P2 = P2(keep,:);  P3 = P3(keep,:);
ci = ci(keep);   cx = cx(keep);   cy = cy(keep);

% 3. Exact gradient of the linear interpolant on each triangle
%    2A is the signed cross product; the sign cancels between numerator and
%    denominator, so orientation does not matter
twoA = (P2(:,1)-P1(:,1)).*(P3(:,2)-P1(:,2)) - (P3(:,1)-P1(:,1)).*(P2(:,2)-P1(:,2));
u1 = uu(tri(:,1)); u2 = uu(tri(:,2)); u3 = uu(tri(:,3));
v1 = vv(tri(:,1)); v2 = vv(tri(:,2)); v3 = vv(tri(:,3));
dudx = (u1.*(P2(:,2)-P3(:,2)) + u2.*(P3(:,2)-P1(:,2)) + u3.*(P1(:,2)-P2(:,2))) ./ twoA;
dvdy = (v1.*(P3(:,1)-P2(:,1)) + v2.*(P1(:,1)-P3(:,1)) + v3.*(P2(:,1)-P1(:,1))) ./ twoA;
% the other two cost nothing -- same interpolant, the other partial. With all four
% the full gradient is available, and a caller that wants a DIRECTIONAL strain
% (eps_rr and eps_tt are n'*E*n for the radial and tangential n) needs the tensor,
% not the trace. See vfield_strain method 'radial_tri'
dudy = (u1.*(P3(:,1)-P2(:,1)) + u2.*(P1(:,1)-P3(:,1)) + u3.*(P2(:,1)-P1(:,1))) ./ twoA;
dvdx = (v1.*(P2(:,2)-P3(:,2)) + v2.*(P3(:,2)-P1(:,2)) + v3.*(P1(:,2)-P2(:,2))) ./ twoA;
div_tri  = dudx + dvdy;
area_tri = abs(twoA) / 2;

% 4. Aggregate to sectors, weighting by area. This is the flux integral over the
%    sector: the interior edges cancel, so only its boundary survives
label_tri  = lab(ci);
div_cells  = NaN(n_labels, 1);
area_cells = zeros(n_labels, 1);
n_tri      = zeros(n_labels, 1);
ok = label_tri > 0 & label_tri <= n_labels & isfinite(div_tri);
if any(ok)
    w = accumarray(label_tri(ok), area_tri(ok),               [n_labels 1]);
    s = accumarray(label_tri(ok), area_tri(ok).*div_tri(ok),  [n_labels 1]);
    n_tri = accumarray(label_tri(ok), 1,                      [n_labels 1]);
    area_cells = w;
    has = w > 0;
    div_cells(has) = s(has) ./ w(has);
end

info = struct('tri', tri, 'xy', xy, 'div_tri', div_tri, 'area_tri', area_tri, ...
    'label_tri', label_tri, 'area_cells', area_cells, 'n_tri_cells', n_tri, ...
    'dropped', dropped, 'max_edge', lim, ...
    'grad_tri', [dudx, dudy, dvdx, dvdy], ...   % nTri x 4, the full gradient
    'cxy_tri',  [cx, cy]);                      % nTri x 2, where each one applies
end
