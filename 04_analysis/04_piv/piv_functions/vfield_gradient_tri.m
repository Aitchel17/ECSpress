function [tri, info] = vfield_gradient_tri(xyuv, opt)
%VFIELD_GRADIENT_TRI  Displacement gradient per triangle of the scattered vectors.
%   Delaunay-triangulates the vectors and takes the exact gradient of the linear
%   interpolant on each triangle. Returns the triangles; it does not aggregate.
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
%   Why the full gradient and not the trace. eps_rr and eps_tt are n'*E*n for the
%   radial and tangential n, so a caller wanting a DIRECTIONAL strain needs the
%   tensor. The two off-diagonal partials come out of the same interpolant and
%   cost nothing. Taking the hoop strain as t'*E*t also removes the need for a
%   radius: u_r/r is the axisymmetric form and wants r from the AXIS, which a
%   distance-from-the-wall map does not have. See CLAUDE_LOG.md
%
%   Why it no longer aggregates. Binning belongs to whoever owns the partition:
%   a label map fixes ring edges to one mask and one crop, which is exactly what
%   stops two recordings from pooling. vfield_polar bins these triangles on a
%   fixed micron scale instead. See CLAUDE_LOG.md
%
% IN   xyuv       ny x nx x 4 float px, (:,:,1:2) = [x y] window centre,
%                 (:,:,3:4) = [u v] displacement. NaN in u or v drops that window.
%                 ny x nx x N x 4 is N consecutive intervals: the displacement
%                 sums over N, the grid does not move
%                 The vectors arrive scattered and stay scattered. The dense
%                 H x W form used to be the input, and the first thing this
%                 function did was find(~isnan) to get these coordinates back --
%                 at pixel resolution, because the stamp had rounded them.
%                 See PIV_PLAN.md 5.2b
%      max_edge   longest triangle edge kept (px). [] = 3x the median shortest
%                 edge, which on a regular vector grid keeps the neighbours and
%                 drops the bridges over holes
%      min_angle  smallest interior angle kept (deg). A sliver has a nearly
%                 singular gradient
%      exclmask   H x W bool, true = never interpolate here. Tested on all three
%                 VERTICES, not the centroid
% OUT  tri        struct of COLUMNS, every field n_tri x 1 (or n_tri x 2/4), so
%                 one selection applies to all of them at once:
%        conn        n_tri x 3 int    vertex indices into xy
%        grad        n_tri x 4 float  [du/dx du/dy dv/dx dv/dy], per px
%        divergence  n_tri x 1 float  grad(:,1) + grad(:,4). Dimensionless
%        area        n_tri x 1 float  px^2
%        cxy         n_tri x 2 float  [x y] px, centroid, where it applies
%        uvc         n_tri x 2 float  [u v] px, the displacement AT that centroid.
%                    Exact, not a resample: the interpolant is linear, so its
%                    value there is the mean of the three vertices
%      info       what was measured about the triangulation
%        xy                n_vec x 2 float  the scattered vectors, [x y] px
%        n_vector          int      vectors that went in
%        n_tri_kept        int      triangles that survived both filters
%        max_edge_measured float px 3x the median shortest edge. ALWAYS the
%                          measured value, never an echo of opt.max_edge -- the
%                          limit that was applied is the caller's own argument
%                          and it still has it. See CLAUDE_LOG.md
%        dropped           struct   counts per filter
%
%   see also VFIELD_POLAR, VFIELD_DIVERGENCE_FD

arguments
    xyuv          {mustBeNumeric, mustBeNonempty}
    opt.max_edge  double = []
    opt.min_angle (1,1) double = 15
    opt.exclmask  logical = []
end

% 0. Setup
% 0.1 Split the four pages. Five dimensions is N consecutive intervals: the
%     displacement sums over them, the grid does not move
if ndims(xyuv) == 5
    x = xyuv(:,:,1,1);
    y = xyuv(:,:,1,2);
    u = sum(xyuv(:,:,:,3), 3);
    v = sum(xyuv(:,:,:,4), 3);
else
    x = xyuv(:,:,1);
    y = xyuv(:,:,2);
    u = xyuv(:,:,3);
    v = xyuv(:,:,4);
end

% 0.2 The vectors, as scattered points. They arrive scattered -- the window
%     centres are the coordinates, at whatever sub-pixel position PIVlab put them
have = ~isnan(u) & ~isnan(v);
xy = [x(have), y(have)];
uu = u(have);
vv = v(have);

% 0.3 Too few to triangulate. Same field set as the normal path, so a caller
%     never has to test which way it returned
if numel(uu) < 3
    tri = empty_tri();
    info = struct('xy', xy, 'n_vector', numel(uu), 'n_tri_kept', 0, ...
        'max_edge_measured', NaN, 'dropped', empty_dropped());
    return
end

% 1. Triangulate
DT = delaunayTriangulation(xy);
conn = DT.ConnectivityList;
P1 = xy(conn(:,1), :);
P2 = xy(conn(:,2), :);
P3 = xy(conn(:,3), :);

% 2. Geometry filters
% 2.1 Edge lengths; the default cut is relative, so it follows the grid step
e1 = vecnorm(P2 - P1, 2, 2);
e2 = vecnorm(P3 - P2, 2, 2);
e3 = vecnorm(P1 - P3, 2, 2);
edge_max = max([e1 e2 e3], [], 2);
edge_min = min([e1 e2 e3], [], 2);
max_edge_measured = 3 * median(edge_min);
lim = opt.max_edge;
if isempty(lim)
    lim = max_edge_measured;
end
keep_edge = edge_max <= lim;

% 2.2 Smallest interior angle, by the law of cosines on the shortest side
ang = @(a,b,c) acosd(min(1, max(-1, (b.^2 + c.^2 - a.^2) ./ (2*b.*c))));
a1 = ang(e1, e3, e2);
a2 = ang(e2, e1, e3);
a3 = ang(e3, e2, e1);
angle_min = min([a1 a2 a3], [], 2);
keep_angle = angle_min >= opt.min_angle;

% 2.3 Centroid, where the triangle's one gradient applies
cx = (P1(:,1) + P2(:,1) + P3(:,1)) / 3;
cy = (P1(:,2) + P2(:,2) + P3(:,2)) / 3;

% 2.4 Mask test on all three VERTICES, not the centroid. A triangle can have its
%     centroid clear of the mask while one corner sits on the boundary, and that
%     corner's window is shared with the masked region, so its vector carries the
%     boundary's bias straight into the gradient
keep_mask = true(size(cx));
if ~isempty(opt.exclmask)
    vx = [P1(:,1) P2(:,1) P3(:,1)];
    vy = [P1(:,2) P2(:,2) P3(:,2)];
    % the mask is the only thing here that lives in frame coordinates, so it is
    % also the only thing that knows the frame size
    vi = sub2ind(size(opt.exclmask), round(vy), round(vx));
    on_mask = opt.exclmask(vi);
    keep_mask = ~any(on_mask, 2);
end

keep = keep_edge & keep_angle & keep_mask;
dropped = struct('long_edge', nnz(~keep_edge), 'sliver', nnz(~keep_angle), ...
                 'masked', nnz(~keep_mask), 'kept', nnz(keep), 'total', numel(keep));

conn = conn(keep, :);
P1 = P1(keep, :);
P2 = P2(keep, :);
P3 = P3(keep, :);
cx = cx(keep);
cy = cy(keep);

% 3. Exact gradient of the linear interpolant on each triangle
%    2A is the signed cross product; the sign cancels between numerator and
%    denominator, so orientation does not matter
twoA = (P2(:,1)-P1(:,1)).*(P3(:,2)-P1(:,2)) - (P3(:,1)-P1(:,1)).*(P2(:,2)-P1(:,2));
u1 = uu(conn(:,1));
u2 = uu(conn(:,2));
u3 = uu(conn(:,3));
v1 = vv(conn(:,1));
v2 = vv(conn(:,2));
v3 = vv(conn(:,3));
dudx = (u1.*(P2(:,2)-P3(:,2)) + u2.*(P3(:,2)-P1(:,2)) + u3.*(P1(:,2)-P2(:,2))) ./ twoA;
dudy = (u1.*(P3(:,1)-P2(:,1)) + u2.*(P1(:,1)-P3(:,1)) + u3.*(P2(:,1)-P1(:,1))) ./ twoA;
dvdx = (v1.*(P2(:,2)-P3(:,2)) + v2.*(P3(:,2)-P1(:,2)) + v3.*(P1(:,2)-P2(:,2))) ./ twoA;
dvdy = (v1.*(P3(:,1)-P2(:,1)) + v2.*(P1(:,1)-P3(:,1)) + v3.*(P2(:,1)-P1(:,1))) ./ twoA;

% the interpolant is linear, so its value at the centroid IS the mean of the
% three vertices -- exact, and at the same point the gradient applies to
uvc = [(u1 + u2 + u3) / 3, (v1 + v2 + v3) / 3];

tri = struct('conn', conn, ...
    'grad', [dudx, dudy, dvdx, dvdy], ...
    'divergence', dudx + dvdy, ...
    'area', abs(twoA) / 2, ...
    'cxy', [cx, cy], ...
    'uvc', uvc);

info = struct('xy', xy, ...
    'n_vector', numel(uu), ...
    'n_tri_kept', nnz(keep), ...
    'max_edge_measured', max_edge_measured, ...
    'dropped', dropped);
end

% ---------------------------------------------------------------------------
function tri = empty_tri()
%EMPTY_TRI  The same fields the normal path returns, with no rows.
tri = struct('conn', zeros(0,3), ...
    'grad', zeros(0,4), ...
    'divergence', zeros(0,1), ...
    'area', zeros(0,1), ...
    'cxy', zeros(0,2), ...
    'uvc', zeros(0,2));
end

% ---------------------------------------------------------------------------
function dropped = empty_dropped()
%EMPTY_DROPPED  The same counter fields, all zero.
dropped = struct('long_edge', 0, 'sliver', 0, 'masked', 0, 'kept', 0, 'total', 0);
end
