function [keep, dropped] = vfield_trifilter(xy, conn, opt)
%VFIELD_TRIFILTER  Which Delaunay cells are fit to differentiate on. A verdict.
%   Three tests, each catching a different failure, and none of them mutates the
%   mesh -- the caller applies the mask in one place. Same shape as
%   vfield_polar.gate_wedge for the same reason: a filter that edited its input
%   would leave no way to see what it took.
%
% IN   xy         nPoint x 2 float px  vertex positions
%      conn       nTri x 3 int         vertex indices into xy
%      max_edge   float px    longest edge kept. [] = 3x the median SHORTEST edge,
%                             which follows the grid step instead of fixing a number
%      min_angle  float deg   smallest interior angle kept, default 15
%      exclmask   H x W bool  true = never interpolate. Tested on all three VERTICES
% OUT  keep     nTri x 1 logical
%      dropped  struct  long_edge / sliver / masked / kept / total counts, plus
%               max_edge_measured (float px), which is reported whether or not
%               max_edge was given
%
%   why  LONG EDGE. Delaunay bridges holes: a cell spanning the masked vessel has
%        no measurement under it at all. Neither the mask test nor the angle test
%        removes it, because its vertices are clear of the mask and its shape is
%        perfectly reasonable
%   why  SLIVER. The gradient divides by twice the signed area, so a thin cell
%        amplifies a sub-pixel error in one vertex without limit -- the same
%        vertex wobble that is negligible on a fat cell is unbounded on a sliver.
%        This is the only reason min_angle exists, and it is why
%        vfield_trigradient carries no guard of its own
%   why  MASK ON THE VERTICES, not the centroid. A cell can have its centroid clear
%        while one corner sits on the boundary, and that corner's interrogation
%        window is shared with the masked region, so its vector carries the
%        boundary's bias straight into the gradient
%   note the three are independent verdicts and the counts overlap. dropped.kept is
%        the size of the intersection, not total minus the three

arguments
    xy            (:,2) double
    conn          (:,3) double
    opt.max_edge  double = []
    opt.min_angle (1,1) double = 15
    opt.exclmask  logical = []
end
if isempty(conn)
    keep = false(0,1);
    dropped = struct('long_edge', 0, 'sliver', 0, 'masked', 0, ...
        'kept', 0, 'total', 0, 'max_edge_measured', NaN);
    return
end

P1 = xy(conn(:,1), :);
P2 = xy(conn(:,2), :);
P3 = xy(conn(:,3), :);

% 1. Edge lengths; the default cut is relative, so it follows the grid step
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

% 2. Smallest interior angle, by the law of cosines on the shortest side
ang = @(a,b,c) acosd(min(1, max(-1, (b.^2 + c.^2 - a.^2) ./ (2*b.*c))));
a1 = ang(e1, e3, e2);
a2 = ang(e2, e1, e3);
a3 = ang(e3, e2, e1);
angle_min = min([a1 a2 a3], [], 2);
keep_angle = angle_min >= opt.min_angle;

% 3. The mask, on all three vertices
keep_mask = true(size(keep_edge));
if ~isempty(opt.exclmask)
    vx = [P1(:,1) P2(:,1) P3(:,1)];
    vy = [P1(:,2) P2(:,2) P3(:,2)];
    % the mask is the only thing here that lives in frame coordinates, so it is
    % also the only thing that knows the frame size
    vi = sub2ind(size(opt.exclmask), round(vy), round(vx));
    keep_mask = ~any(opt.exclmask(vi), 2);
end

keep = keep_edge & keep_angle & keep_mask;
dropped = struct('long_edge', nnz(~keep_edge), ...
                 'sliver',    nnz(~keep_angle), ...
                 'masked',    nnz(~keep_mask), ...
                 'kept',      nnz(keep), ...
                 'total',     numel(keep), ...
                 'max_edge_measured', max_edge_measured);
end
