function grad = vfield_trigradient(xy, conn, uv)
%VFIELD_TRIGRADIENT  The exact gradient of the linear interpolant on each triangle.
%   A linear function on a triangle is fixed by its three vertex values, so the
%   gradient is CONSTANT over the cell and exact -- no fit, no stencil, no window.
%   Four gradient components come out of six numbers, and there is no parameter
%   anywhere in this file that could set the answer.
%
%   The half of the triangle route that depends on what the vectors SAY.
%   vfield_applydelaunay makes the mesh; this differentiates on it.
%
% IN   xy    nPoint x 2 float px  vertex positions, from vfield_applydelaunay
%      conn  nTri x 3 int         vertex indices into xy
%      uv    nPoint x 2 float px  displacement at each vertex, xy's own order
% OUT  grad  nTri x 4 float per px  [du/dx du/dy dv/dx dv/dy]
%
% UNIT  displacement in whatever uv was; the gradient is that per PIXEL, so the
%       trace is dimensionless and does not need a scale
%   why  the four columns are the whole tensor and NOTHING is derived from them
%        here. Divergence is the trace, the two stretches are contractions against
%        an (n, t) pair, rotation is (du/dy - dv/dx)/2 and the shear is the
%        symmetric off-diagonal -- four projections of one object. Computing one of
%        them in this file would make it look like a second measurement, and the
%        caller that picks the direction is the one that should take all four
%   why  the field's VALUE at the centroid is not returned. It is the mean of the
%        three vertices -- exact, for the same reason the gradient is -- but it is
%        the function rather than its derivative, and one caller wants it. Doing it
%        here would make this file answer two questions
%   note the signed cross product 2A appears in numerator and denominator, so the
%        triangle's orientation cancels and no winding convention is needed

arguments
    xy   (:,2) double
    conn (:,3) double
    uv   (:,2) double
end
if size(uv, 1) ~= size(xy, 1)
    error('vfield_trigradient:vertexCount', ...
        'xy has %d vertices and uv has %d.', size(xy, 1), size(uv, 1));
end
if ~isempty(conn) && max(conn(:)) > size(xy, 1)
    error('vfield_trigradient:connRange', ...
        'conn indexes vertex %d; xy has %d.', max(conn(:)), size(xy, 1));
end
if isempty(conn)
    grad = zeros(0,4);
    return
end

P1 = xy(conn(:,1), :);
P2 = xy(conn(:,2), :);
P3 = xy(conn(:,3), :);
twoA = (P2(:,1)-P1(:,1)).*(P3(:,2)-P1(:,2)) - (P3(:,1)-P1(:,1)).*(P2(:,2)-P1(:,2));

u1 = uv(conn(:,1), 1);
u2 = uv(conn(:,2), 1);
u3 = uv(conn(:,3), 1);
v1 = uv(conn(:,1), 2);
v2 = uv(conn(:,2), 2);
v3 = uv(conn(:,3), 2);

dudx = (u1.*(P2(:,2)-P3(:,2)) + u2.*(P3(:,2)-P1(:,2)) + u3.*(P1(:,2)-P2(:,2))) ./ twoA;
dudy = (u1.*(P3(:,1)-P2(:,1)) + u2.*(P1(:,1)-P3(:,1)) + u3.*(P2(:,1)-P1(:,1))) ./ twoA;
dvdx = (v1.*(P2(:,2)-P3(:,2)) + v2.*(P3(:,2)-P1(:,2)) + v3.*(P1(:,2)-P2(:,2))) ./ twoA;
dvdy = (v1.*(P3(:,1)-P2(:,1)) + v2.*(P1(:,1)-P3(:,1)) + v3.*(P2(:,1)-P1(:,1))) ./ twoA;

grad = [dudx, dudy, dvdx, dvdy];
end
