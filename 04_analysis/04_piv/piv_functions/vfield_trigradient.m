function [grad, shapegrad, nodaldisp] = vfield_trigradient(xy, conn, uv)
%VFIELD_TRIGRADIENT  The exact gradient of the linear interpolant on each triangle.
%   A linear function on a triangle is fixed by its three vertex values, so the
%   gradient is CONSTANT over the cell and exact -- no fit, no stencil, no window.
%   Four gradient components come out of six numbers, and there is no parameter
%   anywhere in this file that could set the answer.
%
%   The half of the triangle route that depends on what the vectors SAY.
%   vfield_applydelaunay makes the mesh; this differentiates on it.
%
%   THIS IS THE CST, the constant strain triangle of finite elements. In that
%   vocabulary the whole file is one product:
%
%     J = U * B          J  2 x 2  [du/dx du/dy; dv/dx dv/dy]
%                        U  2 x 3  the three nodal displacements
%                        B  3 x 2  row i = grad N_i, the shape function gradient
%
%     grad N_i = (y_j - y_k,  x_k - x_j) / 2A     j, k = the other two nodes
%
%   READ IT BY INDEX and the point of the product shows:
%
%     J(c,d) = sum over k of  U(c,k) * B(k,d)
%
%       c   which displacement component    u or v         comes from U
%       d   which direction differentiated  d/dx or d/dy   comes from B
%       k   which node                      1, 2, 3        summed away
%
%   J's row says WHAT MOVED and its column says WHICH WAY IT WAS DIFFERENTIATED;
%   the node index that joins them is contracted out. So B is this triangle's
%   discrete differentiation operator -- three nodal scalars in, two derivatives
%   out -- and it is the SAME operator for u and for v, which is why one B serves
%   both. Data and geometry never mix in this form, and that separation is what
%   lets B alone answer the noise question and B'*Sigma*B be the whole of error
%   propagation.
%
%   All three come back, and all three carry the triangle as the FIRST dimension
%   the way grad and nxy already do -- so U and B are indexed (tri, node, xy),
%   not the 2 x 3 and 3 x 2 pages written above. The formula is the maths; the
%   arrays are the maths with the page index moved to the front.
%
%   B is geometry only, so it is the same for any displacement on the same mesh.
%   It is NOT cached: applydelaunay triangulates wherever vectors survived, and
%   piv_validate NaNs a different set every event, so the mesh differs per event.
%   grad N_i points from the opposite edge toward node i and its length is 1/h_i,
%   h_i being that node's height -- which makes 1/h the noise gain of this whole
%   file, and is why vfield_trifilter judges cell shape before measure() runs.
%   A textbook "B matrix" is a DIFFERENT 3 x 6 object: it maps a 6 x 1 nodal
%   vector straight to Voigt strain, so it has already symmetrised and thrown the
%   rotation away, and its third row is the engineering shear gamma = 2*E12.
%
% IN   xy    nPoint x 2 float px  vertex positions, from vfield_applydelaunay
%      conn  nTri x 3 int         vertex indices into xy
%      uv    nPoint x 2 float px  displacement at each vertex, xy's own order
% OUT  grad       nTri x 4 float per px    [du/dx du/dy dv/dx dv/dy], = J
%      shapegrad  SHAPE function GRADient. nTri x 3 x 2 float per px
%                 = B. (:,i,:) is grad N_i. GEOMETRY
%                 ONLY, no displacement enters it. Its row norms are 1/h_i, so
%                 max over i is this triangle's noise gain, and Sigma propagates
%                 into the gradient as B'*Sigma*B. Nobody reads it yet; it is
%                 what PV24 and PV25 need and it costs nothing to hand back
%      nodaldisp  NODAL DISPlacement, the displacement at each of the three
%                 nodes. nTri x 3 x 2 float px
%                 = U. (:,i,:) is node i's [u v].
%                 This one IS uv(conn) reindexed and carries no new information
%                 -- it is returned so the three names in the formula above are
%                 three things a caller can actually hold
%
% UNIT  displacement in whatever uv was; the gradient is that per PIXEL, so the
%       trace is dimensionless and does not need a scale
%   why  the four columns are the whole tensor and NOTHING is derived from them
%        here. Divergence is the trace, the two stretches are contractions against
%        an (n, t) pair, rotation is (dv/dx - du/dy)/2 -- half the curl, + counter-
%        clockwise -- and the shear is the symmetric off-diagonal resolved onto the
%        same (n, t) pair: five projections of one object. Computing one of
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
    grad      = zeros(0, 4);
    shapegrad = zeros(0, 3, 2);
    nodaldisp = zeros(0, 3, 2);
    return
end

P1 = xy(conn(:,1), :);                % node 1 position
P2 = xy(conn(:,2), :);                % node 2
P3 = xy(conn(:,3), :);                % node 3
% determinant of [x y 1; x y 1; x y 1], the system the three nodal values set.
% It equals twice the SIGNED area, and every grad N_i below divides by it
twoA = (P2(:,1)-P1(:,1)).*(P3(:,2)-P1(:,2)) - (P3(:,1)-P1(:,1)).*(P2(:,2)-P1(:,2));

% B. The six numbers below are M's cofactors, written out. Replacing M's column 1
% by the nodal values and expanding gives node i the coefficient C_i1, and column
% 2 gives C_i2; the array's two subscripts are exactly those two indices --
% (:, i, c) is node i, column c replaced. Derivation in GRADIENT_TENSOR_EQUATION
% sections 3 (built from area coordinates) and 10.2 (the cofactor reading).
% Building B here rather than inlining the brackets is what stops the six
% coefficients being written twice, once for u and once for v
shapegrad = zeros(size(conn, 1), 3, 2);
shapegrad(:, 1, 1) = (P2(:,2) - P3(:,2)) ./ twoA;    % grad N_1, d/dx
shapegrad(:, 1, 2) = (P3(:,1) - P2(:,1)) ./ twoA;    % grad N_1, d/dy
shapegrad(:, 2, 1) = (P3(:,2) - P1(:,2)) ./ twoA;    % grad N_2, d/dx
shapegrad(:, 2, 2) = (P1(:,1) - P3(:,1)) ./ twoA;    % grad N_2, d/dy
shapegrad(:, 3, 1) = (P1(:,2) - P2(:,2)) ./ twoA;    % grad N_3, d/dx
shapegrad(:, 3, 2) = (P2(:,1) - P1(:,1)) ./ twoA;    % grad N_3, d/dy

% U. uv gathered onto the triangles, which is all this is
nodaldisp = zeros(size(conn, 1), 3, 2);
nodaldisp(:, 1, :) = uv(conn(:,1), :);               % node 1, [u v]
nodaldisp(:, 2, :) = uv(conn(:,2), :);               % node 2
nodaldisp(:, 3, :) = uv(conn(:,3), :);               % node 3

% J = U * B, summed over the three nodes
dudx = nodaldisp(:,1,1).*shapegrad(:,1,1) + nodaldisp(:,2,1).*shapegrad(:,2,1) + nodaldisp(:,3,1).*shapegrad(:,3,1);
dudy = nodaldisp(:,1,1).*shapegrad(:,1,2) + nodaldisp(:,2,1).*shapegrad(:,2,2) + nodaldisp(:,3,1).*shapegrad(:,3,2);
dvdx = nodaldisp(:,1,2).*shapegrad(:,1,1) + nodaldisp(:,2,2).*shapegrad(:,2,1) + nodaldisp(:,3,2).*shapegrad(:,3,1);
dvdy = nodaldisp(:,1,2).*shapegrad(:,1,2) + nodaldisp(:,2,2).*shapegrad(:,2,2) + nodaldisp(:,3,2).*shapegrad(:,3,2);

grad = [dudx, dudy, dvdx, dvdy];      % J, flattened row-major
end
