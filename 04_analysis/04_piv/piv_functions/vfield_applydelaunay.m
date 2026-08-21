function [mesh, info] = vfield_applydelaunay(xyuv)
%VFIELD_APPLYDELAUNAY  Triangulate the vectors that landed. Every cell, unjudged.
%   The window centres that carried a vector become the vertices, and Delaunay
%   makes the cells. Nothing is filtered here -- vfield_trifilter judges, and the
%   caller applies its verdict. Keeping the two apart is what lets a figure show
%   which cells were thrown away and why.
%
% IN   xyuv  ny x nx x 4 float px  (:,:,1:2) = [x y] window centre,
%                                  (:,:,3:4) = [u v]. Five dims is N consecutive
%                                  intervals: u and v sum over them, the grid
%                                  does not move
% OUT  mesh  struct. xy and uv are one row per VECTOR, the rest one per triangle:
%        xy    nPoint x 2 float px    the window centres that carried a vector
%        uv    nPoint x 2 float px    their displacement, in xy's own order
%        conn  nTri x 3 int           vertex indices into xy
%        cxy   nTri x 2 float px      centroid, where one gradient will apply
%        area  nTri x 1 float px^2    unsigned
%        grid_idx nPoint x 1 int        LINEAR index into the ny x nx grid, xy's
%                                       own order. The way back to everything that
%                                       lives on the grid
%      info  n_vector (int), what the triangulation had to work with
%
% UNIT  px throughout. Nothing here knows about microns
%   why  xy and uv leave together. The gradient indexes uv BY VERTEX, so a point
%        list separated from its values is one reordering away from silently wrong
%   why  grid_idx : logical indexing flattens, so xy alone cannot say WHICH window
%        a vertex came from. Everything that judged that window -- gates(k).mask,
%        the tomasi lambda2, the corr_floor ratio, typevector -- is ny x nx, and
%        without the index none of it can be read back at a triangle. Linear and
%        not [row col] because those arrays take a linear index directly; ind2sub
%        recovers the subscripts when a reader wants them
%   note Delaunay bridges holes, so this output contains cells spanning the masked
%        vessel. That is not a fault to fix here -- it is what vfield_trifilter's
%        long-edge test removes, and it can only be seen because this stage keeps it

arguments
    xyuv {mustBeNumeric, mustBeNonempty}
end

% 0.1 Split the four pages
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
uv = [u(have), v(have)];
grid_idx = find(have);
info = struct('n_vector', size(xy, 1));

% 0.3 Too few to triangulate. Same field set as the normal path, so a caller
%     never has to test which way it returned
if size(xy, 1) < 3
    mesh = struct('xy', xy, 'uv', uv, 'conn', zeros(0,3), ...
        'cxy', zeros(0,2), 'area', zeros(0,1), 'grid_idx', grid_idx);
    return
end

DT = delaunayTriangulation(xy);
conn = DT.ConnectivityList;
P1 = xy(conn(:,1), :);
P2 = xy(conn(:,2), :);
P3 = xy(conn(:,3), :);

% twoA is the signed cross product; the area is its magnitude, and
% vfield_trigradient recovers the sign it needs from the vertices
twoA = (P2(:,1)-P1(:,1)).*(P3(:,2)-P1(:,2)) - (P3(:,1)-P1(:,1)).*(P2(:,2)-P1(:,2));

mesh = struct('xy', xy, ...
    'uv', uv, ...
    'conn', conn, ...
    'cxy', [(P1(:,1) + P2(:,1) + P3(:,1)) / 3, (P1(:,2) + P2(:,2) + P3(:,2)) / 3], ...
    'area', abs(twoA) / 2, ...
    'grid_idx', grid_idx);
end
