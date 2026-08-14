function place = vfield_placetri(cxy, dist, near_idx, wedgemap, pixel2um, bin_edges_um)
%VFIELD_PLACETRI  Where each triangle sits in the wall's own frame.
%   Radius from the traced boundary, which wedge, which radial bin, and the local
%   outward normal. All of it is GEOMETRY: nothing here reads a displacement, so
%   two events over one mask get identical output and the same partition.
%
% IN   cxy           n x 2 float px    the points to place, one row each
%      dist          H x W float px    bwdist from the core wall, first output
%      near_idx      H x W uint32      bwdist's SECOND output, the nearest wall pixel
%      wedgemap      H x W float       1..n_wedge, NaN inside the core
%      pixel2um      float             microns per pixel
%      bin_edges_um  1 x nB+1 float    radial bin edges, MICRONS from the wall
% OUT  place  struct of COLUMNS, n rows:
%        r      n x 1 float  MICRONS from the wall
%        wedge  n x 1 float  1..n_wedge, NaN where the point has none
%        bin    n x 1 float  1..nB, NaN outside the bin range
%        nxy    n x 2 float  [nx ny] outward unit normal at the wall, NaN where
%                            the point sits ON the boundary and has no direction
%
%   why  the normal comes from bwdist's own nearest-pixel index. The direction from
%        the nearest boundary pixel to here IS the outward normal at the wall --
%        no axis, no equivalent circle, and it costs nothing because bwdist already
%        computed it. The alternative, a direction from the mask's centroid, is an
%        AXIS rather than a normal and the two disagree most where the hoop term is
%        largest. See CLAUDE_LOG.md
%   note a point landing exactly ON the boundary has no direction, so nxy is NaN
%        there. Real rather than hypothetical: the vector grid reaches the traced
%        edge

arguments
    cxy          (:,2) double
    dist         double
    near_idx     {mustBeNumeric}
    wedgemap     double
    pixel2um     (1,1) double {mustBePositive}
    bin_edges_um (1,:) double
end
frame_size = size(dist);
if ~isequal(size(near_idx), frame_size) || ~isequal(size(wedgemap), frame_size)
    error('vfield_placetri:mapSize', ...
        'dist is %s; near_idx %s and wedgemap %s must match.', ...
        mat2str(frame_size), mat2str(size(near_idx)), mat2str(size(wedgemap)));
end

% 0. Which pixel each point falls on. Clamped: a centroid can sit half a pixel
%    outside the frame and the maps have no row there
px = min(max(round(cxy(:,1)), 1), frame_size(2));
py = min(max(round(cxy(:,2)), 1), frame_size(1));
ci = sub2ind(frame_size, py, px);

% 1. Radius, wedge, bin
r_um  = dist(ci) * pixel2um;
wedge = wedgemap(ci);
bin   = discretize(r_um, bin_edges_um);

% 2. The outward wall normal, from the nearest boundary pixel to the point itself
q = double(near_idx(ci));
[qy, qx] = ind2sub(frame_size, q);
dx = cxy(:,1) - qx;
dy = cxy(:,2) - qy;
len = hypot(dx, dy);
len(len == 0) = NaN;

place = struct('r', r_um, 'wedge', wedge, 'bin', bin, ...
    'nxy', [dx ./ len, dy ./ len]);
end
