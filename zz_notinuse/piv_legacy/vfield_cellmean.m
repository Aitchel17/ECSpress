function [uv_sparse, info] = vfield_cellmean(uv, cellmap, opt)
%VFIELD_CELLMEAN  One mean vector per partition cell, as a drawable sparse field.
%   Averages the surviving vectors inside each cell of cellmap and returns a map
%   holding that mean at the cell's centroid and NaN everywhere else. showpiv's
%   plot_quiver reads exactly that shape, so no new plotting method is needed and
%   the arrows land on the same axes as the dense field.
%
% IN   uv        H x W x 2 float  dense field, NaN off the vector grid
%      cellmap   H x W int        0 = outside every cell, else the cell label.
%                                 vfield_strain(method='tri') returns one as
%                                 info.sectormap; sector_polar returns the wedges
%      min_n     vectors a cell needs before it gets an arrow, default 3
% OUT  uv_sparse H x W x 2 float  the means at the centroids, NaN elsewhere
%      info      .n_cell     cells that produced an arrow
%                .n_vector   vectors that went into them
%                .per_cell   1 x n_labels int, vectors averaged per cell
%                .mag        1 x n_cell float px, |mean| of the drawn arrows
%                .mag_dense  1 x n float px, |v| of the individual vectors
%
%   why  averaging BEFORE drawing is not the same as drawing fewer arrows. A
%        decimated field shows one arbitrary member of each neighbourhood; this
%        shows what the neighbourhood agreed on, and the disagreement is what
%        gets suppressed
%   caution  the arrow is placed at the centroid of the VECTORS in the cell, not
%        of the cell. With a cell half-covered by the exclusion mask the two are
%        far apart, and the vector centroid is the honest one
%   err  a cell straddling the wall averages inward and outward motion to nearly
%        zero. That is real cancellation, not a bug -- but it is why min_n alone
%        is not enough and the polar wedges exist

arguments
    uv        {mustBeNumeric, mustBeNonempty}
    cellmap   {mustBeNumeric, mustBeNonempty}
    opt.min_n (1,1) double {mustBePositive, mustBeInteger} = 3
end

u = uv(:,:,1);   v = uv(:,:,2);
live = ~isnan(u) & ~isnan(v);
[H, W, ~] = size(uv);
if ~isequal(size(cellmap), [H W])
    error('vfield_cellmean:size', ...
        'cellmap is %s, uv is %d x %d.', mat2str(size(cellmap)), H, W);
end

n_labels  = max(cellmap(:));
uv_sparse = nan(H, W, 2);
per_cell  = zeros(1, n_labels);
mag       = [];
n_vector  = 0;
n_cell    = 0;

[yy, xx] = ndgrid(1:H, 1:W);
for L = 1:n_labels
    sel = live & cellmap == L;
    per_cell(L) = nnz(sel);
    if per_cell(L) < opt.min_n
        continue
    end
    mu = mean(u(sel));   mv = mean(v(sel));
    cy = round(mean(yy(sel)));   cx = round(mean(xx(sel)));
    uv_sparse(cy, cx, 1) = mu;
    uv_sparse(cy, cx, 2) = mv;
    mag(end+1) = hypot(mu, mv);  %#ok<AGROW>
    n_vector = n_vector + per_cell(L);
    n_cell   = n_cell + 1;
end

info = struct( ...
    'n_cell',    n_cell, ...           % int          arrows drawn
    'n_vector',  n_vector, ...         % int          vectors behind them
    'per_cell',  per_cell, ...         % 1 x L int    vectors per cell
    'mag',       mag, ...              % 1 x n float px  |mean vector|
    'mag_dense', hypot(u(live), v(live))');   % 1 x n float px  |individual|
end
