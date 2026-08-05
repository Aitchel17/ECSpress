function step = piv_grid_step(flow)
%PIV_GRID_STEP  Spacing (px) of the sparse vector grid in a PIV field.
%   corr_ensemble returns NaN everywhere except the interrogation-window centres,
%   so the spacing is the gap between the rows and columns that hold vectors. The
%   cell-map builders size their cells from it.
%
% IN   flow  HxWx2 (or HxWxNx2) field, NaN off the vector grid
% OUT  step  grid spacing in px, NaN when fewer than two rows or columns hold vectors

arguments
    flow {mustBeNumeric, mustBeNonempty}
end

% 0. Where the field actually has vectors
if ndims(flow) == 4
    valid = any(~isnan(flow(:,:,:,1)), 3) & any(~isnan(flow(:,:,:,2)), 3);
else
    f     = squeeze(flow);
    valid = ~isnan(f(:,:,1)) & ~isnan(f(:,:,2));
end

% 1. Median gap between occupied rows and columns
vc = find(any(valid, 1));
vr = find(any(valid, 2));
if numel(vc) < 2 || numel(vr) < 2
    step = NaN;
    return
end
step = max(1, round(median([median(diff(vc)), median(diff(vr))])));
end
