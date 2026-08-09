function [div_map, div_mean, div_sum] = vfield_divergence_fd(flow, opt)
%VFIELD_DIVERGENCE_FD  Divergence by finite differences on the sparse vector grid.
%   Not a cell partition: it differentiates the vector grid point by point with
%   MATLAB's divergence(), then tiles the result out by nearest neighbour. Use
%   vfield_divergence with a cell map when a per-region plane fit is wanted.
%   Units per-frame (px displacement / px).
%
% IN   flow      HxWx2 field, last dim 1 = u, 2 = v; HxWxNx2 is summed over time
%      mask      HxW logical ROI for div_mean / div_sum, default [] = whole image
%      exclmask  HxW logical, true = EXCLUDE (PIVlab convention), blanked to NaN
% OUT  div_map   HxW divergence, NaN where excluded or unresolved
%      div_mean  mean over mask
%      div_sum   sum over mask

arguments
    flow         {mustBeNumeric, mustBeNonempty}
    opt.mask     logical = []
    opt.exclmask logical = []
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
% 0.2 ROI mask, default the whole image
if isempty(opt.mask); opt.mask = true(H, W); end
if size(opt.mask, 1) ~= H || size(opt.mask, 2) ~= W
    error('vfield_divergence_fd:maskSize', 'mask must be %dx%d.', H, W);
end
% 0.3 Pixel coordinate grid
[x, y] = meshgrid(1:W, 1:H);

% 1. Divergence map
mask_valid = ~isnan(u) & ~isnan(v);
if ~all(mask_valid, 'all') && any(mask_valid, 'all')
    % 1.1 Rows and columns holding the sparse grid
    valid_rows = find(any(mask_valid, 2));
    valid_cols = find(any(mask_valid, 1));
    % 1.2 Cut out the subgrid
    u_sparse = u(valid_rows, valid_cols);
    v_sparse = v(valid_rows, valid_cols);
    x_sparse = x(valid_rows, valid_cols);
    y_sparse = y(valid_rows, valid_cols);
    % 1.3 Differentiate on the subgrid; spacing comes from x_sparse / y_sparse
    div_sparse = divergence(x_sparse, y_sparse, u_sparse, v_sparse);
    % 1.4 Tile out to the full frame, nearest neighbour
    valid_div = ~isnan(div_sparse);
    if any(valid_div, 'all')
        F = scatteredInterpolant(x_sparse(valid_div), y_sparse(valid_div), ...
            div_sparse(valid_div), 'nearest', 'nearest');
        div_map = F(x, y);
    else
        div_map = NaN(H, W);
    end
else
    % 1.5 Dense field: differentiate directly
    div_map = divergence(x, y, u, v);
end

% 2. Blank the excluded pixels (true = exclude)
if ~isempty(opt.exclmask)
    if ~isequal(size(opt.exclmask), [H, W])
        error('vfield_divergence_fd:maskSize', 'exclmask must be %dx%d.', H, W);
    end
    div_map(opt.exclmask) = NaN;
end

% 3. Mean and sum inside the ROI mask
roi_values = div_map(opt.mask);
div_mean = mean(roi_values, 'omitnan');
div_sum  = sum(roi_values, 'omitnan');
end
