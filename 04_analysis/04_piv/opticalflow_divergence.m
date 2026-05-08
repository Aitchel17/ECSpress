function [div_map, div_mean, div_sum] = opticalflow_divergence(flow, opt)
%OPTICALFLOW_DIVERGENCE Calculates divergence from the optical flow field.
%   and optionally extracts the mean/sum divergence from a specific area (ROI mask).
%
%   USAGE
%     div_map = opticalflow_divergence(flow)
%     [div_map, div_mean] = opticalflow_divergence(flow, mask=BW)
%
%   INPUTS
%     flow       : double H x W x 2 array from opticalflow_wlet
%                  dim 3: 1 = u (horizontal), 2 = v (vertical)
%     mask       : logical H x W binary mask defining the extraction area
%                  (default: [] = calculate for the whole image)
%
%   OUTPUTS
%     div_map    : double H x W array of divergence values
%     div_mean   : double scalar of mean divergence within the mask
%     div_sum    : double scalar of summed divergence within the mask

arguments
    flow          {mustBeNumeric, mustBeNonempty}
    opt.mask      logical = []
end

% Squeeze to handle either H x W x 2 or H x W x 1 x 2 inputs
ndims_flow = ndims(flow);
if ndims_flow == 4
    N = size(flow, 3);
    t_idx = 1:N;
    u = sum(flow(:,:,t_idx,1), 3);   % H x W  (temporal sum)
    v = sum(flow(:,:,t_idx,2), 3);
else
    % Single-frame input: H x W x 2  or  H x W x 1 x 2
    f = squeeze(flow);
    u = f(:,:,1);
    v = f(:,:,2);
end


[H, W] = size(u);

% Generate default mask if not provided
if isempty(opt.mask)
    opt.mask = true(H, W);
end

if size(opt.mask, 1) ~= H || size(opt.mask, 2) ~= W
    error('opticalflow_divergence: Mask size must match the spatial dimensions of the flow array (H x W).');
end

% Create coordinate grid (assuming 1-pixel spacing like PIVlab default)
[x, y] = meshgrid(1:W, 1:H);

% Calculate divergence (PIVlab uses MATLAB's built-in function)
div_map = divergence(x, y, u, v);

% Extract parameters from the area (mask)
% Equivalent to PIVlab's get_mean_of_selection and get_integral_of_selection
roi_values = div_map(opt.mask);
div_mean = mean(roi_values, 'omitnan');
div_sum  = sum(roi_values, 'omitnan');

end
