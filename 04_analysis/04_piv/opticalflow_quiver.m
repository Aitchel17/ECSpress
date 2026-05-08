function opticalflow_quiver(flow, opt)
%OPTICALFLOW_QUIVER  PIVlab-style quiver plot with optional temporal sum.
%
%   USAGE
%     opticalflow_quiver(flow(:,:,k,:))          % single pair
%     opticalflow_quiver(flow, time_range=[s e]) % sum pairs s..e
%
%   INPUTS
%     flow          : H x W x N_pairs x 2  (full 4D from opticalflow_wlet)
%                     or  H x W x 1 x 2 / H x W x 2  (single pair slice)
%     time_range    : [start end] pair indices to sum over (default: all pairs)
%     scale         : vector scale factor (default: 1)
%     block_size    : spatial averaging block (default: 1 = no averaging)
%     uniform       : normalize to unit length (default: false)
%     color_valid   : RGB quiver color (default: [0 1 0] green)
%     linewidth     : arrow line width (default: 0.5)

arguments
    flow             {mustBeNumeric}
    opt.time_range   double  = []   % [start end], default: all pairs
    opt.scale        double  = 1
    opt.block_size   double  = 1
    opt.uniform      logical = false
    opt.color_valid  double  = [0 1 0]
    opt.linewidth    double  = 0.5
end

%% Temporal sum over selected range
%   flow is H x W x N x 2; sum over dim 3 to accumulate displacement
ndims_flow = ndims(flow);
if ndims_flow == 4
    N = size(flow, 3);
    if isempty(opt.time_range)
        t_idx = 1:N;
    else
        t_idx = round(opt.time_range(1)) : round(opt.time_range(2));
        t_idx = t_idx(t_idx >= 1 & t_idx <= N);
    end
    u = sum(flow(:,:,t_idx,1), 3);   % H x W  (temporal sum)
    v = sum(flow(:,:,t_idx,2), 3);
else
    % Single-frame input: H x W x 2  or  H x W x 1 x 2
    f = squeeze(flow);
    u = f(:,:,1);
    v = f(:,:,2);
end

[H, W] = size(u);
[x, y] = meshgrid(1:W, 1:H);

%% Optional: normalize to unit length (direction only, ignore magnitude)
if opt.uniform
    mag = sqrt(u.^2 + v.^2);
    mag(mag == 0) = 1;
    u = u ./ mag;
    v = v ./ mag;
end

%% Spatial block averaging
%   Divides the field into block_size x block_size tiles and averages u,v
%   within each tile. The displayed arrow is placed at the tile center.
s = round(opt.block_size);
if s > 1
    % Trim to integer multiple of block size
    Hr = floor(H/s)*s;
    Wr = floor(W/s)*s;
    u  = u(1:Hr, 1:Wr);
    v  = v(1:Hr, 1:Wr);

    % Block-average rows then columns
    %   reshape to (s, nrow, Wr) → mean over dim1 → (nrow, Wr)
    u = squeeze(mean(reshape(u, s, Hr/s, Wr), 1));
    v = squeeze(mean(reshape(v, s, Hr/s, Wr), 1));
    %   reshape to (s, ncol, nrow) → mean over dim1 → (ncol, nrow) → transpose
    u = squeeze(mean(reshape(u', s, Wr/s, Hr/s), 1))';
    v = squeeze(mean(reshape(v', s, Wr/s, Hr/s), 1))';

    % Coordinates: centre of each block
    xc = (0.5:Wr/s) * s;
    yc = (0.5:Hr/s) * s;
    [x, y] = meshgrid(xc, yc);
end

%% Draw quiver
hold on;
quiver(x, y, u*opt.scale, v*opt.scale, ...
    'Color', opt.color_valid, ...
    'AutoScale', 'off', ...
    'LineWidth', opt.linewidth);
hold off;
axis image;

end