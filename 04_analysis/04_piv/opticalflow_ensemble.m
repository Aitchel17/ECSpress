function [uv_map, corr_map_out] = opticalflow_ensemble(imgstack, opt)
%OPTICALFLOW_ENSEMBLE Run PIVlab ensemble correlation on an in-memory image stack.
%   Accumulates cross-correlation matrices across all consecutive frame pairs
%   and performs multipass window deformation to yield a single, robust,
%   time-averaged velocity field.
%
%   USAGE
%     [uv, corr] = opticalflow_ensemble(imgstack)
%     [uv, corr] = opticalflow_ensemble(imgstack, window_sizes=[64 32; 32 16])
%
%   INPUTS
%     imgstack       : double H x W x N array, values in [0,1]
%     roirect        : [x y w h] ROI in pixel coords (default: full frame)
%     window_sizes   : P x 2 matrix of [window_size, step] per pass, or
%                      P x 1 vector if step = window/2.
%                      (e.g. [64 32; 32 16; 16 8] for 3 passes)
%     subpixfinder   : 1 = 3-point Gaussian (default), 2 = 2D Gaussian
%     mask_auto      : 1 = suppress auto-correlation peak (default: 0)
%     imdeform       : Image deformation interpolation (default: '*spline')
%     repeat         : 1 = repeated/multishift correlation on last pass (default: 0)
%     do_pad         : 1 = linear zero-padded correlation on last pass (default: 1)
%     use_gpu        : 1 = perform FFT/deform operations on GPU (default: 0)
%
%   OUTPUTS
%     uv_map         : H x W x 2 dense velocity map
%                      (:,:,1) = u (horizontal), (:,:,2) = v (vertical)
%                      NaN outside the valid vector grid
%     corr_map_out   : H x W correlation strength map

arguments
    imgstack         {mustBeNumeric, mustBeNonempty}
    opt.window_sizes double  = [64 32; 32 16; 16 8]
    opt.roirect      double  = []
    opt.subpixfinder double  = 1
    opt.mask_auto    double  = 0
    opt.imdeform     char    = '*spline'
    opt.repeat       double  = 0
    opt.do_pad       double  = 1
    opt.use_gpu      logical = false
end

[H, W, N] = size(imgstack);
if N < 2
    error('opticalflow_ensemble: imgstack must have at least 2 frames.');
end

% ROI
roirect = opt.roirect;
if isempty(roirect)
    roirect = [1, 1, W-1, H-1];
end


% Decompose window_sizes into individual parameters for pivensemble
ws = opt.window_sizes;
passes = size(ws, 1);
if size(ws, 2) >= 2
    step_size = ws(1, 2);
else
    step_size = floor(ws(1, 1) / 2);
end
interrogationarea = ws(1, 1);
int2 = ws(min(2, end), 1);
int3 = ws(min(3, end), 1);
int4 = ws(min(4, end), 1);

% Build converted_mask cell array (one mask per image pair, all zeros = no mask)
n_pairs = floor(N / 2);  % pivensemble uses non-overlapping pairs (1:2:N)
roi_h = roirect(4) + 1;
roi_w = roirect(3) + 1;
if isempty(roirect) || all(roirect == [0 0 0 0])
    mask_size = [H, W];
else
    mask_size = [roi_h, roi_w];
end
converted_mask = repmat({zeros(H, W)}, n_pairs, 1);

% Call pivensemble (original PIVlab-style, for debugging)
[xtable, ytable, utable, vtable, typevector, correlation_map] = pivensemble( ...
    imgstack, roirect, converted_mask, interrogationarea, step_size, opt.subpixfinder, ...
    passes, int2, int3, int4, opt.mask_auto, opt.imdeform, opt.repeat, opt.do_pad);

if opt.use_gpu
    xtable = gather(xtable);
    ytable = gather(ytable);
    utable = gather(utable);
    vtable = gather(vtable);
    typevector = gather(typevector);
    correlation_map = gather(correlation_map);
end


% Place sparse PIV vectors into dense H x W grid (no interpolation).
% Vectors are stamped at their rounded grid coordinates; all other
% pixels remain NaN.
uv_map       = NaN(H, W, 2);
corr_map_out = NaN(H, W);

valid = (typevector == 1) & ~isnan(utable) & ~isnan(vtable);
rows = round(double(ytable(valid)));
cols = round(double(xtable(valid)));

% Clamp to image bounds
rows = max(1, min(H, rows));
cols = max(1, min(W, cols));

idx = sub2ind([H, W], rows, cols);
uv_map(idx)         = double(utable(valid));     % u into (:,:,1)
uv_map(idx + H*W)   = double(vtable(valid));     % v into (:,:,2)
corr_map_out(idx)    = double(correlation_map(valid));

end
