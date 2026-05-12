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


% Call the clean ensemble core (no file I/O, no VideoReader)
[xtable, ytable, utable, vtable, typevector, correlation_map] = piv_ensemble( ...
    imgstack, opt.window_sizes, ...
    opt.subpixfinder, roirect, [], opt.imdeform, opt.repeat, opt.mask_auto, opt.do_pad, opt.use_gpu);

if opt.use_gpu
    xtable = gather(xtable);
    ytable = gather(ytable);
    utable = gather(utable);
    vtable = gather(vtable);
    typevector = gather(typevector);
    correlation_map = gather(correlation_map);
end


% Map sparse outputs to dense H x W grid
uv_map       = NaN(H, W, 2);
corr_map_out = NaN(H, W);

valid = (typevector == 1) & ~isnan(utable) & ~isnan(vtable);
if any(valid(:))
    X = double(xtable(valid));
    Y = double(ytable(valid));
    U = double(utable(valid));
    V = double(vtable(valid));
    C = double(correlation_map(valid));

    [Xq, Yq] = meshgrid(1:W, 1:H);

    % scatteredInterpolant with 'none' extrapolation → NaN outside convex hull
    F_u = scatteredInterpolant(X, Y, U, 'linear', 'none');
    F_v = scatteredInterpolant(X, Y, V, 'linear', 'none');
    F_c = scatteredInterpolant(X, Y, C, 'linear', 'none');

    uv_map(:,:,1) = F_u(Xq, Yq);
    uv_map(:,:,2) = F_v(Xq, Yq);
    corr_map_out  = F_c(Xq, Yq);
end

end
