function [u_out, v_out] = piv_validate(utable, vtable, correlation_map, opt)
%PIV_VALIDATE  Reject unreliable PIV vectors on the sparse interrogation grid.
%   Call before stamping the vectors into a dense map, so PIVlab's local-median
%   filter sees real neighbours. Stages 1 and 3 are PIVlab originals and must stay
%   in PIVlab's own order: correlation filter first, then PIVlab_postproc.
%
% IN   utable           sparse u grid from pivensemble
%      vtable           sparse v grid from pivensemble
%      correlation_map  same-size per-window corr2 map from pivensemble
%      corr_thr         corr2 reject threshold; low corr2 = out-of-plane (z)
%                       decorrelation of the window pair, default 0.5
%      peak_width       same-size map from piv_peak_width, [] = skip stage 2
%      min_snr          reject where |disp| < min_snr * peak_width, default 1, so
%                       a displacement must beat its own localisation
%      valid_vel        [umin umax vmin vmax] px limits, [] to skip, default []
%      neigh_thresh     Westerweel-Scarano universal-median threshold, default 2
% OUT  u_out, v_out     validated vectors, rejected ones set to NaN

arguments
    utable          double
    vtable          double
    correlation_map double
    opt.corr_thr     double = 0.5
    opt.peak_width   double = []
    opt.min_snr      double = 1
    opt.valid_vel    double = []
    opt.neigh_thresh double = 2
end

% 1. Correlation-coefficient rejection; arg order is (u, v, threshold, corr2 values)
[u, v] = postproc.PIVlab_correlation_filter(utable, vtable, opt.corr_thr, correlation_map);

% 2. Peak-localisation rejection, which corr2 cannot do: a blurred window still
%    correlates well, but its peak is too wide to place the displacement
if ~isempty(opt.peak_width)
    if ~isequal(size(opt.peak_width), size(u))
        error('piv_validate:widthSize', ...
            'peak_width must be %s to match the vector grid.', mat2str(size(u)));
    end
    weak = hypot(u, v) < opt.min_snr * opt.peak_width;
    u(weak) = NaN;
    v(weak) = NaN;
end

% 3. Spatial outlier rejection (calu = calv = 1, no stdev check, local median on)
[u_out, v_out] = postproc.PIVlab_postproc(u, v, 1, 1, opt.valid_vel, false, [], true, opt.neigh_thresh);

end
