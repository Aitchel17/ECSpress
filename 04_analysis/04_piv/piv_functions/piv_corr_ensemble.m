function corr_planes = piv_corr_ensemble(imgstack, opt)
%PIV_CORR_ENSEMBLE  Ensemble PIV correlation over an in-memory image stack.
%   Sums the cross-correlation planes of every frame pair, then runs multipass
%   window deformation, giving one time-averaged displacement field on the
%   interrogation grid.
%
%   CORRELATION ONLY. No gate, no filter, no rejection of any kind, and one
%   output: what the correlation actually found. Judging a vector needs an order
%   (the field-wide translation has to come out before anything that keys off
%   vector length, PIVlab's correlation filter has to come before its local
%   median), and an engine that gates internally hides that order from whoever
%   owns it. So it stays out here and the caller runs it -- see
%   analysis_pivensemble.gate(), which is one method per gate for the same
%   reason.
%
%   The sparse grid is the result; nothing is stamped into a dense image either,
%   because a vector belongs to a WINDOW and turning that into a pixel is a
%   presentation choice the caller makes.
%
% IN   imgstack       H x W x N double, values in [0,1]. Pairs are (1,2) (3,4)
%                     ..., non-overlapping, so an interleaved stack correlates
%                     each from-frame against its own to-frame
%      window_sizes   Px2 [window step] per pass, or Px1 with step = window/2,
%                     default [64 32; 32 16; 16 8]
%      roirect        [x y w h] ROI in pixel coords, default full frame
%      subpixfinder   1 = 3-point Gaussian (default), 2 = 2D Gaussian
%      mask_auto      1 = suppress the auto-correlation peak, default 0
%      imdeform       deformation interpolant, default '*spline'
%      repeat         1 = repeated/multishift correlation on the last pass, default 0
%      do_pad         1 = zero-padded linear correlation on the last pass, default 1
%      use_gpu        1 = FFT/deform on the GPU, default false
%      exclmask       H x W, true = EXCLUDE that pixel from PIV, default none
%      save_corrmaps  1 = also carry the per-pair planes, default false. They are
%                     IA x IA x ny x nx x nPair single and dwarf everything else,
%                     which is the only reason they are optional
% OUT  corr_planes    one struct, always
%        .xtable/.ytable  ny x nx image coords of each grid vector
%        .utable/.vtable  ny x nx displacement, px, RAW
%        .corr            ny x nx corr2 of each window pair's final cuts
%        .typevector      ny x nx, 1 = valid, 0 = fully masked by exclmask
%        .imsize          [H W] of the stack, so piv_stamp can place the grid
%                         back on the image without being told twice
%        .maps            IA x IA x ny x nx x nPair single, pair p's plane at
%                         (r,c); sum(...,5) = the ensembled plane the peak was
%                         fitted to. [] unless save_corrmaps
%        .corr2           ny x nx x nPair per-pair corr2. [] unless save_corrmaps
%        .pair_frames     nPair x 2 slice indices forming each pair. [] unless
%                         save_corrmaps
%        .pass            which pass the planes came from (always the final one)
%        .per_pass        1 x nPass struct, one row per pass: maps, window, step,
%                         xtable/ytable, utable/vtable (that pass alone) and
%                         utable_total (running). [] unless save_passes

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
    opt.exclmask     {mustBeNumericOrLogical} = []   % H x W, true = EXCLUDE pixel from PIV
    opt.save_corrmaps logical = false  % keep each frame pair's final-pass planes
    opt.save_passes   logical = false  % keep EVERY pass's planes as well. Only pass
                                       % 1 sees an undeformed image2; the rest read
                                       % a residual off an already-shifted one.
                                       % ~20 MB per result, off by default
end

% 0. Setup
% 0.1 Stack size; a pair is the smallest thing that can be correlated
[H, W, N] = size(imgstack);
if N < 2
    error('piv_corr_ensemble:tooFewFrames', 'imgstack must have at least 2 frames.');
end

% 0.2 ROI, defaulting to the full frame
roirect = opt.roirect;
if isempty(roirect)
    roirect = [1, 1, W-1, H-1];
end

% 0.3 Split window_sizes into pivensemble's per-pass arguments
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

% 0.4 One mask per image pair, in PIVlab's polarity: nonzero = pixel EXCLUDED
%     from correlation, zero = analyzed
n_pairs = floor(N / 2);  % pivensemble uses non-overlapping pairs (1:2:N)
if isempty(opt.exclmask)
    frame_mask = zeros(H, W);                % nothing excluded
else
    if ~isequal(size(opt.exclmask), [H, W])
        error('piv_corr_ensemble:maskSize', ...
            'exclmask must be %dx%d (H x W) to match imgstack.', H, W);
    end
    frame_mask = double(opt.exclmask ~= 0);  % logical/numeric -> PIVlab mask
end
converted_mask = repmat({frame_mask}, n_pairs, 1);

% 1. Correlate the stack
[xtable, ytable, utable, vtable, typevector, correlation_map, raw, raw_pass] = pivensemble( ...
    imgstack, roirect, converted_mask, interrogationarea, step_size, opt.subpixfinder, ...
    passes, int2, int3, int4, opt.mask_auto, opt.imdeform, opt.repeat, opt.do_pad, ...
    opt.save_corrmaps, opt.save_passes);

% 2. Bring the tables back from the GPU
if opt.use_gpu
    xtable = gather(xtable);
    ytable = gather(ytable);
    utable = gather(utable);
    vtable = gather(vtable);
    typevector = gather(typevector);
    correlation_map = gather(correlation_map);
end

% 3. The grid, always. Every field is on the same ny x nx, so a caller can gate
%    on it and still line the result up with the planes afterwards
corr_planes = struct( ...
    'xtable',      double(xtable), ...
    'ytable',      double(ytable), ...
    'utable',      double(utable), ...
    'vtable',      double(vtable), ...
    'corr',        double(correlation_map), ...
    'typevector',  typevector, ...
    'imsize',      [H, W], ...
    'maps',        [], ...
    'corr2',       [], ...
    'pair_frames', [], ...
    'pass',        [], ...
    'per_pass',    []);

% 4. Fold the captured planes from window-linear order onto that grid
if opt.save_corrmaps && ~isempty(raw)
    [ny, nx] = size(xtable);
    [s1, s2, ~, np] = size(raw.maps);
    % 4.1 maps(:,:,r,c,p) = pair p at (ytable(r,c), xtable(r,c)); pivensemble's
    %     window index runs x-fastest, hence the nx/ny order
    corr_planes.maps        = permute(reshape(raw.maps, s1, s2, nx, ny, np), [1 2 4 3 5]);
    corr_planes.corr2       = raw.corr2;
    corr_planes.pair_frames = raw.pair_frames;
    corr_planes.pass        = raw.pass;
end

% 5. Every pass on its OWN grid, folded the same way. Not regridded, not cropped
%    and not summed: pass 1 is measured against an unshifted image2 and the later
%    ones against a deformed one, so how to combine them is the caller's question
%    see FINDINGS.md
%    See the note on pass 1 in pivensemble
if opt.save_passes && ~isempty(raw_pass)
    for k = numel(raw_pass) : -1 : 1     % reverse, so the array is sized once
        pk = raw_pass{k};
        [ny_k, nx_k] = size(pk.xtable);
        [q1, q2, ~, np_k] = size(pk.maps);
        pk.maps  = permute(reshape(pk.maps, q1, q2, nx_k, ny_k, np_k), [1 2 4 3 5]);
        pk.reach = (q1 - 1)/2;            % px, the largest displacement it can hold
        per_pass(k) = pk;
    end
    corr_planes.per_pass = per_pass;
end

end
