function player = util_checkoverlay(overlay, opts)
%UTIL_CHECKOVERLAY  Real-time, scrubbable view of a reconstruction_overlay stack.
%
%   Companion to util_checkstack for the 4-channel stack returned by
%   line_fwhm.reconstruction_overlay, which util_checkstack cannot display
%   (it expects [H W 3 T] RGB, whereas the overlay is [H W frames 4]):
%       overlay(:,:,:,1) = ch1 image (boundary pixels blanked)
%       overlay(:,:,:,2) = ch2 image (boundary pixels blanked)
%       overlay(:,:,:,3) = reconstructed BV  boundary
%       overlay(:,:,:,4) = reconstructed PVS boundary
%
%   It composites a grayscale background + red BV + green PVS into an RGB
%   truecolor movie and plays it with implay (play / pause / scrub / loop) for
%   fast real-time viewing.
%
%   player = util_checkoverlay(overlay)                 % base = ch1, 30 fps
%   player = util_checkoverlay(overlay, Name=Value):
%       BaseChannel : image channel used as the grayscale background (1 or 2, default 1)
%       Fps         : playback frame rate (default 30)
%       Clip        : [lo hi] percentiles for background contrast (default [1 99.5])
%       Dilate      : dilate boundary lines by N px for visibility (default 0)
%
%   The RGB movie is uint8 (~8x smaller than the double 4-channel overlay).
%   The same [H W 3 T] RGB also works in util_checkstack if you want the
%   slider-based viewer instead of the player.

arguments
    overlay {mustBeNumeric}
    opts.BaseChannel (1,1) double {mustBeMember(opts.BaseChannel,[1 2])} = 1
    opts.Fps (1,1) double = 30
    opts.Clip (1,2) double = [1 99.5]
    opts.Dilate (1,1) double {mustBeNonnegative} = 0
end

if ndims(overlay) ~= 4 || size(overlay, 4) < 4
    error('util_checkoverlay:badInput', ...
        'overlay must be [H W frames 4] from line_fwhm.reconstruction_overlay.');
end

% --- Grayscale background (uint8) with robust percentile contrast ---
base = double(overlay(:, :, :, opts.BaseChannel));
lo = prctile(base(:), opts.Clip(1));
hi = prctile(base(:), opts.Clip(2));
base = (base - lo) ./ max(hi - lo, eps);
base = uint8(255 * min(max(base, 0), 1));            % [H W T] uint8

% --- Boundary masks ---
bv  = overlay(:, :, :, 3) > 0;
pvs = overlay(:, :, :, 4) > 0;
if opts.Dilate > 0
    se  = strel('disk', round(opts.Dilate));
    bv  = imdilate(bv,  se);
    pvs = imdilate(pvs, se);
end

% --- Compose RGB channels (uint8): BV -> red, PVS -> green ---
R = base; G = base; B = base;
R(bv)  = 255; G(bv)  = 0;   B(bv)  = 0;
R(pvs) = 0;   G(pvs) = 255; B(pvs) = 0;

% implay expects [H W 3 nFrames]
rgb = permute(cat(4, R, G, B), [1 2 4 3]);
player = implay(rgb, opts.Fps);
end
