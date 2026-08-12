function uv_map = piv_stamp(xyuv, imsize, keep)
%PIV_STAMP  Interrogation grid -> dense H x W map. No interpolation.
%   A vector belongs to a WINDOW, not a pixel, so this writes it at the one pixel
%   its window centre sits on and leaves every gap NaN. Nothing is invented
%   between vectors: whatever reads the map sees exactly the samples that were
%   measured and decides for itself how to fill.
%
%   The dense form is a VIEW, not the record. It is self-describing -- the pixel
%   index is the position -- which is why it composes with the H x W masks
%   (coremask, exclmask, a label map) with no bookkeeping. It is also ~19x larger
%   than the grid it came from and rounds the window centre to a pixel, so the
%   sparse form is what gets stored. See PIV_PLAN.md
%
% IN   xyuv    ny x nx x 4 float px  (:,:,1:2) = [x y] window centre,
%                                    (:,:,3:4) = [u v] displacement
%      imsize  1 x 2 int             [H W]. NOT derivable from xyuv -- the grid
%                                    stops half a window short of every edge
%      keep    ny x nx logical       which windows to write. The caller stacks
%                                    every verdict (typevector, the gate masks)
%                                    and applies them here, once
% OUT  uv_map  H x W x 2 float px    (:,:,1) = u horizontal, (:,:,2) = v vertical
%
% UNIT  whatever xyuv(:,:,3:4) was in. The scale belongs to the caller
%   why    keep is an argument rather than NaN in u and v : NaN cannot say WHY a
%          window is absent, and the caller that wants both a gated and an
%          ungated map would have to keep two copies of the field to get them.
%          One field, two masks, two calls
%   err    the predecessor clamped out-of-range centres into the image with
%          max(1, min(H, ...)). A centre outside the frame is a grid that does
%          not belong to this imsize, and folding it onto the border silently
%          stacks vectors on one pixel. It errors now. See CLAUDE_LOG.md

arguments
    xyuv   (:,:,4) double
    imsize (1,2) double {mustBeInteger, mustBePositive}
    keep   (:,:) logical
end
if ~isequal(size(keep), [size(xyuv, 1), size(xyuv, 2)])
    error('piv_stamp:keepShape', ...
        'keep is %s, the grid is %d x %d.', ...
        mat2str(size(keep)), size(xyuv, 1), size(xyuv, 2));
end

H = imsize(1);
W = imsize(2);

% 1. The kept windows, at the pixel their centre rounds to
x = xyuv(:,:,1);
y = xyuv(:,:,2);
u = xyuv(:,:,3);
v = xyuv(:,:,4);
rows = round(y(keep));
cols = round(x(keep));
if any(rows < 1 | rows > H | cols < 1 | cols > W)
    error('piv_stamp:outsideFrame', ...
        'a window centre rounds outside [%d %d]: rows %g..%g, cols %g..%g.', ...
        H, W, min(rows), max(rows), min(cols), max(cols));
end
idx = sub2ind([H, W], rows, cols);
if numel(unique(idx)) ~= numel(idx)
    error('piv_stamp:pixelCollision', ...
        '%d of %d window centres round onto a pixel another one already took.', ...
        numel(idx) - numel(unique(idx)), numel(idx));
end

% 2. One pixel each
uv_map = NaN(H, W, 2);
uv_map(idx)       = u(keep);
uv_map(idx + H*W) = v(keep);
end
