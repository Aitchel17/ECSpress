function [uv_map, corr_map] = vfield_stamp(cp, u, v)
%VFIELD_STAMP  Interrogation grid -> dense H x W maps. No interpolation.
%   A vector belongs to a WINDOW, not a pixel, so this writes it at the one pixel
%   its window centre sits on and leaves every gap NaN. Nothing is invented
%   between vectors: whatever reads the map downstream sees exactly the samples
%   that were measured, and can decide for itself how to fill.
%
%   Takes u and v separately from cp so a caller that has gated the field passes
%   its OWN vectors while still using cp's geometry. cp.utable/vtable are only
%   the default because a raw stamp is worth one line.
%
% IN   cp        corr_planes from piv_corr_ensemble; uses xtable, ytable,
%                typevector and imsize, plus corr for the second output
%      u, v      ny x nx vectors to stamp, default cp.utable / cp.vtable
% OUT  uv_map    H x W x 2, (:,:,1) = u horizontal, (:,:,2) = v vertical
%      corr_map  H x W corr2, at the same pixels and NaN elsewhere

arguments
    cp struct
    u  double = []
    v  double = []
end
if isempty(u); u = cp.utable; end
if isempty(v); v = cp.vtable; end

H = cp.imsize(1);
W = cp.imsize(2);

% 1. Windows PIVlab accepted that still carry a vector
keep = (cp.typevector == 1) & ~isnan(u) & ~isnan(v);
rows = max(1, min(H, round(cp.ytable(keep))));
cols = max(1, min(W, round(cp.xtable(keep))));
idx  = sub2ind([H, W], rows, cols);

% 2. One pixel each
uv_map = NaN(H, W, 2);
uv_map(idx)       = u(keep);
uv_map(idx + H*W) = v(keep);

if nargout > 1
    corr_map = NaN(H, W);
    corr_map(idx) = cp.corr(keep);
end
end
