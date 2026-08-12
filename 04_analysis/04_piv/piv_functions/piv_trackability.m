function [lam2, lam1] = piv_trackability(img, xtable, ytable, window, mask)
%PIV_TRACKABILITY  Shi-Tomasi corner strength per interrogation window.
%   The structure tensor of the image gradients over each window,
%
%       M = [ sum Ix^2   sum IxIy ]
%           [ sum IxIy   sum Iy^2 ]
%
%   whose smaller eigenvalue lam2 says how well the window pins a displacement
%   in BOTH directions. lam2 near zero is an edge or a flat patch: the
%   correlation still peaks, but where it peaks along the degenerate direction
%   is arbitrary. Measured on ev2, lam2 is uncorrelated with how far the window
%   moved (r = -0.06), which is the point -- a stationary structure is trackable
%   and must not be thrown away for being still, which is what a gate on vector
%   length does.
%
%   Read from the image alone, so it is available before any correlation and
%   cannot be confounded by the result it is meant to judge.
%
% IN   img     H x W reference frame, the state the vectors start from
%      xtable  vector-grid x coords
%      ytable  vector-grid y coords
%      window  interrogation window side (px); the same window the last pass used
%      mask    H x W logical, true = pixel excluded, [] = none
% OUT  lam2    size(xtable), smaller eigenvalue. NaN where the window is mostly
%              masked out
%      lam1    larger eigenvalue; lam2/lam1 near 0 marks an edge

arguments
    img    {mustBeNumeric}
    xtable double
    ytable double
    window (1,1) double {mustBePositive}
    mask   = []
end

% 0. Setup
[H, W] = size(img);
if isempty(mask)
    mask = false(H, W);
end
[Ix, Iy] = gradient(double(img));
[ny, nx] = size(xtable);
lam1 = NaN(ny, nx);
lam2 = NaN(ny, nx);
h    = floor(window/2);
need = max(4, round(0.15 * window^2));   % a window mostly masked out has no tensor

% 1. One structure tensor per window
for r = 1:ny
    for c = 1:nx
        yc = round(ytable(r, c));   
        xc = round(xtable(r, c));
        ys = max(1, yc-h) : min(H, yc+h);
        xs = max(1, xc-h) : min(W, xc+h);
        keep = ~mask(ys, xs);
        if nnz(keep) < need
            continue; 
        end
        gx = Ix(ys, xs);   gy = Iy(ys, xs);
        gx = gx(keep);     gy = gy(keep);
        % 1.1 Normalised by the pixel count, so a clipped window is comparable
        M = [gx'*gx, gx'*gy; 
            gx'*gy, gy'*gy] / numel(gx);
        e = sort(eig(M), 'descend');
        lam1(r, c) = e(1);
        lam2(r, c) = e(2);
    end
end
end
