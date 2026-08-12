function [d, ok] = vfield_divergence_plane(xi, yi, ui, vi)
%VFIELD_DIVERGENCE_PLANE  Divergence du/dx+dv/dy from scattered vectors by LSQ plane fit.
%   Fits  u ~ a0+a1*x+a2*y  and  v ~ b0+b1*x+b2*y  to one cell's vectors and
%   returns a1+b2. The kernel vfield_divergence calls once per cell, whatever shape
%   the cell map gave that cell. Units: px displacement / px, per frame.
%
% IN   xi, yi  pixel coordinates of the cell's vectors (any shape, used as column)
%      ui, vi  u, v displacement components at those points (any shape)
% OUT  d       scalar divergence du/dx+dv/dy; NaN when the spread is degenerate
%      ok      true when the design matrix is full rank; false for a collinear or
%              <3-point spread. Callers that skip degenerate cells test ~ok.

    % 1. Design matrix; full rank only for a genuine 2-D stencil
    A  = [ones(numel(xi), 1), xi(:), yi(:)];
    ok = rank(A) >= 3;
    if ~ok
        d = NaN;
        return;
    end

    % 2. Fit each component and add the two in-plane gradients
    cu = A \ ui(:);          % [const; du/dx; du/dy]
    cv = A \ vi(:);          % [const; dv/dx; dv/dy]
    d  = cu(2) + cv(3);
end
