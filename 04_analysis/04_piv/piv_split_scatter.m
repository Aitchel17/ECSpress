function [se, info] = piv_split_scatter(maps, n_splits)
%PIV_SPLIT_SCATTER  Per-window displacement uncertainty, measured not modelled.
%   Splits the ensemble's correlation planes into two disjoint halves, finds the
%   peak of each half separately, and reads the uncertainty off how far apart the
%   two answers land. Repeats over several random splits and pools.
%
%   Why measure it. The obvious candidates for a per-vector error bar are all
%   proxies that turned out to mean something else here: the correlation peak
%   WIDTH carries the spread of displacements inside the window, not a
%   localisation limit; corr2 stays high on a window with nothing to track; and
%   the Lucas-Kanade covariance sigma_n^2 * M^-1 is derived for gradient-based
%   optical flow, which is not the estimator this pipeline uses. Splitting the
%   ensemble asks the estimator itself, so it needs no model of the estimator.
%
%   What it does NOT see: anything the halves get wrong together. Peak locking,
%   decorrelation bias and the majority-vote collapse of a mixed-sign ensemble
%   all move both halves the same way and leave no trace in the scatter. So this
%   is the RANDOM part of the error only, and a floor on the true uncertainty.
%
%   The scaling. Each half holds about half the pairs, so its own scatter is
%   larger than the full ensemble's by the usual root-N. Halves A and B are
%   independent, so var(A - B) = var(A) + var(B), and the standard error of the
%   FULL ensemble follows from the pair counts -- worked out below rather than
%   assumed, since an odd pair count cannot split evenly.
%
% IN   maps      IA x IA x ny x nx x nPair planes, as piv_corr_ensemble saves them
%      n_splits  random half/half partitions to pool, default 10. With few pairs
%                the distinct partitions run out and the extra draws repeat, which
%                is harmless
% OUT  se        ny x nx, STANDARD ERROR of the ensemble displacement, px. NaN
%                where a half had no peak to fit
%      info      n_pairs, n_splits, factor (rms difference -> se), and
%                rms_diff, the raw half-to-half distance before scaling

arguments
    maps     {mustBeNumeric}
    n_splits (1,1) double {mustBePositive} = 10
end

[s, ~, ny, nx, np] = size(maps);
if np < 2
    error('piv_split_scatter:onePair', ...
        'the ensemble holds %d pair(s); a split needs at least 2.', np);
end

% 0. Flatten to (plane cell, window, pair) so a split is one sum
M  = double(reshape(maps, s*s, ny*nx, np));
nA = floor(np/2);
nB = np - nA;
% 0.1 var(A-B) = sigma1^2 (1/nA + 1/nB) and se = sigma1/sqrt(np), so
%     se = rms(A-B) / sqrt(np*(1/nA + 1/nB)). Worked out from the counts because
%     an odd pair count cannot split evenly
factor = 1 / sqrt(np * (1/nA + 1/nB));

% 1. Pool the half-to-half distance over several partitions
acc = zeros(ny*nx, 1);
cnt = zeros(ny*nx, 1);
rng_state = rng(0);                     % same splits every call, so a rerun matches
for it = 1:n_splits
    p  = randperm(np);
    dA = peak_of(M(:, :, p(1:nA)), s);
    dB = peak_of(M(:, :, p(nA+1:end)), s);
    d  = hypot(dA(:,1) - dB(:,1), dA(:,2) - dB(:,2));
    ok = isfinite(d);
    acc(ok) = acc(ok) + d(ok).^2;
    cnt(ok) = cnt(ok) + 1;
end
rng(rng_state);

rms_diff = NaN(ny*nx, 1);
has = cnt > 0;
rms_diff(has) = sqrt(acc(has) ./ cnt(has));
se = reshape(rms_diff * factor, ny, nx);

info = struct('n_pairs', np, 'n_splits', n_splits, 'factor', factor, ...
    'rms_diff', reshape(rms_diff, ny, nx), 'split', [nA nB]);
end

% ---------------------------------------------------------------------------
function d = peak_of(Msub, s)
%PEAK_OF  Subpixel peak of each summed plane, as [dx dy] from the centre.
%   Three-point Gaussian in each axis, the same estimator PIVlab's subpixfinder 1
%   uses, so the scatter describes the vectors the pipeline actually reports.
    E = sum(Msub, 3);
    E = E - min(E, [], 1);
    [~, i] = max(E, [], 1);
    [r, c] = ind2sub([s s], i(:));
    n = numel(r);
    d = NaN(n, 2);
    inner = r > 1 & r < s & c > 1 & c < s;
    if ~any(inner); return; end
    lg  = @(v) log(max(v, realmin));
    idx = @(rr, cc) sub2ind([s s], rr, cc) + (0:n-1)'*s*s;
    rr = r(:);  cc = c(:);
    a = E(idx(max(rr-1,1), cc));  b = E(idx(rr, cc));  g = E(idx(min(rr+1,s), cc));
    dr = (lg(a) - lg(g)) ./ (2*(lg(a) - 2*lg(b) + lg(g)));
    a = E(idx(rr, max(cc-1,1)));  g = E(idx(rr, min(cc+1,s)));
    dc = (lg(a) - lg(g)) ./ (2*(lg(a) - 2*lg(b) + lg(g)));
    good = inner & isfinite(dr) & isfinite(dc);
    d(good, 1) = cc(good) + dc(good) - (floor(s/2) + 1);
    d(good, 2) = rr(good) + dr(good) - (floor(s/2) + 1);
end
