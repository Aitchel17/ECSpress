function [inter, pair_frames] = piv_interleave(S, from, to, halfwin)
%PIV_INTERLEAVE  The two framesets an event spans, interleaved into pairs.
%   Pair p is (from_p, to_p), so a correlation ensemble over this stack averages
%   k measurements of ONE displacement rather than accumulating k steps. The
%   earlier endpoint must be `from`: the pair order sets the sign of uv.
%
% IN   S            H x W x T, the stack. Only its third size is read for clipping
%      from         int, index of the earlier endpoint
%      to           int, index of the later one
%      halfwin      int, frames either side of each endpoint
% OUT  inter        H x W x 2k, (from_1, to_1, from_2, to_2, ...)
%      pair_frames  k x 2 int, [from_idx to_idx] per pair, in original stack
%                   indices. The correlator carries this forward as
%                   cp.pair_frames, so a caller can say which frames a result
%                   came from without recomputing anything
%
%   why  two endpoints and a half-width, not an event row. The arithmetic needs
%        two integers; asking for the row would make this impossible to call on
%        any other pair of frames and would hide what is actually required
%   note k is the SHORTER of the two windows. An endpoint near either end of the
%        recording gets clipped, and pairing across unequal windows would pair
%        frames at different lags
%   note odd slices are the FROM state and even the TO state, which makes this
%        stack a blink comparator when stepped one slice at a time -- wall motion
%        moves the edge, a whole-field drift moves everything

arguments
    S       {mustBeNumeric, mustBeNonempty}
    from    (1,1) double
    to      (1,1) double
    halfwin (1,1) double
end

nT = size(S, 3);
w_from = max(1, from - halfwin):min(nT, from + halfwin);
w_to   = max(1, to   - halfwin):min(nT, to   + halfwin);

% 0. Equal-length index lists for the two framesets
k = min(numel(w_from), numel(w_to));
a = round(linspace(w_from(1), w_from(end), k));
b = round(linspace(w_to(1),   w_to(end),   k));
pair_frames = [a(:), b(:)];

% 1. Interleave into (from_1, to_1, from_2, to_2, ...)
inter = zeros(size(S, 1), size(S, 2), 2 * k, 'like', S);
inter(:, :, 1:2:end) = S(:, :, a);
inter(:, :, 2:2:end) = S(:, :, b);
end
