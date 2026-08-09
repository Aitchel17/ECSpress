function halfwidth = piv_peak_width(maps)
% IN   maps       P x P x ny x nx x npairs correlation planes from piv_corr_ensemble
%                 (corr_planes.maps), one plane per interrogation window per pair
% OUT  halfwidth  ny x nx equivalent half-max radius (px) of the ENSEMBLED peak,
%                 NaN where the window has no usable plane
%
%   How far the displacement is localised. One plane cell is one px of
%   displacement, so this is the uncertainty a vector's length must beat to mean
%   anything. Blur and low contrast keep corr2 high while widening the peak, so
%   this catches a failure the correlation threshold cannot see.

arguments
    maps {mustBeNumeric}
end

% 1. Ensembling the pairs, the same sum the peak was fitted to
E = sum(maps, 5);

% 2.1. Referencing each plane to its own floor
E = E - min(min(E, [], 1), [], 2);
% 2.2. Counting the cells at or above half the peak
pk    = max(max(E, [], 1), [], 2);
area  = squeeze(sum(sum(E >= 0.5 * pk, 1), 2));

% 3. Converting that area to the radius of the equivalent disc
halfwidth = sqrt(area / pi);
halfwidth(squeeze(pk) <= 0 | isnan(squeeze(pk))) = NaN;
end
