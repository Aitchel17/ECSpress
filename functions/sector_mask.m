function [mask, n_cell] = sector_mask(sectormap, rings, wedges)
%SECTOR_MASK  The binary mask of chosen (ring, wedge) cells of a sector_polar map.
%   Reads channels 2 and 3, so it never touches the fused label and there is no
%   ordering to get wrong.
%
% IN   sectormap  H x W x 3 double  from sector_polar
%      rings      index vector into 1..n_rings.   [] or ':' = every ring
%      wedges     index vector into 1..n_sectors. [] or ':' = every wedge
% OUT  mask       H x W bool  true = pixel is in ANY requested cell
%      n_cell     int         cells that actually contributed a pixel
%
%   why   no grid or info argument : the counts are in the map, so asking for them
%         separately would be one more thing a caller could get out of step
%   caution  pass the map you actually USED. vfield_divergence_polar mutates it
%         before measuring -- exclmask blanks pixels (the cell shrinks) and
%         sector_validate blanks whole labels (the cell disappears). One cell's
%         raw and used pixel counts can be far apart, and neither call errors
%         see FINDINGS.md
%   note  a cell that is not in the map gives an all-false mask, not an error, and
%         n_cell counts the ones that did contribute. That is deliberate : the map
%         shrinks under exclmask and sector_validate, so a loop written against
%         the raw partition legitimately asks for rings that are no longer there
%
% EXAMPLE
%   [map, geo] = sector_polar(coremask, 20, 6);
%   band = sector_mask(map, 2:4, []);     % rings 2-4, every wedge
%   one  = sector_mask(map, 3, 2);        % just ring 3, wedge 2

arguments
    sectormap {mustBeNumeric, mustBeNonempty}
    rings   {mustBePositive, mustBeInteger} = []
    wedges  {mustBePositive, mustBeInteger} = []
end

if size(sectormap, 3) ~= 3
    error('sector_mask:channels', ...
        'sectormap must be H x W x 3 from sector_polar; got %s.', ...
        mat2str(size(sectormap)));
end

ring_ch  = sectormap(:,:,2);
wedge_ch = sectormap(:,:,3);

% err  an upper-bound guard here was WRONG and errored in normal use. The bound
%      would have to come from the map, but the map SHRINKS : exclmask and
%      sector_validate blank whole rings, so a loop written against the raw
%      partition asks for ring 10 on a map that now stops at 8. The guard could
%      not tell that from a typo. Out of range now returns an empty mask, and
%      n_cell is how a caller notices
% note  the arguments block still rejects 0, negatives and fractions, which is the
%       typo it CAN catch
if isempty(rings) || isequal(rings, ':')
    rings = 1:max(ring_ch(isfinite(ring_ch)), [], 'all');
end
if isempty(wedges) || isequal(wedges, ':')
    wedges = 1:max(wedge_ch(isfinite(wedge_ch)), [], 'all');
end

mask = ismember(ring_ch, rings) & ismember(wedge_ch, wedges);

n_cell = 0;
for r = rings(:)'
    for w = wedges(:)'
        if any(ring_ch == r & wedge_ch == w, 'all')
            n_cell = n_cell + 1;
        end
    end
end
end
