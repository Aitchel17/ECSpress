function [sectormap, counts] = sector_validate(sectormap, valid, n_labels, min_count, min_span)
% IN   sectormap  HxW label map from sector_block / sector_polar
%      valid      HxW logical, true where the measurement is usable. For a PIV
%                 field: ~isnan(u) & ~isnan(v) & ~exclmask
%      n_labels   how many labels the partition DEFINES, so counts lines up with
%                 it. For sector_block that is info.n_labels, an exact count. For
%                 sector_polar it is n_rings * n_sectors, which the caller builds
%                 from max of the ring channel -- see sector_polar
%                 [] = read it off the map, which comes up SHORT whenever the
%                 highest-numbered cell holds no pixels
%      min_count  samples a sector needs, default 4
%      min_span   distinct rows AND columns its samples must span, default 2
% OUT  sectormap  same map with failing sectors cleared to 0
%      counts     n_labels x 1 samples per sector, before the drop

arguments
    sectormap {mustBeNumeric}
    valid     logical
    n_labels  double = []
    min_count (1,1) double = 4
    min_span  (1,1) double = 2
end

% 0. Setup
lab = double(sectormap);
if ~isequal(size(lab), size(valid))
    error('sector_validate:sizeMismatch', 'valid must match the sectormap frame.');
end
if isempty(n_labels); n_labels = max(lab(:)); end
counts = zeros(n_labels, 1);
[X, Y] = meshgrid(1:size(lab, 2), 1:size(lab, 1));

% 1. Testing every sector, clearing the ones that cannot carry a fit
for k = 1:n_labels
    here = lab == k;
    sel  = valid & here;
    counts(k) = nnz(sel);
    % 1.1. Too few samples, or samples that fall on a line
    if counts(k) < min_count || ...
       numel(unique(X(sel))) < min_span || numel(unique(Y(sel))) < min_span
        sectormap(here) = 0;
    end
end
end
