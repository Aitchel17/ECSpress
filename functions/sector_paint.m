function map = sector_paint(sectormap, values)
% IN   sectormap  HxW label map, 0 = no sector
%      values     n_labels x 1, one value per label
% OUT  map        HxW, each pixel filled with its sector's value. A lookup, not an
%                 interpolation: a sector is one flat value, so the map's
%                 resolution is the sector, not the pixel

arguments
    sectormap {mustBeNumeric}
    values    {mustBeNumeric}
end

% 0. Setup
lab    = double(sectormap);
region = lab > 0;
if max(lab(:)) > numel(values)
    error('sector_paint:tooFewValues', ...
        'sectormap goes up to label %d but only %d values were given.', ...
        max(lab(:)), numel(values));
end

% 1. Looking each pixel's value up by its label
map = NaN(size(lab));
map(region) = values(lab(region));
end
