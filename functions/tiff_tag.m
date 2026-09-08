function tags = tiff_tag(n_frame, n_channel, unit, page_keys, page_step, pixel_density, resolution_unit)
%TIFF_TAG  every tag a page of these TIFFs carries but its own size, the ImageJ block included
%   Caller: io_postsavetiff
%   A copy of tiftag and imagejblock in 00_mdfExtractor's mdf.m; keep the two the same
%
% IN   n_frame          1 x 1 double   pages along the page axis
%      n_channel        1 x 1 double   channels interleaved per page
%      unit             1 x 3 string   the unit of x, y and the page axis
%      page_keys        1 x 2 string   what ImageJ calls the page count and the page
%                                      step: ["slices" "spacing"] for a depth axis,
%                                      ["frames" "finterval"] for a time axis
%      page_step        1 x 1 double   one page in unit(3)
%      pixel_density    1 x 2 double   pixels per resolution_unit, [x y]. TIFF
%                                      stores the reciprocal of a pixel's size,
%                                      so um per pixel arrives here as 10000/um
%      resolution_unit  1 x 1          Tiff.ResolutionUnit.*, or the numeric code
%                                      getTag hands back
% OUT  tags             1 x 1 struct   for Tiff.setTag; the caller adds the size
    tags.ImageDescription = imagejblock(n_frame, n_channel, unit, page_keys, page_step);
    tags.Photometric = Tiff.Photometric.MinIsBlack;
    tags.BitsPerSample = 16;
    tags.SamplesPerPixel = 1;
    tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tags.Compression = Tiff.Compression.None;
    tags.ResolutionUnit = resolution_unit;
    tags.XResolution = pixel_density(1);
    tags.YResolution = pixel_density(2);
    tags.Software = 'MATLAB';
end

function description = imagejblock(n_frame, n_channel, unit, page_keys, page_step)
    %IMAGEJBLOCK  the ImageDescription lines ImageJ reads, as one char row
    %   Caller: tiff_tag
    %
    % IN   n_frame      1 x 1 double   pages along the page axis
    %      n_channel    1 x 1 double   channels
    %      unit         1 x 3 string   the unit of x, y and the page axis
    %      page_keys    1 x 2 string   what ImageJ calls the page count and the page
    %                                  step: ["slices" "spacing"] for a depth axis,
    %                                  ["frames" "finterval"] for a time axis
    %      page_step    1 x 1 double   one page in unit(3)
    % OUT  description  1 x n char
    %
    %   Every line needs its own newline. The page axis is the one thing here a scan
    %   mode decides, and it arrives as page_keys rather than being guessed from unit.
    description = [sprintf('ImageJ=1.54f\n'), ...
                   sprintf('images=%d\n', n_frame * n_channel), ...
                   sprintf('channels=%d\n', n_channel), ...
                   sprintf('%s=%d\n', page_keys(1), n_frame), ...
                   sprintf('unit=%s\n', unit(1)), ...
                   sprintf('yunit=%s\n', unit(2)), ...
                   sprintf('zunit=%s\n', unit(3)), ...
                   sprintf('%s=%g\n', page_keys(2), page_step)];
    if n_channel > 1
        description = [description, sprintf('hyperstack=true\n')];   % ImageJ opens channels x pages as one
        description = [description, sprintf('mode=composite\n')];    % and in colour, not one grey channel at a time
    end
end
