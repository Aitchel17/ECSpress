function io_postsavetiff(zstack, save_path, resolution)
%IO_POSTSAVETIFF  a [ny nx frames channels] array as one ImageJ hyperstack TIFF, a frame's channels side by side
%   Caller: line_fwhm.reconstruction_overlay, stateintegration_test
%
% IN   zstack      H x W x T or H x W x T x C numeric   written as uint16, cast as it is
%      save_path   1 x n char     the .tif
%      resolution  1 x 3 double   [x y] um per pixel, then frames per second
    n_frame = size(zstack, 3);
    n_channel = size(zstack, 4);
    pixel_density = 10000 ./ resolution(1:2);          % pixels per cm
    tags = tiff_tag(n_frame, n_channel, ["um" "um" "sec"], ["frames" "finterval"], 1 / resolution(3), ...
        pixel_density, Tiff.ResolutionUnit.Centimeter);
    tags.ImageLength = size(zstack, 1);
    tags.ImageWidth = size(zstack, 2);
    tags.RowsPerStrip = size(zstack, 1);
    handle = Tiff(save_path, 'w8');
    for k = 1:n_frame
        for channel = 1:n_channel
            frame = uint16(zstack(:, :, k, channel));
            tiff_writepage(handle, tags, (k - 1) * n_channel + channel, frame);
        end
    end
    handle.close();
    fprintf('%s: %d frames x %d channels\n', save_path, n_frame, n_channel);
end
