function analyze_savesimpletif(img3D, filename)
    img3D = uint16(img3D);  % Adjust data type if needed
    t = Tiff(filename, 'w');

    tagstruct.ImageLength = size(img3D, 1);
    tagstruct.ImageWidth = size(img3D, 2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.Software = 'MATLAB';

    for k = 1:size(img3D, 3)
        setTag(t, tagstruct);
        write(t, img3D(:,:,k));
        if k < size(img3D, 3)
            writeDirectory(t);
        end
    end

    close(t);
end
