function [tags_struct, tiff_obj] = io_readtiff_alltags()
    % Open a dialog to select a .tiff file
    [file, path] = uigetfile({'*.tiff;*.tif', 'TIFF Files (*.tiff, *.tif)'}); 
    if isequal(file, 0)
        disp('No file selected.');
        tags_struct = [];
        tiff_obj = [];
        return;
    end
    full_path = fullfile(path, file);

    % Open the TIFF file
    tiff_obj = Tiff(full_path, 'r');

    % List of tags to retrieve
    tags_list = {
        'SubFileType', 'ImageWidth', 'ImageLength', 'BitsPerSample', ...
        'Compression', 'Photometric', 'Thresholding', 'FillOrder', ...
        'DocumentName', 'ImageDescription', 'Make', 'Model', ...
        'StripOffsets', 'Orientation', 'SamplesPerPixel', 'RowsPerStrip', ...
        'StripByteCounts', 'MinSampleValue', 'MaxSampleValue', ...
        'XResolution', 'YResolution', 'PlanarConfiguration', 'PageName', ...
        'XPosition', 'YPosition', 'Group3Options', 'Group4Options', ...
        'ResolutionUnit', 'PageNumber', 'TransferFunction', 'Software', ...
        'DateTime', 'Artist', 'HostComputer', 'WhitePoint', ...
        'PrimaryChromaticities', 'ColorMap', 'HalfToneHints', 'TileWidth', ...
        'TileLength', 'TileOffsets', 'TileByteCounts', 'SubIFD', ...
        'InkSet', 'InkNames', 'NumberOfInks', 'DotRange', 'TargetPrinter', ...
        'ExtraSamples', 'SampleFormat', 'SMinSampleValue', 'SMaxSampleValue', ...
        'YCbCrCoefficients', 'YCbCrSubSampling', 'YCbCrPositioning', ...
        'ReferenceBlackWhite', 'XMP', 'ImageDepth', 'Copyright', ...
        'ModelPixelScaleTag', 'RichTIFFIPTC', 'RPCCoefficientTag', ...
        'ModelTiepointTag', 'ModelTransformationTag', 'Photoshop', ...
        'ICCProfile', 'GeoKeyDirectoryTag', 'GeoDoubleParamsTag', ...
        'GeoASCIIParamsTag', 'SToNits', 'JPEGQuality', 'JPEGColorMode', ...
        'ZipQuality', 'SGILogDataFmt'
    };

    % Initialize structure to hold tags
    tags_struct = struct();

    % Iterate through each tag and try to retrieve its value
    for i = 1:numel(tags_list)
        tag = tags_list{i};
        try
            value = tiff_obj.getTag(tag);
            tags_struct.(tag) = value;
        catch
            % Handle missing or unreadable tags
            tags_struct.(tag) = 'Not available';
        end
    end

    % Check for ImageJ-specific metadata in the ImageDescription tag
    if isfield(tags_struct, 'ImageDescription') && ...
       contains(tags_struct.ImageDescription, 'ImageJ', 'IgnoreCase', true)
        tags_struct.ImageJMetadata = parse_imagej_metadata(tags_struct.ImageDescription);
    else
        tags_struct.ImageJMetadata = 'No ImageJ metadata found.';
    end

    % Iterate through all directories to count slices (optional)
    directory_count = 1;
    while ~tiff_obj.lastDirectory()
        tiff_obj.nextDirectory();
        directory_count = directory_count + 1;
    end
    tags_struct.DirectoryCount = directory_count;

    % Close the TIFF file
    tiff_obj.close();

    % Output the number of directories for debugging
    fprintf('File contains %d directories (image slices).\n', directory_count);
end

function metadata_struct = parse_imagej_metadata(description)
    % Parses ImageJ metadata from the ImageDescription tag
    metadata_lines = splitlines(description);
    metadata_struct = struct();
    for i = 1:numel(metadata_lines)
        line = metadata_lines{i};
        if contains(line, '=')
            [key, value] = strtok(line, '=');
            value = strtrim(value(2:end)); % Remove '=' and trim whitespace
            metadata_struct.(matlab.lang.makeValidName(key)) = value;
        end
    end
end
