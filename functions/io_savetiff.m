function io_savetiff(zstack, info, img_ch)
    % Save a 3D array (zstack) as a multi-page TIFF file with ImageJ-compatible metadata.
    %
    % Input:
    %   - zstack: 3D matrix of the image stack to save
    %   - info: Struct containing metadata such as resolution and file paths.

    % Create the save directory if it does not exist
    save_folder = fullfile(info.mdfPath, info.mdfName(1:end-4));
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    % Construct full file path
    save_path = fullfile(save_folder, [info.mdfName(1:end-4),sprintf('_ch%d.tif',img_ch)]);

    % Extract resolutions
    x_res = str2double(info.objpix(1:end-2)); % Pixel size in X-direction
    y_res = x_res; % Pixel size in Y-direction
    z_res = 1 / info.savefps; % Z-slice spacing

    % Initialize the TIFF file
    t = Tiff(save_path, 'w');

    % Metadata string for ImageJ
    ImgJ_ver = sprintf('ImageJ=1.54f\n');
    num_img = sprintf('images=%d\n',size(zstack,3));
    num_ch = sprintf('channels=1\n');
    num_frames = sprintf('frames=%d\n',size(zstack,3));
    unit = sprintf('unit=um\n');
    zunit = sprintf('zunit=sec\n');
    spacing = sprintf('spacing=%d\n',z_res);

   ImageDescription = [ImgJ_ver,num_img,num_ch,...
        num_frames,unit,zunit,spacing];
   disp(ImageDescription)



        
    % Loop through each slice of the stack
    for i = 1:size(zstack, 3)
        % Extract the i-th frame
        frame = zstack(:, :, i);
        
        % Convert to uint16 if necessary
        if ~isinteger(frame)
            frame = uint16(frame); % Or uint8(frame), depending on your data range
        end

        tagstruct.ImageLength = size(zstack, 1);
        tagstruct.ImageWidth = size(zstack, 2);
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 16; % Adjust if uint8
        tagstruct.SamplesPerPixel = 1;
        tagstruct.RowsPerStrip = size(zstack, 1);
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.ResolutionUnit = 1;
        tagstruct.XResolution = 1/x_res; % Convert microns to cm
        tagstruct.YResolution = 1/y_res; % Convert microns to cm
        tagstruct.Software = 'MATLAB';
        tagstruct.ImageDescription = ImageDescription;
        setTag(t,'StripOffsets',8);

        % Set tags and write the frame
        t.setTag(tagstruct);
        t.write(frame);

        % Append mode for subsequent slices
        if i < size(zstack, 3)
            t.writeDirectory();
        end
    end

    % Close the TIFF file
    t.close();

    % Notify user of the save location
    fprintf('Saved zstack with metadata to %s\n', save_path);
end