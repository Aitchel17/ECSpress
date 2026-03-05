function io_postsavetiff(zstack, save_path, resolution)
    % Save a 4D array [ny, nx, frames, channels] as an ImageJ Hyperstack TIFF.
    % Uses a fully manual binary writer to guarantee correct StripOffsets.
    % This bypasses the MATLAB Tiff class which has a known StripOffsets bug.

    x_res = resolution(1);  % um/pixel
    y_res = resolution(2);  % um/pixel
    z_res = 1/resolution(3);  % seconds/frame (for spacing tag)

    % --- Dimensions ---
    ny          = size(zstack, 1);
    nx          = size(zstack, 2);
    nframes     = size(zstack, 3);
    nchannels   = max(1, size(zstack, 4));
    ntotal      = nframes * nchannels;

    % --- Build ImageJ ImageDescription ---
    desc = sprintf(['ImageJ=1.54f\n' ...
                    'images=%d\n'    ...
                    'channels=%d\n'  ...
                    'frames=%d\n'    ...
                    'hyperstack=true\n' ...
                    'mode=composite\n'  ...
                    'unit=um\n'         ...
                    'pixelWidth=%g\n'   ...
                    'pixelHeight=%g\n'  ...
                    'zunit=sec\n'       ...
                    'finterval=%g\n'],  ...
                    ntotal, nchannels, nframes, x_res, y_res, z_res);
    % Pad to even size per TIFF spec
    if mod(length(desc), 2) ~= 0, desc(end+1) = char(0); end
    desc_bytes = uint8(desc);
    desc_len   = numel(desc_bytes);
    disp(desc);

    % --- Convert to uint16 ---
    zstack = uint16(zstack);
    bpp          = 2;               % bytes per pixel (uint16)
    frame_bytes  = ny * nx * bpp;

    % =================================================================
    % File layout (all offsets are 0-based byte positions):
    %
    %   [0             ]  8-byte TIFF header
    %   [8             ]  ImageDescription string  (desc_len bytes)
    %   [8+desc_len    ]  XResolution rational: uint32×2 (8 bytes)
    %   [8+desc_len+8  ]  YResolution rational: uint32×2 (8 bytes)
    %   [8+desc_len+16 ]  pixel data, ntotal frames back-to-back
    %   [after data    ]  ntotal IFD blocks chained
    %
    % Each IFD = 2 (entry count) + 11*12 (tags) + 4 (next ptr) = 138 bytes
    % =================================================================

    HEADER_SIZE = 8;
    DESC_OFFSET = HEADER_SIZE;                    % 0-based
    % After ImageDescription: 16 bytes for two rationals (XRES, YRES)
    DATA_START  = DESC_OFFSET + desc_len + 16;

    NUM_TAGS    = 12;
    IFD_BYTES   = 2 + NUM_TAGS * 12 + 4;         % 150

    IFD_START   = DATA_START + ntotal * frame_bytes;

    % All StripOffsets and IFD positions known before writing any byte
    strip_offsets = uint32(DATA_START + (0:ntotal-1) .* frame_bytes);
    ifd_positions = uint32(IFD_START  + (0:ntotal-1) .* IFD_BYTES);

    % XResolution/YResolution as proper TIFF RATIONAL (type=5):
    % Since ImageDescription sets unit=um, ImageJ reads XResolution as px/um.
    % So encode: 1/x_res px/um = 1000000 / round(x_res*1000000)
    xres_num = uint32(1000000);
    xres_den = uint32(round(x_res * 1000000));   % e.g. 0.4125*1e6 = 412500
    yres_num = uint32(1000000);
    yres_den = uint32(round(y_res * 1000000));
    
    % File layout for rational values (stored right after ImageDescription)
    XRES_RAT_OFFSET = uint32(DESC_OFFSET + desc_len);      % 8 bytes
    YRES_RAT_OFFSET = uint32(DESC_OFFSET + desc_len + 8);  % 8 bytes

    % --- Write ---
    fid = fopen(save_path, 'wb');
    if fid < 0, error('io_postsavetiff: Cannot open %s for writing', save_path); end

    try
        % Header
        fwrite(fid, uint8([73 73]),      'uint8' );  % 'II' little-endian
        fwrite(fid, uint16(42),          'uint16');  % TIFF magic
        fwrite(fid, ifd_positions(1),    'uint32');  % offset to first IFD

        % ImageDescription string (shared by all IFDs via DESC_OFFSET)
        fwrite(fid, desc_bytes, 'uint8');

        % Rational values for resolution tags (stored right after ImageDescription)
        fwrite(fid, xres_num, 'uint32'); fwrite(fid, xres_den, 'uint32');  % XRes rational
        fwrite(fid, yres_num, 'uint32'); fwrite(fid, yres_den, 'uint32');  % YRes rational

        % Pixel data — all frames, channel-interleaved per frame
        for i = 1:nframes
            for c = 1:nchannels
                if nchannels > 1
                    frame = zstack(:,:,i,c);
                else
                    frame = zstack(:,:,i);
                end
                fwrite(fid, frame.', 'uint16');  % .' = transpose for row-major TIFF scanlines
            end
        end

        % IFDs
        for k = 1:ntotal
            so   = strip_offsets(k);
            next = uint32(0);
            if k < ntotal, next = ifd_positions(k+1); end

            fwrite(fid, uint16(NUM_TAGS), 'uint16');  % entry count

            % Tags in ascending numeric order (TIFF spec requirement)
            wtag(fid, 256, 3, 1,              uint32(nx));         % ImageWidth
            wtag(fid, 257, 3, 1,              uint32(ny));         % ImageLength
            wtag(fid, 258, 3, 1,              uint32(16));         % BitsPerSample
            wtag(fid, 259, 3, 1,              uint32(1));          % Compression=None
            wtag(fid, 262, 3, 1,              uint32(1));          % MinIsBlack
            wtag(fid, 270, 2, uint32(desc_len), uint32(DESC_OFFSET)); % ImageDescription
            wtag(fid, 273, 4, 1,              so);                 % StripOffsets
            wtag(fid, 278, 3, 1,              uint32(ny));         % RowsPerStrip = ny (single strip)
            wtag(fid, 279, 4, 1,              uint32(frame_bytes));% StripByteCounts
            wtag(fid, 282, 5, 1, XRES_RAT_OFFSET);     % XResolution (RATIONAL, px/um)
            wtag(fid, 283, 5, 1, YRES_RAT_OFFSET);     % YResolution (RATIONAL, px/um)
            wtag(fid, 296, 3, 1, uint32(1));             % ResolutionUnit=1 (no unit, defer to ImageJ um)

            fwrite(fid, next, 'uint32');  % next IFD pointer
        end

    catch ME
        fclose(fid);
        rethrow(ME);
    end

    fclose(fid);
    fprintf('Saved %d-ch × %d-fr stack → %s\n', nchannels, nframes, save_path);
end

% Write one 12-byte TIFF IFD entry (little-endian)
function wtag(fid, tag, type, count, value)
    fwrite(fid, uint16(tag),   'uint16');
    fwrite(fid, uint16(type),  'uint16');
    fwrite(fid, uint32(count), 'uint32');
    fwrite(fid, uint32(value), 'uint32');
end