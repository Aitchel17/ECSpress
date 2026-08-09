function band = analyze_interp2(image, vertices, width)
%ANALYZE_INTERP2  interp2-based rotated-band extraction.
%   Drop-in alternative to analyze_affine_rotate: samples a width-wide band
%   along the line vertices(1,:)->vertices(2,:) by bilinear interpolation of the
%   ORIGINAL image (no imwarp), returning a straightened stack.
%
%   Inputs (same as analyze_affine_rotate):
%       image    - 2D image or 3D stack (x,y[,frames])
%       vertices - 2x2 [x1 y1; x2 y2] line endpoints
%       width    - band width (px)
%   Output:
%       band     - [floor(width) x (ceil(len)-2) x frames]
%                  dim1 = across width, dim2 = along line, dim3 = frames.
%                  Same orientation as analyze_affine_rotate (project over dim1
%                  to get the kymograph). Out-of-image samples are NaN.
%
%   vs analyze_affine_rotate: samples only the band points (not a bbox warp),
%   no interpolation on pixels outside the band, no rotation edge artifacts.
%   Column count is matched (2:end-1 trim) for shape parity only.

    lvec    = vertices(2,:) - vertices(1,:);
    llength = norm(lvec);
    u = lvec / max(llength, eps);          % along-line unit vector
    n = [-u(2), u(1)];                     % perpendicular unit vector

    outWidth  = max(ceil(llength), 3);     % along-line samples (match affine)
    outHeight = max(floor(width), 1);      % across-width samples (match affine)

    % pixel-center coordinates matching imwarp's OutputView mapping:
    %   column j -> along-line distance j-0.5 (line origin at vertex 1)
    %   row i    -> perpendicular offset (i-0.5) - outHeight/2 (line at center)
    tcol = (1:outWidth) - 0.5;
    prow = ((1:outHeight) - 0.5) - outHeight/2;
    [T, P] = meshgrid(tcol, prow);         % [outHeight x outWidth]

    p1 = vertices(1,:);
    SX = p1(1) + u(1)*T + n(1)*P;          % sample x (original image coords)
    SY = p1(2) + u(2)*T + n(2)*P;          % sample y

    nframes = size(image, 3);
    band = nan(outHeight, outWidth, nframes);
    for f = 1:nframes
        band(:,:,f) = interp2(double(image(:,:,f)), SX, SY, 'linear', NaN);
    end

    % match analyze_affine_rotate edge trim (drop first/last along-line column)
    band = band(:, 2:end-1, :);
end
