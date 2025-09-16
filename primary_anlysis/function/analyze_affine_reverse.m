function reconstructed = analyze_affine_reverse(rotatedCrop, originalSize, vertices)
% Step 1: Compute rotation properties
lvec = vertices(2,:) - vertices(1,:);
ltheta = atan2(lvec(2), lvec(1));

% Step 2: Compute rectangle corners from forward step
h = size(rotatedCrop, 1);
pvec = [-lvec(2), lvec(1)] / norm(lvec) * (h / 2);
rectVertices = [vertices(1,:) + pvec; vertices(2,:) + pvec;
                vertices(2,:) - pvec; vertices(1,:) - pvec];

% Step 3: Reverse rotation (rotate back to upright orientation)
rot = [cos(ltheta), -sin(ltheta); sin(ltheta), cos(ltheta)];
tform = affinetform2d([rot, [0;0]; 0 0 1]);
Rin = imref2d(size(rotatedCrop), [0 size(rotatedCrop,2)], [-h/2 h/2]);
invRotated = imwarp(rotatedCrop, Rin, tform);

% Step 4: Compute center-based alignment offset
[w, h] = deal(size(invRotated, 2), size(invRotated, 1));
invCenter = [w; h] / 2;
rectCenter = mean(rectVertices, 1)';
offset = round(rectCenter - invCenter);

% Step 5: Place into blank canvas
reconstructed = zeros(originalSize, 'like', rotatedCrop);
yRange = offset(2) + (1:h);
xRange = offset(1) + (1:w);

% Clip to bounds
validY = yRange > 0 & yRange <= originalSize(1);
validX = xRange > 0 & xRange <= originalSize(2);

% Insert rotated image into canvas
reconstructed(yRange(validY), xRange(validX), :) = invRotated(validY, validX, :);
end
