function reconstructed = analyze_affine_reverse(rotatedCrop, originalSize, vertices)
% Step 1: Compute rotation properties
lvec = vertices(2,:)' - vertices(1,:)';
ltheta = atan2(lvec(2), lvec(1));

% Step 2: Compute rectangle corners from forward step
[h, w, ~] = size(rotatedCrop);
pvec = [-lvec(2), lvec(1)] / norm(lvec) * (h / 2);
rectVertices = [vertices(1,:) + pvec; vertices(2,:) + pvec;
                vertices(2,:) - pvec; vertices(1,:) - pvec];

% Step 3: Reverse rotation (rotate back to upright orientation)
rot = [cos(ltheta), -sin(ltheta); sin(ltheta), cos(ltheta)]; 
tform = affinetform2d([rot, [0;0]; 0 0 1]);

% The rotatedCrop was centered on Y=0 (from -h/2 to h/2) and X started at 0
Rin = imref2d([h, w], [0 w], [-h/2, h/2]);
invRotated = imwarp(rotatedCrop, Rin, tform);

% Step 4: Compute center-based alignment offset
[inv_h, inv_w, ~] = size(invRotated);
invCenter = [inv_w; inv_h] / 2;
rectCenter = mean(rectVertices, 1)';
offset = round(rectCenter - invCenter);

% Step 5: Place into blank canvas
reconstructed = zeros(originalSize, 'like', rotatedCrop);
yRange = offset(2) + (1:inv_h);
xRange = offset(1) + (1:inv_w);

% Clip to bounds
validY = yRange > 0 & yRange <= originalSize(1);
validX = xRange > 0 & xRange <= originalSize(2);

% Insert rotated image into canvas
reconstructed(yRange(validY), xRange(validX), :) = invRotated(validY, validX, :);
end
