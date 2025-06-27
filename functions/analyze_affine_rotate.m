function rotatedCrop = analyze_affine_rotate(image, vertices, width)
% 1. line vector property (lvec, llength, ltheta) for rotation
% 2. boundary coordinate, and crop image
% 3. translation origin (tori) for cropped image (lower left vertices)
%   rotated by rotation matrix with ltheta --> affine matrix
% 4. scaling factor, matching final output (one pixel drop for rotation
% boundary blank)
% 5. apply rotation
%1
lvec = vertices(2,:) - vertices(1,:);
llength = norm(lvec);
ltheta = atan2(lvec(2), lvec(1));
%2
pvec = [-lvec(2), lvec(1)] / llength * (width / 2); % 1/2pi rotate unit lvec and scale as half width
rectVertices = [vertices(1,:) + pvec;   vertices(2,:) + pvec;
                vertices(2,:) - pvec;   vertices(1,:) - pvec];
x1 = floor(min(rectVertices(:,1))); y1 = floor(min(rectVertices(:,2)));
x2 = ceil(max(rectVertices(:,1)));  y2 = ceil(max(rectVertices(:,2)));
cropped = image(y1:y2, x1:x2,:);

%3
tori = vertices(1,:) - [x1, y1];
rot = [cos(-ltheta), -sin(-ltheta); sin(-ltheta), cos(-ltheta)];
tshift = -rot * tori';
afmat = [rot, tshift; 0 0 1];
tform = affinetform2d(afmat);

% 4
outputWidth = ceil(llength);
outputHeight = ceil(width - 1);
ymin = -outputHeight / 2;
ymax =  outputHeight / 2;
Rout = imref2d([outputHeight, outputWidth], [0 outputWidth], [ymin ymax]);

% 5
rotatedCrop = imwarp(cropped, tform, 'OutputView', Rout);
rotatedCrop = rotatedCrop(:,2:end-1,:);
end
