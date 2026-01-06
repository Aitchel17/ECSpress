function rgb_img = plot_make_rgb(img_ch1, img_ch2)
% PLOT_MAKE_RGB Creates a 3-channel RGB image from two input images.
%   Applies 5-95% contrast stretching to each channel.
%   Returns an RGB image where R=img_ch1, G=img_ch2, B=0.

arguments
    img_ch1
    img_ch2 = []
end

% If only one input is provided (e.g. 3D stack), split it
if nargin == 1 && size(img_ch1, 3) >= 2
    img_ch2 = double(img_ch1(:,:,2));
    img_ch1 = double(img_ch1(:,:,1));
else
    img_ch1 = double(img_ch1);
    if isempty(img_ch2)
        img_ch2 = zeros(size(img_ch1));
    else
        img_ch2 = double(img_ch2);
    end
end

try
    lim1 = prctile(img_ch1(:), [1 99]);
    lim2 = prctile(img_ch2(:), [1 99]);

    img_ch1_norm = mat2gray(img_ch1, double(lim1));
    img_ch2_norm = mat2gray(img_ch2, double(lim2));
catch
    % Fallback if data is degenerated (e.g. all zeros)
    img_ch1_norm = mat2gray(img_ch1);
    img_ch2_norm = mat2gray(img_ch2);
end

% Generate RGB Image (Ch1=Magenta (R+B), Ch2=Cyan (G+B))
% R = Ch1
% G = Ch2
% B = Ch1 + Ch2
rgb_img = cat(3, img_ch1_norm, img_ch2_norm, img_ch1_norm + img_ch2_norm);
rgb_img(rgb_img > 1) = 1;

end
