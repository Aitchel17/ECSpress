function [rgb_img, lim] = plot_make_rgb(img_ch1, img_ch2, opt)
%PLOT_MAKE_RGB  Two channels merged into one RGB image: ch1 magenta, ch2 green.
%   Each channel is stretched onto its own display range before the merge, the way
%   ImageJ adjusts each channel and only then merges. Without lim1/lim2 that range
%   is this image's own percentiles, which is fine for one picture and wrong for
%   two that have to be compared -- a brightness change between them then reads as
%   if the structure had moved. Pass the same limits to both, or take them off the
%   first call's second output.
%
%   IN   img_ch1  MxN | MxNxC   a stack of two or more planes is split on the third
%                               dimension when img_ch2 is not given
%        img_ch2  MxN           [] with a 2-D ch1 makes the green channel zero
%   opt  lim1     1x2 double    display range for ch1. [] = its own percentiles
%        lim2     1x2 double    display range for ch2
%        prctile  1x2 double    which percentiles those are
%   OUT  rgb_img  MxNx3 double  clipped to 1
%        lim      2x2 double    [lim1; lim2] as used, so a second image can match
    arguments
        img_ch1
        img_ch2 = []
        opt.lim1    (1,:) double = []
        opt.lim2    (1,:) double = []
        opt.prctile (1,2) double = [10 99]
    end

    % the split is decided by img_ch2 being absent rather than by nargin, which
    % counts the name-value pairs too
    if isempty(img_ch2) && size(img_ch1, 3) >= 2
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

    lim1 = opt.lim1;
    if isempty(lim1)
        lim1 = prctile(img_ch1(:), opt.prctile);
    end
    lim2 = opt.lim2;
    if isempty(lim2)
        lim2 = prctile(img_ch2(:), opt.prctile);
    end
    lim = [lim1(:)'; lim2(:)'];

    try
        img_ch1_norm = mat2gray(img_ch1, double(lim1));
        img_ch2_norm = mat2gray(img_ch2, double(lim2));
    catch
        % a flat channel has no spread between the two percentiles and mat2gray
        % refuses a zero-width range
        img_ch1_norm = mat2gray(img_ch1);
        img_ch2_norm = mat2gray(img_ch2);
    end

    % R and B carry ch1, G carries ch2, so ch1 reads magenta and ch2 green
    rgb_img = cat(3, img_ch1_norm, img_ch2_norm, img_ch1_norm);
    rgb_img(rgb_img > 1) = 1;
end
