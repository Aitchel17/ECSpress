function [seg_results,T,sizefilt,sizesmooth,fillcleanopt] = SegmentStack(img,varargin)
%SEGMENTSTACK Segments an entire stack of images, img (MXNXS)
%   Useage: [seg_results,T,sizefilt,sizesmooth,fillcleanopt] = SegmentStack(tiff_stack)
%           [seg_results,T,sizefilt,sizesmooth,fillcleanopt] = SegmentStack(tiff_stack,[T],[sizefilt],[sizesmooth],[fillcleanopt])
%   Returns binary MXNXS arrays with the segmentation (seg_results) and
%   the threshold (T), filter size (sizefilt), size of structuring array 
%   used to smooth (sizesmooth), and binary option indicating whether to
%   clean and filter (fillcleanopt). If no options are supplied, sliderseg
%   opens, which helps the user decide what parameters to use for the
%   segmentation. 
%   ------------Summary of potential values supplied for T:----------------
%   0 < T < 1: basic threshold using value supplied for T to binarize each image
%   T > 1: this occurs when the image(s) supplied is/are a unit16. T is
%   normalized with the same normalization as the entire stack. 
%   T = 0: a separate threshold is calculated for each image. 
%   T < 1: a spatially varying adaptive threshold is calculated for each image using adapatthresh and some other magic (see explanation in step 1 below). 
%   T as a vector supplies a different threshold for every image. If the
%   max value in the vector is larger than 1, T is normalized with the same
%   normalization as the entire stack.
% ------Explanation of the various options------------------
% SegementStack creates a binary segmentation indicating an area of
% interest in a stack of images. The segmentation involves the following
% steps: 
% 1.) Initial binarization: the most basic approach is using the same
% threshold for each image in the stack, and this approach is used if a 
% threshold value between 0 and 1 is supplied. 
% If a 0 is supplied, a separate threshold for each image in the stack is 
    % automatically calculated using the MATLAB function graythresh 
    % (implements Otsu's method). The T that is returned is then a 1 X S 
    % vector with the threshold used for each image.
% If an integer less than 0 supplied, a spatially varying threshold is
    % calculated. However, since the spatially varying threshold is not usually
    % useful by itself, that threshold is masked by the inverse of the
    % graythresh binarized image dilated by a disk structural element of a certain size. The size is
    % specified by the absolute value of the supplied T. The T that is returned is then an struct with two fields: 
    % T.T is an M X N X S array with the spatially varing threshold used for each image. 
    % T.dilatesize is the size of the dilating structure used 
% 2.) After binarizing the image, if the "fill holes, clean" option is
% selected, holes in the image are filled and areas of only one pixel are
% removed. 
% 3.) The resulting segmentation is then smoothed by dilating, filling
% holes, then eroding using a disk shaped structural element of size
% "sizesmooth". It's like an open operation, but it has fills holes
% inbetween the dilating and eroding stage. 
% 4.) The last step is to remove areas that are smaller than a certain
% size (number of pixels), specified by the parameter sizefilt. 
% Written 2021 by Kimberly Boster, kboster@ur.rochester.edu
if nargin>1
    T=varargin{1};
    sizefilt=varargin{2};
    sizesmooth=varargin{3};
    fillcleanopt=varargin{4};
else
    sliderseg(img)
    pause(0.1)
    T=input('What threshold should be used? \n Enter 0 to use a seperate threshold calculated for each image. \n Enter an integer < 0 to use a spatially varying threshold for each image. The absolute value of this number is the size of the dilating element.\n');
    sizefilt=input('What size filter should be used? ');
    sizesmooth=input('What size structuring element should be used for smoothing? ');
    fillcleanopt=input('Fill holes and clean? Enter [1/0] ');
end

if isa(img,'uint16') || isa(img,'uint8')
    normI=double(max(img(:)));
    img=double(img)./normI;
    if max(T)>1
        Tuint16=T;
        T=double(T)/normI;    
    end
elseif isa(img,'logical')
    img=double(img);
end

for s=1:size(img,3)
    im=img(:,:,s);
    if length(T)>1
        [~,BW_smooth]=seg(im,T(s),sizefilt,sizesmooth,fillcleanopt);
    elseif T<0
        se=strel('disk',abs(T));
        T_new=adaptthresh(img(:,:,s)) + ~imdilate(imbinarize(img(:,:,s)),se);
        T_new(T_new>1)=1;
        T_mat(:,:,s)=T_new;
        [~,BW_smooth]=seg(im,T_mat(:,:,s),sizefilt,sizesmooth,fillcleanopt);
    elseif T==0
        T_vect(s)=graythresh(im);
        [~,BW_smooth]=seg(im,T_vect(s),sizefilt,sizesmooth,fillcleanopt);
    else
        [~,BW_smooth]=seg(im,T,sizefilt,sizesmooth,fillcleanopt);
    end
    seg_results(:,:,s)=BW_smooth;
    if mod(s,500)==0
        disp(['Segmenting image ' num2str(s) ' of ' num2str(size(img,3)) '.'])
    end
end
if exist('T_vect','var')
    T=T_vect;
    if exist('normI','var')
        T=T*normI;
    end
end
if exist('T_mat','var')
    dilatesize=T;
    clear T
    T.T=T_mat;
    T.dilatesize=dilatesize;
end
if exist('Tuint16','var')
    T=Tuint16;
end
end

