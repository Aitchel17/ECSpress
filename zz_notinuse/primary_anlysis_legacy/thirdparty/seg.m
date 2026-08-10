function [imbw,BW_final]=seg(im,T,sizefilt,sizesmooth,fillcleanopt)
% SEG segments a single grayscacle (0-1) image im using the supplied parameters.
% T is always a single number between 0 and 1 indicating the threshold. 
% See SegmentStack for a description of the other parameters. Called by
% SegmentStack.m and sliderseg.m.
% kboster@ur.rochester.edu

    imbw=imbinarize(im,T);
    if fillcleanopt==1
%         imbw=padarray(imbw,[0 1],1); % add a border to the left and right to include "holes" on the edges
        imbw=imfill(imbw,'holes');
%         imbw=imbw(:,2:end-1); % trim off pad
        imbw=bwmorph(imbw,'clean');
    end
    if sizesmooth>0
        seD=strel('disk',sizesmooth);
        BW_smooth=imerode(imfill(imdilate(imbw,seD),'holes'),seD); %this is the same as an open, but fills in the holes first.
%         BW_smooth=imerode(imdilate(BW_filt,seD),seD); % this a close
%           operation
%         BW_smooth=imdilate(imerode(BW_filt,seD),seD); % this is an open
%           operation
    else
        BW_smooth=imbw;
    end
    BW_final = bwpropfilt(BW_smooth, 'Area', [sizefilt, numel(im)]);
end

