function [thrstack] = pre_thresholding(inputstack)
%PRE_RESCALE_UINT16 Summary of this function goes here
%   Detailed explanation goes here
    % rescale
    inputstack = gpuArray(inputstack);
    [height, width, ~] = size(inputstack);
    threshold = min(median(inputstack(1+ceil(height/3):end-ceil(height/3),1+ceil(width/3):end-ceil(width/3),:),3),[],'all');
    thrstack = inputstack-threshold;
    thrstack(thrstack<0) = 0;
    thrstack = gather(thrstack);
end

