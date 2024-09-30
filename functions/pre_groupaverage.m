function [outputArg1,outputArg2] = pre_groupaverage(zstack, number_averaging)
%FRAME_AVERAGING Summary of this function goes here
%   Detailed explanation goes here

result = squeeze( ...
    mean( ...
    reshape(zstack,size(data.prestack,1),size(tmp.medstack,2),number_averaging,[]),3) ... % reshape by 
       ... % calculate mean by axis
    ); ... % Remove averaging axis
figure()
sliceViewer(tmp.gpave_stack)

figure()
imshow(mean(data.prestack,3), [0, 1000])
end

