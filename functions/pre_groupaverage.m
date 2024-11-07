function [result] = pre_groupaverage(zstack, number_averaging)
%FRAME_AVERAGING Summary of this function goes here
%   Detailed explanation goes here


zstack = zstack(:,:,1:size(zstack,3)-mod(size(zstack,3),number_averaging));

result = squeeze( ...
    mean( ...
    reshape(zstack,size(zstack,1),size(zstack,2),number_averaging,[]),3) ... % reshape by 
       ... % calculate mean by axis
    ); ... % Remove averaging axis
