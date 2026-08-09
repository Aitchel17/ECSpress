classdef imgstack_reader
%IMGSTACK_READER  VideoReader stand-in serving slices of an in-memory 3-D stack.
%   Wraps an [H x W x N] matrix so PIVlab functions that expect a video file path
%   or a VideoReader object (piv_FFTensemble and friends) can read frame by frame.
%   Frames are indexed along the third dimension.
%
% IN   stack   HxWxN image stack
% OUT  obj     reader; read(obj, idx) returns slice idx

    properties
        imgstack        % HxWxN image stack, frames along dim 3
    end

    methods
        function obj = imgstack_reader(stack)
        %IMGSTACK_READER  Wrap an image stack in a reader.
            obj.imgstack = stack;
        end

        function img = read(obj, idx)
        %READ  Return frame idx as a 2-D image.
            img = obj.imgstack(:,:,idx);
        end
    end
end
