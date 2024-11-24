function [result] = pre_groupaverage(zstack, number_averaging)
    % This function performs group averaging of a 3D z-stack.
    % It averages consecutive frames in the z-stack.
    % The frames are grouped into sets of `number_averaging` frames, 
    % and the average is calculated for each group.
    % Any frames that don't complete a full group are discarded.
    %
    % Inputs:
    %   zstack         - A 3D matrix (height x width x num_frames).
    %   number_averaging - The number of frames to average together.
    %
    % Outputs:
    %   result         - The resulting 3D matrix with the averaged frames.
    
    % Ensure the number of frames is a multiple of the group size
    zstack = zstack(:,:,1:size(zstack,3)-mod(size(zstack,3),number_averaging));

    % Reshape the stack to group the frames
    reshaped_stack = reshape(zstack, size(zstack, 1), size(zstack, 2), number_averaging, []);
    
    % Calculate the mean of each group of frames along the 3rd dimension
    result = squeeze(mean(reshaped_stack, 3)); % Average along the 3rd dimension (grouped frames)
end