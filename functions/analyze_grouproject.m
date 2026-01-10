function [zstack, depricate_frame] = analyze_grouproject(zstack, number_averaging,mode)
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
    arguments
        zstack (:,:,:) {mustBeNumeric}    % 3D numeric array
        number_averaging (1,1) {mustBeInteger, mustBePositive} % positive integer
        mode (1,:) char {mustBeMember(mode, ["mean", "median"])} = "mean" % optional, default "mean"
    end

    % Ensure the number of frames is a multiple of the group size
    if number_averaging == 1
        disp('skip group averaging')
    else
        depricate_frame = mod(size(zstack,3),number_averaging);
        disp(['Group averaging with ' num2str(number_averaging) ' depricate ' num2str(depricate_frame)]);
        zstack = zstack(:,:,1:end-depricate_frame);
    
        % Reshape the stack to group the frames
        zstack = reshape(zstack, size(zstack, 1), size(zstack, 2), number_averaging, []);
        
        % Calculate the mean of each group of frames along the 3rd dimension
        if strcmp(mode,'median')
        zstack = squeeze(median(zstack, 3)); % Average along the 3rd dimension (grouped frames)
        else
        zstack = squeeze(mean(zstack, 3)); % Average along the 3rd dimension (grouped frames)
        end
end