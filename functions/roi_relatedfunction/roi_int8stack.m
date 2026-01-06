function [isRGB, stack] = roi_int8stack(stack)  
% Determine if the stack is RGB
    isRGB = (ndims(stack) == 4 && size(stack, 3) == 3);

    % Normalize the stack for display
    if ~isRGB
        % Grayscale stack normalization
        min_val = min(stack, [], 'all');
        max_val = max(stack, [], 'all');
        if isa(stack, 'double')
            stack = (stack - min_val) / (max_val - min_val) * 256;
            stack = uint8(stack);
        end
    else
        % RGB stack normalization
        for i = 1:3
            channel = stack(:, :, i, :);
            min_val = min(channel, [], 'all');
            max_val = max(channel, [], 'all');
            if isa(channel, 'double')
                stack(:, :, i, :) = (channel - min_val) / (max_val - min_val) * 256;
            end
        end
        stack = uint8(stack);
    end
