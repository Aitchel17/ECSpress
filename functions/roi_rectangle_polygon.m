function [vertices, ref_slice] = roi_rectangle_polygon(stack, roi_type, preexistingvertices)
    % ROI_RECTANGLE_POLYGON - Handles grayscale and RGB stacks for ROI drawing
    %
    % INPUTS:
    %   stack              - 3D grayscale stack or 4D RGB stack.
    %   roi_type           - 'rectangle' or 'polygon'.
    %   preexistingvertices (optional) - Nx2 matrix of vertices for initialization.
    %
    % OUTPUTS:
    %   vertices   - Nx2 matrix of vertices of the ROI.
    %   ref_slice  - Slice number where ROI was drawn.

    % Default for optional argument
    if nargin < 3
        preexistingvertices = [];
    end

    % Determine if the stack is RGB
    isRGB = (ndims(stack) == 4 && size(stack, 3) == 3);

    % Normalize the stack for display
    if ~isRGB
        % Grayscale stack normalization
        min_val = min(stack, [], 'all');
        max_val = max(stack, [], 'all');
        if isa(stack, 'double')
            stack = (stack - min_val) / (max_val - min_val) * 65535;
            stack = uint16(stack);
        end
    else
        % RGB stack normalization
        for i = 1:3
            channel = stack(:, :, i, :);
            min_val = min(channel, [], 'all');
            max_val = max(channel, [], 'all');
            if isa(channel, 'double')
                stack(:, :, i, :) = (channel - min_val) / (max_val - min_val) * 65535;
            end
        end
        stack = uint16(stack);
    end

    % Create the main figure
    fig = uifigure('Name', 'Stack Explorer', 'Position', [100, 100, 600, 400]);

    % Create panels for controls and image display
    imgPanel = uipanel(fig, 'Title', 'Slice Viewer', 'Position', [20, 120, 560, 260]);
    controlPanel = uipanel(fig, 'Title', 'Console', 'Position', [20, 20, 560, 100]);

    % Display the stack using sliceViewer
    if isRGB
        % Extract the first slice for RGB visualization
        sliceData = squeeze(stack(:, :, :, 1)); % First slice of the stack
        hStack = sliceViewer(sliceData, 'Parent', imgPanel);
    else
        % Grayscale stack
        hStack = sliceViewer(stack, 'Parent', imgPanel);
    end

    % Extract the underlying axes object from the sliceViewer
    hAxes = getAxesHandle(hStack);

    % Add a label and sliders for intensity adjustment
    uilabel(controlPanel, 'Text', 'Min Intensity:', 'Position', [20, 60, 100, 20]);
    minIntensitySlider = uislider(controlPanel, ...
        'Position', [130, 65, 400, 3], ...
        'Limits', [0, 65535], ...
        'Value', 0, ...
        'ValueChangedFcn', @(src, ~) updateIntensity());

    uilabel(controlPanel, 'Text', 'Max Intensity:', 'Position', [20, 30, 100, 20]);
    maxIntensitySlider = uislider(controlPanel, ...
        'Position', [130, 35, 400, 3], ...
        'Limits', [0, 65535], ...
        'Value', 65535, ...
        'ValueChangedFcn', @(src, ~) updateIntensity());

    % Function to dynamically update the intensity range
    function updateIntensity()
        % Update the intensity range for the axes
        hAxes.CLim = [minIntensitySlider.Value, maxIntensitySlider.Value];
    end

    % Add instructions label
    uilabel(controlPanel, ...
        'Text', 'Adjust intensity and draw or modify the ROI. Then click Confirm.', ...
        'Position', [20, 20, 500, 20], ...
        'HorizontalAlignment', 'left');

    % Add a Confirm button
    uibutton(controlPanel, ...
        'Text', 'Confirm', ...
        'Position', [480, 10, 70, 30], ...
        'ButtonPushedFcn', @(src, event) uiresume(fig)); % Resume execution when clicked

    % Add Reset button
    uibutton(controlPanel, ...
        'Text', 'Reset', ...
        'Position', [400, 10, 70, 30], ...
        'ButtonPushedFcn', @(~, ~) resetROI());

    % Initialize ROI
    theROI = drawROI();
    resetFlag = false;

    % Main loop to manage ROI interaction
    while true
        uiwait(fig);
        if ~isvalid(fig)
            break; % Exit if the figure is closed
        end
        if resetFlag
            resetFlag = false; % Reset the flag and continue drawing
            continue;
        end
        break; % Exit loop on confirm
    end

    % Extract vertices and slice number

    vertices = round(theROI.Position); % Polygon vertices

    % convert rectangle position to vertices
    if strcmp(roi_type, 'rectangle')
        x1 = vertices(1); 
        y1 = vertices(2);
        x2 = x1 + vertices(3); 
        y2 = y1 + vertices(4);
        vertices = [x1, y1; x2, y1; x2, y2; x1, y2];
    end

    % Extract the current slice from the sliceViewer
    ref_slice = hStack.SliceNumber;

    % Close the figure
    close(fig);

    % Function to draw ROI based on type
    function roi = drawROI()
        if ~isempty(preexistingvertices)
            % Initialize with preexisting vertices
            if strcmp(roi_type, 'rectangle')
                roi = drawrectangle(hAxes, 'Position', preexistingvertices);
            elseif strcmp(roi_type, 'polygon')
                roi = drawpolygon(hAxes, 'Position', preexistingvertices);
            else
                error('Unsupported ROI type: %s. Use "rectangle" or "polygon".', roi_type);
            end
        else
            % Create a new ROI
            if strcmp(roi_type, 'rectangle')
                roi = drawrectangle(hAxes);
            elseif strcmp(roi_type, 'polygon')
                roi = drawpolygon(hAxes);
            else
                error('Unsupported ROI type: %s. Use "rectangle" or "polygon".', roi_type);
            end
        end
    end

    % Function to reset ROI
    function resetROI()
        if isvalid(theROI)
            delete(theROI); % Delete existing ROI
        end
        theROI = drawROI(); % Allow user to redraw
        resetFlag = true; % Set the reset flag
        uiresume(fig); % Resume the UI loop
    end
end
