function [vertices,ref_slice] = roi_rectangle(stack)
    % Create the main figure

    min_val = min(stack,[],'all');
    max_val = max(stack,[],'all');
    stack = (stack - min_val) / (max_val - min_val) * 65535;
    stack = uint16(stack);

    fig = uifigure('Name', 'Stack Explorer', 'Position', [100, 100, 600, 400]);
    
    % Create panels for controls and image display
    imgPanel = uipanel(fig, 'Title', 'Slice Viewer', 'Position', [20, 120, 560, 260]);
    controlPanel = uipanel(fig, 'Title', 'Console', 'Position', [20, 20, 560, 100]);
    
    % Display the stack using sliceViewer
    hStack = sliceViewer(stack, 'Parent', imgPanel);
    
    % Extract the underlying axes object from the sliceViewer
    hAxes = getAxesHandle(hStack);
    
    % Add a label for the intensity range slider
    uilabel(controlPanel, 'Text', 'Intensity Range:', 'Position', [20, 60, 100, 20]);
    
    % Add a range slider for adjusting intensity range
    intensitySlider = uislider(controlPanel, 'range',...
        'Position', [130, 65, 400, 3], ...
        'Limits', [0, 65535], ...
        'Value', [10000, 30000], ...
        'MajorTicks', [], ...
        'Orientation', 'horizontal', ...
        'ValueChangedFcn', @(src, event) updatefig(hAxes, src.Value));
    
    % Add instructions label
    uilabel(controlPanel, ...
        'Text', 'Adjust intensity and draw a rectangle around ROI. Then click Confirm.', ...
        'Position', [20, 20, 500, 20], ...
        'HorizontalAlignment', 'left');
    
    % Add a Confirm button
    uibutton(controlPanel, ...
        'Text', 'Confirm', ...
        'Position', [480, 10, 70, 30], ...
        'ButtonPushedFcn', @(src, event) uiresume(fig)); % Resume execution when clicked
    
    % Draw rectangle on the axes and wait for user confirmation
    instructionLabel = uilabel(controlPanel, ...
        'Text', 'Draw a rectangle around the ROI.', ...
        'Position', [20, 40, 400, 20], ...
        'HorizontalAlignment', 'left');
    thebox = drawrectangle(hAxes); % Allow the user to draw
    
    % Wait for confirmation
    uiwait(fig);
    
    % Get the vertices of the rectangle
    vertices = round(thebox.Vertices);
    % Extract the current slice from the sliceViewer
    ref_slice = hStack.SliceNumber;
    % Close the figure
    close(fig);

    % Function to update intensity range dynamically
    function updatefig(hAxes, range)
        hAxes.CLim = range; % Adjust the display range based on slider values
    end
end

