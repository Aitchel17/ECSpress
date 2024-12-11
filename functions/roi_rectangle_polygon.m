function [vertices,ref_slice] = roi_rectangle_polygon(stack,roi_type)
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
        'Value', [0, 65535], ...
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
    if strcmp(roi_type, 'rectangle')
        vertices = round(theROI.Vertices); % Rectangle vertices
    else
        vertices = round(theROI.Position); % Polygon vertices
    end

    % Extract the current slice from the sliceViewer
    ref_slice = hStack.SliceNumber;
    % Close the figure
    close(fig);

    % Function to update intensity range dynamically
    function updatefig(hAxes, range)
        hAxes.CLim = range; % Adjust display range
    end

    % Function to draw ROI based on type
    function roi = drawROI()
        if strcmp(roi_type, 'rectangle')
            roi = drawrectangle(hAxes); % Draw rectangle
        elseif strcmp(roi_type, 'polygon')
            roi = drawpolygon(hAxes); % Draw polygon
        else
            error('Unsupported ROI type: %s. Use "rectangle" or "polygon".', roi_type);
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

