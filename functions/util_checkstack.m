function [state] = util_checkstack(stack, window_title)
arguments
    stack
    window_title = 'Stack Explorer'
end
% just for inspection purpose

% Create the main figure
stack = double(stack);
min_val = min(stack,[],'all');
max_val = max(stack,[],'all');
stack = (stack - min_val) / (max_val - min_val) * 65535; % Normalize and stretch to min and max pixels.
stack = uint16(stack);
state = false;

fig = uifigure('Name', window_title, 'Position', [100, 100, 600, 400]);

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
    'ButtonPushedFcn', @(src, event) confirm()); % Resume execution when clicked

% Add Reset button
uibutton(controlPanel, ...
    'Text', 'Reject', ...
    'Position', [400, 10, 70, 30], ...
    'ButtonPushedFcn', @(~, ~) uiresume(fig));
uiwait(fig);

% Close the figure
close(fig);

% Function to update intensity range dynamically
    function updatefig(hAxes, range)
        hAxes.CLim = range; % Adjust display range
    end


% Function to reset ROI
    function confirm()
        state = hStack.SliceNumber; % Set the reset flag
        uiresume(fig); % Resume
    end
end
