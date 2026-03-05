function [state] = util_checkstack(stack, window_title)
arguments
    stack
    window_title = 'Stack Explorer'
end
% just for inspection purpose

% Check dimensions to see if it's RGB
sz = size(stack);
if ndims(stack) == 4
    if sz(3) == 3
        is_rgb = true;
        num_slices = sz(4);
    else
        error('4D stack must have 3 channels in the 3rd dimension');
    end
elseif ndims(stack) == 3
    is_rgb = false;
    num_slices = sz(3);
elseif ismatrix(stack)
    is_rgb = false;
    num_slices = 1;
else
    error('Unsupported stack dimensions');
end

% Normalize stack
stack = double(stack);
if is_rgb
    for c = 1:3
        c_stack = stack(:,:,c,:);
        c_min = double(min(c_stack,[],'all'));
        c_max = double(max(c_stack,[],'all'));
        stack(:,:,c,:) = (c_stack - c_min) / max(1e-9, c_max - c_min) * 65535;
    end
else
    min_val = double(min(stack,[],'all'));
    max_val = double(max(stack,[],'all'));
    stack = (stack - min_val) / max(1e-9, max_val - min_val) * 65535; 
end
stack = uint16(stack);

state = false;

% Create the main figure
fig = uifigure('Name', window_title, 'Position', [100, 100, 800, 600]);

% Create panels for controls and image display
sliderPanel = uipanel(fig, 'Title', 'Channels', 'Position', [20, 120, 170, 460]);
imgPanel = uipanel(fig, 'Title', 'Slice Viewer', 'Position', [200, 120, 580, 460]);
controlPanel = uipanel(fig, 'Title', 'Console', 'Position', [20, 20, 760, 90]);

% Axes for displaying the image
ax = uiaxes(imgPanel, 'Position', [10 10 560 410]);
ax.XTick = [];
ax.YTick = [];

current_slice = 1;
curr_img = [];

% Set up sliders in Channel panel
if is_rgb
    climits = repmat([0 65535], 3, 1);
    colors = {'Red', 'Green', 'Blue'};
    for c = 1:3
        y_pos = 400 - (c-1)*100;
        uilabel(sliderPanel, 'Text', [colors{c} ' Range:'], 'Position', [10, y_pos, 150, 20]);
        uislider(sliderPanel, 'range', ...
            'Position', [10, y_pos-20, 150, 3], ...
            'Limits', [0, 65535], ...
            'Value', [0, 65535], ...
            'MajorTicks', [], ...
            'ValueChangingFcn', @(~, event) update_rgb_range(c, event.Value), ...
            'ValueChangedFcn', @(src, ~) update_rgb_range(c, src.Value));
    end
else
    climits = [0 65535];
    uilabel(sliderPanel, 'Text', 'Intensity Range:', 'Position', [10, 400, 150, 20]);
    uislider(sliderPanel, 'range', ...
        'Position', [10, 380, 150, 3], ...
        'Limits', [0, 65535], ...
        'Value', [0, 65535], ...
        'MajorTicks', [], ...
        'ValueChangingFcn', @(~, event) update_gray_range(event.Value), ...
        'ValueChangedFcn', @(src, ~) update_gray_range(src.Value));
end

% Set up controls in Control panel
if num_slices > 1
    uilabel(controlPanel, 'Text', 'Slice:', 'Position', [20, 50, 50, 20]);
    uislider(controlPanel, ...
        'Position', [80, 60, 450, 3], ...
        'Limits', [1, num_slices], ...
        'Value', 1, ...
        'MajorTicks', [], ...
        'ValueChangingFcn', @(~, event) change_slice(round(event.Value)), ...
        'ValueChangedFcn', @(src, ~) change_slice(round(src.Value)));
end

% Add instructions label
uilabel(controlPanel, ...
    'Text', 'Adjust intensity and draw a rectangle around ROI (if needed). Then click Confirm.', ...
    'Position', [20, 15, 500, 20], ...
    'HorizontalAlignment', 'left');

% Add Confirm and Reject buttons
uibutton(controlPanel, ...
    'Text', 'Confirm', ...
    'Position', [660, 50, 80, 30], ...
    'ButtonPushedFcn', @(~, ~) confirm());

uibutton(controlPanel, ...
    'Text', 'Reject', ...
    'Position', [660, 10, 80, 30], ...
    'ButtonPushedFcn', @(~, ~) uiresume(fig));

% Draw the initial image
draw_image();

% Wait until uiresume to be run by confirm() or reject
uiwait(fig);

% Close figure if it still exists
if isvalid(fig)
    close(fig);
end

% == Nested Functions ==

    function draw_image()
        if is_rgb
            img = double(stack(:,:,:,current_slice));
            for c = 1:3
                c_min = climits(c,1);
                c_max = max(climits(c,2), c_min + 1);
                c_img = img(:,:,c);
                c_img = (c_img - c_min) / (c_max - c_min);
                c_img(c_img < 0) = 0;
                c_img(c_img > 1) = 1;
                img(:,:,c) = c_img;
            end
            
            if isempty(curr_img) || ~isvalid(curr_img)
                curr_img = imagesc(ax, img);
            else
                curr_img.CData = img;
            end
        else
            img = stack(:,:,current_slice);
            if isempty(curr_img) || ~isvalid(curr_img)
                curr_img = imagesc(ax, img);
                colormap(ax, 'gray');
                ax.CLim = [climits(1, 1), max(climits(1, 2), climits(1, 1) + 1)];
            else
                curr_img.CData = img;
            end
        end
        axis(ax, 'image');
        title(ax, sprintf('Slice: %d / %d', current_slice, num_slices));
    end

    function update_rgb_range(c, range)
        climits(c, :) = range;
        draw_image();
    end

    function update_gray_range(range)
        climits(1, :) = range;
        if isvalid(ax)
            ax.CLim = [climits(1, 1), max(climits(1, 2), climits(1, 1) + 1)];
        end
    end

    function change_slice(val)
        current_slice = val;
        draw_image();
    end

    function confirm()
        state = current_slice; % Return the confirmed slice
        uiresume(fig);
    end

end
