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

    fig   = uifigure('Name','Stack Explorer','Position',[100 100 900 650]);
    
    outer = uigridlayout(fig,[2 1], ...
            'RowHeight',{'1x','fit'}, ...   % image gets the leftover space
            'RowSpacing',6,'Padding',[8 8 8 8]);
    
    imgPanel      = uipanel(outer,'Title','Slice Viewer');
    controlPanel  = uipanel(outer,'Title','Console');

    g = uigridlayout(controlPanel,[3 1],'RowHeight',{'fit','fit'},'RowSpacing',8,'Padding',[10 10 10 10]);
    % ----- Row 1 : intensity sub-panel --------------------------------
    intensityPanel = uipanel(g,'Title','Intensity');
    intensityPanel.Layout.Row    = 1;   % instead of 'Layout',struct(...)
    intensityPanel.Layout.Column = 1;
    gi = uigridlayout(intensityPanel,[2 2],'ColumnWidth',{100,'1x'},'RowSpacing',5);
    
    uilabel(gi,'Text','Min Intensity:','HorizontalAlignment','left');
    minIntensitySlider = uislider(gi,'Limits',[0 255],'Value',0,...
        'ValueChangedFcn',@(s,~)updateIntensity);
    
    uilabel(gi,'Text','Max Intensity:','HorizontalAlignment','left');
    maxIntensitySlider = uislider(gi,'Limits',[0 255],'Value',255,...
        'ValueChangedFcn',@(s,~)updateIntensity);
    
    % ----- Row 2 : button sub-panel -----------------------------------
    buttonPanel    = uipanel(g,'Title','Actions');
    buttonPanel.Layout.Row    = 3;
    buttonPanel.Layout.Column = 1;    
    gb = uigridlayout(buttonPanel,[1 4], ...
                  'ColumnWidth',{'1x','1x', 120, '1x'}, ... % Reset | Confirm | slider | numeric
                  'ColumnSpacing',10);
    
    uibutton(gb,'Text','Reset',   'ButtonPushedFcn',@(~,~)resetROI);
    uibutton(gb,'Text','Confirm', 'ButtonPushedFcn',@(~,~)uiresume(fig));
    
    % line-width widgets only if we’re drawing a line ROI
    if strcmp(roi_type,'line')
        % 3 a) numeric edit field so user sees the value
        wEdit = uieditfield(gb,'numeric','Limits',[0 30],'Value',10, ...
            'ValueChangedFcn',@syncEditToSlider);
    
        % 3 b) slider – one per GUI, NOT per redraw
        wSlider = uislider(gb,'Limits',[0 30],'Value',10, ...
            'ValueChangedFcn',@updateLineWidth);
    end
    
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

    % Function to dynamically update the intensity range
    function updateIntensity()
        % Update the intensity range for the axes
        hAxes.CLim = [minIntensitySlider.Value, maxIntensitySlider.Value];
    end

  % ---------- shared placeholders (put HIGH in the file, right after you
    % build the UI but BEFORE you call drawROI) ------------------------------
    theROI      = gobjects(0);                % interactive ROI (empty for now)
    overlayLine = gobjects(0);                % translucent overlay
    moveLstn    = event.listener.empty;       % listener handle
    resetFlag   = false;

    % Initialize ROI
    theROI = drawROI();

    % Main loop to manage ROI interaction
    while true
        uiwait(fig); % initial interaction with haxes occurs here, button click will run the cycle
        if ~isvalid(fig)
            break; % Exit if the figure is closed
        end
        if resetFlag % 
            resetFlag = false; % Reset the flag and continue drawing
            continue; % go back to uiwait
        end
        break; % Exit loop on confirm resetFlag
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
            elseif strcmp(roi_type,'line')
                roi = drawline(hAxes,'Color','m', ...
                               'LineWidth',wSlider.Value,'LabelVisible','off');
        
                % translucent overlay
                overlayLine = line(roi.Position(:,1), roi.Position(:,2), ...
                              'Parent',hAxes, ...
                              'Color',[1 0 1 0.3], ...
                              'LineWidth',roi.LineWidth, ...
                              'HitTest','off');
        
            % ----- keep overlay in sync while the ROI moves --------------------
            moveLstn(1) = addlistener(roi,'MovingROI', ...
                             @(~,evt) syncOverlay(evt.CurrentPosition));
            
            moveLstn(2) = addlistener(roi,'ROIMoved',  ...
                             @(~,evt) syncOverlay(evt.CurrentPosition));
            else
                error('Unsupported ROI type: %s. Use "rectangle" or "polygon".', roi_type);
            end
        end
    end


    % Function to reset ROI
    function resetROI
        % remove overlay and its listener first
        if isvalid(moveLstn),    delete(moveLstn);    end
        if isvalid(overlayLine), delete(overlayLine); end
        moveLstn    = event.listener.empty;
        overlayLine = gobjects(0);
    
        % remove the ROI itself
        if isvalid(theROI), delete(theROI); end
    
        % create fresh ROI / overlay
        theROI  = drawROI();
        resetFlag = true;
        uiresume(fig);
    end
    function updateLineWidth(src,~)
        if exist('wEdit','var') && isvalid(wEdit)
            wEdit.Value = src.Value;
        end
        if isvalid(theROI)
            theROI.LineWidth = max(src.Value,0.5);
        end
        if isvalid(overlayLine)
            overlayLine.LineWidth = max(src.Value,0.5);
        end
    end
function syncOverlay(pos)
    if isvalid(overlayLine)
        overlayLine.XData = pos(:,1);
        overlayLine.YData = pos(:,2);
        overlayLine.LineWidth = theROI.LineWidth;
    end
end
    function syncEditToSlider(src,~)
    if isvalid(wSlider)
        wSlider.Value = src.Value;
        updateLineWidth(wSlider);
    end
end

end