function slice3d_polyeditor(imgStack)
      % imgStack: 3D matrix (rows x cols x slices)

    % Initialize Figure
    fig = figure('Name','Custom Slice Viewer with Polygon Editor',...
        'NumberTitle','off');

    % Set initial slice
    currentSlice = 1;

    % Display initial image slice
    ax = axes('Parent',fig);
    imgHandle = imshow(imgStack(:,:,currentSlice),[],'Parent',ax);
    hold(ax,'on');

    % Add slider to change slice
    sliceSlider = uicontrol('Parent',fig,'Style','slider',...
        'Position',[20 20 250 20],'Min',1,'Max',size(imgStack,3),...
        'Value',currentSlice,'SliderStep',[1/(size(imgStack,3)-1),0.1],...
        'Callback',@changeSlice,...
        'Interruptible','off','BusyAction','cancel');
    addlistener(sliceSlider, 'ContinuousValueChange', @changeSlice);

    % Display current slice number
    sliceLabel = uicontrol('Parent',fig,'Style','text',...
        'Position',[2380 20 100 20],'String',sprintf('Slice: %d',currentSlice));

    % Button to launch imcontrast
    uicontrol('Parent',fig,'Style','pushbutton',...
        'Position',[400 20 150 20],'String','Adjust Contrast',...
        'Callback',@(src,event) imcontrast(ax));
    
    % Callback for slice slider
    function changeSlice(src,~)
        currentSlice = round(src.Value);
        updateImage();
        sliceLabel.String = sprintf('Slice: %d',currentSlice);
    end

    % Update displayed image based on slice
    function updateImage()
        imgHandle.CData = imgStack(:,:,currentSlice);
    end
end