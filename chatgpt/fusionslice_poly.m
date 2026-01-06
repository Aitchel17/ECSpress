function fusionslice_poly(imgStack)
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
        'Position',[280 20 100 20],'String',sprintf('Slice: %d',currentSlice));

    % Button to launch imcontrast
    uicontrol('Parent',fig,'Style','pushbutton',...
        'Position',[400 20 150 30],'String','Adjust Contrast',...
        'Callback',@(src,event) imcontrast(ax));

    % Your labeled polygons here
    labeledPolys(1).label = 'Region A';
    labeledPolys(1).polygon = polyshape([10 50 50 10],[10 10 50 50]);

    labeledPolys(2).label = 'Region B';
    labeledPolys(2).polygon = polyshape([60 80 80 60],[60 60 80 80]);

    % Plot polygons initially
    polyHandles = gobjects(numel(labeledPolys),1);
    for i = 1:numel(labeledPolys)
        polyHandles(i) = plot(ax,labeledPolys(i).polygon,'FaceAlpha',0.4);
    end

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