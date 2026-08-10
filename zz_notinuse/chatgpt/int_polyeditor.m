function int_polyeditor(imgStack)


    labeledPolys = [];    
    % Define labeled polygons
    labeledPolys(1).label = 'Constricted_BV';
    labeledPolys(1).polygon = polyshape([50 80 80 50], [30 30 50 50]); % Square
    labeledPolys(1).originalVertices = labeledPolys(1).polygon.Vertices;

    labeledPolys(2).label = 'Dialated_parenchymalborder';
    labeledPolys(2).polygon = polyshape([20 30 30 20], [20 20 30 30]); % Another square
    labeledPolys(2).originalVertices = labeledPolys(2).polygon.Vertices;

    labeledPolys(3).label = 'Surrouding_neuropile';
    labeledPolys(3).polygon = polyshape([40 50 50 40], [10 10 20 20]); % Third square
    labeledPolys(3).originalVertices = labeledPolys(3).polygon.Vertices;

    % Create figure
    fig = figure('Name', 'Interactive Polyshape Editor', 'NumberTitle', 'off', ...
        'KeyPressFcn', @keyPress, 'WindowButtonDownFcn', @mouseClick, ...
        'WindowButtonUpFcn', @releaseVertex, 'WindowButtonMotionFcn', @moveVertex);
    
    hold on;
    ax = gca;
    currentSlice = 1;
    imgHandle = imshow(imgStack(:,:,currentSlice),[],'Parent',ax);
    hold on;
    axis equal;

    % UI Buttons
    selectBtn = uicontrol('Style', 'togglebutton', 'String', 'Select Object', ...
        'Position', [20, 20, 100, 40], 'Callback', @toggleMode, 'Value', 1);
    
    moveBtn = uicontrol('Style', 'togglebutton', 'String', 'Move', ...
        'Position', [130, 20, 80, 40], 'Callback', @toggleMode);
    
    editBtn = uicontrol('Style', 'togglebutton', 'String', 'Edit Vertices', ...
        'Position', [220, 20, 100, 40], 'Callback', @toggleMode);
    
    addDelBtn = uicontrol('Style', 'togglebutton', 'String', 'Add/Delete Vertices', ...
        'Position', [330, 20, 140, 40], 'Callback', @toggleMode);

    % Variables
    selectedIdx = []; % Selected polygon index
    selectedVertex = []; % Selected vertex index
    mode = "select"; % Default mode
    updateDisplay();

    % Additional UI figure for the table only
    tableFig = uifigure('Name','Polygon Labels','Position',[800,100,150,150]);
    labelTable = uitable('Parent',tableFig,...
        'Data',{labeledPolys.label}',...
        'Position',[10,10,130,120],...
        'ColumnName',{'Labels'},...
        'RowName',[],...
        'ColumnWidth',{120},...
        'Tag','polyTable',...
        'CellSelectionCallback',@tableSelectionCallback);

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

    % Toggle Mode Function
    function toggleMode(src, ~)
        % Ensure only one mode is active at a time
        selectBtn.Value = 0;
        moveBtn.Value = 0;
        editBtn.Value = 0;
        addDelBtn.Value = 0;
        src.Value = 1;

        % Set mode based on button pressed
        if src == selectBtn, mode = "select"; end
        if src == moveBtn, mode = "move"; end
        if src == editBtn, mode = "edit"; end
        if src == addDelBtn, mode = "adddelete"; end

        % Clear selections when mode changes
        selectedVertex = [];
        updateDisplay();
    end

    % Mouse Click Callback (Handles Selection, Adding, and Deletion)
    function mouseClick(~, event)
        point = get(ax, 'CurrentPoint');
        xClick = point(1,1);
        yClick = point(1,2);

        if mode == "select" % Select a polygon
            selectedIdx = [];
            for i = 1:numel(labeledPolys)
                if isinterior(labeledPolys(i).polygon, xClick, yClick)
                    selectedIdx = i;
                    updateDisplay();
                    updateTableSelection(); % <-- Insert this line
                    return;
                end
            end

        elseif mode == "edit" % Select vertex for editing
            if isempty(selectedIdx), return; end
            vx = labeledPolys(selectedIdx).polygon.Vertices(:,1);
            vy = labeledPolys(selectedIdx).polygon.Vertices(:,2);
            dists = sqrt((vx - xClick).^2 + (vy - yClick).^2);
            [minDist, minIdx] = min(dists);
            if minDist < 5
                selectedVertex = minIdx;
            end

        elseif mode == "adddelete"
            if isempty(selectedIdx), return; end
            vx = labeledPolys(selectedIdx).polygon.Vertices(:,1);
            vy = labeledPolys(selectedIdx).polygon.Vertices(:,2);

            if strcmp(event.Source.SelectionType, 'alt') % Right-click to delete
                dists = sqrt((vx - xClick).^2 + (vy - yClick).^2);
                [minDist, minIdx] = min(dists);

                if minDist < 0.2 && numel(vx) > 3 % Prevent deleting last 3 vertices
                    vx(minIdx) = [];
                    vy(minIdx) = [];
                    labeledPolys(selectedIdx).polygon = polyshape(vx, vy);
                    updateDisplay();
                end
            else % Left-click to add a vertex
                minDist = inf;
                insertIdx = 1;
                for i = 1:numel(vx)
                    j = mod(i, numel(vx)) + 1;
                    midX = (vx(i) + vx(j)) / 2;
                    midY = (vy(i) + vy(j)) / 2;
                    dist = sqrt((midX - xClick)^2 + (midY - yClick)^2);
                    if dist < minDist
                        minDist = dist;
                        insertIdx = i + 1;
                    end
                end

                vx = [vx(1:insertIdx-1); xClick; vx(insertIdx:end)];
                vy = [vy(1:insertIdx-1); yClick; vy(insertIdx:end)];
                labeledPolys(selectedIdx).polygon = polyshape(vx, vy);
                updateDisplay();
            end
        end
    end
    function moveVertex(~, ~)
        if mode == "edit" && ~isempty(selectedVertex) && ~isempty(selectedIdx)
            point = get(ax, 'CurrentPoint');
            newX = point(1,1);
            newY = point(1,2);
    
            % Current vertices
            currentPoly = labeledPolys(selectedIdx).polygon;
            originalVertexCount = size(currentPoly.Vertices, 1);
    
            vx = currentPoly.Vertices(:,1);
            vy = currentPoly.Vertices(:,2);
    
            % Temporarily move vertex
            tempVx = vx;
            tempVy = vy;
            tempVx(selectedVertex) = newX;
            tempVy(selectedVertex) = newY;
    
            % Attempt to create new polygon
            tempPoly = polyshape(tempVx, tempVy);
    
            % Check if MATLAB altered the polygon
            if tempPoly.NumRegions == 0 || size(tempPoly.Vertices,1) ~= originalVertexCount
                disp('Invalid move detected. Move rejected.');
                return; % Reject the move immediately
            end
    
            % Check for minimum vertex proximity
            minAllowedDist = 1;
            otherVertices = [tempVx, tempVy];
            otherVertices(selectedVertex,:) = [];
            distances = sqrt((otherVertices(:,1)-newX).^2 + (otherVertices(:,2)-newY).^2);
    
            if any(distances < minAllowedDist)
                disp('Vertex too close to another vertex. Move canceled.');
                return; % Reject the move immediately
            end
    
            % If valid, apply move
            labeledPolys(selectedIdx).polygon = tempPoly;
            updateDisplay();
        end
    end

    % Release Vertex Callback
    function releaseVertex(~, ~)
        selectedVertex = [];
    end

    % Update Display
    function updateDisplay()
        delete(findobj(ax, 'Type', 'Polygon'));
        delete(findobj(ax, 'Type', 'Scatter'));
        plotPolygons();
    end
    % Table selection callback (ADD THIS INSIDE int_polyeditor function)
    function tableSelectionCallback(src, event)
        if ~isempty(event.Indices)
            selectedIdx = event.Indices(1); % get selected polygon index
            selectedVertex = []; % clear vertex selection if any
            updateDisplay();     % refresh plot to highlight selected polygon
        end
    end

    function keyPress(~, event)
            if isempty(selectedIdx) || mode ~= "move"
                return;
            end
    
            vx = labeledPolys(selectedIdx).polygon.Vertices(:,1);
            vy = labeledPolys(selectedIdx).polygon.Vertices(:,2);
            step = 0.1;
    
            switch event.Key
                case 'uparrow'
                    vy = vy + step;
                case 'downarrow'
                    vy = vy - step;
                case 'rightarrow'
                    vx = vx + step;
                case 'leftarrow'
                    vx = vx - step;
            end
    
            labeledPolys(selectedIdx).polygon = polyshape(vx, vy);
            updateDisplay();
    end

    % Plot Polygons with Colored Vertices
    function plotPolygons()
        hold on;
        for i = 1:numel(labeledPolys)
            plot(labeledPolys(i).polygon, 'FaceColor', 'c', 'FaceAlpha', 0.1);
        end

        if ~isempty(selectedIdx)
            plot(labeledPolys(selectedIdx).polygon, 'FaceColor', 'm', 'FaceAlpha', 0.2);
            
            vx = labeledPolys(selectedIdx).polygon.Vertices(:,1);
            vy = labeledPolys(selectedIdx).polygon.Vertices(:,2);
            originalVx = labeledPolys(selectedIdx).originalVertices(:,1);
            originalVy = labeledPolys(selectedIdx).originalVertices(:,2);

            isNewVertex = ~ismember([vx vy], [originalVx originalVy], 'rows');

            scatter(vx(~isNewVertex), vy(~isNewVertex), 70, 'k', 'filled', 'MarkerEdgeColor', 'y'); % Black for original
            scatter(vx(isNewVertex), vy(isNewVertex), 70, 'g', 'filled', 'MarkerEdgeColor', 'y'); % Green for new
            if ~isempty(selectedVertex)
                scatter(vx(selectedVertex), vy(selectedVertex), 100, 'b', 'filled', 'MarkerEdgeColor', 'y'); % Selected in blue
            end
        end
    end

function updateTableSelection()
    if isvalid(labelTable)
        labelTable.Data = {labeledPolys.label}';
        if ~isempty(selectedIdx)
            labelTable.Selection = [selectedIdx, 1];

            % Manually focus the table figure briefly
            figure(labelTable.Parent); % Force figure to foreground briefly
            pause(0.05); % Short pause to mimic interaction
        else
            labelTable.Selection = [];
        end
    end
end

end
