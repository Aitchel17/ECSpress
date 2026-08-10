function interactive_polyshape_editor()
    % Define labeled polygons
    labeledPolys(1).label = 'Region A';
    labeledPolys(1).polygon = polyshape([0 1 1 0], [0 0 1 1]); % Square
    labeledPolys(1).originalVertices = labeledPolys(1).polygon.Vertices;

    labeledPolys(2).label = 'Region B';
    labeledPolys(2).polygon = polyshape([2 3 3 2], [2 2 3 3]); % Another square
    labeledPolys(2).originalVertices = labeledPolys(2).polygon.Vertices;

    labeledPolys(3).label = 'Region C';
    labeledPolys(3).polygon = polyshape([4 5 5 4], [1 1 2 2]); % Third square
    labeledPolys(3).originalVertices = labeledPolys(3).polygon.Vertices;

    % Create figure
    fig = figure('Name', 'Interactive Polyshape Editor', 'NumberTitle', 'off', ...
        'KeyPressFcn', @keyPress, 'WindowButtonDownFcn', @mouseClick, ...
        'WindowButtonUpFcn', @releaseVertex, 'WindowButtonMotionFcn', @moveVertex);
    hold on;
    ax = gca;
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
                    return;
                end
            end

        elseif mode == "edit" % Select vertex for editing
            if isempty(selectedIdx), return; end
            vx = labeledPolys(selectedIdx).polygon.Vertices(:,1);
            vy = labeledPolys(selectedIdx).polygon.Vertices(:,2);
            dists = sqrt((vx - xClick).^2 + (vy - yClick).^2);
            [minDist, minIdx] = min(dists);

            if minDist < 0.2
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

    % Move Vertex Callback (Draggable)
    function moveVertex(~, ~)
        if mode == "edit" && ~isempty(selectedVertex) && ~isempty(selectedIdx)
            point = get(ax, 'CurrentPoint');
            vx = labeledPolys(selectedIdx).polygon.Vertices(:,1);
            vy = labeledPolys(selectedIdx).polygon.Vertices(:,2);
            vx(selectedVertex) = point(1,1);
            vy(selectedVertex) = point(1,2);

            labeledPolys(selectedIdx).polygon = polyshape(vx, vy);
            updateDisplay();
        end
    end

    % Release Vertex Callback
    function releaseVertex(~, ~)
        selectedVertex = [];
    end

    % Update Display
    function updateDisplay()
        cla; % Clear plot
        plotPolygons();
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
            plot(labeledPolys(i).polygon, 'FaceColor', 'c', 'FaceAlpha', 0.5);
        end

        if ~isempty(selectedIdx)
            plot(labeledPolys(selectedIdx).polygon, 'FaceColor', 'm', 'FaceAlpha', 0.7);
            
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
end
