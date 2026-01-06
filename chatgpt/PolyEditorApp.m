classdef PolyEditorApp < handle
    properties
        Fig
        Ax
        ImgHandle
        ImgStack
        CurrentSlice = 1
        
        LabeledPolys
        Mode = "select"

        % UI Elements
        ui = struct()
    end

    methods
        function app = PolyEditorApp(imgStack)
            app.ImgStack = imgStack;
            app.initializePolygons();
            app.createUI();
            app.updateDisplay();
        end

        function initializePolygons(app)
            app.ui.SelectedIdx = [];

            app.LabeledPolys(1).label = 'Constricted_BV';
            app.LabeledPolys(1).polygon = polyshape([100 160 160 100], [30 30 50 50]);
            app.LabeledPolys(1).originalVertices = app.LabeledPolys(1).polygon.Vertices;

            app.LabeledPolys(2).label = 'Dialated_parenchymalborder';
            app.LabeledPolys(2).polygon = polyshape([20 30 30 20], [20 20 30 30]);
            app.LabeledPolys(2).originalVertices = app.LabeledPolys(2).polygon.Vertices;
            

            % Add additional polygons similarly...
        end

        function createUI(app)
            % Create figure
            app.Fig = figure('Name', 'Interactive Polyshape Editor', ...
                'NumberTitle', 'off', ...
                'Position',[100, 100, 800, 650], ...
                'KeyPressFcn', @(src,event)app.keyPress(event), ...
                'WindowButtonDownFcn', @(src,event)app.mouseClick(event), ...
                'WindowButtonUpFcn', @(src,event)app.releaseVertex(), ...
                'WindowButtonMotionFcn', @(src,event)app.moveVertex());
        
            % Create axes (make sure it doesn't overlap buttons)
            app.Ax = axes('Parent', app.Fig, 'Units','pixels', 'Position',[50, 120, 700, 500]);
            hold(app.Ax, 'on');
            app.ImgHandle = imshow(app.ImgStack(:,:,app.CurrentSlice), [], 'Parent', app.Ax);
        
            % UI Buttons (placed below axes)
            buttonY = 50; % Y-position for buttons clearly below axes
        
            app.ui.SelectBtn = uicontrol('Style', 'togglebutton', 'String', 'Select Object', ...
                'Position', [50, buttonY, 100, 40], 'Callback', @(src,event)app.toggleMode("select"));
        
            app.ui.MoveBtn = uicontrol('Style', 'togglebutton', 'String', 'Move', ...
                'Position', [160, buttonY, 80, 40], 'Callback', @(src,event)app.toggleMode("move"));
        
            app.ui.EditBtn = uicontrol('Style', 'togglebutton', 'String', 'Edit Vertices', ...
                'Position', [250, buttonY, 100, 40], 'Callback', @(src,event)app.toggleMode("edit"));
        
            app.ui.AddDelBtn = uicontrol('Style', 'togglebutton', 'String', 'Add/Delete Vertices', ...
                'Position', [360, buttonY, 140, 40], 'Callback', @(src,event)app.toggleMode("adddelete"));
        
            % Additional UI elements (slider, labels) placed similarly below the axes
            sliderY = 20;
            app.ui.SliceSlider = uicontrol('Parent', app.Fig, 'Style', 'slider', ...
                'Position', [520, sliderY, 200, 20], ...
                'Min', 1, 'Max', size(app.ImgStack,3), ...
                'Value', app.CurrentSlice, ...
                'SliderStep', [1/(size(app.ImgStack,3)-1), 0.1], ...
                'Callback', @(src,event)app.changeSlice(src));
        
            % Slice Label
            app.ui.SliceLabel = uicontrol('Parent', app.Fig, 'Style', 'text', ...
                'Position', [730, sliderY, 60, 20], ...
                'String', sprintf('Slice: %d', app.CurrentSlice));
        end


        function toggleMode(app, newMode)
                app.Mode = newMode;
                app.ui.SelectBtn.Value = strcmp(newMode, "select");
                app.ui.MoveBtn.Value = strcmp(newMode, "move");
                app.ui.EditBtn.Value = strcmp(newMode, "edit");
                app.ui.AddDelBtn.Value = strcmp(newMode, "adddelete");
                app.ui.SelectedVertex = [];
                app.updateDisplay();
            end
    
        function mouseClick(app, event)
                point = get(app.Ax, 'CurrentPoint');
                        xClick = point(1,1);
            yClick = point(1,2);
    
            if mode == "select" % Select a polygon
                app.ui.LabelTabel = [];
                for i = 1:numel(app.LabeledPolys)
                    if isinterior(app.LabeledPolys(i).polygon, xClick, yClick)
                        app.ui.LabelTabel = i;
                        updateDisplay();
                        updateTableSelection(); % <-- Insert this line
                        return;
                    end
                end
    
            elseif mode == "edit" % Select vertex for editing
                if isempty(app.ui.LabelTabel), return; end
                vx = app.LabeledPolys(app.ui.LabelTabel).polygon.Vertices(:,1);
                vy = app.LabeledPolys(app.ui.LabelTabel).polygon.Vertices(:,2);
                dists = sqrt((vx - xClick).^2 + (vy - yClick).^2);
                [minDist, minIdx] = min(dists);
                if minDist < 5
                    app.ui.SelectedVertex = minIdx;
                end
    
            elseif mode == "adddelete"
                if isempty(app.ui.LabelTabel), return; end
                vx = app.LabeledPolys(app.ui.LabelTabel).polygon.Vertices(:,1);
                vy = app.LabeledPolys(app.ui.LabelTabel).polygon.Vertices(:,2);
    
                if strcmp(event.Source.SelectionType, 'alt') % Right-click to delete
                    dists = sqrt((vx - xClick).^2 + (vy - yClick).^2);
                    [minDist, minIdx] = min(dists);
    
                    if minDist < 0.2 && numel(vx) > 3 % Prevent deleting last 3 vertices
                        vx(minIdx) = [];
                        vy(minIdx) = [];
                        app.LabeledPolys(app.ui.LabelTabel).polygon = polyshape(vx, vy);
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
                    app.LabeledPolys(app.ui.LabelTabel).polygon = polyshape(vx, vy);
                    updateDisplay();
                end
            end
            % handle mouse click logic
        end

        function moveVertex(app)
            % handle vertex moving logic
        end

        function releaseVertex(app)
            app.ui.SelectedVertex = [];
        end

        function keyPress(app, event)
            if isempty(app.ui.SelectedIdx) || mode ~= "move"
                return;
            end
    
            vx = app.LabeledPolys(app.ui.SelectedIdx).polygon.Vertices(:,1);
            vy = app.LabeledPolys(app.ui.SelectedIdx).polygon.Vertices(:,2);
            step = 1;
    
            switch event.Key
                case 'uparrow'
                    vy = vy - step;
                case 'downarrow'
                    vy = vy + step;
                case 'rightarrow'
                    vx = vx + step;
                case 'leftarrow'
                    vx = vx - step;
            end
    
            app.LabeledPolys(app.ui.SelectedIdx).polygon = polyshape(vx, vy);
            updateDisplay();
        end

        function updateDisplay(app)
            cla(app.Ax);
            imshow(app.ImgStack(:,:,app.CurrentSlice),[], 'Parent', app.Ax); hold(app.Ax, 'on');
            app.plotPolygons();
            % Plot polygons, vertices, and selection status
        end

        function plotPolygons(app)
            hold on;
            for i = 1:numel(app.LabeledPolys)
                plot(app.LabeledPolys(i).polygon, 'FaceColor', 'c', 'FaceAlpha', 0.1);
            end
    
            if ~isempty(app.ui.SelectedIdx)
                plot(app.LabeledPolys(app.ui.SelectedIdx).polygon, 'FaceColor', 'm', 'FaceAlpha', 0.2);
                
                vx = app.LabeledPolys(app.ui.SelectedIdx).polygon.Vertices(:,1);
                vy = app.LabeledPolys(app.ui.SelectedIdx).polygon.Vertices(:,2);
                originalVx = app.LabeledPolys(app.ui.SelectedIdx).originalVertices(:,1);
                originalVy = app.LabeledPolys(app.ui.SelectedIdx).originalVertices(:,2);
    
                isNewVertex = ~ismember([vx vy], [originalVx originalVy], 'rows');
    
                scatter(vx(~isNewVertex), vy(~isNewVertex), 70, 'k', 'filled', 'MarkerEdgeColor', 'y'); % Black for original
                scatter(vx(isNewVertex), vy(isNewVertex), 70, 'g', 'filled', 'MarkerEdgeColor', 'y'); % Green for new
                if ~isempty(app.ui.SelectedVertex)
                    scatter(vx(app.ui.SelectedVertex), vy(app.ui.SelectedVertex), 100, 'b', 'filled', 'MarkerEdgeColor', 'y'); % Selected in blue
                end
            end
        end


        % additional methods for managing polygons and GUI state...
    end
end
