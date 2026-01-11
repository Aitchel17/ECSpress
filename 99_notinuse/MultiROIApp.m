classdef MultiROIApp < matlab.apps.AppBase
    properties (Access = public)
        UIFigure       matlab.ui.Figure
        Stacks         cell
        SliceViewers   matlab.ui.control.UIAxes
        ROIContainers  cell
        ConfirmButton  matlab.ui.control.Button
        ResetButton    matlab.ui.control.Button
    end
    
    methods (Access = public)
        function app = MultiROIApp(stack1, stack2)
            if nargin < 2
                stack2 = [];
            end
            
            app.Stacks = {stack1, stack2};
            app.Stacks = app.Stacks(~cellfun('isempty', app.Stacks));
            numStacks = numel(app.Stacks);
            
            % Create figure
            app.UIFigure = uifigure('Name', 'Multi-ROI Viewer', 'Position', [100, 100, 1000, 600]);
            app.SliceViewers = gobjects(1, numStacks);
            app.ROIContainers = cell(1, numStacks);
            
            for i = 1:numStacks
                panel = uipanel(app.UIFigure, 'Title', sprintf('Stack %d', i), 'Position', [(i-1)*500 + 20, 100, 480, 400]);
                app.SliceViewers(i) = sliceViewer(app.Stacks{i}, 'Parent', panel);
                app.ROIContainers{i} = {};  % Store multiple ROIs
            end
            
            % Add control buttons
            app.ConfirmButton = uibutton(app.UIFigure, 'Text', 'Confirm', 'Position', [850, 20, 100, 40], 'ButtonPushedFcn', @(~,~) app.confirmROIs());
            app.ResetButton = uibutton(app.UIFigure, 'Text', 'Reset', 'Position', [730, 20, 100, 40], 'ButtonPushedFcn', @(~,~) app.resetROIs());
            
            % Allow ROI drawing on first stack
            addlistener(app.SliceViewers(1), 'SliceChanged', @(src, ~) app.syncSlices(src));
            app.drawROIs();
        end
        
        function drawROIs(app)
            axesHandle = getAxesHandle(app.SliceViewers(1));
            while true
                roi = drawrectangle(axesHandle);
                if isempty(roi)
                    break;
                end
                app.ROIContainers{1}{end+1} = roi;
                app.syncROIs(roi.Position);
            end
        end
        
        function syncROIs(app, position)
            numStacks = numel(app.Stacks);
            for i = 2:numStacks
                axesHandle = getAxesHandle(app.SliceViewers(i));
                app.ROIContainers{i}{end+1} = drawrectangle(axesHandle, 'Position', position);
            end
        end
        
        function syncSlices(app, src)
            for i = 1:numel(app.SliceViewers)
                if app.SliceViewers(i) ~= src
                    app.SliceViewers(i).SliceNumber = src.SliceNumber;
                end
            end
        end
        
        function confirmROIs(app)
            disp('ROIs Confirmed.');
            for i = 1:numel(app.ROIContainers)
                for j = 1:numel(app.ROIContainers{i})
                    disp(['Stack ', num2str(i), ' ROI ', num2str(j), ': ', mat2str(app.ROIContainers{i}{j}.Position)]);
                end
            end
        end
        
        function resetROIs(app)
            for i = 1:numel(app.ROIContainers)
                for j = 1:numel(app.ROIContainers{i})
                    if isvalid(app.ROIContainers{i}{j})
                        delete(app.ROIContainers{i}{j});
                    end
                end
                app.ROIContainers{i} = {};
            end
            app.drawROIs();
        end
    end
end
