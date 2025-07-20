classdef roiSelector < handle
    properties
        Stack
        ROIType
        PreexistingVertices = []
        IsRGB % bool
        Fig % uifigure object
        HStack % sliceviewer object
        HAxes % sliceviewer axes object
        TheROI % ROI object

        ResetFlag = false 
        MinSlider % contrast slider
        MaxSlider % contrast slider

        OverlayLine % line object for thickness marking 
        MoveListeners = event.listener.empty % line object listener
        WEdit % line width box
        WSlider % line width slider

    end

    methods
        function obj = roiSelector(stack, roi_type, preexistingvertices)
            if nargin > 2
                obj.PreexistingVertices = preexistingvertices;
            end
            obj.Stack = stack;
            obj.ROIType = roi_type;
            obj.IsRGB = ndims(stack) == 4 && size(stack, 3) == 3;
            obj.setupUI();
        end

        function [vertices, ref_slice, mask] = getROI(obj)
            while true
                uiwait(obj.Fig);
                if ~isvalid(obj.Fig)
                    vertices = [];
                    ref_slice = [];
                    return;
                end
                if obj.ResetFlag
                    obj.ResetFlag = false;
                    continue;
                end
                break;
            end

            vertices = obj.getVertices();
            mask = obj.getMask(vertices);
            ref_slice = obj.HStack.SliceNumber;
            close(obj.Fig); % finish and close figure
        end



    end



    methods (Access = private)
        function setupUI(obj)
            obj.normalizeStack();
            obj.Fig = uifigure('Name','Stack Explorer','Position',[100 100 600 950]);
            outer = uigridlayout(obj.Fig);
            outer.RowHeight = {400 ,'1x','1x'};
            outer.ColumnWidth = {'1x'};
            imgPanel = uipanel(outer,'Title','Slice Viewer');
            controlPanel = uipanel(outer,'Title','Console');

            g = uigridlayout(controlPanel,[3 1],'RowHeight',{'fit','fit'},'RowSpacing',8,'Padding',[10 10 10 10]);

            % Intensity Panel
            intensityPanel = uipanel(g,'Title','Intensity'); 
            intensityPanel.Layout.Row = 1;
            gi = uigridlayout(intensityPanel,[2 2],'ColumnWidth',{100,'1x'},'RowSpacing',5);
            uilabel(gi,'Text','Min Intensity:','HorizontalAlignment','left');
            obj.MinSlider = uislider(gi,'Limits',[0 65535],'Value',0, ...
                'ValueChangedFcn', @(s,~) obj.updateIntensity());
            uilabel(gi,'Text','Max Intensity:','HorizontalAlignment','left');
            obj.MaxSlider = uislider(gi,'Limits',[0 65535],'Value',65535, ...
                'ValueChangedFcn', @(s,~) obj.updateIntensity());

            % Button Panel
            buttonPanel = uipanel(g,'Title','Actions'); buttonPanel.Layout.Row = 3;
            gb = uigridlayout(buttonPanel,[1 4], ...
                'ColumnWidth',{'1x','1x', 120, '1x'}, 'ColumnSpacing',10);
            uibutton(gb,'Text','Reset','ButtonPushedFcn', @(~,~) obj.resetROI());
            uibutton(gb,'Text','Confirm','ButtonPushedFcn', @(~,~) uiresume(obj.Fig));

        % Always create WSlider (invisible for non-line)
        obj.WSlider = uislider(gb, 'Limits', [0 30], 'Value', 0, ...
            'ValueChangedFcn', @(src,~) obj.updateLineWidth(src), ...
            'Visible', strcmp(obj.ROIType, 'line'));
        
        % Optional: also create WEdit (only for 'line')
        if strcmp(obj.ROIType,'line')
            obj.WEdit = uieditfield(gb,'numeric','Limits',[0 30],'Value',0,...
                'ValueChangedFcn', @(src,~) obj.syncEditToSlider(src));
            obj.WSlider.Value = 5;  % override to 5 only for line
        end
            % Slice viewer
            if obj.IsRGB
                sliceData = squeeze(obj.Stack(:,:,:,1));
                obj.HStack = sliceViewer(sliceData,'Parent',imgPanel);
            else
                obj.HStack = sliceViewer(obj.Stack,'Parent',imgPanel);
            end

            obj.HAxes = getAxesHandle(obj.HStack);
            obj.TheROI = obj.drawROI();
        end

        function normalizeStack(obj)
            if ~obj.IsRGB
                min_val = min(obj.Stack, [], 'all');
                max_val = max(obj.Stack, [], 'all');
                if isa(obj.Stack, 'double')
                    obj.Stack = uint16((obj.Stack - min_val) / (max_val - min_val) * 65535);
                end
            else
                for i = 1:3
                    ch = obj.Stack(:,:,i,:);
                    min_val = min(ch,[],'all');
                    max_val = max(ch,[],'all');
                    if isa(ch,'double')
                        obj.Stack(:,:,i,:) = uint16((ch - min_val) / (max_val - min_val) * 65535);
                    end
                end
            end
        end

        function updateIntensity(obj)
            obj.HAxes.CLim = [obj.MinSlider.Value, obj.MaxSlider.Value];
        end

        function roi = drawROI(obj)
            if ~isempty(obj.PreexistingVertices)
                switch obj.ROIType
                    case 'rectangle'
                        roi = drawrectangle(obj.HAxes,'Position',obj.PreexistingVertices);
                    case 'polygon'
                        roi = drawpolygon(obj.HAxes,'Position',obj.PreexistingVertices);
                    case 'line'
                        disp('position loaded')
                     roi = drawline(obj.HAxes,'Color','m', ...
                            'LineWidth',obj.WSlider.Value,'LabelVisible','off','Position',obj.PreexistingVertices(1:2,1:2));
                        obj.OverlayLine = line(obj.PreexistingVertices(1:2,1), obj.PreexistingVertices(1:2,2), ...
                            'Parent', obj.HAxes, 'Color', [1 0 1 0.3], ...
                            'LineWidth', roi.LineWidth, 'HitTest','off');
                        obj.MoveListeners(1) = addlistener(roi,'MovingROI', ...
                            @(~,evt)obj.syncOverlay(evt.CurrentPosition));
                        obj.MoveListeners(2) = addlistener(roi,'ROIMoved', ...
                            @(~,evt)obj.syncOverlay(evt.CurrentPosition));
                    otherwise
                        error('Unsupported ROI type.');
                end
            else
                switch obj.ROIType
                    case 'rectangle'
                        roi = drawrectangle(obj.HAxes);
                    case 'polygon'
                        roi = drawpolygon(obj.HAxes);
                    case 'line'
                        roi = drawline(obj.HAxes,'Color','m', ...
                            'LineWidth',obj.WSlider.Value,'LabelVisible','off');
                        obj.OverlayLine = line(roi.Position(:,1), roi.Position(:,2), ...
                            'Parent', obj.HAxes, 'Color', [1 0 1 0.3], ...
                            'LineWidth', roi.LineWidth, 'HitTest','off');
                        obj.MoveListeners(1) = addlistener(roi,'MovingROI', ...
                            @(~,evt)obj.syncOverlay(evt.CurrentPosition));
                        obj.MoveListeners(2) = addlistener(roi,'ROIMoved', ...
                            @(~,evt)obj.syncOverlay(evt.CurrentPosition));
                    otherwise
                        error('Unsupported ROI type.');
                end
            end
        end

        function resetROI(obj)
            if strcmp(obj.ROIType,'line')
                if isvalid(obj.MoveListeners)
                    delete(obj.MoveListeners);
                end
                if isvalid(obj.OverlayLine)
                    delete(obj.OverlayLine);
                end
            end
            if isvalid(obj.TheROI)
                delete(obj.TheROI);
            end
            obj.MoveListeners = event.listener.empty;
            obj.OverlayLine = gobjects(0);
            obj.TheROI = obj.drawROI();
            obj.ResetFlag = true;
            uiresume(obj.Fig);
        end

        function updateLineWidth(obj, src)
            if isvalid(obj.WEdit)
                obj.WEdit.Value = src.Value;
            end
            if isvalid(obj.TheROI)
                obj.TheROI.LineWidth = max(src.Value, 0.5);
            end
            if isvalid(obj.OverlayLine)
                obj.OverlayLine.LineWidth = max(src.Value, 0.5);
            end
        end

        function syncOverlay(obj, pos)
            if isvalid(obj.OverlayLine)
                obj.OverlayLine.XData = pos(:,1);
                obj.OverlayLine.YData = pos(:,2);
                obj.OverlayLine.LineWidth = obj.TheROI.LineWidth;
            end
        end

        function syncEditToSlider(obj, src)
            if isvalid(obj.WSlider)
                obj.WSlider.Value = src.Value;
                obj.updateLineWidth(obj.WSlider);
            end
        end

        function vertices = getVertices(obj)
            vertices = round(obj.TheROI.Position);
            if strcmp(obj.ROIType, 'rectangle')
                x1 = vertices(1); y1 = vertices(2);
                x2 = x1 + vertices(3); y2 = y1 + vertices(4);
                vertices = [x1, y1; x2, y1; x2, y2; x1, y2];
            
            elseif strcmp(obj.ROIType, 'line') 
                vertices = [vertices;[obj.WSlider.Value,-999]];
            end
        end

        function mask = getMask(obj,vertices)
            if strcmp(obj.ROIType, 'line') 
                lwidth = vertices(3,1);
                lvec = vertices(2,:) - vertices(1,:);
                ltheta = atan2(lvec(2),lvec(1));
                llength =norm(lvec);
                pvec = [-lvec(2) lvec(1)]/llength; % perpendicular unit vector
                pvec = (lwidth/2)*pvec; % perpendicular half width vector
                rectVertices = [vertices(1,:)+pvec; vertices(2,:)+pvec; vertices(2,:)-pvec; vertices(1,:) - pvec];
                mask = poly2mask(rectVertices(:,1), rectVertices(:,2), size(obj.Stack(:,:,1),1), size(obj.Stack(:,:,1),2));
            else
                disp('getmask, mask area = ')
                mask = obj.TheROI.createMask;
                disp(sum(mask,'all'))
                uiresume(obj.Fig)
            end
        end
    end
end
