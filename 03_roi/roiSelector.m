classdef roiSelector < handle
    properties
        Stack
        ROIType
        PreexistingVertices = []
        IsRGB % bool
        Fig % uifigure object
        HAxes % axes handle (uiaxes)
        TheROI % ROI object
        ResetFlag = false 
        MinSlider % contrast(min) slider
        MaxSlider % contrast(max) slider
        OverlayLine % line object for thickness marking 
        MoveListeners = event.listener.empty % line object listener
        WEdit % line width box
        WSlider % line width slider
        Displaystack % imshow handle
    end
    properties
        frame = 1;
        minIntensity = 0;
        maxIntensity = 65535;
        fontsize = 14
        nframe = 0;
    end
    properties (Access=private)
        FrameSlider
        FrameLabel
    end

    methods
        function obj = roiSelector(stack, roi_type, preexistingvertices)
            if nargin > 2
                obj.PreexistingVertices = preexistingvertices;
            end
            obj.Stack = stack;
            obj.nframe = size(stack,3);
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
            ref_slice = obj.frame;
            close(obj.Fig); % finish and close figure
        end
    end

    methods (Access = private)
        function setupUI(obj)
            obj.normalizeStack();

            % 0. Main uiFigure
            obj.Fig = uifigure('Name','Stack Explorer','Position',[100 100 600 950]);

            % 0. Main grid
            mainlayout = uigridlayout(obj.Fig, [3,1]); % [이미지, 파라미터, 콘솔, 버튼]
            mainlayout.RowHeight   = {'1x','fit','fit'};
            mainlayout.ColumnWidth = {'1x'};


            % 1) Image panel
            imgPanel = uipanel(mainlayout,'Title','Slice Viewer');
            imgPanel.Layout.Row = 1; imgPanel.Layout.Column = 1;

            imgPanel_layout = uigridlayout(imgPanel, [1,1]);
            obj.HAxes = uiaxes(imgPanel_layout);
            obj.HAxes.Layout.Row = 1; 
            obj.HAxes.Layout.Column = 1;
            obj.HAxes.Toolbar.Visible = 'off';
            disableDefaultInteractivity(obj.HAxes);
            axis(obj.HAxes,'image'); axis(obj.HAxes,'off');

            obj.Displaystack = imshow(obj.Stack(:,:,obj.frame), [], ... % DISPLAY IMAGE
                'Parent', obj.HAxes, 'InitialMagnification','fit'); 
            obj.HAxes.CLim = [obj.minIntensity obj.maxIntensity];

           
            % 3) Console panel
            controlPanel = uipanel(mainlayout,'Title','Console','Scrollable','on');
            controlPanel.Layout.Row = 2;
            controlPanel.Layout.Column = 1;

            controlPanelLayout = uigridlayout(controlPanel,[3 1],'RowSpacing',8,'Padding',[10 10 10 10]);
            controlPanelLayout.RowHeight   = {'fit','fit','fit'};
            controlPanelLayout.ColumnWidth = {'1x'};

            % 3-1) Intensity (Min/Max) – rangeslider 없이 두 개의 uislider로 window 제어
            intensityPanel = uipanel(controlPanelLayout,'Title','Intensity Window (Min/Max)');
            intensityPanel.Layout.Row = 1; intensityPanel.Layout.Column = 1;

            gi = uigridlayout(intensityPanel,[3 2],'ColumnWidth',{110,'1x'},'RowSpacing',6,'Padding',[8 8 8 8]);

            uilabel(gi,'Text','Min Intensity:','HorizontalAlignment','left');
            obj.MinSlider = uislider(gi,'Limits',[0 65535],'Value',obj.minIntensity, ...
                'ValueChangedFcn', @(s,~) obj.updateIntensity());

            uilabel(gi,'Text','Max Intensity:','HorizontalAlignment','left');
            obj.MaxSlider = uislider(gi,'Limits',[0 65535],'Value',obj.maxIntensity, ...
                'ValueChangedFcn', @(s,~) obj.updateIntensity());

            % 3-2) Slice control
            slicesliderPanel = uipanel(controlPanelLayout, "Title", 'Slice control');
            slicesliderPanel.Layout.Row = 2; slicesliderPanel.Layout.Column = 1;

            sliceslider_layout = uigridlayout(slicesliderPanel, [1,2], 'ColumnWidth', {'1x','fit'}, 'RowHeight', {'fit'});
            obj.FrameSlider = uislider(sliceslider_layout, ...
                'Limits', [1, max(1,obj.nframe)], ...
                'Value', obj.frame, ...
                'MajorTicks', round(linspace(1,max(1,obj.nframe),min(6,obj.nframe))), ...
                'ValueChangingFcn', @(src, event) obj.Update_Frame(round(event.Value)), ...
                'ValueChangedFcn',  @(src, event) obj.Update_Frame(round(src.Value)));
            obj.FrameLabel = uilabel(sliceslider_layout, ...
                'Text', sprintf('Frame %d / %d', obj.frame, obj.nframe), ...
                'HorizontalAlignment','left','FontSize',obj.fontsize);

            % 3-3) Line width (for ROI 'line')
            linewidthPanel = uipanel(controlPanelLayout, "Title", 'Line width (ROI: line)');
            linewidthPanel.Layout.Row = 3; linewidthPanel.Layout.Column = 1;

            glw = uigridlayout(linewidthPanel,[1,3],'ColumnWidth',{'fit','1x','fit'},'RowSpacing',6,'Padding',[8 8 8 8]);
            uilabel(glw,'Text','Width:');
            obj.WSlider = uislider(glw, 'Limits', [0 30], 'Value', 5, ...
                'ValueChangedFcn', @(src,~) obj.updateLineWidth(src));
            obj.WEdit = uieditfield(glw,'numeric','Limits',[0 30],'Value',5, ...
                'ValueChangedFcn', @(src,~) obj.syncEditToSlider(src));

            % 4) 버튼 패널
            buttonPanel = uipanel(mainlayout); 
            buttonPanel.Layout.Row = 3; 
            buttonPanel.Layout.Column = 1;
            button_layout = uigridlayout(buttonPanel, [1 2],'ColumnWidth',{'1x','1x'},'Padding',[6 6 6 6]);
            uibutton(button_layout,'Text','Reset','ButtonPushedFcn', @(~,~) obj.resetROI());
            uibutton(button_layout,'Text','Confirm','ButtonPushedFcn', @(~,~) uiresume(obj.Fig));

            % ROI 생성 (sliceViewer 제거 버전)
            obj.TheROI = obj.drawROI();
        end

        function normalizeStack(obj)
            if ~obj.IsRGB
                min_val = min(obj.Stack, [], 'all');
                max_val = max(obj.Stack, [], 'all');
                if isa(obj.Stack, 'double')
                    obj.Stack = uint16((obj.Stack - min_val) / max(1,(max_val - min_val)) * 65535);
                end
            else
                for i = 1:3
                    ch = obj.Stack(:,:,i,:);
                    min_val = min(ch,[],'all');
                    max_val = max(ch,[],'all');
                    if isa(ch,'double')
                        obj.Stack(:,:,i,:) = uint16((ch - min_val) / max(1,(max_val - min_val)) * 65535);
                    end
                end
            end
        end

        function updateIntensity(obj)
            % 두 슬라이더 값 정렬(교차 방지)
            lo = min(obj.MinSlider.Value, obj.MaxSlider.Value);
            hi = max(obj.MinSlider.Value, obj.MaxSlider.Value);
            obj.minIntensity = lo; obj.maxIntensity = hi;
            if isvalid(obj.HAxes)
                obj.HAxes.CLim = [lo hi];
            end
        end

        function roi = drawROI(obj)
            % drawROI는 imshow의 axes(obj.HAxes)에서 동작
            if ~isempty(obj.PreexistingVertices)
                switch obj.ROIType
                    case 'rectangle'
                        roi = drawrectangle(obj.HAxes,'Position',obj.PreexistingVertices);
                    case 'polygon'
                        roi = drawpolygon(obj.HAxes,'Position',obj.PreexistingVertices);
                    case 'line'
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
                if ~isempty(obj.MoveListeners)
                    try
                        delete(obj.MoveListeners); 
                    catch
                    end
                end
                if ~isempty(obj.OverlayLine) && isvalid(obj.OverlayLine)
                    delete(obj.OverlayLine);
                end
            end
            if ~isempty(obj.TheROI) && isvalid(obj.TheROI)
                delete(obj.TheROI);
            end
            obj.MoveListeners = event.listener.empty;
            obj.OverlayLine = gobjects(0);
            obj.TheROI = obj.drawROI();
            obj.ResetFlag = true;
            uiresume(obj.Fig);
        end

        function updateLineWidth(obj, src)
            if ~isempty(obj.WEdit) && isvalid(obj.WEdit)
                obj.WEdit.Value = src.Value;
            end
            if ~isempty(obj.TheROI) && isvalid(obj.TheROI)
                obj.TheROI.LineWidth = max(src.Value, 0.5);
            end
            if ~isempty(obj.OverlayLine) && isvalid(obj.OverlayLine)
                obj.OverlayLine.LineWidth = max(src.Value, 0.5);
            end
        end

        function syncOverlay(obj, pos)
            if ~isempty(obj.OverlayLine) && isvalid(obj.OverlayLine)
                obj.OverlayLine.XData = pos(:,1);
                obj.OverlayLine.YData = pos(:,2);
                if ~isempty(obj.TheROI) && isvalid(obj.TheROI)
                    obj.OverlayLine.LineWidth = obj.TheROI.LineWidth;
                end
            end
        end

        function syncEditToSlider(obj, src)
            if ~isempty(obj.WSlider) && isvalid(obj.WSlider)
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
                vertices = [vertices; [obj.WSlider.Value, -999]];
            end
        end

        function mask = getMask(obj,vertices)
            if strcmp(obj.ROIType, 'line')
                lwidth = vertices(3,1);
                lvec = vertices(2,:) - vertices(1,:);
                llength = norm(lvec) + eps;
                pvec = [-lvec(2) lvec(1)]/llength; % perpendicular unit vector
                pvec = (lwidth/2)*pvec; % perpendicular half width vector
                rectVertices = [vertices(1,:)+pvec; vertices(2,:)+pvec; vertices(2,:)-pvec; vertices(1,:)-pvec];
                mask = poly2mask(rectVertices(:,1), rectVertices(:,2), size(obj.Stack,1), size(obj.Stack,2));
            else
                mask = obj.TheROI.createMask;
                uiresume(obj.Fig)
            end
        end

        % ===== imshow 기반 콜백들 =====

        function Update_Frame(obj, NewFrame)
            NewFrame = max(1, min(obj.nframe, NewFrame));
            if NewFrame == obj.frame, return; end
            obj.frame = NewFrame;

            if isvalid(obj.Displaystack)
                obj.Displaystack.CData = obj.Stack(:,:,obj.frame);
            end

            % 라벨/슬라이더 동기화
            if ~isempty(obj.FrameLabel) && isvalid(obj.FrameLabel)
                obj.FrameLabel.Text = sprintf('Frame %d / %d', obj.frame, obj.nframe);
            end
            if ~isempty(obj.FrameSlider) && isvalid(obj.FrameSlider)
                if obj.FrameSlider.Value ~= obj.frame
                    obj.FrameSlider.Value = obj.frame;
                end
            end
        end
    end
end
