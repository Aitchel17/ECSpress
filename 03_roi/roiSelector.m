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
        Stack2       % optional 2nd stack (e.g. ch2) for the reflect view
        ReflectFig   % reflect figure handle (classic figure)
        ReflectAx    % reflect axes handle
        ReflectImg   % reflect image handle (CData reuse)
        ReflectOverlay = gobjects(0)           % ROI overlay on reflect axes
        ReflectListener = event.listener.empty % ROIMoved -> updateReflect
        ProfileFig   % line-profile figure handle
        ProfileAx    % line-profile axes
        ProfileLine1 % ch1 profile line (red)
        ProfileLine2 % ch2 profile line (green)
        ProfileDrop  % projection-mode dropdown (mean/median/max/min)
        ProfileListeners = event.listener.empty % MovingROI/ROIMoved -> updateProfile
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
        function obj = roiSelector(stack, roi_type, preexistingvertices, stack2)
            if nargin > 2
                obj.PreexistingVertices = preexistingvertices;
            end
            if nargin > 3
                obj.Stack2 = stack2;   % optional co-registered 2nd channel (ch2)
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
            obj.Fig = uifigure('Name','Stack Explorer','Position',[100 100 1000 700]);
            obj.Fig.DeleteFcn = @(~,~) obj.closeChildWindows();  % close reflect/profile with main

            % 0. Main grid
            % Layout: 2 Rows, 2 Cols
            % Col 1: Image (Spans Row 1-2)
            % Col 2: Row 1=Console, Row 2=Button
            mainlayout = uigridlayout(obj.Fig, [2,2]);
            mainlayout.RowHeight   = {'1x','fit'};
            mainlayout.ColumnWidth = {'1x', 350};


            % 1) Image panel (Left, Full Height)
            imgPanel = uipanel(mainlayout,'Title','Slice Viewer');
            imgPanel.Layout.Row = [1 2];
            imgPanel.Layout.Column = 1;

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


            % 3) Console panel (Right, Top)
            controlPanel = uipanel(mainlayout,'Title','Console','Scrollable','on');
            controlPanel.Layout.Row = 1;
            controlPanel.Layout.Column = 2;

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
                'ValueChangedFcn',  @(src, event) obj.onFrameSettled(round(src.Value)));
            obj.FrameLabel = uilabel(sliceslider_layout, ...
                'Text', sprintf('Frame %d / %d', obj.frame, obj.nframe), ...
                'HorizontalAlignment','left','FontSize',obj.fontsize);

            % 3-3) Line width (for ROI 'line')
            linewidthPanel = uipanel(controlPanelLayout, "Title", 'Line width (ROI: line)');
            linewidthPanel.Layout.Row = 3; linewidthPanel.Layout.Column = 1;

            glw = uigridlayout(linewidthPanel,[2,3],'ColumnWidth',{'fit','1x','fit'},'RowSpacing',6,'Padding',[8 8 8 8]);
            uilabel(glw,'Text','Width:');
            obj.WSlider = uislider(glw, 'Limits', [0 30], 'Value', 5, ...
                'ValueChangedFcn', @(src,~) obj.updateLineWidth(src));
            obj.WEdit = uieditfield(glw,'numeric','Limits',[0 30],'Value',5, ...
                'ValueChangedFcn', @(src,~) obj.syncEditToSlider(src));
            uilabel(glw,'Text','Projection:');
            obj.ProfileDrop = uidropdown(glw,'Items',{'mean','median','max','min'},'Value','mean', ...
                'ValueChangedFcn', @(~,~) obj.updateProfile());

            % 4) button pannel (Right, Bottom)
            buttonPanel = uipanel(mainlayout);
            buttonPanel.Layout.Row = 2;
            buttonPanel.Layout.Column = 2;
            button_layout = uigridlayout(buttonPanel, [1 3],'ColumnWidth',{'1x','1x','1x'},'Padding',[6 6 6 6]);
            uibutton(button_layout,'Text','Reset','ButtonPushedFcn', @(~,~) obj.resetROI());
            uibutton(button_layout,'Text','Profile','ButtonPushedFcn', @(~,~) obj.showLineProfile());
            uibutton(button_layout,'Text','Confirm','ButtonPushedFcn', @(~,~) uiresume(obj.Fig));

            % ROI 생성 (sliceViewer 제거 버전)
            obj.TheROI = obj.drawROI();

            % ch2 reflect view (optional 2nd stack): create window + first draw
            if ~isempty(obj.Stack2)
                obj.attachReflect();
            end
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
                        if size(obj.PreexistingVertices,1) == 4 && size(obj.PreexistingVertices,2)==2
                            x = min(obj.PreexistingVertices(:,1));
                            y = min(obj.PreexistingVertices(:,2));
                            w = max(obj.PreexistingVertices(:,1)) - x;
                            h = max(obj.PreexistingVertices(:,2)) - y;
                            roi = drawrectangle(obj.HAxes,'Position',[x y w h]);
                        else
                            roi = drawrectangle(obj.HAxes,'Position',obj.PreexistingVertices);
                        end
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
            if ~isempty(obj.Stack2), obj.attachReflect(); end
            if ~isempty(obj.ProfileFig) && isvalid(obj.ProfileFig)
                obj.attachProfileListener();
                obj.updateProfile();
            end
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
            obj.updateProfile();   % new width -> refresh profile band if open
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
         
        end

        function showLineProfile(obj)
            % Pop up the intensity profile along the current line ROI (current frame).
            if ~strcmp(obj.ROIType, 'line')
                warning('Line profile is only available for line ROIs.');
                return;
            end
            if isempty(obj.TheROI) || ~isvalid(obj.TheROI)
                warning('No line ROI drawn.');
                return;
            end
            obj.ensureProfileFig();
            obj.attachProfileListener();   % live-track while dragging the ROI
            obj.updateProfile();
        end

        function ensureProfileFig(obj)
            % Lazily create the profile window with reusable ch1(red)/ch2(green) lines.
            if ~isempty(obj.ProfileFig) && isvalid(obj.ProfileFig), return; end
            obj.ProfileFig = figure('Name','Line profile','NumberTitle','off');
            obj.ProfileAx  = axes(obj.ProfileFig);
            hold(obj.ProfileAx,'on'); grid(obj.ProfileAx,'on');
            obj.ProfileLine1 = plot(obj.ProfileAx, NaN, NaN, 'r-', 'LineWidth', 1.2); % ch1
            obj.ProfileLine2 = plot(obj.ProfileAx, NaN, NaN, 'g-', 'LineWidth', 1.2); % ch2
            xlabel(obj.ProfileAx, 'Position along line (px)');
            ylabel(obj.ProfileAx, 'Normalized (5-95% capped)');
            ylim(obj.ProfileAx, [-0.02 1.02]);
            legend(obj.ProfileAx, {'ch1','ch2'}, 'Location','best');
        end

        function attachProfileListener(obj)
            % Live-track the profile during ROI drag (MovingROI) and on release (ROIMoved).
            if ~isempty(obj.ProfileListeners)
                delete(obj.ProfileListeners);
            end
            obj.ProfileListeners = event.listener.empty;
            if ~isempty(obj.TheROI) && isvalid(obj.TheROI)
                obj.ProfileListeners(1) = addlistener(obj.TheROI,'MovingROI',@(~,~) obj.updateProfile());
                obj.ProfileListeners(2) = addlistener(obj.TheROI,'ROIMoved', @(~,~) obj.updateProfile());
            end
        end

        function updateProfile(obj)
            % Recompute + redraw ch1/ch2 profiles (cheap: update line data, no clf).
            if isempty(obj.ProfileFig) || ~isvalid(obj.ProfileFig), return; end
            if isempty(obj.TheROI) || ~isvalid(obj.TheROI), return; end
            if ~strcmp(obj.ROIType,'line'), return; end
            if ~isempty(obj.ProfileDrop) && isvalid(obj.ProfileDrop)
                mode = obj.ProfileDrop.Value;
            else
                mode = 'mean';
            end
            [t1, prof1] = obj.bandProfile(obj.Stack(:,:,obj.frame), mode);
            set(obj.ProfileLine1, 'XData', t1, 'YData', prof1);
            if ~isempty(obj.Stack2)
                f2 = min(obj.frame, size(obj.Stack2,3));
                [t2, prof2] = obj.bandProfile(obj.Stack2(:,:,f2), mode);
                set(obj.ProfileLine2, 'XData', t2, 'YData', prof2, 'Visible','on');
            elseif ~isempty(obj.ProfileLine2) && isvalid(obj.ProfileLine2)
                set(obj.ProfileLine2, 'Visible','off');
            end
            title(obj.ProfileAx, sprintf('%s profile @ frame %d', mode, obj.frame));
        end

        function [tline, prof] = bandProfile(obj, img, mode)
            % Width-band profile of img along the current line ROI, projected by mode.
            % Band extraction delegated to analyze_interp2 (drop-in for analyze_affine_rotate).
            v    = obj.TheROI.Position;                 % [x1 y1; x2 y2]
            w    = max(round(obj.WSlider.Value), 1);
            band = analyze_interp2(img, v, w);          % [w x L] straightened band (this frame)
            switch mode
                case 'median', prof = median(band, 1, 'omitnan');
                case 'max',    prof = max(band, [], 1, 'omitnan');
                case 'min',    prof = min(band, [], 1, 'omitnan');
                otherwise,     prof = mean(band, 1, 'omitnan');
            end
            prof = prof(:).';                           % row vector along line
            % cap to [5th, 95th] percentile, then normalize to [0,1]
            lo = prctile(prof, 5);
            hi = prctile(prof, 95);
            if hi > lo
                prof = (min(max(prof, lo), hi) - lo) / (hi - lo);
            end
            tline = 1:numel(prof);                      % position along line (px)
        end

        % ===== ch2 reflect view =====

        function onFrameSettled(obj, NewFrame)
            % Slider RELEASE only (ValueChangedFcn): update ch1 frame + reflect + profile.
            obj.Update_Frame(NewFrame);
            obj.updateReflect();
            obj.updateProfile();
        end

        function attachReflect(obj)
            % Create/refresh reflect window and (re)bind ROIMoved -> updateReflect.
            obj.ensureReflect();
            if ~isempty(obj.ReflectListener)
                delete(obj.ReflectListener);   % no-op if already invalid
            end
            obj.ReflectListener = event.listener.empty;
            if ~isempty(obj.TheROI) && isvalid(obj.TheROI)
                obj.ReflectListener = addlistener(obj.TheROI,'ROIMoved',@(~,~) obj.updateReflect());
            end
            obj.updateReflect();
        end

        function ensureReflect(obj)
            % Lazily create the classic reflect figure (display-only: no drag bug, faster).
            if isempty(obj.Stack2), return; end
            if ~isempty(obj.ReflectFig) && isvalid(obj.ReflectFig), return; end
            f = min(obj.frame, size(obj.Stack2,3));
            obj.ReflectFig = figure('Name','Reflect on ch2','NumberTitle','off');
            obj.ReflectAx  = axes(obj.ReflectFig);
            obj.ReflectImg = imshow(obj.Stack2(:,:,f), [], 'Parent', obj.ReflectAx);
            lo = double(min(obj.Stack2(:))); hi = double(max(obj.Stack2(:)));
            if hi > lo, obj.ReflectAx.CLim = [lo hi]; end
            axis(obj.ReflectAx,'image');
            hold(obj.ReflectAx,'on');
            obj.ReflectOverlay = gobjects(0);
        end

        function updateReflect(obj)
            % Cheap update: swap ch2 CData (current frame) + redraw current ROI overlay.
            if isempty(obj.Stack2), return; end
            obj.ensureReflect();
            if isempty(obj.ReflectImg) || ~isvalid(obj.ReflectImg), return; end
            f = min(obj.frame, size(obj.Stack2,3));
            obj.ReflectImg.CData = obj.Stack2(:,:,f);
            old = obj.ReflectOverlay;
            if ~isempty(old)
                delete(old(isvalid(old)));
            end
            obj.ReflectOverlay = gobjects(0);
            if ~isempty(obj.TheROI) && isvalid(obj.TheROI)
                v = obj.TheROI.Position;
                switch obj.ROIType
                    case 'line'
                        obj.ReflectOverlay = plot(obj.ReflectAx, v(:,1), v(:,2), 'm-', ...
                            'LineWidth', max(obj.WSlider.Value,0.5));
                    case 'rectangle'
                        obj.ReflectOverlay = rectangle(obj.ReflectAx, 'Position', v, ...
                            'EdgeColor','m', 'LineWidth',1.5);
                    case 'polygon'
                        obj.ReflectOverlay = plot(obj.ReflectAx, v([1:end 1],1), v([1:end 1],2), ...
                            'm-', 'LineWidth',1.5);
                end
            end
            title(obj.ReflectAx, sprintf('ch2 @ frame %d', f));
        end

        function closeChildWindows(obj)
            % Close the reflect + profile windows when the main figure is destroyed.
            if ~isempty(obj.ReflectFig) && isvalid(obj.ReflectFig)
                delete(obj.ReflectFig);
            end
            if ~isempty(obj.ProfileFig) && isvalid(obj.ProfileFig)
                delete(obj.ProfileFig);
            end
        end
    end
end
