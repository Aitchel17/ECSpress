savepath = 'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\01_primary_analysis\01_pvscollapse';
vivo_50um1 = mdfExtractLoader();
vivo_50um1.analog = vivo_50um1.loadanalog;
vivo_50um1.stackch1 = vivo_50um1.loadstack("ch1");
vivo_50um1.stackch2 = vivo_50um1.loadstack("ch2");
%%
tmp.gpstack =  analyze_grouproject(vivo_50um1.stackch2,3,"mean");
tmp.medcsf = medfilt3(tmp.gpstack,[1 1 5]);
gausscsf= imgaussfilt(tmp.medcsf,1);
tmp.gpstack = analyze_grouproject(vivo_50um1.stackch1,3,"mean");
tmp.medbv = medfilt3(tmp.gpstack,[1 1 5]);
gaussbv= imgaussfilt(tmp.medbv,1);
%%
clc
roi1 = roi(gaussbv,'pax','rectangle',vivo_50um1);
%%
clc
roi1.modifyroi(gaussbv,'pax')


%%
% class initialization part
obj.Stack = gaussbv;
obj.ROIType = 'rectangle';
obj.IsRGB = ndims(obj.Stack) == 4 && size(obj.Stack, 3) == 3;
% objstupui function called and below is that
obj = normalizeStack(obj); % normalize stack
%% generate window
load mristack
Fig = uifigure('Name','Stack Explorer','Position',[100 100 900 800]);
Fig.Resize = "off";
outer = uigridlayout(Fig);%, 'RowHeight',{'1x','fit'}, 'RowSpacing',6, 'Padding',[8 8 8 8]);
outer.RowHeight = {500,150};
outer.ColumnWidth = {'1x'};
imgPanel = uipanel(outer,'Title','Slice Viewer');
controlPanel = uipanel(outer,'Title','Console');
HStack = sliceViewer(mristack,'Parent',imgPanel);
HAxes = getAxesHandle(HStack);
roi = drawrectangle(HAxes);


%%

%%
delete(controlPanel.Children)
%%
g = uigridlayout(controlPanel);
intensityPanel = uipanel(g,'Title','Intensity'); 
%%
gi = uigridlayout(intensityPanel);
%%
uilabel(gi,'Text','Max Intensity:','HorizontalAlignment','left');



%%
outer.RowHeight = {220,220,'1x'};
outer.ColumnWidth = {700,'1x'};

% Intensity Panel
intensityPanel.Layout.Row = 1;

%%
uilabel(gi,'Text','Min Intensity:','HorizontalAlignment','left');
obj.MinSlider = uislider(gi,'Limits',[0 65535],'Value',0, ...
    'ValueChangedFcn', @(s,~) obj.updateIntensity());

uilabel(gi,'Text','Max Intensity:','HorizontalAlignment','left');
obj.MaxSlider = uislider(gi,'Limits',[0 65535],'Value',65535, ...
    'ValueChangedFcn', @(s,~) obj.updateIntensity());
%%
g = uigridlayout(controlPanel,[3 1],'RowHeight',{'fit','fit','fit'},'RowSpacing',8,'Padding',[10 10 10 10]);

% Intensity Panel
intensityPanel = uipanel(g,'Title','Intensity'); 
intensityPanel.Layout.Row = 1;
gi = uigridlayout(intensityPanel,[2 2],'ColumnWidth',{100,'1x'},'RowSpacing',5);
uilabel(gi,'Text','Min Intensity:','HorizontalAlignment','left');
%%
obj.MinSlider = uislider(gi,'Limits',[0 65535],'Value',0, ...
    'ValueChangedFcn', @(s,~) obj.updateIntensity());
uilabel(gi,'Text','Max Intensity:','HorizontalAlignment','left');
obj.MaxSlider = uislider(gi,'Limits',[0 65535],'Value',65535, ...
    'ValueChangedFcn', @(s,~) obj.updateIntensity());

% Button Panel
buttonPanel = uipanel(g,'Title','Actions'); buttonPanel.Layout.Row = 3;
gb = uigridlayout(buttonPanel,[1 4], ...
    'ColumnWidth',{'1x','1x', '1x', '1x'}, 'ColumnSpacing',10);
uibutton(gb,'Text','Reset','ButtonPushedFcn', @(~,~) obj.resetROI());
uibutton(gb,'Text','Confirm','ButtonPushedFcn', @(~,~) uiresume(obj.Fig));



%%
function obj = normalizeStack(obj)
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