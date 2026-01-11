function hf = imagei(varargin)
% Usage: hf = imagei([x,y,t],im,[seg],[plt_coord],[MovName]);
%
% Given "im", a two-dimensional field that varies over time or depth,
% imagei visualizes the data via an interactive figure in which the user
% can change the slice (time, depth, etc.), play a movie of the data,
% change whether the slice label is displayed as a title or overlaid, and
% save a movie of the data to disk. If "im" is a scalar field like
% temperature, provide it as a three- dimensional rectangular array with
% data arranged for a grid produced by meshgrid.m (not ndgrid.m), with time
% (or depth) varying along the last dimension. If "im" is a vector field
% like velocity, provide it as a cell array of two three-dimensional
% rectangular arrays, one for each component (e.g. {u,v}). If "im" is a
% series of RGB images, provide it as a cell array of three
% three-dimensional rectangular arrays, one for each color channel (e.g.
% {r,g,b}). Optionally specify the spatial and temporal grid using "x",
% "y", and "t", which can either be vectors listing the grid locations and
% times (or depths), or scalars representing the grid sizes and time step
% (or slice thickness). The graphics handle of the resulting figure is
% provided as "hf". To save a movie of the visualization, enter a file name
% ("movie name"; the optional input "MovName" is used by default) and click
% "save...". The duration, downsampling, and frame rate can then be
% adjusted. Conversion from a grid made with ndgrid.m can be done via
% x_mesh = permute(x_nd,[2 1 3]). To overlay a binary segmentation on top
% of the image, optionally provide a cell array seg with the same structure
% as im (either a 1X3 or 1x2 cell containing 3D arrays). To plot coordinates
% on top of the image, provide a list of x,y coordinates plt_coord as an
% N X 2 X S array, where N is the number of points and S is the image number.
% Calls Stack exchange file superSlider, which can be found in
% codelib/boster/FileExchange/. See also slicei.m. 

% Written 14 May 2020 by Doug Kelley, based largely on slicei.m. 
% Updated 19 May 2020 with option to save axes only, otherwise save all but
% controls. 
% Updated 24 May 2020 to make input syntax match image.m. 
% Update 29 Sep 2020 to hide timer and related controls if user provides
% only one frame. 
% Updated Fall 2020 by Kimberly Boster to enable the user to multiply each channel
% by a multiplier, so you can isolate a channel by turning other channels
% off, or oversaturate a channel so that it appears brighter. 
% Updated Fall 2020 to accept an optional input,
% seg, the same structure type as im, but with a binary
% segmentation in each channel. The results of the segmentation are then
% overlaid on im, and the user can adjust the transparency of the overlay
% with text box trans.
% Update 4 Jan 2021 KASB added the capability to plot coordinates, plt_coord
% Update 25 Feb 2021 KASB added sliders for the clim for each channel and a
% separate  segmentation transparency for each channel
% Updated 14 March 2021 by Doug Kelley: Layout all in pixels so controls
% don't resize; renamed variables "slice" instead of "time".
% Updated 16 March 2021 by Doug Kelley: more layout adjustments. 
% Updated 19 March 2021 to handle resizing window with arbitrary units. 
% Updated 5 April 2021 to make 'show curl' text visible.
% Updated 16 April 2021: no longer imposing DataAspectRatio.
% Updated 24 April 2021: Fixed typo (clim_ch23 changed to clim_ch3). Now
% catching non-number inputs to all text boxes expecting numbers. 
% Update 1 June 2021: Default movie playback frame rate is set to match
% mean frame rate in "t" vector, if provided, so playback is in real time. 
% 23 November 2021: Added dialog box to control timing, downsampling, frame
% rate, and file type when saving. 
% 24 Jan 2022 KASB Fixed typo in lines setting inf and nans to zeros
% Updated 27 Mar 2022: HES added section to allow for seg data with <3
% channels. Added toggle to view segmentations as outlines or overlays. 
% Changed trans to a vec to loop through seg channles for plotting
% Updated 02 Mar 2022: KASB added the capability to link the channel limit
% sliders

% -=- Object positions -=-
FigSize = [560 420]; % size of figure
AxesBottomPos = 62; % bottom of axes, if multiple slices, pixels
AxesBottomPosOneSlice = 0; % bottom of axes, if one slice, pixels
AxesLeftPos = 0; % left of axes, if no RGB controls, pixels
AxesLeftPosRGB = 156; % left of axes, if no RGB controls, pixels
SliceLabelPos = [2 5 39 20]; % position of time label, pixels
SliceTextPos = [41 5 40 25]; % position of time text box, pixels
SliceSliderPos = [85 7 120 20]; % position of time slider, pixels
SliceUnitLabelPos = [212 5 35 27]; % position of time unit label, pixels
SliceUnitTextPos = [247 5 45 25]; % position of time unit text box, pixels
MovNameLabelPos = [303 5 42 27]; % position of movie name label, pixels
MovNameTextPos = [345 5 140 25]; % position of movie name text box, pixels
FrameRateLabelPos = [488 5 27 20]; % position of frame rate label, pixels
FrameRateTextPos = [520 5 35 25]; % position of frame rate text box, pixels
QuiverScaleLabelPos = [225 35 45 27]; % position of quiver scale label, pixels
QuiverScaleTextPos = [275 35 35 25]; % position of quiver scale text box, pixels
CurlButtonPos = [315 42 15 15]; % position of curl checkbox, pixels
CurlLabelPos = [330 35 40 27]; % position of overlay label, pixels
OverlayButtonPos = [370 42 15 15]; % position of overlay checkbox, pixels
OverlayLabelPos = [385 35 50 27]; % position of overlay label, pixels
OutlineButtonPos = [300 42 15 15]; % position of overlay checkbox, pixels
OutlineLabelPos = [315 35 50 27]; % position of overlay label, pixels
AxesButtonPos = [440 42 15 15]; % position of axes checkbox, pixels
AxesLabelPos = [455 35 35 27]; % position of axes label, pixels
PlayPauseButtonPos = [495 35 60 25]; % position of play button, pixels
Ch1TransPos = [2 35 50 25];
Ch2TransPos = [54 35 50 25];
Ch3TransPos = [106 35 50 25];
TransLabelPos = [158 35 100 25];
Ch1MultTextPos = [2 62 50 25]; % position of ch1mult text box, pixels
Ch2MultTextPos = [54 62 50 25]; % position of ch2mult text box, pixels
Ch3MultTextPos = [106 62 50 25]; % position of ch3mult text box, pixels
MultLabelPos = [46 40 66 20];
Ch1Slider1ValuePos = [2 89 50 25];
Ch2Slider1ValuePos = [54 89 50 25];
Ch3Slider1ValuePos = [106 89 50 25];
Ch1SliderPos = [17 116 20 125];
Ch2SliderPos = [69 116 20 125];
Ch3SliderPos = [121 116 20 125];
Ch1Slider2ValuePos = [2 243 50 25];
Ch2Slider2ValuePos = [54 243 50 25];
Ch3Slider2ValuePos = [106 243 50 25];
Ch1LabelPos = [2 268 50 30];
Ch2LabelPos = [54 268 50 30];
Ch3LabelPos = [106 268 50 30];
LinkSliderCheckboxPos = [2 300 15 15];
LinkLabelPos = [20 300 50 15];

MovieDialogSize = [400 220];
MoviePromptPos = [10 170 380 40];
MovieSliderPos = [10 140 380 20];
MovieStartTimeLabelPos = [10 110 90 20];
MovieStartTimeValuePos = [102 110 50 20];
MovieEndTimeLabelPos = [248 110 90 20];
MovieEndTimeValuePos = [340 110 50 20];
DownsampleLabelPos = [10 80 160 20];
DownsampleValuePos = [172 80 40 20];
NewFrameRateLabelPos = [215 80 103 20];
NewFrameRateValuePos = [320 80 70 20];
FileTypeLabelPos = [10 50 150 20];
FileTypeMenuPos = [162 50 228 20];
DurationPos = [10 10 220 20];
MovieCancelButtonPos = [260 10 60 30];
MovieOKButtonPos = [330 10 60 30];

% -=- Other Constants -=-
BackgroundColor = 0.94*[1 1 1];
FontSize = 8; % points
MaxQuiverCount = 30; % in each direction
qs=0.9; % initial scaling factor for quivers

% -=- Defaults -=-
x_default = 1; % grid size 1
y_default = 1; % grid size 1
t_default = 1; % time step 1
FrameRate_default = 30; % frames per second
DownSample = 1; % downsample factor, for saving
FileType_default = 'Motion JPEG AVI'; % for saving a movie
ch1mult=1;
ch2mult=1;
ch3mult=1;
trans(1)=0.5;
trans(2)=0.5;
trans(3)=0.5;
linkvar=0;

% -=- Parse inputs and get set up -=-
assert(nargin>0,['Usage: hf = ' mfilename '([x,y,t],im,[seg],[plt_coord],[MovName])'])
args = varargin;
if iscell(args{1}) || ~isvector(args{1}) % looks like (x,y,t) not provided
    im=args{1};
    if numel(args)>1
        for aa=2:numel(args)
            if iscell(args{aa})
                seg=args{aa};
            elseif isa(args{aa},'double')
                plt_coord=args{aa};
            else
                MovName=args{aa};
            end
        end
    end
else % looks like (x,y,t) provided
    x=args{1};
    y=args{2};
    t=args{3};
    im=args{4};
    if numel(args)>4
        for aa=5:numel(args)
            if iscell(args{aa})
                seg=args{aa};
            elseif isa(args{aa},'double') && ~isempty(args{aa})
                plt_coord=args{aa};
            else
                MovName=args{aa};
            end
        end
    end
end % if iscell(args{1}) || ~isvector(args{1})
if ~exist('MovName','var') || isempty(MovName)
    MovName='';
end
isVector=false;
isRGB=false;
if iscell(im)
    Nc=numel(im);
    for jj=1:Nc
        im{jj}=double(im{jj});
    end
    switch Nc
        case 1 % scalar
            im=im{1};
        case 2 % vector components
            isVector=true;
            imComps=im;
            mag=sqrt(im{1}.^2+im{2}.^2); % magnitude
            im=mag;
        case 3 % RGB image
            isRGB=true;
            imComps=im;
            im=imComps{1}; % just to get sz
            AxesLeftPos=AxesLeftPosRGB;
        otherwise
            error('Sorry, im must have 3 or fewer components.')
    end % switch Nc
else % if iscell(im)
    im=double(im);
end % if iscell(im)
sz=size(im);
if numel(sz)==2
    sz(3)=1;
elseif numel(sz)>3 || numel(sz)<2
    error('Sorry, only 2D and 3D arrays are supported.')
end
if ~exist('x','var') || isempty(x)
    x=x_default;
end
if ~exist('y','var') || isempty(y)
    y=y_default;
end
if ~exist('t','var') || isempty(t)
    t=t_default;
else
    FrameRate_default = 1/mean(diff(t));
end
if numel(x)==1
    x = (0:sz(2)-1)*x; % coordinates, meshgrid-style
elseif numel(x)~=sz(2)
    error('Size of x does not match size of im.')
end
if numel(y)==1
    y = (0:sz(1)-1)*y; % coordinates, meshgrid-style
elseif numel(y)~=sz(1)
    error('Size of y does not match size of im.')
end
if numel(t)==1
    t = (0:sz(3)-1)*t;
elseif numel(t)~=sz(3)
    error('Size of t does not match size of im.')
end
if isVector
    if sz(3)>1
        [~,dudy]=gradient(imComps{1},x,y,t);
        [dvdx,~]=gradient(imComps{2},x,y,t);
    else
        [~,dudy]=gradient(imComps{1},x,y);
        [dvdx,~]=gradient(imComps{2},x,y);
    end
    crl=dudy-dvdx; % curl
end
ind = isinf(im(:)) | isnan(im(:));
if any(ind)
    warning('Setting infinities and NaNs to zero.')
    im(ind)=0;
    if isRGB || isVector
        imComps{1}( isinf(imComps{1}) | isnan(imComps{1}) ) = 0;
        imComps{2}( isinf(imComps{2}) | isnan(imComps{2}) ) = 0;
    end
    if isRGB
        imComps{3}( isinf(imComps{3}) | isnan(imComps{3}) ) = 0;
    end
    if isVector
        mag(ind)=0;
        crl(ind)=0;
     end
end
if isRGB
    ch1m=min(imComps{1}(:)); % in case the user wants RGB mode
    ch2m=min(imComps{2}(:));
    ch3m=min(imComps{3}(:));
    ch1r=range(imComps{1}(:));
    ch2r=range(imComps{2}(:));
    ch3r=range(imComps{3}(:));
    if ch1r==0  % avoid NaNs later
        ch1r=1; 
    end
    if ch2r==0  % avoid NaNs later
        ch2r=1; 
    end
    if ch3r==0  % avoid NaNs later
        ch3r=1; 
    end
    clim_ch1=[ch1m ch1m+ch1r];
    clim_ch2=[ch2m ch2m+ch2r];
    clim_ch3=[ch3m ch3m+ch3r];
    
end
if isVector
    qxind=1:ceil(sz(2)/MaxQuiverCount):sz(2); % x indices of quivers
    qyind=1:ceil(sz(1)/MaxQuiverCount):sz(1); % y indices of quivers
    maxlen=max(sqrt( ...
        reshape(imComps{1}(qyind,qxind,:),[],1).^2 ...
        + reshape(imComps{2}(qyind,qxind,:),[],1).^2 )) ... % max magnitude
        / sqrt( mean(diff(x)).^2 + mean(diff(y)).^2 ); % diagonal grid size
end
indt=1; % first slice
t1=[]; % initialize for later

% -=- Adding dummy data to seg if it is not 3 channels -=-
if exist('seg','var')
    if numel(seg)==2
        seg{3}=zeros(size(seg{1}));
    end
end

% -=- Creating segmentation outline data -=-
if isRGB && exist('seg','var')
    for ch=1:size(imComps,2) %normalizing imComps
        im_outline{ch}=imComps{ch}./max(imComps{ch}(:));
    end
    for ch=1:size(seg,2) % creates the outlines using segmentation
        for ss=1:size(seg{ch},3)
            seg_outline{ch}(:,:,ss)=bwperim(seg{ch}(:,:,ss));
        end
    end
    for ch=1:size(seg,2) %makes outline red (1), green (2), or blue (3)
        im_outline{ch}=im_outline{ch}.*~seg_outline{ch}+seg_outline{ch}; 
    end
end %if isRGB

% -=- Plot everything for the first time -=-
hff=figure;
mypos=get(hff,'position');
mypos(3:4)=FigSize;
set(hff,'position',mypos,'name',[num2str(t(indt),'%0.2f') ...
    ' (' num2str(indt) ' of ' num2str(sz(3)) ')'], ...
    'color',BackgroundColor,'SizeChangedFcn',@SizeChanged);
h=axes;
if isRGB
    im1=cat(3, ...
        (squeeze(imComps{1}(:,:,indt))-ch1m)*ch1mult/ch1r, ... % rgb channels, with normalized intensity
        (squeeze(imComps{2}(:,:,indt))-ch2m)*ch2mult/ch2r, ...
        (squeeze(imComps{3}(:,:,indt))-ch3m)*ch3mult/ch3r );
    im1(im1<0)=0;
else
    im1=squeeze(im(:,:,indt));
end
if exist('seg','var') %plots segmentations
    isoutline=0; %plots seg overlay first
    if ~isoutline %seg overlay
        for n=1:numel(seg)
            colormapinput=[0,0,0];colormapinput(n)=1;
            im1=labeloverlay(im1,seg{n}(:,:,indt),'ColorMap',colormapinput,'Trans', trans(n));
        end
    else %seg outline
        for n=1:numel(seg)
            colormapinput=[0,0,0];colormapinput(n)=1;
            im1=labeloverlay(im1,seg_outline{n}(:,:,indt),'ColorMap',colormapinput,'Trans', trans(n));
        end
    end
end %if exist('seg','var')

hi=imagesc(x,y,im1);
set(h,'nextplot','add','units','pixels', ...
    'ydir','reverse','clim',[min(im(:)) max(im(:))], ...
    'outerposition',[AxesLeftPos AxesBottomPos FigSize(1)-AxesLeftPos ...
    FigSize(2)-AxesBottomPos])
if exist('plt_coord','var')
    hp=plot(plt_coord(:,1,indt),plt_coord(:,2,indt),'.y'); % Yellow to blue
end
if isVector % plot quiver arrows for vector field
    [q1x,q1y]=meshgrid(x(qxind),y(qyind));
    hq=quiver(q1x(:),q1y(:), ...
        reshape(imComps{1}(qyind,qxind,indt),[],1)/maxlen*qs, ...
        reshape(imComps{2}(qyind,qxind,indt),[],1)/maxlen*qs, ...
        'k','autoscale','off');
end
ht=title(num2str(t(indt),'%0.2f'));
ht(2)=text(0,0,num2str(t(indt),'%0.2f'));
set(ht(2),'units','points','position',[5 5 0], ...
    'verticalalignment','bottom','horizontalalignment','left', ...
    'fontsize',FontSize,'fontweight','bold','backgroundcolor','k', ...
    'color','w','visible','off')
if sz(3)<=1 % data includes 1 slice or less
    set(ht(1),'visible','off')
end
if ~isRGB
    hcb=colorbar;
    hcb.Label.String='magnitude';
end

% -=- Make controls -=-
ha = annotation('rectangle',[0 0 1 1]);
set(ha,'units','pixels','position',[AxesLeftPos AxesBottomPos ...
    FigSize(1)-AxesLeftPos FigSize(2)-AxesBottomPos],'color','r', ...
    'linewidth',5,'facecolor','k','facealpha',0.2,'visible','off');

if isRGB
    Ch1MultText = uicontrol('Style','edit','units','pixels', ...
        'Position',Ch1MultTextPos,'String',num2str(ch1mult), ...
        'fontsize',FontSize,'Callback',@inputFromch1mult);
    Ch2MultText = uicontrol('Style','edit','units','pixels', ...
        'Position',Ch2MultTextPos,'String',num2str(ch2mult), ...
        'fontsize',FontSize,'Callback',@inputFromch2mult);
    Ch3MultText = uicontrol('Style','edit','units','pixels', ...
        'Position',Ch3MultTextPos,'String',num2str(ch3mult), ...
        'fontsize',FontSize,'Callback',@inputFromch3mult);
    uicontrol('Style','text','units','pixels', ... % multiplier label
        'Position',MultLabelPos,'FontSize',FontSize, ...
        'FontWeight','bold','String','multipliers', ...
        'ForegroundColor','k','horizontalalignment','center', ...
        'backgroundcolor',BackgroundColor);
    Ch1Slider=superSlider(gcf,'position', ...
        Ch1SliderPos,'value', [0 1], 'numslides',2,'max',1,'min',0, ...
        'stepsize',.05,'Callback',@inputFromCh1Slider);
    set(Ch1Slider,'units','pixels')
    set(Ch1Slider,'position',Ch1SliderPos)
    Ch2Slider=superSlider(gcf,'position', ...
        Ch2SliderPos,'value', [0 1], 'numslides',2,'max',1,'min',0, ...
        'stepsize',.05,'Callback',@inputFromCh2Slider);
    set(Ch2Slider,'units','pixels')
    set(Ch2Slider,'position',Ch2SliderPos)
    Ch3Slider=superSlider(gcf,'position', ...
        Ch3SliderPos,'value', [0 1], 'numslides',2,'max',1,'min',0, ...
        'stepsize',.05,'Callback',@inputFromCh3Slider);
    set(Ch3Slider,'units','pixels')
    set(Ch3Slider,'position',Ch3SliderPos)
    Ch1Slider1Value = uicontrol('Style','edit','units','pixels', ...
        'Position',Ch1Slider1ValuePos,'String',num2str(clim_ch1(1)), ...
        'fontsize',FontSize,'Callback',@inputFromCh1Slider1Value); 
    Ch2Slider1Value = uicontrol('Style','edit','units','pixels', ...
        'Position',Ch2Slider1ValuePos,'String',num2str(clim_ch2(1)'), ...
        'fontsize',FontSize,'Callback',@inputFromCh2Slider1Value);
    Ch3Slider1Value = uicontrol('Style','edit','units','pixels', ...
        'Position',Ch3Slider1ValuePos,'String',num2str(clim_ch3(1)), ...
        'fontsize',FontSize,'Callback',@inputFromCh3Slider1Value);
    Ch1Slider2Value = uicontrol('Style','edit','units','pixels', ...
        'Position',Ch1Slider2ValuePos,'String',num2str(clim_ch1(2)), ...
        'fontsize',FontSize,'Callback',@inputFromCh1Slider2Value);
    Ch2Slider2Value = uicontrol('Style','edit','units','pixels', ...
        'Position',Ch2Slider2ValuePos,'String',num2str(clim_ch2(2)), ...
        'fontsize',FontSize,'Callback',@inputFromCh2Slider2Value);
    Ch3Slider2Value = uicontrol('Style','edit','units','pixels', ...
        'Position',Ch3Slider2ValuePos,'String',num2str(clim_ch3(2)), ...
        'fontsize',FontSize,'Callback',@inputFromCh3Slider2Value);
    uicontrol('Style','text','units','pixels', ... % Channel 1 label
        'Position',Ch1LabelPos,'FontSize',FontSize,'FontWeight','bold', ...
        'String','Ch 1 (red)','ForegroundColor','k', ...
        'horizontalalignment','center','backgroundcolor',BackgroundColor);
    uicontrol('Style','text','units','pixels', ... % Channel 1 label
        'Position',Ch2LabelPos,'FontSize',FontSize,'FontWeight','bold', ...
        'String','Ch 2 (green)','ForegroundColor','k', ...
        'horizontalalignment','center','backgroundcolor',BackgroundColor);
    uicontrol('Style','text','units','pixels', ... % Channel 1 label
        'Position',Ch3LabelPos,'FontSize',FontSize,'FontWeight','bold', ...
        'String','Ch 3 (blue)','ForegroundColor','k', ...
        'horizontalalignment','center','backgroundcolor',BackgroundColor);
    LinkSlidersCheckbox = uicontrol('style','checkbox', ...
            'units','pixels','position',LinkSliderCheckboxPos, ...
            'string','','value',0, ...
            'backgroundcolor',BackgroundColor,'callback',@inputFromLinkSliderCheckbox);
    uicontrol('Style','text','units','pixels', ... % link label
            'Position',LinkLabelPos,'FontSize',FontSize,'FontWeight','bold', ...
            'String',{'link sliders','dont link'},'ForegroundColor','k', ...
            'horizontalalignment','left','backgroundcolor',BackgroundColor);
    if exist('seg','var') 
        Ch1TransText = uicontrol('Style','edit','units','pixels', ...
            'Position',Ch1TransPos,'String',num2str(trans(1)), ...
            'fontsize',FontSize,'Callback',@inputFromCh1TransText);
        Ch2TransText = uicontrol('Style','edit','units','pixels', ...
            'Position',Ch2TransPos,'String',(trans(2)), ...
            'fontsize',FontSize,'Callback',@inputFromCh2TransText);
        Ch3TransText = uicontrol('Style','edit','units','pixels', ...
            'Position',Ch3TransPos,'String',(trans(3)), ...
            'fontsize',FontSize,'Callback',@inputFromCh3TransText);
        uicontrol('Style','text','units','pixels', ... % transparency label
            'Position',TransLabelPos,'FontSize',FontSize, ...
            'FontWeight','bold','String','seg transparency', ...
            'ForegroundColor','k','horizontalalignment','left', ...
            'backgroundcolor',BackgroundColor);
        OutlineButton = uicontrol('style','checkbox', ...
            'units','pixels','position',OutlineButtonPos, ...
            'string','','value',0, ...
            'backgroundcolor',BackgroundColor,'callback',@OutlineToggle);
        uicontrol('Style','text','units','pixels', ... % outline label
            'Position',OutlineLabelPos,'FontSize',FontSize,'FontWeight','bold', ...
            'String',{'seg','outline'},'ForegroundColor','k', ...
            'horizontalalignment','left','backgroundcolor',BackgroundColor);
    end
end
if sz(3)>1 % data includes more than 1 slice
    OverlayButton = uicontrol('style','checkbox', ...
        'units','pixels','position',OverlayButtonPos, ...
        'string','','value',0, ...
        'backgroundcolor',BackgroundColor,'callback',@OverlayToggle);
    uicontrol('Style','text','units','pixels', ... % overlay label
        'Position',OverlayLabelPos,'FontSize',FontSize,'FontWeight','bold', ...
        'String','overlay label','ForegroundColor','k', ...
        'horizontalalignment','left','backgroundcolor',BackgroundColor);
    uicontrol('Style','text','units','pixels', ... % Slice label
        'Position',SliceLabelPos,'FontSize',FontSize,'FontWeight','bold', ...
        'String','slice:','ForegroundColor','k', ...
        'horizontalalignment','right','backgroundcolor',BackgroundColor);
    SliceText = uicontrol('Style','edit','units','pixels', ...
        'Position',SliceTextPos,'String',num2str(t(1),'%.4g'), ...
        'fontsize',FontSize,'Callback',@inputFromSliceText);
    SliceSlider = uicontrol('Style','slider','units','pixels', ...
        'Position',SliceSliderPos,'Min',t(1),'Max',t(sz(3)), ...
        'Value',t(1),'SliderStep',[1/sz(3) 1/sz(3)+eps], ...
        'Callback',@inputFromSliceSlider);
    uicontrol('Style','text','units','pixels', ... % Slice unit label
        'Position',SliceUnitLabelPos,'FontSize',FontSize, ...
        'FontWeight','bold','String','slice unit:','ForegroundColor','k', ...
        'horizontalalignment','right','backgroundcolor',BackgroundColor);
    SliceUnitText = uicontrol('Style','edit','units','pixels', ...
        'Position',SliceUnitTextPos,'String','', ...
        'fontsize',FontSize,'Callback',@inputFromSliceUnitText);
    uicontrol('Style','text','units','pixels', ... % Movie name label
        'Position',MovNameLabelPos,'FontSize',FontSize, ...
        'FontWeight','bold','String','movie name:', ...
        'ForegroundColor','k','horizontalalignment','right', ...
        'backgroundcolor',BackgroundColor);
    MovNameText = uicontrol('Style','edit','units','pixels', ...
        'fontsize',FontSize,'Position',MovNameTextPos,'String',MovName, ...
        'callback',@inputFromMovNameText);
    uicontrol('Style','text','units','pixels', ... % Frame rate label
        'Position',FrameRateLabelPos,'FontSize',FontSize, ...
        'FontWeight','bold','String','fps:','ForegroundColor','k', ...
        'horizontalalignment','right','backgroundcolor',BackgroundColor);
    FrameRateText = uicontrol('Style','edit','units','pixels', ...
        'fontsize',FontSize,'Position',FrameRateTextPos, ...
        'String',FrameRate_default);
    PlayPauseButton = uicontrol('style','pushbutton', ...
        'units','pixels','position',PlayPauseButtonPos, ...
        'string','play','callback',@PlayPause);
    if ~isempty(get(MovNameText,'string'))
        set(PlayPauseButton,'string','save...')
    end
    AxesButton = uicontrol('style','checkbox', ...
        'units','pixels','position',AxesButtonPos, ...
        'string','save axes only','value',0, ...
        'backgroundcolor',BackgroundColor,'callback',@AxesToggle);
    uicontrol('Style','text','units','pixels', ... % axes button label
        'Position',AxesLabelPos,'FontSize',FontSize,'FontWeight','bold', ...
        'String','axes only','ForegroundColor','k', ...
        'horizontalalignment','left','backgroundcolor',BackgroundColor);
else % if sz(3)>1
    AxesBottomPos = AxesBottomPosOneSlice;
    SizeChanged();
end % if sz(3)>1
if isVector
    uicontrol('Style','text','units','pixels', ...
        'Position',QuiverScaleLabelPos,'FontSize',FontSize,'FontWeight','bold', ...
        'String','quiver scale: ','ForegroundColor','k', ...
        'horizontalalignment','right','backgroundcolor',BackgroundColor);
    QuiverScaleText = uicontrol('Style','edit','units','pixels', ...
        'Position',QuiverScaleTextPos,'String',num2str(qs,'%.3g'), ...
        'fontsize',FontSize,'Callback',@inputFromQuiverScaleText);
    CurlButton = uicontrol('style','checkbox', ...
        'units','pixels','position',CurlButtonPos, ...
        'string','','value',0, ...
        'backgroundcolor',BackgroundColor,'callback',@CurlToggle);
    uicontrol('Style','text','units','pixels', ... % curl button label
        'Position',CurlLabelPos,'FontSize',FontSize,'FontWeight','bold', ...
        'String','show curl','ForegroundColor','k', ...
        'horizontalalignment','left','backgroundcolor',BackgroundColor);
end % if isVector
if nargout>0
    hf=hff;
end

% -=- Callback: When the user resizes the window -=-
function SizeChanged(~,~)
    old_ax_units = get(h,'units');
    old_fig_units = get(hff,'units');
    set(hff,'units','pixels');
    hf_pos = get(hff,'position');
    set(h,'units','pixels','outerposition',[AxesLeftPos AxesBottomPos ...
        hf_pos(3)-AxesLeftPos hf_pos(4)-AxesBottomPos]);
    set(h,'units',old_ax_units);
    set(hff,'units',old_fig_units);
end % function SizeChanged(~,~)

% -=- Callback: When the user changes time unit -=-
    function inputFromSliceUnitText(~,~)
        t1=get(SliceSlider,'value');
        changeSlice();
    end % function inputFromSliceUnitText(~,~)

% -=- Callback: When the user changes mult -=-
    function inputFromch1mult(~,~)
        ch1multtmp=str2double(get(Ch1MultText,'string'));
        if isnan(ch1multtmp) % user entered a non-number; replace with current value
            set(Ch1MultText,'string',num2str(ch1mult));
            return
        end
        ch1mult=ch1multtmp;
        updateImagesAndQuivers()
    end 
    function inputFromch2mult(~,~)
        ch2multtmp=str2double(get(Ch2MultText,'string'));
        if isnan(ch2multtmp) % user entered a non-number; replace with current value
            set(Ch2MultText,'string',num2str(ch2mult));
            return
        end
        ch2mult=ch2multtmp;
        updateImagesAndQuivers()
    end 
    function inputFromch3mult(~,~)
        ch3multtmp=str2double(get(Ch3MultText,'string'));
        if isnan(ch3multtmp) % user entered a non-number; replace with current value
            set(Ch3MultText,'string',num2str(ch3mult));
            return
        end
        ch3mult=ch3multtmp;
        updateImagesAndQuivers()
    end 

% -=- Callback: When the user changes trans -=-
    function inputFromCh1TransText(~,~)
        trans1tmp=str2double(get(Ch1TransText,'string'));
        if isnan(trans1tmp) % user entered a non-number; replace with current value
            set(Ch1TransText,'string',num2str(trans(1)));
            return
        end
        trans(1)=trans1tmp;
        updateImagesAndQuivers()
    end 
    function inputFromCh2TransText(~,~)
        trans2tmp=str2double(get(Ch2TransText,'string'));
        if isnan(trans2tmp) % user entered a non-number; replace with current value
            set(Ch2TransText,'string',num2str(trans(2)));
            return
        end
        trans(2)=trans2tmp;
        updateImagesAndQuivers()
    end 
    function inputFromCh3TransText(~,~)
        trans3tmp=str2double(get(Ch3TransText,'string'));
        if isnan(trans3tmp) % user entered a non-number; replace with current value
            set(Ch3TransText,'string',num2str(trans(3)));
            return
        end
        trans(3)=trans3tmp;
        updateImagesAndQuivers()
    end 

% -=- Callback: When the user changes Channel limit Sliders -=-
    function inputFromCh1Slider(~,~)
        ch1lim=get(Ch1Slider, 'UserData');
        clim_ch1=ch1lim(1,:)*ch1r+ch1m;  
        set(Ch1Slider1Value,'String', num2str(clim_ch1(1),'%0.5g'))
        set(Ch1Slider2Value,'String', num2str(clim_ch1(2),'%0.5g'))
        if linkvar
            clim_ch1(1)=clim_ch1(2);
            set(Ch1Slider1Value,'String', num2str(clim_ch1(1),'%0.5g'))
            set(Ch1Slider2Value,'String', num2str(clim_ch1(2),'%0.5g'))
            % get rid of slider 1 if sliders are linked
            allSlides = get(Ch1Slider, 'Children');
            location = get(allSlides(1), 'Position');
            location(2)=(-.5); 
            set(allSlides(1),'Position',location)
        end
        updateImagesAndQuivers() 
    end 
    function inputFromCh2Slider(~,~)
        ch2lim=get(Ch2Slider, 'UserData');
        clim_ch2=ch2lim(1,:)*ch2r+ch2m;
        set(Ch2Slider1Value,'String', num2str(clim_ch2(1),'%0.5g'))
        set(Ch2Slider2Value,'String', num2str(clim_ch2(2),'%0.5g'))
        if linkvar
            clim_ch2(1)=clim_ch2(2);
            set(Ch2Slider1Value,'String', num2str(clim_ch2(1),'%0.5g'))
            set(Ch2Slider2Value,'String', num2str(clim_ch2(2),'%0.5g'))
            % get rid of slider 1 
            allSlides = get(Ch2Slider, 'Children');
            location = get(allSlides(1), 'Position');
            location(2)=(-.5); 
            set(allSlides(1),'Position',location)
        end
        updateImagesAndQuivers()
    end 
    function inputFromCh3Slider(~,~)
        ch3lim=get(Ch3Slider, 'UserData');
        clim_ch3=ch3lim(1,:)*ch3r+ch3m;
        set(Ch3Slider1Value,'String', num2str(clim_ch3(1),'%0.5g'))
        set(Ch3Slider2Value,'String', num2str(clim_ch3(2),'%0.5g'))
        if linkvar
            clim_ch3(1)=clim_ch3(2);
            set(Ch3Slider1Value,'String', num2str(clim_ch3(1),'%0.5g'))
            set(Ch3Slider2Value,'String', num2str(clim_ch3(2),'%0.5g'))
            % get rid of slider 1 
            allSlides = get(Ch3Slider, 'Children');
            location = get(allSlides(1), 'Position');
            location(2)=(-.5); 
            set(allSlides(1),'Position',location)
        end
        updateImagesAndQuivers()
    end 

% -=- Callback: When the user changes Channel limit Value Boxes -=-
    function inputFromCh1Slider1Value(~,~)
        limtmp=str2double(get(Ch1Slider1Value,'string')); 
        if isnan(limtmp)
            set(Ch1Slider1Value,'string',num2str(clim_ch1(1),'%0.5g'))
            return
        end
        clim_ch1(1)=limtmp;
        if linkvar
           clim_ch1(2)=clim_ch1(1);
           set(Ch1Slider2Value,'string',num2str(clim_ch1(2),'%0.5g'))
        end
        updateImagesAndQuivers()
        
        % update slider location
        allSlides = get(Ch1Slider, 'Children');
        location = get(allSlides(1), 'Position');
        location(2)=((clim_ch1(1)-ch1m)/ch1r)*(1-2*location(4)); % scale since the maximum location is 1-2*location(4) 
        set(allSlides(1),'Position',location)   
        if linkvar
           location = get(allSlides(2), 'Position');
           location(2)=((clim_ch1(2)-ch1m)/ch1r)*(1-2*location(4)); % scale since the maximum location is 1-2*location(4) 
           set(allSlides(2),'Position',location)  
        end
    end 
    function inputFromCh2Slider1Value(~,~)
        limtmp=str2double(get(Ch2Slider1Value,'string')); 
        if isnan(limtmp)
            set(Ch2Slider1Value,'string',num2str(clim_ch2(1),'%0.5g'))
            return
        end
        clim_ch2(1)=limtmp;
        if linkvar
           clim_ch2(2)=clim_ch2(1);
           set(Ch2Slider2Value,'string',num2str(clim_ch2(2),'%0.5g'))
        end
        updateImagesAndQuivers() 
        
        % update slider location
        allSlides = get(Ch2Slider, 'Children');
        location = get(allSlides(1), 'Position');
        location(2)=((clim_ch2(1)-ch2m)/ch2r)*(1-2*location(4)); % scale since the maximum location is 1-2*location(4) 
        set(allSlides(1),'Position',location) 
        if linkvar
           location = get(allSlides(2), 'Position');
           location(2)=((clim_ch2(2)-ch2m)/ch2r)*(1-2*location(4)); % scale since the maximum location is 1-2*location(4) 
           set(allSlides(2),'Position',location)  
        end
    end 
    function inputFromCh3Slider1Value(~,~)
        limtmp=str2double(get(Ch3Slider1Value,'string')); 
        if isnan(limtmp)
            set(Ch3Slider1Value,'string',num2str(clim_ch3(1),'%0.5g'))
            return
        end
        clim_ch3(1)=limtmp;
        if linkvar
           clim_ch3(2)=clim_ch3(1);
           set(Ch3Slider2Value,'string',num2str(clim_ch3(2),'%0.5g'))
        end
        updateImagesAndQuivers()
                        
        % update slider location
        allSlides = get(Ch3Slider, 'Children');
        location = get(allSlides(1), 'Position');
        location(2)=((clim_ch3(1)-ch3m)/ch3r)*(1-2*location(4)); % scale since the maximum location is 1-2*location(4) 
        set(allSlides(1),'Position',location) 
        if linkvar
           location = get(allSlides(2), 'Position');
           location(2)=((clim_ch3(2)-ch3m)/ch3r)*(1-2*location(4)); % scale since the maximum location is 1-2*location(4) 
           set(allSlides(2),'Position',location)  
        end
    end 
    function inputFromCh1Slider2Value(~,~)
        limtmp=str2double(get(Ch1Slider2Value,'string')); 
        if isnan(limtmp)
            set(Ch1Slider2Value,'string',num2str(clim_ch1(2),'%0.5g'))
            return
        end
        clim_ch1(2)=limtmp;
        if linkvar
           clim_ch1(1)=clim_ch1(2);
           set(Ch1Slider1Value,'string',num2str(clim_ch1(1),'%0.5g'))
        end
        updateImagesAndQuivers()
        
        % update slider location
        allSlides = get(Ch1Slider, 'Children');
        location = get(allSlides(2), 'Position');
        location(2)=((clim_ch1(2)-ch1m)/ch1r)*(1-2*location(4))+location(4); % scale since the minimum location is location(4) and the max is 1-location(4)
        set(allSlides(2),'Position',location)
        if linkvar
            location = get(allSlides(1), 'Position');
            location(2)=((clim_ch1(1)-ch1m)/ch1r)*(1-2*location(4))+location(4); % scale since the minimum location is location(4) and the max is 1-location(4)
            set(allSlides(1),'Position',location)
        end
    end 
    function inputFromCh2Slider2Value(~,~)
        limtmp=str2double(get(Ch2Slider2Value,'string')); 
        if isnan(limtmp)
            set(Ch2Slider2Value,'string',num2str(clim_ch2(2),'%0.5g'))
            return
        end
        clim_ch2(2)=limtmp;
        if linkvar
           clim_ch2(1)=clim_ch2(2);
           set(Ch2Slider1Value,'string',num2str(clim_ch2(1),'%0.5g'))
        end
        updateImagesAndQuivers()
        
        % update slider location
        allSlides = get(Ch2Slider, 'Children');
        location = get(allSlides(2), 'Position');
        location(2)=((clim_ch2(2)-ch2m)/ch2r)*(1-2*location(4))+location(4); % scale since the minimum location is location(4) and the max is 1-location(4)
        set(allSlides(2),'Position',location)
        if linkvar
            location = get(allSlides(1), 'Position');
            location(2)=((clim_ch2(1)-ch2m)/ch2r)*(1-2*location(4))+location(4); % scale since the minimum location is location(4) and the max is 1-location(4)
            set(allSlides(1),'Position',location)
        end
    end 
    function inputFromCh3Slider2Value(~,~)
        limtmp=str2double(get(Ch3Slider2Value,'string')); 
        if isnan(limtmp)
            set(Ch3Slider2Value,'string',num2str(clim_ch3(2),'%0.5g'))
            return
        end
        clim_ch3(2)=limtmp;
        if linkvar
           clim_ch3(1)=clim_ch3(2);
           set(Ch3Slider1Value,'string',num2str(clim_ch3(1),'%0.5g'))
        end
        updateImagesAndQuivers()
        
        % update slider location
        allSlides = get(Ch3Slider, 'Children');
        location = get(allSlides(2), 'Position');
        location(2)=((clim_ch3(2)-ch3m)/ch3r)*(1-2*location(4))+location(4); % scale since the minimum location is location(4) and the max is 1-location(4)
        set(allSlides(2),'Position',location)
        if linkvar
            location = get(allSlides(1), 'Position');
            location(2)=((clim_ch3(1)-ch3m)/ch3r)*(1-2*location(4))+location(4); % scale since the minimum location is location(4) and the max is 1-location(4)
            set(allSlides(1),'Position',location)
        end
    end 
% -=- Callback: When user links or unlinks channel sliders
    function inputFromLinkSliderCheckbox(~,~)
        linkvar=get(LinkSlidersCheckbox,'value');
    end % function inputFromSliceSlider(~,~)
% -=- Callback: When the user changes time using the slider bar -=-
    function inputFromSliceSlider(~,~)
        t1=get(SliceSlider,'value');
        changeSlice();
    end % function inputFromSliceSlider(~,~)

% -=- Callback: When the user changes time using the text box -=-
    function inputFromSliceText(~,~)
        t1=str2double(get(SliceText,'string'));
        if isnan(t1) % user entered a non-number; replace with current time
            set(SliceText,'string',num2str(t(indt),'%.4g'));
            return
        end
        changeSlice();
    end % function inputFromSliceText(~,~)

% -=- Callback: When the user changes quiver scale using the text box -=-
    function inputFromQuiverScaleText(~,~)
        qs1=str2double(get(QuiverScaleText,'string'));
        if isnan(qs1) % user entered a non-number; replace with current quiver scale
            set(QuiverScaleText,'string',num2str(qs,'%.3g'));
        else
            qs=qs1;
            updateImagesAndQuivers();
        end
    end % function inputFromQuiverScaleText(~,~)

% -=- Callback: When user clicks PlayPauseButton -=-
    function PlayPause(~,~)
        if strcmpi(get(PlayPauseButton,'string'),'play') || ...
            strcmpi(get(PlayPauseButton,'string'),'save...') % could implement "Save" as a separate function, but setup and cleanup duplicate "Play"
            Play()
        else
            Pause()
        end
    end % function PlayPause(~,~)

% -=- Callback: When the user toggles overlay checkbox -=-
    function OverlayToggle(~,~)
        if get(OverlayButton,'value')
            set(ht(1),'visible','off')
            set(ht(2),'visible','on')
        else
            set(ht(1),'visible','on')
            set(ht(2),'visible','off')
        end
    end % function OverlayToggle(~,~)

% -=- Callback: When the user toggles outline checkbox -=-
    function OutlineToggle(~,~)
        if get(OutlineButton,'value')
            isoutline=1;
            updateImagesAndQuivers
        else
            isoutline=0;
            updateImagesAndQuivers
        end
    end % function OverlayToggle(~,~)

% -=- Callback: When user changes movie name -=-
    function inputFromMovNameText(~,~)
        if isempty(get(MovNameText,'string'))
            set(PlayPauseButton,'string','play')
        else
            set(PlayPauseButton,'string','save...')
        end
    end

% -=- Callback: When the user toggles curl checkbox -=-
    function CurlToggle(~,~)
        if get(CurlButton,'value')
            im=crl;
            hcb.Label.String='curl';
            caxis auto
        else
            im=mag;
            hcb.Label.String='magnitude';
            caxis auto
        end
        updateImagesAndQuivers()
    end % function CurlToggle(~,~)

% -=- Callback: When the user toggles axes checkbox -=-
    function AxesToggle(~,~)
        if get(AxesButton,'value')
            li = GetLayoutInformation(h);
            set(ha,'position',li.PlotBox);
        else
            old_units = get(hff,'units');
            set(hff,'units','pixels');
            hf_pos = get(hff,'position');
            set(ha,'position',...
                [AxesLeftPos AxesBottomPos hf_pos(3)-AxesLeftPos ...
                hf_pos(4)-AxesBottomPos]);
            set(hff,'units',old_units);
        end
        set(ha,'visible','on')
        pause(1)
        set(ha,'visible','off')
    end % function AxesToggle(~,~)

% -=- Start playback -=-
    function Play(~)
        set([SliceText MovNameText FrameRateText OverlayButton ...
            SliceUnitText AxesButton],'enable','off'); % disable controls
        if isVector
            set([QuiverScaleText CurlButton],'enable','off'); % disable controls
        end
        set(SliceSlider,'callback',''); % disable slider while movie plays
        FrameRate=str2double(get(FrameRateText,'string'));
        if isnan(FrameRate)
            FrameRate = FrameRate_default;
            set(FrameRateText,'string',num2str(FrameRate))
        end
        MovName=get(MovNameText,'string');
        if ~isempty(MovName) % save playback to disk
            set(PlayPauseButton,'string','pause','enable','off')
            indt0=indt; % save the time for later
            [~,~,ext]=fileparts(MovName);
            if ~strcmpi(ext,'.avi')
                MovName=[MovName '.avi']; % append extension if necessary
            end
            old_units = get(hff,'units');
            MovieTimes = [t(1) t(end)];

            % -=- Make dialog box and controls -=-
            figPos = get(gcf,'position');
            MovieDialog = dialog('position', ...
                [figPos(1)+figPos(3)/2-MovieDialogSize(1)/2 ...
                    figPos(2)+figPos(4)/2-MovieDialogSize(2)/2 ...
                    MovieDialogSize], ...
                'name',['Saving ' MovName ' ...']);
            uicontrol(MovieDialog,'style','text','units','pixels', ...
                'position',MoviePromptPos, ...
                'foregroundcolor','k', ...
                'backgroundcolor',BackgroundColor,'string', ...
                'Please specify the times, downsampling, and frame rate.')
            MovieDialogSlider = superSlider(MovieDialog, ...
                'position',MovieSliderPos, ...
                'value',[t(1) t(end)],'numslides',2, ...
                'max',t(end),'min',t(1),'stepsize',range(t)/20, ...
                'callback',@inputFromMovieDialogSlider);
            set(MovieDialogSlider,'units','pixels')
            set(MovieDialogSlider,'position',MovieSliderPos)
            uicontrol(MovieDialog,'style','text','units','pixels', ...
                'position',MovieStartTimeLabelPos,'foregroundcolor','k', ...
                'backgroundcolor',BackgroundColor, ...
                'string','initial slice:', ...
                'horizontalalignment','right')
            MovieStartTimeValue = uicontrol(MovieDialog,'style','edit', ...
                'units','pixels','position',MovieStartTimeValuePos, ...
                'string',num2str(t(1),'%0.2f'), ...
                'callback',@inputFromMovieStartTimeValue);
            uicontrol(MovieDialog,'style','text','units','pixels', ...
                'position',MovieEndTimeLabelPos,'foregroundcolor','k', ...
                'backgroundcolor',BackgroundColor, ...
                'string','final slice:', ...
                'horizontalalignment','right')
            MovieEndTimeValue = uicontrol(MovieDialog,'style','edit', ...
                'units','pixels','position',MovieEndTimeValuePos, ...
                'string',num2str(t(end),'%0.2f'), ...
                'callback',@inputFromMovieEndTimeValue);
            uicontrol(MovieDialog,'style','text','units','pixels', ...
                'position',DownsampleLabelPos,'foregroundcolor','k', ...
                'backgroundcolor',BackgroundColor, ...
                'string','downsample factor:', ...
                'horizontalalignment','right')
            DownsampleValue = uicontrol(MovieDialog,'style','edit', ...
                'units','pixels','position',DownsampleValuePos, ...
                'string',num2str(DownSample), ...
                'callback',@inputFromDownsampleValue);
            uicontrol(MovieDialog,'style','text','units','pixels', ...
                'position',NewFrameRateLabelPos,'foregroundcolor','k', ...
                'backgroundcolor',BackgroundColor, ...
                'string','frame rate (Hz):', ...
                'horizontalalignment','right')
            NewFrameRateValue = uicontrol(MovieDialog,'style','edit', ...
                'units','pixels','position',NewFrameRateValuePos, ...
                'string',num2str(FrameRate,'%0.3f'), ...
                'callback',@inputFromNewFrameRateValue);
            uicontrol(MovieDialog,'style','text','units','pixels', ...
                'position',FileTypeLabelPos,'foregroundcolor','k', ...
                'backgroundcolor',BackgroundColor, ...
                'string','file type:','horizontalalignment','right')
            if strcmp(computer,'GLNXA64') % Linux / Unix
                FileTypeList = {'Motion JPEG AVI','Archival'};
            else % Max / Windows
                FileTypeList = {'Motion JPEG AVI','MPEG-4','Archival'};
            end
            FileType = FileType_default;
            FileTypeMenu = uicontrol(MovieDialog,'style','popupmenu', ...
                'units','pixels','position', ...
                FileTypeMenuPos,'foregroundcolor','k', ...
                'backgroundcolor',BackgroundColor,'string', ...
                FileTypeList,'horizontalalignment','left', ...
                'callback',@inputFromFileTypeMenu);
            DurationMessage = uicontrol(MovieDialog,'style','text', ...
                'units','pixels','position', ...
                DurationPos,'foregroundcolor','k', ...
                'backgroundcolor',BackgroundColor,'string', ...
                ['new duration: ' num2str(range(t),'%0.2f') ' s'], ...
                'horizontalalignment','left');
            uicontrol(MovieDialog,'style','pushbutton', ...
                'units','pixels','position',MovieCancelButtonPos, ...
                'string','Cancel','callback',@inputFromMovieCancelButton);
            uicontrol(MovieDialog,'style','pushbutton', ...
                'units','pixels','position',MovieOKButtonPos, ...
                'string','OK','callback',@inputFromMovieOKButton);
            waitfor(MovieDialog);

            set(hff,'units',old_units);
            t1=t(indt0);
            set([SliceText MovNameText FrameRateText OverlayButton ...
                SliceUnitText AxesButton PlayPauseButton],'enable','on'); % disable controls
            if isempty(get(MovNameText,'string'))
                set(PlayPauseButton,'string','play')
            else
                set(PlayPauseButton,'string','save...')
            end
            if isVector
                set([QuiverScaleText CurlButton],'enable','on'); % disable controls
            end
            set(SliceSlider,'callback',''); % disable slider while movie plays
            changeSlice() % go back to the same time

        else % if ~isempty(MovName) - just play, don't save
            set(PlayPauseButton,'string','pause')
            while strcmpi(get(PlayPauseButton,'string'),'pause')
                indt=indt+1;
                if indt>sz(3)
                    indt=1; % loop it
                end
                t1=t(indt);
                changeSlice();
                pause(1/FrameRate);
            end % while strcmpi(get(PlayPauseButton,'string'),'pause')
        end % if ~isempty(MovName)

% -=- Callback: When the user changes movie dialog slider -=-
    function inputFromMovieDialogSlider(~,~)
        limtmp = get(MovieDialogSlider,'userdata');
        [~,indt1]=min(abs(limtmp(1,1)-t));
        [~,indt2]=min(abs(limtmp(1,2)-t));
        MovieTimes = t([indt1 indt2]);
        set(MovieStartTimeValue,'String', num2str(MovieTimes(1),'%0.2f'))
        set(MovieEndTimeValue,'String', num2str(MovieTimes(2),'%0.2f'))
        UpdateDuration();
    end 

% -=- Callback: When the user changes file type -=-
    function inputFromFileTypeMenu(~,~)
        FileType = FileTypeList{get(FileTypeMenu,'value')};
        [~,MovName,~]=fileparts(MovName);
        switch FileType
            case 'Motion JPEG AVI'
                MovName = [MovName '.avi'];
            case 'MPEG-4'
                MovName = [MovName '.mp4'];
            case 'Archival'
                MovName = [MovName '.mj2'];
        end
        set(MovieDialog,'name',['Saving ' MovName ' ...']);
    end

% -=- Callback: When the user changes movie start/end time value boxes -=-
    function inputFromMovieStartTimeValue(~,~)
        limtmp=str2double(get(MovieStartTimeValue,'string')); 
        if isnan(limtmp)
            set(MovieStartTimeValue,'string',num2str(MovieTimes(1),'%0.2f'))
            return
        elseif limtmp<t(1)
            MovieTimes(1) = t(1);
        elseif limtmp>MovieTimes(2)
            MovieTimes(1) = MovieTimes(2);
        else
            [~,indt]=min(abs(limtmp-t));
            MovieTimes(1) = t(indt);
        end
        set(MovieStartTimeValue,'string',num2str(MovieTimes(1),'%0.2f'))
        allSlides = get(MovieDialogSlider, 'Children');
        loc = get(allSlides(1),'Position');
        loc(1) = (MovieTimes(1) - t(1)) ...
            * (1-loc(3)) / range(t);
        set(allSlides(1),'position',loc)
        UpdateDuration();
    end 
    function inputFromMovieEndTimeValue(~,~)
        limtmp=str2double(get(MovieEndTimeValue,'string')); 
        if isnan(limtmp)
            set(MovieEndTimeValue,'string',num2str(MovieTimes(2),'%0.2f'))
            return
        elseif limtmp>t(end)
            MovieTimes(2) = t(end);
        elseif limtmp<MovieTimes(1)
            MovieTimes(2) = MovieTimes(1);
        else
            [~,indt]=min(abs(limtmp-t));
            MovieTimes(2) = t(indt);
        end
        set(MovieEndTimeValue,'string',num2str(MovieTimes(2),'%0.2f'))
        allSlides = get(MovieDialogSlider, 'Children');
        loc = get(allSlides(2),'Position');
        loc(1) = (MovieTimes(2) - t(1)) ...
            * (1-loc(3)) / range(t);
        set(allSlides(2),'position',loc)
        UpdateDuration();
    end 

% -=- Callback: When the user changes new frame rate -=-
    function inputFromNewFrameRateValue(~,~)
        tmp = str2double(get(NewFrameRateValue,'string')); 
        if isnan(tmp) || tmp<0
            set(NewFrameRateValue,'string',num2str(FrameRate,'%0.3f')); 
            return
        end
        FrameRate = tmp;
        UpdateDuration();
    end

% -=- Callback: When the user changes downsample rate -=-
    function inputFromDownsampleValue(~,~)
        tmp = str2double(get(DownsampleValue,'string')); 
        if isnan(tmp) || tmp<=0.5 || isinf(tmp)
            set(DownsampleValue,'string',num2str(DownSample)); 
            return
        end
        DownSample = round(tmp);
        set(DownsampleValue,'string',num2str(DownSample)); 
        UpdateDuration();
    end

% -=- Update the duration of the movie to be saved -=-
    function UpdateDuration(~,~)
        set(DurationMessage,'string', ...
            ['duration: ' ...
                num2str(diff(MovieTimes)/DownSample/FrameRate, ...
                '%0.2f') ' s'])
    end

% -=- Callback: When the user clicks OK on the movie dialog -=-
    function inputFromMovieOKButton(~,~)
        close(MovieDialog)
        vid = VideoWriter(MovName,FileType);
        vid.FrameRate=FrameRate;
        open(vid);
        set(hff,'units','pixels');
        ind = find(t>=MovieTimes(1),1,'first') : ...
            DownSample : ...
            find(t<MovieTimes(2),1,'last');
        TimeList = t(ind);
        for jj=1:numel(TimeList)
            t1=TimeList(jj);
            changeSlice()
            pause(1/FrameRate);
            hf_pos = get(hff,'position');
            if get(AxesButton,'value') % axes only
                snap=getframe(ha);
            else
                snap=getframe(hff,[AxesLeftPos AxesBottomPos ...
                    hf_pos(3)-AxesLeftPos hf_pos(4)-AxesBottomPos]);
            end
            writeVideo(vid,snap.cdata);
        end % for jj=1:numel(TimeList)
        close(vid);
        msgbox(['Saved ' MovName '.'])
    end

% -=- Callback: When the user clicks Cancel on the movie dialog -=-
    function inputFromMovieCancelButton(~,~)
        close(MovieDialog)
    end

    end % function Play(~)

% -=- Pause playback -=-
    function Pause(~)
        if isempty(get(MovNameText,'string'))
            set(PlayPauseButton,'string','play')
        else
            set(PlayPauseButton,'string','save...')
        end
        set([SliceText MovNameText FrameRateText OverlayButton ...
            SliceUnitText AxesButton],'enable','on'); % re-enable controls
        if isVector
            set([QuiverScaleText CurlButton],'enable','on'); % re-enable controls
        end
        set(SliceSlider,'callback',@inputFromSliceSlider); % re-enable slider
        t1=t(indt);
        changeSlice() % just to clean up
        return
    end % function Pause(~)

% -=- Update controls and images when time changes -=-
    function changeSlice(~)
        [~,indt]=min(abs(t1-t));
        indt=min([max([indt 1]) sz(3)]); % don't allow times outside the range
        set(SliceText,'string',num2str(t(indt),'%.4g')) % update the text box
        set(SliceSlider,'Value',t(indt)); % update slider
        set(hff,'name',[num2str(t(indt),'%0.2f') ' ' ...
            get(SliceUnitText,'string') ' (' num2str(indt) ' of ' ...
            num2str(sz(3)) ')']);
        set(ht,'string',[num2str(t(indt),'%0.2f') ' ' ...
            get(SliceUnitText,'string')])
        updateImagesAndQuivers();
    end % function changeSlice(~)

% -=- Update images and quivers -=-
    function updateImagesAndQuivers(~)
        if isRGB 
            im1=cat(3, ...
                (squeeze(imComps{1}(:,:,indt))-clim_ch1(1))*ch1mult/diff(clim_ch1), ... % rgb channels, with normalized intensity
                (squeeze(imComps{2}(:,:,indt))-clim_ch2(1))*ch2mult/diff(clim_ch2), ...
                (squeeze(imComps{3}(:,:,indt))-clim_ch3(1))*ch3mult/diff(clim_ch3) );
            im1(im1<0)=0;
            if diff(clim_ch1)==0 || diff(clim_ch2)==0 || diff(clim_ch3)==0
                if diff(clim_ch1)==0
                    ch1=(squeeze(imComps{1}(:,:,indt)) ...
                        - clim_ch1(1))*ch1mult;
                else
                    ch1=(squeeze(imComps{1}(:,:,indt)) ...
                        - clim_ch1(1))*ch1mult/diff(clim_ch1);
                end
                if diff(clim_ch2)==0 
                    ch2=(squeeze(imComps{2}(:,:,indt)) ...
                        - clim_ch2(1))*ch2mult;
                else
                    ch2=(squeeze(imComps{2}(:,:,indt)) ...
                        - clim_ch2(1))*ch2mult/diff(clim_ch2);
                end
                if diff(clim_ch3)==0 
                    ch3=(squeeze(imComps{3}(:,:,indt)) ...
                        - clim_ch3(1))*ch3mult;
                else
                    ch3=(squeeze(imComps{3}(:,:,indt)) ...
                        - clim_ch3(1))*ch3mult/diff(clim_ch3);
                end
                im1=cat(3, ch1, ch2, ch3);
                im1(im1<0)=0;
            end
        else
            im1=squeeze(im(:,:,indt));
        end
        if exist('seg','var') %plots segmentations
            if ~isoutline %seg overlay
                for n=1:numel(seg)
                    colormapinput=[0,0,0];colormapinput(n)=1;
                    im1=labeloverlay(im1,seg{n}(:,:,indt),'ColorMap',colormapinput,'Trans', trans(n));
                end
            else %seg outline
                for n=1:numel(seg)
                    colormapinput=[0,0,0];colormapinput(n)=1;
                    im1=labeloverlay(im1,seg_outline{n}(:,:,indt),'ColorMap',colormapinput,'Trans', trans(n));
                end
            end
        end %if exist('seg','var')
        set(hi,'cdata',im1) % update image
        if exist('plt_coord','var')
            set(hp, 'xdata', plt_coord(:,1,indt), 'ydata', plt_coord(:,2,indt));
        end
        if isVector
            set(hq, ... % update quivers
                'udata', ...
                reshape(imComps{1}(qyind,qxind,indt),[],1)/maxlen*qs, ...
                'vdata', ...
                reshape(imComps{2}(qyind,qxind,indt),[],1)/maxlen*qs);
        end % if isVector

    end % function updateImagesAndQuivers(~)

end % function imagei
