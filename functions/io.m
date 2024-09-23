%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prerequisite: MCSX 
%                       >> https://www.sutter.com/MICROSCOPES/mcs.html

% FUNCTION NAME:    io_mdf
%
% DESCRIPTION:      Sutter .mdf Data input - output function
% INPUT:            .mdf
%
% NOTES:            1. If wavelength and objective length 
%
% WRITTEN BY:       C. Hyunseok Lee 2024-09-14
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[path.mdfName, path.folderPath] = uigetfile('*.mdf'); % select file by UI
path.mdfPath = [path.folderPath, path.mdfName];
mobj = actxserver('MCSX.Data'); % Create Component Object Model (COM)
mobj.invoke('OpenMCSFile',path.mdfPath); % Using COM open .mdf file

%%
% General Info
info.User = mobj.ReadParameter('Created by');
info.Date = mobj.ReadParameter('Created on');

%% two photon scanning microscope info
tpsm_info.scanmode      = mobj.ReadParameter('Scan Mode');
tpsm_info.pclock        = mobj.ReadParameter('Pixel Clock');
tpsm_info.yoffset       = mobj.ReadParameter('Y Frame Offset');

% Objective info
tpsm_info.objname       = mobj.ReadParameter('Objective');
tpsm_info.objpix        = mobj.ReadParameter('Microns per Pixel');
tpsm_info.objx          = mobj.ReadParameter('X Position');
tpsm_info.objy          = mobj.ReadParameter('Y Position');
tpsm_info.objz          = mobj.ReadParameter('Z Position');

% Scan info
tpsm_info.zoom          = mobj.ReadParameter('Magnification');
tpsm_info.tpex          = strcat(mobj.ReadParameter('Laser Wavelength (nm)'),' nm');
tpsm_info.fbit          = mobj.ReadParameter('Frame Bit Depth');
tpsm_info.fdur          = mobj.ReadParameter('Frame Duration (s)');
tpsm_info.fcount        = str2double(mobj.ReadParameter('Frame Count'));
tpsm_info.fint          = mobj.ReadParameter('Frame Interval (ms)');
tpsm_info.fps           = 1/str2double(tpsm_info.fdur(1:end-1)); % Hz
tpsm_info.fh            = mobj.ReadParameter('Frame Height');
tpsm_info.fw            = str2double(mobj.ReadParameter('Frame Width'));
tpsm_info.lpower        = mobj.ReadParameter('Laser intensity');

%% Behavicor camera info
bcam_info.enabled = mobj.ReadParameter('Behavior Video Enabled');


%% to generate 1 Hz summary video
 
for frame = 10:tpsm_info.fps:tpsm_info.fcount
    double(invoke(mObj, 'ReadFrame', 2, 4000))'


if strcmp(info.scanmode,'XY Movie')
    print('XY movie loaded')
end

tmp.fps = 1/str2double(tpsm_info.fdur(1:end-1)); % Hz


%%
img = double(invoke(mObj, 'ReadFrame', 2, 4000))'; % double(invoke(mObj, 'ReadFrame', ch, frameIdx))';

figure(1)
imshow(img,[0 500])
%%
framesNum = length(1:mdfInfo.totalFrames());
%%
h = waitbar(0,'Buffering Frames...');
%%
bufImg = double(invoke(mObj, 'ReadFrame', 4, 1))';
%%
imSize = size(bufImg, 1);
imgs = zeros(imSize, imSize, framesNum);
for i = startIdx:50:endIdx
	% disp(i)
	curIdx = (i-startIdx)/framesNum;
	waitbar(curIdx, h, ...
		sprintf(['Loading:' '%d /' num2str(endIdx-startIdx+1)], i-startIdx));
	bufImg = mcsxReadFrame(mObj, ch, i);
	imgs(:,:,1+i-startIdx)=bufImg;
end
delete(h)

