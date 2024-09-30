clear, clc
%%

%make mcsx obj, get general (info), two photon scanning microscope imaging (info_ tpsm), and imaging mode specific (info_mode)  
[info, info_tpsm, info_mode, mobj] = io_initmdf();
%%
param.start         = 0;       % [sec]
param.duration      = 600;    % [sec]
param.pshift        = 1;        % [pixel]
param.medfilt       = [3 3 1]; % 3d median filter voxel size [Y X Z]
%%
param.framestart  = round(1+param.start*info_tpsm.fps); % fps*sec
param.frameduration = round(param.duration*info_tpsm.fps); % fps*sec

% Preprcocessing (Padding removal -> post pixel shift correction -> Trim)
% Find padded region caused by sinusoidal correction
[tmp.xpadStart,tmp.xpadEnd,tmp.meanyz] = pre_findpadding(mobj,2,info_tpsm,param.framestart);
% Load raw stack without padding
tmp.rawstack = io_readframes(mobj,2,[tmp.xpadStart,tmp.xpadEnd],[1,info_tpsm.fh],[param.framestart,param.framestart+param.frameduration]);
% pixelshift correction
[tmp.pcorstack, data.avepcor] = pre_pshiftcorr(tmp.rawstack,param.pshift);
%
figure('Name','After pixel correction');
imshow(data.avepcor(:,10:end-param.pshift), [0,2000]);

%
% remove aberrant region at the side and caused by post pixel correction
data.prestack = tmp.pcorstack(:,10:end-param.pshift,:);

clear('tmp')
figure('Name','raw imagestack');
sliceViewer(data.prestack)
% Median filter
tmp.medstack = medfilt3(data.prestack,[3 3 1]);
figure('Name','medfilt imagestack')
sliceViewer(tmp.medstack)

%%
param.gpave = 19;
tmp.gpave_stack = squeeze( ...
    mean( ...
    reshape(tmp.medstack,size(data.prestack,1),size(tmp.medstack,2),param.gpave,[]),3) ... % reshape by 
       ... % calculate mean by axis
    );
figure()
figure('Name','gpave imagestack')
sliceViewer(tmp.gpave_stack)
%%
info.savefolder = uigetdir(); % select file by UI
%%
info.savePath = [info.savefolder,'\',info.mdfName(1:end-4)];
mkdir(info.savePath)
%%
info.savename = ['medfilt',erase(int2str(param.medfilt),' '),'gpave',int2str(param.gpave),'.tif'];
%%
imwrite(tmp.gpave_stack,[info.savePath,'\',info.savename])
%%
tifobj = Tiff([info.savePath,'\',info.savename],"w");


%%
write(tifobj,tmp.gpave_stack(:,:,1))