clear, clc
%make mcsx obj, get general (info), two photon scanning microscope imaging (info_ tpsm), and imaging mode specific (info_mode)  
[info, analog, mobj] = io_initmdf();

analogfig = analogtofig(info,analog);

%%
%Image loading parameter
param.start         = 300;       % [sec]
param.duration      = 1200;    % [sec]

param.frameend = param.totalframe + param.framestart;
param.groupz = 10;
info.refchannel = 1;

info.savefps = info.fps/param.groupz; 
param.framestart  = round(1+param.start*info.fps); % fps*sec
param.totalframe = round(param.duration*info.fps); % fps*sec

if param.frameend > info.fcount % if calculated frame end exceed end of frame, load from start to the end
    param.frameend = info.fcount;
end

%%
zstack = io_readframes(mobj,info.refchannel,[param.framestart, param.frameend]); % read frame from start to end (start+duration)

% Preprcocessing (Padding removal -> post pixel shift correction -> Trim -> Non Negative)
[tmp.xpadStart,tmp.xpadEnd] = pre_findpadding(zstack); % Find padded region caused by sinusoidal correction
zstack = zstack(:,tmp.xpadStart:tmp.xpadEnd,:);
zstack = zstack - min(zstack,[],'all');
%%
% motion correction by dft registration
dft_stack = pre_groupaverage(zstack, param.groupz); % denoise by group averaging 
dft_stack = medfilt3(dft_stack,[3,3,5]); % denoise by 3d median filter
drift_table = pre_estimatemotion(dft_stack); % using dft_registration, get drift table [error,diffphase,net_row_shift,net_col_shift]
%% using pixel shift information register the zstack
[interpdrifttable, corrected_z] = pre_pshiftcorr(zstack,drift_table);
gc_stack = pre_groupaverage(corrected_z,param.groupz);
gc_stack = medfilt3(gc_stack,[3,3,3]);
figure()
sliceViewer(gc_stack)
%% save
io_savetiff(gc_stack, info, info.refchannel)

%% process other channel

% clear memory
clearvars('corrected_z')
clearvars('zstack')
if info.refchannel == 1
    img_channel = 2;
else 
    img_channel = 1;
end
 
zstack = io_readframes(mobj,img_channel,[param.framestart, param.frameend]); % read frame from start to end (start+duration)
zstack = zstack(:,tmp.xpadStart:tmp.xpadEnd,:);
zstack = zstack - min(zstack,[],'all');

[interpdrifttable, corrected_z] = pre_pshiftcorr(zstack,drift_table);
gc_stack = pre_groupaverage(corrected_z,param.groupz);
gc_stack = medfilt3(gc_stack,[3,3,3]);

%%
figure()
sliceViewer(gc_stack)
%% save
io_savetiff(gc_stack, info, img_channel)
%%


