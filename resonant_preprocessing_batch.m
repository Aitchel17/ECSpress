clear, clc
% make mcsx obj, get general (info), two photon scanning microscope imaging (info_ tpsm), and imaging mode specific (info_mode)  
[info, analog, mobj] = io_initmdf();

% Image loading parameter

info.start         = 0;       % [sec], start frame
info.duration      = 2000;    % [sec], if duration is -1 or exceeding end of frame, read to end of frame
info.groupz = 10;
info.refchannel = 1;


info.savefps = info.fps/info.groupz; 
info.totalframe = round(info.duration*info.fps); % fps*sec
info.framestart  = round(1+info.start*info.fps); % fps*sec
info.frameend = info.totalframe + info.framestart;
if info.frameend > info.fcount % if calculated frame end exceed end of frame, load from start to the end
    info.frameend = info.fcount;
elseif info.frameend == -1 %
    info.frameend = info.fcount;
end
% demo processing
% Demo Image processing


% Create the save directory if it does not exist
save_folder = fullfile(info.mdfPath, info.mdfName(1:end-4));
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end



demo.fend = round((info.fcount - info.framestart)/20);
demo.stack = io_readframes(mobj,info.refchannel,[info.framestart, demo.fend]); % read frame from start to end (start+duration)
[tmp.xpadStart,tmp.xpadEnd] = pre_findpadding(demo.stack); % Find padded region caused by sinusoidal correction
info.xshift = pre_pshiftexplorer(demo.stack);
demo.stack = pre_pshiftcorrection(demo.stack,info.xshift);
demo.stack = pre_groupaverage(demo.stack(:,tmp.xpadStart:tmp.xpadEnd,:), info.groupz);
demo.stack = medfilt3(demo.stack,[3,3,5]);
[info.motionvertices, info.refslice] = roi_rectangle(demo.stack);
demo.drift_table = pre_estimatemotion(demo.stack,info.refslice,info.motionvertices);
[demo.ip_Drifttable, demo.correctedstack] = pre_applymotion(demo.stack,demo.drift_table);
%figure('Name','')
%sliceViewer(demo.correctedstack)
%figure()
%sliceViewer(demo.stack)

% Analog signal processing
% Ball processing
    tmp.numplot = 0;
    if isfield(analog,'raw_Ball') % if Ball channel detected
        % convert mdf metadata from string to double
        tmp.analogfreq = str2double(info.analogfreq(1:end-2)); % Hz
        tmp.analogresolution = str2double(info.analogresolution(1:end-4)); % bit
        tmp.Ballinputrange = str2double(info.Ballinputrange(2:end-1)); % V
        % make parameter array
        tmp.Ball_parameter = [tmp.analogfreq, tmp.analogresolution, tmp.Ballinputrange]; %[Hz, bit, V]
        % downsample input=[rawdata,parameterlist,downsamplingfactor,1d medianfilter],output=[xaxis, yaxis]
        analog.ds_ball = analog_preprocessing(analog.raw_Ball,tmp.Ball_parameter,10,3); %
        tmp.numplot = tmp.numplot +1;
    end
% EMG processing
    if isfield(analog,'raw_EMG')
        % mdf metadata to double
        tmp.analogfreq = str2double(info.analogfreq(1:end-2)); % Hz

        tmp.analogresolution = str2double(info.analogresolution(1:end-4)); % bit
        tmp.EMGinputrange = str2double(info.EMGinputrange(2:end-1)); % V
        tmp.EMG_parameter = [tmp.analogfreq, tmp.analogresolution,tmp.EMGinputrange];
        analog.ds_EMG = analog_preprocessing(analog.raw_EMG,tmp.EMG_parameter,10,3); % [xaxis, yaxis]
    end
    
% ECoG processing
    if isfield(analog,'raw_ECoG')
        % mdf metadata to double
        tmp.analogfreq = str2double(info.analogfreq(1:end-2)); % Hz
        % mutitaper spectrum
        analog.ecog_spectrum = analog_ecogspectrum(tmp.analogfreq,analog.raw_ECoG); % calculate ecogspectrum  .taxis (sec), .faxis(Hz), .spectral power(DB)
    end

% Plot analog channel result
[analog_fig, analog_axes] = analog_plot(analog);

% Real image processing
zstack = io_readframes(mobj,info.refchannel,[info.framestart, info.frameend]); % read frame from start to end (start+duration)
% Preprcocessing (Padding removal -> post pixel shift correction -> Trim -> Non Negative)
disp('Padding removal')
zstack = zstack(:,tmp.xpadStart:tmp.xpadEnd,:);
% xshiftcorrection
zstack = pre_pshiftcorrection(zstack,info.xshift);
% non negative
disp('min subtraction for non negative array')
zstack = zstack - min(zstack,[],'all');
%
% generate drift table by dft registration
dft_stack = pre_groupaverage(zstack, info.groupz); % denoise by group averaging
disp('3D median filtering')
dft_stack = medfilt3(dft_stack,[3,3,5]); % denoise by 3d median filter
drift_table = pre_estimatemotion(dft_stack,info.refslice,info.motionvertices); % using dft_registration, get drift table [error,diffphase,net_row_shift,net_col_shift]

% apply drift table, correct motion
[applied_drifttable, zstack] = pre_applymotion(zstack,drift_table);
gc_stack = pre_groupaverage(zstack,info.groupz);

figure('Name','Corrected reference channel')
sliceViewer(gc_stack)
% save
io_savetiff(gc_stack, save_folder,info, info.refchannel)

% process other channel
% clear memory
clearvars('corrected_z')
clearvars('zstack')
if info.refchannel == 1
    img_channel = 2;
else 
    img_channel = 1;
end
 
zstack = io_readframes(mobj,img_channel,[info.framestart, info.frameend]); % read frame from start to end (start+duration)
% Padding removal
disp('Padding removal')
zstack = zstack(:,tmp.xpadStart:tmp.xpadEnd,:);
% xshiftcorrection
zstack = pre_pshiftcorrection(zstack,info.xshift);
% non negative
disp('min subtraction for non negative array')
zstack = zstack - min(zstack,[],'all');
% Motion correction using 
[interpdrifttable, zstack] = pre_applymotion(zstack,drift_table);
gc_stack = pre_groupaverage(zstack,info.groupz);

%
figure('Name','Corrected following channel')
sliceViewer(gc_stack)
% save
io_savetiff(gc_stack, save_folder, info, img_channel);

io_saveinfo(info,save_folder);
io_saveanalog(analog,save_folder,info);
%%
max(demo.stack,[],"all")
%%
in16array = uint16(doublearray);


