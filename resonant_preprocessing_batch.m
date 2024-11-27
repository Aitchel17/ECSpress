clear, clc
% make mcsx obj, get general (info), two photon scanning microscope imaging (info_ tpsm), and imaging mode specific (info_mode)  
[info, analog, mobj] = io_initmdf();

%% Image loading parameter
param.start         = 0;       % [sec], start frame
param.duration      = 2000;    % [sec], if duration exceed end of frame, read entire frame
param.groupz = 10;
info.refchannel = 1;

%% Analog signal processing
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
analog_plot(analog);


%%
param.totalframe = round(param.duration*info.fps); % fps*sec
param.framestart  = round(1+param.start*info.fps); % fps*sec
param.frameend = param.totalframe + param.framestart;
if param.frameend > info.fcount % if calculated frame end exceed end of frame, load from start to the end
    param.frameend = info.fcount;
elseif param.frameend == -1 %
    param.frameend = info.fcount;
end

info.savefps = info.fps/param.groupz; 

%%
demostack = io_readframes(mobj,info.refchannel,[param.framestart, param.frameend]); % read frame from start to end (start+duration)
[tmp.xpadStart,tmp.xpadEnd] = pre_findpadding(demostack); % Find padded region caused by sinusoidal correction

zstack = io_readframes(mobj,info.refchannel,[param.framestart, param.frameend]); % read frame from start to end (start+duration)
% Preprcocessing (Padding removal -> post pixel shift correction -> Trim -> Non Negative)
zstack = zstack(:,tmp.xpadStart:tmp.xpadEnd,:);
zstack = zstack - min(zstack,[],'all');





% motion correction by dft registration
dft_stack = pre_groupaverage(zstack, param.groupz); % denoise by group averaging
disp('3D median filtering')
dft_stack = medfilt3(dft_stack,[3,3,5]); % denoise by 3d median filter
drift_table = pre_estimatemotion(dft_stack); % using dft_registration, get drift table [error,diffphase,net_row_shift,net_col_shift]
% using pixel shift information register the zstack
[~, corrected_z] = pre_applymotion(zstack,drift_table);
gc_stack = pre_groupaverage(corrected_z,param.groupz);
disp('3D median filtering')
gc_stack = medfilt3(gc_stack,[3,3,3]);
figure()
sliceViewer(gc_stack)
% save
io_savetiff(gc_stack, info, info.refchannel)

% process other channel

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

[interpdrifttable, corrected_z] = pre_applymotion(zstack,drift_table);
gc_stack = pre_groupaverage(corrected_z,param.groupz);
gc_stack = medfilt3(gc_stack,[3,3,3]);

%
figure()
sliceViewer(gc_stack)
% save
io_savetiff(gc_stack, info, img_channel)

%%
infoFields = fieldnames(info);

infoValues = struct2cell(info);
infoTable = table(infoFields, infoValues, 'VariableNames', {'Field', 'Value'});


