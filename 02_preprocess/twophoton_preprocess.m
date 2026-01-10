function downsampled_2p = twophoton_preprocess(mdfstruct)
%TWOPHOTON_PREPROCESS Summary of this function goes here
% 1. Downsample images to 5 fps
% 2. Median filter to suppress noise

%   Detailed explanation goes here
disp('twophoton_preprocess running')
downsampled_2p.targetfrequency = 5;
downsampled_2p.medianframes = 5;

% Transfer metadata
if isfield(mdfstruct, 'img_param')
    downsampled_2p.startframe = mdfstruct.img_param.imgstartframe;
    downsampled_2p.endframe = mdfstruct.img_param.imgendframe;
end

disp('twophoton_preprocess ch2 downsample')
[preprocessed_ch2, outfps_ch2] = ...
    analyze_resample(mdfstruct.stackch2, mdfstruct.img_param.save_fps, downsampled_2p.targetfrequency);
disp('twophoton_preprocess ch1 downsample')
[preprocessed_ch1, ~] = ...
    analyze_resample(mdfstruct.stackch1, mdfstruct.img_param.save_fps, downsampled_2p.targetfrequency);
downsampled_2p.outfps = outfps_ch2;

disp('twophoton_preprocess ch1 medfilt')
downsampled_2p.ch1 = medfilt3(preprocessed_ch1,[1 1 downsampled_2p.medianframes]);
disp('twophoton_preprocess ch2 medfilt')
downsampled_2p.ch2 = medfilt3(preprocessed_ch2,[1 1 downsampled_2p.medianframes]);


% Calculate and attach metadata for downsampled data
if isfield(mdfstruct, 'img_param')
    downsampled_2p.pixel2um = mdfstruct.img_param.pixel2um;
    downsampled_2p.fps = outfps_ch2;
    % Recalculate time axis matching the new frame count
    n_frames = size(downsampled_2p.ch1, 3);
    downsampled_2p.t_axis = linspace(mdfstruct.img_param.imgstarttime, mdfstruct.img_param.imgendtime, n_frames);
end

end

