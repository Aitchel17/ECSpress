function downsampled_2p = twophoton_preprocess(session_instance)
%TWOPHOTON_PREPROCESS Summary of this function goes here
% 1. Downsample images to 5 fps
% 2. Median filter to suppress noise

%   Detailed explanation goes here
disp('twophoton_preprocess running')
downsampled_2p.targetfrequency = 3;
downsampled_2p.medianframes = 11;

% Transfer metadata
if isfield(session_instance, 'img_param')
    downsampled_2p.startframe = session_instance.img_param.imgstartframe;
    downsampled_2p.endframe = session_instance.img_param.imgendframe;
end

disp('twophoton_preprocess ch2 downsample')
[preprocessed_ch2, outfps_ch2] = ...
analyze_resample(session_instance.stackch2, session_instance.img_param.save_fps, downsampled_2p.targetfrequency);
disp('twophoton_preprocess ch1 downsample')
[preprocessed_ch1, ~] = ...
analyze_resample(session_instance.stackch1, session_instance.img_param.save_fps, downsampled_2p.targetfrequency);
downsampled_2p.outfps = outfps_ch2;

disp('twophoton_preprocess ch1 medfilt')
downsampled_2p.ch1 = medfilt3(preprocessed_ch1,[1 1 downsampled_2p.medianframes]);
disp('twophoton_preprocess ch2 medfilt')
downsampled_2p.ch2 = medfilt3(preprocessed_ch2,[1 1 downsampled_2p.medianframes]);


% Calculate and attach metadata for downsampled data
if isprop(session_instance, 'img_param')
    downsampled_2p.pixel2um = session_instance.img_param.pixel2um;
    downsampled_2p.fps = outfps_ch2;
    % Recalculate time axis matching the new frame count
    n_frames = size(downsampled_2p.ch1, 3);
    downsampled_2p.t_axis = linspace(session_instance.img_param.imgstarttime, session_instance.img_param.imgendtime, n_frames);
end

end

