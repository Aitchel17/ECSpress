function downsampled_2p = twophoton_preprocess(mdfstruct)
%TWOPHOTON_PREPROCESS Summary of this function goes here
    % 1. Downsample images to 5 fps
    % 2. Median filter to suppress noise
    
%   Detailed explanation goes here
    disp('twophoton_preprocess running')
    downsampled_2p.frequency = 5;
    downsampled_2p.medianframes = 5;

    disp('twophoton_preprocess ch2 downsample')
    [preprocessed_ch2, outfps_ch2] = ...
        analyze_resample(mdfstruct.stackch2, mdfstruct.img_param.save_fps, downsampled_2p.frequency);
    disp('twophoton_preprocess ch1 downsample')
    [preprocessed_ch1, outfps_ch1] = ...
        analyze_resample(mdfstruct.stackch1, mdfstruct.img_param.save_fps, downsampled_2p.frequency);

    disp('twophoton_preprocess ch1 medfilt')
    downsampled_2p.ch1 = medfilt3(preprocessed_ch1,[1 1 downsampled_2p.medianframes]);
    disp('twophoton_preprocess ch2 medfilt')
    downsampled_2p.ch2 = medfilt3(preprocessed_ch2,[1 1 downsampled_2p.medianframes]);

end

