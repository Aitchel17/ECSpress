function imgstack_out = opticalflow_preprocess(imgstack, opt)
% This function just mimic PIVlab Image settings tab > Image pre-processing tab
% The main function is in PIVlab folder > preproc folder > PIVlab_preproc.m
% Input:  image stack (double)
% Output: image stack (double)

% 1. Normalize image from 0~1
% 2. PIV lab preprocessing to do 
    % 2.1 CLAHE (require window size)
    % 2.2 high-pass (require maginitude)
    % 2.3 intensity capping (true false), intensity above median+2std capped by median+2std
    % 2.4 Wiener (require size)
    % 2.5 Low pass filter automatically applied with Wiener filter using wiener filter size and half of wiener filter size as sigma


arguments
    imgstack   {mustBeNumeric, mustBeNonempty}
    opt.roirect        double  = []
    opt.prctile_low        double  = 0.05
    opt.prctile_high        double  = 0.95
    opt.clahe          logical = true
    opt.clahe_size     double  = 64
    opt.highpass       logical = false
    opt.highpass_size  double  = 15
    opt.intens_cap     logical = false
    opt.wiener         logical = true
    opt.wiener_size    double  = 3
end


n_frames   = size(imgstack, 3);
imgstack_out = imgstack;          % pre-allocate output (same class)

%% Main loop
for k = 1:n_frames

    frame = imgstack(:,:,k);
    norm_frame = (frame-min(frame(:)))/max(frame(:));
    frame_proc = preproc.PIVlab_preproc( ...
        norm_frame, opt.roirect, ...
        opt.clahe, opt.clahe_size, ...
        opt.highpass, opt.highpass_size, ...
        opt.intens_cap, ...
        opt.wiener, opt.wiener_size, ...
        opt.prctile_low, opt.prctile_high);

    imgstack_out(:,:,k) = frame_proc;
    if mod(k,10)==1
        fprintf('  frame %d / %d  (%.0f%%)\n', k, n_frames, 100*k/n_frames);
    end
end

end
