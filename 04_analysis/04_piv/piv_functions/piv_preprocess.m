function imgstack_out = piv_preprocess(imgstack, opt)
%PIV_PREPROCESS  Run PIVlab image pre-processing over a whole stack.
%   Frame-by-frame wrapper around preproc.PIVlab_preproc, mirroring PIVlab's
%   Image settings > Image pre-processing tab. Each frame is scaled to 0~1 first.
%   With roirect set, only that rectangle is filtered; the rest passes through.
%
% IN   imgstack       H x W x T numeric stack
%      roirect        [x y width height] filter region, default [] (whole frame)
%      prctile_low    imadjust lower limit on the 0~1 frame, default 0.05
%      prctile_high   imadjust upper limit, default 0.95
%      clahe          CLAHE on/off, default true
%      clahe_size     CLAHE tile size (px), default 64
%      highpass       high-pass on/off, default false
%      highpass_size  high-pass kernel size (px), default 15
%      intens_cap     cap intensity above median+2*std at median+2*std, default false
%      wiener         Wiener denoise on/off; also applies the low-pass, default true
%      wiener_size    Wiener window (px); low-pass sigma is half of it, default 3
% OUT  imgstack_out   same size and class as imgstack, filtered

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


% 0. Setup
n_frames   = size(imgstack, 3);
imgstack_out = imgstack;          % pre-allocate output (same class)

% 1. Filter one frame at a time
for k = 1:n_frames

    % 1.1 Scale the frame to 0~1
    frame = imgstack(:,:,k);
    norm_frame = (frame-min(frame(:)))/max(frame(:));
    % 1.2 Run the PIVlab filter chain
    frame_proc = preproc.PIVlab_preproc( ...
        norm_frame, opt.roirect, ...
        opt.clahe, opt.clahe_size, ...
        opt.highpass, opt.highpass_size, ...
        opt.intens_cap, ...
        opt.wiener, opt.wiener_size, ...
        opt.prctile_low, opt.prctile_high);

    imgstack_out(:,:,k) = frame_proc;
    % 1.3 Progress line every tenth frame
    if mod(k,10)==1
        fprintf('  frame %d / %d  (%.0f%%)\n', k, n_frames, 100*k/n_frames);
    end
end

end
