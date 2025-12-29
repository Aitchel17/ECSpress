function [resampled_stack, resampled_fps] = analyze_resample(stack, inFPS, outFPS)
%%
    % Input data: 3d image stack with (x,y,frames)
    % Input parameters: recording fps 'inFPS' and output fps 'outFPS'
    % Output: resampled image stack
    % Output parameter: actual resampled_fps

    % Steps
    % 1. Collapse data into (pixels, frames) for GPU processing
    % 2. Determine upsample and downsample multiple at 1e-6 accyuracy
    % 3. Resample
    % 4. Reconstruct

    sz = size(stack);
    t_pixels = reshape(stack,[],sz(3)).';
    t_pixels = gpuArray(t_pixels);
    [p,q] = rat(outFPS/inFPS, 1e-6);
    resampled = resample(double(t_pixels),p,q);
    out_tdim = size(resampled,1); 
    resampled_stack = reshape(resampled.',[sz(1), sz(2) out_tdim]);
    resampled_stack = gather(resampled_stack);
    resampled_fps = inFPS* (p/q);
end