function flow = opticalflow_wlet(imgstack, opt)
%OPTICALFLOW_WLET  Run PIVlab wavelet optical flow (wOFV) on an image stack.
%   Applies wOFV to each consecutive frame pair (k, k+1) in imgstack.
%   imgstack should already be preprocessed (e.g. via opticalflow_preprocess).
%
%   USAGE
%     result = opticalflow_wlet(imgstack)
%     result = opticalflow_wlet(imgstack, pyramid_levels=3, smoothness=50)
%
%   INPUTS
%     imgstack         : double H x W x N array, values in [0,1]
%     pyramid_levels   : int    - number of pyramid levels (default: 3)
%                        higher → captures larger displacements, slower
%     smoothness       : scalar - regularization strength, GUI scale [0,100]
%                        (default: 50 → eta = 10^0 = 1)
%                        lower = less smooth / more detail
%                        higher = smoother / less noise
%     roirect          : [x y w h] - ROI to analyze (default: full frame)
%
%   OUTPUT
%     flow : double H x W x N_pairs x 2
%            flow(:,:,k,1) = u (horizontal velocity) for pair k
%            flow(:,:,k,2) = v (vertical velocity)   for pair k

arguments
    imgstack         {mustBeNumeric, mustBeNonempty}
    opt.pyramid_levels double  = 3
    opt.smoothness     double  = 50    % GUI scale [0,100]
    opt.roirect        double  = []

end

n_frames  = size(imgstack, 3);
n_pairs   = n_frames - 1;

if n_pairs < 1
    error('opticalflow_wlet: imgstack must have at least 2 frames.');
end

%% Convert GUI smoothness [0,100] → eta (same formula PIVlab uses)
%  etaUnscaled = 50 → eta = 10^(50*0.1 - 5) = 10^0 = 1
eta = 10^(opt.smoothness * 0.1 - 5);

%% Build roirect (full frame if empty)
[H, W, ~] = size(imgstack);
roirect = opt.roirect;
if isempty(roirect)
    roirect = [1, 1, W-1, H-1];
end

%% Pre-compute filter matrices from first frame (all frames same size)
img0     = imgstack(:,:,1);
roi_img  = img0(roirect(2):roirect(2)+roirect(4)-1, ...
                roirect(1):roirect(1)+roirect(3)-1);
PatchSize = 2^floor(log2(min(size(roi_img))));
Fmats     = wOFV.getFmatPyramid(PatchSize, opt.pyramid_levels);

%% No mask (all pixels valid)
mask = false(H, W);   % PIVlab convention: true = exclude

%% Median filter settings (off by default, matching PIVlab default 'Off')
MedFiltFlag = false;
MedFiltSize = [3, 3];

vartheta = ones(size(roi_img));   % uniform spatial weighting

%% Main loop - RunMain_Parallel_DatasetProc uses parfor over spatial patches
%  Open a parpool before calling this function to activate parallel workers.
%% Pre-allocate output: H x W x N_pairs x 2  (dim 4: 1=u, 2=v)
flow = zeros(H-1, W-1, n_pairs, 2);

for k = 1:n_pairs

    fprintf('  wOFV pair %d/%d\n', k, n_pairs);

    [~, ~, u, v, ~] = wOFV.RunMain_Parallel_DatasetProc( ...
        imgstack(:,:,k), imgstack(:,:,k+1), ...
        mask, roirect, ...
        eta, vartheta, ...
        MedFiltFlag, MedFiltSize, ...
        opt.pyramid_levels, ...
        Fmats, PatchSize);

    flow(:,:,k,1) = u;
    flow(:,:,k,2) = v;
end

fprintf('opticalflow_wlet: done. %d pairs processed.\n', n_pairs);

end