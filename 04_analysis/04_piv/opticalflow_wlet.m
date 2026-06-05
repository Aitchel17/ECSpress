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
%% Median filter settings (off by default, matching PIVlab default 'Off')
MedFiltFlag = false;
MedFiltSize = [3, 3];
%% Pre-allocate output: H x W x N_pairs x 2  (dim 4: 1=u, 2=v)
flow = NaN(H, W, n_pairs, 2);
% Get roirect coordinates for mapping back
rx = roirect(1);
ry = roirect(2);
rw = roirect(3);
rh = roirect(4);

%% single patch size to prevent grid/aperture artifacts
M = 2^ceil(log2(max(rh, rw)));
PatchSize = M;
Fmats = wOFV.getFmatPyramid(PatchSize, opt.pyramid_levels);

% Determine the ideal extraction box centered around the ROI to grab actual context
cx = round(rx + rw/2);
cy = round(ry + rh/2);

nx = round(cx - M/2);
ny = round(cy - M/2);

% Clamp to image boundaries
nx = max(1, nx);
ny = max(1, ny);
if nx + M - 1 > W
    nx = W - M + 1;
    nx = max(1, nx);
end
if ny + M - 1 > H
    ny = H - M + 1;
    ny = max(1, ny);
end

ex = min(W, nx + M - 1);
ey = min(H, ny + M - 1);

actual_w = ex - nx + 1;
actual_h = ey - ny + 1;

% Offset of the user's ROI inside the M x M block
offset_x = rx - nx + 1;
offset_y = ry - ny + 1;

%% Pre-allocate padded inputs
mask_pad = true(M, M); % true means exclude in PIVlab
mask_pad(1:actual_h, 1:actual_w) = false; % Include ALL actual image context
roirect_pad = [1, 1, M, M];
vartheta_pad = ones(M, M);

for k = 1:n_pairs

    fprintf('  wOFV pair %d/%d\n', k, n_pairs);

    % Extract actual context around ROI to prevent boundary regularizer artifacts
    I0_actual = imgstack(ny:ey, nx:ex, k);
    I1_actual = imgstack(ny:ey, nx:ex, k+1);
    
    % Pad only if M > image dimensions ('symmetric' prevents sharp black edges)
    I0_pad = padarray(I0_actual, [M-actual_h, M-actual_w], 'symmetric', 'post');
    I1_pad = padarray(I1_actual, [M-actual_h, M-actual_w], 'symmetric', 'post');

    [~, ~, u_pad, v_pad, ~] = wOFV.RunMain_DatasetProc( ...
        I0_pad, I1_pad, ...
        mask_pad, roirect_pad, ...
        eta, vartheta_pad, ...
        MedFiltFlag, MedFiltSize, ...
        opt.pyramid_levels, ...
        Fmats, PatchSize);

    % Extract the specific ROI region back from the padded result
    u = u_pad(offset_y : offset_y + rh - 1, offset_x : offset_x + rw - 1);
    v = v_pad(offset_y : offset_y + rh - 1, offset_x : offset_x + rw - 1);

    % Map cropped u, v back into full-frame coordinates
    flow(ry : ry+rh-1, rx : rx+rw-1, k, 1) = u;
    flow(ry : ry+rh-1, rx : rx+rw-1, k, 2) = v;
end

fprintf('opticalflow_wlet: done. %d pairs processed.\n', n_pairs);

end