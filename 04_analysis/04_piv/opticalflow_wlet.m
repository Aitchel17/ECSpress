function flow = opticalflow_wlet(imgstack, opt)
%OPTICALFLOW_WLET  Run PIVlab wavelet optical flow (wOFV) on an image stack.
%   Applies wOFV to each consecutive frame pair (k, k+1). imgstack must already be
%   preprocessed (e.g. by piv_preprocess). Flow is px/frame, unscaled.
%
% IN   imgstack        H x W x N double, values in [0,1]
%      pyramid_levels  pyramid depth; higher = larger displacements, slower, default 3
%      smoothness      PIVlab GUI scale [0,100]; lower = more detail, default 50
%      roirect         [x y w h] region to analyse, default full frame
% OUT  flow            H x W x (N-1) x 2 double; (:,:,k,1) = u horizontal,
%                      (:,:,k,2) = v vertical, for pair k

arguments
    imgstack         {mustBeNumeric, mustBeNonempty}
    opt.pyramid_levels double  = 3
    opt.smoothness     double  = 50    % GUI scale [0,100]
    opt.roirect        double  = []
end

% 0. Setup
% 0.1 Number of frame pairs
n_frames  = size(imgstack, 3);
n_pairs   = n_frames - 1;
if n_pairs < 1
    error('opticalflow_wlet: imgstack must have at least 2 frames.');
end
% 0.2 GUI smoothness to regularizer eta, PIVlab formula: 50 -> 10^0 = 1
eta = 10^(opt.smoothness * 0.1 - 5);
% 0.3 ROI rect, full frame if empty
[H, W, ~] = size(imgstack);
roirect = opt.roirect;
if isempty(roirect)
    roirect = [1, 1, W-1, H-1];
end
% 0.4 Median filter settings, off as in PIVlab
MedFiltFlag = false;
MedFiltSize = [3, 3];
% 0.5 Output, dim 4: 1 = u, 2 = v
flow = NaN(H, W, n_pairs, 2);
% 0.6 ROI corner and size
rx = roirect(1);
ry = roirect(2);
rw = roirect(3);
rh = roirect(4);

% 1. Build the M x M analysis block around the ROI
% 1.1 One power-of-two patch size covering the whole ROI
M = 2^ceil(log2(max(rh, rw)));
PatchSize = M;
Fmats = wOFV.getFmatPyramid(PatchSize, opt.pyramid_levels);
% 1.2 Block centred on the ROI centre
cx = round(rx + rw/2);
cy = round(ry + rh/2);
nx = round(cx - M/2);
ny = round(cy - M/2);
% 1.3 Clamp the block inside the image
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
% 1.4 Image extent actually available inside the block
ex = min(W, nx + M - 1);
ey = min(H, ny + M - 1);
actual_w = ex - nx + 1;
actual_h = ey - ny + 1;
% 1.5 Offset of the ROI inside the block
offset_x = rx - nx + 1;
offset_y = ry - ny + 1;

% 2. Padded inputs shared by every pair
mask_pad = true(M, M);                      % true = exclude, PIVlab convention
mask_pad(1:actual_h, 1:actual_w) = false;
roirect_pad = [1, 1, M, M];
vartheta_pad = ones(M, M);

% 3. Solve the flow field for each consecutive frame pair
for k = 1:n_pairs

    fprintf('  wOFV pair %d/%d\n', k, n_pairs);

    % 3.1 Image context around the ROI
    I0_actual = imgstack(ny:ey, nx:ex, k);
    I1_actual = imgstack(ny:ey, nx:ex, k+1);

    % 3.2 Fill the block out to M x M
    I0_pad = padarray(I0_actual, [M-actual_h, M-actual_w], 'symmetric', 'post');
    I1_pad = padarray(I1_actual, [M-actual_h, M-actual_w], 'symmetric', 'post');

    % 3.3 Run wOFV on the block
    [~, ~, u_pad, v_pad, ~] = wOFV.RunMain_DatasetProc( ...
        I0_pad, I1_pad, ...
        mask_pad, roirect_pad, ...
        eta, vartheta_pad, ...
        MedFiltFlag, MedFiltSize, ...
        opt.pyramid_levels, ...
        Fmats, PatchSize);

    % 3.4 Crop the ROI back out of the block
    u = u_pad(offset_y : offset_y + rh - 1, offset_x : offset_x + rw - 1);
    v = v_pad(offset_y : offset_y + rh - 1, offset_x : offset_x + rw - 1);

    % 3.5 Place the ROI in full-frame coordinates
    flow(ry : ry+rh-1, rx : rx+rw-1, k, 1) = u;
    flow(ry : ry+rh-1, rx : rx+rw-1, k, 2) = v;
end

fprintf('opticalflow_wlet: done. %d pairs processed.\n', n_pairs);

end
