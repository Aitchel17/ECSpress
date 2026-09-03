function [tirs_result, sinogram_stack] = analyze_radon(hold_stack,the_angles,rtd_threshold)
%ANALYZE_RADON  Vessel width per angle per frame, by thresholding in Radon space.
%   Every frame is projected at each of the_angles. Each projection is scaled to
%   0~1 on its own, so a frame getting brighter or dimmer does not move the
%   boundaries -- only the shape of the profile does. The two boundaries are then
%   the outermost bins on either side of the peak that fall below rtd_threshold,
%   which makes the width a full-width-at-that-fraction: 0.5 is FWHM.
%
%   The frames go onto the GPU on entry and both outputs come back on the CPU,
%   so the caller never sees a gpuArray.
%
% IN   hold_stack     H x W x T numeric   the frames to measure
%      the_angles     1 x A double deg    projection angles, as radon takes them
%      rtd_threshold  1 x 1 double        fraction of the normalised peak, 0~1
% OUT  tirs_result    1 x 1 struct
%        radon_size                  1 x 3 double  size of the projection stack
%        idx_maxloc                  A x T double  peak bin, per angle per frame
%        idx_uploc                   A x T double  boundary above the peak
%        idx_downloc                 A x T double  boundary below it
%        diameter                    A x T double  downloc - uploc, in BINS
%        center_loc                  A x T double  downloc + uploc -- the SUM,
%                                                  so twice the midpoint
%        median_diameter             A x 1 double  over frames
%        normalized_diameterchange   A x T double  diameter / median - 1
%        var_normdiameter            A x 1 double  over frames
%        var_diameter                A x 1 double  over frames
%      sinogram_stack  P x A x T double  the projections themselves, scaled
%                                        the same way the boundaries were found
%                                        on. Only built when asked for -- it is
%                                        the one array here big enough to matter,
%                                        and re-thresholding is what it is for
%
%   UNIT  bins of the Radon projection, not microns. One bin is one pixel only
%         along the projection direction; converting needs objpix and the angle
%
% 0. Put stack to memory of gpu
hold_stack = gpuArray(hold_stack);

% 3. Normalize stack
radon_stack = zeros([size(radon(hold_stack(:,:,1),the_angles)),size(hold_stack,3)]); % 3.1
disp('Radon transform start')
for idx = 1:size(hold_stack,3)
    radon_stack(:,:,idx)=radon(hold_stack(:,:,idx),the_angles); %3.2
end
disp('Radon transform end')

% 4. Normalize radon stack
radon_stack = radon_stack-min(radon_stack,[],1); % 4.1
radon_stack = radon_stack./max(radon_stack,[],1); % 4.2

% 5. Make maximum location matrix
disp('Radon thresholding start')
[~, maxlocarray] = max(radon_stack,[],1); % 5.1
maxlocarray = squeeze(maxlocarray);
sz = size(radon_stack);

% We need temporary mask for calculation, but we won't store full tirs yet.
mask = false(sz); % 5.2
[row_idx, ~, ~] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3)); % 5.3
maxlocarray3d = repmat(reshape(maxlocarray, [1, sz(2), sz(3)]), [sz(1), 1, 1]); % 5.4

% 6. upper processing
mask(row_idx <= maxlocarray3d) = 1; % 6.1
clearvars maxlocarray3d
radon_thr = radon_stack<rtd_threshold; % 6.2
upperboundary_idx = row_idx .* radon_thr.* mask; % 6.3
upperboundary_idx = max(upperboundary_idx, [], 1); % 6.4
upperboundary_idx = squeeze(upperboundary_idx);

% 7. bottom processing
mask = ~mask; % 7.1
bottomboundary_idx = row_idx.*radon_thr.*mask; % 7.2
bottomboundary_idx(bottomboundary_idx == 0) = Inf; % 7.3
bottomboundary_idx = min(bottomboundary_idx, [], 1); % 7.4
bottomboundary_idx = squeeze(bottomboundary_idx);

% Clear large variables we don't need for the full stack anymore
clearvars row_idx radon_thr mask

disp('Radon thresholding end')

% Store Statistics results
tirs_result.radon_size = size(radon_stack); % Store size for reconstruction
tirs_result.idx_maxloc = maxlocarray;
tirs_result.idx_uploc = upperboundary_idx;
tirs_result.idx_downloc = bottomboundary_idx;
tirs_result.diameter = bottomboundary_idx-upperboundary_idx;
tirs_result.center_loc = bottomboundary_idx+upperboundary_idx;
tirs_result.median_diameter = median(tirs_result.diameter,2);
tirs_result.normalized_diameterchange = tirs_result.diameter./tirs_result.median_diameter;
tirs_result.normalized_diameterchange = tirs_result.normalized_diameterchange-1;
tirs_result.var_normdiameter = var(tirs_result.normalized_diameterchange,0,2);
tirs_result.var_diameter = var(tirs_result.diameter,0,2);


disp('Inverse radon transform end')

% Gather all data to CPU
tirs_result = gather_radon_result(tirs_result);
if nargout >= 2
    sinogram_stack = gather(radon_stack);
end

end

function result = gather_radon_result(result)
% GATHER_RADON_RESULT Bring all GPU arrays to CPU for saving.
fields = fieldnames(result);
for i = 1:numel(fields)
    fname = fields{i};
    val = result.(fname);
    if isa(val, 'gpuArray')
        result.(fname) = gather(val);
    end
end
end
