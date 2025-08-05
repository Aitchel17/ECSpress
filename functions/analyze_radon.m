function radon_result = analyze_radon(hold_stack)

% 1. Set parameter
    % 1.1 the_agnles: Determine range and interval of projection angle 
    % 1.2 rtd_threshod: Threshold after normalization, 0.5 equal full width half max of each angle 
% 2. gpuArray
    % Bring memory to gpu memory
% 3. Radon transform
    % 3.1 preallocate memory for radonimages
    % 3.2 radon transform (gpu)
% 4. Column wise radon transform before thresholding
    % 4.1 column wise minimum subtraction -> minimum to 0
    % 4.2 column wise maximum division -> maximum to 1
    % radon_stack now normalized
% 5. Make maximum location matrix
    % 5.1 column wise find location of maximum, row dimension become 1
    % 5.2 initialize zero mask for radon_stack
    % 5.3 make row index mask
    % 5.4 copy maxlocarray to row direction to match the size of radon_stack
% 6. Upper processing
    % 6.1 create the mask that fill upper position than max to be 1
    % 6.2 radon thresholding
    % 6.3 apply both radon threshold and upper mask
    % 6.4 bottom column coordination deleted by masked threshold 
    % 6.5 find lowest position of upper area using max()
            % remind: matrix idx start from 1, coordination system y inverted
% 7. Down thresholding
    % 7.1 Invert the mask
    % 7.2 apply both radon threshold and lower mask
    % 7.3 deleted index become 9999 to avoid interuption while min()  
    % 7.4 find highest boundary of lower area using min()
% 8. Final mask processing
    % 8.1 variable mask reset now its final mask
    % 8.2 same function as 5.4 but now its upper boundary 
    % 8.3 convert all lower area of upper boundary to be 1
    % 8.4 same as 5.4 but lower boundary
    % 8.5 convert all lower area of lower boundary back to 0
% 9. Inverse radon transform
radon_result = struct();

% 1. Set parameter
the_angles=1:1:180;
rtd_threshold = 0.5;

% 2. Put stack to memory of gpu
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
mask = false(sz); % 5.2
radon_result.tirs = zeros(sz);
[row_idx, ~, ~] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3)); % 5.3
maxlocarray3d = repmat(reshape(maxlocarray, [1, sz(2), sz(3)]), [sz(1), 1, 1]); % 5.4

% 6. upper processing
mask(row_idx <= maxlocarray3d) = 1; % 6.1
radon_result.tirs(row_idx == maxlocarray3d) = 2;
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

% 8. Final mask processing
mask = false(sz); % 8.1
uplocarray3d = repmat(reshape(upperboundary_idx, [1, sz(2), sz(3)]), [sz(1), 1, 1]); % 8.2
mask(row_idx >= uplocarray3d) =1; % 8.3
radon_result.tirs(row_idx == uplocarray3d) = 1;
clearvars uplocarray3d
downlocarray3d = repmat(reshape(bottomboundary_idx, [1, sz(2), sz(3)]), [sz(1), 1, 1]); % 8.4
mask(row_idx > downlocarray3d) =0; % 8.5
radon_result.tirs(row_idx == downlocarray3d) = 3;
clearvars downlocarray3d
disp('Radon thresholding end')
% mask = mask; % inspection purpose

radon_result.idx_maxloc = maxlocarray;
radon_result.idx_uploc = upperboundary_idx;
radon_result.idx_downloc = bottomboundary_idx;

% 9. Inverse radon thresholding 
disp('Inverse radon transform start')
radon_result.irtd = zeros([size(iradon(double(mask(:,:,1)),(the_angles),'linear','Hamming',1,size(mask,1))),size(mask,3)]);
mask = gpuArray(mask);
for f = 1:size(radon_result.irtd,3)
    radon_result.irtd(:,:,f)=iradon(mask(:,:,f),(the_angles),'linear','Hamming',1,size(mask,1));
end
disp('Inverse radon transform end')
end
