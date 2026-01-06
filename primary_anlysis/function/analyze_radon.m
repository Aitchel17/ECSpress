function radon_result = analyze_radon(hold_stack)

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
radon_result.radon_size = size(radon_stack); % Store size for reconstruction
radon_result.idx_maxloc = maxlocarray;
radon_result.idx_uploc = upperboundary_idx;
radon_result.idx_downloc = bottomboundary_idx;
radon_result.diameter = bottomboundary_idx-upperboundary_idx;
radon_result.median_diameter = median(radon_result.diameter,2);
radon_result.normalized_diameterchange = radon_result.diameter./radon_result.median_diameter;
radon_result.normalized_diameterchange = radon_result.normalized_diameterchange-1;
radon_result.var_normdiameter = var(radon_result.normalized_diameterchange,0,2);
radon_result.var_diameter = var(radon_result.diameter,0,2);


disp('Inverse radon transform end')

% Gather all data to CPU
radon_result = gather_radon_result(radon_result);

end

function radon_result = gather_radon_result(radon_result)
% GATHER_RADON_RESULT Bring all GPU arrays to CPU for saving.
fields = fieldnames(radon_result);
for i = 1:numel(fields)
    fname = fields{i};
    val = radon_result.(fname);
    if isa(val, 'gpuArray')
        radon_result.(fname) = gather(val);
    end
end
end
