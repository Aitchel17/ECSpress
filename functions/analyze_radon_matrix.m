function [outputArg1,outputArg2] = analyze_radon_matrix(inputArg1,inputArg2)
[vertices,~] = roi_rectangle(finalFilteredStack);
param.xy = round(vertices);


% Extract the selected region from the stack
hold_stack = finalFilteredStack(param.xy(1,2):param.xy(3,2), param.xy(1,1):param.xy(3,1), :);
hold_stack = gpuArray(hold_stack);
submean_stack = hold_stack - mean(hold_stack,[1,2]);
radon_stack = zeros([size(radon(submean_stack(:,:,1),tmp.the_angles)),size(submean_stack,3)]);
hold_stack = gather(hold_stack);
% matrix calculation, radon transform and normalization
disp('Radon transform start')
for idx = 1:size(submean_stack,3)
    radon_stack(:,:,idx)=radon(submean_stack(:,:,idx),tmp.the_angles);
end
disp('Radon transform end')

submean_stack = gather(submean_stack);
radon_subtract = radon_stack-min(radon_stack,[],1);
radon_norm = radon_subtract./max(radon_subtract,[],1);

% clear gpu memory
clear radon_stack radon_subtract;

% Radon thresholding
[~, maxlocarray] = max(radon_norm,[],1);
maxlocarray = squeeze(maxlocarray);
maxlocarray = gather(maxlocarray);

sz = size(radon_norm); % Make column maximum mask
% Initialize the mask array with zeros
mask = zeros(sz);
[row_idx, ~, ~] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
% Expand maxlocarray to match the size of radon_norm
maxlocarray3d = repmat(reshape(maxlocarray, [1, sz(2), sz(3)]), [sz(1), 1, 1]);
% Create the mask by comparing row indices with maxlocarray3d
mask(row_idx <= maxlocarray3d) = 1;

% until here I confirmed by directly subtract with original result.
% Thresholding... from here deviate
% upper processing
radon_thr = radon_norm<param.rtd_threshold;
upper_radon_thr = radon_thr .* mask;
upperboundary_idx = row_idx .* upper_radon_thr;
upperboundary_idx = max(upperboundary_idx, [], 1);
upperboundary_idx = squeeze(upperboundary_idx);
% bottom processing
mask = ~mask;
bottom_radon_thr = radon_thr.*mask;
bottomboundary_idx = row_idx .* bottom_radon_thr;
bottomboundary_idx(bottomboundary_idx == 0) = 9999;
bottomboundary_idx = min(bottomboundary_idx, [], 1);
bottomboundary_idx = squeeze(bottomboundary_idx);
% deviation occured

% Final mask processing
uplocarray3d = repmat(reshape(upperboundary_idx, [1, sz(2), sz(3)]), [sz(1), 1, 1]);
downlocarray3d = repmat(reshape(bottomboundary_idx, [1, sz(2), sz(3)]), [sz(1), 1, 1]);

rtd_mask = zeros(sz);
rtd_mask(row_idx >= uplocarray3d) =1;
rtd_mask(row_idx > downlocarray3d) =0;

tirs = radon_norm.*rtd_mask;

irtd_norm = zeros([size(iradon(double(tirs(:,:,1)),(tmp.the_angles),'linear','Hamming',size(tirs,2))),size(tirs,3)]);

for idx = 1:size(radon_norm,3)
    tirs_slice = double(tirs(:,:,idx)>param.rtd_threshold*max(tirs(:,:,idx)));
    irtd_norm(:,:,idx)=iradon(tirs_slice,(tmp.the_angles),'linear','Hamming',size(tirs,2));
end



radoncontours = cell(size(irtd_norm,3),1); %countour lines obtained of the vessel lumen using TiRS method

for idx = 1:size(irtd_norm,3)
    disp(idx)
    radoncontours{idx} = contour(irtd_norm(:,:,idx),[0.2 0.2]*max(irtd_norm(:)),'r','LineWidth',2);
end


radon_norm = gather(radon_norm);
irtd_norm = gather(irtd_norm);
tirs = gather(tirs);
toc
end
end

