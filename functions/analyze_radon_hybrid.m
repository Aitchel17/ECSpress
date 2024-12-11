function [area,radoncontours] = analyze_radon_hybrid(inputStack)
% Set roi
[vertices,~] = roi_rectangle(inputStack);
vertices = round(vertices);
param.the_angles=gpuArray(1:1:180);
rtd_threshold = 0.5;
irtd_threshold = 0.2;


% Extract the selected region from the stack
hold_stack = inputStack(vertices(1,2):vertices(3,2), vertices(1,1):vertices(3,1), :);
hold_stack = gpuArray(hold_stack);
submean_stack = hold_stack - mean(hold_stack,[1,2]);
radon_stack = zeros([size(radon(submean_stack(:,:,1),param.the_angles)),size(submean_stack,3)]);
hold_stack = gather(hold_stack);

% matrix calculation, radon transform and normalization
disp('Radon transform start')
for idx = 1:size(submean_stack,3)
    radon_stack(:,:,idx)=radon(submean_stack(:,:,idx),param.the_angles);
end
disp('Radon transform end')

submean_stack = gather(submean_stack);
radon_subtract = radon_stack-min(radon_stack,[],1);
radon_norm = radon_subtract./max(radon_subtract,[],1);

% clear gpu memory
clear radon_stack radon_subtract;

% % Radon thresholding
% [~, maxlocarray] = max(radon_norm,[],1);
% maxlocarray = squeeze(maxlocarray);
% maxlocarray = gather(maxlocarray);
% sz = size(radon_norm); % Make column maximum mask
% 
% % Initialize the mask array with zeros
% mask = zeros(sz);
% [row_idx, ~, ~] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
% % Expand maxlocarray to match the size of radon_norm
% maxlocarray3d = repmat(reshape(maxlocarray, [1, sz(2), sz(3)]), [sz(1), 1, 1]);
% % Create the mask by comparing row indices with maxlocarray3d
% mask(row_idx <= maxlocarray3d) = 1;
% 
% % until here I confirmed by directly subtract with original result.
% % Thresholding... from here deviate
% % upper processing
% radon_thr = radon_norm<param.rtd_threshold;
% upper_radon_thr = radon_thr .* mask;
% upperboundary_idx = row_idx .* upper_radon_thr;
% upperboundary_idx = max(upperboundary_idx, [], 1);
% upperboundary_idx = squeeze(upperboundary_idx);
% % bottom processing
% mask = ~mask;
% bottom_radon_thr = radon_thr.*mask;
% bottomboundary_idx = row_idx .* bottom_radon_thr;
% bottomboundary_idx(bottomboundary_idx == 0) = 9999;
% bottomboundary_idx = min(bottomboundary_idx, [], 1);
% bottomboundary_idx = squeeze(bottomboundary_idx);
% % deviation occured
% 
% % Final mask processing
% uplocarray3d = repmat(reshape(upperboundary_idx, [1, sz(2), sz(3)]), [sz(1), 1, 1]);
% downlocarray3d = repmat(reshape(bottomboundary_idx, [1, sz(2), sz(3)]), [sz(1), 1, 1]);
% rtd_mask = zeros(sz);
% rtd_mask(row_idx >= uplocarray3d) = 1;
% rtd_mask(row_idx > downlocarray3d) = 0;
% 
% tirs = radon_norm.*rtd_mask; 

% 
% % inverse radon tranform
% irtd_norm = zeros([size(iradon(double(radon_norm(:,:,1)),(tmp.the_angles),'linear','Hamming',size(radon_norm,2))),size(radon_norm,3)]);
% for idx = 1:size(radon_norm,3)
%     tirs_slice = double(radon_norm(:,:,idx)>param.rtd_threshold*max(tirs(:,:,idx)));
%     irtd_norm(:,:,idx)=iradon(tirs_slice,(tmp.the_angles),'linear','Hamming',size(radon_norm,2));
% end
disp('Radon thresholding start')

% Hybrid thresholding in radon space
irtd_stack = zeros([size(radon_norm,2),size(radon_norm,2),size(radon_norm,3)]);

% find the threshold crossings in Radon Space
for f = 1:size(radon_norm,3)
    disp(f)
    radon_hold_image = radon_norm(:,:,f);
    for k = 1:length(param.the_angles)
        [maxpoint(k),maxpointlocation(k)]=max(radon_hold_image(:,k));
        [~,min_edge(k)]=max(find(radon_hold_image(1:maxpointlocation(k),k)<rtd_threshold));
        [~,max_edge(k)]=max(find(radon_hold_image(maxpointlocation(k)+1:end,k)>rtd_threshold));
        
        % set all other pixels to 0
        radon_hold_image(1:min_edge(k),k)=0;
        radon_hold_image((maxpointlocation(k)+max_edge(k)):end,k)=0;
    end
    radon_norm(:,:,f) = radon_hold_image;
    tirs_slice = double(radon_hold_image>rtd_threshold*max(radon_hold_image));
    irtd_norm = iradon(tirs_slice,(param.the_angles),'linear','Hamming',size(radon_norm,2));
    irtd_stack(:,:,idx)= irtd_norm;
    [cc,l] = bwboundaries(irtd_norm>irtd_threshold*max(irtd_norm(:)));
    numPixels = cellfun(@length,cc);
    [~,idx] = max(numPixels);
    area_filled=regionprops(l,'FilledArea','Image','FilledImage');
    area(f)=length(find(area_filled(idx).FilledImage));
    radoncontours{f}=contour(irtd_norm(1:end,1:end),[irtd_threshold irtd_threshold]*max(irtd_norm(:)),'r', 'LineWidth', 2);
end
end


