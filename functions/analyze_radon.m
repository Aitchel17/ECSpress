function [result] = analyze_radon(input_Stack)

% Parameter for threshold in radon space
the_angles=1:1:180;
rtd_threshold = 0.5;
irtd_threshold = 0.2;

[result.vertices, ~] = roi_rectangle_polygon(input_Stack,'rectangle');

result.vertices = round(result.vertices);

% Extract the selected region from the stack
result.hold_stack = input_Stack(result.vertices(1,2):result.vertices(3,2), result.vertices(1,1):result.vertices(3,1), :);
result.hold_stack = gpuArray(result.hold_stack);
submean_stack = result.hold_stack - mean(result.hold_stack,[1,2]);
radon_stack = zeros([size(radon(submean_stack(:,:,1),the_angles)),size(submean_stack,3)]);
result.hold_stack = gather(result.hold_stack);
% matrix calculation, radon transform and normalization
disp('Radon transform start')
for idx = 1:size(submean_stack,3)
    radon_stack(:,:,idx)=radon(submean_stack(:,:,idx),the_angles);
end
disp('Radon transform end')

radon_subtract = radon_stack-min(radon_stack,[],1);
radon_norm = radon_subtract./max(radon_subtract,[],1);

% clear gpu memory
clearvars submean_stack radon_stack radon_subtract;


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
disp('Radon thresholding start')
% upper processing
radon_thr = radon_norm<rtd_threshold;
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
rtd_mask = zeros(sz);
uplocarray3d = repmat(reshape(upperboundary_idx, [1, sz(2), sz(3)]), [sz(1), 1, 1]);
rtd_mask(row_idx >= uplocarray3d) =1;
clearvars uplocarray3d
downlocarray3d = repmat(reshape(bottomboundary_idx, [1, sz(2), sz(3)]), [sz(1), 1, 1]);
rtd_mask(row_idx > downlocarray3d) =0;
clearvars downlocarray3d

result.tirs = radon_norm.*rtd_mask;
clearvars rtd_mask radon_norm
result.tirs = gpuArray(result.tirs);
result.irtd_norm = zeros([size(iradon(double(result.tirs(:,:,1)),(the_angles),'linear','Hamming',size(result.tirs,2))),size(result.tirs,3)]);

disp('Radon thresholding end')

disp('Inverse radon transform start')

%
for f = 1:size(result.irtd_norm,3)
    tirs_slice = double(result.tirs(:,:,f)>rtd_threshold*max(result.tirs(:,:,f)));
    result.irtd_norm(:,:,f)=iradon(tirs_slice,(the_angles),'linear','Hamming',size(result.tirs,2));
end

result.tirs = gather(result.tirs);
result.irtd_norm = gather(result.irtd_norm);

disp('Inverse radon transform end')
disp('Area calculation and contour start')

result.area = zeros([size(result.irtd_norm,3),1]);
result.radoncontours=cell(size(result.irtd_norm,3),1); %countour lines obtained of the vessel lumen using TiRS method
for f = 1:size(result.irtd_norm,3)
    [cc,l]=bwboundaries(result.irtd_norm(:,:,f)>irtd_threshold*max(result.irtd_norm(:,:,f),[],'all'));
    numPixels = cellfun(@length,cc);
    [~,idx] = max(numPixels);
    area_filled=regionprops(l,'FilledArea','Image','FilledImage');
    result.area(f)=length(find(area_filled(idx).FilledImage));
    result.radoncontours{idx} = contour(result.irtd_norm(:,:,idx),[irtd_threshold irtd_threshold]*max(result.irtd_norm(:,:,idx),[],'all'),'r', 'LineWidth', 2);
end
disp('Area calculation and contour end')


end
% Original ploting step
% figure()
% for idx = 1:size(result.irtd_norm,3)
% subplot(111)
%     % hold_image = result.hold_stack(:,:,idx);
%     % imagesc(hold_image,'XData',[1:size(hold_image,2)]+90-size(hold_image,2)/2,'YData',[1:size(hold_image,1)]+90-size(hold_image,1)/2)
%     % axis([min([1:size(hold_image,2)]+90-size(hold_image,2)/2) max([1:size(hold_image,2)]+90-size(hold_image,2)/2)...
%     % min([1:size(hold_image,1)]+90-size(hold_image,1)/2) max([1:size(hold_image,1)]+90-size(hold_image,1)/2)])
%     % axis equal
%     % axis xy
%     % hold on
%  disp(idx)
%     %subplot(111)
%  radoncontours{idx}=contour(result.irtd_norm(:,:,idx),[irtd_threshold irtd_threshold]*max(result.irtd_norm(:,:,idx),[],'all'),'r', 'LineWidth', 2);
%  %pause(0.01)
% end
% %%
% numPixels_contours = cellfun(@length,radoncontours);
% %%
% figure()
% plot(numPixels_contours)
