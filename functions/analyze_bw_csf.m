function [bw_csf,outputArg2] = analyze_bw_csf(target_stack,thr_refstack,sigma)
%ANALYZE_BW_CSF Summary of this function goes here
 % 1. Using constricted blood vessel mask, calculate mean + n*sigma for
 % threshold of each slice and apply.  bw_csf
 % 2. 

bw_csf = false(size(target_stack,3));       % Preallocate binary output

%   Detailed explanation goes here
for sli = 1:size(target_stack,3)
    slice = target_stack(:,:,sli);
    % Flatten and remove NaNs
    data = slice(~isnan(thr_refstack(:,:,sli)));
    % % Trim top 10%
    % upper = prctile(data, 50);
    % trimmed_data = data(data <= upper);
    trimmed_data = data;
    % Fit normal distribution
    pd = fitdist(trimmed_data, 'Normal');
    threshold = pd.mu + sigma*pd.sigma;
    % Apply threshold to current slice
    bw_csf(:,:,sli) = slice > threshold;
end
se = strel('square', 3);   % 2D structuring element
bw3D_clean = false(size(bw_csf));  % Preallocate

for z = 1:size(bw_csf, 3)
    % Apply 2D operations per slice
    slice = bw_csf(:,:,z);
    slice = bwareaopen(slice,2); % clean the pixels below 2
    imclose(slice, se);  % closing + opening
    slice_clean = imclose(slice,se);
    bw3D_clean(:,:,z) = slice_clean;
end
end

