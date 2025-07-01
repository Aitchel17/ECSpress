function [idx, kymographmask] = analyze_csfoutter(kymograph, bv_upidx, bv_downidx, threshold)
% calculate_csf_boundaries - calculates upper and lower CSF boundaries based on given inputs.
%
% Inputs:
%   c                  - normalized CSF intensity profile (2D array)
%   upperboundary_idx  - upper boundary indices from blood vessel
%   bottomboundary_idx - bottom boundary indices from blood vessel
%   bv_centeridx       - blood vessel center indices
%   threshold          - intensity threshold for boundary detection
%   thresholded        - initial thresholded binary mask
%
% Outputs:
%   up_csf_boundary    - upper CSF boundary indices
%   down_csf_boundary  - lower CSF boundary indices
%   cb_thresholded     - updated binary mask with CSF boundaries

% output struct
idx = []; % 1D row index per slice
kymographmask = []; % 2D Mask corresponding to kymograph or masked kymograph
sz = size(kymograph);
[row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 2D row index grid
bv_centeridx = ceil((bv_downidx+bv_upidx)/2);

% Top process (CSF external boundary using vessel upper boundary)
edgebv_upkymograph = kymograph;
edgebv_upkymograph(row_idx_grid > bv_upidx) = NaN;

% Outer boundary detection (upper)
up_offset = prctile(edgebv_upkymograph, 10, 1);
mask_up = edgebv_upkymograph < up_offset;
[~, upoffsetloc] = max(mask_up .* row_idx_grid, [], 1);

% Upper FWHM calculation (CSF internal boundary using vessel center)
centbv_upkymograph = kymograph;
centbv_upkymograph(row_idx_grid > bv_centeridx | row_idx_grid < upoffsetloc) = NaN;
centbv_upkymograph = centbv_upkymograph - up_offset; % don't care about negative value will not effect on thresholding
idx.maxup = max(centbv_upkymograph, [], 1);
centbv_upkymograph = centbv_upkymograph ./ max(centbv_upkymograph, [], 1);
thr_cup = centbv_upkymograph > threshold;
csf_thridx = row_idx_grid .* thr_cup;
csf_thridx(csf_thridx == 0) = Inf;
upperupboundary_idx = min(csf_thridx, [], 1);

% Bottom process 
% CSF external boundary using vessel bottom boundary
edgebv_downkymograph = kymograph; % initialize bottom kymograph
edgebv_downkymograph(row_idx_grid < bv_downidx) = NaN; % 
down_offset = prctile(edgebv_downkymograph, 10, 1); % lower 10% of down kymograph 
mask_bot = edgebv_downkymograph < down_offset; % make mask of lower 10% intensity at the bottom side
mask_bot = mask_bot .* row_idx_grid;  % make row index array
mask_bot(mask_bot == 0) = Inf; % to find top most lower 10% intensity position, maked region changed to Inf 
[~, botoffsetloc] = min(mask_bot, [], 1); % find top most lower 10% intenisty position

centbv_downkymograph = kymograph;  % 1.1 Initialize array
centbv_downkymograph(row_idx_grid < bv_centeridx | row_idx_grid > botoffsetloc) = NaN;  % 1.2 Remove upper portion of image from center of vessel
kymographmask.testdown = centbv_downkymograph;
centbv_downkymograph = centbv_downkymograph - down_offset;
idx.maxdown = max(centbv_downkymograph, [], 1);
centbv_downkymograph = centbv_downkymograph ./ max(centbv_downkymograph, [], 1);
thr_cdown = centbv_downkymograph > threshold;
csf_thrdownidx = thr_cdown .* row_idx_grid;
lowerlowboundary_idx = max(csf_thrdownidx, [], 1);
mask = false(sz); % 8.1
thresholded = false(sz); % 8.1

uplocarray2d = repmat(upperupboundary_idx, [sz(1), 1]); % 8.2
mask(row_idx_grid >= uplocarray2d) = 1; % 8.3
mask(row_idx_grid > bv_upidx) =0; % 8.5
mask(row_idx_grid >= uplocarray2d) = 1; % 8.3
mask(row_idx_grid > bv_downidx) =1; % 8.5
lowerlocarray2d = repmat(lowerlowboundary_idx, [sz(1), 1]); % 8.2
mask(row_idx_grid >= lowerlocarray2d) = 0; % 8.3
thresholded(row_idx_grid == lowerlocarray2d) = 1;
thresholded(row_idx_grid == uplocarray2d) = 1;

idx.upboundary = upperupboundary_idx;
idx.downboundary = lowerlowboundary_idx;

kymographmask.boundary = thresholded;
kymographmask.mask = mask;

kymographmask.up = edgebv_upkymograph;
kymographmask.down = edgebv_downkymograph;

end

% Debug code
% figure()
% imagesc(c_up)
% caxis([-2 1])

