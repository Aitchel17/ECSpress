function [idx, kymographmask] = analyze_csfoutter(kymograph, bv_upidx, bv_downidx, threshold, offset)
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
    up_offset = prctile(edgebv_upkymograph, offset, 1);
    mask_up = edgebv_upkymograph < up_offset;
    [~, upoffsetloc] = max(mask_up .* row_idx_grid, [], 1);
    
    % Upper FWHM calculation (CSF internal boundary using vessel center)
    centbv_upkymograph = kymograph;
    centbv_upkymograph(row_idx_grid > bv_centeridx | row_idx_grid < upoffsetloc) = NaN;
    centbv_upkymograph = centbv_upkymograph - up_offset; % don't care about negative value will not effect on thresholding
    maxup = max(centbv_upkymograph, [], 1);
    centbv_upkymograph = centbv_upkymograph ./ maxup;
    thr_cup = centbv_upkymograph > threshold;
    csf_thridx = row_idx_grid .* thr_cup;
    csf_thridx(csf_thridx == 0) = Inf;
    upperupboundary_idx = min(csf_thridx, [], 1);
    
    % Bottom process 
    % CSF external boundary using vessel bottom boundary
    edgebv_downkymograph = kymograph; % initialize bottom kymograph
    edgebv_downkymograph(row_idx_grid < bv_downidx) = NaN; % 
    down_offset = prctile(edgebv_downkymograph, offset, 1); % lower 10% of down kymograph 
    mask_bot = edgebv_downkymograph < down_offset; % make mask of lower 10% intensity at the bottom side
    mask_bot = mask_bot .* row_idx_grid;  % make row index array
    mask_bot(mask_bot == 0) = Inf; % to find top most lower 10% intensity position, maked region changed to Inf 
    [~, downoffsetloc] = min(mask_bot, [], 1); % find top most lower 10% intenisty position
    
    centbv_downkymograph = kymograph;  % 1.1 Initialize array
    centbv_downkymograph(row_idx_grid < bv_centeridx | row_idx_grid > downoffsetloc) = NaN;  % 1.2 Remove upper portion of image from center of vessel
    kymographmask.csf_down = centbv_downkymograph;
    centbv_downkymograph = centbv_downkymograph - down_offset;
    maxdown = max(centbv_downkymograph, [], 1);
    centbv_downkymograph = centbv_downkymograph ./ maxdown;
    thr_cdown = centbv_downkymograph > threshold;
    csf_thrdownidx = thr_cdown .* row_idx_grid;
    lowerlowboundary_idx = max(csf_thrdownidx, [], 1);
    mask = false(sz); % 8.1
    
    
    uplocarray2d = repmat(upperupboundary_idx, [sz(1), 1]); % 8.2
    mask(row_idx_grid >= uplocarray2d) = 1; % 8.3
    mask(row_idx_grid > bv_upidx) =0; % 8.5
    mask(row_idx_grid >= uplocarray2d) = 1; % 8.3
    mask(row_idx_grid > bv_downidx) =1; % 8.5
    lowerlocarray2d = repmat(lowerlowboundary_idx, [sz(1), 1]); % 8.2
    mask(row_idx_grid >= lowerlocarray2d) = 0; % 8.3
    upline = false(sz); % 8.1
    upline(row_idx_grid == uplocarray2d) = 1;
    downline = false(sz); % 8.1
    downline(row_idx_grid == lowerlocarray2d) = 1;
    
    % output
    idx.pvs_upmax = maxup;
    idx.pvs_downmax = maxdown;
    idx.pvs_upperboundary = upperupboundary_idx;
    idx.pvs_lowerboundary = lowerlowboundary_idx;
    idx.pvs_upshadeloc = upoffsetloc;
    idx.pvs_downshadeloc = downoffsetloc;

    kymographmask.pvs_upline = upline;
    kymographmask.pvs_downline = downline;
    kymographmask.pvs_bwin = mask;
    kymographmask.pvs_up = edgebv_upkymograph; % upper csf kymograph 
    kymographmask.pvs_down = edgebv_downkymograph; % lower csf kymograph 

end

% Debug code
% figure()
% imagesc(c_up)
% caxis([-2 1])

