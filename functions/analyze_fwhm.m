function [idx, kymograph_mask] = analyze_fwhm(kymograph,threshold,offset)
arguments
    kymograph (:,:) {mustBeNumeric}
    threshold (1,1)  double {mustBeGreaterThanOrEqual(threshold,0), mustBeLessThanOrEqual(threshold,1)}
    offset (1,1) double {mustBeGreaterThanOrEqual(offset,0), mustBeLessThanOrEqual(offset,100)}
end

    [~, max_idx] = max(kymograph,[],1); % 5.1 
    sz = size(kymograph);

    [row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 5.3
    maxlocarray2d = repmat(max_idx, [sz(1), 1]); % 5.4
    % 6. upper processing

    upkymograph = kymograph;
    upkymograph(row_idx_grid > maxlocarray2d) = NaN;
    upkymograph = upkymograph - prctile(upkymograph,offset,1,"Method","approximate");
    upkymograph = upkymograph./max(upkymograph,[],1);
    upkymograph_thr = upkymograph < threshold; % 6.2
    upperboundary_idx = row_idx_grid .* (upkymograph_thr & row_idx_grid <= maxlocarray2d); % 6.3
    upperboundary_idx = max(upperboundary_idx, [], 1); % 6.4
    upperboundary_idx = squeeze(upperboundary_idx);
    % 7. bottom proessing
    downkymograph = kymograph;
    downkymograph(row_idx_grid < maxlocarray2d) = NaN;
    downkymograph = downkymograph - prctile(downkymograph,offset,1,"Method","approximate");
    downkymograph = downkymograph./max(downkymograph,[],1);
    downkymograph_thr = downkymograph<threshold; % 6.2
    lowerboundary_idx = row_idx_grid.*(downkymograph_thr & row_idx_grid >= maxlocarray2d); % 7.2
    lowerboundary_idx(lowerboundary_idx == 0) = Inf; % 7.3
    lowerboundary_idx = min(lowerboundary_idx, [], 1); % 7.4
    lowerboundary_idx = squeeze(lowerboundary_idx);
    % 8. Final mask processing
    uplocarray2d = repmat(upperboundary_idx, [sz(1), 1]); % 8.2
    downlocarray2d = repmat(lowerboundary_idx, [sz(1), 1]); % 8.4
    mask = row_idx_grid >= uplocarray2d & row_idx_grid <= downlocarray2d;

    maxline = false(sz);
    maxline(row_idx_grid == maxlocarray2d) = 1;
    maxline(row_idx_grid == uplocarray2d | row_idx_grid == downlocarray2d) = 1;

    upline = false(sz);
    upline(row_idx_grid == uplocarray2d) = 1;
    
    downline = false(sz);
    downline(row_idx_grid == downlocarray2d) = 1;

    midline = false(sz);
    midline(row_idx_grid == ceil((downlocarray2d+uplocarray2d)/2)) = 1;

    % 9. return
    idx =[];
    idx.bv_upperboundary = upperboundary_idx;
    idx.bv_lowerboundary = lowerboundary_idx;
    idx.max_idx = max_idx;
    
    kymograph_mask = [];
    kymograph_mask.rowidx = row_idx_grid;
    kymograph_mask.bv_bwin = mask;
    kymograph_mask.bv_up = upkymograph;
    kymograph_mask.bv_down = downkymograph;
    kymograph_mask.bv_maxline = maxline; 
    kymograph_mask.bv_upline = upline;
    kymograph_mask.bv_downline = downline;
    kymograph_mask.bv_midline = midline;
end

