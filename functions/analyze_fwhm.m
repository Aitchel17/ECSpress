function [idx, kymograph_mask] = analyze_fwhm(kymograph,threshold,offset)
arguments
    kymograph (:,:) {mustBeNumeric}
    threshold (1,1)  double {mustBeGreaterThanOrEqual(threshold,0), mustBeLessThanOrEqual(threshold,1)}
    offset (1,1) double {mustBeGreaterThanOrEqual(offset,0), mustBeLessThanOrEqual(offset,100)}
end

    [~, max_idx] = max(kymograph,[],1); % 5.1 
    sz = size(kymograph);
    thresholded = false(sz);
    [row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 5.3
    maxlocarray2d = repmat(max_idx, [sz(1), 1]); % 5.4
    % 6. upper processing
    mask = false(sz); % 5.2
    mask(row_idx_grid <= maxlocarray2d) = 1; % 6.1
    thresholded(row_idx_grid == maxlocarray2d) = 1;
    upkymograph = kymograph;
    upkymograph(~mask) = NaN;
    upkymograph = upkymograph - prctile(upkymograph,offset,1,"Method","approximate");
    upkymograph = upkymograph./max(upkymograph,[],1);
    upkymograph_thr = upkymograph<threshold; % 6.2
    upperboundary_idx = row_idx_grid .* upkymograph_thr.* mask; % 6.3
    upperboundary_idx = max(upperboundary_idx, [], 1); % 6.4
    upperboundary_idx = squeeze(upperboundary_idx);
    % 7. bottom proessing
    mask = false(sz); % 5.2
    mask(row_idx_grid >= maxlocarray2d) = 1; % 6.1; % 7.1
    downkymograph = kymograph;
    downkymograph(~mask) = NaN;
    downkymograph = downkymograph - prctile(downkymograph,offset,1,"Method","approximate");
    downkymograph = downkymograph./max(downkymograph,[],1);
    downkymograph_thr = downkymograph<threshold; % 6.2
    lowerboundary_idx = row_idx_grid.*downkymograph_thr.*mask; % 7.2
    lowerboundary_idx(lowerboundary_idx == 0) = 9999; % 7.3
    lowerboundary_idx = min(lowerboundary_idx, [], 1); % 7.4
    lowerboundary_idx = squeeze(lowerboundary_idx);
    % 8. Final mask processing
    mask = false(sz); % 8.1
    uplocarray2d = repmat(upperboundary_idx, [sz(1), 1]); % 8.2
    mask(row_idx_grid >= uplocarray2d) =1; % 8.3
    thresholded(row_idx_grid == uplocarray2d) = 1;
    downlocarray2d = repmat(lowerboundary_idx, [sz(1), 1]); % 8.4
    mask(row_idx_grid > downlocarray2d) =0; % 8.5
    thresholded(row_idx_grid == downlocarray2d) = 1; % 
    
    % 9. saving
    idx =[];
    idx.upperboundary_idx = upperboundary_idx;
    idx.lowerboundary_idx = lowerboundary_idx;
    idx.max_idx = max_idx;

    kymograph_mask = [];
    kymograph_mask.rowidx = row_idx_grid;
    kymograph_mask.inside = mask;
    kymograph_mask.boundary = thresholded; 

end

