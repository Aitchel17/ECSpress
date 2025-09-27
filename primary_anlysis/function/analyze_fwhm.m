function [idx, kymograph_mask] = analyze_fwhm(kymograph,threshold,offset)
arguments
    kymograph (:,:) {mustBeNumeric}
    threshold (1,1)  double {mustBeGreaterThanOrEqual(threshold,0), mustBeLessThanOrEqual(threshold,1)}
    offset (1,1) double {mustBeGreaterThanOrEqual(offset,0), mustBeLessThanOrEqual(offset,100)}
end

    [max_value, max_idx] = max(kymograph,[],1); % 5.1 
    sz = size(kymograph);
    norm_kymograph = kymograph./max_value;
    [row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 5.3

    filtered_max_idx = medfilt1(max_idx,31);
    %devidation from midline
    % dev_mid = sz(1)/2 - filtered_max_idx;
    % mid_idx = (sz(1)-dev_mid)./2 ;
    maxlocarray2d = repmat(filtered_max_idx, [sz(1), 1]); % 5.4 
    upkymograph = norm_kymograph;
    upkymograph(row_idx_grid > maxlocarray2d) = NaN;
    up_offset = prctile(upkymograph,offset,1,"Method","approximate");
    upoffsetloc =  upkymograph<up_offset;
    upoffsetloc = max(row_idx_grid.*upoffsetloc,[],1);
    upoffsetloc = medfilt1(upoffsetloc,31);
    upkymograph(row_idx_grid < upoffsetloc) = NaN;
    upkymograph = upkymograph - up_offset; % below offset become negative
    upkymograph = upkymograph./max(upkymograph,[],1); % Max to be 1
    upkymograph_thr = upkymograph > threshold; % 6.2 thresholding
    upperboundary_idx = row_idx_grid .* (upkymograph_thr & row_idx_grid <= maxlocarray2d); % 6.3
    upperboundary_idx(upperboundary_idx==0) = Inf;
    upperboundary_idx = min(upperboundary_idx, [], 1); % 6.4
    upperboundary_idx = squeeze(upperboundary_idx);
%%
    % 7. bottom proessing
    downkymograph = kymograph;
    downkymograph(row_idx_grid < maxlocarray2d) = NaN;
    %%
    down_offset = prctile(downkymograph,offset,1,"Method","approximate");
    downoffsetloc =  downkymograph<down_offset;

    %%
    downoffsetloc = row_idx_grid.*downoffsetloc;
    %%
    downoffsetloc(downoffsetloc == 0) = Inf;
    %%
    downoffsetloc = min(downoffsetloc,[],1);
    downoffsetloc = medfilt1(downoffsetloc,31);

    %%
    downkymograph(row_idx_grid > downoffsetloc) = NaN;
    %%
    downkymograph = downkymograph - down_offset;
    downkymograph = downkymograph./max(downkymograph,[],1);
    %%
    downkymograph_thr = downkymograph>threshold; % 6.2
    %%
    lowerboundary_idx = row_idx_grid.*downkymograph_thr; % 7.2
    %%
    %%lowerboundary_idx(lowerboundary_idx == 0) = Inf; % 7.3
    lowerboundary_idx = max(lowerboundary_idx, [], 1); % 7.4
    lowerboundary_idx = squeeze(lowerboundary_idx);
    %%
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
    idx.upperboundary = upperboundary_idx;
    idx.lowerboundary = lowerboundary_idx;
    idx.max_idx = filtered_max_idx;
    idx.downoffsetloc = downoffsetloc;
    idx.upoffsetloc = upoffsetloc;
    kymograph_mask = [];
    kymograph_mask.rowidx = row_idx_grid;
    kymograph_mask.bwin = mask;
    kymograph_mask.up = upkymograph;
    kymograph_mask.down = downkymograph;
    kymograph_mask.maxline = maxline; 
    kymograph_mask.upline = upline;
    kymograph_mask.downline = downline;
    kymograph_mask.midline = midline;
end

