function [top5_framecell,top5legnth] = analyze_prctileframecluster(bv_diameter_area_array, percentile_range)
%ANALYZE_PRCTILEFRAMECLUSTER Summary of this function goes here
%   Detailed explanation goes here
    lower_percentile = prctile(bv_diameter_area_array,percentile_range(1));
    upper_percentile = prctile(bv_diameter_area_array,percentile_range(2));
    thresholded_areaidx = find(bv_diameter_area_array > lower_percentile & bv_diameter_area_array < upper_percentile);

    binary_differential = diff(thresholded_areaidx);
    binary_differential = binary_differential < 3; % binarize with neighboring frames with 2 frame tolerance

    neighboring_idxs = bwconncomp(binary_differential); 
    length_idxs = cellfun(@numel,neighboring_idxs.PixelIdxList);
    [top5legnth, top5idxs] = maxk(length_idxs(2:end-1),5); % remove beginning and end part
    top5idxs = top5idxs+1;
    top5_idxlists = neighboring_idxs.PixelIdxList(top5idxs);
    top5_framecell = cellfun(@(idx) thresholded_areaidx(idx), top5_idxlists, 'UniformOutput', false);
end

