function result = heatmap_postprocessing(bin_counts, x_edge, y_edge, x_origin, y_origin)
%HEATMAP_POSTPROCESSING  The heatmap's main blob, put on the caller's axis.
%   histcounts2's three outputs as they come; everything returns in the unit the
%   edges were given in. Bin in pixels, scale the edges -- see CLAUDE_LOG.md
%
%   IN   bin_counts  nx x ny double  histcounts2's counts, x down ROWS
%        x_edge      1 x nx+1 double bin edges, in the unit wanted out
%        y_edge      1 x ny+1 double
%        x_origin    1x1 double      where zero goes, same unit as the edges. The
%                                    caller picks it, so two maps can share one
%        y_origin    1x1 double
%   OUT  result     1x1 struct   xy_counts_clean   ny' x nx'  NaN outside the blob,
%                                                  all-NaN rows and columns trimmed
%                                x_baseceneters    1 x nx'  x_centers minus x_origin
%                                y_baseceneters    1 x ny'
%                                x_mode            1x1  most-occupied column, on the
%                                                  shifted axis
%                                y_mode            1x1
%                                mode_curve        1 x nx'  each column's most-occupied
%                                                  row, read in y_baseceneters
%                                column_std        1x1  spread WITHIN a column

% Into this file's own layout: y down rows, and one centre per bin
xy_counts = bin_counts';
x_centers = (x_edge(1:end-1) + x_edge(2:end)) / 2;
y_centers = (y_edge(1:end-1) + y_edge(2:end)) / 2;

% 1. Spread inside one column, on the counts as given (before the mask trims them)
tmp.column_total = sum(xy_counts, 1);
tmp.cumulative = cumsum(xy_counts, 1);
tmp.past_half = tmp.cumulative >= tmp.column_total / 2;
[~, tmp.median_row] = max(tmp.past_half, [], 1);          % first row past the half
tmp.column_median = y_centers(tmp.median_row);            % 1 x nx
tmp.residual = y_centers(:) - tmp.column_median;          % ny x nx, expanded on purpose
tmp.weight_total = sum(xy_counts, 'all');
tmp.residual_mean = sum(xy_counts .* tmp.residual, 'all') / tmp.weight_total;
tmp.centred = tmp.residual - tmp.residual_mean;
result.column_std = sqrt(sum(xy_counts .* tmp.centred.^2, 'all') / tmp.weight_total);

% 2. Extract main population of heatmap
tmp.bwmask = xy_counts > 0;
tmp.cc = bwconncomp(tmp.bwmask);
tmp.stats = regionprops(tmp.cc);
tmp.area = [tmp.stats.Area];
[~,tmp.maxareaidx] = max(tmp.area);
tmp.cc.PixelIdxList{tmp.maxareaidx};
tmp.bwmask = zeros(size(tmp.bwmask));
tmp.bwmask(tmp.cc.PixelIdxList{tmp.maxareaidx}) = 1;

tmp.xy_maskedcounts = xy_counts;
tmp.xy_maskedcounts(~tmp.bwmask) = NaN;

% 3. Trim the all-NaN columns, then rows
result.xy_counts_clean = tmp.xy_maskedcounts;
tmp.allnan_x = all(isnan(result.xy_counts_clean),1);
result.xy_counts_clean = result.xy_counts_clean(:,~tmp.allnan_x);
tmp.x_centers_clean = x_centers(~tmp.allnan_x);
tmp.allnan_y = all(isnan(result.xy_counts_clean),2);
result.xy_counts_clean = result.xy_counts_clean(~tmp.allnan_y,:);
tmp.y_centers_clean = y_centers(~tmp.allnan_y);

% 4. Onto the caller's anchor, and where this map's own peak sits on it
result.x_baseceneters = tmp.x_centers_clean - x_origin;
result.y_baseceneters = tmp.y_centers_clean - y_origin;
tmp.count_by_bv = sum(result.xy_counts_clean,1,'omitmissing');
[~,tmp.maxbvloc]= max(tmp.count_by_bv);
result.x_mode = result.x_baseceneters(tmp.maxbvloc);
tmp.count_by_y = sum(result.xy_counts_clean,2,'omitmissing');
[~,tmp.maxyloc]= max(tmp.count_by_y);
result.y_mode = result.y_baseceneters(tmp.maxyloc);


% 5. Mode of PVS for each x point
[~,tmp.modepvslocs] = max(result.xy_counts_clean,[],1);
result.modepvs = result.y_baseceneters(tmp.modepvslocs);

% 6. The mode curve's slope, whole and either side of zero
% Zero on x is wherever the caller anchored, so x < 0 is constricted against that
% anchor and x > 0 dilated. Which columns land on which side therefore moves with
% the anchor, and so do these two slopes. The 2 given to fit_bothsides is not a
% threshold -- it is the fewest points a straight line has. How many columns a side
% NEEDS is the caller's judgement, and n_constricted / n_dilated is what it judges.
% robustfit needs more points than parameters AND a scale to weight by. A
% microarousal is a few seconds long and its map can come out two or three columns
% wide, which has neither -- it errored on uarousal before this was here. NaN and
% not a crash: how short is too short is the caller judgement, and n_constricted +
% n_dilated is what it judges on.
if numel(result.x_baseceneters) >= 4
    tmp.mode_fit = robustfit(result.x_baseceneters, result.modepvs);
    result.mode_slope = tmp.mode_fit(2);
else
    result.mode_slope = NaN;
end
[tmp.slope_low, tmp.slope_high, tmp.side_info] = fit_bothsides(result.x_baseceneters, ...
    result.modepvs, 2);
result.slope_constricted = tmp.slope_low;
result.slope_dilated = tmp.slope_high;
result.n_constricted = tmp.side_info.n_low;
result.n_dilated = tmp.side_info.n_high;
result.bend = tmp.side_info.bend;
end

