function result = heatmap_postprocessing(bin_counts, x_edge, y_edge, x_origin, y_origin)
%HEATMAP_POSTPROCESSING  The heatmap's main blob, put on the caller's axis and measured.
%   Takes histcounts2's three outputs AS THEY COME, so a caller passes them straight
%   through. The transpose and the edges-to-centres step both live here: taking the
%   edges in histcounts2's layout and the counts in some other would be a contract
%   that reads right and is wrong.
%
%   Everything comes back in the unit the EDGES were given in, so scale them on the
%   way in and nothing downstream needs scaling. The binning itself must still be in
%   pixels -- both series are differences of row indices, and a micrometre-wide bin
%   lets an edge fall between two allowed values and rules empty lines across.
%
%   IN   bin_counts  nx x ny double  histcounts2's counts, x down ROWS
%        x_edge      1 x nx+1 double bin edges, in whatever unit is wanted out
%        y_edge      1 x ny+1 double
%        x_origin    1x1 double      where zero goes, same unit as the edges. The
%                                    caller owns this: two maps of one vessel share
%                                    an anchor only if something outside both picks
%                                    it, and without that their axes cannot be read
%                                    against each other at all
%        y_origin    1x1 double
%   OUT  result     1x1 struct   xy_counts_clean   ny' x nx'  NaN outside the blob,
%                                                  all-NaN rows and columns trimmed
%                                x_baseceneters    1 x nx'  x_centers minus x_origin
%                                y_baseceneters    1 x ny'
%                                x_mode            1x1  the most-occupied column, on
%                                                  the shifted axis. Where THIS map
%                                                  rests, against the shared zero
%                                y_mode            1x1
%                                modepvs           1 x nx'  each column's most-occupied
%                                                  row, read in y_baseceneters
%                                column_std        1x1  spread WITHIN a column
%                                mode_slope        1x1  robustfit over the whole curve
%                                slope_constricted 1x1  fitted where x < 0
%                                slope_dilated     1x1  where x > 0, on the same span
%                                n_constricted     1x1  columns that fit was made on
%                                n_dilated         1x1
%                                bend              1x1  dilated minus constricted

% Into this file's own layout: y down rows, and one centre per bin
xy_counts = bin_counts';
x_centers = (x_edge(1:end-1) + x_edge(2:end)) / 2;
y_centers = (y_edge(1:end-1) + y_edge(2:end)) / 2;

% 1. Spread inside one column, before anything is thrown away
% Each column is the distribution of y at one x, so subtracting the column's own
% median leaves only the scatter that x does not account for. Read off the counts
% AS GIVEN, and it has to stay that way: the mask below deletes scattered bins, and
% a session that scatters loses the most, so measured afterwards this number would
% be smallest exactly where the scatter is worst.
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
tmp.xy_maskedcounts(~tmp.bwmask) = NaN; % 마스크 적용
% Dev purpose, check data
% figure()
% imshow(tmp.xy_maskedcounts)

% 3. NaN removal
result.xy_counts_clean = tmp.xy_maskedcounts;
% NaN removal for x axis
tmp.allnan_x = all(isnan(result.xy_counts_clean),1);
result.xy_counts_clean = result.xy_counts_clean(:,~tmp.allnan_x);
tmp.x_centers_clean = x_centers(~tmp.allnan_x);
% NaN removal for y axis
tmp.allnan_y = all(isnan(result.xy_counts_clean),2);
result.xy_counts_clean = result.xy_counts_clean(~tmp.allnan_y,:);
tmp.y_centers_clean = y_centers(~tmp.allnan_y);

% 4. Onto the caller's anchor, and where this map's own peak sits on it
% The shift uses the given origin, NOT this map's peak. A map re-origined on its own
% peak always reads as resting at zero, so two states of one vessel would both sit at
% zero and the difference between them -- which is the thing being asked about --
% would be gone. x_mode is that difference: the peak measured on the shared axis.
result.x_baseceneters = tmp.x_centers_clean - x_origin;
result.y_baseceneters = tmp.y_centers_clean - y_origin;
tmp.baseline_bv = sum(result.xy_counts_clean,1,'omitmissing');
[~,tmp.maxbvloc]= max(tmp.baseline_bv);
result.x_mode = result.x_baseceneters(tmp.maxbvloc);
tmp.baseline_pvs = sum(result.xy_counts_clean,2,'omitmissing');
[~,tmp.maxpvsloc]= max(tmp.baseline_pvs);
result.y_mode = result.y_baseceneters(tmp.maxpvsloc);


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

