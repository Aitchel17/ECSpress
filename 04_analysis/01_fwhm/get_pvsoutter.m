function [idx, kymographmask] = get_pvsoutter(kymograph, bv_upidx, bv_downidx, min_prominence, sg_ord, sg_frame, min_diprows, medfilt_cols)
% calculate_csf_boundaries - calculates upper and lower CSF boundaries based on given inputs.
%
% Inputs:
%   kymograph          - normalized CSF intensity profile (2D array)
%   bv_upidx  - upper boundary indices from blood vessel
%   bv_downidx - bottom boundary indices from blood vessel
%   min_prominence  1x1 float 0-1  least prominence of a band peak or a shade dip, and the least
%                                  drop that counts as shade, on the rescaled kymograph
%   sg_ord          1x1 count      Savitzky-Golay order along a column profile
%   sg_frame        1x1 count      Savitzky-Golay frame (odd); the band-peak search skips a shorter column
%   min_diprows     1x1 count      least rows outward of the band peak for a dip search
%   medfilt_cols    1x1 count      median window across columns
%   (line_fwhm.param.pvsedge_* holds all five, so every paxfwhm.mat records them)
%
% Outputs:
%   idx.up_csf_boundary    - upper CSF boundary indices
%   idx.down_csf_boundary  - lower CSF boundary indices
%   kymographmask.cb_thresholded     - updated binary mask with CSF boundaries
%
% A column with no band peak outward of the vessel midline, or no shade beyond the peak, gets
% NaN edge boundaries; a brighter region beyond the band no longer captures the anchor, and the
% shade is the first dip beyond the peak rather than a percentile mask. The band peak may sit
% inside the vessel wall; the search runs from the peak either way, and an edge that lands inside
% the wall is set to the wall afterwards (thickness 0) rather than by cropping the search.
% Both sides run as the upper side: the lower side on the kymograph flipped upside down
% (flip_kymograph), its results flipped back. A cleaned boundary sits between rows and flips as
% n_row + 2 - b, a row as n_row + 1 - r.

%%
fig = figure('Name','Check kymograph'); % 0.1.1 initialize figure
ax = axes(fig); % 0.1.2 initialize axis
imagesc(ax,kymograph) % 0.1.3 show original kymograph for interactive crop
%%
% 0. output struct
idx = []; % 1D row index per slice
kymographmask = []; % 2D Mask corresponding to kymograph or masked kymograph

% 1. size and grid matrix generation
sz = size(kymograph);
[row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 2D row index grid
% bv_centeridx = ceil((bv_downidx+bv_upidx)/2); % center point to
% divide kymograph  251013. CSF boundary restraint


% 1. manual scope restiction
while true
    cla(ax)
    imagesc(ax,kymograph)
    upcropval = input('Get upper crop idx: '); % 0.2.1 Get input
    upcropkymograph = kymograph;

    upcropkymograph(1:upcropval,:) = NaN;
    imagesc(ax,upcropkymograph)
    flag = input('looks ok? type "y":  ', 's'); % 0.3.5 check
    if strcmp(flag,'y') % 0.3.5 check
        break
    end
end

while true
    cla(ax)
    imagesc(ax,kymograph)
    downcropval = input('Get lower crop idx: '); % 0.2.1 Get input
    cropkymograph = upcropkymograph;
    cropkymograph(downcropval:end,:) = NaN;
    imagesc(ax,cropkymograph)
    flag = input('looks ok? type "y":  ', 's'); % 0.3.5 check
    if strcmp(flag,'y') % 0.3.5 check
        break
    end
end

% 2. Rescale from cropped region
cropkymograph = cropkymograph - min(cropkymograph, [], 1);
cropkymograph = cropkymograph ./ max(cropkymograph, [], 1);

% 3. Bisect kymograph by center of vessel
bv_mididx = round((bv_upidx+bv_downidx)/2);  % 3.1 Center loc array

% 4. Each side as the upper side : the lower side runs on the kymograph flipped upside down, then flips back
n_row = sz(1);
edge_idx = nan(2, sz(2));       % 2 x W  outer PVS boundary, row 1 upper side, row 2 lower side
cent_idx = nan(2, sz(2));       % 2 x W  inner PVS boundary
edge_rawidx = nan(2, sz(2));    % 2 x W  outer half-max row before cleaning
cent_rawidx = nan(2, sz(2));    % 2 x W  inner half-max row before cleaning
for side_k = 1:2
    % 4.1 The side's kymograph, midline and wall
    if side_k == 1
        side_kgph = cropkymograph;
        side_mididx = bv_mididx;
        side_wallidx = bv_upidx;
    else
        [side_kgph, side_mididx] = flip_kymograph(cropkymograph, bv_mididx);
        side_wallidx = n_row + 2 - bv_downidx;   % compared with a cleaned boundary, which sits between rows
    end

    % 4.2 Find first peak of ECS signal from vessel center
    search_kgph = side_kgph; % peak search starts at the midline : the dark lumen keeps a band flush to the wall an interior maximum
    search_kgph(row_idx_grid > side_mididx) = NaN;
    maxidx = find_bandpeak(search_kgph, true, min_prominence, sg_ord, sg_frame, medfilt_cols); % first prominent peak outward of the midline; NaN = no band in that column

    % 4.3 Find outer side pvs wall
    % 4.3.1 Apply max peak
    edge_kgph = side_kgph;
    edge_kgph(row_idx_grid > maxidx) = NaN;
    edge_kgph(:, isnan(maxidx)) = NaN;   % no band : the whole chain leaves the column NaN
    % 4.3.2 Find first shadow region or minimum dip region
    [edge_confined, minidx] = find_banddip(edge_kgph, true, min_prominence, sg_ord, sg_frame, min_diprows, medfilt_cols);
    % 4.3.3 Find half max from max to first shadow region
    [edge_thresholded, side_edgerawidx] = halfmax_confinedkymograph(edge_confined, true);
    side_edgeidx = clean_thresholdedkymograph(edge_thresholded, true, ax);
    % 4.3.4 An edge that landed inside the vessel wall is a PVS thinner than the half-max resolves : it sits on the wall, thickness 0
    is_insidewall = side_edgeidx > side_wallidx;
    side_edgeidx(is_insidewall) = side_wallidx(is_insidewall);
    fprintf('   side %d : PVS edge inside the wall in %d of %d columns, set to the wall\n', ...
        side_k, nnz(is_insidewall), numel(side_edgeidx));

    % 4.4 Find center side pvs wall, between the peak and the midline
    cent_kgph = search_kgph;
    cent_kgph(row_idx_grid < maxidx) = NaN;
    cent_kgph(:, isnan(maxidx)) = NaN;   % no band peak : without this the cap silently vanishes and the cent edge escapes past the shade
    [cent_thresholded, side_centrawidx] = halfmax_confinedkymograph(cent_kgph, false);
    side_centidx = clean_thresholdedkymograph(cent_thresholded, false, ax);

    %% debug : the side's scope, as the upper side sees it
    imagesc(ax, side_kgph)
    hold(ax, 'on')
    plot(ax, side_mididx, 'g')
    plot(ax, maxidx, 'r')
    plot(ax, minidx, 'm')
    plot(ax, side_edgeidx, 'w')
    plot(ax, side_centidx, 'c')
    hold(ax, 'off')

    % 4.5 Back to the kymograph's own orientation
    if side_k == 2
        side_rawrows = [side_edgerawidx; side_centrawidx];
        [~, side_rawrows] = flip_kymograph(side_kgph, side_rawrows);   % half-max rows
        side_edgerawidx = side_rawrows(1, :);
        side_centrawidx = side_rawrows(2, :);
        side_edgeidx = n_row + 2 - side_edgeidx;   % cleaned boundaries sit between rows
        side_centidx = n_row + 2 - side_centidx;
    end
    edge_idx(side_k, :) = side_edgeidx;
    cent_idx(side_k, :) = side_centidx;
    edge_rawidx(side_k, :) = side_edgerawidx;
    cent_rawidx(side_k, :) = side_centrawidx;
end

%% output
idx.pvsupedge_idx = edge_idx(1, :);
idx.pvsupcent_idx = cent_idx(1, :);
idx.pvsdownedge_idx = edge_idx(2, :);
idx.pvsdowncent_idx = cent_idx(2, :);

idx.pvsupedge_rawidx   = edge_rawidx(1, :);
idx.pvsupcent_rawidx   = cent_rawidx(1, :);
idx.pvsdownedge_rawidx = edge_rawidx(2, :);
idx.pvsdowncent_rawidx = cent_rawidx(2, :);

% 5. The masks : one row per edge, and the cropped kymograph between the edge and the vessel wall
uplocarray2d = repmat(idx.pvsupedge_idx, [sz(1), 1]); % 8.2
upline = false(sz);
upline(row_idx_grid == uplocarray2d) = 1;

lowlocarray2d = repmat(idx.pvsdownedge_idx, [sz(1), 1]); % 8.2
downline = false(sz);
downline(row_idx_grid == lowlocarray2d) = 1;

kymographmask.pvs_upline = upline;
kymographmask.pvs_downline = downline;

is_upbetween = (row_idx_grid - idx.pvsupedge_idx) .* (row_idx_grid - bv_upidx) <= 0;   % NaN edge : nothing
kymographmask.pvs_up = cropkymograph; % upper csf kymograph
kymographmask.pvs_up(~is_upbetween) = NaN;
is_downbetween = (row_idx_grid - idx.pvsdownedge_idx) .* (row_idx_grid - bv_downidx) <= 0;
kymographmask.pvs_down = cropkymograph; % lower csf kymograph
kymographmask.pvs_down(~is_downbetween) = NaN;

%% the four boundaries over the cropped kymograph
imagesc(ax, cropkymograph)
hold(ax, 'on')
plot(ax, idx.pvsupedge_idx, 'r')
plot(ax, idx.pvsupcent_idx, 'g')
plot(ax, idx.pvsdownedge_idx, 'r')
plot(ax, idx.pvsdowncent_idx, 'g')
hold(ax, 'off')

end

function [kgph_flipped, rows_flipped] = flip_kymograph(kgph, rows)
%FLIP_KYMOGRAPH  The kymograph upside down, and row indices on it carried onto the flipped rows.
%   Applied twice it gives back what it was given. This is the rule for a ROW; a boundary between
%   rows, as clean_thresholdedkymograph returns, flips as n_row + 2 - b instead.
% IN   kgph          H x W float  the kymograph
%      rows          K x W float  row indices on kgph, one line per index
% OUT  kgph_flipped  H x W float
%      rows_flipped  K x W float  n_row + 1 - rows
kgph_flipped = flipud(kgph);
rows_flipped = size(kgph, 1) + 1 - rows;
end

function maxidx = find_bandpeak(cropkymograph, upmod, min_prominence, sg_ord, sg_frame, medfilt_cols)
%FIND_BANDPEAK  Per column: the first prominent peak outward of the vessel midline, NaN if none.
%   The search starts at the midline, not the wall: the dark lumen is the baseline, so a band
%   flush to the wall is still an interior maximum. The half-max region stays wall-outward.

% IN   cropkymograph   H x W float  one side's kymograph, NaN inward of the midline
%      upmod           1x1 bool     true = upper side, outward is upward
%      min_prominence  1x1 float    least peak prominence, on the 0-1 rescaled kymograph
%      sg_ord          1x1 count    Savitzky-Golay order
%      sg_frame        1x1 count    Savitzky-Golay frame (odd); a column with fewer valid rows is skipped
%      medfilt_cols    1x1 count    median window across columns
% OUT  maxidx          1 x W float  band-peak row per column; NaN = no band
sz = size(cropkymograph);
maxidx = nan(1, sz(2));
for c = 1:sz(2)
    valid_rows = find(~isnan(cropkymograph(:, c)));
    if numel(valid_rows) < sg_frame
        continue
    end
    if upmod
        outward_rows = (valid_rows(end):-1:valid_rows(1))';   % wall -> outward (upward)
    else
        outward_rows = (valid_rows(1):valid_rows(end))';      % wall -> outward (downward)
    end
    profile = cropkymograph(outward_rows, c);
    if any(isnan(profile))
        profile = fillmissing(profile, 'linear');
    end
    profile = sgolayfilt(profile, sg_ord, sg_frame);
    [~, peak_locs] = findpeaks(profile, 'MinPeakProminence', min_prominence);
    if ~isempty(peak_locs)
        maxidx(c) = outward_rows(peak_locs(1));
    end
end
maxidx = round(medfilt1(maxidx, medfilt_cols, [], 2, 'omitnan', 'truncate'));   % bridges flicker up to half the window
end

function [crop_mask,minidx] = find_banddip(maxcropkymograph, upmod, min_prominence, sg_ord, sg_frame, min_diprows, medfilt_cols)
%FIND_BANDDIP  Per column: the first prominent dip outward of the band peak with nothing clearly
%   darker beyond it, and the kymograph cropped to [peak, dip]. A dip the profile later undercuts
%   by min_prominence is a wiggle inside the band, not the shade; judging the dip against the
%   column's own floor keeps a shallow shade, where an absolute fraction of the peak lost whole
%   sessions. With no prominent dip the deepest point stands in. A column whose profile never
%   falls min_prominence below the band peak has no shade and stays NaN.
%   was findshadow, then an absolute shade_frac gate; see CLAUDE_LOG.md
%
% IN   maxcropkymograph  H x W float  one side's kymograph, NaN inward of the band peak
%      upmod             1x1 bool     true = upper side, outward is upward
%      min_prominence    1x1 float    least dip prominence and least shade drop, on the 0-1 rescaled kymograph
%      sg_ord            1x1 count    Savitzky-Golay order
%      sg_frame          1x1 count    Savitzky-Golay frame (odd), shrunk to fit a shorter column
%      min_diprows       1x1 count    least valid rows for a dip search
%      medfilt_cols      1x1 count    median window across columns
% OUT  crop_mask         H x W float  the kymograph, NaN beyond the dip; all NaN where minidx is NaN
%      minidx            1 x W float  dip row per column; NaN where the column holds nothing or no shade
sz = size(maxcropkymograph);
[row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 2D row index grid
minidx = nan(1, sz(2));
is_noshade = false(1, sz(2));   % 1 x W bool  profile never falls min_prominence below the band peak
for c = 1:sz(2)
    valid_rows = find(~isnan(maxcropkymograph(:, c)));
    if numel(valid_rows) < min_diprows
        continue
    end
    if upmod
        outward_rows = (valid_rows(end):-1:valid_rows(1))';   % peak -> outward (upward)
    else
        outward_rows = (valid_rows(1):valid_rows(end))';      % peak -> outward (downward)
    end
    profile = maxcropkymograph(outward_rows, c);
    if any(isnan(profile))
        profile = fillmissing(profile, 'linear');
    end
    frame = min(sg_frame, numel(profile));
    if mod(frame, 2) == 0
        frame = frame - 1;
    end
    if frame > sg_ord
        profile = sgolayfilt(profile, sg_ord, frame);
    end
    shade_drop = profile(1) - min(profile);
    if shade_drop < min_prominence
        is_noshade(c) = true;     % no contrast beyond the band : nothing to threshold on
        continue
    end
    [~, dip_locs] = findpeaks(-profile, 'MinPeakProminence', min_prominence);
    dip_loc = [];
    for dip_candidate = dip_locs(:)'
        beyond_min = min(profile(dip_candidate:end));
        if beyond_min < profile(dip_candidate) - min_prominence
            continue              % a wiggle : something clearly darker lies beyond
        end
        dip_loc = dip_candidate;
        break
    end
    if isempty(dip_loc)
        [~, dip_loc] = min(profile);      % plain decay : the deepest point is the shade
    end
    if dip_loc < 2
        continue
    end
    minidx(c) = outward_rows(dip_loc);
end
minidx = round(medfilt1(minidx, medfilt_cols, [], 2, 'omitnan', 'truncate'));   % bridges flicker up to half the window
minidx(is_noshade) = NaN;   % after the filter : a frame with no shade is not bridged from its neighbours
crop_mask = maxcropkymograph;
if upmod
    crop_mask(row_idx_grid < minidx) = NaN;
else
    crop_mask(row_idx_grid > minidx) = NaN;
end
crop_mask(:, isnan(minidx)) = NaN;   % no dip : the column stays NaN instead of an unclosed window
end

function [thresholded_kgph,thresholded_rawidx] = halfmax_confinedkymograph(confined_kymograph,upmod)
sz = size(confined_kymograph);
[row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 2D row index grid
normconfined_kgph = confined_kymograph -min(confined_kymograph,[],1);
zerocolumn = sum(normconfined_kgph,1,'omitmissing');
zerocolumn = zerocolumn ==0 ;
normconfined_kgph(:,zerocolumn) = confined_kymograph(:,zerocolumn);
normconfined_kgph = normconfined_kgph./max(normconfined_kgph,[],1);
%%
zerocolumn = sum(normconfined_kgph,1,'omitmissing');
zerocolumn = zerocolumn ==0 ;
%%

thresholded_kgph = normconfined_kgph>=0.5;


thresholded_rawidx = thresholded_kgph.*row_idx_grid;
if upmod
    thresholded_rawidx(thresholded_rawidx==0) = Inf;
    thresholded_rawidx = min(thresholded_rawidx,[],1);
else
    thresholded_rawidx = max(thresholded_rawidx,[],1);
end
is_nocrossing = ~any(thresholded_kgph, 1);   % 1 x W bool  nothing reaches half max : no row to report
thresholded_rawidx(is_nocrossing) = NaN;
%%

end

function [sz1d, sz1dind,sizematrix] = binary2columnsize(binary_matrix)
sz = size(binary_matrix);
transit = abs(diff(binary_matrix,1,1));
transit = [true(1,sz(2));transit];
%%
columnszind = cumsum(transit,1);
%%
sz1dind = columnszind + (0:sz(2)-1)*max(columnszind,[],'all'); % make each column to be unique 1d count and reconstruction
sz1dind = sz1dind(:);
sz1d = accumarray(sz1dind, 1); % count idx = size
%%
sizematrix = sz1d(sz1dind); % duplicate size back to correspoinding idx
sizematrix = reshape(sizematrix,sz(1), sz(2));
end

