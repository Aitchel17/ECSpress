function [idx, kymographmask] = get_invertedpvsoutter(kymograph, bv_upidx, bv_downidx)
% calculate_csf_boundaries - calculates upper and lower CSF boundaries based on given inputs.
%
% Inputs:
%   kymograph          - normalized CSF intensity profile (2D array)
%   bv_upidx  - upper boundary indices from blood vessel
%   bv_downidx - bottom boundary indices from blood vessel
%
% Outputs:
%   idx.up_csf_boundary    - upper CSF boundary indices
%   idx.down_csf_boundary  - lower CSF boundary indices
%   kymographmask.cb_thresholded     - updated binary mask with CSF boundaries

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

%% Restrict scope

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



%% 2.Separation process (CSF external boundary using vessel upper boundary)
% 2.1 Upper vessel kymograph generation (above upper vessel boundary)
bv_mididx = round((bv_upidx+bv_downidx)/2);

% upper vessel edge side kymograph
edgebv_upkymograph = cropkymograph;
edgebv_upkymograph(row_idx_grid > bv_upidx) = NaN; % upper vessel edge kymograph
% Debugging plot
% test_Fig = figure();
% test_ax = axes(test_Fig);
% figure(test_Fig)
% cla(test_ax)
% imagesc(test_ax,downedge_confined)
% hold on
% plot(test_ax,downedge_minidx,'r')
% plot(test_ax,downedge_rawidx,'g')

% findmax latched onto the bright vessel-wall bleed (contamination is severe here).
% Two explicit steps instead: (1) find the dark-PVS dip + the tissue max beyond it,
% (2) rise-onset edge where the signal first departs the dark floor on the outer
% half [dip -> max].
[up_dipidx, up_maxidx]      = find_dip_max(edgebv_upkymograph, true);
[upedge_idx, upedge_rawidx] = rise_onset_edge(edgebv_upkymograph, up_dipidx, up_maxidx, true);

%% Debugging
figure(fig)
imagesc(ax,cropkymograph)
hold on
plot(ax,upedge_idx,'g')

%% lower vessel edge side kymograph
edgebv_downkymograph = cropkymograph; % lower vessel edge kymograph
edgebv_downkymograph(row_idx_grid < bv_downidx) = NaN;
% Same two explicit steps, downward.
[dn_dipidx, dn_maxidx]          = find_dip_max(edgebv_downkymograph, false);
[downedge_idx, downedge_rawidx] = rise_onset_edge(edgebv_downkymograph, dn_dipidx, dn_maxidx, false);
%%
figure(fig)
cla(ax)
imagesc(ax,cropkymograph)
hold on
plot(ax,upedge_idx,'r')
plot(ax,downedge_idx,'r')

%% output
idx.pvsupedge_idx = upedge_idx;
idx.pvsdownedge_idx = downedge_idx;

idx.pvsupedge_rawidx   = upedge_rawidx;
idx.pvsdownedge_rawidx = downedge_rawidx;

%%
uplocarray2d = repmat(idx.pvsupedge_idx, [sz(1), 1]); % 8.2
upline = false(sz);
upline(row_idx_grid == uplocarray2d) = 1;

lowlocarray2d = repmat(idx.pvsdownedge_idx, [sz(1), 1]); % 8.2
downline = false(sz);
downline(row_idx_grid == lowlocarray2d) = 1;
%%
kymographmask.pvs_upline = upline;
kymographmask.pvs_downline = downline;
%%
imagesc(ax,kymographmask.pvs_upline+kymographmask.pvs_downline)
%%
kymographmask.pvs_up = edgebv_upkymograph; % upper csf kymograph
kymographmask.pvs_up(row_idx_grid<upedge_idx) = NaN;
kymographmask.pvs_down = edgebv_downkymograph; % lower csf kymograph
kymographmask.pvs_down(row_idx_grid>downedge_rawidx) = NaN; % lower csf kymograph

imagesc(ax,sum(cat(3,kymographmask.pvs_up,kymographmask.pvs_down),3,'omitnan'));

end

function [edge_idx, edge_raw] = rise_onset_edge(edgebv, dip_idx, max_idx, upmod)
%RISE_ONSET_EDGE  PVS edge = where the signal first departs the dark floor.
%   Given the per-column dark-PVS dip and tissue-max anchors (dip_idx, max_idx from
%   find_dip_max, which the caller runs explicitly), walk the OUTER half [dip -> max]
%   outward and take the first crossing of:
%       floor = min(outer-half),  peak = max(outer-half)
%       edge  = first outward crossing of  floor + frac*(peak - floor).
%   The dark PVS is normalized down to a flat 0 floor, so a half-max / Gaussian sits
%   mid-rise and spreads too far outward; a low frac catches where the signal LEAVES
%   the dark floor (the true outer PVS boundary) instead. floor/peak are taken from
%   the raw outer half (no offset subtraction). frac tunes how far out the edge sits
%   (0.5 = half-max, lower = nearer the rise onset). Cross-column median smoothed.
    [H, W]   = size(edgebv);
    edge_raw = nan(1, W);
    frac = 0.2;                                              % departure level fraction (tune)
    for c = 1:W
        if isnan(dip_idx(c)) || isnan(max_idx(c))
            continue; 
        end
        d = round(dip_idx(c));  
        m = round(max_idx(c));
        if d == m 
            continue; 
        end
        % dip -> max, outward
        if upmod
            rr = (d:-1:m).';
        else
            rr = (d:m).';
        end

        rr = rr(rr >= 1 & rr <= H);
        y  = edgebv(rr, c);

        if any(isnan(y))
            y = fillmissing(y, 'linear');
        end

        if numel(y) < 4
            continue
        end
        fr = min(5, numel(y));
        if mod(fr, 2) == 0
            fr = fr - 1;
        end
        if fr >= 3
            y = sgolayfilt(y, 2, fr);
        end
        flo = min(y);
        pk  = max(y);
        if pk - flo < 1e-6
            continue
        end
        L = flo + frac * (pk - flo);
        k = find(y >= L, 1, 'first');                        % first departure from the floor
        if isempty(k)
            continue
        end
        if k < 2
            edge_raw(c) = rr(1);
        else
            edge_raw(c) = rr(k-1) + (L - y(k-1)) / (y(k) - y(k-1)) * (rr(k) - rr(k-1));
        end
    end
    edge_raw = medfilt1(edge_raw, 11, [], 2, 'omitnan', 'truncate');
    edge_idx = round(edge_raw);
end

function [dip_idx, max_idx] = find_dip_max(edgebv, upmod)
%FIND_DIP_MAX  Per column: approximate dark-PVS dip (min) + the max just beyond it.
%   Marches the sgolay-smoothed profile from the vessel wall outward. The dip is the
%   smoothed minimum (the dark PVS valley -- robust to the bright wall bleed that
%   fooled findmax); the max is the largest value outward of the dip (the tissue
%   beyond the PVS). Both are approximate anchors; halfmax refines the actual edge.
%   Returns row indices (1 x W), median-smoothed across columns; NaN where undefined.
    [~, W]  = size(edgebv);
    dip_idx = nan(1, W);  
    max_idx = nan(1, W);
    for c = 1:W
        rr = find(~isnan(edgebv(:, c)));
        % 1. Don't restrict if it goes below 7 pixels
        if numel(rr) < 7 
            continue; 
        end
        % 2. flip direction if its upside outside decrease idx
        if upmod 
            rows = (rr(end):-1:rr(1)).'; 
        else
            rows = (rr(1):rr(end)).';
        end  % wall -> outward
        % 3. 
        p = edgebv(rows, c);
        if any(isnan(p))
            p = fillmissing(p, 'linear');
        end
        fr = min(9, numel(p));  
        if mod(fr, 2) == 0
            fr = fr - 1; 
        end
        if fr < 3
            continue
        end
        ps = sgolayfilt(p, 2, fr);
        [~, dloc] = min(ps);                 % dark PVS dip (approx valley centre)
        if dloc >= numel(ps)
            continue 
        end
        [~, mrel] = max(ps(dloc:end));       % max just beyond the dip (outward tissue)
        mloc = dloc + mrel - 1;
        if mloc <= dloc
            continue
        end
        dip_idx(c) = rows(dloc);
        max_idx(c) = rows(mloc);
    end
    dip_idx = medfilt1(dip_idx, 11, [], 2, 'omitnan', 'truncate');
    max_idx = medfilt1(max_idx, 11, [], 2, 'omitnan', 'truncate');
end

