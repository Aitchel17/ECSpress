function [up_csf_boundary, down_csf_boundary] = analyze_csfoutter(csf_kymograph, bv_upidx, bv_downidx, threshold)
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

sz = size(csf_kymograph);
[row_idx, ~] = ndgrid(1:sz(1), 1:sz(2));

% Top process (CSF external boundary using vessel upper boundary)
upcsf = row_idx < bv_upidx;
c_up = csf_kymograph;
c_up(~upcsf) = NaN;

% Outer boundary detection (upper)
up_csfoffset = prctile(c_up, 10, 1);
mask_up = c_up < up_csfoffset;
[~, upoffsetloc] = max(mask_up .* row_idx, [], 1);

% Upper FWHM calculation (CSF internal boundary using vessel center)

bv_centeridx = ceil((bv_downidx+bv_upidx)/2);


upcsf = row_idx < bv_centeridx;
c_up = csf_kymograph;
c_up(~upcsf) = NaN;
c_up(row_idx < upoffsetloc) = NaN;

c_up = c_up - up_csfoffset;
c_up = c_up ./ max(c_up, [], 1);

thr_cup = c_up > threshold;
csf_thridx = row_idx .* thr_cup;
csf_thridx(csf_thridx == 0) = Inf;
up_csf_boundary = min(csf_thridx, [], 1);

% Bottom process (CSF external boundary using vessel bottom boundary)
botcsf = row_idx > bv_downidx;
c_bot = csf_kymograph;
c_bot(~botcsf) = NaN;

% Outer boundary detection (bottom)
bot_csfoffset = prctile(c_bot, 10, 1);
mask_bot = c_bot < bot_csfoffset;
mask_bot_idx = mask_bot .* row_idx;
mask_bot_idx(mask_bot_idx == 0) = Inf;
[~, botoffsetloc] = min(mask_bot_idx, [], 1);

% Lower FWHM calculation (CSF internal boundary using vessel center)
downcsf = row_idx > bv_centeridx;
c_down = csf_kymograph;
c_down(~downcsf) = NaN;
c_down(row_idx > botoffsetloc) = NaN;

c_down = c_down - bot_csfoffset;
c_down = c_down ./ max(c_down, [], 1);

thr_cdown = c_down > threshold;
csf_thrdownidx = thr_cdown .* row_idx;
down_csf_boundary = max(csf_thrdownidx, [], 1);

%%
%
thresholded(row_idx == down_csf_boundary) = 1;
downlocarray2d = repmat(bottomboundary_idx, [sz(1), 1]); % 8.4
mask(row_idx > downlocarray2d) =0; % 8.5
thresholded(row_idx == up_csf_boundary) = 1;



end