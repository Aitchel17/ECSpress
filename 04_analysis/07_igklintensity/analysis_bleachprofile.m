%ANALYSIS_BLEACHPROFILE  Bleaching decay of intensity bands, one session.
%
%   Reads the imaging beam as the probe rather than as a nuisance. Two-photon
%   excitation is confined to the focal volume, so a fixed imaging plane digs a
%   bleached slab about one axial PSF thick through otherwise unbleached tissue.
%   What is immobile inside that slab burns off and does not come back; what
%   exchanges is resupplied from outside the slab. The two therefore separate in
%   TIME, not in brightness, and the decay curve is the separation.
%
%   Every comparison is between intensity bands INSIDE one region, so no ruler is
%   hung off a moving boundary and no volume fraction crosses a comparison. The
%   masks are fixed at the constricted state and never track: a fixed mask is
%   sometimes wrong per frame, but it is wrong the same way throughout, which is
%   all a trend needs. See CLAUDE_LOG.md for what tracking masks did instead.

%% settings
param.box_px        = 10;           % int      sector_block box side, outside the PVS
param.window_s      = 600;          % float    only the early, steep part is read
param.band_pct      = [0 25; 30 60; 75 100];
                                    % 3x2 pct  low / middle / high intensity bands,
                                    %          as rank ranges over a region's pixels
param.band_name     = ["low", "mid", "high"];
param.decay_model   = "double";     % str      "single" or "double" exponential plus a floor
param.diam_tol_pct  = 5;            % float pct  frames are read only while the FWHM
                                    %          diameter sits within this of its median
param.min_block_px  = 40;           % int      blocks smaller than this are dropped
param.max_dist_um   = 30;           % float    outermost block layer, from the PVS wall
param.n_sectors     = 24;           % int      polar_frame sector count, matches clusterpolar
param.contour_bv    = "constrictedBV_contour";
param.contour_pvs   = "constrictedPVS_contour";

dirs.session = 'G:\tmp\00_igkl\hql090\251016_hql090_sleep\HQL090_sleep251016_005';
dirs.scratch = ['C:\Users\hql5715\AppData\Local\Temp\claude\g--03-program\' ...
    'aba78c46-71f9-4044-865e-82b968e72537\scratchpad'];

%% load, following analysis_main
% ECSSession inherits from mdfExtractLoader, which lives in the extractor tree, so
% both roots have to be on the path. util_ecspath adds exactly the three ECSpress
% roots and nothing above them -- other folders under the program root carry
% copies of these same file names, and the two trees already share one,
% pre_thresholding, whose two versions differ.
addpath('g:\03_program\01_ecspress\functions');   % where util_ecspath lives
util_ecspath;                                     % three roots, minus zz_notinuse
session = ECSSession(dirs.session);
session = session.load_primary_results;
session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
twophoton_processed = twophoton_preprocess(session);


%% which channel carries which label
% recordings may swap ch1/ch2, so the map comes off pax_fwhm when it is there
name_lumen = 'ch1';
name_pvs = 'ch2';
if isfield(session.pax_fwhm.param, 'channel_lumen')
    name_lumen = char(session.pax_fwhm.param.channel_lumen);
end
if isfield(session.pax_fwhm.param, 'channel_pvs')
    name_pvs = char(session.pax_fwhm.param.channel_pvs);
end
fprintf('lumen channel %s, pvs channel %s\n', name_lumen, name_pvs);

fs = twophoton_processed.fps;
pixel2um = twophoton_processed.pixel2um;
n_keep = min(size(twophoton_processed.(name_pvs), 3), round(param.window_s * fs));
stack_pvs = single(twophoton_processed.(name_pvs)(:,:,1:n_keep));
stack_lumen = single(twophoton_processed.(name_lumen)(:,:,1:n_keep));
clear twophoton_processed
t_axis = (0:n_keep-1)' / fs;
[n_y, n_x, n_frame] = size(stack_pvs);
fprintf('%d x %d x %d frames kept, %.1f s at %.3f Hz, %.4f um/px\n', ...
    n_y, n_x, n_frame, t_axis(end), fs, pixel2um);

%% the three regions, from the constricted polar contours
polarcluster = session.polarcluster;
roilist = session.roilist;
mask_cbv = roilist.getmask(char(param.contour_bv));
props = regionprops(mask_cbv, 'Centroid');
center = props(1).Centroid;
frame = polar_frame([n_y n_x], center, 'n_sectors', param.n_sectors);

% polar_theta is stored in RADIANS while angle_map is in degrees; querying one
% with the other silently extrapolates, because makima extrapolates by default
theta_bin = rad2deg(polarcluster.polar_theta(:).');
r_bv_bin = pick_profile(polarcluster.polar_profiles, param.contour_bv);
r_pvs_bin = pick_profile(polarcluster.polar_profiles, param.contour_pvs);
r_bv_map = wrap_makima(theta_bin, r_bv_bin, frame.angle_map);
r_pvs_map = wrap_makima(theta_bin, r_pvs_bin, frame.angle_map);

% the crop carries black borders where registration had no data; a block landing
% there is all zeros and would read as a region that never bleaches
valid = polarcluster.constricted_medianimg > 0;

mask.lumen = frame.radius_map < r_bv_map & valid;
mask.pvs = frame.radius_map >= r_bv_map & frame.radius_map < r_pvs_map & valid;
core = frame.radius_map < r_pvs_map;
fprintf('measured angles: bv %d/360, pvs %d/360\n', ...
    sum(isfinite(r_bv_bin)), sum(isfinite(r_pvs_bin)));
fprintf('lumen %d px, pvs %d px, core %d px, valid %d px\n', ...
    nnz(mask.lumen), nnz(mask.pvs), nnz(core), nnz(valid));
assert(nnz(mask.lumen & mask.pvs) == 0, 'lumen and PVS masks overlap');
assert(all(wrap_makima(theta_bin, r_bv_bin, theta_bin) < ...
    wrap_makima(theta_bin, r_pvs_bin, theta_bin)), 'the fitted annulus inverts at some angle');

%% the outside, tiled into equal-area boxes growing away from the PVS wall
[sectormap, block_info] = sector_block(core, param.box_px);
sectormap(~valid) = 0;
n_block = block_info.n_labels;
block_layer = zeros(1, n_block);
for layer = 1:numel(block_info.n_per)
    from_label = block_info.offset(layer) + 1;
    to_label = block_info.offset(layer) + block_info.n_per(layer);
    block_layer(from_label:to_label) = layer;
end
% layers keep growing to the frame corners, well past any tissue the PVS could
% supply, so the sweep is cut at a stated distance rather than at the crop edge
max_layer = floor(param.max_dist_um / (param.box_px * pixel2um));
fprintf('%d blocks in %d layers, box side %.2f um, layers kept %d\n', ...
    n_block, numel(block_info.n_per), param.box_px * pixel2um, max_layer);

%% band traces
% the bands are RANK ranges recomputed every frame, not a pixel cohort: drift and
% axial excursion mean a given pixel is not the same tissue at both ends of the
% window, so only the population is followed
stack_flat = reshape(stack_pvs, n_y * n_x, n_frame);
stack_lumen_flat = reshape(stack_lumen, n_y * n_x, n_frame);
clear stack_pvs stack_lumen

info.session = string(dirs.session);
info.fs = fs;
info.pixel2um = pixel2um;
info.t_axis = t_axis;
info.center = center;
info.n_px = struct('lumen', nnz(mask.lumen), 'pvs', nnz(mask.pvs));

info.pvs = band_trace(stack_flat, find(mask.pvs), param.band_pct);
info.lumen = band_trace(stack_flat, find(mask.lumen), param.band_pct);
info.lumen_control = band_trace(stack_lumen_flat, find(mask.lumen), param.band_pct);

block_band = nan(n_block, size(param.band_pct,1), n_frame);
block_npx = zeros(1, n_block);
for b = 1:n_block
    if block_layer(b) > max_layer
        continue
    end
    idx = find(sectormap == b);
    block_npx(b) = numel(idx);
    if block_npx(b) < param.min_block_px
        continue
    end
    block_band(b,:,:) = band_trace(stack_flat, idx, param.band_pct);
end
info.block = block_band;
info.block_layer = block_layer;
info.block_npx = block_npx;
info.block_dist_um = (block_layer - 0.5) * param.box_px * pixel2um;
info.band_name = param.band_name;
info.band_pct = param.band_pct;

fprintf('%d blocks kept of %d (>= %d px)\n', ...
    sum(block_npx >= param.min_block_px), n_block, param.min_block_px);

%% frames read only while the vessel sits at its own median diameter
% The masks are fixed, so a dilating wall walks into the PVS mask and drags its
% darker ranks down. Fitting the whole trace and calling the residual "vessel" put
% that drag inside the fit; dropping the frames instead takes it out of the sample.
% The ranks are recomputed per frame, so gating afterwards is the same as gating
% first, and keeping the full trace leaves the discarded frames visible.
bv_px = session.pax_fwhm.thickness.bv(:).';
n_shared = min(numel(bv_px), n_frame);
bv_px = bv_px(1:n_shared);
median_bv_px = median(bv_px);
frame_keep = false(1, n_frame);
frame_keep(1:n_shared) = abs(bv_px / median_bv_px - 1) <= param.diam_tol_pct / 100;
info.bv_px = bv_px;
info.median_bv_px = median_bv_px;
info.frame_keep = frame_keep;
fprintf('diameter median %.1f px = %.2f um; %d of %d frames within %.0f%%\n', ...
    median_bv_px, median_bv_px * pixel2um, sum(frame_keep), n_frame, param.diam_tol_pct);
fprintf('   whole-record median %.1f px, so the window is %+.1f%% off it\n', ...
    median(session.pax_fwhm.thickness.bv), ...
    100 * (median_bv_px / median(session.pax_fwhm.thickness.bv) - 1));
assert(sum(frame_keep) > 30, 'too few frames survive the diameter gate to fit anything');

%% the monotone part
% Fitted on the gated frames only. Every amplitude and time constant goes through
% exp(), so all of them are positive and the model is a sum of decaying
% exponentials on a positive floor -- monotone decreasing by construction.
% The amplitudes are the point: a fast one is what burns off and does not come
% back, a slow one is what is being resupplied, and the floor is what the beam
% never reaches inside this window.
n_band = size(param.band_pct, 1);
t_fit = t_axis(frame_keep);
[info.pvs_fit, info.pvs_coef] = fit_bands(t_fit, info.pvs(:,frame_keep), param.decay_model);
[info.lumen_fit, info.lumen_coef] = fit_bands(t_fit, info.lumen(:,frame_keep), param.decay_model);
[info.lumen_control_fit, info.lumen_control_coef] = ...
    fit_bands(t_fit, info.lumen_control(:,frame_keep), param.decay_model);
info.t_fit = t_fit;

fprintf('\n%-16s %-6s %10s %10s %10s %10s %10s %8s\n', ...
    'region', 'band', 'A_fast', 'tau_fast', 'A_slow', 'tau_slow', 'floor', 'R2');
report_fit('PVS', param.band_name, info.pvs(:,frame_keep), info.pvs_fit, info.pvs_coef);
report_fit('lumen', param.band_name, info.lumen(:,frame_keep), info.lumen_fit, info.lumen_coef);
report_fit('lumen ctrl', param.band_name, info.lumen_control(:,frame_keep), ...
    info.lumen_control_fit, info.lumen_control_coef);

%% figure
param.block_show    = [300 700];    % 1x2 int  label range averaged in the last panel

clee = color_lee;
band_color = [clee.clist.blue; clee.clist.darkgreen; clee.clist.red];
close(findall(groot, 'Type', 'figure', 'Name', 'bleachprofile'));
fig_bleach = figure('Name', 'bleachprofile', 'Color', 'w', 'Position', [40 40 1500 900]);

ax = subplot(3, 3, 1);
imagesc(ax, polarcluster.constricted_medianimg);
colormap(ax, gray);
axis(ax, 'image');
hold(ax, 'on');
visboundaries(ax, mask.lumen, 'Color', clee.clist.magenta, 'LineWidth', 1, 'EnhanceVisibility', false);
visboundaries(ax, mask.pvs, 'Color', clee.clist.teal, 'LineWidth', 1, 'EnhanceVisibility', false);
title(ax, 'constricted median, lumen and PVS masks');

ax = subplot(3, 3, 2);
show_block = double(sectormap);
show_block(sectormap == 0) = NaN;
imagesc(ax, show_block, 'AlphaData', ~isnan(show_block));
axis(ax, 'image');
colormap(ax, parula);
colorbar(ax);
title(ax, sprintf('sector\\_block, %d boxes', n_block));

ax = subplot(3, 3, 3);
hold(ax, 'on');
plot(ax, theta_bin, r_bv_bin, '.', 'Color', clee.clist.magenta, 'MarkerSize', 4);
plot(ax, theta_bin, r_pvs_bin, '.', 'Color', clee.clist.teal, 'MarkerSize', 4);
plot(ax, theta_bin, wrap_makima(theta_bin, r_bv_bin, theta_bin), '-', 'Color', clee.clist.magenta);
plot(ax, theta_bin, wrap_makima(theta_bin, r_pvs_bin, theta_bin), '-', 'Color', clee.clist.teal);
xlabel(ax, 'angle (deg)');
ylabel(ax, 'r (px)');
title(ax, 'contours, dots measured and lines makima');
grid(ax, 'on');

ax = subplot(3, 3, 4);
plot_bands(ax, t_axis, info.pvs, band_color, param.band_name);
title(ax, 'PVS');
ylabel(ax, '/ first 30 s');

ax = subplot(3, 3, 5);
plot_bands(ax, t_axis, info.lumen, band_color, param.band_name);
title(ax, 'lumen, tracer channel');

ax = subplot(3, 3, 6);
plot_bands(ax, t_axis, info.lumen_control, band_color, param.band_name);
title(ax, 'lumen, other channel (PMT onset control)');

ax = subplot(3, 3, [7 8 9]);
show = false(1, n_block);
show(max(1, param.block_show(1)):min(n_block, param.block_show(2))) = true;
show = show & block_npx >= param.min_block_px & block_layer <= max_layer;
band_show = squeeze(mean(block_band(show,:,:), 1, 'omitnan'));
plot_bands(ax, t_axis, band_show, band_color, param.band_name);
layer_show = block_layer(show);
dist_show = info.block_dist_um(show);
title(ax, sprintf('blocks %d-%d averaged: %d kept, layers %d-%d, %.1f-%.1f um from the PVS wall', ...
    param.block_show(1), param.block_show(2), sum(show), ...
    min(layer_show), max(layer_show), min(dist_show), max(dist_show)));
ylabel(ax, '/ first 30 s');
info.band_show = band_show;
info.block_show_kept = find(show);

exportgraphics(fig_bleach, fullfile(dirs.scratch, 'bleachprofile.png'), 'Resolution', 130);

%% figure, the monotone part on its own
% PVS against the block average, both divided by their own fitted value at t = 0,
% so the axis is a fraction of each region's own starting level and the two sit on
% one scale whatever their absolute brightness. The raw trace is kept faint behind
% each fit: the gap between them is what the vessel did, and it is not small.
[block_fit, info.block_show_coef] = fit_bands(t_fit, band_show(:,frame_keep), param.decay_model);
info.block_show_fit = block_fit;
report_fit('block', param.band_name, band_show(:,frame_keep), block_fit, info.block_show_coef);

close(findall(groot, 'Type', 'figure', 'Name', 'bleachmonotone'));
fig_monotone = figure('Name', 'bleachmonotone', 'Color', 'w', 'Position', [80 80 1200 820]);

ax = subplot(3, 1, 1);
hold(ax, 'on');
plot(ax, t_axis, bv_px(1:n_shared) * pixel2um, '-', 'Color', clee.clist.gray, 'LineWidth', 0.6);
plot(ax, t_fit, bv_px(frame_keep(1:n_shared)) * pixel2um, '.', ...
    'Color', clee.clist.black, 'MarkerSize', 4);
yline(ax, median_bv_px * pixel2um, '-', 'Color', clee.clist.red);
yline(ax, median_bv_px * pixel2um * (1 + param.diam_tol_pct/100), '--', 'Color', clee.clist.red);
yline(ax, median_bv_px * pixel2um * (1 - param.diam_tol_pct/100), '--', 'Color', clee.clist.red);
ylabel(ax, 'diameter (um)');
grid(ax, 'on');
title(ax, sprintf('%d of %d frames within %.0f%% of the median diameter', ...
    sum(frame_keep), n_frame, param.diam_tol_pct));

ax = subplot(3, 1, [2 3]);
hold(ax, 'on');
handle_cmp = gobjects(1, 2 * n_band);
label_cmp = strings(1, 2 * n_band);
for b = 1:n_band
    y_raw = info.pvs(b,frame_keep) / info.pvs_fit(b,1);
    y_fit = info.pvs_fit(b,:) / info.pvs_fit(b,1);
    plot(ax, t_fit, y_raw, '.', 'Color', [band_color(b,:) 0.25], 'MarkerSize', 3);
    handle_cmp(b) = plot(ax, t_fit, y_fit, '-', 'Color', band_color(b,:), 'LineWidth', 2.2);
    label_cmp(b) = "PVS " + param.band_name(b);
end
for b = 1:n_band
    y_raw = band_show(b,frame_keep) / block_fit(b,1);
    y_fit = block_fit(b,:) / block_fit(b,1);
    plot(ax, t_fit, y_raw, '.', 'Color', [band_color(b,:) 0.25], 'MarkerSize', 3);
    handle_cmp(n_band + b) = plot(ax, t_fit, y_fit, '--', 'Color', band_color(b,:), 'LineWidth', 1.8);
    label_cmp(n_band + b) = "block " + param.band_name(b);
end
xlabel(ax, 's');
ylabel(ax, '/ fitted value at t = 0');
grid(ax, 'on');
legend(ax, handle_cmp, cellstr(label_cmp), 'Location', 'eastoutside', 'AutoUpdate', 'off');
title(ax, sprintf('%s fit on the gated frames, PVS solid and block dashed', param.decay_model));

exportgraphics(fig_monotone, fullfile(dirs.scratch, 'bleachmonotone.png'), 'Resolution', 130);
save(fullfile(dirs.scratch, 'bleachprofile.mat'), 'info', 'param', 'mask', 'sectormap', '-v7.3');
fprintf('saved to %s\n', dirs.scratch);

%% ---------------------------------------------------------------- helpers
function r = pick_profile(polar_profiles, roi_label)
    % the 1x360 radius at each angle bin, NaN where the edge detector found nothing
    hit = strcmp({polar_profiles.roi_label}, char(roi_label));
    r = polar_profiles(hit).binned_r(:).';
end

function value = wrap_makima(theta_bin, r_bin, query_deg)
    % makima across the measured angles, padded by one turn each way so the
    % interpolation closes on itself rather than extrapolating at 0 and 360
    ok = isfinite(r_bin);
    theta_pad = [theta_bin(ok) - 360, theta_bin(ok), theta_bin(ok) + 360];
    r_pad = [r_bin(ok), r_bin(ok), r_bin(ok)];
    value = interp1(theta_pad, r_pad, query_deg, 'makima');
end

function out = band_trace(stack_flat, pixel_idx, band_pct)
    % out  B x T float   mean of each rank band of the region's pixels, per frame
    values = double(stack_flat(pixel_idx, :));
    sorted = sort(values, 1);
    n_px = size(sorted, 1);
    n_band = size(band_pct, 1);
    out = nan(n_band, size(sorted, 2));
    for b = 1:n_band
        from_rank = max(1, floor(band_pct(b,1) / 100 * n_px) + 1);
        to_rank = min(n_px, ceil(band_pct(b,2) / 100 * n_px));
        out(b,:) = mean(sorted(from_rank:to_rank, :), 1);
    end
end

function [fitted, coef] = fit_bands(t_axis, band, model)
    % fitted  B x T float   the monotone part of each band
    % coef    B x 5 float   A_fast, tau_fast, A_slow, tau_slow, floor.
    %                       a single-exponential fit leaves the slow pair NaN
    n_band = size(band, 1);
    fitted = nan(size(band));
    coef = nan(n_band, 5);
    for b = 1:n_band
        [fitted(b,:), coef(b,:)] = fit_decay(t_axis(:), band(b,:).', model);
    end
end

function [y_hat, coef] = fit_decay(t, y, model)
    % every parameter is fitted through exp(), so amplitudes, time constants and
    % the floor are all positive and the model cannot rise
    level = median(y(1:min(numel(y), 30)));
    span = max(t);
    if model == "single"
        start_p = log([0.3*level, span/6, 0.7*level]);
        shape = @(p, tt) exp(p(1)) * exp(-tt/exp(p(2))) + exp(p(3));
    else
        start_p = log([0.2*level, span/20, 0.15*level, span/2, 0.65*level]);
        shape = @(p, tt) exp(p(1)) * exp(-tt/exp(p(2))) ...
            + exp(p(3)) * exp(-tt/exp(p(4))) + exp(p(5));
    end
    cost = @(p) sum((shape(p, t) - y).^2);
    opt = optimset('MaxFunEvals', 5e4, 'MaxIter', 5e4, 'Display', 'off', 'TolX', 1e-8);
    best_p = fminsearch(cost, start_p, opt);
    y_hat = shape(best_p, t).';
    raw = exp(best_p);
    coef = nan(1, 5);
    if model == "single"
        coef([1 2 5]) = raw;
    else
        % order the two exponentials so the faster one is always first
        if raw(2) <= raw(4)
            coef = [raw(1) raw(2) raw(3) raw(4) raw(5)];
        else
            coef = [raw(3) raw(4) raw(1) raw(2) raw(5)];
        end
    end
end

function report_fit(region_name, band_name, band, fitted, coef)
    for b = 1:size(band, 1)
        residual = band(b,:) - fitted(b,:);
        r2 = 1 - var(residual) / var(band(b,:));
        fprintf('%-16s %-6s %10.1f %10.1f %10.1f %10.1f %10.1f %8.3f\n', ...
            region_name, band_name(b), coef(b,1), coef(b,2), coef(b,3), ...
            coef(b,4), coef(b,5), r2);
    end
end

function plot_bands(ax, t_axis, band, band_color, band_name)
    hold(ax, 'on');
    early = t_axis < 30;
    handle = gobjects(1, size(band,1));
    for b = 1:size(band, 1)
        y = band(b,:) / median(band(b, early));
        handle(b) = plot(ax, t_axis, y, 'Color', band_color(b,:), 'LineWidth', 1.5);
    end
    xlabel(ax, 's');
    grid(ax, 'on');
    legend(ax, handle, cellstr(band_name), 'Location', 'southwest', 'AutoUpdate', 'off');
end
