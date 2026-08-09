%PIV_EVENT_CONTROL_FIGURES  One event set against its lag-matched controls.
%   Correlates two transitions and two controls, then draws every view of the same
%   four fields: displacement, cell-mean vectors, Delaunay divergence, polar
%   divergence, radial strain, and the blink GIFs. Every panel of a given kind
%   shares one colour limit and one arrow scale, which is the only reason the
%   panels can be compared at all.
%
%   Assumes the base workspace already holds what piv_integration_testbed cells
%   0-4 leave: session, tp, t_axis, events, det, img_state, sleep_integrate.
%
%   why      controls are matched on LAG, not picked for flatness : PIV error
%            accumulates with the number of frames between the endpoints, so a
%            25 s transition's reference is a 25 s control. A flat window is a
%            different control answering a different question -- see the blink
%            pair in section 7, which IS picked for flatness
%   see FINDINGS.md
%   err      the naive sqrt(n_triangle) understates the divergence error by about
%            sqrt(6) : Delaunay reuses every vertex in ~6 triangles, so the
%            independent count is the VECTOR count over 3

clc
addpath(genpath('g:\03_program\01_ecspress'));

param.rem_index   = 2;      % which rem bout to bracket
param.box_px      = 40;     % cell side for the mean-vector partition
param.arrow_dense = 5;      % showpiv scale for the full field
param.arrow_cell  = 25;     % showpiv scale for the cell means
param.ring_w      = 20;     % polar ring width, px -- the geometry the scoring picked
param.n_sect      = 6;      % polar wedge count
param.band        = 2:4;    % rings read; ring 1 is 108 vectors over 6 wedges, sign flips
param.min_n       = 3;      % vectors a cell needs before it gets an arrow
param.gif_navg    = 15;     % frames averaged per blink panel

dirs.gif_dir      = '';     % char   '' = skip the blink GIFs, else an existing folder

exclmask = ~session.roilist.getmask('piv_include') | ...
            session.roilist.getmask('manual_dilated_pvs');
coremask = session.roilist.getmask('manual_dilated_pvs');

%% 1. The two transitions bracketing the chosen REM bout
rem_bout = img_state.state_idx.rem(param.rem_index, :);
row_nr = find(strcmp({events.state}, 'nrem2rem') & ...
    [events.to] >= rem_bout(1) - 60 & [events.from] <= rem_bout(2), 1);
row_ra = find(strcmp({events.state}, 'rem2awake') & ...
    [events.from] >= rem_bout(1) & [events.from] <= rem_bout(2) + 200, 1);

%% 2. Lag-matched controls
is_control = ~ismember({events.pol}, {'dilation', 'constriction'});
control_row = find(is_control);
control_lag = [events(control_row).to] - [events(control_row).from];
nearest_lag = @(target) control_row(find(abs(control_lag - target) == ...
    min(abs(control_lag - target)), 1));
row_c1 = nearest_lag(events(row_nr).to - events(row_nr).from);
row_c2 = nearest_lag(events(row_ra).to - events(row_ra).from);
if row_c2 == row_c1                     % never correlate the same window twice
    distance = abs(control_lag - (events(row_ra).to - events(row_ra).from));
    distance(control_row == row_c1) = Inf;
    row_c2 = control_row(find(distance == min(distance), 1));
end

panel_row = [row_nr, row_ra, row_c1, row_c2];
panel_tag = ["nrem2rem", "rem2awake", "control-1", "control-2"];
fprintf('%-11s %-12s %6s %6s %7s %9s\n', 'panel','state','from','to','lag(f)','dD(um)');
for k = 1:4
    e = events(panel_row(k));
    fprintf('%-11s %-12s %6d %6d %7d %+9.2f\n', panel_tag(k), e.state, ...
        e.from, e.to, e.to - e.from, e.diameter_change);
end

%% 3. Correlate and gate. endpoint only : consecutive costs the same again
ensemble = cell(1, 4);
for k = 1:4
    e = events(panel_row(k));
    fprintf('\n--- %s ---\n', panel_tag(k));
    pe = analysis_pivensemble(tp.ch1, [e.from, e.to], tp.fps, tp.pixel2um, halfwin = 2);
    pe.param.pivlab_mask = exclmask;
    pe.param.verbose = false;
    pe.correlate("endpoint");
    pe.gate("endpoint");
    ensemble{k} = pe;
    fprintf('%s: %d vectors kept\n', panel_tag(k), nnz(~isnan(pe.endpoint.uv(:,:,1))));
end

%% 4. Displacement, full field
for k = 1:4
    P = showpiv(sprintf('%s displacement', panel_tag(k)), param.arrow_dense);
    P.update_figsize([7 5]);
    P.plot_background(ensemble{k}.endpoint.stack(:,:,1));
    P.plot_quiver(ensemble{k}.endpoint.uv);
    P.put_xaxistitle('x (px)');  P.put_yaxistitle('y (px)');
end

%% 5. Mean vector per box cell
% why  averaging before drawing is not decimation : a decimated field shows one
%      arbitrary member of each neighbourhood, this shows what it agreed on
cellmean = cell(1, 4);
fprintf('\n%-11s %8s %9s %10s %12s %12s %7s\n', 'panel', 'arrows', 'vectors', ...
    'per cell', 'mean|cell|', 'mean|dense|', 'ratio');
for k = 1:4
    [~, ~, box_info] = vfield_strain(ensemble{k}.endpoint.uv, method = 'tri', ...
        exclmask = exclmask, box_px = param.box_px);
    [uv_cell, cm_info] = vfield_cellmean(ensemble{k}.endpoint.uv, ...
        box_info.sectormap, min_n = param.min_n);
    cellmean{k} = uv_cell;
    fprintf('%-11s %8d %9d %10.1f %12.3f %12.3f %7.2f\n', panel_tag(k), ...
        cm_info.n_cell, cm_info.n_vector, cm_info.n_vector/max(cm_info.n_cell,1), ...
        mean(cm_info.mag), mean(cm_info.mag_dense), ...
        mean(cm_info.mag)/mean(cm_info.mag_dense));

    M = showpiv(sprintf('%s cell-mean uv box %d x%d', panel_tag(k), param.box_px, param.arrow_cell), ...
        param.arrow_cell);
    M.update_figsize([7 5]);
    M.plot_background(ensemble{k}.endpoint.stack(:,:,1));
    M.plot_quiver(uv_cell);
    M.put_xaxistitle('x (px)');  M.put_yaxistitle('y (px)');
end

%% 6. Divergence and radial strain, triangles then polar wedges
[part, geo] = sector_polar(coremask, param.ring_w, param.n_sect);
wedges  = part(:,:,1);                       % the fused label, for the tri route
ring_ch = part(:,:,2);
n_rings = max(ring_ch(isfinite(ring_ch)), [], 'all');
tri_info = cell(1,4);  polar_cell = cell(1,4);  polar_map = cell(1,4);
radial_cell = cell(1,4);  radial_info = cell(1,4);
abs_tri = [];  abs_polar = [];  abs_eps = [];
for k = 1:4
    [~, ~, ti] = vfield_strain(ensemble{k}.endpoint.uv, method = 'tri', ...
        exclmask = exclmask);
    tri_info{k} = ti;
    abs_tri = [abs_tri; abs(ti.div_tri(:))];                       %#ok<AGROW>

    [dmap, dcells] = vfield_strain(ensemble{k}.endpoint.uv, method = 'tri', ...
        sectormap = wedges, exclmask = exclmask);
    % err  the 'tri' branch sizes its output by max(sectormap(:)), which falls
    %      short of the grid count whenever the outermost wedge holds no pixels.
    %      Pad before reshaping or the ring index is wrong
    padded = nan(1, n_rings * param.n_sect);
    padded(1:numel(dcells)) = dcells;
    polar_cell{k} = reshape(padded, n_rings, param.n_sect);
    polar_map{k}  = dmap;
    abs_polar = [abs_polar; abs(dcells(~isnan(dcells)))];          %#ok<AGROW>

    [~, ecells, ei] = vfield_strain(ensemble{k}.endpoint.uv, method = 'radial_tri', ...
        coremask = coremask, ring_width = param.ring_w, n_sectors = param.n_sect, ...
        exclmask = exclmask);
    radial_cell{k} = ecells;  radial_info{k} = ei;
    abs_eps = [abs_eps; abs(ecells(:))];                           %#ok<AGROW>
end
lim_tri   = prctile(abs_tri, 95);
lim_polar = prctile(abs_polar, 95);
lim_eps   = prctile(abs_eps(~isnan(abs_eps)), 95);
fprintf('\nshared limits | triangle %.5f | polar %.5f | eps_rr %.5f\n', ...
    lim_tri, lim_polar, lim_eps);

fprintf('\n%-11s %11s %11s %11s %11s %8s\n', ...
    'panel', 'div tri', 'div polar24', 'eps ring1', 'eps rings2-4', 'sigma');
for k = [3 4 1 2]                          % controls first so events draw on top
    d  = tri_info{k}.div_tri(~isnan(tri_info{k}.div_tri));
    pc = polar_cell{k};   ec = radial_cell{k};
    band = param.band(param.band <= size(pc,1));
    band_v = reshape(pc(band,:), 1, []);  band_v = band_v(~isnan(band_v));
    fprintf('%-11s %+11.5f %+11.5f %+11.5f %+11.5f %8.1f\n', panel_tag(k), ...
        median(d), median(band_v), median(ec(1,:), 'omitnan'), ...
        median(reshape(ec(band,:),1,[]), 'omitnan'), ...
        abs(mean(band_v)) / (std(band_v)/sqrt(numel(band_v))));

    T = showpiv(sprintf('%s divergence per triangle', panel_tag(k)), param.arrow_dense);
    T.update_figsize([7 5]);
    T.plot_background(ensemble{k}.endpoint.stack(:,:,1));
    T.plot_tri(tri_info{k}.xy, tri_info{k}.tri, tri_info{k}.div_tri, lim_tri);
    T.put_xaxistitle('x (px)');  T.put_yaxistitle('y (px)');

    G = showpiv(sprintf('%s divergence polar %dpx x %d', panel_tag(k), param.ring_w, param.n_sect), ...
        param.arrow_dense);
    G.update_figsize([7 5]);
    G.plot_background(ensemble{k}.endpoint.stack(:,:,1));
    G.plot_overlay(polar_map{k}, lim_polar);
    G.plot_quiver(cellmean{k});
    G.put_xaxistitle('x (px)');  G.put_yaxistitle('y (px)');

    E = showpiv(sprintf('%s eps_rr ring%d x %d', panel_tag(k), param.ring_w, param.n_sect), ...
        param.arrow_dense);
    E.update_figsize([7 5]);
    E.plot_background(ensemble{k}.endpoint.stack(:,:,1));
    E.plot_tri(radial_info{k}.xy, radial_info{k}.tri, radial_info{k}.eps_rr_tri, lim_eps);
    E.plot_quiver(ensemble{k}.endpoint.uv);
    E.put_xaxistitle('x (px)');  E.put_yaxistitle('y (px)');
end

%% 6.1 The wall every polar measure is referenced to
% caution  bwboundaries returns [row col]; plot takes x first, so the columns are
%          swapped. The axes y already runs downward because plot_background
%          laid the frame down with image()
wall = bwboundaries(coremask);
for f = findobj(groot, 'Type', 'figure')'
    if ~any(contains(string(f.Name), ["divergence", "eps_rr", "cell-mean"])); continue; end
    ax = findobj(f, 'Type', 'axes');
    if isempty(ax); continue; end
    delete(findobj(ax(1), 'Tag', 'pvs_wall'));      % do not stack on a re-run
    for b = 1:numel(wall)
        plot(ax(1), wall{b}(:,2), wall{b}(:,1), 'w-', 'LineWidth', 1.8, ...
            'Tag', 'pvs_wall', 'HandleVisibility', 'off');
        plot(ax(1), wall{b}(:,2), wall{b}(:,1), 'k-', 'LineWidth', 0.7, ...
            'Tag', 'pvs_wall', 'HandleVisibility', 'off');
    end
end

%% 7. Blink GIFs
% why  a DIFFERENT control from section 2 : this pair is picked for flatness and
%      caliber separation, because a blink shows a diameter difference and a
%      lag-matched control is not flat, it is net-zero
% err  `dirs.gif_dir == ""` compares a char array element by element and returns a
%      vector, so a non-empty path made the if fire on its first character
if isempty(dirs.gif_dir)
    return
end
if ~isfolder(dirs.gif_dir)
    error('piv_event_control_figures:gifdir', ...
        'dirs.gif_dir does not exist: %s', dirs.gif_dir);
end
n_frame = size(tp.ch1, 3);
frameset = @(centre) min(max(1, centre - floor(param.gif_navg/2)), n_frame - param.gif_navg + 1) ...
    + (0:param.gif_navg-1);

% 7.1 the transition's own two endpoints
e = events(row_nr);
piv_blink_gif(tp.ch1, tp.ch2, {frameset(e.from), frameset(e.to)}, ["", ""], ...
    fullfile(dirs.gif_dir, 'blink_nrem2rem.gif'), ...
    pixel2um = tp.pixel2um, scalebar_um = 10);

% 7.2 the flattest two controls that differ most in caliber
% err  gating on flatness first and duration second returns nothing : the flattest
%      windows ARE the shortest. Duration is structural, so it goes first
param.min_dur_s = 10;
ctl_stats = struct('row', {}, 'state', {}, 'from', {}, 'to', {}, ...
    'dur_s', {}, 'mean_um', {}, 'sd', {});
for k = 1:numel(control_row)
    ce = events(control_row(k));
    segment = det.dsg(ce.from:ce.to);
    ctl_stats(k) = struct('row', control_row(k), 'state', string(ce.state), ...
        'from', ce.from, 'to', ce.to, 'dur_s', t_axis(ce.to) - t_axis(ce.from), ...
        'mean_um', mean(segment), 'sd', std(segment));
end
long_enough = find([ctl_stats.dur_s] >= param.min_dur_s);
flat = long_enough([ctl_stats(long_enough).sd] <= ...
    prctile([ctl_stats(long_enough).sd], 50));
best = struct('separation', -Inf, 'a', NaN, 'b', NaN);
for i = 1:numel(flat)
    for j = i+1:numel(flat)
        gap = abs(ctl_stats(flat(i)).mean_um - ctl_stats(flat(j)).mean_um);
        if gap > best.separation
            best = struct('separation', gap, 'a', flat(i), 'b', flat(j));
        end
    end
end
if ctl_stats(best.a).mean_um > ctl_stats(best.b).mean_um   % smaller caliber first
    swap = best.a;  best.a = best.b;  best.b = swap;
end
A = ctl_stats(best.a);  B = ctl_stats(best.b);
fprintf('\nblink pair | %s %.2f um (sd %.3f) vs %s %.2f um (sd %.3f) | %.2f um apart\n', ...
    A.state, A.mean_um, A.sd, B.state, B.mean_um, B.sd, best.separation);
piv_blink_gif(tp.ch1, tp.ch2, ...
    {frameset(round((A.from+A.to)/2)), frameset(round((B.from+B.to)/2))}, ["", ""], ...
    fullfile(dirs.gif_dir, 'blink_control_pair.gif'), ...
    pixel2um = tp.pixel2um, scalebar_um = 10);
