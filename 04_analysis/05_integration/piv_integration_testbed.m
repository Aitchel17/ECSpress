%PIV_INTEGRATION_TESTBED 
%% 0. Load session, stack and sleep state

% 0.1 Paths : the repo roots, never genpath(pwd) -- a second copy of this tree
%     once sat in ECSpress_py, see CLAUDE.md
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
close all
% 0.2 Session directory
% sessiondir = 'G:\tmp\01_igkltdt\hql104\260607_hql104_sleep\HQL104_sleep260607_006';
sessiondir ='E:\02_egfptdt\hql099\260601_hql99_sleep\HQL99_sleep260601_005';

%% 0.3 Load session and both stacks
session = ECSSession(sessiondir);
session = session.load_primary_results;
session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
% 0.4 FPS matching + preprocessing
twophoton_processed = twophoton_preprocess(session);

% 0.6 Load sleep state and put the bouts on the imaging frame axis
score = load(fullfile(sessiondir, 'peripheral', 'sleep_score.mat'));
sleep_integrate = state_integration(score);
img_state = state_image(sleep_integrate);
img_state.get_state_indices(twophoton_processed.t_axis, twophoton_processed.outfps);


%% 1. Event detection on the FWHM diameter
% analysis_event owns this stage. 

% The partition is analysis_event.param -- six transitions and the three states
event_det = analysis_event(session.pax_fwhm.thickness.bv, session.pax_fwhm.t_axis);
event_det.arousaltransition(img_state.state_idx);

% The controls
event_det.control(img_state.state_idx);

% Both traces on the object are PX; the figure layer converts
disp(struct2table(rmfield(event_det.eventlist, {'search_from','search_to'})));

%% 1.1 Review the picked events
% The gate before PIV: a mis-picked endpoint shows up here. The figure layer is
% where px becomes um, so pixel2um is passed here and nowhere upstream
event_review(event_det.eventlist, event_det.diameter, event_det.smooth_diameter, ...
    event_det.info.P.sg_win_s, twophoton_processed.t_axis, ...
    'pixel2um',    twophoton_processed.pixel2um, ...
    'sleep_score', sleep_integrate.sleep_score, ...
    'state_idx',   img_state.state_idx, ...
    'states',      event_det.param.states);
%% 2. Pick the event
% events is chronological, so ask by label + id rather than by index
ev = find(strcmp({event_det.eventlist.state}, 'nrem2rem'), 1);  % first nrem2rem = id 2
fprintf('ev%d  %s  %s  dD %+.2f um  rise %.1f s  frames %d -> %d\n', ...
    ev, event_det.eventlist(ev).state, event_det.eventlist(ev).pol, ...
    event_det.eventlist(ev).diameter_change * twophoton_processed.pixel2um, ...
    event_det.eventlist(ev).rise_s, event_det.eventlist(ev).from, event_det.eventlist(ev).to);

%% 3. Make the ensemble object
% The object takes the frames it will correlate and nothing else. piv_getframescope
% cuts them, and refuses an event too close to either end of the recording
halfwin = 2;                                          % frames either side of from / to
[stack_span, span] = piv_getframescope(twophoton_processed.ch1, ...
    event_det.eventlist(ev).from, event_det.eventlist(ev).to, halfwin);
piv_ensemble = analysis_pivensemble(stack_span, halfwin, ...
    twophoton_processed.fps, twophoton_processed.pixel2um);

% 3.1 PIV masks

if ~has_roi(session.roilist, 'piv_include')
    session.roilist.addormodifyroi(twophoton_processed.ch1, 'piv_include', ...
        'rectangle', twophoton_processed.ch2);
    session.roilist.save2disk(session.dir_struct.primary_analysis);
end
if ~has_roi(session.roilist, 'manual_dilated_pvs')
    session.roilist.addormodifyroi( ...
        twophoton_processed.(session.pax_fwhm.param.channel_pvs), ...
        'manual_dilated_pvs', 'polygon', ...
        twophoton_processed.(session.pax_fwhm.param.channel_lumen));
    session.roilist.save2disk(session.dir_struct.primary_analysis);
end

coremask = session.roilist.getmask('manual_dilated_pvs');   % H x W bool, the wall
piv_ensemble.param.pivlab_mask = ~session.roilist.getmask('piv_include') | coremask;

% 4. Correlate. The expensive half; endpoint gives the displacement, consecutive
%    says whether it accumulated. Settings are param, set here
piv_ensemble.param.pivlab_windows     = [30 15; 20 10];
piv_ensemble.param.corr_floor_medfrac = [];           % corr_floor off
piv_ensemble.correlate("endpoint");
piv_ensemble.gate("endpoint");

%% 5. Confirmation of the displacement field
% vectors over the FROM frame, at showpiv's fixed scale
P = showpiv(sprintf('ev%d displacement', ev), 5);
P.update_figsize([7 5]);
P.plot_background(piv_ensemble.endpoint.stack(:,:,1));
P.plot_quiver(piv_ensemble.endpoint.uv);
P.put_xaxistitle('x (px)');  P.put_yaxistitle('y (px)');

%% 5.1 Blink comparator
% endpoint.stack is already interleaved (from_1, to_1, from_2, ...), so stepping
% one slice at a time flips between the two states
figure('Name', sprintf('ev%d blink  odd = FROM, even = TO', ev));
sliceViewer(piv_ensemble.endpoint.stack);

%% 7. One correlation plane
% get_corrplane correlates the result's stack again and returns the final-pass
% planes; the object keeps none. Centre cell = zero displacement, the offset is the
% final-pass residual. For consecutive, t picks the interval:
%   piv_ensemble.plane_at("consecutive", maps, 170, 86, 1:5)    the first five summed
maps = piv_ensemble.get_corrplane("endpoint");
[plane, info] = piv_ensemble.plane_at("endpoint", maps, 170, 86);
figure('Name', sprintf('%s plane at (%.0f,%.0f)', info.mode, info.x, info.y));
imagesc(plane - min(plane(:)));   axis image;   colorbar
title(sprintf('%s, %d pair(s) | uv (%+.2f, %+.2f) px | valid %d', ...
    info.mode, numel(info.pairs), info.uv(1), info.uv(2), info.valid));

%% 8.1 What the gates were reading, one square per interrogation window
% lambda2 and peak/rms are the two numbers analysis_pivensemble judges a window by,

for q = 1:numel(piv_ensemble.endpoint.gates)
    g = piv_ensemble.endpoint.gates(q);
    if ~g.info.on
        continue
    end
    switch g.name
        case 'tomasi'
            val = min(g.info.lambda2, [], 3);   % the weaker end; the gate reads that
            axis_name = '\lambda_2, weaker end';
        case 'corr_floor'
            val = g.info.peak_floor;
            axis_name = 'peak / rms';
        otherwise
            continue                            % neighbour keeps no per-window number
    end

    G = showpiv(sprintf('ev%d %s', ev, g.name), 5);
    G.update_figsize([7 5]);
    G.plot_background(piv_ensemble.endpoint.stack(:,:,1));
    G.plot_grid(piv_ensemble.endpoint.planes.xtable, ...
        piv_ensemble.endpoint.planes.ytable, val, [0 2*g.info.cut]);
    G.put_xaxistitle('x (px)');  G.put_yaxistitle('y (px)');
    title(G.ax, sprintf('%s   cut %.4g = %.2f x median   %d of %d rejected', ...
        axis_name, g.info.cut, g.info.cut / g.info.median, ...
        nnz(g.mask), numel(g.mask)));
end


%% 8.2 The mesh the divergence theorem will run on
% One object for both this cell and 8.3, so the field is triangulated ONCE.

grid_xyuv = piv_ensemble.endpoint.xyuv;   % gated, on the interrogation grid
vp = vfield_polar(grid_xyuv, coremask, twophoton_processed.pixel2um, ...
    exclmask = piv_ensemble.param.pivlab_mask, gate_name = piv_ensemble.endpoint.gate_name);
vp.applydelaunay();   % every cell Delaunay makes, none judged
vp.trifilter();       % long edges, slivers, cells on the mask -- a verdict
vp.placetri();        % each triangle's radius, wedge, bin and wall normal
vp.measure();         % the gradient, resolved against that normal

T = showpiv(sprintf('ev%d divergence per triangle', ev), 5);
T.update_figsize([7 5]);
T.plot_background(piv_ensemble.endpoint.stack(:,:,1));
% tri.conn, not mesh.conn : tri.* is the kept subset, mesh.conn every cell Delaunay made
T.plot_tri(vp.mesh.xy, vp.tri.conn, vp.tri.divergence, ...
    prctile(abs(vp.tri.divergence), 95));
T.put_xaxistitle('x (px)');  T.put_yaxistitle('y (px)');


%% 8.3 Every other quantity those same triangles carry, one figure each
frame_from = piv_ensemble.endpoint.stack(:,:,1);
clee       = color_lee;
disp_mag   = hypot(vp.tri.disp_radial, vp.tri.disp_tangential);
mag_lim    = prctile(disp_mag, 98);

% name | one value per triangle | colour limit ([] = symmetric on the data) | label
tri_panel = { ...
    'disp_radial',     vp.tri.disp_radial,     [],          'disp radial  (+ outward, \mum)'; ...
    'disp_tangential', vp.tri.disp_tangential, [],          'disp tangential  (+ CCW, \mum)'; ...
    'disp_magnitude',  disp_mag,               [0 mag_lim], 'disp magnitude  (\mum)'; ...
    'strain_radial',   vp.tri.strain_radial,   [],          'strain radial  n''En'; ...
    'strain_hoop',     vp.tri.strain_hoop,     [],          'strain hoop  t''Et'};

for q = 1:size(tri_panel, 1)
    name  = tri_panel{q, 1};
    val   = tri_panel{q, 2};
    lim   = tri_panel{q, 3};
    label = tri_panel{q, 4};
    fig_name = sprintf('ev%d %s tri', ev, name);

    Q = showpiv(fig_name, 5);
    Q.update_figsize([7 5]);
    Q.alpha = 0.9;                      % triangles tile with no gaps; 0.55 washes out
    if isempty(lim)
        lim = prctile(abs(val), 98);
    else
        % a magnitude has no negative half; sequential map
        Q.cmap = clee.gradient.inferno;
    end
    Q.plot_background(frame_from);
    Q.plot_tri(vp.mesh.xy, vp.tri.conn, val, lim);
    Q.plot_pivwall(coremask);        % radial and hoop are defined against it
    Q.put_xaxistitle('x (px)');  Q.put_yaxistitle('y (px)');
    title(Q.ax, label);
    % Q.save2svg(fig_dir);   one panel, one SVG; the destination is not this script's
end
fprintf('%d of %d triangles, %d panels | wedges %d, radial bins %d\n', ...
    size(vp.tri.conn, 1), size(vp.mesh.conn, 1), size(tri_panel, 1), ...
    vp.n_wedge, vp.n_bin);

%% 8.4 What the gates took
% red = uv_ungated minus uv, the gates' whole verdict
uv_removed = piv_ensemble.endpoint.uv_ungated;
uv_removed(~isnan(piv_ensemble.endpoint.uv)) = NaN;

V = showpiv(sprintf('ev%d gate verdict', ev), 5);
V.update_figsize([7 5]);
V.plot_background(piv_ensemble.endpoint.stack(:,:,1));
V.plot_quiver(uv_removed, 'red',  'red');
V.plot_quiver(piv_ensemble.endpoint.uv, 'cyan', 'blue');
V.put_xaxistitle('x (px)');  V.put_yaxistitle('y (px)');
fprintf('blue kept %d | red removed %d\n', ...
    nnz(~isnan(piv_ensemble.endpoint.uv(:,:,1))), nnz(~isnan(uv_removed(:,:,1))));

%% 9. Why the correlation plane cannot pick the vectors -- the evidence figure
% Writes piv_validation_corrgates.svg next to the source. Needs piv_ensemble
% correlated and gated, which cells 4 and 4.1 have done
piv_validation_testbed


% ======================================================================
%  BATCH PATH, cells 10-14. A DIFFERENT API from everything above.
%  caution  cell 10 is its cell 0. param.windows / exclmask / halfwin / piv_ch are
%           built here and nowhere else; cells 11-14 will not run without it
% ======================================================================

%% 10. BATCH PATH -- masks and the settings
% piv_include        : analyze ONLY inside this region (true = keep)
% manual_dilated_pvs : vessel + PVS at maximum dilation, traced by setup_rois. It is
%                      dropped from the PIV field and is the blob the polar rings
%                      are measured from.
piv_ch  = 'ch1';                                     % PIV channel -- set per dataset
halfwin = round(0.5 * twophoton_processed.outfps);   % frames either side of from/to
param.windows = [40 20; 20 10; 12 6];                % [window step] per pass, coarse to fine

% 10.1 piv_include and the traced wall are cell 3.1's; this path only reads them
if ~has_roi(session.roilist, 'piv_include') || ~has_roi(session.roilist, 'manual_dilated_pvs')
    error('piv_integration_testbed:noRois', ...
        'run cell 3.1 first -- it draws piv_include and manual_dilated_pvs.');
end

% 10.2 Combine into the two masks the rest of the pipeline uses
excl = session.roilist.getmask('manual_dilated_pvs');
try
    incl = session.roilist.getmask('piv_include');
catch; incl = [];
end
if isempty(incl) || ~any(incl(:)); incl = true(size(excl)); end
exclmask = ~incl | excl;                             % PIVlab convention: true = excluded

% 10.3 Check the traced boundary against the per-frame FWHM; eps is an outer-to-outer
%      diameter, so the PVS outer wall sits at eps/2
pvs_outer_r = double(session.pax_fwhm.thickness.eps(:).') / 2;      % px, per frame
mask_r      = sqrt(nnz(excl) / pi);
fprintf(['coremask r %.1f px | PVS outer r p50 %.1f p95 %.1f | ' ...
         'min margin %.1f px, %d frames uncovered\n'], ...
    mask_r, median(pvs_outer_r), prctile(pvs_outer_r, 95), ...
    min(mask_r - pvs_outer_r), nnz(pvs_outer_r > mask_r));

% 10.4 The traced boundary IS the ring origin; it is not grown. see FINDINGS.md
coremask = excl;

%% 11. Run PIV over every event
% One analysis_pivensemble per event, the object cell 3 makes; review cell 4 first.
% THE CONTROLS GO THROUGH PIV TOO : pick = {} is all events. Narrowing to the
% transitions is a tuning shortcut, and difference() has nothing to subtract then
pick    = {};
ev_list = 1:numel(event_det.eventlist);
if ~isempty(pick)
    ev_list = find(ismember({event_det.eventlist.state}, pick));
end
fprintf('PIV on %d of %d events | %d controls\n', numel(ev_list), ...
    numel(event_det.eventlist), ...
    nnz(strcmp({event_det.eventlist(ev_list).pol}, 'none')));

piv_run = piv_run_events(event_det.eventlist, twophoton_processed.(piv_ch), ...
    twophoton_processed.fps, twophoton_processed.pixel2um, ...
    halfwin = halfwin, sel = ev_list, windows = param.windows, exclmask = exclmask);

%% 12. Perivascular divergence, every event
% test_label narrows to one transition while tuning; '' = every event. A label
% drops the controls too, so leave it '' for any number that gets quoted
test_label = '';
sel = 1:numel(piv_run);
if ~isempty(test_label)
    sel = find(strcmp({piv_run.state}, test_label));
end
piv_sel = piv_run(sel);
fprintf('divergence on %d of %d events (%s) | %d controls\n', ...
    numel(sel), numel(piv_run), strjoin(unique({piv_sel.state}), ', '), ...
    nnz(strcmp({piv_sel.pol}, 'none')));

% every row carries gate_name; vfield_profile.append refuses to pool rows gated differently
[vprofile, vfields] = piv_polar_events(piv_sel, coremask, twophoton_processed.pixel2um, ...
    exclmask = exclmask, n_wedge = 12, bin_edges_um = 0:1.5:40, ...
    min_tri_wedge = 10, verbose = true);

%% 13. Batch figures -- shared decisions first
% the quiver scale and the colour limits are the same across every panel;
% computed here once and passed in
fig_style = struct('scale', 5, 'alpha', 0.55, 'cmap', 'turbo', ...
                   'head', 5, 'linew', 1, 'arrow', 'green');
lim_dr = vprofile.peak('disp_radial');     % one limit per quantity, over ALL rows
lim_dv = vprofile.peak('divergence');      % pass a per-event limit and each map
lim_vo = vprofile.peak('volume_out');      % self-normalises -- showpiv says so too

% the two per-event numbers the checks below read : one value per event over all its cells
mean_dr = arrayfun(@(r) mean(r.cells.disp_radial(:), 'omitnan'), vprofile.rows)';
mean_dv = arrayfun(@(r) mean(r.cells.divergence(:),  'omitnan'), vprofile.rows)';
pol_all = string({piv_sel.pol});
dD_all  = [piv_sel.diameter_change];

%% 13.1 Sign check -- the two polarities must land in opposite half-planes
% disp_radial is a CHECK (the wall has to follow dD); divergence is a RESULT and
% the count is reported, not asserted -- see DIVERGENCE_RESULT.md
co = color_lee.polarity(pol_all);
figure('Name', 'sign check', 'Position', [200 200 950 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact');
for q = {'disp_radial', mean_dr, 'wall motion follows \DeltaD'; ...
         'divergence',  mean_dv, 'divergence follows polarity'}'
    nexttile;
    hold on
    grid on
    scatter(gca, dD_all, q{2}, 70, co, 'filled');
    xline(gca, 0, 'k--');
    yline(gca, 0, 'k--');
    xlabel(gca, '\DeltaD (\mum)');
    ylabel(gca, showpiv.axis_label(q{1}));
    [n_ok, n_tot] = event_validatesign(pol_all, q{2});
    title(gca, sprintf('%s   (%d/%d consistent)', q{3}, n_ok, n_tot));
end

%% 13.2 Radial profile -- ONE wedge set for every row, fixed before the curves are built
quantity = 'volume_out';
allrows  = true(1, vprofile.n_row);
[bin_range, wedge_ok] = vprofile.widest_common(allrows, quantity);
if any(isnan(bin_range))
    warning('piv_integration_testbed:noCommonSupport', ...
        'no wedge spans any interval for every row; nothing to plot');
else
    curves = vprofile.curve_rows(allrows, bin_range, quantity = quantity);
    r      = vprofile.bin_center_um;
    figure('Name', 'radial profile', 'Position', [120 80 1000 460]);
    hold on
    grid on
    % the quiet rows are on the same axes : the pipeline's own floor
    for k = 1:vprofile.n_row
        plot(r, curves(k, :), '-o', 'Color', co(k, :), 'LineWidth', 1.1, ...
            'MarkerSize', 4, 'HandleVisibility', 'off');
    end
    yline(0, 'k--', 'HandleVisibility', 'off');
    xlabel('\mum from the vessel wall');
    ylabel(showpiv.axis_label(quantity));
    title(sprintf(['red = dilation, blue = constriction, grey = quiet   |   ' ...
        '%d/%d wedges over %.1f-%.1f \mum'], nnz(wedge_ok), vprofile.n_wedge, ...
        r(bin_range(1)), r(bin_range(2))));
end

%% 13.3 One event, one quantity, one figure
% Nothing tiles; the panels are arranged in Illustrator. The maps are PAINTED by
% paint_cells, not measured, and nothing downstream reads them
k = 1;
R = piv_sel(k);
bg   = mat2gray(twophoton_processed.(piv_ch)(:, :, R.from));   % the EARLIER state
keep = ~isnan(R.xyuv(:,:,3)) & ~isnan(R.xyuv(:,:,4));
uv   = piv_stamp(R.xyuv, R.planes.imsize, keep);        % drawing only
cap  = [sprintf('ev%d  %s bout %d  %s  %.1f s  ', R.id, R.state, R.bout, ...
        R.pol, R.rise_s) '\DeltaD=' sprintf('%+.2f', R.diameter_change) ' \mum'];

for q = {'disp_radial', lim_dr; 'divergence', lim_dv; 'volume_out', lim_vo}'
    which = q{1};
    lim   = q{2};
    P = showpiv(sprintf('%s ev%d %s %s', which, R.id, R.state, R.pol), fig_style.scale);
    P.alpha      = fig_style.alpha;
    P.linew      = fig_style.linew;
    P.head_len   = fig_style.head;
    P.head_width = fig_style.head;
    P.cmap       = fig_style.cmap;
    % one colour for shaft and head : there is one population here
    P.fill_color = fig_style.arrow;
    P.line_color = fig_style.arrow;
    % background, overlay, wall, arrows : in that order
    P.plot_background(bg);
    P.plot_overlay(vfields{k}.paint_cells(vprofile.rows(k).cells.(which)), lim);
    P.plot_pivwall(coremask);
    P.plot_quiver(uv);
    title(P.ax, [showpiv.axis_label(which) '   ' cap]);
end

%% 13.4 One event's diameter trace, zoomed
% the same event_review call as cell 1.1, then xlim
k = 1;
ax_tr = event_review(event_det.eventlist, event_det.diameter, ...
    event_det.smooth_diameter, event_det.info.P.sg_win_s, twophoton_processed.t_axis, ...
    'pixel2um',    twophoton_processed.pixel2um, ...
    'sleep_score', sleep_integrate.sleep_score, ...
    'state_idx',   img_state.state_idx, ...
    'states',      event_det.param.states);
w  = [event_det.eventlist(piv_sel(k).id).start_f, event_det.eventlist(piv_sel(k).id).end_f];
rg = max(1, w(1) - 60):min(numel(event_det.diameter), w(2) + 60);
xlim(ax_tr, twophoton_processed.t_axis([rg(1) rg(end)]));

%% 13.5 Raw frames PIV used (blink comparator)
% odd slice = FROM, even = TO. piv_interleave is the function the class ran on its
% span, so these are the frames PIV used, before the preprocessing
[inter, pair_frames] = piv_interleave(twophoton_processed.(piv_ch), ...
    piv_sel(1).from, piv_sel(1).to, halfwin);
figure('Name', sprintf('ev%d %s %s  PIV frames (odd=from, even=to)', ...
    piv_sel(1).id, piv_sel(1).state, piv_sel(1).pol));
sliceViewer(inter);
fprintf('pair_frames matches the stored run: %d\n', ...
    isequal(pair_frames, piv_sel(1).pair_frames));

%% 14. Assemble piv_events / piv_result  (not saved -- inspect only)
% 14.1 Polar frame from the pax axis and the vessel centroid
paxv      = session.roilist.getvertices('pax');
pax_angle = atan2d(paxv(2,2) - paxv(1,2), paxv(2,1) - paxv(1,1));
[vy, vx]  = find(coremask);  center = [mean(vx), mean(vy)];
azframe   = polar_frame(size(coremask), center, 'pax_angle', pax_angle, 'n_sectors', 8);

% 14.2 Detection record
piv_events = struct('source_ref', 'fwhm', ...
    'events', event_det.eventlist, 'focus', {event_det.focus}, ...
    'detect', struct('params', det.P, 'sign_ok', det.sign_ok, ...
                     'clipped', det.clipped, 'dropped', det.dropped));

% 14.3 PIV record
piv_result = struct();
piv_result.meta    = struct('pixel2um', twophoton_processed.pixel2um, 'outfps', twophoton_processed.outfps, 'halfwin', halfwin, ...
    'channel', piv_ch);
piv_result.frame   = struct('center', center, 'pax_angle', pax_angle, 'n_sectors', 8, ...
    'ring_edges', polar{1}.info.ring_edges);
% ring_width comes from whoever set it (cell 12); sector_polar does not echo its arguments
piv_result.byevent = piv_sel;
piv_result.polar   = polar;

fprintf('\npiv_sel: %d events. piv_events / piv_result assembled (not saved).\n', numel(piv_sel));

%% ---- local functions ----
function tf = has_roi(rlist, name)
%HAS_ROI  True if the roilist already carries this ROI (getmask errors if not).
    try
        m  = rlist.getmask(name);
        tf = ~isempty(m) && any(m(:));
    catch
        tf = false;
    end
end
