%PIV_INTEGRATION_TESTBED  Diameter-event-driven PIV, run cell by cell.
%   FWHM diameter -> merge bouts -> pick excursions -> review -> per-event
%   ensemble PIV -> gating -> strain. Cell 0 loads the session itself; run
%   analysis_main through FWHM and setup_rois first, because this reads
%   .pax_fwhm and the 'pax' ROI. 'piv_include' and 'manual_dilated_pvs' are drawn
%   here when the roilist has not got them; re-tracing is a command you type.
%
%   TWO PATHS, and they do not share an API. Know which one you are in.
%     cells 0-9   ONE event, through the analysis_pivensemble CLASS. Current.
%                 Everything written since 260806 targets this.
%     cells 10-14 EVERY event, through piv_run_events / piv_polar_events, and
%                 then plain functions. Older, function-based, still runs.
%   caution  cell 10 is the batch path's cell 0 : it builds S, param.windows, exclmask,
%            halfwin and piv_ch, which cells 11-14 need and cells 0-9 never make.
%            Running cell 11 straight after cell 9 fails on an undefined S
%   note     the class path re-preprocesses inside analysis_pivensemble, so S and
%            the class are two independent copies of the same frames. That is
%            wasteful, not wrong -- collapsing them is the obvious next cleanup
%
%   note  this file replaced the 260803 testbed on 260808. The settled window
%         ladder and coremask measurements that used to live in this header moved
%         to analysis_pivensemble, beside the settings they justify. The old
%         header's "OPEN: ring 1 v_r runs below rings 2-3" is now ANSWERED --
%         small n, 108 vectors over 6 wedges -- and is recorded in vfield_strain
%
%% 0. Load session, stack and sleep state

% 0.1 Paths
% The repo roots, not pwd. genpath(pwd) from g:\03_program pulled in a second copy
% of this whole tree that used to sit in ECSpress_py, and MATLAB resolves by
% filename, so which copy ran depended on path order. See CLAUDE.md.
addpath('g:\03_program\01_ecspress\functions');   % where util_ecspath lives
util_ecspath;                                     % three roots, minus zz_notinuse
close all
% 0.2 Session directory
% sessiondir = 'G:\tmp\01_igkltdt\hql104\260607_hql104_sleep\HQL104_sleep260607_006';
sessiondir ='E:\02_egfptdt\hql099\260601_hql99_sleep\HQL99_sleep260601_005';

% 0.3 Load session and both stacks
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
% analysis_event owns this stage. It is handed a trace and the bout tables and
% knows nothing about line_fwhm or about sleep scoring, so it can be re-run on a
% saved session without rebuilding either.
%   .bv is the ABSOLUTE caliber. .bvchanges is median-subtracted and would band an
%   amplitude with no caliber in it.
%   PX throughout, as every analysis class here stores; the figure layer converts.

% The partition is analysis_event.param -- six transitions and the three states the
% controls search, defaulted in the class. It is not repeated here because control()
% requires the SAME one arousaltransition was given, and a table written at two call
% sites cannot show whether the two agree. Override with event_det.param.merge for a
% different design; event_det.focus is derived from its column 3 and is never set
event_det = analysis_event(session.pax_fwhm.thickness.bv, session.pax_fwhm.t_axis);
event_det.arousaltransition(img_state.state_idx);

% The controls, into the same list. A row is (state, pol): a transition is
% ('nrem2rem','dilation'), a stable stretch is ('nrem','none') -- two axes, so a
% within-state event later is just ('rem','dilation') and needs no new category.
% Several lengths per bout because most of the PIV error accumulates with LAG,
% and range vs rise_s across these rows is how that gets measured. See FINDINGS.md
%   a long stable window is the FLATTEST available, not flat -- its range can be
%   wider than the median transition, which is what the per-row range field is for
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
% Whole recording plus the two endpoints; the margin, the cut and the filtering
% are the object's business. Only the frames it correlates survive
piv_ensemble = analysis_pivensemble(twophoton_processed.ch1, ...
    [event_det.eventlist(ev).from, event_det.eventlist(ev).to], ...
    twophoton_processed.fps, twophoton_processed.pixel2um, halfwin = 2);

% 3.1 PIV masks
% Analyse inside piv_include and outside the PVS, which is one mask in PIVlab's
% polarity. Both ROIs are drawn when the roilist has not got them, and neither is
% reopened on a re-run -- addormodifyroi opens an editor, and a cell that does that
% every time it runs cannot be re-run.
%   TO RE-TRACE the wall, paste this at the prompt and re-run the cell. It is
%   written out rather than using a variable, so it works whatever has been run:
%     session.roilist.modifyroi( ...
%         twophoton_processed.(session.pax_fwhm.param.channel_pvs), ...
%         'manual_dilated_pvs', ...
%         twophoton_processed.(session.pax_fwhm.param.channel_lumen))
%     session.roilist.save2disk(session.dir_struct.primary_analysis)
%   Deliberate, and it leaves nothing behind. A redraw flag would have to be set,
%   run, and then remembered back to false
%   caution  the channel is not ch1/ch2 by fiat -- the two can be swapped between
%            datasets, and pax_fwhm.param is the only thing that says which is
%            which. setup_rois draws this ROI on the PVS channel with the lumen
%            channel as the second view, and this matches that
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

% 4. Correlate
% The expensive half, tens of seconds. endpoint gives the displacement,
% consecutive says whether it accumulated; name one to skip the other's cost
piv_ensemble.correlate("endpoint", window_sizes = [30 15; 20 10]);

piv_ensemble.gate(corr_floor = false); % gate

%% 5. Confirmation of the displacement field
% Vectors over the FROM frame. showpiv fixes the scale instead of autoscaling,
% so one event can be held against another
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
% Planes are captured before any gating, so a window the gates rejected still
% has one. Centre cell = zero displacement; the offset is the final-pass residual.
% Naming t switches to consecutive and picks the frame interval:
%   piv_ensemble.plane_at(170, 86, 7)      one interval
%   piv_ensemble.plane_at(170, 86, 1:5)    the first five summed
[plane, info] = piv_ensemble.plane_at(170, 86);
figure('Name', sprintf('%s plane at (%.0f,%.0f)', info.which, info.x, info.y));
imagesc(plane - min(plane(:)));   axis image;   colorbar
title(sprintf('%s, %d pair(s) | uv (%+.2f, %+.2f) px | valid %d', ...
    info.which, numel(info.pairs), info.uv(1), info.uv(2), info.valid));

%% 8.1 What the gates were reading, one square per interrogation window
% lambda2 and peak/rms are the two numbers analysis_pivensemble judges a window by,
% and both are PER WINDOW. plot_grid is the layer for that: it draws the window at
% its own size, where plot_overlay would paint it as a pixel and plot_contour would
% interpolate between windows -- either one claiming a resolution never measured.
%   the colour axis runs 0..2*cut on a diverging map, so the THRESHOLD is the white
%   line: below it the gate rejects, above it the gate keeps. That is the question
%   this figure exists to answer, and a 98th-percentile axis would not answer it
%   caution  the squares are the STEP, not the interrogation area. Windows overlap
%            2:1 here, so neighbours touch while the data behind them share half
%            their pixels
%   note     a gate whose SWITCH is off still measures: info.on follows the
%            threshold in param, not the switch. So corr_floor is drawn even though
%            cell 4.1 leaves it out, and this is what it would have cut
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
% One object for both this cell and 8.3, so the field is triangulated ONCE. Two
% calls, because the two halves answer different questions: applydelaunay depends
% on WHERE the vectors are, measure on what they say.
%   Each triangle IS the cell whose gradient gets computed, so painting it flat
%   draws the partition the number came from and nothing invented. Slivers, and the
%   bridges Delaunay throws across the masked vessel, are dropped before the
%   gradient runs; info.dropped counts them
grid_xyuv = piv_ensemble.endpoint.xyuv;   % gated, on the interrogation grid
vp = vfield_polar(grid_xyuv, coremask, twophoton_processed.pixel2um, ...
    exclmask = piv_ensemble.param.pivlab_mask, gated = true);
vp.applydelaunay();   % every cell Delaunay makes, none judged
vp.trifilter();       % long edges, slivers, cells on the mask -- a verdict
vp.placetri();        % each triangle's radius, wedge, bin and wall normal
vp.measure();         % the gradient, resolved against that normal

T = showpiv(sprintf('ev%d divergence per triangle', ev), 5);
T.update_figsize([7 5]);
T.plot_background(piv_ensemble.endpoint.stack(:,:,1));
% tri.conn, not mesh.conn : measure() is the one place trifilter's verdict is
% applied, so every tri.* column is the KEPT subset while mesh.conn is every cell
% Delaunay made. Passing the two together is a size mismatch plot_tri refuses
T.plot_tri(vp.mesh.xy, vp.tri.conn, vp.tri.divergence, ...
    prctile(abs(vp.tri.divergence), 95));
T.put_xaxistitle('x (px)');  T.put_yaxistitle('y (px)');


%% 8.3 Every other quantity those same triangles carry, one figure each
% The same object, already measured in 8.2 -- vfield_polar put those triangles into
% polar coordinates about the wall, and that is what turns four gradient components
% into a radial and a hoop part. Five more numbers on one partition, not a second
% measurement. Divergence is not repeated here; 8.2 has it.
%
% One object is one panel, and nothing tiles: the panels are arranged outside
% MATLAB. Every signed panel takes the same colour rule, symmetric on the 98th
% percentile, because a panel that self-normalises cannot be held against another
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
        % a magnitude has no negative half, so a diverging map would spend half
        % its range on values that cannot occur
        Q.cmap = clee.gradient.inferno;
    end
    Q.plot_background(frame_from);
    Q.plot_tri(vp.mesh.xy, vp.tri.conn, val, lim);
    Q.plot_pivwall(coremask);        % radial and hoop are defined against it
    Q.put_xaxistitle('x (px)');  Q.put_yaxistitle('y (px)');
    title(Q.ax, label);
    % Q.save2svg(fig_dir);   one panel, one SVG. Left off because it writes five
    %                        files per event and the destination is not this
    %                        script's to pick
end
fprintf('%d of %d triangles, %d panels | wedges %d, radial bins %d\n', ...
    size(vp.tri.conn, 1), size(vp.mesh.conn, 1), size(tri_panel, 1), ...
    vp.n_wedge, vp.n_bin);

%% 8.4 What the gates took
% Two calls, two colours. uv_ungated is every vector after the common mode came
% out; uv is what survived, so the difference is the three gates' whole verdict
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
%  caution  cell 10 is its cell 0. S / param.windows / exclmask / halfwin / piv_ch are
%           built here and nowhere else; cells 11-14 will not run without it
% ======================================================================

%% 10. BATCH PATH -- masks and the preprocessed stack
% piv_include        : analyze ONLY inside this region (true = keep)
% manual_dilated_pvs : vessel + PVS at maximum dilation, traced by setup_rois. It is
%                      dropped from the PIV field and is the blob the polar rings
%                      are measured from.
piv_ch  = 'ch1';                                     % PIV channel -- set per dataset
halfwin = round(0.5 * twophoton_processed.outfps);   % frames either side of from/to
param.windows = [40 20; 20 10; 12 6];                      % [window step] per pass, coarse to fine
S       = piv_preprocess(twophoton_processed.(piv_ch));       % motion-preserving preprocess

% 5.1 piv_include and the traced wall are cell 3.1's; this path only reads them
if ~has_roi(session.roilist, 'piv_include') || ~has_roi(session.roilist, 'manual_dilated_pvs')
    error('piv_integration_testbed:noRois', ...
        'run cell 3.1 first -- it draws piv_include and manual_dilated_pvs.');
end

% 5.2 Combine into the two masks the rest of the pipeline uses
excl = session.roilist.getmask('manual_dilated_pvs');
try
    incl = session.roilist.getmask('piv_include');
catch; incl = [];
end
if isempty(incl) || ~any(incl(:)); incl = true(size(excl)); end
exclmask = ~incl | excl;                             % PIVlab convention: true = excluded

% 5.3 Check the traced boundary against the per-frame FWHM. eps = totalpvs + bv is
%     an outer-to-outer DIAMETER, so the PVS outer wall sits at eps/2
pvs_outer_r = double(session.pax_fwhm.thickness.eps(:).') / 2;      % px, per frame
mask_r      = sqrt(nnz(excl) / pi);
fprintf(['coremask r %.1f px | PVS outer r p50 %.1f p95 %.1f | ' ...
         'min margin %.1f px, %d frames uncovered\n'], ...
    mask_r, median(pvs_outer_r), prctile(pvs_outer_r, 95), ...
    min(mask_r - pvs_outer_r), nnz(pvs_outer_r > mask_r));

% 5.4 The traced boundary IS the ring origin: it was drawn on the maximum-dilation
%     frame, so growing it would push ring 1 into tissue the PVS never reaches.
%     Measured: a 6 px grow moved ring 1 v_r by 2% and cost 3 cells
coremask = excl;

%% 11. Run PIV over every event
% ~20 s per event, so review cell 4 first
% THE CONTROLS GO THROUGH PIV TOO. event_det.focus is the transition labels only,
% and running that alone leaves the profile with zero pol=="none" rows -- at which
% point difference() returns a table of NaN without erroring, and every published
% number from this pipeline was a difference against those rows. {} = all of them.
%   caution  it is ~2x the events and ~2x the PIV time. Narrowing to the
%            transitions is a tuning shortcut, not a setting to leave behind
pick    = {};
ev_list = 1:numel(event_det.eventlist);
if ~isempty(pick)
    ev_list = find(ismember({event_det.eventlist.state}, pick));
end
fprintf('PIV on %d of %d events | %d controls\n', numel(ev_list), ...
    numel(event_det.eventlist), ...
    nnz(strcmp({event_det.eventlist(ev_list).pol}, 'none')));

piv_run = piv_run_events(event_det.eventlist, S, 'halfwin', halfwin, 'sel', ev_list, ...
    'piv_opt', {'window_sizes', param.windows, 'repeat', 1, 'do_pad', 1, ...
                'exclmask', exclmask, 'use_gpu', true});

%% 12. Perivascular divergence, every event
% test_label narrows to one transition while tuning; '' = every event in piv_run
%   caution  a label here drops the CONTROLS along with everything else, and the
%            profile then has nothing to difference against. Leave it '' for any
%            number that gets quoted
test_label = '';
sel = 1:numel(piv_run);
if ~isempty(test_label)
    sel = find(strcmp({piv_run.state}, test_label));
end
piv_sel = piv_run(sel);
fprintf('divergence on %d of %d events (%s) | %d controls\n', ...
    numel(sel), numel(piv_run), strjoin(unique({piv_sel.state}), ', '), ...
    nnz(strcmp({piv_sel.pol}, 'none')));

% gated = anything was applied to the field after the correlation. TRUE here
% because piv_run_events calls piv_validate at ens_interleave:97 -- PIVlab's
% outlier removal, NOT the tomasi / corr_floor gates that analysis_pivensemble
% applies. The flag is coarse on purpose: what it protects is pooling rows that
% were filtered differently, and vfield_profile.append refuses a mixed set.
% Whether either kind of filtering belongs here is a separate question with a
% measured answer -- see CLAUDE_LOG.md, the gates remove the larger vectors
[vprofile, vfields] = piv_polar_events(piv_sel, coremask, twophoton_processed.pixel2um, ...
    exclmask = exclmask, n_wedge = 12, bin_edges_um = 0:1.5:40, ...
    min_tri_wedge = 10, gated = true, verbose = true);

%% 13. Batch figures -- shared decisions first
% Two things must be the same across every panel or the panels cannot be read
% against each other: the quiver scale and the colour limits. They are computed
% HERE, once, and passed in. There is no figure class holding them any more --
% a class that owns a shared setting hides it, and the setting is the part a
% reader has to see
fig_style = struct('scale', 5, 'alpha', 0.55, 'cmap', 'turbo', ...
                   'head', 5, 'linew', 1, 'arrow', 'green');
lim_dr = vprofile.peak('disp_radial');     % one limit per quantity, over ALL rows
lim_dv = vprofile.peak('divergence');      % pass a per-event limit and each map
lim_vo = vprofile.peak('volume_out');      % self-normalises -- showpiv says so too

% the two per-event numbers the checks below read. One value per event, every
% cell of that event one sample; the wedge-membership question does not arise
% here the way it does for a radial profile
mean_dr = arrayfun(@(r) mean(r.cells.disp_radial(:), 'omitnan'), vprofile.rows)';
mean_dv = arrayfun(@(r) mean(r.cells.divergence(:),  'omitnan'), vprofile.rows)';
pol_all = string({piv_sel.pol});
dD_all  = [piv_sel.diameter_change];

%% 13.1 Sign check -- the two polarities must land in opposite half-planes
% One quantity per axis. disp_radial is a CHECK (the wall has to follow dD or
% something is wrong upstream); divergence is a RESULT and the count is
% reported, not asserted -- see DIVERGENCE_RESULT.md.
% Drawn here rather than in a function: it is a scatter with two reference
% lines, and wrapping that would hide the one thing worth sharing, which is the
% colour convention -- that is color_lee.polarity
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

%% 13.2 Radial profile -- ONE wedge set for every row
% The wedge set is fixed BEFORE the curves are built. A mean over whatever
% wedges reach each ring averages more of them near the wall than far from it,
% and the slope that comes out is partly the change in membership
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
    % the quiet rows are on the SAME axes. They are the only thing that says
    % whether a curve is a measurement or the pipeline's own floor
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
% Nothing tiles. The deliverable is a representative event per transition, and
% the panels are placed next to each other in Illustrator -- a figure that bakes
% in an arrangement is a figure that has to be taken apart again.
% The maps are PAINTED, not measured: vfield_polar bins triangles and keeps no
% dense map; paint_cells spreads a per-cell value back over the frame so it can
% be looked at where it came from. Nothing downstream reads a painted map
k = 1;
R = piv_sel(k);
bg   = mat2gray(S(:, :, R.from));                       % the EARLIER state
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
    % one colour for shaft and head: showpiv's two-colour default marks a second
    % population, and there is only one here
    P.fill_color = fig_style.arrow;
    P.line_color = fig_style.arrow;
    % the wall goes down AFTER the overlay and BEFORE the arrows -- the overlay
    % is opaque enough to swallow it, and the arrows stay on top of both
    P.plot_background(bg);
    P.plot_overlay(vfields{k}.paint_cells(vprofile.rows(k).cells.(which)), lim);
    P.plot_pivwall(coremask);
    P.plot_quiver(uv);
    title(P.ax, [showpiv.axis_label(which) '   ' cap]);
end

%% 13.4 One event's diameter trace, zoomed
% event_review OWNS this picture -- raw trace, the smoothed one the picking saw,
% the state bands and the from/to markers. A second copy lived on the old figure
% class purely to zoom to one event, and a zoom is xlim, not a function.
% Same call as cell 1.1 so the two figures are the same picture
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
% Odd slice = FROM, even = TO, so stepping one slice at a time separates wall
% motion from a whole-field drift. piv_interleave is the SAME function the
% correlator ran, not a copy of it, so these really are the frames PIV used
[inter, pair_frames] = piv_interleave(S, piv_sel(1).from, piv_sel(1).to, halfwin);
figure('Name', sprintf('ev%d %s %s  PIV frames (odd=from, even=to)', ...
    piv_sel(1).id, piv_sel(1).state, piv_sel(1).pol));
sliceViewer(inter);
fprintf('pair_frames matches the stored run: %d\n', ...
    isequal(pair_frames, piv_sel(1).planes.pair_frames));

%% 14. Assemble piv_events / piv_result  (not saved -- inspect only)
% 9.1 Polar frame from the pax axis and the vessel centroid
paxv      = session.roilist.getvertices('pax');
pax_angle = atan2d(paxv(2,2) - paxv(1,2), paxv(2,1) - paxv(1,1));
[vy, vx]  = find(coremask);  center = [mean(vx), mean(vy)];
azframe   = polar_frame(size(coremask), center, 'pax_angle', pax_angle, 'n_sectors', 8);

% 9.2 Detection record
piv_events = struct('source_ref', 'fwhm', ...
    'events', event_det.eventlist, 'focus', {event_det.focus}, ...
    'detect', struct('params', det.P, 'sign_ok', det.sign_ok, ...
                     'clipped', det.clipped, 'dropped', det.dropped));

% 9.3 PIV record
piv_result = struct();
piv_result.meta    = struct('pixel2um', twophoton_processed.pixel2um, 'outfps', twophoton_processed.outfps, 'halfwin', halfwin, ...
    'channel', piv_ch);
piv_result.frame   = struct('center', center, 'pax_angle', pax_angle, 'n_sectors', 8, ...
    'ring_edges', polar{1}.info.ring_edges);
% note  ring_width used to be read back off info. sector_polar stopped echoing its
%       own arguments 260808, so it comes from whoever set it -- here, cell 12
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
