%PIV_INTEGRATION_TESTBED  Diameter-event-driven PIV, run cell by cell.
%   FWHM diameter -> merge bouts -> pick excursions -> review -> per-event
%   ensemble PIV -> gating -> strain. Cell 0 loads the session itself; run
%   analysis_main through FWHM and setup_rois first, because this reads
%   .pax_fwhm and the 'pax' / 'manual_dilated_pvs' ROIs. Only 'piv_include' is
%   drawn here.
%
%   TWO PATHS, and they do not share an API. Know which one you are in.
%     cells 0-9   ONE event, through the analysis_pivensemble CLASS. Current.
%                 Everything written since 260806 targets this.
%     cells 10-14 EVERY event, through piv_run_events / piv_polar_events /
%                 piv_figure. Older, function-based, still runs.
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
% The repo root, not pwd. genpath(pwd) from g:\03_program pulled in a second copy
% of this whole tree that used to sit in ECSpress_py, and MATLAB resolves by
% filename, so which copy ran depended on path order. See CLAUDE.md.
addpath(genpath('g:\03_program\01_ecspress'));
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

% 0.5 Shorthands used throughout
tp     = twophoton_processed;
fps    = tp.outfps;
p2u    = tp.pixel2um;
t_axis = tp.t_axis;

% 0.6 Load sleep state and put the bouts on the imaging frame axis
sleep_integrate = state_integration(sessiondir);
img_state = state_image(sleep_integrate);
img_state.get_state_indices(t_axis, fps);     

%% 1. Event detection on the FWHM diameter
% analysis_event owns this stage. It is handed a trace and the bout tables and
% knows nothing about line_fwhm or about sleep scoring, so it can be re-run on a
% saved session without rebuilding either.
%   .bv is the ABSOLUTE caliber. .bvchanges is median-subtracted and would band an
%   amplitude with no caliber in it.
%   PX throughout, as every analysis class here stores; the figure layer converts.

param.merge = { ...
%    earlier(A)    later(B)      label             polarity
    'nrem',       'rem',        'nrem2rem',       'dilation'     ; ...
    'uarousal',   'nrem',       'uarousal2nrem',  'dilation'     ; ...
    'roughawake', 'nrem',       'awake2nrem',     'dilation'     ; ...
    'rem',        'roughawake', 'rem2awake',      'constriction' ; ...
    'nrem',       'roughawake', 'nrem2awake',     'constriction' ; ...
    'nrem',       'uarousal',   'nrem2uarousal',  'constriction' };

param.focus  = param.merge(:, 3)';                         % labels to run PIV on
param.states = {'roughnrem', 'rem', 'roughawake'};   % the merge partition: the states a
                                               % bout is scored into, so also the
                                               % ones a stable window lives inside

caliber   = session.pax_fwhm.thickness.bv;
event_det = analysis_event(double(caliber(:).'), session.pax_fwhm.t_axis);
event_det.arousaltransition(img_state.state_idx, param.merge);

% The controls, into the same list. A row is (state, pol): a transition is
% ('nrem2rem','dilation'), a stable stretch is ('nrem','none') -- two axes, so a
% within-state event later is just ('rem','dilation') and needs no new category.
% Several lengths per bout because most of the PIV error accumulates with LAG,
% and range vs rise_s across these rows is how that gets measured. See FINDINGS.md
%   a long stable window is the FLATTEST available, not flat -- its range can be
%   wider than the median transition, which is what the per-row range field is for
event_det.control(img_state.state_idx, states = param.states);

% Shorthands the PIV cells below read. Both traces are PX
events = event_det.eventlist;
dtrace = event_det.diameter;
dsg    = event_det.smooth_diameter;   % picking trace, what every figure marks on
disp(struct2table(rmfield(events, {'search_from','search_to'})));

%% 1.1 Review the picked events
% The gate before PIV: a mis-picked endpoint shows up here. The figure layer is
% where px becomes um, so pixel2um is passed here and nowhere upstream
event_review(events, dtrace, dsg, event_det.info.P.sg_win_s, t_axis, ...
    'pixel2um',    p2u, ...
    'sleep_score', sleep_integrate.sleep_score, ...
    'state_idx',   img_state.state_idx, ...
    'states',      param.states);
%% 2. Pick the event
% events is chronological, so ask by label + id rather than by index
ev = find(strcmp({events.state}, 'nrem2rem'), 1);  % first nrem2rem = id 2
fprintf('ev%d  %s  %s  dD %+.2f um  rise %.1f s  frames %d -> %d\n', ...
    ev, events(ev).state, events(ev).pol, events(ev).diameter_change * p2u, ...
    events(ev).rise_s, events(ev).from, events(ev).to);

%% 3. Make the ensemble object
% Whole recording plus the two endpoints; the margin, the cut and the filtering
% are the object's business. Only the frames it correlates survive
piv_ensemble = analysis_pivensemble(twophoton_processed.ch1, ...
    [events(ev).from, events(ev).to], ...
    twophoton_processed.fps, twophoton_processed.pixel2um, halfwin = 2);

%% 3.1 PIV masks
% piv_include is drawn here when the roilist has not got it; the PVS boundary
% comes from setup_rois. Analyse inside piv_include and outside the PVS, which is
% one mask in PIVlab's polarity
if ~any(strcmp(session.roilist.list(), 'piv_include'))
    session.roilist.addormodifyroi(twophoton_processed.ch1, 'piv_include', ...
        'rectangle', twophoton_processed.ch2);
    session.roilist.save2disk(session.dir_struct.primary_analysis);
end

piv_ensemble.param.pivlab_mask = ~session.roilist.getmask('piv_include') | ...
                                  session.roilist.getmask('manual_dilated_pvs');

%% 4. Correlate
% The expensive half, tens of seconds. endpoint gives the displacement,
% consecutive says whether it accumulated; name one to skip the other's cost
piv_ensemble.correlate();                % or correlate("endpoint")

%% 4.1 Gate
% Cheap and repeatable, so this cell can be re-run on its own after changing
% param or switching a gate off. gate() always restarts from the correlation, so
% calls do not stack. It prints one line per gate; when a line looks wrong,
% piv_ensemble.endpoint.gates is a 1x3 table of them --
%   struct2table(rmfield(piv_ensemble.endpoint.gates, {'mask','info'}))
% and gates(k).info holds the thresholds that gate was given and what it measured
%   piv_ensemble.gate(tomasi = false)     lambda2 gate off, threshold kept
%   piv_ensemble.gate(neighbour = true)   put the local median back for one call
%   piv_ensemble.gate("endpoint")         gate one result only
%   piv_ensemble.param.verbose = false    silence
% note  neighbour runs OFF : stage NaN in the table, but gates(3).mask still holds
%       the verdict it would have carried out
piv_ensemble.gate();

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

%% 6. Strain on polar wedges
% eps_rr, eps_tt and their trace all come out of one call now; info carries the
% components whatever method is asked for. Near the PVS wall the trace is
% POSITIVE -- tangential stretch beats radial compression there, which is the
% regime Causemann et al. (PNAS 2026) attribute to a stiff PVS. Beyond ~30 px it
% flips to mild compression. Read the two bands separately, the pooled number
% averages the sign away
coremask = session.roilist.getmask('manual_dilated_pvs');
[~, ~, pol2] = vfield_divergence_polar(piv_ensemble.endpoint.uv, coremask, ...
    ring_width = 10, n_sectors = 12, min_valid = 1, method = 'divergence');
near = false(size(pol2.eps_rr));  near(1:3,:) = true;
kk   = ~isnan(pol2.eps_rr) & ~isnan(pol2.eps_tt);
for band = {near & kk, ~near & kk, kk}
    m = band{1};
    fprintf(['%-10s n %3d | eps_rr %+.5f | eps_tt %+.5f | div %+.5f | ' ...
             'eps_tt/|eps_rr| %5.2f | expanding %3.0f%%\n'], ...
        "", nnz(m), median(pol2.eps_rr(m)), median(pol2.eps_tt(m)), ...
        median(pol2.divergence(m)), ...
        median(pol2.eps_tt(m))/abs(median(pol2.eps_rr(m))), ...
        100*mean(pol2.divergence(m) > 0));
end

%% 6.1 Radial strain map
% Wedges exist to isolate the RADIAL term. Each pixel is projected onto the
% radial direction first, then averaged per wedge, then differenced between rings
% along the same wedge -- so tangential stretch never enters. That is the
% 'radial' method; 'plane' fits du/dx+dv/dy per wedge and mixes the two, which is
% what the square boxes are for. r = 0 is the PVS wall, not the centroid
coremask = session.roilist.getmask('manual_dilated_pvs');
% min_valid is a wedge's vector count. Four is the rank the PLANE fit needs; the
% radial method only averages a scalar, so it can go lower, but the CENTRED
% DIFFERENCE that follows cannot -- measured, 40% of its gradient cells lean on a
% wedge holding one or two vectors, and those cells came out at +0.002 against
% -0.008 for the rest. 'radial_fit' fits the vectors in a radial window instead
% and carries its own conditioning, so it does not need this knob
[dvrdr_map, dvrdr, pol] = vfield_divergence_polar(piv_ensemble.endpoint.uv, ...
    coremask, ring_width = 10, n_sectors = 12, min_valid = 5, method = 'radial');

fprintf('%d rings x %d wedges | %d wedges carried a v_r | %d carried a gradient\n', ...
    pol.n_rings, pol.n_sectors, nnz(~isnan(pol.vr_cells)), nnz(~isnan(dvrdr)));
fprintf('mean v_r %+.3f px | dv_r/dr p50 %+.5f | compressed (dv_r/dr<0) %.0f%%\n', ...
    mean(pol.vr_cells(:),'omitnan'), median(dvrdr(~isnan(dvrdr))), ...
    100*mean(dvrdr(~isnan(dvrdr)) < 0));

fprintf('\n ring  r from wall   v_r(px)   dv_r/dr    wedges\n');
for k = 1:pol.n_rings
    if all(isnan(pol.vr_cells(k,:))); continue; end
    fprintf('  %2d   %3.0f-%3.0f px    %+7.3f  %+9.5f   %2d\n', k, ...
        pol.ring_edges(k), pol.ring_edges(k+1), mean(pol.vr_cells(k,:),'omitnan'), ...
        mean(dvrdr(k,:),'omitnan'), nnz(~isnan(pol.vr_cells(k,:))));
end

D = showpiv(sprintf('ev%d dv_r/dr', ev), 5);
D.update_figsize([7 5]);
D.plot_background(piv_ensemble.endpoint.stack(:,:,1));
D.plot_overlay(dvrdr_map, prctile(abs(dvrdr(~isnan(dvrdr))), 95));
D.plot_quiver(piv_ensemble.endpoint.uv);
D.put_xaxistitle('x (px)');  D.put_yaxistitle('y (px)');

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

%% 8. Divergence three ways, from the same u,v map
% vfield_strain is a router: one input, one output shape, and only method=
% changes. fd differentiates point by point, tri integrates over Delaunay
% triangles (the divergence theorem, so interior edges cancel), the polar wedges
% isolate the RADIAL component. Agreement in sign is the cross-validation; the
% spread is how much vector noise each route lets through
% see FINDINGS.md
% err       'radial_fit' is eps_rr, ONE strain component, not the trace. Holding
%           it against fd compares two different quantities
coremask = session.roilist.getmask('manual_dilated_pvs');
for route = ["fd" "tri" "plane" "radial_fit"]
    [strain_map, strain_cells, strain_info] = vfield_strain( ...
        piv_ensemble.endpoint.uv, method = char(route), coremask = coremask, ...
        exclmask = piv_ensemble.param.pivlab_mask);
    strain_finite = strain_map(~isnan(strain_map));
    fprintf('%-11s cells %4d | p50 %+.5f | IQR %.5f | %s\n', route, ...
        numel(strain_cells), median(strain_finite), iqr(strain_finite), strain_info.unit);
end

%% 8.1 The window grid, coloured
% plot_grid draws one flat square PER INTERROGATION WINDOW, at the step size.
% plot_overlay would paint the same numbers per pixel and claim a resolution the
% measurement never had; plot_contour would interpolate between windows
[div_map, ~, div_info] = vfield_strain(piv_ensemble.endpoint.uv, method = 'fd', ...
    exclmask = piv_ensemble.param.pivlab_mask);
grid_planes = piv_ensemble.endpoint.planes;
grid_live   = ~isnan(piv_ensemble.endpoint.uv_grid(:,:,1));
grid_value  = div_map(sub2ind(size(div_map), ...
    round(grid_planes.ytable(grid_live)), round(grid_planes.xtable(grid_live))));

G = showpiv(sprintf('ev%d divergence per window', ev), 5);
G.update_figsize([7 5]);
G.plot_background(piv_ensemble.endpoint.stack(:,:,1));
G.plot_grid(grid_planes.xtable(grid_live), grid_planes.ytable(grid_live), ...
    grid_value, prctile(abs(grid_value), 95));
G.plot_quiver(piv_ensemble.endpoint.uv);
G.put_xaxistitle('x (px)');  G.put_yaxistitle('y (px)');
fprintf('%d windows drawn | %s | mean %+.5f\n', nnz(grid_live), div_info.unit, div_info.mean);

%% 8.2 The triangles the divergence theorem ran on
% Each triangle IS the cell whose gradient was computed, so painting it flat
% shows the partition and nothing invented. Slivers and the bridges Delaunay
% throws across the masked vessel are already dropped -- info.dropped counts them
[~, ~, tri_info] = vfield_strain(piv_ensemble.endpoint.uv, method = 'tri', ...
    exclmask = piv_ensemble.param.pivlab_mask);
T = showpiv(sprintf('ev%d divergence per triangle', ev), 5);
T.update_figsize([7 5]);
T.plot_background(piv_ensemble.endpoint.stack(:,:,1));
T.plot_tri(tri_info.xy, tri_info.tri, tri_info.div_tri, ...
    prctile(abs(tri_info.div_tri), 95));
T.put_xaxistitle('x (px)');  T.put_yaxistitle('y (px)');
fprintf('%d of %d triangles kept | long edge %d, sliver %d, masked %d\n', ...
    tri_info.dropped.kept, tri_info.dropped.total, tri_info.dropped.long_edge, ...
    tri_info.dropped.sliver, tri_info.dropped.masked);

%% 8.3 What the gates took
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
halfwin = round(0.5 * fps);                          % frames either side of from/to
param.windows = [40 20; 20 10; 12 6];                      % [window step] per pass, coarse to fine
S       = piv_preprocess(tp.(piv_ch));       % motion-preserving preprocess

% 5.1 piv_include is the only ROI drawn here; the vessel comes from setup_rois
redraw_rois = false;
if redraw_rois || ~has_roi(session.roilist, 'piv_include')
    session.roilist.addormodifyroi(twophoton_processed.ch1, 'piv_include', 'rectangle', twophoton_processed.ch2);
    session.roilist.save2disk(session.dir_struct.primary_analysis);
    fprintf('piv_include drawn and saved to %s\n', session.dir_struct.primary_analysis);
else
    fprintf('piv_include already present -- reusing (set redraw_rois=true to redraw)\n');
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
pick    = param.focus;                                     % subset of the labels, or {} for all
ev_list = 1:numel(events);
if ~isempty(pick); ev_list = find(ismember({events.state}, pick)); end

piv_run = piv_run_events(events, S, 'halfwin', halfwin, 'sel', ev_list, ...
    'piv_opt', {'window_sizes', param.windows, 'repeat', 1, 'do_pad', 1, ...
                'exclmask', exclmask, 'use_gpu', true});

%% 12. Perivascular divergence, every event
% test_label narrows to one transition while tuning; '' = every event in piv_run
test_label = 'nrem2rem';
sel = 1:numel(piv_run);
if ~isempty(test_label); sel = find(strcmp({piv_run.state}, test_label)); end
piv_sel = piv_run(sel);
fprintf('divergence on %d of %d events (%s)\n', numel(sel), numel(piv_run), ...
    strjoin(unique({piv_sel.state}), ', '));

% gated = anything was applied to the field after the correlation. TRUE here
% because piv_run_events calls piv_validate at ens_interleave:97 -- PIVlab's
% outlier removal, NOT the tomasi / corr_floor gates that analysis_pivensemble
% applies. The flag is coarse on purpose: what it protects is pooling rows that
% were filtered differently, and vfield_profile.append refuses a mixed set.
% Whether either kind of filtering belongs here is a separate question with a
% measured answer -- see CLAUDE_LOG.md, the gates remove the larger vectors
[vprofile, vfields] = piv_polar_events(piv_sel, coremask, p2u, ...
    exclmask = exclmask, n_wedge = 12, bin_edges_um = 0:1.5:40, ...
    min_tri_wedge = 10, gated = true, verbose = true);

%% 13. Batch figures
% ctx carries everything the figures need; piv_figure holds the shared quiver scale
% and colour limits
ctx = struct('S', S, 'coremask', coremask, 'exclmask', exclmask, 't_axis', t_axis, ...
             'dtrace', dtrace, 'dsg', dsg, 'fps', fps, 'p2u', p2u, ...
             'halfwin', halfwin, 'events', events, 'state_idx', img_state.state_idx);
F = piv_figure(piv_sel, vfields, vprofile, ctx);
T = F.summary();

%% 13.1 Sign and radial-profile overview
F.signcheck();                 % dD vs measured direction -- the polarities must split
F.radprofile();                % volume_out vs distance, ONE wedge set for every row

%% 13.2 One event in full
F.panel(1);                    % disp_radial | divergence | volume_out | diameter

%% 13.3 All events of one polarity
F.tiles('dilation');           % shared colour limits
F.tiles('constriction');

%% 13.4 Raw frames PIV used (blink comparator)
F.slices(1);

%% 14. Assemble piv_events / piv_result  (not saved -- inspect only)
% 9.1 Polar frame from the pax axis and the vessel centroid
paxv      = session.roilist.getvertices('pax');
pax_angle = atan2d(paxv(2,2) - paxv(1,2), paxv(2,1) - paxv(1,1));
[vy, vx]  = find(coremask);  center = [mean(vx), mean(vy)];
azframe   = polar_frame(size(coremask), center, 'pax_angle', pax_angle, 'n_sectors', 8);

% 9.2 Detection record
piv_events = struct('source_ref', 'fwhm', 'events', events, 'focus', {param.focus}, ...
    'detect', struct('params', det.P, 'sign_ok', det.sign_ok, ...
                     'clipped', det.clipped, 'dropped', det.dropped));

% 9.3 PIV record
piv_result = struct();
piv_result.meta    = struct('pixel2um', p2u, 'outfps', fps, 'halfwin', halfwin, ...
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
