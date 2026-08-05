%PIV_INTEGRATION_TESTBED  Diameter-event-driven PIV pipeline, run cell by cell.
%   Stages:
%     FWHM diameter -> merge bouts -> pick excursions -> review
%       -> per-event ensemble PIV -> per-event polar divergence -> assemble
%
%   Run analysis_main through FWHM and setup_rois first: this reads .pax_fwhm and
%   the 'pax' / 'manual_dilated_pvs' ROIs. Only 'piv_include' is drawn here.
%   Cell 0 loads the session itself.
%
%   SETTLED 260803, all measured on this session:
%     window_sizes is [window step] PER PASS. A single [40 20] row is one coarse
%       pass and no refinement; the 3-pass ladder below moved displacement over
%       its own resolution limit (disp/peak-halfwidth 0.35 -> 1.12, 149 -> 1537
%       vectors, |disp| 0.687 -> 0.926 because a 40 px window under-reads it).
%     Ensembling more pairs does NOT fix a wide peak: 5 -> 53 pairs moved the
%       half-width 1.99 -> 1.93 px, so the width is window geometry, not frame
%       noise. Event span caps the pair count anyway (framesets must not overlap).
%     corr2 cannot see this failure -- blur keeps it high (0.73 at the wall) while
%       the peak spreads, hence the min_snr gate in piv_validate.
%     coremask is NOT grown: it is already the maximum-dilation envelope, covers
%       the FWHM PVS wall on all 5462 frames (min margin 0.9 px), and a 6 px grow
%       changed ring-1 v_r by 2% while costing cells.
%     32 events: v_r sign matches the prescribed polarity 32/32, dv_r/dr 27/32
%       (the 5 are dilations sitting at dv_r/dr ~ 0).
%   OPEN: ring 1 v_r runs below rings 2-3 in BOTH polarities. Not window straddle
%     (survives the grow), not coverage (synthetic pure-radial field is exact at
%     17% coverage), not rigid translation or rotation (v_theta is 15-20% of v_r).

%% 0. Load session, stack and sleep state

% 0.1 Paths
addpath(genpath(pwd));
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

%% 1. Diameter source  (FWHM)
% Absolute BV diameter: thickness.bv = lower - upper BV boundary, px -> um.
% Not .bvchanges (median-subtracted) -- the amplitude bands need true caliber.
dtrace = double(session.pax_fwhm.thickness.bv(:).') * p2u;      % um
if numel(dtrace) ~= numel(t_axis)
    error('piv_testbed:lengthMismatch', ...
        'dtrace has %d samples but t_axis has %d -- they must be the same series.', ...
        numel(dtrace), numel(t_axis));
end
figure; plot(t_axis, dtrace); grid on
xlabel('t (s)'); ylabel('diameter (\mum)'); title('BV diameter (FWHM)');

%% 2. Merge adjacent bouts into transition pairs

MERGE = { ...
%    earlier(A)    later(B)      label             polarity
    'nrem',       'rem',        'nrem2rem',       'dilation'     ; ...
    'uarousal',   'nrem',       'uarousal2nrem',  'dilation'     ; ...
    'roughawake', 'nrem',       'awake2nrem',     'dilation'     ; ...
    'rem',        'roughawake', 'rem2awake',      'constriction' ; ...
    'nrem',       'roughawake', 'nrem2awake',     'constriction' ; ...
    'nrem',       'uarousal',   'nrem2uarousal',  'constriction' };

pairs = piv_merge_bouts(img_state.state_idx, MERGE, t_axis, 30);   % 30 s each side
fprintf('merged %d bout pairs\n', numel(pairs));

%% 3. Pick the excursion endpoints
%   sg_win_s = 3 s picking window, peaktol = 0.80 of the way to the window extreme
[events, det] = piv_pick_excursions(dtrace, t_axis, pairs, 3, 0.80);
fprintf('picked %d events (dropped: %s)\n', numel(events), ...
    strjoin(cellfun(@(f) sprintf('%s %d', f, det.dropped.(f)), ...
                    fieldnames(det.dropped)', 'uni', 0), ', '));
if det.sign_ok
    fprintf('sign check : all %d events match their prescribed polarity\n', numel(events));
else
    warning('piv_testbed:signMismatch', 'events %s contradict their polarity', ...
        mat2str(det.sign_bad));
end
if ~isempty(det.clipped)
    fprintf('edge check : %d event(s) hit a search-window edge (dD is a lower bound): %s\n', ...
        numel(det.clipped), mat2str(det.clipped));
end
dsg = det.dsg;              % picking trace, also what every figure marks on
disp(struct2table(rmfield(events, {'search_from','search_to'})));

FOCUS  = MERGE(:, 3)';                         % labels to run PIV on
STATES = {'roughnrem', 'rem', 'roughawake'};   % merge partition, for the review figure

%% 4. Review the picked events
piv_review_events(events, det, dtrace, t_axis, ...
    'sleep_score', sleep_integrate.sleep_score, ...
    'state_idx',   img_state.state_idx, ...
    'states',      STATES);

%% 5. PIV masks + preprocessed stack
% piv_include        : analyze ONLY inside this region (true = keep)
% manual_dilated_pvs : vessel + PVS at maximum dilation, traced by setup_rois. It is
%                      dropped from the PIV field and is the blob the polar rings
%                      are measured from.
piv_ch  = 'ch1';                                     % PIV channel -- set per dataset
halfwin = round(0.5 * fps);                          % frames either side of from/to
WINDOWS = [40 20; 20 10; 12 6];                      % [window step] per pass, coarse to fine
S       = opticalflow_preprocess(tp.(piv_ch));       % motion-preserving preprocess

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

%% 6. Run PIV
% ~20 s per event, so review cell 4 first
pick    = FOCUS;                                     % subset of the labels, or {} for all
ev_list = 1:numel(events);
if ~isempty(pick); ev_list = find(ismember({events.state}, pick)); end

piv_run = piv_run_events(events, S, 'halfwin', halfwin, 'sel', ev_list, ...
    'piv_opt', {'window_sizes', WINDOWS, 'repeat', 1, 'do_pad', 1, ...
                'exclmask', exclmask, 'use_gpu', true});

%% 7. Perivascular divergence
% test_label narrows to one transition while tuning; '' = every event in piv_run
test_label = 'nrem2rem';
sel = 1:numel(piv_run);
if ~isempty(test_label); sel = find(strcmp({piv_run.state}, test_label)); end
piv_sel = piv_run(sel);
fprintf('divergence on %d of %d events (%s)\n', numel(sel), numel(piv_run), ...
    strjoin(unique({piv_sel.state}), ', '));

polar = piv_polar_events(piv_sel, coremask, 'exclmask', exclmask, ...
    'n_sectors', 8, 'min_valid', 4, 'box_step', 3, 'verbose', true);

%% 8. Figures
% ctx carries everything the figures need; piv_figure holds the shared quiver scale
% and colour limits
ctx = struct('S', S, 'coremask', coremask, 'exclmask', exclmask, 't_axis', t_axis, ...
             'dtrace', dtrace, 'dsg', dsg, 'fps', fps, 'p2u', p2u, ...
             'halfwin', halfwin, 'events', events, 'state_idx', img_state.state_idx);
F = piv_figure(piv_sel, polar, ctx);
T = F.summary();

%% 8.1 Sign and radial-profile overview
F.signcheck();                 % dD vs measured direction -- the polarities must split
F.radprofile();                % dv_r/dr vs distance from the wall, all events

%% 8.2 One event in full
F.panel(1);                    % v_r | dv_r/dr | plain div | diameter

%% 8.3 All events of one polarity
F.tiles('dilation');           % shared colour limits
F.tiles('constriction');

%% 8.4 Raw frames PIV used (blink comparator)
F.slices(1);

%% 9. Assemble piv_events / piv_result  (not saved -- inspect only)
% 9.1 Polar frame from the pax axis and the vessel centroid
paxv      = session.roilist.getvertices('pax');
pax_angle = atan2d(paxv(2,2) - paxv(1,2), paxv(2,1) - paxv(1,1));
[vy, vx]  = find(coremask);  center = [mean(vx), mean(vy)];
azframe   = polar_frame(size(coremask), center, 'pax_angle', pax_angle, 'n_sectors', 8);

% 9.2 Detection record
piv_events = struct('source_ref', 'fwhm', 'events', events, 'focus', {FOCUS}, ...
    'detect', struct('params', det.P, 'sign_ok', det.sign_ok, ...
                     'clipped', det.clipped, 'dropped', det.dropped));

% 9.3 PIV record
piv_result = struct();
piv_result.meta    = struct('pixel2um', p2u, 'outfps', fps, 'halfwin', halfwin, ...
    'channel', piv_ch);
piv_result.frame   = struct('center', center, 'pax_angle', pax_angle, 'n_sectors', 8, ...
    'ring_width', polar{1}.info.ring_width, 'ring_edges', polar{1}.info.ring_edges);
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
