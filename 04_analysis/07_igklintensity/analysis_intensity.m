% ANALYSIS_INTENSITY  IgKL brightness inside the PVS.
%
% The PVS is read through a per-frame mask: on the thicker side, the rows between
% the vessel wall and the PVS outer edge that FWHM already found. Read as a MEAN
% -- the sum over the mask mostly re-measures how thick the PVS is, which FWHM has
% already said; the mean is the label concentration inside it, which it has not.
%
% An event runs from one trough of the diameter, through the peak, to the NEXT
% trough. Both ends are stationary points, so they match in phase.
%
% Pooling brightness across recordings would be wrong -- laser power and PMT gain
% differ. What is pooled is not brightness: every number is a ratio between two
% masks read in the same frames, or a fraction of an event's own onset level, so
% the gain has cancelled inside each one before it leaves its session.
%
% THE UNIT IS THE VESSEL. Sessions collapse to their vessel first, and the
% across-vessel numbers are taken on those. Depth belongs to the SESSION, not the
% vessel -- the same arteriole was imaged at several depths -- so the filter keeps
% shallow SESSIONS and a vessel contributes only the shallow ones.
%
% twophoton_preprocess ran medfilt3 over param.preproc_medianframes frames per
% fixed pixel before these kymographs existed. Anything faster than that in time
% is the filter, not the tissue.

clear
addpath('g:\03_program\01_ecspress\functions');   % where util_ecspath lives
util_ecspath;                                     % three roots, minus zz_notinuse

dirs = util_centraldirs();
dirs.dirtable = fullfile(dirs.secondary_root, 'merged_igkl_igkltdt', 'merged_igkl_igkltdt_dirtable.xlsx');
% a session temp folder was written down here and then deleted with the session.
% tempdir is the same kind of place and always exists
dirs.scratch = tempdir;

param.vessel_prefix        = "PA";  % str    penetrating arterioles
param.depth_thr_um         = 70;    % float  tableManager.m sets the L1 / L2 boundary here
param.paren_gap_px         = 3;     % int    gap outside the PVS edge before the reference band
param.paren_width_px       = 10;    % int    reference band width
param.preproc_medianframes = 11;    % int    twophoton_preprocess.medianframes, for the caveat
param.smooth_frames        = 9;     % int    movmean on the diameter before its landmarks are read
param.min_separation       = 12;    % int    minimum frames between troughs and peaks
param.min_prominence_pct   = 10;    % float  a trough or peak has to stand this far out of its
                                    %        surroundings, as a percentage of the vessel's own
                                    %        median diameter. A fixed um threshold is a
                                    %        different fraction of every vessel
param.max_duration_s       = 30;    % float  longer than this and the event is dropped
param.edge_frames          = 2;     % int    frames averaged at each end to read a level
param.half_columns         = 40;    % int    columns per half of the stretched event axis
param.margin_frac          = 0.20;  % float  margin outside each end, in stretched units
param.min_event_session    = 5;     % int    fewer events than this and the session is left out
param.pvsedge_band_px      = [-20 -11; -10 -1; 1 10; 11 20];
                                    % Bx2 int  bands measured from the PVS-parenchyma
                                    %        boundary in rows, negative going INWARD into
                                    %        the PVS and positive OUTWARD into parenchyma.
                                    %        All four share one anchor, so they slide
                                    %        together and their differences are not one
                                    %        anchor moving against another. Inward bands
                                    %        are clipped at the vessel wall.
                                    %        The offsets are PIXELS, so one band is a
                                    %        different physical distance in every session,
                                    %        and a fixed length inside a PVS whose width
                                    %        moves does not stay at one distance. The
                                    %        outward bands are in parenchyma and do not
                                    %        have this problem -- see CLAUDE_LOG.md
param.radial_grid_um       = 0:0.2:8;  % 1xU float  distance outside the wall, in um.
                                    %        NOT in PVS widths. A ruler tied to the PVS
                                    %        boundary puts a feature at that boundary
                                    %        whichever way it is handled -- scaled with the
                                    %        instantaneous thickness it moves with the very
                                    %        edge being tested, frozen at the onset it sits
                                    %        on the steepest part of the profile while the
                                    %        edge slides underneath. An absolute axis cannot
                                    %        find the boundary, and answers only the question
                                    %        that is left: does the response fall off with
                                    %        distance from the blood, as attenuation would

% --- which sessions ----------------------------------------------------------------
dir_table = readtable(dirs.dirtable, 'TextType', 'string');
vessel_tag = fillmissing(dir_table.VesselID, 'constant', "BLANK");
vessel_tag = upper(regexprep(vessel_tag, '[^A-Za-z0-9]', ''));
vessel_key = upper(dir_table.MouseID) + "_" + vessel_tag;
depth_um = str2double(extractBefore(fillmissing(dir_table.Depth, 'constant', "") + "um", "um"));

keep = dir_table.Primary_paxFWHM == "paxfwhm.mat" ...
     & startsWith(vessel_tag, param.vessel_prefix) ...
     & isfinite(depth_um) & depth_um <= param.depth_thr_um;
session_dir = dir_table.Directory(keep);
session_vessel = vessel_key(keep);

% The dirtable is a scan and folders have been renamed since. Missing ones are
% counted here rather than left to fail one at a time
on_disk = isfolder(session_dir);
if ~all(on_disk)
    fprintf('%d of %d directories no longer exist, skipped\n', nnz(~on_disk), numel(on_disk));
end
session_dir = session_dir(on_disk);
session_vessel = session_vessel(on_disk);
fprintf('%d sessions, %d vessels, %d mice, depth <= %.0f um\n\n', numel(session_dir), ...
    numel(unique(session_vessel)), numel(unique(extractBefore(session_vessel, "_"))), ...
    param.depth_thr_um);

% --- one row per session -------------------------------------------------------------
frac_edge = 1 + param.margin_frac;
frac_axis = linspace(-frac_edge, frac_edge, 2 * param.half_columns + 1);
peak_col = param.half_columns + 1;
n_session = numel(session_dir);
row = struct('vessel', {}, 'session', {}, 'n_event', {}, 'amp_um', {}, 'duration_s', {}, ...
    'gap_um', {}, 'ratio_static', {}, 'thickness_corr', {}, 'pvs_um', {}, ...
    'peak_moving', {}, 'peak_band', {}, 'band_um', {}, 'band_clipped', {}, ...
    'peak_pvs', {}, 'peak_paren', {}, ...
    'peak_diff', {}, 'end_pvs', {}, 'end_paren', {}, 'end_diff', {}, ...
    'curve_pvs', {}, 'curve_paren', {}, 'curve_dD', {}, 'radial', {}, 'curve_band', {});

failed = struct('session', {}, 'message', {});
for s = 1:n_session
    % one unreadable session must not cost the other seventy-odd. Failures are
    % collected and named at the end rather than passed over quietly
    try
        this = session_result(session_dir(s), ...
            paren_gap_px = param.paren_gap_px, paren_width_px = param.paren_width_px, ...
            smooth_frames = param.smooth_frames, min_separation = param.min_separation, ...
            min_prominence_pct = param.min_prominence_pct, max_duration_s = param.max_duration_s, ...
            edge_frames = param.edge_frames, frac_axis = frac_axis, ...
            radial_grid_um = param.radial_grid_um, pvsedge_band_px = param.pvsedge_band_px);
    catch failure
        failed(end+1).session = session_dir(s);        %#ok<SAGROW>
        failed(end).message = string(failure.message);
        fprintf('%3d/%d  %-32s FAILED  %s\n', s, n_session, session_dir(s), failure.message);
        continue
    end
    if this.n_event < param.min_event_session
        fprintf('%3d/%d  %-32s %-14s  %d events -- left out\n', s, n_session, ...
            this.session, session_vessel(s), this.n_event);
        continue
    end
    row(end+1) = mergefields(this, session_vessel(s), peak_col);    %#ok<SAGROW>
    fprintf('%3d/%d  %-32s %-14s  n %3d  amp %.2f  dur %5.1f  diff %+.5f\n', s, n_session, ...
        this.session, session_vessel(s), row(end).n_event, row(end).amp_um, ...
        row(end).duration_s, row(end).peak_diff);
end

% --- sessions collapse to their vessel ------------------------------------------------
[vessel_name, ~, vessel_group] = unique([row.vessel]);
n_vessel = numel(vessel_name);
scalar_names = ["amp_um", "duration_s", "gap_um", "ratio_static", "thickness_corr", "pvs_um", ...
    "peak_moving", ...
    "peak_pvs", "peak_paren", "peak_diff", "end_pvs", "end_paren", "end_diff"];
vessel = struct();
for name = scalar_names
    per_vessel = nan(1, n_vessel);
    for v = 1:n_vessel
        per_vessel(v) = median([row(vessel_group == v).(name)], 'omitnan');
    end
    vessel.(name) = per_vessel;
end
curve_pvs = nan(n_vessel, numel(frac_axis));
curve_paren = nan(n_vessel, numel(frac_axis));
curve_dD = nan(n_vessel, numel(frac_axis));
curve_radial = nan(n_vessel, numel(param.radial_grid_um));
n_band = size(param.pvsedge_band_px, 1);
curve_band = nan(n_band, numel(frac_axis), n_vessel);
peak_band = nan(n_band, n_vessel);
for v = 1:n_vessel
    curve_pvs(v,:) = median(vertcat(row(vessel_group == v).curve_pvs), 1, 'omitnan');
    curve_paren(v,:) = median(vertcat(row(vessel_group == v).curve_paren), 1, 'omitnan');
    curve_dD(v,:) = median(vertcat(row(vessel_group == v).curve_dD), 1, 'omitnan');
    curve_radial(v,:) = median(vertcat(row(vessel_group == v).radial), 1, 'omitnan');
    curve_band(:,:,v) = median(cat(3, row(vessel_group == v).curve_band), 3, 'omitnan');
    peak_band(:,v) = median(vertcat(row(vessel_group == v).peak_band), 1, 'omitnan');
end

if ~isempty(failed)
    fprintf('\n%d sessions failed to read:\n', numel(failed));
    for f = 1:numel(failed)
        fprintf('   %s  --  %s\n', failed(f).session, failed(f).message);
    end
end
fprintf('\n%d sessions kept, %d vessels\n', numel(row), n_vessel);
fprintf('events per session  median %d,  amplitude median %.2f um,  duration median %.1f s\n', ...
    round(median([row.n_event])), median(vessel.amp_um, 'omitnan'), ...
    median(vessel.duration_s, 'omitnan'));
fprintf('PVS mask / reference  median %.4f  (%d of %d vessels finite)\n', ...
    median(vessel.ratio_static, 'omitnan'), nnz(isfinite(vessel.ratio_static)), n_vessel);
fprintf('mask mean vs PVS thickness  median r %+.3f\n', median(vessel.thickness_corr));

fprintf('\n=== vessel level, n %d ===\n', n_vessel);
fprintf('%-12s %12s %12s %10s %8s\n', '', 'median', 'IQR', 'negative', 'p');
report_row('peak PVS', vessel.peak_pvs);
report_row('peak ref', vessel.peak_paren);
report_row('peak diff', vessel.peak_diff);
report_row('end PVS', vessel.end_pvs);
report_row('end ref', vessel.end_paren);
report_row('end diff', vessel.end_diff);
fprintf('\n=== bands hung off the PVS-parenchyma boundary, inward and outward ===\n');
band_median = median(cat(3, row.band_um), 3);
clip_median = median(vertcat(row.band_clipped), 1);
fprintf('PVS is %.2f um at the vessel median. Negative is inward, into the PVS\n', ...
    median(vessel.pvs_um, 'omitnan'));
band_label = strings(1, n_band);
for b = 1:n_band
    band_label(b) = sprintf('%+.1f..%+.1f um', band_median(b,1), band_median(b,2));
    fprintf('   band %d  %-18s  clipped at the wall in %.1f%% of frames\n', ...
        b, band_label(b), 100 * clip_median(b));
end
fprintf('   pixel size differs between sessions, so a px offset is not one distance\n');
for b = 1:n_band
    report_row(sprintf('band %d', b), peak_band(b,:));
end
report_row('inner - outer', peak_band(1,:) - peak_band(n_band,:));
fprintf('   p is a signed-rank test against zero, ACROSS VESSELS\n');

% --- one figure ------------------------------------------------------------------------
clee = color_lee();
fig_pool = make_fig('intensity_pool');
fig_pool.xy_sizeinch = [5.6 3.4];
fig_pool.update_figsize(fig_pool.xy_sizeinch);
fig_pool.hold_axis(true);
% All four bands hang off the SAME edge, the PVS-parenchyma boundary, so the gaps
% between them cannot be one anchor moving against another. The PVS mask and its
% edge-anchored reference are not drawn: those two are anchored to different things
% -- the wall and the PVS edge -- which advance at different rates, and the gap
% between them takes its sign from that
band_color = [clee.clist.red; clee.clist.orange; clee.clist.teal; clee.clist.blue];
h_band = gobjects(1, n_band);
for b = 1:n_band
    h_band(b) = plot(fig_pool.ax, frac_axis, median(curve_band(b,:,:), 3, 'omitnan'), ...
        'Color', band_color(b,:), 'LineWidth', 2);
end
yline(fig_pool.ax, 0, 'k--', 'LineWidth', 1);
xline(fig_pool.ax, 0, 'k-', 'LineWidth', 1);
xlim(fig_pool.ax, [-frac_edge frac_edge]);
xticks(fig_pool.ax, [-1 -0.5 0 0.5 1]);
xticklabels(fig_pool.ax, {'trough', '', 'peak', '', 'next trough'});
fig_pool.put_xaxistitle('event time, stretched');
fig_pool.put_yaxistitle('fraction of onset level');
yyaxis(fig_pool.ax, 'right');
h_dD = plot(fig_pool.ax, frac_axis, median(curve_dD, 1, 'omitnan'), ...
    'Color', clee.clist.black, 'LineWidth', 1);
ylabel(fig_pool.ax, 'diameter from onset (um)');
fig_pool.ax.YAxis(2).Color = clee.clist.black;
yyaxis(fig_pool.ax, 'left');
legend([h_band h_dD], [cellstr(band_label), {'diameter'}], ...
    'Location', 'southwest', 'Box', 'off', 'AutoUpdate', 'off');
print(fig_pool.fig, fullfile(dirs.scratch, 'intensity_pool.png'), '-dpng', '-r150');

% --- the radial profile of the peak response ------------------------------------------
% Distance is in um from the wall and owes nothing to the PVS boundary. A response
% that decays monotonically from the wall is attenuation falling off with distance
% from the blood. A response that stays flat over the first few um and then drops is
% not, whatever else it turns out to be. The median PVS width is drawn for reference
% only -- the axis was not built from it
fig_r = make_fig('intensity_radial');
fig_r.xy_sizeinch = [5.2 3.4];
fig_r.update_figsize(fig_r.xy_sizeinch);
fig_r.hold_axis(true);
for v = 1:n_vessel
    plot(fig_r.ax, param.radial_grid_um, curve_radial(v,:), 'Color', [clee.clist.red 0.2], ...
        'LineWidth', 0.5);
end
plot(fig_r.ax, param.radial_grid_um, median(curve_radial, 1, 'omitnan'), ...
    'Color', clee.clist.red, 'LineWidth', 2);
yline(fig_r.ax, 0, 'k--', 'LineWidth', 1);
xline(fig_r.ax, median(vessel.pvs_um, 'omitnan'), 'k:', 'LineWidth', 1);
xlim(fig_r.ax, [0 max(param.radial_grid_um)]);
fig_r.put_xaxistitle('distance outside the wall (um)');
fig_r.put_yaxistitle('response at the peak, fraction of onset');
print(fig_r.fig, fullfile(dirs.scratch, 'intensity_radial.png'), '-dpng', '-r150');

show_u = [0 0.5 1 1.5 2 3 4 5 6 8];
[~, col_u] = min(abs(param.radial_grid_um(:) - show_u), [], 1);
fprintf('\n=== radial profile of the peak response, vessel median ===\n');
fprintf('%-12s', 'um');
fprintf('%9.2f', show_u);
fprintf('\n%-12s', 'response');
fprintf('%+9.5f', median(curve_radial(:, col_u), 1, 'omitnan'));
fprintf('\n%-12s', 'negative');
fprintf('%6d/%-2d', [sum(curve_radial(:, col_u) < 0, 1); repmat(n_vessel, 1, numel(col_u))]);
fprintf('\n   monotone decay from the wall is attenuation falling off with distance\n');

save(fullfile(dirs.scratch, 'intensity_pool.mat'), 'row', 'vessel', 'vessel_name', 'peak_band', 'curve_band', 'band_label', ...
    'curve_pvs', 'curve_paren', 'curve_dD', 'curve_radial', 'frac_axis', 'peak_col', 'param');
fprintf('\nsaved to %s\n', dirs.scratch);

% ---------------------------------------------------------------------------
function report_row(label, value)
    ok = isfinite(value);
    fprintf('%-12s %+12.5f %12s %6d/%-3d %8.4f\n', label, median(value(ok)), ...
        sprintf('%+.5f..%+.5f', prctile(value(ok),25), prctile(value(ok),75)), ...
        nnz(value(ok) < 0), nnz(ok), signrank(value(ok)));
end

function out = mergefields(this, vessel_key, peak_col)
    out.vessel = vessel_key;
    out.session = this.session;
    out.n_event = this.n_event;
    out.amp_um = median(this.amp_um);
    out.duration_s = median(this.duration_s);
    out.gap_um = median(this.gap_um);
    out.ratio_static = this.ratio_static;
    out.thickness_corr = this.thickness_corr;
    out.pvs_um = this.pvs_um;
    out.peak_moving = this.peak_moving;
    out.peak_band = this.peak_band;
    out.band_um = this.band_um;
    out.band_clipped = this.band_clipped;
    out.curve_band = cell2mat(cellfun(@(m) mean(m, 1, 'omitnan'), ...
        this.event_band(:), 'UniformOutput', false));
    out.peak_pvs = mean(this.event_pvs(:, peak_col), 'omitnan');
    out.peak_paren = mean(this.event_paren(:, peak_col), 'omitnan');
    out.peak_diff = out.peak_pvs - out.peak_paren;
    out.end_pvs = median(this.end_pvs);
    out.end_paren = median(this.end_paren);
    out.end_diff = out.end_pvs - out.end_paren;
    out.curve_pvs = mean(this.event_pvs, 1, 'omitnan');
    out.curve_paren = mean(this.event_paren, 1, 'omitnan');
    out.curve_dD = mean(this.event_dD, 1, 'omitnan');
    out.radial = this.radial;
end

function out = session_result(session_dir, opt)
    arguments
        session_dir             (1,1) string
        opt.paren_gap_px        (1,1) double
        opt.paren_width_px      (1,1) double
        opt.smooth_frames       (1,1) double
        opt.min_separation      (1,1) double
        opt.min_prominence_pct  (1,1) double
        opt.max_duration_s      (1,1) double
        opt.edge_frames         (1,1) double
        opt.frac_axis           (1,:) double
        opt.pvsedge_band_px     (:,2) double
        opt.radial_grid_um      (1,:) double
    end

    [~, stem, ext] = fileparts(session_dir);
    out.session = string([char(stem) char(ext)]);

    loaded = load(fullfile(session_dir, 'primary_analysis', 'paxfwhm.mat'));
    pax_fwhm = loaded.line_fwhm;
    clear loaded

    info_text = string(fileread(fullfile(session_dir, out.session + "_info.txt")));
    objpix_token = regexp(info_text, 'objpix,([^\r\n]+)', 'tokens', 'once');
    pixel2um = str2double(extractBefore(strtrim(objpix_token{1}) + " ", " "));

    % param.fs was added to line_fwhm after some of these files were written, so
    % it comes off the time axis when the field is absent
    if isfield(pax_fwhm.param, 'fs')
        fs = pax_fwhm.param.fs;
    else
        fs = 1 / median(diff(pax_fwhm.t_axis));
    end
    kgph = double(pax_fwhm.kymograph.kgph_pvs);
    bv_um = pax_fwhm.thickness.bv(:).' * pixel2um;

    % the mask, on the side line_fwhm already found to be thicker. Outward is a
    % falling row index above the vessel and a rising one below
    if pax_fwhm.up_thicker
        wall_row = round(pax_fwhm.idx.clean_upperBVboundary(:).');
        pvs_row = round(pax_fwhm.idx.clean_pvsupedge_idx(:).');
        outward = -1;
    else
        wall_row = round(pax_fwhm.idx.clean_lowerBVboundary(:).');
        pvs_row = round(pax_fwhm.idx.clean_pvsdownedge_idx(:).');
        outward = 1;
    end
    pvs_thickness = outward * (pvs_row - wall_row);
    trace_pvs = mask_mean(kgph, wall_row + outward, pvs_row);
    paren_from = pvs_row + outward * opt.paren_gap_px;
    paren_to = paren_from + outward * (opt.paren_width_px - 1);
    trace_paren = mask_mean(kgph, paren_from, paren_to);

    % bands hung off the PVS-parenchyma boundary, inward and outward. Neither follows
    % the tissue -- nothing here does -- but they slide the same way, so the difference
    % between them is not contaminated by one anchor moving against the other. A
    % PVS-edge band next to a wall band is exactly that contamination: the PVS edge
    % advances by only about half of what the wall does.
    % A COMMON ANCHOR IS NOT ENOUGH. The band LENGTH is fixed in pixels while the PVS
    % width moves, so an inward band walks across the PVS intensity crest and the
    % innermost one is not a measurement at a fixed distance. band_clipped is the
    % fraction of frames it hit the wall; read it before quoting an inward band.
    % see CLAUDE_LOG.md
    n_band = size(opt.pvsedge_band_px, 1);
    trace_band = cell(1, n_band);
    out.band_clipped = nan(1, n_band);
    for b = 1:n_band
        from_row = pvs_row + outward * opt.pvsedge_band_px(b,1);
        to_row = pvs_row + outward * opt.pvsedge_band_px(b,2);
        % an inward band must not cross the wall into the lumen. Where it would, it
        % is clipped at the wall and the fraction of frames that happened to is kept
        if outward > 0
            clipped = from_row < wall_row + 1;
            from_row = max(from_row, wall_row + 1);
            to_row = max(to_row, wall_row + 1);
        else
            clipped = from_row > wall_row - 1;
            from_row = min(from_row, wall_row - 1);
            to_row = min(to_row, wall_row - 1);
        end
        out.band_clipped(b) = mean(clipped);
        trace_band{b} = mask_mean(kgph, from_row, to_row);
    end
    out.band_um = opt.pvsedge_band_px * pixel2um;
    out.ratio_static = mean(trace_pvs) / mean(trace_paren);
    out.thickness_corr = corr(trace_pvs(:), pvs_thickness(:));
    out.pvs_um = median(pvs_thickness) * pixel2um;

    n_u = numel(opt.radial_grid_um);

    % events, trough to trough
    bv_smooth = movmean(bv_um, opt.smooth_frames);
    prominence_um = opt.min_prominence_pct / 100 * median(bv_um);
    trough_idx = find(islocalmin(bv_smooth, 'MinSeparation', opt.min_separation, ...
        'MinProminence', prominence_um));
    peak_idx = find(islocalmax(bv_smooth, 'MinSeparation', opt.min_separation, ...
        'MinProminence', prominence_um));
    onset = [];
    peak = [];
    back = [];
    for j = 1 : numel(trough_idx) - 1
        between = peak_idx(peak_idx > trough_idx(j) & peak_idx < trough_idx(j+1));
        if isempty(between)
            continue
        end
        if (trough_idx(j+1) - trough_idx(j)) / fs > opt.max_duration_s
            continue
        end
        onset(end+1) = trough_idx(j);            %#ok<AGROW>
        peak(end+1) = between(1);                %#ok<AGROW>
        back(end+1) = trough_idx(j+1);           %#ok<AGROW>
    end
    out.n_event = numel(onset);
    out.amp_um = bv_smooth(peak) - bv_smooth(onset);
    out.gap_um = bv_smooth(back) - bv_smooth(onset);
    out.duration_s = (back - onset) / fs;

    % each event stretched onto the common axis, as a fraction of its own onset
    n_col = numel(opt.frac_axis);
    out.event_pvs = nan(out.n_event, n_col);
    out.event_paren = nan(out.n_event, n_col);
    out.event_dD = nan(out.n_event, n_col);
    out.end_pvs = nan(1, out.n_event);
    out.end_paren = nan(1, out.n_event);
    % the radial profile needs only the two ends and the peak, not the whole curve
    radial_event = nan(out.n_event, n_u);
    peak_event_pvs = nan(1, out.n_event);
    peak_event_paren = nan(1, out.n_event);
    peak_event_band = nan(out.n_event, n_band);
    out.event_band = cell(1, n_band);
    for b = 1:n_band
        out.event_band{b} = nan(out.n_event, numel(opt.frac_axis));
    end
    margin_frac = max(opt.frac_axis) - 1;
    for k = 1:out.n_event
        rise_seg = onset(k) : peak(k);
        fall_seg = peak(k) : back(k);
        if numel(rise_seg) < 2 || numel(fall_seg) < 2
            continue
        end
        step_rise = 1 / (numel(rise_seg) - 1);
        step_fall = 1 / (numel(fall_seg) - 1);
        before_seg = max(1, onset(k) - round(margin_frac / step_rise)) : onset(k) - 1;
        after_seg = back(k) + 1 : min(numel(bv_um), back(k) + round(margin_frac / step_fall));
        frac_rise = linspace(-1, 0, numel(rise_seg));
        frac_fall = linspace(0, 1, numel(fall_seg));
        frac_before = -1 - (numel(before_seg):-1:1) * step_rise;
        frac_after = 1 + (1:numel(after_seg)) * step_fall;
        frac_event = [frac_before, frac_rise, frac_fall(2:end), frac_after];
        segment = [before_seg, rise_seg, fall_seg(2:end), after_seg];

        onset_range = onset(k) : min(onset(k) + opt.edge_frames - 1, peak(k));
        back_range = max(back(k) - opt.edge_frames + 1, peak(k)) : back(k);
        base_pvs = mean(trace_pvs(onset_range));
        base_paren = mean(trace_paren(onset_range));
        out.event_pvs(k,:) = interp1(frac_event, trace_pvs(segment) / base_pvs - 1, opt.frac_axis);
        out.event_paren(k,:) = interp1(frac_event, trace_paren(segment) / base_paren - 1, opt.frac_axis);
        out.event_dD(k,:) = interp1(frac_event, bv_um(segment) - bv_um(onset(k)), opt.frac_axis);
        out.end_pvs(k) = mean(trace_pvs(back_range)) / base_pvs - 1;
        out.end_paren(k) = mean(trace_paren(back_range)) / base_paren - 1;

        % the response read at a series of distances outside the wall, so it can be
        % read against distance instead of in two bands. Distance is in um and owes
        % nothing to the PVS boundary -- see param.radial_grid_um for why a ruler
        % tied to that boundary produces a feature there whichever way it is handled
        peak_range = max(peak(k) - opt.edge_frames + 1, onset(k)) : ...
            min(peak(k) + opt.edge_frames - 1, back(k));

        % the same three bands read at the peak, so the wall-anchored reference can
        % be put beside the edge-anchored one on identical events
        peak_event_pvs(k) = mean(trace_pvs(peak_range)) / base_pvs - 1;
        peak_event_paren(k) = mean(trace_paren(peak_range)) / base_paren - 1;
        for b = 1:n_band
            base_band = mean(trace_band{b}(onset_range));
            peak_event_band(k,b) = mean(trace_band{b}(peak_range)) / base_band - 1;
            out.event_band{b}(k,:) = interp1(frac_event, ...
                trace_band{b}(segment) / base_band - 1, opt.frac_axis);
        end

        for u = 1:n_u
            offset_px = outward * round(opt.radial_grid_um(u) / pixel2um);
            level_onset = row_mean(kgph, wall_row(onset_range) + offset_px, onset_range);
            level_peak = row_mean(kgph, wall_row(peak_range) + offset_px, peak_range);
            radial_event(k,u) = level_peak / level_onset - 1;
        end
    end
    out.radial = mean(radial_event, 1, 'omitnan');
    out.peak_moving = mean(peak_event_paren, 'omitnan');
    out.peak_band = mean(peak_event_band, 1, 'omitnan');
    out.peak_mask = mean(peak_event_pvs, 'omitnan');
end

function level = row_mean(kgph, read_row, frame_idx)
    % One row per frame, averaged over the frames given. The row moves with the
    % wall; frames whose row falls outside the crop drop out
    [n_row, ~] = size(kgph);
    inside = read_row >= 1 & read_row <= n_row;
    if ~any(inside)
        level = NaN;
        return
    end
    linear_idx = sub2ind(size(kgph), read_row(inside), frame_idx(inside));
    level = mean(kgph(linear_idx));
end

function trace = mask_mean(kgph, row_from, row_to)
    % Mean over the rows between two moving boundaries, one frame at a time.
    % Either order is allowed, so the same call serves both sides of the vessel
    [n_row, n_frame] = size(kgph);
    trace = nan(1, n_frame);
    for t = 1:n_frame
        first_row = max(1, min(row_from(t), row_to(t)));
        last_row = min(n_row, max(row_from(t), row_to(t)));
        trace(t) = mean(kgph(first_row:last_row, t));
    end
end
