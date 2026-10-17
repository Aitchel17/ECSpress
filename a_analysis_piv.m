%A_ANALYSIS_PIV  The second main : one session's diameter transitions through PIV, as pivtransition.mat and pivpolar.mat
%   Runs on a_analysis_pvs's products (paxfwhm.mat, roilist.mat), the sleep score and a_check_pivquality's pass.
%   Inputs - stacks - ROIs - transitions - mask - PIV - save - QC - polar.
clc, clear
addpath('g:/03_program/01_ecspress/09_dirstruct');    % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

%% settings
dirs.session = 'K:\tmp\02_egfptdt\hql099\260601_hql99_sleep\HQL99_sleep260601_005';   % char  the session folder; tdTomato / eGFP cohorts only, igkl has no texture to track
param.halfwin_s     = 0.5;                    % frames either side of an endpoint, as seconds (s)
param.windows       = [40 20; 20 10; 12 6];   % [window step] per pass, coarse to fine (px)
param.roi_exclude   = "manual_dilated_pvs";   % ROI dropped from the field : vessel + PVS at max dilation, drawn in cell 3 (str)
param.roi_include   = "piv_include";          % ROI the field is confined to, drawn in cell 3 (str)
param.n_wedge       = 12;                     % angular bins around the vessel centroid (count)
param.bin_edges_um  = 0:1.5:40;               % radial bin edges, distance from the wall (um)
param.min_tri_wedge = 10;                     % triangles a wedge needs before it counts (count)

%% 1. Inputs : load session files (FWHM and sleep score)
session = ECSSession(dirs.session);
session = session.load_primary_results;
load(fullfile(session.dir_struct.primary_analysis, 'pivquality_pass.mat'), 'pivquality');   % a session not passed by a_check_pivquality stops here
score = load(fullfile(dirs.session, 'peripheral', 'sleep_score.mat'));
sleep_bout = bout_sleep(score.behavState, score.statecodes, score.binwidth_sec);
channel = char(session.pax_fwhm.param.channel_pvs);        % the tissue channel, off the FWHM record

%% 2. Load stack
session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
twophoton_processed = twophoton_preprocess(session);

%% 3. ROIs : the wall and the PIV region drawn or confirmed
stack = twophoton_processed.(channel);
% the wall at maximum dilation, on setup_rois's window around the widest frame, as a_analysis_pvs 4.1.5 draws it
setup_rois(session.roilist, twophoton_processed, session.pax_fwhm, session.dir_struct.primary_analysis, ...
    labels = param.roi_exclude);
% the region PIV is confined to
session.roilist.addormodifyroi(stack, char(param.roi_include), 'rectangle');
session.roilist.save2disk(session.dir_struct.primary_analysis);

%% 4. Transitions on the diameter, analysis_event's default merge table. control() is not called : transitions only
state_idx = struct();
for table_name = string(fieldnames(sleep_bout))'
    state_idx.(table_name) = timetable2frame(twophoton_processed.t_axis, sleep_bout.(table_name));   % each bout table on the stack's frames
end
if numel(session.pax_fwhm.t_axis) ~= size(twophoton_processed.ch1, 3)
    error('a_analysis_piv:frameAxis', ...
        'paxfwhm.mat has %d frames, the stack %d -- the two were not made from the same extraction', ...
        numel(session.pax_fwhm.t_axis), size(twophoton_processed.ch1, 3));
end
event_det = analysis_event(session.pax_fwhm.thickness.bv, session.pax_fwhm.t_axis);
event_det.arousaltransition(state_idx);

%% 5. The mask. PIVlab convention : true = pixel excluded
excl = session.roilist.getmask(char(param.roi_exclude));
incl = session.roilist.getmask(char(param.roi_include));
exclmask = ~incl | excl;

%% 6. PIV, one transition at a time : a refused window drops one row, not the session
% the margin at the recording ends (piv_getframescope) and the frameset overlap
% (analysis_pivensemble) both refuse by error
halfwin = round(param.halfwin_s * twophoton_processed.outfps);
n_event = numel(event_det.eventlist);
rows    = cell(1, n_event);
refused = struct('id', {}, 'state', {}, 'message', {});
for e = 1:n_event
    event = event_det.eventlist(e);
    fprintf('[%d/%d] %-14s %-12s %5.1f s ', e, n_event, event.state, event.pol, event.rise_s);
    try
        rows{e} = piv_run_events(event_det.eventlist, stack, twophoton_processed.fps, ...
            twophoton_processed.pixel2um, halfwin = halfwin, sel = e, ...
            windows = param.windows, exclmask = exclmask, verbose = false);
        fprintf('| %d vectors\n', nnz(~isnan(rows{e}.xyuv(:, :, 3))));
    catch err
        refused(end+1) = struct('id', e, 'state', event.state, 'message', string(err.message)); %#ok<SAGROW>
        fprintf('| refused : %s\n', err.message);
    end
end
% one row per transition that went through, in eventlist order
made = ~cellfun(@isempty, rows);
pivtransition = [rows{made}];
fprintf('%d transitions -> %d rows, %d refused\n', n_event, numel(pivtransition), numel(refused));

%% 7. Save : info holds what was measured or resolved, param the settings
info = struct();
info.fps         = twophoton_processed.fps;         % Hz, after resampling
info.pixel2um    = twophoton_processed.pixel2um;    % um per px
info.n_frame     = size(stack, 3);
info.channel     = string(channel);                 % which stack PIV ran on
info.halfwin     = halfwin;                         % frames, param.halfwin_s at this fps
info.eventlist   = event_det.eventlist;             % every transition detected, refused ones included
info.event_param = event_det.param;                 % the partition the detection ran on
info.event_info  = event_det.info;                  % what the picking measured
info.preprocess  = rmfield(twophoton_processed, {'ch1', 'ch2'});   % the resampling, minus the stacks
info.refused     = refused;

% beside paxfwhm.mat and roilist.mat, where the primary products of a session live
out_path = fullfile(session.dir_struct.primary_analysis, 'pivtransition.mat');
if isfile(out_path)
    backup_stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    backup_path  = fullfile(session.dir_struct.primary_analysis, ['pivtransition_backup_' backup_stamp '.mat']);
    movefile(out_path, backup_path);
    fprintf('previous file kept as %s\n', backup_path);
end
save(out_path, 'pivtransition', 'param', 'info');
listing = dir(out_path);
fprintf('%d rows, %.1f MB -> %s\n', numel(pivtransition), listing.bytes / 1e6, out_path);

%% 8. QC : per state, its largest and its median |dD| row, the field beside the FROM/TO blink
row_state = string({pivtransition.state});
abs_dD = abs([pivtransition.diameter_change]);
frame_hw = size(stack, [1 2]);
qc_panel_in = [5.2, 5.2 * frame_hw(1) / frame_hw(2) + 0.8];   % one panel, inches; the top 0.8 is the title
for state = unique(row_state, 'stable')
    idx_state = find(row_state == state);
    [~, order] = sort(abs_dD(idx_state));
    show_idx = unique(idx_state(order([end, ceil(end / 2)])), 'stable');   % largest, then the median (the lower middle)
    show_tag = ["largest |dD|", "median |dD|"];
    if isscalar(show_idx)
        show_tag = "only row";
    end
    board = fusion_figure('grid', [numel(show_idx) 2], 'width', 2 * qc_panel_in(1), ...
        'height', numel(show_idx) * qc_panel_in(2));
    board.fig.Name = char("PIV QC  " + state + "  arrows x10, every other vector");
    for k = 1:numel(show_idx)
        row = pivtransition(show_idx(k));
        from_img = mean(stack(:, :, row.pair_frames(:, 1)), 3);   % the frame sets PIV averaged
        to_img = mean(stack(:, :, row.pair_frames(:, 2)), 3);
        n_vector = nnz(~isnan(row.xyuv(:, :, 3)));
        head_line = sprintf('%s  ev%d %s %s', show_tag(k), row.id, row.state, row.pol);
        value_line = sprintf('dD %+.2f um  rise %.1f s  %d vectors, mean of %d frames each', ...
            row.diameter_change * twophoton_processed.pixel2um, row.rise_s, n_vector, row.nfr);
        mode_line = sprintf('common mode removed %+.2f %+.2f px', row.common_mode(1), row.common_mode(2));
        [quiver_fig, blink_fig] = vfield_qcfig(from_img, to_img, row.xyuv);
        title(quiver_fig.ax, {head_line, value_line, mode_line}, 'Interpreter', 'none', 'FontSize', 7);
        title(blink_fig.ax, {head_line, value_line, mode_line}, 'Interpreter', 'none', 'FontSize', 7);
        board.add(quiver_fig.ax, 'at', [k 1]);
        board.add(blink_fig.ax, 'at', [k 2]);
        close([quiver_fig.fig, blink_fig.fig]);
    end
    plot_blink(board.fig, 'blink_image');
end

%% 9. Polar : every row's field read against the wall at maximum dilation, saved as pivpolar.mat
% rebuilt from the rows in seconds, so an overwrite keeps no backup
[polar_profile, polar_fields] = piv_polar_events(pivtransition, excl, twophoton_processed.pixel2um, ...
    n_wedge = param.n_wedge, bin_edges_um = param.bin_edges_um, exclmask = exclmask, ...
    min_tri_wedge = param.min_tri_wedge);
mesh_rows = cellfun(@(field) field.info, polar_fields, 'UniformOutput', false);
polar_info = struct();
polar_info.mesh = [mesh_rows{:}];   % 1 x N struct, per row what vfield_polar measured : vectors, triangles kept and placed
polar_param = rmfield(param, {'halfwin_s', 'windows'});
pivpolar = polar_profile.to_struct();
polar_file = struct('pivpolar', pivpolar, 'param', polar_param, 'info', polar_info);
polar_path = fullfile(session.dir_struct.primary_analysis, 'pivpolar.mat');
save(polar_path, '-struct', 'polar_file');
fprintf('%d rows, %d bins x %d wedges -> %s\n', polar_profile.n_row, polar_profile.n_bin, polar_profile.n_wedge, polar_path);

%% 10. Polar figures : wall motion against dD, then divergence and radial strain against distance from the wall
% dD is stored in px and drawn in um. All three figures read the same rows, bins and wedges : the widest span
% common to every row, so no point or curve averages over a set of cells its neighbours do not share
dD_um = [pivtransition.diameter_change] * twophoton_processed.pixel2um;
pol_rows = string({pivtransition.pol});
row_rgb = color_lee.polarity(pol_rows);
all_rows = true(1, polar_profile.n_row);
[bin_range, wedge_ok] = polar_profile.widest_common(all_rows, 'divergence');
r_um = polar_profile.bin_center_um;
span_line = sprintf('%d of %d wedges, %.1f-%.1f \\mum from the wall', nnz(wedge_ok), polar_profile.n_wedge, ...
    r_um(bin_range(1)), r_um(bin_range(2)));

% one point per row : its wall-normal displacement averaged over that common span
disp_curves = polar_profile.curve_rows(all_rows, bin_range, quantity = 'disp_radial');
span_bins = bin_range(1):bin_range(2);
mean_disp = mean(disp_curves(:, span_bins), 2).';
% the sign check : dilation has to push the field outward, constriction pull it inward
[n_follow, n_signed] = event_validatesign(pol_rows, mean_disp);
disp_fig = make_fig('PIV polar  mean wall displacement against dD', 'normal');
disp_fig.update_figsize([6 4]);
hold(disp_fig.ax, 'on');
scatter(disp_fig.ax, dD_um, mean_disp, 50, row_rgb, 'filled');
xline(disp_fig.ax, 0, 'k--');
yline(disp_fig.ax, 0, 'k--');
disp_fig.put_xaxistitle('\DeltaD (\mum)');
disp_label = strrep(showpiv.axis_label('disp_radial'), '_', '\_');   % the names are TeX, their underscores are not
disp_fig.put_yaxistitle(disp_label);
title(disp_fig.ax, {sprintf('%d of %d rows move the way their polarity says', n_follow, n_signed), span_line});

for quantity = {'divergence', 'strain_radial'}
    curves = polar_profile.curve_rows(all_rows, bin_range, quantity = quantity{1});
    profile_fig = make_fig(['PIV polar  ' quantity{1} ' against distance from the wall'], 'normal');
    profile_fig.update_figsize([6 4]);
    hold(profile_fig.ax, 'on');
    for k = 1:polar_profile.n_row
        plot(profile_fig.ax, r_um, curves(k, :), '-o', 'Color', row_rgb(k, :), 'MarkerSize', 3);
    end
    yline(profile_fig.ax, 0, 'k--');
    profile_fig.put_xaxistitle('distance from the wall (\mum)');
    quantity_label = strrep(showpiv.axis_label(quantity{1}), '_', '\_');
    profile_fig.put_yaxistitle(quantity_label);
    title(profile_fig.ax, {'red = dilation, blue = constriction', span_line});
end

%% Opt. One row on screen : its displacement field over the FROM frame, and the frames PIV used
R = pivtransition(1);
keep = ~isnan(R.xyuv(:, :, 3)) & ~isnan(R.xyuv(:, :, 4));
keep(2:2:end, :) = false;                             % drawing only : every other grid row and column
keep(:, 2:2:end) = false;
uv_map = piv_stamp(R.xyuv, R.planes.imsize, keep);
P = showpiv(sprintf('ev%d %s %s', R.id, R.state, R.pol), 10);
P.plot_background(stack(:, :, R.from));
P.plot_quiver(uv_map);
inter = piv_interleave(stack, R.from, R.to, halfwin);    % odd slice = FROM, even = TO
figure('Name', sprintf('ev%d blink  odd = FROM, even = TO', R.id));
sliceViewer(inter);
