sessiondir = 'G:\tmp\00_igkl\hql090\251016_hql090_sleep\HQL090_sleep251016_006';
session = ECSSession(sessiondir);
session = session.load_primary_results();
%%
session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
%% Two photon analysis time axis calculation
session.pax_fwhm.t_axis = linspace(session.img_param.imgstarttime,session.img_param.imgendtime,numel(session.pax_fwhm.idx.upperBVboundary));
sleep_score = load(fullfile(sessiondir,"peripheral","sleep_score.mat"));
% Calculate
% Find indices for all REM epochs
% Calculate State Analysis (Indices, Bouts, Sliced Data)
rem.pax_fwhm = add_state_analysis(session.pax_fwhm, sleep_score.REMTimes);
%%
nrem.pax_fwhm = add_state_analysis(session.pax_fwhm, sleep_score.NREMTimes, sleep_score.uArousalTimes);

%%
awake.pax_fwhm = add_state_analysis(session.pax_fwhm, sleep_score.AwakeTimes);
drowsy.pax_fwhm = add_state_analysis(session.pax_fwhm, sleep_score.DrowsyTimes);

% Image frame location (kept separate for now as it uses img_param.taxis)
rem.img.stateidx = statebin2frame(sleep_score.REMTimes, session.img_param.taxis);
nrem.img.stateidx = statebin2frame(sleep_score.NREMTimes, sleep_score.uArousalTimes, session.img_param.taxis);
awake.img.stateidx = statebin2frame(sleep_score.AwakeTimes, session.img_param.taxis);
drowsy.img.stateidx = statebin2frame(sleep_score.DrowsyTimes, session.img_param.taxis);

% bouts loc of image
nrem.img.bouts = stateidx2bouts(nrem.img.stateidx);
awake.img.bouts = stateidx2bouts(awake.img.stateidx);
rem.img.bouts = stateidx2bouts(rem.img.stateidx);
drowsy.img.bouts = stateidx2bouts(drowsy.img.stateidx);

%% Spectrogram Calculation
% Using the bouts from the structured analysis
rem.thickness.bv_spectrogramlist = get_spectboutarray(rem.pax_fwhm.bouts, session.pax_fwhm.thickness.bvchanges, session.pax_fwhm.param.fs, 0.1);
awake.thickness.bv_spectrogramlist = get_spectboutarray(awake.pax_fwhm.bouts, session.pax_fwhm.thickness.bvchanges, session.pax_fwhm.param.fs, 0.1);
nrem.thickness.bv_spectrogramlist = get_spectboutarray(nrem.pax_fwhm.bouts, session.pax_fwhm.thickness.bvchanges, session.pax_fwhm.param.fs, 0.1);
drowsy.thickness.bv_spectrogramlist = get_spectboutarray(drowsy.pax_fwhm.bouts, session.pax_fwhm.thickness.bvchanges, session.pax_fwhm.param.fs, 0.1);
%%
% Average Spectrograms (handle varying frequency axes)
nrem.thickness.bv_spectral_summary = get_summaryspectrogram(nrem.thickness.bv_spectrogramlist,session.pax_fwhm.param.fs);
rem.thickness.bv_spectral_summary = get_summaryspectrogram(rem.thickness.bv_spectrogramlist,session.pax_fwhm.param.fs);
awake.thickness.bv_spectral_summary = get_summaryspectrogram(awake.thickness.bv_spectrogramlist,session.pax_fwhm.param.fs);
drowsy.thickness.bv_spectral_summary = get_summaryspectrogram(drowsy.thickness.bv_spectrogramlist,session.pax_fwhm.param.fs);

%%
figure()
hold on
% Plot individual (log-log)
target = awake.thickness;
for i = 1:numel(target.bv_spectrogramlist)
    loglog(target.bv_spectrogramlist(i).F,target.bv_spectrogramlist(i).S, 'Color', [0.8 0.6 0.6]) % light grey
end
% Plot average
loglog(target.bv_spectral_summary.F, target.bv_spectral_summary.S, 'r', 'LineWidth', 2)
xlabel('Frequency (Hz)'); ylabel('Power');
title('Average spectral density');

% Plot individual (log-log)

target = nrem.thickness;
for i = 1:numel(target.bv_spectrogramlist)
    loglog(target.bv_spectrogramlist(i).F,target.bv_spectrogramlist(i).S, 'Color', [0.6 0.8 0.6]) % light grey
end
% Plot average
loglog(target.bv_spectral_summary.F, target.bv_spectral_summary.S, 'g', 'LineWidth', 2)
set(gca, 'XScale', 'linear', 'YScale', 'linear')
xlabel('Frequency (Hz)'); ylabel('Power');
title('Average spectral density');

%% mtspectrum based analysis of fwhm thickness

%% spectral analysis

focus = drowsy.fwhmloc_bouts;
cnt = 1;
for i = 1:numel(focus)
    loc = focus{i};
    if numel(loc)>eventlen_thr
        disp(i)
        % Preprocessing: Sgolay filter + Detrend
        bv_signal = session.pax_fwhm.thickness.bvchanges(loc);
        % bv_signal = detrend(sgolayfilt(bv_signal, 3, 5));
        bv_spectrogram = get_spectrogram(session.pax_fwhm.param.fs,bv_signal);
        bv_spectrogram.boutidx = i;
        bv_spectrogramlist(cnt) = bv_spectrogram;
        cnt = cnt + 1;
    end
end


%%
figure(name='nrem_pvs')
hold on
clc
for i = 1:numel(focus)
    loc = focus{i};
    if numel(loc)>eventlen_thr
        disp(i)
        % Preprocessing: Sgolay filter + Detrend
        pvs_signal = session.pax_fwhm.thickness.pvschanges_total(loc);
        % bv_signal = detrend(sgolayfilt(bv_signal, 3, 5));
        bv_spectrogram = get_spectrogram(session.pax_fwhm.param.fs,pvs_signal);
        loglog(bv_spectrogram.F,bv_spectrogram.S)
    end
end

%%

% Calculate PSD for each NREM bout
nrem_psd = struct('bv', {}, 'pvs', {});
fs = session.pax_fwhm.param.fs; % Calculate fs from actual time axis

for i = 1:numel(fwhmloc.nrem_bouts)
    loc = fwhmloc.nrem_bouts{i};

    % Preprocessing: Sgolay filter + Detrend
    bv_signal = session.pax_fwhm.thickness.bvchanges(loc);
    bv_signal = detrend(sgolayfilt(bv_signal, 3, 11));

    pvs_signal = session.pax_fwhm.thickness.pvschanges_dynamic(loc);
    pvs_signal = detrend(sgolayfilt(pvs_signal, 3, 11));

    % BV power spectrum
    [pxx_bv, f] = pwelch(bv_signal-mean(bv_signal), [], [], [], fs);
    nrem_psd(i).bv.Pxx = pxx_bv;
    nrem_psd(i).bv.F = f;

    % PVS power spectrum
    [pxx_pvs, f] = pwelch(pvs_signal-mean(pvs_signal), [], [], [], fs);
    nrem_psd(i).pvs.Pxx = pxx_pvs;
    nrem_psd(i).pvs.F = f;
end


%%

fig = figure('Name','Thickness plot');
t = tiledlayout(fig,3,1);
loc = fwhmloc.nrem;
ax.idx = nexttile(t, 1);
hold on
imagesc(ax.idx, t_axis, 1:size(session.pax_fwhm.kymograph.kgph_pvs_processed,1), session.pax_fwhm.kymograph.kgph_pvs_processed)
plot(ax.idx,t_axis,session.pax_fwhm.idx.clean_lowerBVboundary,'r')
plot(ax.idx,t_axis,session.pax_fwhm.idx.clean_pvsdownedge_idx,'b')
plot(ax.idx,t_axis,session.pax_fwhm.idx.clean_upperBVboundary,'r')
plot(ax.idx,t_axis,session.pax_fwhm.idx.clean_pvsupedge_idx,'b')


ax.upidx = nexttile(t, 3);
plot(ax.upidx,t_axis(loc),session.pax_fwhm.idx.clean_upperBVboundary(loc),'r')
hold on
plot(ax.upidx,t_axis(loc),session.pax_fwhm.idx.clean_pvsdownedge_idx(loc),'b')

ax.thickness = nexttile(t, 2);
plot(ax.thickness,t_axis(loc),session.pax_fwhm.thickness.bvchanges(loc),'r')
hold on
plot(ax.thickness,t_axis(loc),session.pax_fwhm.thickness.pvschanges_dynamic(loc),'b')

all_axes = struct2cell(ax);
all_axes = [all_axes{:}]; % Convert cell array to object array
linkaxes(all_axes, 'x');



%%

channel = ["bv","csf"];
state = ["rem","nrem","awake"];

%%
fig = figure();
fig.Units = 'inches';
fig.Position = [21 1 18 7]; % [left bottom width height]
t = tiledlayout(numel(channel),numel(state), 'TileSpacing', 'none', 'Padding', 'none');

img_handles = gobjects(numel(channel),numel(state));
%
ax = struct();
for ch =1:numel(channel)
    for st = 1:numel(state)
        fieldname =  strcat(channel(ch),"_",state(st));
        % Use nexttile to get the next tile in the layout
        ax.(fieldname) = nexttile(t);
        current_locs = imgloc.(state(st));
        if isempty(current_locs)
            title(ax.(fieldname), [fieldname " (Empty)"], 'Interpreter', 'none');
            continue;
        end

        if ch == 1
            img_data = mean(session.stackch1(:,:,current_locs), 3);
        elseif ch == 2
            img_data = mean(session.stackch2(:,:,current_locs), 3);
        end

        % Contrast adjustment: 10th to 90th percentile
        clim_val = prctile(img_data(:), [10 99]);
        img_handles(ch, st) = imagesc(ax.(fieldname), img_data, clim_val);
        axis(ax.(fieldname), 'image'); % Fix aspect ratio
        title(ax.(fieldname), fieldname, 'Interpreter', 'none');
    end
end
%% Sync all axes for simultaneous panning/zooming
all_axes = struct2cell(ax);
all_axes = [all_axes{:}]; % Convert cell array to object array
linkaxes(all_axes, 'xy');
%%
%% Interactive Marker Control (Click on Figure)
fprintf('--- Interactive Marker Mode ---\n');
fprintf('Click on any image to mark position on all axes.\n');
fprintf('Press Enter (key) to stop.\n');

marker_handles = [];
while true
    try
        [x_val, y_val, button] = ginput(1);
    catch
        break; % Handle figure closure
    end

    % Break if Enter (empty button) or any key press (button > 3 typically)
    % Mouse clicks are 1 (Left), 2 (Middle), 3 (Right)
    if isempty(button) || ~ismember(button, [1, 2, 3])
        fprintf('Exiting marker mode.\n');
        break;
    end

    % Remove old markers
    delete(marker_handles);
    marker_handles = gobjects(0);

    % Add new markers to all axes
    for i = 1:numel(all_axes)
        try
            curr_ax = all_axes(i);
            hold(curr_ax, 'on');
            h = plot(curr_ax, x_val, y_val, 'rx', 'MarkerSize', 8, 'LineWidth', 1);
            marker_handles(end+1) = h;
        catch
            % Handle closed axes
        end
    end
    drawnow;
end

%%

imagesc(ans)
axis image
%%




%%
session.pax_fwhm.kymograph.kgph_lumen_processed(:,drowsy_locs)


%%
clee = color_lee;
%%
thickness_fig = make_fig('Thickness figure');
%%
thickness_fig.hold_axis(true)
thickness_fig.resolution = session.img_param.pixel2um;
thickness_fig.update_figsize([15,3])
thickness_fig.monitor_xyinch = [21, 5];
thickness_fig.reset_axis
thickness_fig.loc.x = t_axis;
xlim(thickness_fig.ax,[0 max(t_axis)])
thickness_fig.plot_line(session.pax_fwhm.thickness.bvchanges,'r')
thickness_fig.plot_line(session.pax_fwhm.thickness.pvschanges_total, clee.clist.darkgreen)
thickness_fig.put_xaxistitle('Time (sec)')
thickness_fig.put_yaxistitle('Thickness change (um)')
%%
kymograph_fig = make_fig('Kymograph figure');
kymograph_fig.resolution = session.img_param.pixel2um;
kymograph_fig.update_figsize([15,3])
kymograph_fig.monitor_xyinch = [21, 5];
kymograph_fig.reset_axis
kymograph_fig.loc.x = t_axis;

kymograph_fig.plot_kymograph(session.pax_fwhm.kymograph.kgph_pvs_processed(:,fwhmloc.rem),t_axis(fwhmloc.rem))
xlim(kymograph_fig.ax,[min(t_axis(fwhmloc.rem)) max(t_axis(fwhmloc.rem))])


%%
scope = fwhmloc.nrem;
kymograph_fig = make_fig('Kymograph figure');
kymograph_fig.resolution = session.img_param.pixel2um;
kymograph_fig.update_figsize([15,3])
kymograph_fig.monitor_xyinch = [21, 5];
kymograph_fig.reset_axis
kymograph_fig.loc.x = t_axis;

kymograph_fig.plot_kymograph(session.pax_fwhm.kymograph.kgph_pvs_processed(:,scope),t_axis(scope))
xlim(kymograph_fig.ax,[min(t_axis(scope)) max(t_axis(scope))])







