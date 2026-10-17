%B_REVIEWSESSION  One session folder, whatever it holds : the FWHM and cluster figures, and the traces
%   stacked on one time axis over the sleep states or the puff x locomotion bouts.
%   A FWHM or cluster figure whose data the folder lacks is skipped with one line; every trace
%   panel is drawn, empty where there is nothing to draw. Comment a panel's block out to drop it.

% Path setup
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

% Directory setup
% sessiondir = 'G:\tmp\01_igkltdt\hql104\260607_hql104_sleep\HQL104_sleep260607_006';
sessiondir = 'K:\tmp\02_egfptdt\hql108\261002_hql108_whiskerb\hql108_whiskerb261002_007';


% 1. Load : the primary results, the peripheral products and the tissue motion
session = ECSSession(sessiondir);
session = session.load_primary_results;
session = session.load_peripheral;
session = session.load_motion;
%%

%% 2. FWHM panels
if ~isempty(session.pax_fwhm)
    pax_fig = analysis_pax_makefig(session.pax_fwhm, session.pax_fwhm.t_axis,...
        session.img_param.pixel2um, session.dir_struct.figures_fwhm);
end

%% 3. the cluster panels. No stack is read: load_primary_results already brought
% polarcluster.mat in, and t_axis and pixel2um stand in for twophoton_processed the
% same way the pax call above does it
if ~isempty(session.polarcluster) && ~isempty(session.pax_fwhm)
    analysis_cluster_makefig(session.polarcluster, session.roilist, session.pax_fwhm,...
        session.pax_fwhm.t_axis, session.img_param.pixel2um,...
        session.dir_struct.figures_polarcluster);
end

%% 4. the polar profiles. Only for a session that went through
% analysis_clusterpolar.transform -- this reads pax_angle, polar_theta,
% polar_profiles and polar_binstart_deg, and a polarcluster.mat saved before that
% step carries none of the four
if ~isempty(session.polarcluster)
    analysis_polar_makefig(session.polarcluster, session.dir_struct.figures_polarcluster);
end

%% 5. the awake bouts with the FWHM, for a whisker session that has both
if isempty(session.sleep_score) && ~isempty(session.primary_analog) && ~isempty(session.pax_fwhm)
    air_on = session.primary_analog.rawdata.raw_Air_puff1 > 0;
    awake_bout = bout_puffmove(air_on, session.primary_analog.samplingfrequency, ...
        session.primary_analog.ball.resampled_ball, session.primary_analog.ball.rs_taxis, ...
        session.primary_analog.ball.ds_fps);
    awake_fig = analysis_awake_makefig(awake_bout, session.pax_fwhm.t_axis, ...
        session.pax_fwhm.thickness.bv, session.pax_fwhm.thickness.totalpvs, ...
        session.img_param.pixel2um, session.primary_analog.ball.rs_taxis, ...
        session.primary_analog.ball.resampled_ball, session.dir_struct.figures_fwhm);
end

%% 6. the traces on one time axis : one window per panel, top to bottom, empty where the session has nothing
% 6.1 what goes behind them : the scored states, or the puff x locomotion bouts off the air channel
has_analog = ~isempty(session.primary_analog);
air_on = [];
if has_analog && isfield(session.primary_analog.rawdata, 'raw_Air_puff1')
    air_on = session.primary_analog.rawdata.raw_Air_puff1 > 0;
end
awake_bout = struct();
if ~isempty(session.sleep_score)
    bout_kind = "sleep";
elseif any(air_on)
    bout_kind = "awake";
    awake_bout = bout_puffmove(air_on, session.primary_analog.samplingfrequency, ...
        session.primary_analog.ball.resampled_ball, session.primary_analog.ball.rs_taxis, ...
        session.primary_analog.ball.ds_fps);
else
    bout_kind = "none";
end
if ~has_analog
    fprintf('no analysis_analog.mat : ball, EMG, ECoG and the puff bouts need peripheral_batch first\n');
end
trace_fig = struct();

% 6.2 state : the bout strip
trace_fig.state = state_window(bout_kind, session.sleep_score, awake_bout);

% 6.3 ball : the signed velocity as saved, backward rolls below zero
t = [];
y = [];
if has_analog && isfield(session.primary_analog.ball, 'resampled_ball')
    t = session.primary_analog.ball.rs_taxis(:);
    y = session.primary_analog.ball.resampled_ball(:);
end
trace_fig.ball = trace_window("ball", t, y, 'Ball velocity (cm/s)');

% 6.4 rostral : tissue displacement toward the image top
t = [];
y = [];
if ~isempty(session.motion)
    t = session.motion.t_s;
    y = session.motion.rostral_um;
end
trace_fig.rostral = trace_window("rostral", t, y, 'rostral (\mum)');

% 6.5 lateral : tissue displacement toward the image right
t = [];
y = [];
if ~isempty(session.motion)
    t = session.motion.t_s;
    y = session.motion.lateral_um;
end
trace_fig.lateral = trace_window("lateral", t, y, 'lateral (\mum)');

% 6.6 bv : the lumen FWHM
t = [];
y = [];
if ~isempty(session.pax_fwhm)
    t = session.pax_fwhm.t_axis(:);
    y = session.pax_fwhm.thickness.bv(:) * session.img_param.pixel2um;
end
trace_fig.bv = trace_window("bv", t, y, 'BV (\mum)');

% 6.7 pvs : both perivascular sides
t = [];
y = [];
if ~isempty(session.pax_fwhm)
    t = session.pax_fwhm.t_axis(:);
    y = session.pax_fwhm.thickness.totalpvs(:) * session.img_param.pixel2um;
end
trace_fig.pvs = trace_window("pvs", t, y, 'PVS, both (\mum)');

% 6.8 pupil
t = [];
y = [];
if ~isempty(session.camera_session) && isfield(session.camera_session.eye, 'diameter')
    t = session.camera_session.eye.t_axis(:);
    y = session.camera_session.eye.diameter(:);
end
trace_fig.pupil = trace_window("pupil", t, y, 'Pupil (px)');

% 6.9 emg : resampled power
t = [];
y = [];
if has_analog && isfield(session.primary_analog.emg, 'resampled_power')
    t = session.primary_analog.emg.rs_taxis(:);
    y = session.primary_analog.emg.resampled_power(:);
end
trace_fig.emg = trace_window("emg", t, y, 'EMG power');

% 6.10 ecog : the spectrogram
ecog_spectrum = [];
if has_analog && isfield(session.primary_analog.ecog, 'ecogspectrum')
    ecog_spectrum = session.primary_analog.ecog.ecogspectrum;
end
trace_fig.ecog = spectrum_window(ecog_spectrum);

%% 6.11 one time range over every window, and the bouts behind the traces
panel_names = string(fieldnames(trace_fig))';
t_end = 0;
for name = panel_names
    t_end = max(t_end, data_end(trace_fig.(name).ax));
end
for name = panel_names
    ax = trace_fig.(name).ax;
    if t_end > 0
        xlim(ax, [0 t_end]);
    end
    if ismember(name, ["ecog", "state"])
        continue
    end
    switch bout_kind
        case "sleep"
            sleep_patch = plot_sleep_patches(ax, session.sleep_score);
            uistack(sleep_patch, 'bottom');
        case "awake"
            plot_awake_patches(ax, awake_bout);
    end
end

%% 6.12 the board, in the order the windows were made
[~, session_name] = fileparts(sessiondir);
title(trace_fig.(panel_names(1)).ax, session_name, 'Interpreter', 'none', 'FontWeight', 'normal');
% the strip takes one grid row, a trace two
row_span = 2 * ones(size(panel_names));
row_span(panel_names == "state") = 1;
row_start = cumsum([1, row_span(1:end-1)]);
board = fusion_figure('grid', [sum(row_span) 1], 'width', 9, ...
    'height', 0.8 * sum(row_span), 'fontsize', 9);
for k = 1:numel(panel_names)
    board.add(trace_fig.(panel_names(k)).ax, 'at', [row_start(k) 1], 'span', [row_span(k) 1]);
end
board.link_x(1);

%% ---------------------------------------------------------------- helpers
function f = trace_window(name, t, y, y_label)
%TRACE_WINDOW  One make_fig window holding one trace against time.
%   Caller: b_reviewsession
% IN   name     1x1 str        window name
%      t        N x 1 float s  times of the trace
%      y        N x 1 float    the trace
%      y_label  1 x n char
% OUT  f        make_fig
    f = make_fig(char("review_" + name));
    f.update_figsize([9 2]);
    f.reset_axis();
    hold(f.ax, 'on');
    plot(f.ax, t, y, '-', 'Color', color_lee.clist.black, 'LineWidth', 0.7);
    f.put_yaxistitle(y_label);
    f.put_xaxistitle('Time (sec)');
end

function f = state_window(bout_kind, sleep_score, awake_bout)
%STATE_WINDOW  The bouts as a strip of their own : the sleep states in one row, or puff and locomotion / rest in two.
%   Caller: b_reviewsession
% IN   bout_kind    1x1 str      "sleep" | "awake"
%      sleep_score  1x1 struct   sleep_score.mat contents; read when bout_kind is "sleep"
%      awake_bout   1x1 struct   bout_puffmove's time_table; read when bout_kind is "awake"
% OUT  f            make_fig
    f = make_fig('review_state');
    f.update_figsize([9 1.2]);
    f.reset_axis();
    switch bout_kind
        case "sleep"
            plot_sleep_patches(f.ax, sleep_score, 'YLim', [0 1]);
            ylim(f.ax, [0 1]);
            set(f.ax, 'YTick', []);
        case "awake"
            patch_list = plot_awake_patches(f.ax, awake_bout, rows = true);
            ylim(f.ax, [0.4 2.6]);
            set(f.ax, 'YTick', 1:2, 'YTickLabel', {'locomotion / rest', 'puff'});
            legend(f.ax, patch_list, 'Orientation', 'horizontal', 'Location', 'northoutside', ...
                'Box', 'off', 'Interpreter', 'none');
    end
    f.put_xaxistitle('Time (sec)');
end

function f = spectrum_window(spectrum)
%SPECTRUM_WINDOW  The ECoG spectrogram in one make_fig window, the way main_analog draws it.
%   Caller: b_reviewsession
% IN   spectrum  1x1 struct | []   analysis_analog ecog.ecogspectrum : T, F and baseline_S
%                                   (sleep-scored) or log_norm_spectrum, frequency x time;
%                                   [] = an empty window
% OUT  f         make_fig
    f = make_fig('review_ecog');
    f.update_figsize([9 2]);
    f.reset_axis();
    f.put_yaxistitle('Freq (Hz)');
    f.put_xaxistitle('Time (sec)');
    if isempty(spectrum)
        return
    end
    power = spectrum.log_norm_spectrum;
    if isfield(spectrum, 'baseline_S')
        power = spectrum.baseline_S;
    end
    surface(f.ax, spectrum.T, spectrum.F, zeros(size(power)), power, 'LineStyle', 'none');
    set(f.ax, 'YScale', 'log');
    view(f.ax, 2);
end

function t_end = data_end(ax)
%DATA_END  The largest x any line, patch or surface on the axes reaches; 0 when it holds none.
%   Caller: b_reviewsession
% IN   ax     1x1 axes
% OUT  t_end  1x1 float s
    t_end = 0;
    for child = findobj(ax, '-property', 'XData')'
        x_data = child.XData(:);
        if ~isempty(x_data)
            t_end = max(t_end, max(x_data));
        end
    end
end
