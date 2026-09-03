%FIG_FWHM_EVENTRASTER  Lumen diameter around each scored event, one row per vessel.
%   A check on the DETECTION, not a result : every row should turn at x = 0. A row of
%   one colour caught no excursion; one that turns early or late is a window off the
%   event. Colour is dbv, diameter against the event's own pre-event level, computed
%   by tablegeneration_fwhmrelation. Nothing saved unless block 4 is uncommented.

%% 1. Read
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');
dirs_ecspath;
clee = color_lee();

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
param.fontsize = 11;
% name and label side by side, matched by position : the table stores an identifier, an axis names the thing
param.state       = ["an_trans",      "na_trans",      "nr_trans",    "ra_trans",     "uarousal"];
param.state_label = ["awake to NREM", "NREM to awake", "NREM to REM", "REM to awake", "microarousal"];
param.clim = 1.5;                  % um   colour saturates here, symmetric about 0
param.sample_rate = 3;             % Hz   from param.fs of the state analysis

% the ordering window matches how long the response lasts; a whole-window mean
% buries a short one. Both are the events' own scale, see FINDINGS.md
param.window_bout_s = 10;          % s    about the median microarousal
param.window_trans_s = 25;         % s    the transition half window

dirs = dirs_central();
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);
relation = load(fullfile(dirs.save_dir, 'fwhmrelation.mat')).save_content;
state = relation.state;
vessel = unique(state.vessel_key);
n_vessel = numel(vessel);
fprintf('%d state rows over %d vessels\n', height(state), n_vessel);

%% 2. One row order for every panel
% taken from ONE state, so a vessel keeps its row across the five panels
param.order_pick = 1;              % which of param.state the row order comes from

order_state = param.state(param.order_pick);
order_label = param.state_label(param.order_pick);
order_trace = stack_by_vessel(state(state.state == order_state, :), vessel);
order_t = time_axis(size(order_trace, 2), param.sample_rate);
on_window = order_t > 0 & order_t <= param.window_trans_s;
order_amplitude = mean(order_trace(:, on_window), 2, 'omitnan');
% sort leaves NaN last, so a vessel with no such event sits at the bottom of every one
[~, row_order] = sort(order_amplitude);
fprintf('ordered on %s : %d of %d vessels have it, %+.2f .. %+.2f um\n', ...
    order_label, sum(isfinite(order_amplitude)), n_vessel, ...
    min(order_amplitude), max(order_amplitude));

%% 3. Draw, one figure per state
% not filtered by any verdict : a window that caught nothing is what is looked for.
% A vessel absent from a state stays as a blank row, so the rows do not shift
for k = 1:numel(param.state)
    this_state = param.state(k);
    this_label = param.state_label(k);
    trace = stack_by_vessel(state(state.state == this_state, :), vessel);
    if all(isnan(trace), 'all')
        fprintf('%-14s no events\n', this_label);
        continue
    end
    trace = trace(row_order, :);
    t_axis = time_axis(size(trace, 2), param.sample_rate);
    if endsWith(this_state, "_trans")
        window_s = param.window_trans_s;
    else
        window_s = param.window_bout_s;
    end
    after = mean(trace(:, t_axis > 0 & t_axis <= window_s), 2, 'omitnan');

    fig_raster = figure('Color', 'w', 'Position', [80 + 30*k, 60, 560, 520], ...
        'Name', 'fig_fwhm_eventraster_' + this_state);
    ax = axes(fig_raster); %#ok<LAXES>
    imagesc(ax, t_axis, 1:n_vessel, trace, 'AlphaData', ~isnan(trace));
    clim(ax, [-param.clim param.clim])
    colormap(ax, clee.gradient.hilo)   % blue-white-red, made for reading against zero
    xline(ax, 0, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
    bar_handle = colorbar(ax);
    bar_handle.Label.String = 'lumen diameter, um from the pre-event level';
    xlabel(ax, 'time from the event, s')
    ylabel(ax, "vessel, ordered by the " + order_label + " excursion")
    set(ax, 'FontSize', param.fontsize, 'Box', 'on', 'YDir', 'normal', 'Color', 'w')
    title(ax, sprintf("%s  |  %d of %d vessels have it", this_label, ...
        sum(isfinite(after)), n_vessel), 'FontWeight', 'normal')
    fprintf('%-14s %2d vessels | first %2d s after %+.2f .. %+.2f um, median %+.2f\n', ...
        this_label, sum(isfinite(after)), window_s, ...
        min(after, [], 'omitnan'), max(after, [], 'omitnan'), median(after, 'omitnan'));
end

%% 4. Save, when a panel above is worth keeping
% save_figure(fig_raster, dirs.save_dir, "fig_fwhm_eventraster_" + this_state)

%% ---------------------------------------------------------------- helpers
function trace = stack_by_vessel(state_rows, vessel)
%STACK_BY_VESSEL  Each vessel's event-triggered trace, sessions averaged within it.
%   A mean of per-session means, weighted by neither event count nor session count.
%   IN   state_rows  table     the rows of ONE state
%        vessel      Vx1 str   every vessel in the set, in the order rows come out
%   OUT  trace       V x T double   NaN in the row of a vessel with no such event
    trace = nan(numel(vessel), 1);
    if isempty(state_rows)
        return
    end
    stacked = cell2mat(state_rows.event_dbv);
    [~, vessel_of] = ismember(state_rows.vessel_key, vessel);
    trace = nan(numel(vessel), size(stacked, 2));
    for v = 1:numel(vessel)
        on_vessel = vessel_of == v;
        if ~any(on_vessel)
            continue
        end
        trace(v, :) = mean(stacked(on_vessel, :), 1, 'omitnan');
    end
end

function t_axis = time_axis(n_sample, sample_rate)
%TIME_AXIS  Seconds either side of the event, read off the trace's own length.
    half = (n_sample - 1) / 2;
    t_axis = (-half:half) / sample_rate;
end

function save_figure(fig, save_dir, fig_name) %#ok<DEFNU>
% R2023b exportgraphics cannot write svg, so the vector copy goes through print.
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    print(fig, fullfile(save_dir, fig_name + ".svg"), '-dsvg', '-vector');
    exportgraphics(fig, fullfile(save_dir, fig_name + ".png"), 'Resolution', 200);
    fprintf('wrote %s\n', fullfile(save_dir, fig_name));
end
