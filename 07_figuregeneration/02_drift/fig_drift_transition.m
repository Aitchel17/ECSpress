%FIG_DRIFT_TRANSITION  Drift trace around each state transition, one line per vessel,
%   both components. The drift counterpart of transition_fig. Reads
%   drift_transition.mat and draws one figure per normalisation: rows are the two
%   components, columns the four transitions.
%
%   What differs from transition_fig, all settled in DRIFT_FIGURE_PLAN.md: two
%   components instead of three series, no awake-normalised arm, colour by the state
%   being entered, y matched across every tile, and the bold line and its band are
%   both over vessels -- transition_fig drew a vessel-level band around the
%   mouse-level mean.
clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
addpath('g:/03_program/01_ecspress/09_dirstruct');    % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
clee = color_lee();
param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
dirs = dirs_central();
exp_path = fullfile(dirs.secondary_root, param.dataset);
table_path = fullfile(exp_path, "drift_transition.mat");
load_struct = load(table_path);
transition_table = load_struct.save_content;

% The stored trace is center +/- transition_window/2 and the window is not in the
% table; 25 is state_integration.param.transition_window. see CLAUDE_LOG.md
param.transition_window = 25;      % s
t_axis.dpoints = size(transition_table.abs_numeric.filtered_table.data{1});
t_axis.tpoints = linspace(-param.transition_window/2, param.transition_window/2, t_axis.dpoints(2));

%% Filters : which rows the figure is about
filt.component_drawn = ["drift_row", "drift_col"];   % one row of tiles each (str)
filt.transition_drawn = ["an_trans", "na_trans", "nr_trans", "ra_trans"];   % one column each

% Keyed by the name, not by position, so dropping one from filt cannot slide the
% rest along by one. +row is image-down, which is caudal, and the axis reads
% rostral-positive, so the sign flips here and nowhere else.
param.component_label = struct('drift_row', "rostral − caudal (um)", ...
    'drift_col', "lateral − medial (um)");
param.component_sign = struct('drift_row', -1, 'drift_col', 1);
label_map = struct( ...
    'an_trans', struct('pre',"Awake",'post',"NREM",  'enters', "nrem"), ...
    'na_trans', struct('pre',"NREM", 'post',"Awake", 'enters', "awake"), ...
    'nr_trans', struct('pre',"NREM", 'post',"REM",   'enters', "rem"), ...
    'ra_trans', struct('pre',"REM",  'post',"Awake", 'enters', "awake"));

param.norm_drawn = ["abs_data", "awakesubtract_data"];
param.norm_label = struct('abs_data', "absolute", 'awakesubtract_data', "awake-subtracted");

param.fontsize        = 14;
param.tile_widthinch  = 2.2;
param.tile_heightinch = 2.4;

%% One figure per normalisation, rows = components, columns = transitions
n_component = numel(filt.component_drawn);
n_transition = numel(filt.transition_drawn);
for ni = 1:numel(param.norm_drawn)
    norm_field = param.norm_drawn(ni);
    vessel_ave = transition_table.(norm_field).Date_ave;
    fig = figure("Name", "Drift around transitions – " + param.norm_label.(norm_field));
    monitor_xyinch = [10 2];
    xy_sizeinch = [param.tile_widthinch*n_transition param.tile_heightinch*n_component];
    set(fig, 'Units','inches', ...
        "Position",[monitor_xyinch(1) monitor_xyinch(2) xy_sizeinch(1) xy_sizeinch(2)])
    tl = tiledlayout(fig, n_component, n_transition, "TileSpacing", "compact");

    all_axes = [];
    for ci = 1:n_component
        component_name = filt.component_drawn(ci);
        for ti = 1:n_transition
            transition_name = filt.transition_drawn(ti);
            ax = nexttile(tl);
            on_tile = vessel_ave.state_name == transition_name & vessel_ave.DataType == component_name;
            trace_mat = param.component_sign.(component_name) ...
                * cell2mat(vessel_ave.data(on_tile))';   % T x n_vessel
            labels = label_map.(transition_name);
            plot_tile(ax, t_axis.tpoints, trace_mat, state_color(clee, labels.enters), labels, param.fontsize);
            if ti == 1
                ylabel(ax, param.component_label.(component_name))
            end
            all_axes = [all_axes, ax]; %#ok<AGROW>
        end
    end
    xlabel(tl, "time from transition (s)")
    match_ylim(all_axes);
end

%% ── Local functions ────────────────────────────────────────────────────

function cfg = state_color(clee, state_name)
    cfg.bold  = clee.clist.(state_name);
    cfg.faint = 1 - 0.55 * (1 - cfg.bold);
end

function plot_tile(ax, taxis, trace_mat, cfg, labels, fontsize)
    % one line per vessel, the mean over vessels, and its 95 % band -- same rung
    n_vessel = size(trace_mat, 2);
    mean_trace = mean(trace_mat, 2, 'omitnan');
    t_val = tinv(0.975, n_vessel - 1);
    ci95 = t_val * std(trace_mat, 0, 2, 'omitnan') / sqrt(n_vessel);

    cla(ax)
    hold(ax, "on")
    % the band is an OPAQUE light grey drawn FIRST : print -dsvg cannot express
    % per-object alpha, so nothing in this axes is transparent
    patch(ax, [taxis, fliplr(taxis)], [(mean_trace + ci95)', fliplr((mean_trace - ci95)')], ...
        [0.9 0.9 0.9], 'EdgeColor', 'none');
    plot(ax, taxis, trace_mat, "Color", cfg.faint);
    plot(ax, taxis, mean_trace, "Color", cfg.bold, "LineWidth", 2);
    xlim(ax, [taxis(1) taxis(end)])
    xline(ax, 0, "--")
    yline(ax, 0, '--k', 'Alpha', 0.3, 'LineWidth', 1);
    text(ax, 0.48, 0.98, labels.pre,  'Units','normalized', ...
        "HorizontalAlignment", "right", "VerticalAlignment", "top", "FontSize", fontsize-2)
    text(ax, 0.52, 0.98, labels.post, 'Units','normalized', ...
        "HorizontalAlignment", "left",  "VerticalAlignment", "top", "FontSize", fontsize-2)
    text(ax, 0.98, 0.02, "n=" + n_vessel, 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', fontsize - 4);
    xticks(ax, [-10 0 10])
    set(ax,'Box','off', ...
                'TickDir','out', ...
                'LineWidth',1, ...
                'FontName','Arial', ...
                'FontSize',fontsize, ...
                'Layer','top')
end

function match_ylim(axes_array)
    ylims = [];
    for i = 1:numel(axes_array)
        ylims = [ylims; ylim(axes_array(i))]; %#ok<AGROW>
    end
    new_yl = [min(ylims(:,1)), max(ylims(:,2))];
    for i = 1:numel(axes_array)
        ylim(axes_array(i), new_yl);
    end
end
