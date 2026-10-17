%FIG_DRIFT_PREPOST  Drift across each transition, post minus pre, one dot per vessel,
%   both components. The drift counterpart of transition_prepost. Reads
%   drift_transition.mat and draws one figure: tiles are the two components, items
%   the four transitions.
%
%   What differs from transition_prepost, all settled in DRIFT_FIGURE_PLAN.md: two
%   components instead of three series, the absolute arm only, colour by the state
%   being entered, y matched across the two tiles, seeded jitter. Its second
%   grouping (by series) has no counterpart: with two series it is the same ten
%   marks rearranged.
clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
addpath('g:/03_program/01_ecspress/09_dirstruct');    % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
clee = color_lee();
param.dataset = 'merged';   % char  the dataset whose tables this reads
dirs = dirs_central();
exp_path = fullfile(dirs.secondary_root, param.dataset);
table_path = fullfile(exp_path, "drift_transition.mat");
load_struct = load(table_path);
transition_table = load_struct.save_content;

%% Filters : which rows the figure is about
filt.component_drawn = ["drift_row", "drift_col"];   % one tile each (str)
filt.transition_drawn = ["an_trans", "na_trans", "nr_trans", "ra_trans"];   % one item each

% Keyed by the name, not by position, so dropping one from filt cannot slide the
% rest along by one. +row is image-down, which is caudal, and the axis reads
% rostral-positive, so the sign flips here and nowhere else.
param.component_label = struct('drift_row', "rostral motion", ...
    'drift_col', "lateral motion");
param.component_sign = struct('drift_row', -1, 'drift_col', 1);
param.transition_label = struct('an_trans', "Awake-NREM", 'na_trans', "NREM-Awake", ...
    'nr_trans', "NREM-REM", 'ra_trans', "REM-Awake");
param.transition_enters = struct('an_trans', "nrem", 'na_trans', "awake", ...
    'nr_trans', "rem", 'ra_trans', "awake");

% post minus pre cancels any per-session constant, so the awake-subtracted arm
% is the absolute one again and is not drawn
param.norm_drawn = "abs_numeric";
param.norm_label = struct('abs_numeric', "absolute");

param.fontsize        = 14;
param.tile_widthinch  = 2.8;
param.jitter_seed     = 1;      % count  so the jitter redraws the same
param.match_ylim      = true;   % both tiles are um on orthogonal axes of one vector

%% One figure per normalisation, tiles = components, items = transitions
for ni = 1:numel(param.norm_drawn)
    norm_field = param.norm_drawn(ni);
    vessel_ave = transition_table.(norm_field).Date_ave;
    spec = struct();
    spec.title       = "Drift post − pre by transition – " + param.norm_label.(norm_field);
    spec.ylabel      = "post − pre (um)";
    spec.tile_layout = [1 numel(filt.component_drawn)];
    spec.fontsize    = param.fontsize;
    spec.tile_width  = param.tile_widthinch;
    spec.match_ylim  = param.match_ylim;
    spec.tiles = struct('title', {}, 'items', {});
    for ci = 1:numel(filt.component_drawn)
        component_name = filt.component_drawn(ci);
        spec.tiles(ci).title = param.component_label.(component_name);
        for ti = 1:numel(filt.transition_drawn)
            transition_name = filt.transition_drawn(ti);
            on_item = vessel_ave.state_name == transition_name & vessel_ave.DataType == component_name;
            rows = vessel_ave(on_item, :);
            spec.tiles(ci).items(ti).values = param.component_sign.(component_name) ...
                * (rows.post_median - rows.pre_median)';
            spec.tiles(ci).items(ti).color  = state_color(clee, param.transition_enters.(transition_name));
            spec.tiles(ci).items(ti).label  = param.transition_label.(transition_name);
        end
    end
    rng(param.jitter_seed)
    render_scatter_figure(spec);
end

%% ── Local functions ────────────────────────────────────────────────────

function cfg = state_color(clee, state_name)
    cfg.bold  = clee.clist.(state_name);
    cfg.faint = 1 - 0.55 * (1 - cfg.bold);
end

function render_scatter_figure(spec)
    fig = figure("Name", spec.title);
    monitor_xyinch = [10 2];
    tile_aspect = 2;            % float   height / width, matched to transition_prepost
    xy_sizeinch = [spec.tile_width*spec.tile_layout(2) spec.tile_width*tile_aspect];
    set(fig, 'Units','inches', ...
        "Position",[monitor_xyinch(1) monitor_xyinch(2) xy_sizeinch(1) xy_sizeinch(2)])

    tl = tiledlayout(fig, 1, 1, "TileSpacing", "compact");
    tl.GridSize = spec.tile_layout;

    all_axes = [];
    for i = 1:numel(spec.tiles)
        ax = axes(Parent=tl);
        ax.Layout.Tile = i;
        plot_scatter_tile(ax, spec.tiles(i), spec.fontsize);
        all_axes = [all_axes, ax]; %#ok<AGROW>
    end
    ylabel(all_axes(1), spec.ylabel)
    if spec.match_ylim
        match_ylim(all_axes);
    end
end

function plot_scatter_tile(ax, tile, fontsize)
    hold(ax, 'on');
    n_items = numel(tile.items);
    x_positions = 1:n_items;
    x_labels = strings(1, n_items);
    for i = 1:n_items
        item = tile.items(i);
        changes = item.values(~isnan(item.values));
        n_vessel = numel(changes);

        jitter_amount = 0.2;
        x_jitter = x_positions(i) + (rand(1, n_vessel) - 0.5) * jitter_amount;
        % no MarkerFaceAlpha : print -dsvg cannot express per-object alpha and the
        % bold dot washes out. The faint colour is already light enough
        scatter(ax, x_jitter, changes, 25, item.color.faint, 'filled', ...
            'MarkerEdgeColor', 'none');

        % 95 % interval over vessels, t distribution
        c_mean = mean(changes);
        t_val = tinv(0.975, n_vessel - 1);
        c_err = t_val * std(changes) / sqrt(n_vessel);
        errorbar(ax, x_positions(i), c_mean, c_err, ...
                 'Color', item.color.bold, 'LineWidth', 2, 'CapSize', 6);
        scatter(ax, x_positions(i), c_mean, 60, item.color.bold, 'filled', ...
                 'MarkerEdgeColor', 'w', 'LineWidth', 1);
        x_norm = (x_positions(i) - 0.5) / n_items;     % xlim is [0.5, n_items + 0.5]
        text(ax, x_norm, 0.02, "n=" + n_vessel, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontsize - 4);

        x_labels(i) = item.label;
    end
    yline(ax, 0, '--k', 'Alpha', 0.3, 'LineWidth', 1);

    title(ax, tile.title, 'FontWeight', 'normal');
    xticks(ax, x_positions);
    % the arrow forms carry the label in a third of the width. NREM before REM, or
    % "NREM" becomes "NR"
    short_labels = replace(x_labels,     "Awake", "A");
    short_labels = replace(short_labels, "NREM",  "N");
    short_labels = replace(short_labels, "REM",   "R");
    short_labels = replace(short_labels, "-",     char(8594));
    xticklabels(ax, short_labels);
    if n_items >= 4
        xtickangle(ax, 30);
    else
        xtickangle(ax, 0);
    end
    xlim(ax, [0.5, n_items + 0.5]);
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
