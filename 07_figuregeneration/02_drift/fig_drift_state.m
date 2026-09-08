%FIG_DRIFT_STATE  Drift level by sleep state, one dot per vessel, both components.
%   The drift counterpart of tonic_state. Reads drift_summary.mat and draws one
%   figure per normalisation: absolute (um about the session's own trend line) and
%   awake-subtracted. Tiles are the two components, items the states.
%
%   What differs from tonic_state, all settled in DRIFT_FIGURE_PLAN.md: two
%   components instead of three series, no awake-normalised arm, state colours from
%   color_lee.clist, the reference line at 0, one render, seeded jitter, no alpha.
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
table_path = fullfile(exp_path, "drift_summary.mat");
load_struct = load(table_path);
summary_table = load_struct.save_content;

%% Shared definitions
% +row is image-down, which is caudal. The axis reads rostral-positive, so row is
% negated here and nowhere else.
components       = ["drift_row", "drift_col"];
component_labels = ["rostral − caudal", "lateral − medial"];
component_sign   = [-1, 1];

norm_fields = ["abs_numeric", "awakesubtract_numeric"];
norm_labels = ["absolute", "awake-subtracted"];
% awake is a structural zero once subtracted, so it is drawn only on the absolute
norm_states = {["awake", "drowsy", "nrem", "rem"], ["drowsy", "nrem", "rem"]};
state_labels = struct('awake', "Awake", 'drowsy', "Drowsy", 'nrem', "NREM", 'rem', "REM");

param.fontsize        = 14;
param.tile_widthinch  = 2.8;
param.jitter_seed     = 1;      % count  so the jitter redraws the same
param.match_ylim      = true;   % both tiles are um on orthogonal axes of one vector

%% One figure per normalisation
for ni = 1:numel(norm_fields)
    vessel_ave = summary_table.(norm_fields(ni)).Date_ave;
    states = norm_states{ni};
    spec = struct();
    spec.title       = "Drift level by state – " + norm_labels(ni);
    spec.ylabel      = "level (um)";
    spec.tile_layout = [1 numel(components)];
    spec.fontsize    = param.fontsize;
    spec.tile_width  = param.tile_widthinch;
    spec.match_ylim  = param.match_ylim;
    spec.tiles = struct('title', {}, 'items', {});
    for ci = 1:numel(components)
        spec.tiles(ci).title = component_labels(ci);
        for si = 1:numel(states)
            on_item = vessel_ave.state_name == states(si) & vessel_ave.DataType == components(ci);
            spec.tiles(ci).items(si).values = component_sign(ci) * vessel_ave.Q2Q3_mean(on_item)';
            spec.tiles(ci).items(si).color  = state_color(clee, states(si));
            spec.tiles(ci).items(si).label  = state_labels.(states(si));
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
        values = item.values(~isnan(item.values));
        n_vessel = numel(values);

        jitter_amount = 0.2;
        x_jitter = x_positions(i) + (rand(1, n_vessel) - 0.5) * jitter_amount;
        % no MarkerFaceAlpha : print -dsvg cannot express per-object alpha and the
        % bold dot washes out. The faint colour is already light enough
        scatter(ax, x_jitter, values, 25, item.color.faint, 'filled', ...
            'MarkerEdgeColor', 'none');

        % 95 % interval over vessels, t distribution
        c_mean = mean(values);
        t_val = tinv(0.975, n_vessel - 1);
        c_err = t_val * std(values) / sqrt(n_vessel);
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
    xticklabels(ax, x_labels);
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
