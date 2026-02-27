clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
clee = color_lee();
exp_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl';
table_path = fullfile(exp_path, "transition.mat");
load_struct = load(table_path);
transition_table = load_struct.save_content;

%% Shared definitions
transitions  = ["an_trans", "na_trans", "nr_trans", "ra_trans"];
transition_labels = ["Awake-NREM", "NREM-Awake", "NREM-REM", "REM-Awake"];
data_types = ["thickness_bv", "thickness_totalpvs", "thickness_eps"];
vessel_names = ["BV", "PVS", "EPS"];

figconfig = struct();
figconfig.bv.faint  = clee.lch(80, 130,  10);
figconfig.bv.bold   = clee.lch(45, 150,  15);
figconfig.pvs.faint = clee.lch(80, 130, 110);
figconfig.pvs.bold  = clee.lch(45, 150, 115);
figconfig.eps.faint = clee.lch(80, 130,  70);
figconfig.eps.bold  = clee.lch(45, 120,  75);
vessel_colors = {figconfig.bv, figconfig.pvs, figconfig.eps};
norm_fields = ["abs_numeric", "awakenorm_numeric", "awakesubtract_numeric"];

%% Extract specific configuration data
% Using awakenorm as per original script
sess_ave  = transition_table.(norm_fields(3)).Date_ave;
mouse_ave = transition_table.(norm_fields(3)).MouseID_ave;
plotdata  = build_plotdata(sess_ave, mouse_ave, transitions, data_types);

%% ── Tiled Figure 1: Grouped by Transition (1 figure, 4 tiles) ─────────────
spec_trans = struct();
spec_trans.title = "Pre-Post Changes grouped by Transition";
spec_trans.tile_layout = [1 4];
for ti = 1:numel(transitions)
    spec_trans.tiles(ti).title = transition_labels(ti);
    for vi = 1:numel(data_types)
        spec_trans.tiles(ti).items(vi).transition = transitions(ti);
        spec_trans.tiles(ti).items(vi).datatype   = data_types(vi);
        spec_trans.tiles(ti).items(vi).color      = vessel_colors{vi};
        spec_trans.tiles(ti).items(vi).label      = vessel_names(vi);
    end
end
render_scatter_figure(spec_trans, plotdata);

%% ── Tiled Figure 2: Grouped by Vessel Type (1 figure, 3 tiles) ────────────
spec_vess = struct();
spec_vess.title = "Pre-Post Changes grouped by Vessel Type";
spec_vess.tile_layout = [1 3];
for vi = 1:numel(data_types)
    spec_vess.tiles(vi).title = vessel_names(vi);
    for ti = 1:numel(transitions)
        spec_vess.tiles(vi).items(ti).transition = transitions(ti);
        spec_vess.tiles(vi).items(ti).datatype   = data_types(vi);
        % Consistent vessel identity colors 
        spec_vess.tiles(vi).items(ti).color      = vessel_colors{vi};
        spec_vess.tiles(vi).items(ti).label      = transition_labels(ti);
    end
end
render_scatter_figure(spec_vess, plotdata);


%% ── Local functions ────────────────────────────────────────────────────

function plotdata = build_plotdata(sess_tbl, mouse_tbl, transitions, datatypes)
    for ti = 1:numel(transitions)
        tr = transitions(ti);
        for di = 1:numel(datatypes)
            dt = datatypes(di);
            % field name: replace non-identifier chars
            dtfield = strrep(dt, "thickness_", "");
            plotdata.(tr).(dtfield).sess  = ...
                sess_tbl(sess_tbl.state_name == tr & sess_tbl.DataType == dt, :);
            plotdata.(tr).(dtfield).mouse = ...
                mouse_tbl(mouse_tbl.state_name == tr & mouse_tbl.DataType == dt, :);
        end
    end
end

function render_scatter_figure(spec, plotdata)
    fig = figure("Name", spec.title);
    monitor_xyinch = [10 2];
    xy_sizeinch = [2.5*spec.tile_layout(2) 3];
    set(fig, 'Units','inches', ...
        "Position",[monitor_xyinch(1) monitor_xyinch(2) xy_sizeinch(1) xy_sizeinch(2)])

    tl = tiledlayout(fig, 1, 1, "TileSpacing", "compact");
    tl.GridSize = spec.tile_layout;

    an_axes = [];
    na_axes = [];
    nr_axes = [];
    ra_axes = [];
    all_axes = [];

    for i = 1:numel(spec.tiles)
        ax = axes(Parent=tl);
        ax.Layout.Tile = i;
        tile = spec.tiles(i);
        plot_scatter_tile(ax, tile, plotdata);
        all_axes = [all_axes, ax];

        % Categorize axes by transition if the tile represents a single transition
        trans_in_tile = unique([tile.items.transition]);
        if numel(trans_in_tile) == 1
            switch trans_in_tile
                case "an_trans", an_axes = [an_axes, ax];
                case "na_trans", na_axes = [na_axes, ax];
                case "nr_trans", nr_axes = [nr_axes, ax];
                case "ra_trans", ra_axes = [ra_axes, ax];
            end
        end
    end

    % Match ylims based on figure configuration
    if ~isempty(an_axes) || ~isempty(na_axes) || ~isempty(nr_axes) || ~isempty(ra_axes)
        % For figures grouped by transition, match the pairs
        match_ylim([an_axes, na_axes]);
        match_ylim([nr_axes, ra_axes]);
    else
        % For figures grouped by vessel (where all transitions share the same tile),
        % match all vessel tiles so magnitudes can be compared directly
        match_ylim(all_axes);
    end
end

function plot_scatter_tile(ax, tile, plotdata)
    hold(ax, 'on');
    
    n_items = numel(tile.items);
    x_positions = 1:n_items;
    x_labels = strings(1, n_items);
    
    for i = 1:n_items
        item = tile.items(i);
        dtfield = strrep(item.datatype, "thickness_", "");
        d = plotdata.(item.transition).(dtfield);
        
        % Calculate Post-Pre change from session level data
        % (adjust column names based on your precise table structure)
        changes = d.sess.post_median - d.sess.pre_median;
        
        % Drop NaNs just in case
        changes(isnan(changes)) = [];
        
        % Add horizontal jitter to scatter points to prevent overlap
        jitter_amount = 0.2;
        x_jitter = x_positions(i) + (rand(size(changes)) - 0.5) * jitter_amount;
        
        % Plot individual data points (faint color for individual points)
        scatter(ax, x_jitter, changes, 25, item.color.faint, 'filled', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
        
        % Calculate statistics
        c_mean = mean(changes);
        
        % 95% Confidence interval over Session changes using t-distribution
        cn = numel(changes);
        t_val = tinv(0.975, cn - 1);
        % Or use standard deviation: c_err = std(changes); 
        c_err = t_val * std(changes) / sqrt(cn); 
        
        % Plot Error bar & Mean dot
        errorbar(ax, x_positions(i), c_mean, c_err, ...
                 'Color', item.color.bold, 'LineWidth', 2, 'CapSize', 6);
        scatter(ax, x_positions(i), c_mean, 60, item.color.bold, 'filled', ...
                 'MarkerEdgeColor', 'w', 'LineWidth', 1);

        x_labels(i) = item.label;
    end
    
    % Draw 0-line
    yline(ax, 0, '--k', 'Alpha', 0.3, 'LineWidth', 1);

    % Configure Limits & Appearance
    title(ax, tile.title, 'FontWeight', 'normal');
    xticks(ax, x_positions);
    xticklabels(ax, x_labels);
    xtickangle(ax, 25);
    xlim(ax, [0.5, n_items + 0.5]);

    set(ax,'Box','off', ...
                'TickDir','out', ...
                'LineWidth',1, ...
                'FontName','Arial', ...
                'FontSize',11, ...
                'Layer','top')
end

function match_ylim(axes_array)
    if ~isempty(axes_array)
        ylims = [];
        for i = 1:numel(axes_array)
            ylims = [ylims; ylim(axes_array(i))];
        end
        new_yl = [min(ylims(:,1)), max(ylims(:,2))];
        for i = 1:numel(axes_array)
            ylim(axes_array(i), new_yl);
        end
    end
end
