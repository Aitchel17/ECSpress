clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
clee = color_lee();
exp_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl';
table_path = fullfile(exp_path, "state_summary.mat");
load_struct = load(table_path);
statsummary_table = load_struct.save_content;


%% Shared definitions
states  = ["awake", "drowsy", "nrem", "rem"];
states_labels = ["Awake", "Drowsy", "NREM", "REM"];
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

%% Define labels explicitly
norm_fields = ["abs_numeric", "awakenorm_numeric", "awakesubtract_numeric"];

%%
sess_ave  = statsummary_table.(norm_fields(3)).Date_ave;
mouse_ave = statsummary_table.(norm_fields(3)).MouseID_ave;
plotdata  = build_plotdata(sess_ave, mouse_ave, states, data_types);
%% ── Tiled Figure 1: Grouped by Transition (1 figure, 4 tiles) ─────────────
spec_states = struct();
spec_states.title = "Q2Q3 mean Grouped by states";
spec_states.tile_layout = [1 4];
%%
for ti = 1:numel(states)
    spec_states.tiles(ti).title = states_labels(ti);
    for vi = 1:numel(data_types)
        spec_states.tiles(ti).items(vi).states = states(ti);
        spec_states.tiles(ti).items(vi).datatype   = data_types(vi);
        spec_states.tiles(ti).items(vi).color      = vessel_colors{vi};
        spec_states.tiles(ti).items(vi).label      = vessel_names(vi);
    end
end
%%
render_scatter_figure(spec_states, plotdata);



%%
function plotdata = build_plotdata(sess_tbl, mouse_tbl, states, datatypes)
    for si = 1:numel(states)
        sr = states(si);
        for di = 1:numel(datatypes)
            dt = datatypes(di);
            % field name: replace non-identifier chars
            dtfield = strrep(dt, "thickness_", "");
            plotdata.(sr).(dtfield).sess  = ...
                sess_tbl(sess_tbl.state_name == sr & sess_tbl.DataType == dt, :);
            plotdata.(sr).(dtfield).mouse = ...
                mouse_tbl(mouse_tbl.state_name == sr & mouse_tbl.DataType == dt, :);
        end
    end
end

%%
function render_scatter_figure(spec, plotdata)
    fig = figure("Name", spec.title);
    monitor_xyinch = [10 2];
    xy_sizeinch = [2.5*spec.tile_layout(2) 3];
    set(fig, 'Units','inches', ...
        "Position",[monitor_xyinch(1) monitor_xyinch(2) xy_sizeinch(1) xy_sizeinch(2)])

    tl = tiledlayout(fig, 1, 1, "TileSpacing", "compact");
    tl.GridSize = spec.tile_layout;



    for i = 1:numel(spec.tiles)
        ax = axes(Parent=tl);
        ax.Layout.Tile = i;
        tile = spec.tiles(i);
        plot_scatter_tile(ax, tile, plotdata);
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
        d = plotdata.(item.states).(dtfield);
        
        % Calculate Post-Pre change from session level data
        % (adjust column names based on your precise table structure)
        q2q3mean = d.sess.Q2Q3_mean';
        
        % Drop NaNs just in case
        q2q3mean(isnan(q2q3mean)) = [];
        
        % Add horizontal jitter to scatter points to prevent overlap
        jitter_amount = 0.2;
        x_jitter = x_positions(i) + (rand(size(q2q3mean)) - 0.5) * jitter_amount;
        
        % Plot individual data points (faint color for individual points)
        scatter(ax, x_jitter, q2q3mean, 25, item.color.faint, 'filled', ...
            'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
        
        % Calculate statistics
        c_mean = mean(q2q3mean);
        
        % 95% Confidence interval over Session changes using t-distribution
        cn = numel(q2q3mean);
        t_val = tinv(0.975, cn - 1);
        % Or use standard deviation: c_err = std(changes); 
        c_err = t_val * std(q2q3mean) / sqrt(cn); 
        
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