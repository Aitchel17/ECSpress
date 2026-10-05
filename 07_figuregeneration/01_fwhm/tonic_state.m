clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
addpath('g:/03_program/01_ecspress/09_dirstruct');    % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
clee = color_lee();
% integrate_analysisresult sets these in the environment, which survives the clear
% above; run this file on its own and the defaults below apply.
%   this used to be pinned to 00_igkl while transition_fig drew merged_igkl_igkltdt,
%   so the two figures were of different datasets. see CLAUDE_LOG.md
param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
dirs = dirs_central();
exp_path = fullfile(dirs.secondary_root, param.dataset);
table_path = fullfile(exp_path, "state_summary.mat");
load_struct = load(table_path);
statsummary_table = load_struct.save_content;


%% Filters : which rows the figure is about
filt.state_drawn  = ["drowsy", "nrem", "rem"];    % one item per state, in this order (str)
filt.series_drawn = ["thickness_bv", "thickness_totalpvs", "thickness_eps"];   % one tile each

%% Picture
% Labels and colours are keyed by the name they belong to, not by position, so
% dropping one from filt cannot slide the rest along by one.
param.state_label = struct('awake', "Awake", 'drowsy', "Drowsy", ...
    'nrem', "NREM", 'rem', "REM");
param.series_label = struct('thickness_bv', "BV", 'thickness_totalpvs', "PVS", ...
    'thickness_eps', "EPS");
param.series_colour = struct( ...
    'thickness_bv',       struct('faint', clee.lch(80, 130,  10), 'bold', clee.lch(45, 150,  15)), ...
    'thickness_totalpvs', struct('faint', clee.lch(80, 130, 110), 'bold', clee.lch(45, 150, 115)), ...
    'thickness_eps',      struct('faint', clee.lch(80, 130,  70), 'bold', clee.lch(45, 120,  75)));

param.norm_drawn = ["abs_numeric", "awakenorm_numeric", "awakesubtract_numeric"];
param.norm_label = struct('abs_numeric', "absolute", ...
    'awakenorm_numeric', "awake-normalized", 'awakesubtract_numeric', "awake-subtracted");
% The reference line is where a state that did nothing would sit, and each arm puts
% that somewhere else: awake divides to exactly 1 and subtracts to exactly 0. The
% absolute arm has no such fixed point -- micrometres of lumen are not read against
% anything -- so it gets no line. NaN means none.
param.norm_refline = struct('abs_numeric', NaN, ...
    'awakenorm_numeric', 1, 'awakesubtract_numeric', 0);

%% One figure per normalisation, one tile per series, one item per state
% A first spec was built here and then overwritten by an identical block, and the
% surviving one was rendered twice, so the file opened two identical windows.
for ni = 1:numel(param.norm_drawn)
    norm_field = param.norm_drawn(ni);
    sess_ave  = statsummary_table.(norm_field).Date_ave;
    mouse_ave = statsummary_table.(norm_field).MouseID_ave;
    plotdata  = build_plotdata(sess_ave, mouse_ave, filt.state_drawn, filt.series_drawn);

    spec_states = struct();
    spec_states.title = "Q2Q3 mean by state – " + param.norm_label.(norm_field);
    spec_states.tile_layout = [1 numel(filt.series_drawn)];
    spec_states.refline = param.norm_refline.(norm_field);
    spec_states.tiles = struct('title', {}, 'items', {});
    for vi = 1:numel(filt.series_drawn)
        series_name = filt.series_drawn(vi);
        spec_states.tiles(vi).title = param.series_label.(series_name);
        for ti = 1:numel(filt.state_drawn)
            state_name = filt.state_drawn(ti);
            spec_states.tiles(vi).items(ti).states   = state_name;
            spec_states.tiles(vi).items(ti).datatype = series_name;
            spec_states.tiles(vi).items(ti).color    = param.series_colour.(series_name);
            spec_states.tiles(vi).items(ti).label    = param.state_label.(state_name);
        end
    end
    render_scatter_figure(spec_states, plotdata);
end



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
    xy_sizeinch = [2*spec.tile_layout(2) 4];
    set(fig, 'Units','inches', ...
        "Position",[monitor_xyinch(1) monitor_xyinch(2) xy_sizeinch(1) xy_sizeinch(2)])

    tl = tiledlayout(fig, 1, 1, "TileSpacing", "compact");
    tl.GridSize = spec.tile_layout;



    for i = 1:numel(spec.tiles)
        ax = axes(Parent=tl);
        ax.Layout.Tile = i;
        tile = spec.tiles(i);
        plot_scatter_tile(ax, tile, plotdata, spec.refline);
    end

end

function plot_scatter_tile(ax, tile, plotdata, refline)
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
                 'Color', item.color.bold, 'LineWidth', 1, 'CapSize', 3);
        scatter(ax, x_positions(i), c_mean, 60, item.color.bold, 'filled', ...
                 'MarkerEdgeColor', 'w', 'LineWidth', 1);

        x_labels(i) = item.label;
    end
    
    % where a state that did nothing would sit, if the arm has such a place
    if isfinite(refline)
        yline(ax, refline, '--k', 'Alpha', 0.3, 'LineWidth', 1);
    end

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
