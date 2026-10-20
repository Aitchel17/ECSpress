clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
addpath('g:/03_program/01_ecspress/09_dirstruct');    % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
clee = color_lee();
% Which dataset to draw. tablegeneration_paxfwhm writes transition.mat beside the
% dirtable it was pointed at, so this must name the same folder.
param.dataset = 'merged';   % char  the dataset whose tables this reads
dirs = dirs_central();
exp_path = fullfile(dirs.secondary_root, param.dataset);
table_path = fullfile(exp_path, "transition.mat");
load_struct = load(table_path);
transition_table = load_struct.save_content;

%% Filters : which rows the figure is about
filt.transition_drawn = ["an_trans", "na_trans", "nr_trans", "ra_trans"];   % (str)
filt.series_drawn = ["thickness_bv", "thickness_totalpvs", "thickness_eps"];

% Labels and colours are keyed by the name they belong to, not by position, so
% dropping one from filt cannot slide the rest along by one.
param.transition_label = struct('an_trans', "Awake-NREM", 'na_trans', "NREM-Awake", ...
    'nr_trans', "NREM-REM", 'ra_trans', "REM-Awake");
param.series_label = struct('thickness_bv', "BV", 'thickness_totalpvs', "PVS", ...
    'thickness_eps', "EPS");

% Figure sizing. Matched to transition_fig so the two read as one panel in a
% slide -- see the note in render_scatter_figure for why the height moves and
% not the width.
%   see FINDINGS.md
param.fontsize       = 14;
param.tile_widthinch = 2.8;

% why  a shared y range is set by whichever vessel moves most, so BV's -5.2..3.9
%      forced EPS onto a -6..4 axis where its own data only spans -3.2..3.1 and
%      the spread became unreadable. Each tile autoscales instead
% caution  magnitudes are then NOT comparable BY EYE across tiles -- read the axis
param.match_ylim = false;

param.series_colour = struct( ...
    'thickness_bv',       struct('faint', clee.lch(80, 130,  10), 'bold', clee.lch(45, 150,  15)), ...
    'thickness_totalpvs', struct('faint', clee.lch(80, 130, 110), 'bold', clee.lch(45, 150, 115)), ...
    'thickness_eps',      struct('faint', clee.lch(80, 130,  70), 'bold', clee.lch(45, 120,  75)));

%% Extract specific configuration data
% Every normalisation, not one. Absolute and relative have both been wanted from
% the start; transition_fig draws the same set so the two stay comparable.
param.norm_drawn = ["abs_numeric", "awakenorm_numeric", "awakesubtract_numeric"];
param.norm_label = struct('abs_numeric', "absolute", ...
    'awakenorm_numeric', "awake-normalized", 'awakesubtract_numeric', "awake-subtracted");

for ni = 1:numel(param.norm_drawn)
norm_field = param.norm_drawn(ni);
norm_name = param.norm_label.(norm_field);
sess_ave  = transition_table.(norm_field).Date_ave;
mouse_ave = transition_table.(norm_field).MouseID_ave;
plotdata  = build_plotdata(sess_ave, mouse_ave, filt.transition_drawn, filt.series_drawn);

%% ── Tiled Figure 1: Grouped by Transition (1 figure, 4 tiles) ─────────────
spec_trans = struct();
spec_trans.title = "Pre-Post by Transition – " + norm_name;
spec_trans.tile_layout = [1 numel(filt.transition_drawn)];
spec_trans.fontsize    = param.fontsize;
spec_trans.tile_width  = param.tile_widthinch;
spec_trans.match_ylim  = param.match_ylim;
spec_trans.tiles = struct('title', {}, 'items', {});
for ti = 1:numel(filt.transition_drawn)
    transition_name = filt.transition_drawn(ti);
    spec_trans.tiles(ti).title = param.transition_label.(transition_name);
    for vi = 1:numel(filt.series_drawn)
        series_name = filt.series_drawn(vi);
        spec_trans.tiles(ti).items(vi).transition = transition_name;
        spec_trans.tiles(ti).items(vi).datatype   = series_name;
        spec_trans.tiles(ti).items(vi).color      = param.series_colour.(series_name);
        spec_trans.tiles(ti).items(vi).label      = param.series_label.(series_name);
    end
end
render_scatter_figure(spec_trans, plotdata);

%% ── Tiled Figure 2: Grouped by Vessel Type (1 figure, 3 tiles) ────────────
spec_vess = struct();
spec_vess.title = "Pre-Post by Vessel Type – " + norm_name;
spec_vess.tile_layout = [1 numel(filt.series_drawn)];
spec_vess.fontsize    = param.fontsize;
spec_vess.tile_width  = param.tile_widthinch;
spec_vess.match_ylim  = param.match_ylim;
spec_vess.tiles = struct('title', {}, 'items', {});
for vi = 1:numel(filt.series_drawn)
    series_name = filt.series_drawn(vi);
    spec_vess.tiles(vi).title = param.series_label.(series_name);
    for ti = 1:numel(filt.transition_drawn)
        transition_name = filt.transition_drawn(ti);
        spec_vess.tiles(vi).items(ti).transition = transition_name;
        spec_vess.tiles(vi).items(ti).datatype   = series_name;
        % Consistent vessel identity colors
        spec_vess.tiles(vi).items(ti).color      = param.series_colour.(series_name);
        spec_vess.tiles(vi).items(ti).label      = param.transition_label.(transition_name);
    end
end
render_scatter_figure(spec_vess, plotdata);
end   % ni, the normalisation loop


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
    % why  transition_fig draws 1.5 x 3 inch tiles (h/w = 2). At 2.5 x 3 the
    %      prepost tiles were flatter than the trace tiles beside them, and the
    %      4-tile figure came out 10 x 3 -- too wide to read the spread
    % note  width per tile stays 2.5 : 3-4 jittered categories with angled labels
    %       do not fit in 1.5, so the height moves instead
    tile_aspect = 2;            % float   height / width, matched to transition_fig
    xy_sizeinch = [spec.tile_width*spec.tile_layout(2) spec.tile_width*tile_aspect];
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
        plot_scatter_tile(ax, tile, plotdata, spec.fontsize);
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
    if ~spec.match_ylim
        return      % every tile keeps its own autoscale
    end
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

function plot_scatter_tile(ax, tile, plotdata, fontsize)
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
        % err  MarkerFaceAlpha 0.6 was the same trap as the CI patch in
        %      transition_fig : print -dsvg cannot express per-object alpha, so
        %      the whole axes goes through a group opacity and the bold mean dot
        %      washes out in the deck. The faint colour is already light enough
        scatter(ax, x_jitter, changes, 25, item.color.faint, 'filled', ...
            'MarkerEdgeColor', 'none');
        
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
    % caution  the full "Awake-NREM" labels overlap once the font passes ~13 pt.
    %          The arrow forms carry the same information in a third of the width,
    %          which is what lets the angle go back to 0. NREM must be replaced
    %          before REM or "NREM" becomes "NR"
    short_labels = replace(x_labels,     "Awake", "A");
    short_labels = replace(short_labels, "NREM",  "N");
    short_labels = replace(short_labels, "REM",   "R");
    short_labels = replace(short_labels, "-",     char(8594));   % right arrow
    xticklabels(ax, short_labels);
    % see FINDINGS.md
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
