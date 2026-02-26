clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
clee = color_lee();
exp_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl';
table_path = fullfile(exp_path, "transition.mat");
load_struct = load(table_path);
transition_table = load_struct.save_content;
%% taxis calculation
t_axis.duration = transition_table.abs_numeric.filtered_table.bout_duration(1);
t_axis.dpoints = size(transition_table.abs_numeric.filtered_table.data{1});
t_axis.tpoints = linspace(-t_axis.duration/2,t_axis.duration/2, t_axis.dpoints(2));
%% Color config
figconfig = struct();

figconfig.bv.faint  = clee.lch(80, 130,  10);
figconfig.bv.bold   = clee.lch(45, 150,  15);

figconfig.pvs.faint = clee.lch(80, 130, 110);
figconfig.pvs.bold  = clee.lch(45, 150, 115);

figconfig.eps.faint = clee.lch(80, 130,  70);
figconfig.eps.bold  = clee.lch(45, 120,  75);

%% Shared definitions
transitions  = ["an_trans", "na_trans", "nr_trans", "ra_trans"];
vessel_types = ["thickness_bv", "thickness_totalpvs", "thickness_eps"];
vessel_names = ["BV", "PVS", "EPS"];
vessel_colors = {figconfig.bv, figconfig.pvs, figconfig.eps};

label_map = struct( ...
    'an_trans', struct('pre',"Awake",'post',"NREM"), ...
    'na_trans', struct('pre',"NREM", 'post',"Awake"), ...
    'nr_trans', struct('pre',"NREM", 'post',"REM"), ...
    'ra_trans', struct('pre',"REM",  'post',"Awake") ...
);

%% ── Generate figures for each normalization type ───────────────────────

norm_fields = ["abs_data", "awakenorm_data", "awakesubtract_data"];
norm_labels = ["absolute", "awake-normalized", "awake-subtracted"];

for ni = 1:numel(norm_fields)
    sess_ave  = transition_table.(norm_fields(ni)).Date_ave;
    mouse_ave = transition_table.(norm_fields(ni)).MouseID_ave;
    plotdata  = build_plotdata(sess_ave, mouse_ave, transitions, vessel_types);

    % ── Per-vessel figures (4 transition tiles each) ───────────────────
    for vi = 1:numel(vessel_types)
        spec = struct();
        spec.title = vessel_names(vi) + " – " + norm_labels(ni);
        spec.tile_layout = [1 4];
        for ti = 1:numel(transitions)
            spec.tiles(ti).transition = transitions(ti);
            spec.tiles(ti).datatype   = vessel_types(vi);
            spec.tiles(ti).color      = vessel_colors{vi};
        end
        render_figure(spec, plotdata, t_axis.tpoints, label_map);
    end
end


    % % ── Per-transition figures (3 vessel-type tiles each) ──────────────
    % for ti = 1:numel(transitions)
    %     spec = struct();
    %     spec.title = label_map.(transitions(ti)).pre + "-" ...
    %                + label_map.(transitions(ti)).post + " – " + norm_labels(ni);
    %     spec.tile_layout = [1 numel(vessel_types)];
    %     for vi = 1:numel(vessel_types)
    %         spec.tiles(vi).transition = transitions(ti);
    %         spec.tiles(vi).datatype   = vessel_types(vi);
    %         spec.tiles(vi).color      = vessel_colors{vi};
    %     end
    %     render_figure(spec, plotdata, t_axis.tpoints, label_map);
    % end



%% ── Local functions ────────────────────────────────────────────────────

function plotdata = build_plotdata(sess_tbl, mouse_tbl, transitions, datatypes)
% BUILD_PLOTDATA  Index all data by (transition, datatype).
%   plotdata.(trans).(dtype).sess  = filtered session-level table rows
%   plotdata.(trans).(dtype).mouse = filtered mouse-level table rows
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


function render_figure(spec, plotdata, taxis, label_map)
% RENDER_FIGURE  Create a figure from a spec struct.
%   spec.title        – figure window title
%   spec.tile_layout  – [rows cols]
%   spec.tiles(i)     – struct with .transition, .datatype, .color

    fig = figure("Name", spec.title);
    monitor_xyinch = [10 2];
    xy_sizeinch = [1.5*spec.tile_layout(2) 3];
    set(fig, 'Units','inches', ...
        "Position",[monitor_xyinch(1) monitor_xyinch(2) xy_sizeinch(1) xy_sizeinch(2)])

    tl = tiledlayout(fig, spec.tile_layout(1), spec.tile_layout(2), ...
        "TileSpacing", "compact");

    an_axes = [];
    na_axes = [];
    nr_axes = [];
    ra_axes = [];

    for i = 1:numel(spec.tiles)
        ax = nexttile(tl);
        tile = spec.tiles(i);

        % Resolve datatype field name (strip "thickness_" prefix)
        dtfield = strrep(tile.datatype, "thickness_", "");
        d = plotdata.(tile.transition).(dtfield);

        labels = label_map.(tile.transition);
        plot_tile(ax, taxis, d.sess, d.mouse, tile.color, labels);

        switch tile.transition
            case "an_trans"
                an_axes = [an_axes, ax];
            case "na_trans"
                na_axes = [na_axes, ax];
            case "nr_trans"
                nr_axes = [nr_axes, ax];
            case "ra_trans"
                ra_axes = [ra_axes, ax];
        end
    end

    match_ylim(an_axes, na_axes);
    match_ylim(nr_axes, ra_axes);
end


function plot_tile(ax, taxis, sess_data, mouse_data, cfg, labels)
% PLOT_TILE  Draw session traces + mouse-average on a single axes.
    pervesselmat  = cell2mat(sess_data.data)';
    mpervesselmat = cell2mat(mouse_data.data)';
    
    % Recalculate 95% CI from sess_data using t-distribution
    n_samples = size(pervesselmat, 2);
    t_val = tinv(0.975, n_samples - 1);
    mci95mat = t_val * std(pervesselmat, 0, 2, 'omitnan') / sqrt(n_samples);

    % Auto-scale ylim to data range with 5% padding, unless cfg.ylim is set
    if isfield(cfg, 'ylim') && ~isempty(cfg.ylim)
        yl = cfg.ylim;
    else
        data_min = min(pervesselmat(:), [], 'omitnan');
        data_max = max(pervesselmat(:), [], 'omitnan');
        pad = (data_max - data_min) * 0.05;
        yl = [data_min - pad, data_max + pad];
    end

    cla(ax)
    plot(ax, taxis, pervesselmat,  "Color", cfg.faint);
    hold(ax, "on")

    % Plot 95% confidence interval patch
    ci_upper = mpervesselmat + mci95mat;
    ci_lower = mpervesselmat - mci95mat;
    patch(ax, [taxis, fliplr(taxis)], [ci_upper', fliplr(ci_lower')], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    plot(ax, taxis, mpervesselmat, "Color", cfg.bold, "LineWidth", 2);
    xlim(ax, [-25 25])
    ylim(ax, yl)
    xline(ax, 0, "--")
    text(ax, -1, yl(2)*0.99, labels.pre,  "HorizontalAlignment", "right", "FontSize", 11)
    text(ax,  1, yl(2)*0.99, labels.post, "HorizontalAlignment", "left",  "FontSize", 11)
    ylabel(ax, "")
    xlabel(ax, "")

    set(ax,'Box','off', ...
                'TickDir','out', ...
                'LineWidth',1, ...
                'FontName','Arial', ...
                'FontSize',11, ...
                'Layer','top')
end

function match_ylim(axes1, axes2)
    if ~isempty(axes1) && ~isempty(axes2)
        ylims = [];
        all_axes = [axes1, axes2];
        for i = 1:numel(all_axes)
            ylims = [ylims; ylim(all_axes(i))];
        end
        new_yl = [min(ylims(:,1)), max(ylims(:,2))];
        for i = 1:numel(all_axes)
            ylim(all_axes(i), new_yl);
        end
    end
end
