clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
clee = color_lee();
% Which dataset to draw. tablegeneration_main writes transition.mat beside the
% dirtable it was pointed at, so this must name the same folder.
% integrate_analysisresult sets these in the environment, which survives the clear
% above; run this file on its own and the defaults below apply.
param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';      % or '00_igkl' / '01_igkltdt'
end
dirs.secondary_root = getenv('ECSPRESS_ROOT');
if isempty(dirs.secondary_root)
    dirs.secondary_root = ['E:\OneDrive - The Pennsylvania State University\' ...
        '2023ecspress\02_secondary_analysis'];
end
exp_path = fullfile(dirs.secondary_root, param.dataset);
table_path = fullfile(exp_path, "transition.mat");
load_struct = load(table_path);
transition_table = load_struct.save_content;

% taxis calculation
% err  bout_duration is the BOUT's length, (end-start)/fs -- NOT the window. It
%      reads 50 s because get_transition takes +/-transition_window around the
%      transition, and 41 or 46 s where the recording ran out. Building the axis
%      from bout_duration(1) therefore drew whatever row 1 happened to be
% why  the cut is state_linefwhm line 203, center_i +/- onset with
%      onset = transition_window*fs/2 -- HALF the window again. So the trace spans
%      transition_window seconds total, +/- half of it
% caution  25 is state_integration.m line 203. If that changes, this must too:
%          the window is not stored anywhere in the transition table
param.transition_window = 25;      % s, state_integration.param.transition_window
t_axis.duration = param.transition_window;
t_axis.dpoints = size(transition_table.abs_numeric.filtered_table.data{1});
t_axis.tpoints = linspace(-t_axis.duration/2,t_axis.duration/2, t_axis.dpoints(2));

% Color config
figconfig = struct();

figconfig.bv.faint  = clee.lch(80, 130,  10);
figconfig.bv.bold   = clee.lch(45, 150,  15);

figconfig.pvs.faint = clee.lch(80, 130, 110);
figconfig.pvs.bold  = clee.lch(45, 150, 115);

figconfig.eps.faint = clee.lch(80, 130,  70);
figconfig.eps.bold  = clee.lch(45, 120,  75);

% Figure sizing. These land in a slide at roughly half their drawn width, so a
% point size that reads on screen is half that once pasted.
%   see FINDINGS.md
%   why       the tile grew with the font : 16 pt tick labels collide inside a
%             1.5 in tile, and thinning the ticks alone was not enough
%   see FINDINGS.md
param.fontsize       = 14;
param.tile_widthinch = 2.2;
param.tile_heightinch = 3.4;

% Shared definitions
transitions  = ["an_trans", "na_trans", "nr_trans", "ra_trans"];
data_types = ["thickness_bv", "thickness_totalpvs", "thickness_eps"];
vessel_names = ["BV", "PVS", "EPS"];
vessel_colors = {figconfig.bv, figconfig.pvs, figconfig.eps};

label_map = struct( ...
    'an_trans', struct('pre',"Awake",'post',"NREM"), ...
    'na_trans', struct('pre',"NREM", 'post',"Awake"), ...
    'nr_trans', struct('pre',"NREM", 'post',"REM"), ...
    'ra_trans', struct('pre',"REM",  'post',"Awake") ...
);

% ── Generate figures for each normalization type ───────────────────────

norm_fields = ["abs_data", "awakenorm_data", "awakesubtract_data"];
norm_labels = ["absolute", "awake-normalized", "awake-subtracted"];
% Every normalisation, not one. Absolute and relative have both been wanted from
% the start, and which one reads better depends on the panel.
param.norms = [1 2 3];

%% ── Per-vessel figures (4 transition tiles each) ───────────────────
for ni = param.norms
    sess_ave  = transition_table.(norm_fields(ni)).Date_ave;
    mouse_ave = transition_table.(norm_fields(ni)).MouseID_ave;
    plotdata  = build_plotdata(sess_ave, mouse_ave, transitions, data_types);

    spec = struct();
    spec.tile_layout = [1 4];
    spec.fontsize    = param.fontsize;
    spec.tile_inch   = [param.tile_widthinch param.tile_heightinch];
    for vi = 1:numel(data_types)
        spec.title = vessel_names(vi) + " – " + norm_labels(ni);
        spec.tiles = struct('transition', {}, 'datatype', {}, 'color', {});
        for ti = 1:numel(transitions)
            spec.tiles(ti).transition = transitions(ti);
            spec.tiles(ti).datatype   = data_types(vi);
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
    xy_sizeinch = [spec.tile_inch(1)*spec.tile_layout(2) spec.tile_inch(2)];
    set(fig, 'Units','inches', ...
        "Position",[monitor_xyinch(1) monitor_xyinch(2) xy_sizeinch(1) xy_sizeinch(2)])

    tl = tiledlayout(fig, 1, 1, "TileSpacing", "compact");
    tl.GridSize = spec.tile_layout;

    an_axes = [];
    na_axes = [];
    nr_axes = [];
    ra_axes = [];

    for i = 1:numel(spec.tiles)
        ax = axes(Parent=tl);
        ax.Layout.Tile = i;
        tile = spec.tiles(i);

        % Resolve datatype field name (strip "thickness_" prefix)
        dtfield = strrep(tile.datatype, "thickness_", "");
        d = plotdata.(tile.transition).(dtfield);

        labels = label_map.(tile.transition);
        plot_tile(ax, taxis, d.sess, d.mouse, tile.color, labels, spec.fontsize);

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


function plot_tile(ax, taxis, sess_data, mouse_data, cfg, labels, fontsize)
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
    hold(ax, "on")

    % err  the CI patch used to carry FaceAlpha 0.5 and sit BETWEEN the traces and
    %      the mean. print -dsvg -vector cannot express per-object alpha, so it
    %      wraps the axes content in a group opacity and the mean line came out
    %      washed in the deck even though the PNG looked right
    % why  nothing in this axes is transparent now. The patch is an OPAQUE light
    %      grey drawn FIRST, so the session traces and the mean both sit on top of
    %      it and stay at full strength
    % note 0.9 is grey 0.8 at 50% over white -- the band reads the same as before
    ci_upper = mpervesselmat + mci95mat;
    ci_lower = mpervesselmat - mci95mat;
    patch(ax, [taxis, fliplr(taxis)], [ci_upper', fliplr(ci_lower')], ...
        [0.9 0.9 0.9], 'EdgeColor', 'none');

    plot(ax, taxis, pervesselmat, "Color", cfg.faint);
    plot(ax, taxis, mpervesselmat, "Color", cfg.bold, "LineWidth", 2);
    % err  this was hardcoded [-25 25] while the trace spans +/-transition_window/2
    %      = +/-12.5 s, so the data filled only the middle half of every tile
    xlim(ax, [taxis(1) taxis(end)])
    ylim(ax, yl)
    xline(ax, 0, "--")
    % why  normalized units : yl(2)*0.99 put the text at 0.99 of the DATA value, so
    %      on an axis running 6.5-20.5 it landed a quarter of the way down, and on
    %      a subtracted axis crossing zero it landed off-tile
    text(ax, 0.48, 0.98, labels.pre,  'Units','normalized', ...
        "HorizontalAlignment", "right", "VerticalAlignment", "top", "FontSize", fontsize-2)
    text(ax, 0.52, 0.98, labels.post, 'Units','normalized', ...
        "HorizontalAlignment", "left",  "VerticalAlignment", "top", "FontSize", fontsize-2)
    ylabel(ax, "")
    xlabel(ax, "")

    % caution  at 14 pt the default 5-6 x ticks collide and MATLAB silently rotates
    %          them. Three is what fits. y ticks are left to MATLAB : forcing
    %          linspace(yl) gives ragged decimals instead of round numbers
    xticks(ax, [-10 0 10])

    set(ax,'Box','off', ...
                'TickDir','out', ...
                'LineWidth',1, ...
                'FontName','Arial', ...
                'FontSize',fontsize, ...
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
