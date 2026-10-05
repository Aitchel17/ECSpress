%FIG_FWHM_STATEDIAMETER  Lumen diameter by arousal state, one vessel per observation.
%   Reads state_summary.mat and overlays the states twice: a kernel density and the
%   cumulative fraction. Nothing is measured here: Q2Q3_mean is the middle half of
%   each bout, already folded to one value per vessel by tablegeneration_main.
%
%   Neither panel bins. About thirty vessels carry a state, so a histogram of them
%   swings between one and ten per bin on Poisson noise alone and the eye reads the
%   swings as structure; the per-bin counts go to the console instead.
%
%   The x axis is absolute micrometres, not the deviation the state heatmap uses --
%   fig_fwhm_statepool centres every session on its own median, so calibre cannot be
%   read off that figure at all.
clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
addpath('g:/03_program/01_ecspress/09_dirstruct');    % where dirs_ecspath lives
dirs_ecspath;                                         % three roots, minus zz_notinuse
clee = color_lee();
param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
dirs = dirs_central();
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);

%% Filters : which rows the figure is about
filt.data_type = "thickness_bv";                          % str  the lumen series
filt.states = ["awake", "drowsy", "nrem", "rem"];         % str  drawn, in this order

%% Picture
% keyed by the state's own name, so dropping one from filt.states cannot slide the
% legend along by one
param.state_label = struct('awake', "Awake", 'drowsy', "Drowsy", ...
    'nrem', "NREM", 'rem', "REM");
param.bin_width = 2;              % um   the console counts only, nothing is drawn on bins
param.fontsize = 12;
param.linewidth = 1.5;

%% What is drawn
table_path = fullfile(dirs.save_dir, "state_summary.mat");
vessel_ave = load(table_path).save_content.abs_numeric.Date_ave;
on_series = vessel_ave.DataType == filt.data_type;
in_series = vessel_ave(on_series, :);
fprintf('%s\n   %s : %d vessel rows\n', table_path, filt.data_type, height(in_series));

n_state = numel(filt.states);
diameter = cell(n_state, 1);
for k = 1:n_state
    on_state = in_series.state_name == filt.states(k);
    values = in_series.Q2Q3_mean(on_state);
    diameter{k} = values(isfinite(values));
end

% one set of edges for every state, or the bars are not comparable
all_diameter = vertcat(diameter{:});
bin_edges = floor(min(all_diameter)) : param.bin_width : ceil(max(all_diameter)) + param.bin_width;

%% The board
fig_diameter = figure('Color', 'w', 'Position', [60 60 1080 440], ...
    'Name', 'fig_fwhm_statediameter');
layout_diameter = tiledlayout(fig_diameter, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
axis_hist = nexttile(layout_diameter);
axis_ecdf = nexttile(layout_diameter);
hold(axis_hist, 'on')
hold(axis_ecdf, 'on')

legend_entry = strings(1, n_state);
for k = 1:n_state
    state_colour = clee.clist.(filt.states(k));
    % stairs, not filled bars : print -dsvg cannot express per-object alpha, so four
    % overlapping filled histograms come out of the deck as one solid block
    histogram(axis_hist, diameter{k}, bin_edges, ...
        'DisplayStyle', 'stairs', 'EdgeColor', state_colour, ...
        'LineWidth', param.linewidth);
    % the same vessels with no bin to choose : every vessel is one step up
    sorted_diameter = sort(diameter{k});
    cumulative = (1:numel(sorted_diameter)) / numel(sorted_diameter);
    stairs(axis_ecdf, sorted_diameter, cumulative, 'Color', state_colour, ...
        'LineWidth', param.linewidth);
    legend_entry(k) = sprintf('%s  n=%d, median %.1f um', ...
        param.state_label.(filt.states(k)), numel(diameter{k}), median(diameter{k}));
    fprintf('   %-7s %2d vessels | %.1f to %.1f um, median %.1f\n', filt.states(k), ...
        numel(diameter{k}), min(diameter{k}), max(diameter{k}), median(diameter{k}));
end

legend(axis_hist, legend_entry, 'Location', 'northeast', 'Box', 'off', ...
    'FontSize', param.fontsize - 1);
ylabel(axis_hist, 'vessels');
ylabel(axis_ecdf, 'fraction of vessels at or below');
ylim(axis_ecdf, [0 1]);
yline(axis_ecdf, 0.5, '--k', 'Alpha', 0.3, 'LineWidth', 1);
xlabel(layout_diameter, 'lumen diameter (um)');
for ax = [axis_hist, axis_ecdf]
    xlim(ax, [bin_edges(1), bin_edges(end)]);
    set(ax, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, ...
        'FontName', 'Arial', 'FontSize', param.fontsize, 'Layer', 'top');
end

save_figure(fig_diameter, dirs.save_dir, 'fig_fwhm_statediameter');

function save_figure(fig, save_dir, fig_name)
    % R2023b exportgraphics cannot write svg, so the vector copy goes through
    % print, which is what make_fig.save2svg does. see CLAUDE.md
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    print(fig, fullfile(save_dir, fig_name + ".svg"), '-dsvg', '-vector');
    exportgraphics(fig, fullfile(save_dir, fig_name + ".png"), 'Resolution', 200);
    fprintf('wrote %s.svg and .png to %s\n', fig_name, save_dir);
end
