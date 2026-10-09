%FIG_FWHM_STATEDEVIATION  Where the lumen sits within its own session, by arousal state.
%   Reads fwhmrelation.mat and presses each state's pooled heatmap down onto its
%   diameter axis, so every frame of every recording is counted. Nothing is measured
%   here; pool_state_by_vessel already normalised each session's map to sum one and
%   averaged within vessel, so a long recording weighs the same as a short one.
%
%   The x axis is micrometres FROM THAT SESSION'S OWN MEDIAN, not absolute calibre --
%   the maps are stacked on a shared zero. A 20 um vessel and a 10 um vessel both sit
%   at 0 when they are at rest.
clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
addpath('g:/03_program/01_ecspress/09_dirstruct');    % where dirs_ecspath lives
dirs_ecspath;                                         % three roots, minus zz_notinuse
clee = color_lee();
param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged';
end
dirs = dirs_central();
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);

%% Picture
% keyed by the state's own name, so a state entering or leaving state_pool cannot
% slide the legend along by one
param.state_label = struct('awake', "Awake", 'drowsy', "Drowsy", ...
    'nrem', "NREM", 'rem', "REM");
param.fontsize = 12;
param.linewidth = 1.5;

%% What is drawn
relation_path = fullfile(dirs.save_dir, "fwhmrelation.mat");
relation = load(relation_path).save_content;
state_pool = relation.state_pool;
grid_x = relation.vessel_pool.grid_x;
n_state = height(state_pool);
fprintf('%s\n   state_pool : %s\n', relation_path, strjoin(state_pool.state', ', '));

% the columns a state did not reach are NaN in its map, and they stay out of the
% marginal rather than counting as zero occupancy
marginal = cell(n_state, 1);
for k = 1:n_state
    column_weight = sum(state_pool.map{k}, 1, 'omitnan');
    column_weight(all(isnan(state_pool.map{k}), 1)) = NaN;
    marginal{k} = column_weight;
end

%% The board
fig_deviation = figure('Color', 'w', 'Position', [60 60 640 460], ...
    'Name', 'fig_fwhm_statedeviation');
axis_deviation = axes(fig_deviation);
hold(axis_deviation, 'on')

legend_entry = strings(1, n_state);
for k = 1:n_state
    plot(axis_deviation, grid_x, marginal{k}, 'Color', clee.clist.(state_pool.state(k)), ...
        'LineWidth', param.linewidth);
    drawn = isfinite(marginal{k});
    centre_of_mass = sum(grid_x(drawn) .* marginal{k}(drawn)) / sum(marginal{k}(drawn));
    legend_entry(k) = sprintf('%s  %d vessels, mean %+.2f um', ...
        param.state_label.(state_pool.state(k)), state_pool.n_vessel(k), centre_of_mass);
    fprintf('   %-7s %2d vessels | drawn %+.1f to %+.1f um | weighted mean %+.2f\n', ...
        state_pool.state(k), state_pool.n_vessel(k), min(grid_x(drawn)), ...
        max(grid_x(drawn)), centre_of_mass);
end

xline(axis_deviation, 0, '--k', 'Alpha', 0.3, 'LineWidth', 1);
legend(axis_deviation, legend_entry, 'Location', 'northeast', 'Box', 'off', ...
    'FontSize', param.fontsize - 1);
xlabel(axis_deviation, 'lumen diameter, um from this session''s median');
ylabel(axis_deviation, 'occupancy, fraction of frames per column');
set(axis_deviation, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, ...
    'FontName', 'Arial', 'FontSize', param.fontsize, 'Layer', 'top');

save_figure(fig_deviation, dirs.save_dir, 'fig_fwhm_statedeviation');

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
