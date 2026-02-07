function fig = plot_transition(state_obj, field_name, state_name, options)
% PLOT_TRANSITION Plot before/after statistics from transition table
%
% Inputs:
%   state_obj   - state_linefwhm object with transition data
%   field_name  - Name of the field in transition property (e.g., 'BV_thickness')
%   state_name  - State to filter (e.g., 'WAKE_to_NREM_trans')
%   options     - (optional) struct with plotting options
%                 .plot_type: 'mean_median', 'individual', 'both' (default: 'both')
%                 .show_errorbars: true/false (default: true)
%
% Example:
%   fig = plot_transition(paxfwhm_state, 'BV_thickness', 'NREM_to_REM_trans');

if nargin < 4
    options = struct();
end

% Default options
if ~isfield(options, 'plot_type')
    options.plot_type = 'both';
end
if ~isfield(options, 'show_errorbars')
    options.show_errorbars = true;
end

% Get filtered transition table
trans_table = state_obj.get_filtered_table('transition', field_name, state_name);

if isempty(trans_table)
    warning('No data found for state: %s', state_name);
    fig = [];
    return;
end

nbouts = height(trans_table);

% Create figure
fig = figure('Position', [100, 100, 1200, 800]);

switch options.plot_type
    case 'mean_median'
        % Single plot showing mean/median comparison
        plot_summary_stats(trans_table, state_name, field_name, options);

    case 'individual'
        % Multiple subplots for individual bouts
        plot_individual_bouts(trans_table, state_name, field_name);

    case 'both'
        % Combined view
        subplot(2, 1, 1)
        plot_summary_stats(trans_table, state_name, field_name, options);

        subplot(2, 1, 2)
        plot_bout_changes(trans_table, state_name, field_name);
end

sgtitle(sprintf('%s - %s Transition Analysis', field_name, state_name), ...
    'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
end

function plot_summary_stats(trans_table, ~, ~, options)
% Plot aggregate before/after statistics

% Extract before and after data
before_means = trans_table.before_mean;
after_means = trans_table.after_mean;
before_medians = trans_table.before_median;
after_medians = trans_table.after_median;

% Calculate statistics across all bouts
stats_data = [
    mean(before_means, 'omitnan'), mean(after_means, 'omitnan');
    median(before_medians, 'omitnan'), median(after_medians, 'omitnan')
    ];

stats_err = [
    std(before_means, 'omitnan'), std(after_means, 'omitnan');
    std(before_medians, 'omitnan'), std(after_medians, 'omitnan')
    ];

% Bar plot
x_pos = [1, 2; 3.5, 4.5];
colors_before = [0.3, 0.6, 0.9];  % Blue
colors_after = [0.9, 0.4, 0.3];   % Red

hold on;
for i = 1:2  % Mean, Median
    bar(x_pos(i, 1), stats_data(i, 1), 0.4, 'FaceColor', colors_before);
    bar(x_pos(i, 2), stats_data(i, 2), 0.4, 'FaceColor', colors_after);

    if options.show_errorbars
        errorbar(x_pos(i, 1), stats_data(i, 1), stats_err(i, 1), ...
            'k.', 'LineWidth', 1.5, 'CapSize', 10);
        errorbar(x_pos(i, 2), stats_data(i, 2), stats_err(i, 2), ...
            'k.', 'LineWidth', 1.5, 'CapSize', 10);
    end
end
hold off;

% Formatting
xticks([mean(x_pos(1, :)), mean(x_pos(2, :))]);
xticklabels({'Mean', 'Median'});
ylabel('Value');
title(sprintf('Before vs After Transition (n=%d bouts)', height(trans_table)));
legend({'Before', 'After'}, 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);
end

function plot_bout_changes(trans_table, ~, ~)
% Plot individual bout changes (before -> after)

before_means = trans_table.before_mean;
after_means = trans_table.after_mean;
nbouts = length(before_means);

% Calculate change
delta = after_means - before_means;
percent_change = (delta ./ before_means) * 100;

% Create paired plot
hold on;
for i = 1:nbouts
    plot([1, 2], [before_means(i), after_means(i)], ...
        'o-', 'Color', [0.5, 0.5, 0.5, 0.4], 'LineWidth', 1, 'MarkerSize', 6);
end

% Add mean line
plot([1, 2], [mean(before_means, 'omitnan'), mean(after_means, 'omitnan')], ...
    'ro-', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;

xlim([0.5, 2.5]);
xticks([1, 2]);
xticklabels({'Before', 'After'});
ylabel('Mean Value');
title(sprintf('Individual Bout Changes (Δ = %.3f ± %.3f)', ...
    mean(delta, 'omitnan'), std(delta, 'omitnan')));
grid on;
set(gca, 'FontSize', 11);

% Add text annotation
text(2.2, mean(after_means, 'omitnan'), ...
    sprintf('Δ%%: %.1f%%', mean(percent_change, 'omitnan')), ...
    'FontSize', 10, 'FontWeight', 'bold');
end

function plot_individual_bouts(trans_table, ~, ~)
% Plot each bout's raw data before and after transition

nbouts = height(trans_table);
ncols = min(4, nbouts);
nrows = ceil(nbouts / ncols);

for i = 1:nbouts
    subplot(nrows, ncols, i);

    before_data = trans_table.before_data{i};
    after_data = trans_table.after_data{i};

    % Create time axis
    fs = 1 / trans_table.before_duration(i) * length(before_data);
    t_before = (0:length(before_data)-1) / fs;
    t_after = (0:length(after_data)-1) / fs + t_before(end);

    hold on;
    plot(t_before, before_data, 'b-', 'LineWidth', 1.5);
    plot(t_after, after_data, 'r-', 'LineWidth', 1.5);

    % Mark transition point
    xline(t_before(end), 'k--', 'LineWidth', 2, 'Label', 'Transition');
    hold off;

    title(sprintf('Bout %d', i));
    xlabel('Time (s)');
    ylabel('Value');
    grid on;

    if i == 1
        legend({'Before', 'After'}, 'Location', 'best');
    end
end
end
