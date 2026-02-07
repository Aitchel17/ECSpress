%% Example: Plot Transition Analysis
% This script demonstrates how to use plot_transition to visualize
% before/after transition statistics from state_linefwhm objects

% Assuming you have already run the transition analysis:
% paxfwhm_state.get_transitionsummary('BV_thickness', session.pax_fwhm.thickness.bv)

%% Example 1: Default plot (both summary and individual changes)
fig1 = plot_transition(paxfwhm_state, 'BV_thickness', 'NREM_to_REM_trans');

%% Example 2: Summary statistics only
options = struct();
options.plot_type = 'mean_median';
options.show_errorbars = true;

fig2 = plot_transition(paxfwhm_state, 'PVStotal_thickness', 'WAKE_to_NREM_trans', options);

%% Example 3: Individual bout changes only
options.plot_type = 'both';
fig3 = plot_transition(paxfwhm_state, 'PVSdynamic_thickness', 'REM_to_WAKE_trans', options);

%% Example 4: Plot all available transitions for a field
% Get all transition states
field_name = 'BV_thickness';
trans_table_full = paxfwhm_state.transition.(field_name);
unique_states = unique(trans_table_full.state_name);

% Plot each transition type
for i = 1:length(unique_states)
    state_name = char(unique_states(i));
    fig = plot_transition(paxfwhm_state, field_name, state_name);

    % Optional: save figure
    % saveas(fig, fullfile(savepath, sprintf('%s_%s.png', field_name, state_name)));
end

%% Example 5: Comparing multiple parameters for the same transition
transition_state = 'NREM_to_REM_trans';
parameters = {'BV_thickness', 'PVStotal_thickness', 'PVSdynamic_thickness', 'Extraparenchyma_thickness'};

figure('Position', [50, 50, 1600, 1000]);
for i = 1:length(parameters)
    subplot(2, 2, i);

    % Get filtered data
    trans_table = paxfwhm_state.get_filtered_table('transition', parameters{i}, transition_state);

    if ~isempty(trans_table)
        % Plot before vs after means
        before_means = trans_table.before_mean;
        after_means = trans_table.after_mean;

        hold on;
        for j = 1:length(before_means)
            plot([1, 2], [before_means(j), after_means(j)], ...
                'o-', 'Color', [0.5, 0.5, 0.5, 0.4], 'LineWidth', 1);
        end

        % Mean line
        plot([1, 2], [mean(before_means, 'omitnan'), mean(after_means, 'omitnan')], ...
            'ro-', 'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        hold off;

        xlim([0.5, 2.5]);
        xticks([1, 2]);
        xticklabels({'Before', 'After'});
        ylabel('Mean Value');
        title(parameters{i}, 'Interpreter', 'none');
        grid on;
    end
end
sgtitle(sprintf('Transition Comparison: %s', transition_state), ...
    'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
