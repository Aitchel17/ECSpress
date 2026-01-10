function [fig, figAxes] = analog_plot(analog)
    %ANALOG_PLOT Plots the analog signals and returns the figure and axes handles.

    % Make figure
    fig = figure('Name', 'Analog signal', 'NumberTitle', 'off');
    % Define the number of subplots
    num_subplots = length(fieldnames(analog)) / 2;
    current_subplot = 0; % Counter for subplots
    
    % Initialize axes handle storage
    figAxes = gobjects(0); % Create empty graphics object array

    % Ball plot
    if isfield(analog, 'raw_Ball')
        current_subplot = current_subplot + 1;
        ax = subplot(num_subplots, 1, current_subplot);
        plot(analog.ds_ball(1, :) / 60, analog.ds_ball(2, :));
        ylabel('Ball (V)');
        xlim([0 30]);
        xlabel('min');
        set(ax, 'FontSize', 14);
        figAxes(end + 1) = ax; % Append current axes to figAxes
    end

    % EMG plot
    if isfield(analog, 'raw_EMG')
        current_subplot = current_subplot + 1;
        ax = subplot(num_subplots, 1, current_subplot);
        plot(analog.ds_EMG(1, :) / 60, analog.ds_EMG(2, :));
        xlim([0 30]);
        ylabel('EMG (V)');
        xlabel('min');
        set(ax, 'FontSize', 14);
        figAxes(end + 1) = ax; % Append current axes to figAxes
    end

    % ECoG plot
    if isfield(analog, 'raw_ECoG')
        current_subplot = current_subplot + 1;
        ax = subplot(num_subplots, 1, current_subplot);
        surface(ax, analog.ecog_spectrum.t_axis / 60, analog.ecog_spectrum.f_axis, ...
            zeros(size(analog.ecog_spectrum.log_norm_spectrum)), ...
            analog.ecog_spectrum.log_norm_spectrum, 'LineStyle', 'none');
        set(ax, 'YScale', 'log');
        xlim([0 30]);
        ylabel('Freq (Hz)');
        xlabel('min');
        set(ax, 'FontSize', 14);
        figAxes(end + 1) = ax; % Append current axes to figAxes
    end
end
