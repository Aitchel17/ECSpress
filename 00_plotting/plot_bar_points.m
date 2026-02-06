function plot_bar_points(data_cell, labels, varargin)
%PLOT_BAR_POINTS Create a bar plot with individual data points and SEM
%   plot_bar_points({data1, data2}, {'Group1', 'Group2'})
%
%   Inputs:
%       data_cell: Cell array of vectors (e.g., {vecA, vecB})
%       labels: Cell array of strings for x-axis labels
%
%   Optional Name-Value Parameters:
%       'Colors', matrix : Nx3 RGB matrix for bar colors
%       'YLabel', str    : Label for Y-axis
%       'Title', str     : Plot title
%       'Ax', axes handle: Target axes

p = inputParser;
addRequired(p, 'data_cell', @iscell);
addRequired(p, 'labels', @iscell);
addParameter(p, 'Colors', lines(length(data_cell)), @isnumeric);
addParameter(p, 'YLabel', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Title', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Ax', [], @(x) isempty(x) || isgraphics(x, 'axes'));

parse(p, data_cell, labels, varargin{:});

colors = p.Results.Colors;

if isempty(p.Results.Ax)
    figure('Color', 'w');
    ax = gca;
else
    ax = p.Results.Ax;
end

hold(ax, 'on');

nGroups = length(data_cell);

for i = 1:nGroups
    dat = data_cell{i};
    dat = dat(~isnan(dat)); % Remove NaNs

    if isempty(dat)
        continue;
    end

    mu = mean(dat);

    % Bar
    bar(ax, i, mu, 'FaceColor', colors(i,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'BarWidth', 0.6);

    % Error Bar
    errorbar(ax, i, mu, std(dat), 'k', 'LineWidth', 1.5, 'CapSize', 10);

    % Scatter Points (Jittered)
    x_scatter = i + (rand(size(dat)) - 0.5) * 0.25;
    scatter(ax, x_scatter, dat, 25, 'k', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
end

set(ax, 'XTick', 1:nGroups, 'XTickLabel', labels, 'FontSize', 12);
ylabel(ax, p.Results.YLabel, 'FontSize', 12, 'Interpreter', 'none');
title(ax, p.Results.Title, 'FontSize', 14, 'Interpreter', 'none');
grid(ax, 'on');
box(ax, 'off');
end
