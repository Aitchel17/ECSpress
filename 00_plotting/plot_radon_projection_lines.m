function plot_radon_projection_lines(ax, center, radius, angles, color)
%PLOT_RADON_PROJECTION_LINES Plots lines corresponding to Radon projection angles.
%   ax: Axes handle to plot on.
%   center: [y, x] coordinates of the image center (1-based).
%   radius: Radius of the lines to draw.
%   angles: Vector of angles in degrees (0 to 180).
%   color: Color of the lines (char or [r g b]).

if nargin < 5
    color = 'r';
end

hold(ax, 'on');
for theta = angles
    % Radon theta 0 corresponds to Vertical Integration (Horizontal Projection Axis)
    % theta 90 corresponds to Horizontal Integration (Vertical Projection Axis)
    % We plot the line passing through center at geometric angle theta.
    % Image Y-axis is down, so geometric angle uses -theta.

    rad_theta = deg2rad(-theta);

    x = center(2) + radius * cos(rad_theta) * [-1, 1];
    y = center(1) + radius * sin(rad_theta) * [-1, 1];

    plot(ax, x, y, 'Color', color, 'LineStyle', '-');

    % Label position (at the positive end)
    text(ax, x(2), y(2), sprintf('%d', theta), ...
        'Color', 'y', 'FontSize', 8, 'HorizontalAlignment', 'center');
end
end
