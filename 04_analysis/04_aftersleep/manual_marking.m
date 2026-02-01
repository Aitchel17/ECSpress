fprintf('--- Interactive Marker Mode ---\n');
fprintf('Click on any image to mark position on all axes.\n');
fprintf('Press Enter (key) to stop.\n');

marker_handles = [];
while true
    try
        [x_val, y_val, button] = ginput(1);
    catch
        break; % Handle figure closure
    end

    % Break if Enter (empty button) or any key press (button > 3 typically)
    % Mouse clicks are 1 (Left), 2 (Middle), 3 (Right)
    if isempty(button) || ~ismember(button, [1, 2, 3])
        fprintf('Exiting marker mode.\n');
        break;
    end

    % Remove old markers
    delete(marker_handles);
    marker_handles = gobjects(0);

    % Add new markers to all axes
    for i = 1:numel(all_axes)
        try
            curr_ax = all_axes(i);
            hold(curr_ax, 'on');
            h = plot(curr_ax, x_val, y_val, 'rx', 'MarkerSize', 8, 'LineWidth', 1);
            marker_handles(end+1) = h;
        catch
            % Handle closed axes
        end
    end
    drawnow;
end