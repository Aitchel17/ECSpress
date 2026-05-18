function timetable = logic2timetable(taxis,logic,min_duration)
%LOGIC2TIMETABLE Summary of this function goes here
%   Detailed explanation goes here
    diff_binarymove = diff([0; logic(:); 0]);
    rise_idx = find(diff_binarymove == 1);
    fall_idx = find(diff_binarymove == -1);
    % Calculate duration of each bout in seconds
    bout_durations = taxis(fall_idx(:) - 1) - taxis(rise_idx(:));
    valid_bouts = bout_durations >= min_duration;
    rise_idx = rise_idx(valid_bouts);
    fall_idx = fall_idx(valid_bouts);
    
    % Create Nx2 array of event start and end times from taxis
    % (fall_idx - 1 gets the exact last sample of the movement bout)
    timetable = [taxis(rise_idx(:))', taxis(fall_idx(:) - 1)'];
end

