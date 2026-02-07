function idx_table = timetable2frame(t_axis,timetable)
%TIME2FRAME Convert time windows to frame indices
%   idx_table = timearr2frame(t_axis, timetable)
%   timetable: [n x 2] matrix of start and end times
%   t_axis: time axis vector
%
%   This function ensures that all returned index windows have the same length,
%   determined by the duration of the first event in timetable.
%   Events that fall outside the bounds of t_axis are excluded.

% Calculate sampling frequency
fs = 1 / mean(diff(t_axis));
n_time = size(timetable,1);



% Collect valid indices
idx_list = [];

for eventidx = 1:n_time
    % Find closest start index
    [~, start_idx] = min(abs(t_axis - timetable(eventidx,1)));
    [~, end_idx] = min(abs(t_axis - timetable(eventidx,2)));
    idx_list = [idx_list; start_idx, end_idx];
end

idx_table = idx_list;
end
