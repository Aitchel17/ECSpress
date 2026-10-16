function idx_table = timetable2frame(t_axis, timetable)
%TIMETABLE2FRAME  Windows in seconds as [first last] frames on one time axis, half-open.
%   Caller: state_1d.get_state_indices, analysis_main_piv, piv_integration_testbed, stateintegration_test
%   A frame is inside a window when start <= t < end, so two bouts that meet at one
%   boundary time share no frame. A window that holds no frame of t_axis -- a bout the
%   scoring ran past the recording -- is dropped.
%
% IN   t_axis     1 x T float s   frame times
%      timetable  N x 2 float s   [start end] per window
% OUT  idx_table  M x 2 int       [first last] frame of each window that holds one, M <= N
    n_window = size(timetable, 1);
    idx_table = zeros(0, 2);
    for window_idx = 1:n_window
        % 0. make bool array
        in_window = t_axis >= timetable(window_idx, 1) & t_axis < timetable(window_idx, 2);
        % 1. find first true
        first_frame = find(in_window, 1, 'first');
        if isempty(first_frame) % if first true not exist pass
            continue
        end
        % 2. find last true
        last_frame = find(in_window, 1, 'last');
        % 3. add to idx table
        idx_table(end+1, :) = [first_frame, last_frame]; %#ok<AGROW>
    end
end
