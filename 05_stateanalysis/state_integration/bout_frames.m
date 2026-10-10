function state_idx = bout_frames(t_axis, time_table)
%BOUT_FRAMES  Every bout table of a time_table as frame windows on one time axis.
%   Caller: state_1d.get_state_indices, analysis_main_piv, piv_integration_testbed, stateintegration_test
%
% IN   t_axis      1 x T float s   the measurement's time axis
%      time_table  1x1 struct      table name -> N x 2 [start end] bouts in seconds,
%                                  as bout_sleep or bout_puffmove builds it
% OUT  state_idx   1x1 struct      the same names -> M x 2 int frame windows, as
%                                  timetable2frame cuts them on t_axis
    state_idx = struct();
    table_names = fieldnames(time_table);
    for k = 1:numel(table_names)
        table_name = table_names{k};
        state_idx.(table_name) = timetable2frame(t_axis, time_table.(table_name));
    end
end
