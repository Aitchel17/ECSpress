classdef sleep_integration
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %   This class is designed tobe initialize time_table and hold indices
    %   corresponding to time table using time_axis of target analysis
    %   Each
    properties
        sleep_score
        binary_bin
        time_table
        param
        info         
    end

    methods
        function obj = sleep_integration(sleep_score)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.param.transition_window = 25;
            obj.param.bigchunk_windowsize = 300;
            obj.sleep_score = sleep_score;
            obj.param.bigchunk_awake_weight = struct('w_awake', 10, 'w_nrem', -1, 'w_rem', -1000, 'w_drowsy', 0);
            obj.param.bigchunk_nrem_weight = struct('w_awake', -1, 'w_nrem', 10, 'w_rem', -1000, 'w_drowsy', 0);
            obj.param.bigchunk_rem_weight = struct('w_awake', -1, 'w_nrem', 0, 'w_rem', 10, 'w_drowsy', 0);
            obj.param.bigchunk_drowsy_weight = struct('w_awake', 0, 'w_nrem', -1, 'w_rem', -1000, 'w_drowsy', 10);

            % Binary array generation
            obj.binary_bin.awake = sleep_score.behavState == sleep_score.statecodes.NotSleep;
            obj.binary_bin.nrem = sleep_score.behavState == sleep_score.statecodes.NREM;
            obj.binary_bin.rem = sleep_score.behavState == sleep_score.statecodes.REM;
            obj.binary_bin.drowsy = sleep_score.behavState == sleep_score.statecodes.Drowsy;
            [obj.time_table.long_nrem, obj.info.bigchunk_nrem_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_nrem_weight,...
                                        obj.param.bigchunk_windowsize, sleep_score.binwidth_sec);
            [obj.time_table.long_awake, obj.info.bigchunk_awake_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_awake_weight,...
                                        obj.param.bigchunk_windowsize, sleep_score.binwidth_sec);
            [obj.time_table.long_drowsy, obj.info.bigchunk_drowsy_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_drowsy_weight,...
                                        obj.param.bigchunk_windowsize, sleep_score.binwidth_sec);
            [obj.time_table.long_rem, obj.info.bigchunk_rem_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_rem_weight,...
                                        obj.param.bigchunk_windowsize, sleep_score.binwidth_sec);
            % Time table generation
            % for chunk analysis
            obj.time_table.roughawake = statebin2timetable(obj.sleep_score.AwakeTimes, obj.sleep_score.DrowsyTimes);
            obj.time_table.roughnrem = statebin2timetable(obj.sleep_score.NREMTimes, obj.sleep_score.uArousalTimes);
            obj.time_table.rem = statebin2timetable(obj.sleep_score.REMTimes);
            % not actively using
            obj.time_table.awake = statebin2timetable(obj.sleep_score.AwakeTimes);
            obj.time_table.drowsy = statebin2timetable(obj.sleep_score.DrowsyTimes);
            obj.time_table.uarousal = statebin2timetable(obj.sleep_score.uArousalTimes);
            obj.time_table.nrem = statebin2timetable(obj.sleep_score.NREMTimes);
            % transition
            obj.time_table.nr_trans = get_transition(obj.param.transition_window,obj.time_table.roughnrem, obj.time_table.rem);
            obj.time_table.na_trans = get_transition(obj.param.transition_window,obj.time_table.roughnrem, obj.time_table.roughawake);
            obj.time_table.ra_trans = get_transition(obj.param.transition_window,obj.time_table.rem, obj.time_table.roughawake);
            obj.time_table.an_trans = get_transition(obj.param.transition_window,obj.time_table.roughawake, obj.time_table.roughnrem);
        end

        function state_idx = add_taxis(obj,t_axis)
            state_idx = struct();
            ttable_names = fieldnames(obj.time_table);
            for tname_idx = 1:numel(ttable_names)
                ttable_name = ttable_names{tname_idx};
                disp(ttable_name)
                state_idx.(ttable_name) = timetable2frame(t_axis,obj.time_table.(ttable_name));
            end
        end
    end
end

