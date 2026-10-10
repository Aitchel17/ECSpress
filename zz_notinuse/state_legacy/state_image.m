classdef state_image < handle
    %STATE_IMAGE  Sleep-state indices on an image stack's frame axis.
    %   Caller: analysis_main_piv, piv_integration_testbed, stateintegration_test
    %   Built on a state_integration's bout tables.

    properties
        time_table  % 1x1 struct   table name -> N x 2 [start end] bouts in seconds, state_integration's
        t_axis      % Time axis for the specific data being analyzed
        state_idx
        param        % Name of this analysis instance
        imgch1
        imgch2
    end



    methods
        function obj = state_image(bout_table)
            % bout_table  1x1 struct   state_integration's time_table
            obj.time_table = bout_table;
        end

        function obj = get_state_indices(obj, t_axis,fs)
            obj.state_idx = bout_frames(t_axis, obj.time_table);   % every bout table on the stack's frames
            obj.t_axis = t_axis;
            obj.param.fs = fs;
        end
    end
end