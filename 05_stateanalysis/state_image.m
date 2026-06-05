classdef state_image < handle
    %STATEANALYSIS_ABSTRACT Base class for sleep state analysis.
    %   Coordinated by sleep_integration object.

    properties
        sleep_obj   % Handle to sleep_integration object (Temporal Source of Truth)
        t_axis      % Time axis for the specific data being analyzed
        state_idx
        param        % Name of this analysis instance
        imgch1
        imgch2
    end



    methods
        function obj = state_image(sleep_obj)
            % Constructor
            obj.sleep_obj = sleep_obj;
        end

        function obj = get_state_indices(obj, t_axis,fs)
            % from
            obj.state_idx =  obj.sleep_obj.add_taxis(t_axis);
            obj.t_axis = t_axis;
            obj.param.fs = fs;
        end
    end
end