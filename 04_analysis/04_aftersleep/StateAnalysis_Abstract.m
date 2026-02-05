classdef (Abstract) StateAnalysis_Abstract < handle
    %STATEANALYSIS_ABSTRACT Base class for sleep state analysis.
    %   Coordinated by sleep_integration object.

    properties
        sleep_obj   % Handle to sleep_integration object (Temporal Source of Truth)
        t_axis      % Time axis for the specific data being analyzed
        data        % The raw data to analyze (Dimension varies by subclass)
        results     % Struct to store analysis results
        name        % Name of this analysis instance
        fs          % Sampling frequency
    end

    methods
        function obj = StateAnalysis_Abstract(sleep_obj, data, t_axis, name)
            % Constructor
            obj.sleep_obj = sleep_obj;
            obj.data = data;
            obj.t_axis = t_axis;
            obj.name = name;

            % Calculate FS if not provided
            if length(t_axis) > 1
                obj.fs = 1 / mean(diff(t_axis));
            else
                obj.fs = 1;
            end

            obj.results = struct();
        end

        function obj = get_state_indices(obj, t_axis,fs)
            obj.sleep_obj.add_taxis(t_axis,fs)
        end
    end

    methods (Abstract)
        % Subclasses must implement the actual analysis logic
        run_analysis(obj)
        duration_filter(obj,sec)
        save(obj)
    end
end
