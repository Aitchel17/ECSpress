classdef Vessel < handle
    % VESSEL Represents a single vessel across multiple imaging sessions.
    % Aggregates data from Primary, State, and other analyses for longitudinal averaging.
    
    properties
        VesselID        % String Identifier (e.g., 'V1')
        MouseID         % String Identifier (e.g., 'M1')
        Depth           % Numeric Depth (um)
        DepthLayer      % 'L1' or 'L2'
        
        % Aggregated Data Storage (Cell arrays of structs/tables)
        SessionIDs      % Array of session IDs
        Dates           % Array of dates
        
        PrimaryData     % struct array of primary analysis structs
        StateData       % struct array of state analysis structs
        
        % Computed Averages
        AverageData     % Struct containing averaged metrics
    end
    
    methods
        function obj = Vessel(vessel_id, mouse_id, depth_layer)
            % Constructor
            obj.VesselID = string(vessel_id);
            obj.MouseID = string(mouse_id);
            obj.DepthLayer = string(depth_layer);
            % Initialize storage
            obj.SessionIDs = [];
            obj.Dates = [];
            obj.PrimaryData = [];
            obj.StateData = [];
            obj.AverageData = struct();
        end
        
        function addSession(obj, session_id, date, depth,primary_struct, state_struct)
            % Adds data from a single session to this vessel
            obj.SessionIDs(end+1) = string(session_id);
            obj.Dates(end+1) = string(date);
            obj.Depth(end+1)  = depth;
            
            % Add Primary Data (if exists)
            if isempty(obj.PrimaryData)
                obj.PrimaryData = primary_struct;
            else
                obj.PrimaryData(end+1) = primary_struct;
            end
            
            % Add State Data (if exists)
            if isempty(obj.StateData)
                obj.StateData = state_struct;
            else
                obj.StateData(end+1) = state_struct;
            end


        end
        
        function count = getSessionCount(obj)
            count = numel(obj.SessionIDs);
        end
        
        function computeAverage(obj)
            % Aggregates and averages data across sessions
            
            if isempty(obj.StateData)
                return;
            end
            
            % Collect all transition tables from all sessions
            all_transition_tables = struct();
            
            for i = 1:numel(obj.StateData)
                session_data = obj.StateData(i);
                
                % Check if transition field exists
                if isfield(session_data, 'transition')
                    % Get all fields within transition
                    trans_fields = fieldnames(session_data.transition);
                    
                    for f = 1:numel(trans_fields)
                        field_name = trans_fields{f};
                        trans_table = session_data.transition.(field_name);
                        
                        % Convert struct to table if needed
                        if isstruct(trans_table)
                            trans_table = struct2table(trans_table);
                        end
                        
                        % Concatenate tables across sessions
                        if ~isfield(all_transition_tables, field_name)
                            all_transition_tables.(field_name) = trans_table;
                        else
                            all_transition_tables.(field_name) = [all_transition_tables.(field_name); trans_table];
                        end
                    end
                end
            end
            
            % Average each transition table using the helper function
            if ~isempty(fieldnames(all_transition_tables))
                averaged_fields = fieldnames(all_transition_tables);
                
                for f = 1:numel(averaged_fields)
                    field_name = averaged_fields{f};
                    merged_table = all_transition_tables.(field_name);
                    
                    % Apply averaging function
                    averaged_table = average_transition_table(merged_table);
                    
                    % Store in AverageData
                    obj.AverageData.(field_name) = averaged_table;
                    
                    fprintf('Vessel %s: Averaged %s (%d states, %d total traces).\n', ...
                        obj.VesselID, field_name, height(averaged_table), sum(averaged_table.n_traces));
                end
            end
        end
    end
end
