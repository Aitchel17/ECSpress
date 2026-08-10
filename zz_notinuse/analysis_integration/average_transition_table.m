function averaged_table = average_transition_table(transition_table)
    % AVERAGE_TRANSITION_TABLE Averages data traces for each unique state
    %
    % Input:
    %   transition_table - Table with columns including 'state_name' and 'data'
    %
    % Output:
    %   averaged_table - Table with one row per unique state, containing:
    %       - state_name: The state identifier
    %       - mean_data: The averaged trace (1xN double)
    %       - all_data: Cell array containing all individual traces (NxM double)
    %       - n_traces: Number of traces averaged
    %       - Other columns from the first occurrence (metadata)
    
    if isempty(transition_table)
        averaged_table = table();
        return;
    end
    
    % Get unique state names
    unique_states = unique(transition_table.state_name);
    n_states = numel(unique_states);
    
    % Initialize output struct array
    result_struct = struct();
    
    for i = 1:n_states
        state_name = unique_states(i);
        
        % Filter for this state
        state_logic = transition_table.state_name == state_name;
        state_rows = transition_table(state_logic, :);
        
        % Extract and concatenate all data traces
        all_traces = [];
        for k = 1:height(state_rows)
            trace = state_rows.data{k};
            % Ensure row vector
            if iscolumn(trace)
                trace = trace';
            end
            all_traces = [all_traces; trace];
        end
        
        % Compute mean
        if ~isempty(all_traces)
            mean_trace = mean(all_traces, 1, 'omitnan');
        else
            mean_trace = [];
        end
        
        % Preserve metadata from first occurrence
        first_row = state_rows(1, :);
        
        % Build result struct
        result_struct(i).state_name = state_name;
        result_struct(i).mean_data = mean_trace;
        result_struct(i).all_data = {all_traces};
        result_struct(i).n_traces = size(all_traces, 1);
        
        % Copy other columns (excluding 'data' which we've processed)
        other_cols = setdiff(first_row.Properties.VariableNames, {'data', 'state_name'});
        for col = other_cols
            col_name = col{1};
            result_struct(i).(col_name) = first_row.(col_name);
        end
    end
    
    % Convert to table
    averaged_table = struct2table(result_struct);
end
