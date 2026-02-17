function joined_table = average_table(targettable, key_cols, numeric_colnames, data_colname)
%GET_VESSELAVGTABLE Summary of this function goes here
%   Detailed explanation goes here
% Identify numeric columns to average (exclude metadata)

    
    % Find numeric columns (excluding 'data' which needs special handling)
    numeric_summary = groupsummary(targettable, key_cols, "mean", numeric_colnames);



% Handle 'data' column separately (average time series)
    averaged_data = table();
    if ~isempty(data_colname) && ismember(data_colname, targettable.Properties.VariableNames)
        % Group by key_cols
        unique_keys = unique(targettable(:, key_cols));
        
        for i = 1:height(unique_keys)
            key_row = unique_keys(i, :);
            
            % Filter for this key combination
            logic = ismember(targettable(:, key_cols), key_row);
            subset = targettable(logic, :);
            
            if height(subset) > 0
                % Concatenate all data traces
                all_traces = [];
                raw_data = subset.(data_colname);
                
                for k = 1:numel(raw_data)
                    trace = raw_data{k};
                    if iscolumn(trace), trace = trace'; end
                    all_traces = [all_traces; trace];
                end
                
                % Compute mean
                mean_trace = mean(all_traces, 1, 'omitnan');
                
                % Store result
                % Start with the key columns
                row = key_row;
                row.mean_data = {mean_trace};
                row.n_traces = height(subset);
                
                averaged_data = [averaged_data; row];
            end
        end
    end
    
% Join the tables if both exist
    if ~isempty(averaged_data)
        % Perform inner join on key columns
        joined_table = innerjoin(numeric_summary, averaged_data, ...
            'Keys', key_cols);
        disp('Joined scalar summary and time-series data.');
    else
        % If no time-series data, just return the numeric summary
        joined_table = numeric_summary;
        disp('Returned scalar summary only.');
    end
end

