function blown_table = explode_nest(struct_array, key_table)
    % Explode nested struct array and bring things into table
    %
    % Inputs:
    %   struct_array - Nx1 struct array where each element contains
    %                             transition tables (e.g., thickness_bv, thickness_pvs)
    %   key_table - Table with N rows containing metadata (MouseID, Date, VesselID, Depth)
    %
    % Output:
    %   blown_table - Long-format table with all transition data and metadata

    
    % Get field names (e.g., 'thickness_bv', 'thickness_pvs', etc.)
    struct_fnames = fieldnames(struct_array(1));
    
    % Initialize cell array to collect all rows
    aftermath_rows = {};
    
    % Loop through each session
    for session_idx = 1:numel(struct_array)
        % Get metadata for this session
        session_keyrow = key_table(session_idx, :);
        
        % Loop through each transition field
        for field_idx = 1:numel(struct_fnames)
            field_name = struct_fnames{field_idx};
            target_cell = struct_array(session_idx).(field_name);
            
            % Convert to table if struct
            if isstruct(target_cell)
                fragment_table = struct2table(target_cell);
            else istable(target_cell);
                
                fragment_table = target_cell;
            end
            
            % Add metadata columns to each row
            n_rows = height(fragment_table);
            fragment_table.SessionIndex = repmat(session_idx, n_rows, 1);
            fragment_table.DataType = repmat(string(field_name), n_rows, 1);
            
            % Add key columns
            for var = session_keyrow.Properties.VariableNames
                var_name = var{1};
                fragment_table.(var_name) = repmat(session_keyrow.(var_name), n_rows, 1);
            end
            
            % Reorder: move key columns to the beginning
            n_keys = length(session_keyrow.Properties.VariableNames);
            fragment_table = fragment_table(:, circshift(fragment_table.Properties.VariableNames, n_keys, 2));

            % Collect this table
            aftermath_rows{end+1} = fragment_table;
        end
    end
    
    % Concatenate all rows
    blown_table = vertcat(aftermath_rows{:});
end
