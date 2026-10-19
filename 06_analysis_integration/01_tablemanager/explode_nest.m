function blown_table = explode_nest(session_nest, key_row)
%EXPLODE_NEST  Convert nested struct of a session into a long table,
%   Caller: bout_paxfwhm, bout_drift
%
%   IN   session_nest  1x1 struct    one session's nested data column
%        key_row       1-row table   keys for the session
%                             
%   OUT  blown_table   long table: key columns, the series' own, then DataType

    % 0. Check nested struct to be empty
    if isempty(session_nest)
        blown_table = [];   % no series filled this nest; the caller's is_filled drops it
        return
    end

    % 1. Find key and data names
    key_names = key_row.Properties.VariableNames;
    n_keys = numel(key_names);
    series_names = fieldnames(session_nest);

    % 2. Explode series
    fragment_list = {};
    for field_idx = 1:numel(series_names)
        series_name = series_names{field_idx};
        fragment_table = session_nest.(series_name);
        % 0. check data exist
        if isempty(fragment_table)
            continue
        end
        % 1. Count number of duplication (n bouts)
        n_rows = height(fragment_table);
        % 2. Append Data type (eg. pvs_thickness ...)
        fragment_table.DataType = repmat(string(series_name), n_rows, 1);
        % 3. Append Keys
        for key_idx = 1:n_keys
            key_name = key_names{key_idx};
            fragment_table.(key_name) = repmat(key_row.(key_name), n_rows, 1);
        end

        % 4. Bring keys to front
        column_order = circshift(fragment_table.Properties.VariableNames, n_keys, 2);
        fragment_table = fragment_table(:, column_order);
        fragment_list{end+1} = fragment_table; %#ok<AGROW>
    end

    blown_table = vertcat(fragment_list{:});
end
