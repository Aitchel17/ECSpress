function blown_table = explode_nest(nest_cell, key_table)
%EXPLODE_NEST  One nest, one session per cell, flattened into a long table.
%   Each session's nest is a struct whose FIELDS are the nine series the state
%   analysis was run over (thickness_bv, displacement_downpvs, ...) and whose
%   VALUES are that series' result table. This drops the field name into a
%   DataType column and stacks every session's every series into one table.
%
%   IN   nest_cell  Nx1 cell   nest_cell{k} = the struct for session k
%        key_table  N-row table, one row per session, the columns that say
%                   WHICH session it is. Whatever they are, they are copied onto
%                   every row this makes
%   OUT  blown_table  long table: key columns, DataType, then the series' own
%
%   The identity of a row is the key columns and nothing else. An earlier version
%   also stamped the loop counter as SessionIndex, and a reader downstream then
%   had to index a second table with it to find out what the session WAS.
%   see CLAUDE_LOG.md

    n_session = numel(nest_cell);
    if height(key_table) ~= n_session
        error('explode_nest:keyMismatch', ...
            '%d nests but %d key rows', n_session, height(key_table));
    end
    key_names = key_table.Properties.VariableNames;
    n_keys = numel(key_names);

    fragment_list = {};
    for session_idx = 1:n_session
        session_nest = nest_cell{session_idx};
        if isempty(session_nest)
            continue
        end
        session_keyrow = key_table(session_idx, :);
        % per session, not off the first one: a session saved before the state
        % analysis gained a series has fewer fields, and reading the first
        % session's list would ask every other session for a field it may not
        % have. see CLAUDE_LOG.md
        series_names = fieldnames(session_nest);

        for field_idx = 1:numel(series_names)
            series_name = series_names{field_idx};
            series_content = session_nest.(series_name);
            if isempty(series_content)
                continue
            end
            if isstruct(series_content)
                fragment_table = struct2table(series_content);
            else
                fragment_table = series_content;
            end

            n_rows = height(fragment_table);
            fragment_table.DataType = repmat(string(series_name), n_rows, 1);
            for key_idx = 1:n_keys
                key_name = key_names{key_idx};
                fragment_table.(key_name) = repmat(session_keyrow.(key_name), n_rows, 1);
            end

            % the keys were appended, so rotate them to the front
            column_order = circshift(fragment_table.Properties.VariableNames, n_keys, 2);
            fragment_table = fragment_table(:, column_order);
            fragment_list{end+1} = fragment_table; %#ok<AGROW>
        end
    end

    blown_table = vertcat(fragment_list{:});
end
