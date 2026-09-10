function subtable = flatten_state(central, key_columns, nest_names)
%FLATTEN_STATE  A centralized state table's nests, flattened into one long table each.
%   Caller: tablegeneration_main, tablegeneration_drift
%   Each row of `central` holds one session's state analysis in a `data` struct whose
%   fields are the nests (state_summary, transition, ...). This pulls one nest out of
%   every session and stacks them, once per nest. session_key is stamped here, once
%   per session, rather than rebuilt per flattened row by each caller.
%
%   IN   central      N-row table   centralized_<x>_state, with a `data` cell column
%        key_columns  1 x K str     which of its columns say WHICH session a row is
%        nest_names   1 x M str     the nests to flatten
%   OUT  subtable     1x1 struct    one field per nest, each a long table:
%                                   key columns, session_key, DataType, then the series'

    key_table = central(:, key_columns);
    key_table.session_key = util_sessionkey(central);
    n_session = height(central);

    subtable = struct();
    for name_idx = 1:numel(nest_names)
        nest_name = nest_names(name_idx);
        nest_cell = cell(n_session, 1);
        for session_idx = 1:n_session
            session_row = central.data{session_idx};
            if isfield(session_row, nest_name)
                nest_cell{session_idx} = session_row.(nest_name);
            end
        end
        subtable.(nest_name) = explode_nest(nest_cell, key_table);
    end
end
