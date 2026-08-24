function keys = util_sessionkey(row_table)

%   IN   row_table  N-row table carrying the four columns below
%   OUT  keys       Nx1 str   "hql090|251012|sleep|005"

key_columns = ["MouseID", "Date", "SessionType", "SessionID"];

missing = setdiff(key_columns, string(row_table.Properties.VariableNames));
if ~isempty(missing)
    error('util_sessionkey:noColumn', ...
        'the table has no %s column', strjoin(missing, ', '));
end

n_row = height(row_table);
parts = strings(n_row, numel(key_columns));
for c = 1:numel(key_columns)
    parts(:, c) = string(row_table.(key_columns(c)));
end
keys = strings(n_row, 1);
for k = 1:n_row
    keys(k) = strjoin(parts(k, :), "|");
end
end
