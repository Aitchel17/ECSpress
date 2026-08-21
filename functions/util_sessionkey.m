function keys = util_sessionkey(row_table)
%UTIL_SESSIONKEY  The string that says which recording a table row is.
%   One recording is one mouse, one date, one kind of session and one session
%   number. Every stage of the secondary chain looks a row up by this and by
%   nothing else, so both halves of it -- WHICH columns, and how they are joined
%   into a string -- belong in one place. They were copied into three files and a
%   change to one of them would have made two stages disagree about which rows are
%   the same recording without either of them failing.
%
%   IN   row_table  N-row table carrying the four columns below
%   OUT  keys       Nx1 str   "hql090|251012|sleep|005"
%
%   NOT a caller's choice, which is why there is no way to pass a different set.
%
%   caution  the key does NOT separate a session folder from a renamed copy of
%        itself, and it does not separate two recordings whose folder name defeats
%        the parser: HQL072_250325_005_PA01_1 and HQL072_250325_007_PA03_1 both
%        come out as hql072|250325|unknown|1 because the parser takes the trailing
%        _1 as the session number. merge_dirtable reduces a repeated key to the row
%        whose folder still exists and NAMES any key that still points at more than
%        one live folder rather than dropping either. see CLAUDE_LOG.md
%   note  Directory is the only thing that is unique per row, which is why
%        merge_dirtable dedups on it first and on this second

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
