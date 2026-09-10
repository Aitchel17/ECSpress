function baseline = extract_awakebaseline(bout_idx_ave, key_cols)
%EXTRACT_AWAKEBASELINE  Each session's awake Q2Q3 mean, the level the other states are read against.
%   Caller: tablegeneration_main, tablegeneration_drift
%
%   IN   bout_idx_ave  table        one row per (session x state x DataType)
%        key_cols      1 x K str    the columns the caller will join on
%   OUT  baseline      table        key_cols, then awake_boutq2q3mean

    awake_rows = bout_idx_ave(bout_idx_ave.state_name == "awake", :);
    baseline = awake_rows(:, [key_cols, "Q2Q3_mean"]);
    baseline.Properties.VariableNames("Q2Q3_mean") = "awake_boutq2q3mean";
end
