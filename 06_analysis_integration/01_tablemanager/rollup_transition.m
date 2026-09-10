function analyzer = rollup_transition(analyzer)
%ROLLUP_TRANSITION  The four folds bout -> session -> vessel -> mouse, for the numbers and the traces.
%   Caller: tablegeneration_main, tablegeneration_drift
%   match_tracelength resamples the traces to one length before they are averaged.
%
%   IN   analyzer  1x1 tableAnalyzer   filtered and scaled
%   OUT  analyzer  1x1 tableAnalyzer   numeric_tables and data_tables both gain
%                                      filtered_table, bout_idx_ave, Date_ave,
%                                      VesselID_ave, MouseID_ave

    analyzer.numeric_tables.filtered_table = analyzer.filtered_table;
    analyzer.get_numericsummary("bout_idx","filtered_table")   % intra session
    analyzer.get_numericsummary("Date","bout_idx_ave")         % inter session
    analyzer.get_numericsummary("VesselID","Date_ave")         % vessel
    analyzer.get_numericsummary("MouseID","VesselID_ave")      % mouse
    analyzer.get_datasummary("bout_idx","filtered_table")
    analyzer.get_datasummary("Date","bout_idx_ave")
    analyzer.get_datasummary("VesselID","Date_ave")
    analyzer.get_datasummary("MouseID","VesselID_ave")
end
