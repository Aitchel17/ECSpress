function analyzer = rollup_summary(analyzer)
%ROLLUP_SUMMARY  Q2Q3 mean and percentiles per bout, then the four folds bout -> session -> vessel -> mouse.
%   Caller: tablegeneration_main, tablegeneration_drift
%   Each fold takes the table the fold before it wrote.
%
%   IN   analyzer  1x1 tableAnalyzer   filtered and scaled
%   OUT  analyzer  1x1 tableAnalyzer   numeric_tables gains filtered_table,
%                                      bout_idx_ave, Date_ave, VesselID_ave, MouseID_ave

    analyzer.meanFrom2("raw_data","Q2Q3_mean",0.25,0.75)
    analyzer.addPrctilecol("raw_data","prctile_95", 95);
    analyzer.addPrctilecol("raw_data","prctile_5",5);
    analyzer.numeric_tables.filtered_table = analyzer.filtered_table;
    analyzer.get_numericsummary("bout_idx","filtered_table")   % intra session
    analyzer.get_numericsummary("Date","bout_idx_ave")         % inter session
    analyzer.get_numericsummary("VesselID","Date_ave")         % vessel
    analyzer.get_numericsummary("MouseID","VesselID_ave")      % mouse
end
