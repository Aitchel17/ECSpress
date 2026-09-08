%TABLEGENERATION_DRIFT  Centralized drift state analysis -> the two tables the drift
%   figures read. Reads centralized_drift_state.mat, flattens its two nests, applies
%   the same filters as tablegeneration_main plus one -- only recordings the FWHM
%   pipeline has judged -- and writes drift_summary.mat and drift_transition.mat
%   beside the dataset, with param and filt saved beside each table.
%
%   A copy of tablegeneration_main with three differences: the input file, the two
%   output names, and no awake-normalised arm -- a drift baseline is a position, not
%   a magnitude, so a ratio against it has no meaning. The two files are meant to be
%   folded into one function once this one has a figure on it; see CLAUDE_LOG.md
%
%   The first version of this file built its own bout table by hand with a local
%   quiet-flank baseline. That is superseded: the series is now flattened by
%   detrend in centralize_driftsession and goes through state_1d like a thickness.
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';      % or '00_igkl' / '01_igkltdt'
end
dirs = dirs_central();
dirs.out = fullfile(dirs.secondary_root, param.dataset);

param.nest_names = ["state_summary", "transition"];
param.key_columns = ["MouseID", "Date", "SessionType", "SessionID", "VesselID", ...
    "NumericDepth", "NumericResolution"];

%% Filters : which rows the tables are about
filt.depth_thr    = 70;                                   % um   L1 / L2-3 boundary
filt.min_bout_s   = 30;                                   % s    shorter is not a state
filt.state_names  = ["awake", "drowsy", "nrem", "rem"];   % str  the scored states drawn
filt.require_fwhm = true;                                 % bool only recordings the FWHM pipeline has judged

%% Flatten the centralized state analysis
central_path = fullfile(dirs.central, 'centralized_drift_state.mat');
if ~isfile(central_path)
    error('tablegeneration_drift:noCentral', ...
        'run centralize_driftstate first; %s is missing', central_path);
end
loaded_central = load(central_path);
central = loaded_central.central;
clear loaded_central
fprintf('centralized drift state : %d sessions\n   %s\n', height(central), central_path);

key_table = central(:, param.key_columns);
fwhm_state = load(fullfile(dirs.central, 'centralized_paxfwhm_state.mat')).central;
key_fwhm = util_sessionkey(fwhm_state);
clear fwhm_state
n_session = height(central);
subtable = struct();
for name_idx = 1:numel(param.nest_names)
    nest_name = param.nest_names(name_idx);
    nest_cell = cell(n_session, 1);
    for session_idx = 1:n_session
        session_row = central.data{session_idx};
        if isfield(session_row, nest_name)
            nest_cell{session_idx} = session_row.(nest_name);
        end
    end
    subtable.(nest_name) = explode_nest(nest_cell, key_table);
    fprintf('   %-14s %d rows\n', nest_name, height(subtable.(nest_name)));
end
clear central

%% State summary analysis
analysis_table = subtable.state_summary;
vessel_id = string(analysis_table.VesselID);
is_artery = contains(vessel_id, "PA", 'IgnoreCase', true);
is_layer1 = analysis_table.NumericDepth < filt.depth_thr;
is_longbout = analysis_table.bout_duration > filt.min_bout_s;
is_scored = ismember(analysis_table.state_name, filt.state_names);
has_fwhm = ismember(util_sessionkey(analysis_table), key_fwhm) | ~filt.require_fwhm;
keep_summary = is_artery & is_layer1 & is_longbout & is_scored & has_fwhm;
fprintf('state_summary : %d of %d rows kept\n', sum(keep_summary), numel(keep_summary));
summary_filtered = analysis_table(keep_summary, :);

%% Absolute analysis, px -> um
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
summarystatAnalyzer = tableAnalyzer(summary_filtered, struct());
summarystatAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames,"division");
summarystatAnalyzer = common_summarystat_analysis(summarystatAnalyzer);
result.drift_summary.abs_numeric = summarystatAnalyzer.numeric_tables;

%% Awake subtraction
tmp.awakenumtable = result.drift_summary.abs_numeric.bout_idx_ave(result.drift_summary.abs_numeric.bout_idx_ave.state_name=="awake",:);
tmp.key_cols = ["MouseID","Date","VesselID","DataType"];
tmp.awake_baseline = tmp.awakenumtable(:, [tmp.key_cols,"Q2Q3_mean"]);
tmp.awake_baseline.Properties.VariableNames("Q2Q3_mean") = "awake_boutq2q3mean";

summarystat_awake = outerjoin(result.drift_summary.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
awakesubtract_summarystatAnalyzer = tableAnalyzer(summarystat_awake, struct());
awakesubtract_summarystatAnalyzer.scale_table("awake_boutq2q3mean",data_colnames,numeric_colnames,"subtraction");
awakesubtract_summarystatAnalyzer = common_summarystat_analysis(awakesubtract_summarystatAnalyzer);
result.drift_summary.awakesubtract_numeric = awakesubtract_summarystatAnalyzer.numeric_tables;

%% Transition analysis, sleep recordings only
analysis_table = subtable.transition;
vessel_id = string(analysis_table.VesselID);
is_artery = contains(vessel_id, "PA", 'IgnoreCase', true);
is_layer1 = analysis_table.NumericDepth < filt.depth_thr;
session_type = string(analysis_table.SessionType);
is_sleep = session_type == "sleep";
has_fwhm = ismember(util_sessionkey(analysis_table), key_fwhm) | ~filt.require_fwhm;
keep_transition = is_artery & is_layer1 & is_sleep & has_fwhm;
dropped_type = unique(session_type(~is_sleep));
fprintf('SessionType : sleep %d rows | dropped %d (%s)\n', ...
    nnz(is_sleep), nnz(~is_sleep), strjoin(dropped_type(:).', ','));
fprintf('transition : %d of %d rows kept\n', sum(keep_transition), numel(keep_transition));
transition_filtered = analysis_table(keep_transition, :);

data_colnames = {"data"};
numeric_colnames = {'pre_mean','pre_median','pre_q1','pre_q3', 'pre_var',...
                'post_mean','post_median','post_q1','post_q3', 'post_var'};
transitionAnalyzer = tableAnalyzer(transition_filtered, struct());
transitionAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames,"division");
transitionAnalyzer = common_transition_analysis(transitionAnalyzer);
result.drift_transition.abs_numeric = transitionAnalyzer.numeric_tables;
result.drift_transition.abs_data = transitionAnalyzer.data_tables;

transition_awake = outerjoin(result.drift_transition.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
awakesubtracted_transitionAnalyzer = tableAnalyzer(transition_awake, struct());
awakesubtracted_transitionAnalyzer.scale_table("awake_boutq2q3mean", data_colnames, numeric_colnames,"subtraction");
awakesubtracted_transitionAnalyzer = common_transition_analysis(awakesubtracted_transitionAnalyzer);
result.drift_transition.awakesubtract_numeric = awakesubtracted_transitionAnalyzer.numeric_tables;
result.drift_transition.awakesubtract_data = awakesubtracted_transitionAnalyzer.data_tables;

%% saving mechanism
if ~isfolder(dirs.out)
    mkdir(dirs.out)
end
tmp.savefieldnames = fieldnames(result);
tmp.savefilenames = strcat(string(tmp.savefieldnames),'.mat');
tmp.savepaths = fullfile(dirs.out, string(tmp.savefilenames));

for save_idx = 1:numel(tmp.savefieldnames)
    save_content = result.(tmp.savefieldnames{save_idx});
    save(tmp.savepaths{save_idx}, "save_content", "param", "filt")
    fprintf('wrote %s\n', tmp.savepaths{save_idx});
end


function tableAnalyzer = common_summarystat_analysis(tableAnalyzer)
    tableAnalyzer.meanFrom2("raw_data","Q2Q3_mean",0.25,0.75)
    tableAnalyzer.addPrctilecol("raw_data","prctile_95", 95);
    tableAnalyzer.addPrctilecol("raw_data","prctile_5",5);
    tableAnalyzer.numeric_tables.filtered_table = tableAnalyzer.filtered_table;
    tableAnalyzer.get_numericsummary("bout_idx","filtered_table") % intra session
    tableAnalyzer.get_numericsummary("Date","bout_idx_ave") % inter session
    tableAnalyzer.get_numericsummary("VesselID","Date_ave") % Vessel
    tableAnalyzer.get_numericsummary("MouseID","VesselID_ave") % mouse
end


function tableAnalyzer = common_transition_analysis(tableAnalyzer)
    tableAnalyzer.numeric_tables.filtered_table = tableAnalyzer.filtered_table;
    tableAnalyzer.get_numericsummary("bout_idx","filtered_table") % intra session bout average
    tableAnalyzer.get_numericsummary("Date","bout_idx_ave") % inter session average
    tableAnalyzer.get_numericsummary("VesselID","Date_ave") % within mouse vessel average
    tableAnalyzer.get_numericsummary("MouseID","VesselID_ave") % mouse average
    tableAnalyzer.get_datasummary("bout_idx","filtered_table")
    tableAnalyzer.get_datasummary("Date","bout_idx_ave")
    tableAnalyzer.get_datasummary("VesselID","Date_ave")
    tableAnalyzer.get_datasummary("MouseID","VesselID_ave")
end
