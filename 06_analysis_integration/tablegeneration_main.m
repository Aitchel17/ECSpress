%TABLEGENERATION_MAIN  Centralized state analysis -> the two tables the figures read.
%   Reads centralized_paxfwhm_state.mat, flattens the two nests it works on, applies
%   the analysis filters, and writes state_summary.mat and transition.mat beside the
%   dataset.
%
%   It used to walk the dirtable, open one paxfwhm_state.mat per session folder and
%   rebuild mtable_FWHMsleep.mat on the way. That is centralize_primary and
%   centralize_state now, which happen once and are incremental; this file no longer
%   touches the data tree, the sheet, or the mtable.
%
%   The FILTERS below are the analysis decision -- artery only, layer 1, bouts over
%   30 s, four scored states, sleep recordings for the transitions. They stay in this
%   file rather than moving up into integrate_analysisresult, and they are the reason
%   nothing caches what the flattening produces: change one and every table below is
%   a different table.
%
%   integrate_analysisresult exports the two environment variables and they survive
%   the clear below. Run this file on its own and the defaults apply.
%     caution  the OUTPUT lands in the dataset folder. Pointing this at another
%              dataset writes a different pair of files
clc, clear
addpath('g:\03_program\01_ecspress\functions');   % where util_ecspath lives
util_ecspath;                                     % three roots, minus zz_notinuse

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';      % or '00_igkl' / '01_igkltdt'
end
dirs = util_centraldirs();
dirs.central = fullfile(dirs.secondary_root, 'centralized');
dirs.out = fullfile(dirs.secondary_root, param.dataset);

% The two nests this file works on. centralize_state carries five; band_decomposition,
% powerdensity and peak_trough are not opened here and stay nested where they are.
param.nest_names = ["state_summary", "transition"];

% What a row of the flattened tables carries before DataType and the series' own
% columns. Every one is read below: the first four say which recording, VesselID and
% NumericDepth are filtered on, NumericResolution is the pixel-to-micrometre scale.
param.key_columns = ["MouseID", "Date", "SessionType", "SessionID", "VesselID", ...
    "NumericDepth", "NumericResolution"];

param.depth_thr  = 70;                                    % um  L1 / L2-3 boundary
param.min_bout_s = 30;                                    % s   shorter is not a state
param.state_names = ["awake", "drowsy", "nrem", "rem"];

%% Flatten the centralized state analysis
central_path = fullfile(dirs.central, 'centralized_paxfwhm_state.mat');
if ~isfile(central_path)
    error('tablegeneration:noCentral', ...
        'run centralize_state first; %s is missing', central_path);
end
loaded_central = load(central_path);
central = loaded_central.central;
clear loaded_central
fprintf('centralized state : %d sessions\n   %s\n', height(central), central_path);

key_table = central(:, param.key_columns);
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
is_layer1 = analysis_table.NumericDepth < param.depth_thr;
is_longbout = analysis_table.bout_duration > param.min_bout_s;
is_scored = ismember(analysis_table.state_name, param.state_names);
keep_summary = is_artery & is_layer1 & is_longbout & is_scored;
fprintf('state_summary : %d of %d rows kept\n', sum(keep_summary), numel(keep_summary));
summary_filtered = analysis_table(keep_summary, :);

%% Absolute analysis
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
summarystatAnalyzer = tableAnalyzer(summary_filtered, struct());
summarystatAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames,"division");
summarystatAnalyzer = common_summarystat_analysis(summarystatAnalyzer);
result.state_summary.abs_numeric = summarystatAnalyzer.numeric_tables;


% Awake normalization/subtraction
tmp.awakenumtable = result.state_summary.abs_numeric.bout_idx_ave(result.state_summary.abs_numeric.bout_idx_ave.state_name=="awake",:);

% 0. Extract the denominator baseline
tmp.key_cols = ["MouseID","Date","VesselID","DataType"];
tmp.awake_baseline = tmp.awakenumtable(:, [tmp.key_cols,"Q2Q3_mean"]);
tmp.awake_baseline.Properties.VariableNames("Q2Q3_mean") = "awake_boutq2q3mean";
tmp.awake_baseline.awake_scale = 1./tmp.awake_baseline.awake_boutq2q3mean;

summarystat_awakenorm = outerjoin(result.state_summary.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
% 1. get awake normalized table
awakescaled_summarystatAnalyzer = tableAnalyzer(summarystat_awakenorm, struct());
awakescaled_summarystatAnalyzer.scale_table("awake_scale",data_colnames,numeric_colnames,"division");
awakescaled_summarystatAnalyzer = common_summarystat_analysis(awakescaled_summarystatAnalyzer);
result.state_summary.awakenorm_numeric = awakescaled_summarystatAnalyzer.numeric_tables;

% 2. get awake subtracted table
awakesubtract_summarystatAnalyzer = tableAnalyzer(summarystat_awakenorm, struct());
awakesubtract_summarystatAnalyzer.scale_table("awake_boutq2q3mean",data_colnames,numeric_colnames,"subtraction");
awakesubtract_summarystatAnalyzer = common_summarystat_analysis(awakesubtract_summarystatAnalyzer);
result.state_summary.awakesubtract_numeric = awakesubtract_summarystatAnalyzer.numeric_tables;



%% transition analysis
analysis_table = subtable.transition;
vessel_id = string(analysis_table.VesselID);
is_artery = contains(vessel_id, "PA", 'IgnoreCase', true);
is_layer1 = analysis_table.NumericDepth < param.depth_thr;
% Sleep only.
%   err  without this a whiskerb session rides along : hql090 251020 sets
%        transition_window to 20 s where every sleep session uses 50, so its traces
%        are a different WINDOW and not just a different sampling of one
% SessionType is a column of this table now. It used to be reached by indexing the
% dirtable with SessionIndex, explode_nest's loop counter, which needed a check that
% the counter still meant the row it stood for. see CLAUDE_LOG.md
session_type = string(analysis_table.SessionType);
is_sleep = session_type == "sleep";
keep_transition = is_artery & is_layer1 & is_sleep;
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
result.transition.abs_numeric = transitionAnalyzer.numeric_tables;
result.transition.abs_data = transitionAnalyzer.data_tables;

% Normalize transition data by awake Q2Q3_mean
transition_awakenorm = outerjoin(result.transition.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
% Apply normalization using awake_scale
awakescaled_transitionAnalyzer = tableAnalyzer(transition_awakenorm, struct());
awakescaled_transitionAnalyzer.scale_table("awake_scale", data_colnames, numeric_colnames,"division");
% Re-calculate numeric and data summaries
awakescaled_transitionAnalyzer = common_transition_analysis(awakescaled_transitionAnalyzer);
result.transition.awakenorm_numeric = awakescaled_transitionAnalyzer.numeric_tables;
result.transition.awakenorm_data = awakescaled_transitionAnalyzer.data_tables;

awakesubtracted_transitionAnalyzer = tableAnalyzer(transition_awakenorm, struct());
awakesubtracted_transitionAnalyzer.scale_table("awake_boutq2q3mean", data_colnames, numeric_colnames,"subtraction");
awakesubtracted_transitionAnalyzer = common_transition_analysis(awakesubtracted_transitionAnalyzer);
result.transition.awakesubtract_numeric = awakesubtracted_transitionAnalyzer.numeric_tables;
result.transition.awakesubtract_data = awakesubtracted_transitionAnalyzer.data_tables;


%% saving mechanism
if ~isfolder(dirs.out)
    mkdir(dirs.out)
end
tmp.savefieldnames = fieldnames(result);
tmp.savefilenames = strcat(string(tmp.savefieldnames),'.mat');
tmp.savepaths = fullfile(dirs.out, string(tmp.savefilenames));

for save_idx = 1:numel(tmp.savefieldnames)
    save_content = result.(tmp.savefieldnames{save_idx});
    save(tmp.savepaths{save_idx}, "save_content")
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
