function info = tablegeneration_paxfwhm(central_dir, out_dir, vessel_prefix, depth_thr, min_bout_s, state_names)
%TABLEGENERATION_PAXFWHM  Centralized FWHM state analysis -> the two tables the figures read.
%   Reads the two long tables it works on out of centralized_paxfwhm_state.mat, applies
%   the caller's row filters, and writes state_summary.mat and transition.mat into
%   out_dir, replacing the copies there. Nothing is read from the data tree, the sheet
%   or the mtable.
%   Caller: e_tablegeneration
%
%   The row filters are the analysis decision and e_tablegeneration sets them all.
%   Change one and every table below is a different table, which is why nothing caches
%   them. Transitions are sleep recordings only, see step 5.
%
% IN   central_dir    1 x 1 str     where the centralized_*.mat files live
%      out_dir        1 x 1 str     the dataset folder the two tables are written to
%      vessel_prefix  1 x 1 str     VesselID prefix kept ("PA" = arteries)
%      depth_thr      1 x 1 um      kept shallower than this
%      min_bout_s     1 x 1 s       bouts kept longer than this
%      state_names    1 x S str     the scored states kept
% OUT  info           1 x 1 struct  n_summary, n_transition, out_path

arguments
    central_dir (1,1) string
    out_dir (1,1) string
    vessel_prefix (1,1) string
    depth_thr (1,1) double
    min_bout_s (1,1) double
    state_names (1,:) string
end

% 1. check the state file exists. Of its five long tables this reads two; every row
%    carries bout_paxfwhm's carried_columns, session_key and DataType
paxfwhmstate_path = fullfile(central_dir, 'centralized_paxfwhm_state.mat');
if ~isfile(paxfwhmstate_path)
    error('tablegeneration_paxfwhm:noCentral', ...
        'run bout_paxfwhm first; %s is missing', paxfwhmstate_path);
end
fprintf('centralized state : %s\n', paxfwhmstate_path);

% 2. filter state_summary
state_summary = load(paxfwhmstate_path, 'state_summary').state_summary;
vessel_id = string(state_summary.VesselID);
is_artery = startsWith(vessel_id, vessel_prefix, 'IgnoreCase', true);
is_layer1 = state_summary.NumericDepth < depth_thr;
is_longbout = state_summary.bout_duration > min_bout_s;
is_scored = ismember(state_summary.state_name, state_names);
keep_summary = is_artery & is_layer1 & is_longbout & is_scored;
fprintf('state_summary : %d of %d rows kept\n', sum(keep_summary), numel(keep_summary));
summary_filtered = state_summary(keep_summary, :);

% 3. scale state_summary px -> um and roll it up
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
summarystatAnalyzer = tableAnalyzer(summary_filtered, struct());
summarystatAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames,"division");
summarystatAnalyzer = rollup_summary(summarystatAnalyzer);
result.state_summary.abs_numeric = summarystatAnalyzer.numeric_tables;

% 4. awake baseline per vessel and series
tmp.key_cols = ["MouseID","Date","VesselID","DataType"];
tmp.awake_baseline = extract_awakebaseline(result.state_summary.abs_numeric.bout_idx_ave, tmp.key_cols);
tmp.awake_baseline.awake_scale = 1 ./ tmp.awake_baseline.awake_boutq2q3mean;

summarystat_awakenorm = outerjoin(result.state_summary.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
% 4.1 awake-normalized state_summary
awakescaled_summarystatAnalyzer = tableAnalyzer(summarystat_awakenorm, struct());
awakescaled_summarystatAnalyzer.scale_table("awake_scale",data_colnames,numeric_colnames,"division");
awakescaled_summarystatAnalyzer = rollup_summary(awakescaled_summarystatAnalyzer);
result.state_summary.awakenorm_numeric = awakescaled_summarystatAnalyzer.numeric_tables;

% 4.2 awake-subtracted state_summary
awakesubtract_summarystatAnalyzer = tableAnalyzer(summarystat_awakenorm, struct());
awakesubtract_summarystatAnalyzer.scale_table("awake_boutq2q3mean",data_colnames,numeric_colnames,"subtraction");
awakesubtract_summarystatAnalyzer = rollup_summary(awakesubtract_summarystatAnalyzer);
result.state_summary.awakesubtract_numeric = awakesubtract_summarystatAnalyzer.numeric_tables;

% 5. filter transition, sleep recordings only
% whiskerb sessions use a 20 s transition window, sleep sessions 50 s; see CLAUDE_LOG.md
transition = load(paxfwhmstate_path, 'transition').transition;
vessel_id = string(transition.VesselID);
is_artery = startsWith(vessel_id, vessel_prefix, 'IgnoreCase', true);
is_layer1 = transition.NumericDepth < depth_thr;
session_type = string(transition.SessionType);
is_sleep = session_type == "sleep";
keep_transition = is_artery & is_layer1 & is_sleep;
dropped_type = unique(session_type(~is_sleep));
fprintf('SessionType : sleep %d rows | dropped %d (%s)\n', ...
    nnz(is_sleep), nnz(~is_sleep), strjoin(dropped_type(:).', ','));
fprintf('transition : %d of %d rows kept\n', sum(keep_transition), numel(keep_transition));
transition_filtered = transition(keep_transition, :);

% 6. scale transition px -> um and roll it up
data_colnames = {"data"};
numeric_colnames = {'pre_mean','pre_median','pre_q1','pre_q3', 'pre_var',...
                'post_mean','post_median','post_q1','post_q3', 'post_var'};
transitionAnalyzer = tableAnalyzer(transition_filtered, struct());
transitionAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames,"division");
transitionAnalyzer = rollup_transition(transitionAnalyzer);
result.transition.abs_numeric = transitionAnalyzer.numeric_tables;
result.transition.abs_data = transitionAnalyzer.data_tables;

% 7. awake-normalized and awake-subtracted transition, against the step 4 baseline
transition_awakenorm = outerjoin(result.transition.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
% 7.1 awake-normalized transition
awakescaled_transitionAnalyzer = tableAnalyzer(transition_awakenorm, struct());
awakescaled_transitionAnalyzer.scale_table("awake_scale", data_colnames, numeric_colnames,"division");
awakescaled_transitionAnalyzer = rollup_transition(awakescaled_transitionAnalyzer);
result.transition.awakenorm_numeric = awakescaled_transitionAnalyzer.numeric_tables;
result.transition.awakenorm_data = awakescaled_transitionAnalyzer.data_tables;

% 7.2 awake-subtracted transition
awakesubtracted_transitionAnalyzer = tableAnalyzer(transition_awakenorm, struct());
awakesubtracted_transitionAnalyzer.scale_table("awake_boutq2q3mean", data_colnames, numeric_colnames,"subtraction");
awakesubtracted_transitionAnalyzer = rollup_transition(awakesubtracted_transitionAnalyzer);
result.transition.awakesubtract_numeric = awakesubtracted_transitionAnalyzer.numeric_tables;
result.transition.awakesubtract_data = awakesubtracted_transitionAnalyzer.data_tables;

% 8. save one file per table
if ~isfolder(out_dir)
    mkdir(out_dir)
end
tmp.savefieldnames = fieldnames(result);
tmp.savefilenames = strcat(string(tmp.savefieldnames),'.mat');
tmp.savepaths = fullfile(out_dir, string(tmp.savefilenames));

for save_idx = 1:numel(tmp.savefieldnames)
    save_content = result.(tmp.savefieldnames{save_idx});
    save(tmp.savepaths{save_idx}, "save_content")
    fprintf('wrote %s\n', tmp.savepaths{save_idx});
end
info = struct('n_summary', sum(keep_summary), 'n_transition', sum(keep_transition), ...
    'out_path', tmp.savepaths);
end
