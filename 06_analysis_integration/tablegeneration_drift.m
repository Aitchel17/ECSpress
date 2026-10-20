function info = tablegeneration_drift(central_dir, out_dir, vessel_prefix, depth_thr, min_bout_s, state_names, require_fwhm)
%TABLEGENERATION_DRIFT  Centralized drift state analysis -> the two tables the drift figures read.
%   Reads its two long tables out of centralized_drift_state.mat, applies the caller's
%   row filters -- tablegeneration_paxfwhm's plus require_fwhm -- and writes
%   drift_summary.mat and drift_transition.mat into out_dir, replacing the copies there,
%   with filt saved beside each table. e_tablegeneration sets every filter.
%   Caller: e_tablegeneration
%
%   A copy of tablegeneration_paxfwhm with three differences: the input file, the two
%   output names, and no awake-normalised arm -- a drift baseline is a position, not
%   a magnitude, so a ratio against it has no meaning. Folding the two into one
%   function waits until the radon and PIV tables arrive.
%
% IN   central_dir    1 x 1 str     where the centralized_*.mat files live
%      out_dir        1 x 1 str     the dataset folder the two tables are written to
%      vessel_prefix  1 x 1 str     VesselID prefix kept ("PA" = arteries)
%      depth_thr      1 x 1 um      kept shallower than this
%      min_bout_s     1 x 1 s       bouts kept longer than this
%      state_names    1 x S str     the scored states kept
%      require_fwhm   1 x 1 bool    only recordings centralized_paxfwhm_state holds
% OUT  info           1 x 1 struct  n_summary, n_transition, out_path

arguments
    central_dir (1,1) string
    out_dir (1,1) string
    vessel_prefix (1,1) string
    depth_thr (1,1) double
    min_bout_s (1,1) double
    state_names (1,:) string
    require_fwhm (1,1) logical
end

% 1. check the drift state file exists, and read the sessions the FWHM pipeline judged.
%    Every row of its two long tables carries bout_drift's carried_columns, session_key
%    and DataType
driftstate_path = fullfile(central_dir, 'centralized_drift_state.mat');
if ~isfile(driftstate_path)
    error('tablegeneration_drift:noCentral', ...
        'run bout_drift first; %s is missing', driftstate_path);
end
fprintf('centralized drift state : %s\n', driftstate_path);

paxfwhmstate_path = fullfile(central_dir, 'centralized_paxfwhm_state.mat');
fwhm_state = load(paxfwhmstate_path, 'central').central;
key_fwhm = util_sessionkey(fwhm_state);
clear fwhm_state

% 2. filter state_summary
state_summary = load(driftstate_path, 'state_summary').state_summary;
vessel_id = string(state_summary.VesselID);
is_artery = startsWith(vessel_id, vessel_prefix, 'IgnoreCase', true);
is_layer1 = state_summary.NumericDepth < depth_thr;
is_longbout = state_summary.bout_duration > min_bout_s;
is_scored = ismember(state_summary.state_name, state_names);
has_fwhm = ismember(state_summary.session_key, key_fwhm) | ~require_fwhm;
keep_summary = is_artery & is_layer1 & is_longbout & is_scored & has_fwhm;
fprintf('state_summary : %d of %d rows kept\n', sum(keep_summary), numel(keep_summary));
summary_filtered = state_summary(keep_summary, :);

% 3. scale state_summary px -> um and roll it up
data_colnames = {"raw_data"};
numeric_colnames = {'raw_mean','raw_median','raw_q1','raw_q3', 'raw_var'};
summarystatAnalyzer = tableAnalyzer(summary_filtered, struct());
summarystatAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames,"division");
summarystatAnalyzer = rollup_summary(summarystatAnalyzer);
result.drift_summary.abs_numeric = summarystatAnalyzer.numeric_tables;

% 4. awake-subtracted state_summary, against the awake baseline per vessel and series
tmp.key_cols = ["MouseID","Date","VesselID","DataType"];
tmp.awake_baseline = extract_awakebaseline(result.drift_summary.abs_numeric.bout_idx_ave, tmp.key_cols);

summarystat_awake = outerjoin(result.drift_summary.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
awakesubtract_summarystatAnalyzer = tableAnalyzer(summarystat_awake, struct());
awakesubtract_summarystatAnalyzer.scale_table("awake_boutq2q3mean",data_colnames,numeric_colnames,"subtraction");
awakesubtract_summarystatAnalyzer = rollup_summary(awakesubtract_summarystatAnalyzer);
result.drift_summary.awakesubtract_numeric = awakesubtract_summarystatAnalyzer.numeric_tables;

% 5. filter transition, sleep recordings only
transition = load(driftstate_path, 'transition').transition;
vessel_id = string(transition.VesselID);
is_artery = startsWith(vessel_id, vessel_prefix, 'IgnoreCase', true);
is_layer1 = transition.NumericDepth < depth_thr;
session_type = string(transition.SessionType);
is_sleep = session_type == "sleep";
has_fwhm = ismember(transition.session_key, key_fwhm) | ~require_fwhm;
keep_transition = is_artery & is_layer1 & is_sleep & has_fwhm;
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
result.drift_transition.abs_numeric = transitionAnalyzer.numeric_tables;
result.drift_transition.abs_data = transitionAnalyzer.data_tables;

% 7. awake-subtracted transition, against the step 4 baseline
transition_awake = outerjoin(result.drift_transition.abs_numeric.filtered_table, tmp.awake_baseline, ...
    'MergeKeys', true, ...
    'LeftKeys', tmp.key_cols, ...
    'RightKeys', tmp.key_cols);
awakesubtracted_transitionAnalyzer = tableAnalyzer(transition_awake, struct());
awakesubtracted_transitionAnalyzer.scale_table("awake_boutq2q3mean", data_colnames, numeric_colnames,"subtraction");
awakesubtracted_transitionAnalyzer = rollup_transition(awakesubtracted_transitionAnalyzer);
result.drift_transition.awakesubtract_numeric = awakesubtracted_transitionAnalyzer.numeric_tables;
result.drift_transition.awakesubtract_data = awakesubtracted_transitionAnalyzer.data_tables;

% 8. save one file per table, the filters that made it beside it
filt = struct('vessel_prefix', vessel_prefix, 'depth_thr', depth_thr, 'min_bout_s', min_bout_s, ...
    'state_names', state_names, 'require_fwhm', require_fwhm);
if ~isfolder(out_dir)
    mkdir(out_dir)
end
tmp.savefieldnames = fieldnames(result);
tmp.savefilenames = strcat(string(tmp.savefieldnames),'.mat');
tmp.savepaths = fullfile(out_dir, string(tmp.savefilenames));

for save_idx = 1:numel(tmp.savefieldnames)
    save_content = result.(tmp.savefieldnames{save_idx});
    save(tmp.savepaths{save_idx}, "save_content", "filt")
    fprintf('wrote %s\n', tmp.savepaths{save_idx});
end
info = struct('n_summary', sum(keep_summary), 'n_transition', sum(keep_transition), ...
    'out_path', tmp.savepaths);
end
