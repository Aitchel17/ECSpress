function [row, failure] = centralize_driftsession(k, disp_row, disp_col, drift_fps, ...
    score, primary_analog, stamp, old)
%CENTRALIZE_DRIFTSESSION  One session's row of centralized_drift_state.
%   Caller: centralize_driftstate
%   The two drift components go through state_1d the way the nine FWHM series do
%   in centralize_statesession, so the tables downstream read both alike. Each is
%   first flattened by one straight line over the whole session, the way
%   analog_ecogspectrum flattens the ECoG.
%   IN   k               1x1 double        which row of the drift table this is
%        disp_row        N x 1 double      tissue displacement along rows, px
%        disp_col        N x 1 double      along columns, px
%        drift_fps       1x1 double        Hz
%        score           1x1 struct | []   the sleep score, if it has one
%        primary_analog  1x1 struct | []   the whisker scoring instead
%        stamp           1x1 double        what this row is rebuilt against
%        old             0|1 row table     what the last run wrote for this session
%   OUT  row             0x0 | 1x1 struct  drift_i, data, input_stamp, reused
%        failure         1x1 str           "" unless the analysis threw
    row = struct('drift_i', {}, 'data', {}, 'input_stamp', {}, 'reused', {});
    failure = "";

    unmoved = ~isempty(old) && old.input_stamp == stamp;
    if ~unmoved
        try
            t_drift = (0:numel(disp_row) - 1) / drift_fps;
            state_integrate = state_integration(score, primary_analog);
            state_integrate.trim_to_duration(t_drift(end));
            drift_state = state_1d(state_integrate);
            drift_state.get_state_indices(t_drift, drift_fps);

            row_flat = detrend(disp_row)';    % 1 x N, the orientation the FWHM series carry
            col_flat = detrend(disp_col)';
            drift_state.get_summary("drift_row", row_flat);
            drift_state.get_transitionsummary("drift_row", row_flat);
            drift_state.get_summary("drift_col", col_flat);
            drift_state.get_transitionsummary("drift_col", col_flat);

            row(1).drift_i = k;
            row(1).data = strip_state(drift_state);
            row(1).input_stamp = stamp;
            row(1).reused = false;
            return
        catch err
            failure = string(err.message);
        end
    end
    if isempty(old)
        return
    end
    row(1).drift_i = k;
    row(1).data = old.data{1};
    row(1).input_stamp = stamp;
    row(1).reused = true;
end

function out = strip_state(drift_state)
    % the two result structs the tables nest on, the axes, and the pieces of the
    % state_integration anything outside the class reads
    out = struct();
    for f = ["state_summary", "transition", "t_axis", "state_idx", "param"]
        out.(f) = drift_state.(f);
    end
    out.sleep_obj = struct();
    if ~isempty(drift_state.sleep_obj)
        out.sleep_obj.param = drift_state.sleep_obj.param;
        out.sleep_obj.time_table = drift_state.sleep_obj.time_table;
        out.sleep_obj.info = drift_state.sleep_obj.info;
    end
end
