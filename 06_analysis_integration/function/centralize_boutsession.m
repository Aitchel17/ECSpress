function [row, failure] = centralize_boutsession(k, score, primary_analog, stamp, old)
%CENTRALIZE_BOUTSESSION  One session's row of centralized_bout.
%   Caller: centralize_bout
%   IN   k               1x1 double        which row of session_table this is
%        score           1x1 struct | []   the sleep score, if it has one
%        primary_analog  1x1 struct | []   the whisker scoring instead
%        stamp           1x1 double        what this row is rebuilt against
%        old             0|1 row table     what the last run wrote for this session
%   OUT  row             0x0 | 1x1 struct  the six fields are declared below
%        failure         1x1 str           "" unless bout_sleep or bout_puffmove threw
    row = struct('session_i', {}, 'scoring', {}, 'data', {}, 'input_stamp', {}, ...
        'computed', {}, 'reused', {});
    failure = "";

    unmoved = ~isempty(old) && old.input_stamp == stamp;
    if ~unmoved
        try
            if ~isempty(score)
                scoring = "sleep";
                [time_table, param] = bout_sleep(score.behavState, score.statecodes, score.binwidth_sec);
            else
                scoring = "analog";
                [time_table, param] = bout_puffmove(primary_analog.air.on, primary_analog.air.fs, ...
                    primary_analog.ball.resampled_ball, primary_analog.ball.rs_taxis, primary_analog.ball.ds_fps);
            end
            row(1).session_i = k;
            row(1).scoring = scoring;
            row(1).data = struct('time_table', time_table, 'param', param);
            row(1).input_stamp = stamp;
            row(1).computed = convertTo(datetime('now'), 'datenum');   % the clock source_modified is on
            row(1).reused = false;
            return
        catch err
            % a scoring that throws is no reason to drop what the last run
            % computed, so the old row stands and the failure is named
            failure = string(err.message);
        end
    end
    if isempty(old)
        return
    end
    row(1).session_i = k;
    row(1).scoring = string(old.scoring);
    row(1).data = old.data{1};
    row(1).input_stamp = stamp;
    row(1).computed = old.computed;
    row(1).reused = true;
end
