function [time_table, param] = bout_sleep(behav_state, state_codes, binwidth_s)
%BOUT_SLEEP  Generate timetable of each bout including individual sleep state, rough category, and transition
%   Caller: centralize_bout, analysis_main_piv, piv_integration_testbed, stateintegration_test
%
% IN   behav_state  B x 1 double    state code per scored bin, sleep_score.mat's behavState
%      state_codes  1x1 struct      state name -> code, sleep_score.mat's statecodes
%      binwidth_s   1 x 1 float s   bin width, sleep_score.mat's binwidth_sec
% OUT  time_table  1x1 struct   table name -> N x 2 [start end] bouts in seconds
%      param       1x1 struct   transition_window (s)

%% 1. parameter
    param.transition_window = 25;   %sec

    time_table = struct();

%% 2. Time table generation
    % rough chunk analysis
    time_table.roughawake = code_bouts(behav_state, [state_codes.NotSleep, state_codes.Drowsy], binwidth_s);
    time_table.roughnrem = code_bouts(behav_state, [state_codes.NREM, state_codes.uarousal], binwidth_s);
    time_table.rem = code_bouts(behav_state, state_codes.REM, binwidth_s);

    % individual state
    time_table.awake = code_bouts(behav_state, state_codes.NotSleep, binwidth_s);
    time_table.drowsy = code_bouts(behav_state, state_codes.Drowsy, binwidth_s);
    time_table.uarousal = code_bouts(behav_state, state_codes.uarousal, binwidth_s);
    time_table.nrem = code_bouts(behav_state, state_codes.NREM, binwidth_s);

%% 3. transition
    time_table.nr_trans = get_transition(param.transition_window, time_table.roughnrem, time_table.rem);
    time_table.na_trans = get_transition(param.transition_window, time_table.roughnrem, time_table.roughawake);
    time_table.ra_trans = get_transition(param.transition_window, time_table.rem, time_table.roughawake);
    time_table.an_trans = get_transition(param.transition_window, time_table.roughawake, time_table.roughnrem);

end

function bouts = code_bouts(behav_state, codes, binwidth_s)
%CODE_BOUTS  The runs of bins whose code is one of codes, as N x 2 [start end] seconds.
%   The scoring GUI's grid : bin k spans [(k-1) binwidth_s, k binwidth_s].
    is_in = ismember(behav_state(:), codes);
    edge = diff([0; is_in; 0]);
    first_bin = find(edge == 1);
    last_bin = find(edge == -1) - 1;
    bouts = [(first_bin - 1) * binwidth_s, last_bin * binwidth_s];
end
