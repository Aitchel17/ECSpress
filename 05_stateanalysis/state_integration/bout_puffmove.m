function [time_table, param] = bout_puffmove(air_on, air_fs, speed, taxis, fs)
%BOUT_PUFFMOVE  The puffs crossed with the movement, as bout tables.
%   Caller: centralize_boutsession, review_session
%   The movement block reads the ball; split_puffs and zhang_loco cross it with the puffs and
%   never see the sensor, so another rig (a cylinder) replaces only that block.
%
% IN   air_on  T x 1 bool       the air-puff channel, on (> 0) per sample
%      air_fs  1 x 1 float Hz   its sample rate
%      speed   S float cm/s     the ball speed, signed; analysis_analog's ball.resampled_ball
%      taxis   S float s        its time axis, ball.rs_taxis
%      fs      1 x 1 float Hz   its rate, ball.ds_fps
% OUT  time_table      1x1 struct   table name -> N x 2 [start end] bouts in seconds :
%        dur<X>                      every puff of that duration
%        puff_go_dur<X>              a puff to a mouse not moving in [onset - puff_quiet_s,
%                                    onset - puff_quietend_s] that starts moving during the
%                                    puff. Movement is Zhang's |acceleration|, the stretches
%                                    loco_* are built from
%        puff_nogo_dur<X>            the same, still through the puff
%        puff_moving_dur<X>          a puff to a mouse already moving before it
%        loco_puff / loco_nopuff     Zhang's locomotion events, by whether the onset falls
%                                    during a puff
%        loco_rest_puff / _nopuff    Zhang's rests, cut at the puffs; _nopuff is true rest
%      param           1x1 struct   the settings below, and per duration dur<X> each puff's
%                                   Duration, FrequencyHz, DutyCycle and rounded_duration

    %% Criteria
    % Puff
    param.puff_gap_s = 10;          % pulse onsets further apart than this are separate puffs (s), 02_othersignal get_airtable
    param.puff_quiet_s = 3;        % puff criterion: still for this long before onset, else puff_moving (s), Zhang 2026
    param.puff_quietend_s = 0.1;    % that window ends this long before onset; filtfilt smears an evoked start earlier (s)

    % Locomotion preprocessing
    param.loco_lowpass_hz = 10;     % low-pass on the speed before the derivative (Hz), Zhang 2026
    param.loco_accel_thresh = 3;    % |acceleration| at or above it is moving (cm/s^2), Zhang 2026

    % Locomotion time
    param.loco_quiet_s = 3;         % minimum rest before a locomotion event (s), Zhang 2026
    param.loco_min_s = 5;           % minimum locomotion event duration (s), Zhang 2026
    param.loco_join_s = 1;          % moving stretches closer than this are joined (s), Zhang 2022
    % Rest time
    param.loco_rest_pad_s = 4;      % rest starts this long after the last movement (s), Zhang 2026
    param.loco_rest_min_s = 10;     % minimum rest duration, after the pad (s), Zhang 2026

    %% the puffs, off the air channel : one table per rounded duration
    atable = group_pulses(air_on, air_fs, param.puff_gap_s);
    atable.('rounded_duration')= round(atable.Duration);
    unique_dur = unique(atable.rounded_duration);
    time_table = struct();
    for idx = 1:numel(unique_dur)
        stim_dur = unique_dur(idx); % duration double 3,5,...
        table_loc = atable.rounded_duration == unique_dur(idx);
        time_table.(strcat("dur",string(stim_dur))) = atable{table_loc,1:2};
        param.(strcat("dur",string(stim_dur))) = atable(table_loc,3:end);
    end

    %% the movement, off the ball : Zhang's moving stretches, on the acceleration ; both directions count
    is_moving = gate_accel(speed, fs, param.loco_lowpass_hz, param.loco_accel_thresh);
    moving_stretch = logic2timetable(taxis, is_moving, 0);
    moving_table = bout_joingaps(moving_stretch, param.loco_join_s);

    %% puff-anchored : each duration's puffs by the movement before and during them
    all_puff = [];   % P x 2 float s  every puff, in unique_dur order
    for idx = 1:numel(unique_dur)
        dur_field = strcat("dur", string(unique_dur(idx)));
        puff_table = time_table.(dur_field);
        [puff_go, puff_nogo, puff_moving] = split_puffs(puff_table, moving_table, ...
            param.puff_quiet_s, param.puff_quietend_s);
        time_table.(strcat("puff_go_", dur_field)) = puff_go;
        time_table.(strcat("puff_nogo_", dur_field)) = puff_nogo;
        time_table.(strcat("puff_moving_", dur_field)) = puff_moving;
        all_puff = [all_puff; puff_table]; %#ok<AGROW>
    end

    %% Zhang's locomotion events and rests, split by puff
    record_span = [taxis(1), taxis(end)];
    [time_table.loco_puff, time_table.loco_nopuff, time_table.loco_rest_nopuff, time_table.loco_rest_puff] = ...
        zhang_loco(moving_table, record_span, all_puff, ...
        param.loco_quiet_s, param.loco_min_s, param.loco_rest_pad_s, param.loco_rest_min_s);
end

function puff_table = group_pulses(air_on, air_fs, gap_s)
%GROUP_PULSES  The puffs on the air channel: pulses whose onsets lie within gap_s are one puff.
%   02_othersignal's analysis_analog.get_airtable, on the sample clock (k - 1) / air_fs.
% IN   air_on      T x 1 bool         the air channel, on per sample
%      air_fs      1 x 1 float Hz
%      gap_s       1 x 1 float s      param.puff_gap_s
% OUT  puff_table  N x 5 table        StartTime, EndTime, Duration (s), FrequencyHz, DutyCycle;
%                                     the last two the mean over the pulses that have a next
%                                     one in the same puff, NaN for a single pulse
    edge = diff(air_on(:));
    rise_idx = find(edge == 1);    % the last off sample before each pulse
    fall_idx = find(edge == -1);   % the last on sample of each pulse
    rise_s = (rise_idx - 1) / air_fs;
    fall_s = fall_idx / air_fs;    % the first off sample after it
    onset_gap = diff(rise_s);
    is_first = [true; onset_gap > gap_s];
    is_last = [onset_gap > gap_s; true];
    start_s = rise_s(is_first);
    end_s = fall_s(is_last);
    % pulse k against pulse k + 1, kept where both are in one puff
    pulse_duty = (fall_idx(1:end-1) - rise_idx(1:end-1)) ./ diff(rise_idx);
    has_next = ~is_last(1:end-1);
    puff_id = cumsum(is_first);
    pulse_puff = puff_id(has_next);
    n_puff = numel(start_s);
    mean_gap = accumarray(pulse_puff, onset_gap(has_next), [n_puff 1], @mean, NaN);
    mean_duty = accumarray(pulse_puff, pulse_duty(has_next), [n_puff 1], @mean, NaN);
    puff_table = table(start_s, end_s, end_s - start_s, 1 ./ mean_gap, mean_duty, ...
        'VariableNames', {'StartTime', 'EndTime', 'Duration', 'FrequencyHz', 'DutyCycle'});
end

function is_moving = gate_accel(speed, fs, lowpass_hz, accel_thresh)
%GATE_ACCEL  Zhang's moving samples: the low-passed speed's |acceleration| at or over accel_thresh.
% IN   speed         T float cm/s         the ball speed, signed
%      fs            1 x 1 float Hz       its rate
%      lowpass_hz    1 x 1 float Hz       on the speed, before the derivative
%      accel_thresh  1 x 1 float cm/s^2
% OUT  is_moving     1 x T bool
    [lp_zero, lp_pole, lp_gain] = butter(5, lowpass_hz / (fs / 2), 'low');
    [lp_sos, lp_g] = zp2sos(lp_zero, lp_pole, lp_gain);
    speed_lp = filtfilt(lp_sos, lp_g, speed(:)');
    is_moving = abs([0, diff(speed_lp)] * fs) >= accel_thresh;
end

function [puff_go, puff_nogo, puff_moving] = split_puffs(puff_table, moving_table, quiet_s, quietend_s)
%SPLIT_PUFFS  One duration's puffs, by whether a moving stretch falls before them or during them.
% IN   puff_table    N x 2 float s   one duration's puffs, [onset end]
%      moving_table  L x 2 float s   the moving stretches
%      quiet_s       1 x 1 float s   param.puff_quiet_s
%      quietend_s    1 x 1 float s   param.puff_quietend_s
% OUT  puff_go       K x 2 float s   to a mouse not moving before, moving during the puff
%      puff_nogo     J x 2 float s   to a mouse not moving before, with none
%      puff_moving   I x 2 float s   to a mouse already moving before
    quiet_end = puff_table(:, 1) - quietend_s;   % N x 1 float s
    quiet_table = [quiet_end, quiet_end];   % N x 2 float s  the window's end as a point
    quiet_span = quiet_s - quietend_s;
    is_moving_before = bout_hasoverlap(quiet_table, moving_table, quiet_span, 0);
    is_go = bout_hasoverlap(puff_table, moving_table, 0, 0);
    keep_go = ~is_moving_before & is_go;
    keep_nogo = ~is_moving_before & ~is_go;
    puff_go = puff_table(keep_go, :);
    puff_nogo = puff_table(keep_nogo, :);
    puff_moving = puff_table(is_moving_before, :);
end

function [loco_puff, loco_nopuff, loco_rest_nopuff, loco_rest_puff] = zhang_loco(moving_table, record_span, all_puff, ...
    quiet_s, min_s, rest_pad_s, rest_min_s)
%ZHANG_LOCO  Zhang's locomotion events and rests on the moving stretches, split by puff.
%   The events and rests PNAS 2026 aligns to.
% IN   moving_table  L x 2 float s   the moving stretches
%      record_span   1 x 2 float s   the record's first and last time
%      all_puff      P x 2 float s   every puff, [onset end]
%      quiet_s       1 x 1 float s   param.loco_quiet_s
%      min_s         1 x 1 float s   param.loco_min_s
%      rest_pad_s    1 x 1 float s   param.loco_rest_pad_s
%      rest_min_s    1 x 1 float s   param.loco_rest_min_s
% OUT  loco_puff, loco_nopuff             events, by whether the onset falls during a puff
%      loco_rest_nopuff, loco_rest_puff   rests outside every puff (Winder's true rest), and the piece a puff landed on
    quiet_before = moving_table(:, 1) - [record_span(1); moving_table(1:end-1, 2)];
    long_enough = moving_table(:, 2) - moving_table(:, 1) >= min_s;
    loco_event = moving_table(quiet_before >= quiet_s & long_enough, :);
    % an event whose onset falls during a puff is the puff's; the rest are spontaneous
    onset_in_puff = false(size(loco_event, 1), 1);
    for event_idx = 1:size(loco_event, 1)
        onset_in_puff(event_idx) = any(all_puff(:, 1) <= loco_event(event_idx, 1) & loco_event(event_idx, 1) <= all_puff(:, 2));
    end
    loco_puff = loco_event(onset_in_puff, :);
    loco_nopuff = loco_event(~onset_in_puff, :);
    % a rest : from the pad after one moving stretch to the start of the next, the record's ends included
    rest_start = [record_span(1); moving_table(:, 2) + rest_pad_s];
    rest_end = [moving_table(:, 1); record_span(2)];
    rest_table_loco = [rest_start, rest_end];
    rest_table_loco = rest_table_loco(rest_end - rest_start > rest_min_s, :);
    % a rest the puff landed on is the puff's; what lies outside every puff is Winder's true rest
    [loco_rest_in_window, loco_rest_nopuff] = ...
        bout_cutwindows(rest_table_loco, all_puff, rest_min_s);
    loco_rest_puff = loco_rest_in_window(bout_beginsat(loco_rest_in_window, all_puff(:, 1)), :);
end
