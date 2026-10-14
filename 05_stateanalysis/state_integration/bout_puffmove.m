function [time_table, param] = bout_puffmove(air_on, air_fs, speed, taxis, fs)
%BOUT_PUFFMOVE  The puffs crossed with the movement, as bout tables.
%   Caller: centralize_bout, review_session
%   The movement block reads the ball; split_puffs and zhang_loco cross it with the puffs and
%   never see the sensor, so another rig (a cylinder) replaces only that block.
%
% IN   air_on  1 x T bool       the air-puff channel, on (> 0) per sample
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
%        loco_restpuff / _nopuff     Zhang's rests, cut at the puffs; _nopuff is true rest
%      param           1x1 struct   the settings below, and per duration dur<X> each puff's
%                                   Duration, FrequencyHz, DutyCycle

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

    % 1. make airtable
    airtable = make_airtable(air_on, air_fs, param.puff_gap_s);
    % 2. make duration field
    rounded_puffdur = round(airtable.Duration);   % P x 1 float s  the nominal puff length, 3, 5, ...
    unique_puffdur = unique(rounded_puffdur);
    dur_field = "dur" + string(unique_puffdur);   % P x 1 str  "dur3", "dur5", ... : each duration's field name

    % 3. make timetables for each puff-movement or movement-puff
    time_table = struct();
    for idx = 1:numel(dur_field)
        table_loc = rounded_puffdur == unique_puffdur(idx);
        time_table.(dur_field(idx)) = airtable{table_loc,1:2};
        param.(dur_field(idx)) = airtable(table_loc,3:end);
    end

    % 4. make move timetable 
    is_moving = gate_accel(speed, fs, param.loco_lowpass_hz, param.loco_accel_thresh);
    moving_stretch = logic2timetable(taxis, is_moving, 0);
    moving_table = bout_joingaps(moving_stretch, param.loco_join_s);

    % 5.1 make puff-move table: move with puff or stay with puff
    all_puff = [];   % P x 2 float s  every puff, in unique_puffdur order
    for idx = 1:numel(dur_field)
        puff_table = time_table.(dur_field(idx));
        [puff_go, puff_nogo, puff_moving] = split_puffs(puff_table, moving_table, ...
            param.puff_quiet_s, param.puff_quietend_s);
        time_table.("puff_go_" + dur_field(idx)) = puff_go;
        time_table.("puff_nogo_" + dur_field(idx)) = puff_nogo;
        time_table.("puff_moving_" + dur_field(idx)) = puff_moving;
        all_puff = [all_puff; puff_table]; %#ok<AGROW>
    end

    % 5.2 make move-puff table (Zhang 2026 PNAS)
    record_span = [taxis(1), taxis(end)];
    [time_table.loco_puff, time_table.loco_nopuff, time_table.loco_restnopuff, time_table.loco_restpuff] = ...
        zhang_loco(moving_table, record_span, all_puff, ...
        param.loco_quiet_s, param.loco_min_s, param.loco_rest_pad_s, param.loco_rest_min_s);
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

function timetable = logic2timetable(taxis,logic,min_duration)
%LOGIC2TIMETABLE  The runs of true in logic, as [first last] times on taxis, at least min_duration long.
% IN   taxis         1 x S float s   the time axis
%      logic         S bool          per sample
%      min_duration  1 x 1 float s   a shorter run is dropped
% OUT  timetable     N x 2 float s   [first, last] sample time of each run
    diff_binarymove = diff([0; logic(:); 0]);
    rise_idx = find(diff_binarymove == 1);
    fall_idx = find(diff_binarymove == -1);
    % Calculate duration of each bout in seconds
    bout_durations = taxis(fall_idx(:) - 1) - taxis(rise_idx(:));
    valid_bouts = bout_durations >= min_duration;
    rise_idx = rise_idx(valid_bouts);
    fall_idx = fall_idx(valid_bouts);

    % Create Nx2 array of event start and end times from taxis
    % (fall_idx - 1 gets the exact last sample of the movement bout)
    timetable = [taxis(rise_idx(:))', taxis(fall_idx(:) - 1)'];
end

function joined = bout_joingaps(bout_table, gap_s)
%BOUT_JOINGAPS  Bouts closer than gap_s merged into one.
% IN   bout_table  N x 2 float s   [start, end], sorted by start, non-overlapping
%      gap_s       1 x 1 float s   a gap shorter than this joins its two neighbours
% OUT  joined      M x 2 float s   M <= N
    joined = zeros(0, 2);
    for k = 1:size(bout_table, 1)
        if ~isempty(joined) && bout_table(k, 1) - joined(end, 2) < gap_s
            joined(end, 2) = max(joined(end, 2), bout_table(k, 2));
        else
            joined(end + 1, :) = bout_table(k, :); %#ok<AGROW>
        end
    end
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

function [loco_puff, loco_nopuff, loco_restnopuff, loco_restpuff] = zhang_loco(moving_table, record_span, all_puff, ...
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
%      loco_restnopuff, loco_restpuff   rests outside every puff (Winder's true rest), and the piece a puff landed on
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
    [loco_rest_in_window, loco_restnopuff] = ...
        bout_cutwindows(rest_table_loco, all_puff, rest_min_s);
    loco_restpuff = loco_rest_in_window(bout_beginsat(loco_rest_in_window, all_puff(:, 1)), :);
end
