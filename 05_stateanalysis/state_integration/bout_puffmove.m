function [time_table, param] = bout_puffmove(air_on, air_fs, speed, taxis, fs)
%BOUT_PUFFMOVE  The puffs crossed with the movement, as bout tables.
%   Caller: centralize_bout, b_reviewsession
%   The movement block reads the ball; categorize_puff, categorize_locomotion and categorize_rest cross it with the puffs and
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
%        puff_gogo_dur<X>            a puff to a mouse already moving before it, moving during it
%        loco_puff / loco_nopuff     Zhang's locomotion events, by whether the onset falls
%                                    during a puff
%        rest_puff / _nopuff     Zhang's rests, cut at the puffs; _nopuff is true rest
%      param           1x1 struct   the settings below, and per duration dur<X> each puff's
%                                   Duration, FrequencyHz, DutyCycle

    %% Criteria
    % Puff
    param.puff_gap_s = 10;          % pulse onsets further apart than this are separate puffs (s), 02_othersignal get_airtable
    param.puff_quiet_s = 3;        % puff criterion: still for this long before onset, else puff_gogo, or dropped when still through the puff (s), Zhang 2026
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
        [puff_go, puff_nogo, puff_gogo] = categorize_puff(puff_table, moving_table, ...
            param.puff_quiet_s, param.puff_quietend_s);
        time_table.("puff_go_" + dur_field(idx)) = puff_go;
        time_table.("puff_nogo_" + dur_field(idx)) = puff_nogo;
        time_table.("puff_gogo_" + dur_field(idx)) = puff_gogo;
        all_puff = [all_puff; puff_table]; %#ok<AGROW>
    end

    % 5.2 make move-puff table (Zhang 2026 PNAS)
    [time_table.loco_puff, time_table.loco_nopuff] = categorize_locomotion(moving_table, all_puff, ...
        param.loco_quiet_s, param.loco_min_s);
    % 5.3 make rest tables, cut at the puffs
    [time_table.rest_nopuff, time_table.rest_puff] = categorize_rest(moving_table, all_puff, ...
        param.loco_rest_pad_s, param.loco_rest_min_s);
end

function is_moving = gate_accel(speed, fs, lowpass_hz, accel_thresh)
%GATE_ACCEL  Zhang's moving samples: the low-passed speed's |acceleration| at or over accel_thresh.
% IN   speed         T float cm/s         the ball speed, signed
%      fs            1 x 1 float Hz       its rate
%      lowpass_hz    1 x 1 float Hz       on the speed, before the derivative
%      accel_thresh  1 x 1 float cm/s^2
% OUT  is_moving     1 x T bool
    % 1. Low pass
    [lp_zero, lp_pole, lp_gain] = butter(5, lowpass_hz / (fs / 2), 'low'); 
    [lp_sos, lp_g] = zp2sos(lp_zero, lp_pole, lp_gain); 
    speed_lp = filtfilt(lp_sos, lp_g, speed(:)');
    % 2. Calculate acceleration
    accel = gradient(speed_lp) * fs;
    % 3. Threshold to make gate
    is_moving = abs(accel) >= accel_thresh;
end

function timetable = logic2timetable(taxis,logic,min_duration)
%LOGIC2TIMETABLE  The runs of true in logic, as [first last] times on taxis, at least min_duration long.
% IN   taxis         1 x S float s   the time axis
%      logic         S bool          per sample
%      min_duration  1 x 1 float s   a shorter run is dropped
% OUT  timetable     N x 2 float s   [first, last] sample time of each run
    % 1. zero-pad to find start from 1 and end with one
    diff_binarymove = diff([0; logic(:); 0]);
    % 2. Detect rising and falling edge
    rise_idx = find(diff_binarymove == 1);
    fall_idx = find(diff_binarymove == -1);
    % 3. Calculate gap between rise and fall edge
    bout_durations = taxis(fall_idx(:) - 1) - taxis(rise_idx(:));
    % 4. Pulse purification
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
    if isempty(bout_table)
        joined = zeros(0, 2);
        return
    end
    joined = bout_table(1, :);
    for k = 2:size(bout_table, 1)
        gap = bout_table(k, 1) - joined(end, 2);   % s  from the last joined bout's end
        is_shorter = gap < gap_s;
        if is_shorter
            joined(end, 2) = max(joined(end, 2), bout_table(k, 2));
        else
            joined(end + 1, :) = bout_table(k, :); %#ok<AGROW>
        end
    end
end

function [puff_go, puff_nogo, puff_gogo] = categorize_puff(puff_table, moving_table, quiet_s, quietend_s)
%CATEGORIZE_PUFF  One duration's puffs, by whether a moving stretch falls before them or during them.
%   A puff to a mouse moving before it and still through it is dropped: when it stopped is unknown.
% IN   puff_table    N x 2 float s   one duration's puffs, [onset end]
%      moving_table  L x 2 float s   the moving stretches
%      quiet_s       1 x 1 float s   param.puff_quiet_s
%      quietend_s    1 x 1 float s   param.puff_quietend_s
% OUT  puff_go       K x 2 float s   to a mouse not moving before, moving during the puff
%      puff_nogo     J x 2 float s   to a mouse not moving before, with none
%      puff_gogo     I x 2 float s   to a mouse moving before, moving during the puff
    % 1. make buffer time table - quiet before puff
    quiet_table = [puff_table(:, 1) - quiet_s, puff_table(:, 1) - quietend_s];   % N x 2 float s  [onset - quiet_s, onset - quietend_s]
    % 2. gate1 does the mouse move before puff 
    is_moving_before = bout_hasoverlap(quiet_table, moving_table, 0, 0);
    % 3. gate2 does the mouse run during puff
    is_go = bout_hasoverlap(puff_table, moving_table, 0, 0);
    % 4. gate1xgate2 does the mouse not move and moved during puff or not 
    is_puffgo = ~is_moving_before & is_go;
    is_puffnogo = ~is_moving_before & ~is_go;
    is_puffgogo = is_moving_before & is_go;
    % out
    puff_go = puff_table(is_puffgo, :);
    puff_nogo = puff_table(is_puffnogo, :);
    puff_gogo = puff_table(is_puffgogo, :);
end

function [loco_puff, loco_nopuff] = categorize_locomotion(moving_table, all_puff, quiet_s, min_s)
%CATEGORIZE_LOCOMOTION  Zhang's locomotion events on the moving stretches, by whether the onset falls during a puff.
%   The events PNAS 2026 aligns to. The first stretch is never one: nothing before it to measure quiet from.
% IN   moving_table  L x 2 float s   the moving stretches
%      all_puff      P x 2 float s   every puff, [onset end]
%      quiet_s       1 x 1 float s   param.loco_quiet_s
%      min_s         1 x 1 float s   param.loco_min_s
% OUT  loco_puff     K x 2 float s   events whose onset falls during a puff
%      loco_nopuff   J x 2 float s   the rest, spontaneous
    % 1. make buffer timetable - quiet time before move
    candidate = moving_table(2:end, :);   % (L-1) x 2  drop first movement with unknown prior state
    quiet_before = candidate(:, 1) - moving_table(1:end-1, 2);   % (L-1) x 1  gap since the stretch before
    % 2. Duration gate
    long_enough = candidate(:, 2) - candidate(:, 1) >= min_s;
    % 3. Filter with duration and quiet
    loco_event = candidate(quiet_before >= quiet_s & long_enough, :);
    % 4. Categorize overlap with puff
    onset_in_puff = false(size(loco_event, 1), 1);
    for event_idx = 1:size(loco_event, 1)
        onset_in_puff(event_idx) = any(all_puff(:, 1) <= loco_event(event_idx, 1) & loco_event(event_idx, 1) <= all_puff(:, 2));
    end
    loco_puff = loco_event(onset_in_puff, :);
    loco_nopuff = loco_event(~onset_in_puff, :);
end

function [rest_nopuff, rest_puff] = categorize_rest(moving_table, all_puff, pad_s, min_s)
%CATEGORIZE_REST  Zhang's rests between the moving stretches, cut at the puffs.
%   The rests PNAS 2026 aligns to. Only the gaps between two stretches: the record's ends are not rests.
% IN   moving_table     L x 2 float s   the moving stretches
%      all_puff         P x 2 float s   every puff, [onset end]
%      pad_s            1 x 1 float s   param.loco_rest_pad_s
%      min_s            1 x 1 float s   param.loco_rest_min_s
% OUT  rest_nopuff  Q x 2 float s   rest outside every puff (Winder's true rest)
%      rest_puff    R x 2 float s   the piece of a rest a puff landed on
    % 1. make rest time table 
    rest_start = moving_table(1:end-1, 2) + pad_s;   % (L-1) x 1  rest start after stop move + pad, drop last moveend
    rest_end = moving_table(2:end, 1);               % (L-1) x 1  rest end with start move
    rest_table = [rest_start, rest_end];
    rest_table = rest_table(rest_end - rest_start > min_s, :);

    % 2. rest table touch puff or not touching puff
    [rest_puff, rest_nopuff] = ...
        bout_cutwindows(rest_table, all_puff, min_s);
    % 3. rest touch puff that 
    rest2puff_gap = abs(rest_puff(:, 1) - all_puff(:, 1)');   % N x P float s  each piece's start against each puff onset
    is_at_onset = rest2puff_gap < 1e-3;                            % N x P bool     within a millisecond
    begins_at_puff = any(is_at_onset, 2);                      % N x 1 bool     at any puff's onset
    rest_puff = rest_puff(begins_at_puff, :);
end

function hit = bout_hasoverlap(bout_table, other_table, pre_s, post_s)
%BOUT_HASOVERLAP  Per row of bout_table, whether any row of other_table touches [start - pre_s, end + post_s].
%
% IN   bout_table   N x 2 float s   [start, end] bouts to label
%      other_table  M x 2 float s   [start, end] bouts to look for
%      pre_s        1 x 1 float s   window extension before each bout
%      post_s       1 x 1 float s   window extension after each bout
% OUT  hit          N x 1 bool
    hit = false(size(bout_table, 1), 1);
    if isempty(bout_table) || isempty(other_table)
        return
    end
    window_start = bout_table(:, 1) - pre_s;
    window_end = bout_table(:, 2) + post_s;
    for k = 1:size(bout_table, 1)
        hit(k) = any(other_table(:, 1) <= window_end(k) & other_table(:, 2) >= window_start(k));
    end
end

function [inside, outside] = bout_cutwindows(bout_table, window_table, min_s)
%BOUT_CUTWINDOWS  Each bout split at the windows: the pieces inside any window, and the pieces outside all.
%
% IN   bout_table    N x 2 float s   [start, end] bouts to cut
%      window_table  M x 2 float s   [start, end] windows; may overlap one another
%      min_s         1 x 1 float s   a piece shorter than this is dropped
% OUT  inside        P x 2 float s   pieces that lie in a window
%      outside       Q x 2 float s   pieces that lie in no window
    inside = zeros(0, 2);
    outside = zeros(0, 2);
    if isempty(bout_table)
        return
    end
    if isempty(window_table)
        outside = bout_table;
        return
    end
    % windows merged into disjoint spans, so each bout is walked once
    sorted_window = sortrows(window_table, 1);
    span = sorted_window(1, :);
    for k = 2:size(sorted_window, 1)
        if sorted_window(k, 1) <= span(end, 2)
            span(end, 2) = max(span(end, 2), sorted_window(k, 2));
        else
            span(end+1, :) = sorted_window(k, :); %#ok<AGROW>
        end
    end
    for k = 1:size(bout_table, 1)
        bout_start = bout_table(k, 1);
        bout_end = bout_table(k, 2);
        cursor = bout_start;
        touching = span(span(:, 2) > bout_start & span(:, 1) < bout_end, :);
        for j = 1:size(touching, 1)
            piece_start = max(touching(j, 1), bout_start);
            piece_end = min(touching(j, 2), bout_end);
            if piece_start - cursor >= min_s
                outside(end+1, :) = [cursor, piece_start]; %#ok<AGROW>
            end
            if piece_end - piece_start >= min_s
                inside(end+1, :) = [piece_start, piece_end]; %#ok<AGROW>
            end
            cursor = piece_end;
        end
        if bout_end - cursor >= min_s
            outside(end+1, :) = [cursor, bout_end]; %#ok<AGROW>
        end
    end
end
