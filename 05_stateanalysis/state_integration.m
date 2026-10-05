classdef state_integration < handle
    %STATE_INTEGRATION  One recording's bout tables: the sleep score into state bouts, or the puff x locomotion bouts off the analog.
    %   Caller: piv_integration_testbed, piv_transition_session, centralize_driftsession,
    %           centralize_statesession, review_session, stateintegration_test
    %   This class is designed tobe initialize time_table and hold indices
    %   corresponding to time table using time_axis of target analysis
    %   Each
    properties
        sleep_score
        binary_bin
        time_table
        param
        info
        dir_struct
    end

    methods

        function obj = state_integration(score, primary_analog)
            %STATE_INTEGRATION  Time tables for one recording, from its scoring.
            % IN   score           struct   contents of sleep_score.mat, [] if none
            %      primary_analog  struct   contents of analysis_analog.mat, [] if none
            %
            %   Takes the contents, not a session path. Which kind of recording this
            %   is comes from WHICH ARGUMENT IS FILLED, not from whether the folder
            %   name contains "sleep" or "whis", and nothing is read or created on
            %   disk -- the old constructor made <session>/state_analysis just by
            %   being called. see CLAUDE_LOG.md
            arguments
                score = []
                primary_analog = []
            end

            if ~isempty(score)
                obj.sleep_integration(score);
            end
            if ~isempty(primary_analog)
                obj.awake_integration(primary_analog);
            end
            if isempty(score) && isempty(primary_analog)
                error('state_integration:noInput', ...
                    'neither a sleep score nor an analog analysis was given');
            end
        end


        function awake_integration(obj, primary_analog)
            % primary_analog  struct   the primary_analog field of analysis_analog.mat,
            %                          already loaded. Only airtable and ball are read
            %
            % Bout tables, two anchors on one window of puff_margin_s around every puff:
            %   dur<X>                      every puff of that duration, as the airtable has it
            %   puff_go_dur<X>              a puff with a run anywhere in its window
            %   puff_nogo_dur<X>            a puff with none
            %   go_puff / go_nopuff         the runs cut at the puff windows: the pieces inside
            %                               one, and the pieces outside every one
            %   nogo_puff / nogo_nopuff     the rests likewise; a rest's window is the puff and
            %                               its after-margin only
            %   <from>_<to>_trans           rest into a run or a puff, and back; the suffix is
            %                               what get_transitionsummary selects on
            %   loco_puff / loco_nopuff     Zhang's locomotion events, by whether the onset falls in a
            %                               puff's window. see CLAUDE_LOG.md
            %   loco_rest_puff / _nopuff    Zhang's rests, cut at the puff windows; _nopuff is true rest
            %% settings, the sleep side's way : what the ball and the puffs are cut with
            obj.param.ball_smooth_s = 1;        % Gaussian window on |speed| before the threshold (s)
            obj.param.ball_thresh = 0.02;       % running above it, on the smoothed speed and at the raw edges (cm/s)
            obj.param.ball_gap_s = 3;           % a pause shorter than this joins the runs around it (s)
            obj.param.ball_run_thresh = 0.1;    % a bout has to reach this somewhere to be a run, else it is a twitch (cm/s)
            obj.param.bout_min_s = 5;           % shortest run, rest, or piece of one that is kept (s)
            obj.param.puff_margin_s = 10;       % window either side of a puff (s)
            % Zhang, Garborg et al. PNAS 2026 : locomotion on the acceleration, see CLAUDE_LOG.md
            obj.param.loco_lowpass_hz = 10;     % on the speed, before the derivative (Hz)
            obj.param.loco_accel_thresh = 0.23; % |acceleration| above it is moving (cm/s^2); their 3 through our conversion, provisional
            obj.param.loco_join_s = 1;          % moving stretches closer than this are one (s), Zhang 2022
            obj.param.loco_quiet_s = 3;         % an event needs this much not moving before it (s)
            obj.param.loco_min_s = 5;           % and lasts at least this long (s)
            obj.param.loco_rest_pad_s = 4;      % a rest starts this long after the last moving stretch (s)
            obj.param.loco_rest_min_s = 10;     % and lasts longer than this (s)
            if ~isfield(obj.param, 'transition_window')
                obj.param.transition_window = 10;   % half window of a transition (s); sleep_integration's stays if it ran first
            end
            %% Whisker stimulation time table
            atable = primary_analog.airtable;
            atable.('rounded_duration')= round(atable.Duration);
            unique_dur = unique(atable.rounded_duration);
            for idx = 1:numel(unique_dur)
                stim_dur = unique_dur(idx); % duration double 3,5,...
                table_loc = atable.rounded_duration == unique_dur(idx);
                obj.time_table.(strcat("dur",string(stim_dur))) = atable{table_loc,1:2};
                obj.param.(strcat("dur",string(stim_dur))) = atable(table_loc,3:end);
            end
            %%
            fs = primary_analog.ball.ds_fps;
            taxis = primary_analog.ball.rs_taxis;
            forward_speed = max(primary_analog.ball.resampled_ball, 0);   % a backward roll is not running
            %% 1 second smoothing
            smoothingKernel = gausswin(obj.param.ball_smooth_s*fs)/sum(gausswin(obj.param.ball_smooth_s*fs));
            ballSmooth = conv(forward_speed, smoothingKernel, 'same');
            %% Awake Integration: Continuous Movement Detection
            % 1. Detect binary movement from the SMOOTHED ball data first
            binary_smooth = ballSmooth > obj.param.ball_thresh;
            binary_move   = get_binaryball(forward_speed, binary_smooth, obj.param.ball_gap_s * fs, obj.param.ball_thresh);
            % a run has to run somewhere : a bout whose smoothed speed never reaches
            % ball_run_thresh is a twitch, and goes back to not moving
            move_edge = diff([0, binary_move(:)', 0]);
            move_rise = find(move_edge == 1);
            move_fall = find(move_edge == -1) - 1;
            for move_idx = 1:numel(move_rise)
                if max(ballSmooth(move_rise(move_idx):move_fall(move_idx))) < obj.param.ball_run_thresh
                    binary_move(move_rise(move_idx):move_fall(move_idx)) = false;
                end
            end

            % Convert binary to Nx2 [start, end] time table
            move_table = logic2timetable(taxis, binary_move, obj.param.bout_min_s);
            rest_table = logic2timetable(taxis, ~binary_move, obj.param.bout_min_s);

            %% Puff-anchored: each puff by whether a run falls in its window
            all_puff = [];
            for idx = 1:numel(unique_dur)
                dur_field = strcat("dur", string(unique_dur(idx)));
                puff_table = obj.time_table.(dur_field);
                is_go = bout_hasoverlap(puff_table, move_table, obj.param.puff_margin_s, obj.param.puff_margin_s);
                obj.time_table.(strcat("puff_go_", dur_field)) = puff_table(is_go, :);
                obj.time_table.(strcat("puff_nogo_", dur_field)) = puff_table(~is_go, :);
                all_puff = [all_puff; puff_table]; %#ok<AGROW>
            end

            %% Locomotion-anchored: each run and each rest cut at the puff windows
            % a whole-bout label made one 30-min run go_puff for nine puffs; see CLAUDE_LOG.md
            run_window = [all_puff(:, 1) - obj.param.puff_margin_s, all_puff(:, 2) + obj.param.puff_margin_s];
            rest_window = [all_puff(:, 1), all_puff(:, 2) + obj.param.puff_margin_s];
            [obj.time_table.go_puff, obj.time_table.go_nopuff] = ...
                bout_cutwindows(move_table, run_window, obj.param.bout_min_s);
            [obj.time_table.nogo_puff, obj.time_table.nogo_nopuff] = ...
                bout_cutwindows(rest_table, rest_window, obj.param.bout_min_s);

            %% Zhang's locomotion, on the acceleration : the events and the rests PNAS 2026 aligns to, split by puff
            % both directions count, as there; the speed-defined go / nogo tables above are a different question
            [lp_zero, lp_pole, lp_gain] = butter(5, obj.param.loco_lowpass_hz / (fs / 2), 'low');
            [lp_sos, lp_g] = zp2sos(lp_zero, lp_pole, lp_gain);
            speed_lp = filtfilt(lp_sos, lp_g, primary_analog.ball.resampled_ball(:)');
            is_moving = abs([0, diff(speed_lp)] * fs) >= obj.param.loco_accel_thresh;
            moving_table = bout_joingaps(logic2timetable(taxis, is_moving, 0), obj.param.loco_join_s);
            quiet_before = moving_table(:, 1) - [taxis(1); moving_table(1:end-1, 2)];
            long_enough = moving_table(:, 2) - moving_table(:, 1) >= obj.param.loco_min_s;
            loco_event = moving_table(quiet_before >= obj.param.loco_quiet_s & long_enough, :);
            % an event whose onset falls in a puff's window is the puff's; the rest are spontaneous
            onset_in_puff = false(size(loco_event, 1), 1);
            for event_idx = 1:size(loco_event, 1)
                onset_in_puff(event_idx) = any(run_window(:, 1) <= loco_event(event_idx, 1) & loco_event(event_idx, 1) <= run_window(:, 2));
            end
            obj.time_table.loco_puff = loco_event(onset_in_puff, :);
            obj.time_table.loco_nopuff = loco_event(~onset_in_puff, :);
            % a rest : from the pad after one moving stretch to the start of the next, the record's ends included
            rest_start = [taxis(1); moving_table(:, 2) + obj.param.loco_rest_pad_s];
            rest_end = [moving_table(:, 1); taxis(end)];
            rest_table_loco = [rest_start, rest_end];
            rest_table_loco = rest_table_loco(rest_end - rest_start > obj.param.loco_rest_min_s, :);
            % a rest a puff reaches is the puff's; what lies outside every puff window is Winder's true rest
            [obj.time_table.loco_rest_puff, obj.time_table.loco_rest_nopuff] = ...
                bout_cutwindows(rest_table_loco, rest_window, obj.param.loco_rest_min_s);

            %% Transitions out of and back into the plain rest
            trans_win = obj.param.transition_window;
            obj.time_table.nogo_nopuff_go_nopuff_trans = ...
                get_transition(trans_win, obj.time_table.nogo_nopuff, obj.time_table.go_nopuff, 0.5);
            obj.time_table.go_nopuff_nogo_nopuff_trans = ...
                get_transition(trans_win, obj.time_table.go_nopuff, obj.time_table.nogo_nopuff, 0.5);
            for idx = 1:numel(unique_dur)
                dur_field = strcat("dur", string(unique_dur(idx)));
                for puff_kind = ["puff_go_", "puff_nogo_"]
                    puff_field = strcat(puff_kind, dur_field);
                    obj.time_table.(strcat("nogo_nopuff_", puff_field, "_trans")) = ...
                        get_transition(trans_win, obj.time_table.nogo_nopuff, obj.time_table.(puff_field), 0.5);
                    obj.time_table.(strcat(puff_field, "_nogo_nopuff_trans")) = ...
                        get_transition(trans_win, obj.time_table.(puff_field), obj.time_table.nogo_nopuff, 0.5);
                end
            end

        end

        function sleep_integration(obj, score)
            % score  struct   the contents of sleep_score.mat, already loaded
            obj.sleep_score = score;
            obj.param.transition_window = 25;
            obj.param.bigchunk_windowsize = 300;            obj.param.bigchunk_awake_weight = struct('w_awake', 10, 'w_nrem', -1, 'w_rem', -1000, 'w_drowsy', 0);
            obj.param.bigchunk_nrem_weight = struct('w_awake', -1, 'w_nrem', 10, 'w_rem', -1000, 'w_drowsy', 0);
            obj.param.bigchunk_rem_weight = struct('w_awake', -1, 'w_nrem', 0, 'w_rem', 10, 'w_drowsy', 0);
            obj.param.bigchunk_drowsy_weight = struct('w_awake', 0, 'w_nrem', -1, 'w_rem', -1000, 'w_drowsy', 10);

            % Binary array generation
            obj.binary_bin.awake = score.behavState == score.statecodes.NotSleep;
            obj.binary_bin.nrem = score.behavState == score.statecodes.NREM;
            obj.binary_bin.rem = score.behavState == score.statecodes.REM;
            obj.binary_bin.drowsy = score.behavState == score.statecodes.Drowsy;
            [obj.time_table.long_nrem, obj.info.bigchunk_nrem_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_nrem_weight,...
                obj.param.bigchunk_windowsize, score.binwidth_sec);
            [obj.time_table.long_awake, obj.info.bigchunk_awake_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_awake_weight,...
                obj.param.bigchunk_windowsize, score.binwidth_sec);
            [obj.time_table.long_drowsy, obj.info.bigchunk_drowsy_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_drowsy_weight,...
                obj.param.bigchunk_windowsize, score.binwidth_sec);
            [obj.time_table.long_rem, obj.info.bigchunk_rem_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_rem_weight,...
                obj.param.bigchunk_windowsize, score.binwidth_sec);
            % Time table generation
            % for chunk analysis
            obj.time_table.roughawake = statebin2timetable(obj.sleep_score.AwakeTimes, obj.sleep_score.DrowsyTimes);
            obj.time_table.roughnrem = statebin2timetable(obj.sleep_score.NREMTimes, obj.sleep_score.uArousalTimes);
            obj.time_table.rem = statebin2timetable(obj.sleep_score.REMTimes);
            % not actively using
            obj.time_table.awake = statebin2timetable(obj.sleep_score.AwakeTimes);
            obj.time_table.drowsy = statebin2timetable(obj.sleep_score.DrowsyTimes);
            obj.time_table.uarousal = statebin2timetable(obj.sleep_score.uArousalTimes);
            obj.time_table.nrem = statebin2timetable(obj.sleep_score.NREMTimes);
            % transition
            obj.time_table.nr_trans = get_transition(obj.param.transition_window,obj.time_table.roughnrem, obj.time_table.rem);
            obj.time_table.na_trans = get_transition(obj.param.transition_window,obj.time_table.roughnrem, obj.time_table.roughawake);
            obj.time_table.ra_trans = get_transition(obj.param.transition_window,obj.time_table.rem, obj.time_table.roughawake);
            obj.time_table.an_trans = get_transition(obj.param.transition_window,obj.time_table.roughawake, obj.time_table.roughnrem);
        end

        function state_idx = add_taxis(obj,t_axis)
            % from time axis (sec) get
            state_idx = struct();
            ttable_names = fieldnames(obj.time_table);
            for tname_idx = 1:numel(ttable_names)
                ttable_name = ttable_names{tname_idx};
                disp(ttable_name)
                state_idx.(ttable_name) = timetable2frame(t_axis,obj.time_table.(ttable_name));
            end
        end

        function trim_to_duration(obj, max_time_sec)
            % TRIM_TO_DURATION  이미징이 중간에 중단된 경우, sleep_score(peripheral)이
            % 더 길 수 있으므로 짧은 쪽(이미징 끝시간)에 맞게 time_table을 자른다.
            % transition_window 1개 분량을 safety margin으로 추가로 자름 ->
            % get_transitionsummary에서 center ± onset 윈도우가 절대 넘치지 않음.
            %
            % Usage:
            %   state_integrate.trim_to_duration(pax_fwhm.t_axis(end))
            %
            % Input:
            %   max_time_sec - 유지할 최대 시간 (초), 보통 pax_fwhm.t_axis(end)

            % Safety margin: transition_window 1개 (없으면 보수적으로 30초)
            if isfield(obj.param, 'transition_window')
                safety_margin = obj.param.transition_window;
            else
                safety_margin = 30;
            end
            effective_max = max(0, max_time_sec - safety_margin);

            fprintf('[trim_to_duration] imaging=%.1f sec, margin=%.1f sec -> clip to %.1f sec\n', ...
                max_time_sec, safety_margin, effective_max);

            %% 1. Clip binary_bin arrays (bin 단위 -> 초 변환 필요)
            if ~isempty(obj.sleep_score) && isfield(obj.sleep_score, 'binwidth_sec')
                max_bin = floor(effective_max / obj.sleep_score.binwidth_sec);
                bin_fields = fieldnames(obj.binary_bin);
                for i = 1:numel(bin_fields)
                    arr = obj.binary_bin.(bin_fields{i});
                    if numel(arr) > max_bin
                        obj.binary_bin.(bin_fields{i}) = arr(1:max_bin);
                    end
                end
                fprintf('  binary_bin clipped to %d bins (%.1f sec)\n', max_bin, max_bin * obj.sleep_score.binwidth_sec);
            end

            %% 2. Clip time_table Nx2 matrices (초 단위)
            if isempty(obj.time_table)
                return;
            end
            tfields = fieldnames(obj.time_table);
            for i = 1:numel(tfields)
                tbl = obj.time_table.(tfields{i});
                if isempty(tbl)
                    continue;
                end
                % start > effective_max 인 row 제거
                valid_rows = tbl(:, 1) < effective_max;
                tbl = tbl(valid_rows, :);
                % end > effective_max 인 row의 end를 effective_max로 클립
                tbl(tbl(:, 2) > effective_max, 2) = effective_max;
                obj.time_table.(tfields{i}) = tbl;
            end
            fprintf('  time_table fields clipped: %s\n', strjoin(tfields', ', '));
        end

    end
end
