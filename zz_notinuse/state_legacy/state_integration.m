classdef state_integration < handle
    %STATE_INTEGRATION  One recording's bout tables: the sleep score into state bouts, or the puff x locomotion bouts off the analog.
    %   Caller: centralize_boutsession, piv_integration_testbed, analysis_main_piv,
    %           review_session, stateintegration_test
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
            % Bout tables, two anchors on one window per puff, the puff itself and puff_post_s after it:
            %   dur<X>                      every puff of that duration, as the airtable has it
            %   puff_go_dur<X>              a puff to a mouse not moving in [onset - puff_quiet_s,
            %                               onset - puff_quietend_sec],
            %                               with movement anywhere in its window. Movement is Zhang's
            %                               |acceleration|, the same stretches loco_* are built from
            %   puff_nogo_dur<X>            the same, with none
            %   puff_moving_dur<X>          a puff to a mouse already moving in that window
            %   go_puff / go_nopuff         a run begun in a window is the puff's whole, until the mouse
            %                               stops; a run already going when the puff came is cut at the
            %                               window edges and only the window is the puff's. see CLAUDE_LOG.md
            %   nogo_puff / nogo_nopuff     the rests cut at the window edges: _puff is the in-window
            %                               piece of a rest the puff landed on, _nopuff the pieces outside
            %                               every window. A rest begun inside a window is in neither
            %   <from>_<to>_trans           rest into a run or a puff, and back; the suffix is
            %                               what get_transitionsummary selects on
            %   loco_puff / loco_nopuff     Zhang's locomotion events, by whether the onset falls in a
            %                               puff's window. see CLAUDE_LOG.md
            %   loco_rest_puff / _nopuff    Zhang's rests, cut at the puff windows; _nopuff is true rest
            %% settings, the sleep side's way : what the ball and the puffs are cut with
            obj.param.ball_smooth_s = 1;        % Gaussian window on |speed| before the threshold (s)
            obj.param.ball_thresh = 0.26;       % running above it, on the smoothed speed and at the raw edges (cm/s)
            obj.param.ball_gap_s = 3;           % a pause shorter than this joins the runs around it (s)
            obj.param.ball_run_thresh = 1.31;   % a bout has to reach this somewhere to be a run, else it is a twitch (cm/s)
            obj.param.bout_min_s = 5;           % shortest run, rest, or piece of one that is kept (s)
            obj.param.puff_post_s = 10;         % a puff's window runs from its onset to this long after its end (s)
            obj.param.puff_quiet_s = 3;         % a puff is a trial only if no moving stretch touches this long before its onset (s), Zhang 2026
            obj.param.puff_quietend_sec = 0.1;  % that window ends this long before the onset (s); the zero-phase filters pull an evoked onset back, see CLAUDE_LOG.md
            % Zhang, Garborg et al. PNAS 2026 : locomotion on the acceleration, see CLAUDE_LOG.md
            obj.param.loco_lowpass_hz = 10;     % on the speed, before the derivative (Hz)
            obj.param.loco_accel_thresh = 3;    % |acceleration| above it is moving (cm/s^2), Zhang 2022; the three thresholds moved x13.09 when the ETACH2 scale was corrected, see CLAUDE_LOG.md
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

            %% Zhang's moving stretches, on the acceleration : both directions count, as there
            [lp_zero, lp_pole, lp_gain] = butter(5, obj.param.loco_lowpass_hz / (fs / 2), 'low');
            [lp_sos, lp_g] = zp2sos(lp_zero, lp_pole, lp_gain);
            speed_lp = filtfilt(lp_sos, lp_g, primary_analog.ball.resampled_ball(:)');
            is_moving = abs([0, diff(speed_lp)] * fs) >= obj.param.loco_accel_thresh;
            moving_table = bout_joingaps(logic2timetable(taxis, is_moving, 0), obj.param.loco_join_s);

            %% Puff-anchored: a puff to a still mouse, by whether a moving stretch falls in its window
            all_puff = [];
            for idx = 1:numel(unique_dur)
                dur_field = strcat("dur", string(unique_dur(idx)));
                puff_table = obj.time_table.(dur_field);
                quiet_end = puff_table(:, 1) - obj.param.puff_quietend_sec;   % N x 1 float s
                quiet_table = [quiet_end, quiet_end];   % N x 2 float s  the window's end as a point
                quiet_span = obj.param.puff_quiet_s - obj.param.puff_quietend_sec;
                is_moving_before = bout_hasoverlap(quiet_table, moving_table, quiet_span, 0);
                is_go = bout_hasoverlap(puff_table, moving_table, 0, obj.param.puff_post_s);
                keep_go = ~is_moving_before & is_go;
                keep_nogo = ~is_moving_before & ~is_go;
                obj.time_table.(strcat("puff_go_", dur_field)) = puff_table(keep_go, :);
                obj.time_table.(strcat("puff_nogo_", dur_field)) = puff_table(keep_nogo, :);
                obj.time_table.(strcat("puff_moving_", dur_field)) = puff_table(is_moving_before, :);
                all_puff = [all_puff; puff_table]; %#ok<AGROW>
            end

            %% Locomotion-anchored: a run begun in a window is the puff's whole, any other is cut at the window; rests are cut
            % a whole-bout label made one 30-min run go_puff for nine puffs, and a margin before the
            % puff handed spontaneous runs to it; see CLAUDE_LOG.md
            puff_window = [all_puff(:, 1), all_puff(:, 2) + obj.param.puff_post_s];
            [obj.time_table.go_puff, obj.time_table.go_nopuff] = ...
                bout_splitat(move_table, puff_window, obj.param.bout_min_s);
            % a rest the puff landed on begins, once cut, at the puff onset; a rest begun inside a
            % window, after a puff-run, is neither the puff's nor true rest. see CLAUDE_LOG.md
            [rest_in_window, obj.time_table.nogo_nopuff] = bout_cutwindows(rest_table, puff_window, obj.param.bout_min_s);
            obj.time_table.nogo_puff = rest_in_window(bout_beginsat(rest_in_window, all_puff(:, 1)), :);

            %% Zhang's locomotion, on the moving stretches : the events and the rests PNAS 2026 aligns to, split by puff
            % the speed-defined go / nogo tables above are a different question
            quiet_before = moving_table(:, 1) - [taxis(1); moving_table(1:end-1, 2)];
            long_enough = moving_table(:, 2) - moving_table(:, 1) >= obj.param.loco_min_s;
            loco_event = moving_table(quiet_before >= obj.param.loco_quiet_s & long_enough, :);
            % an event whose onset falls in a puff's window is the puff's; the rest are spontaneous
            onset_in_puff = false(size(loco_event, 1), 1);
            for event_idx = 1:size(loco_event, 1)
                onset_in_puff(event_idx) = any(puff_window(:, 1) <= loco_event(event_idx, 1) & loco_event(event_idx, 1) <= puff_window(:, 2));
            end
            obj.time_table.loco_puff = loco_event(onset_in_puff, :);
            obj.time_table.loco_nopuff = loco_event(~onset_in_puff, :);
            % a rest : from the pad after one moving stretch to the start of the next, the record's ends included
            rest_start = [taxis(1); moving_table(:, 2) + obj.param.loco_rest_pad_s];
            rest_end = [moving_table(:, 1); taxis(end)];
            rest_table_loco = [rest_start, rest_end];
            rest_table_loco = rest_table_loco(rest_end - rest_start > obj.param.loco_rest_min_s, :);
            % a rest the puff landed on is the puff's; what lies outside every puff window is Winder's true rest
            [loco_rest_in_window, obj.time_table.loco_rest_nopuff] = ...
                bout_cutwindows(rest_table_loco, puff_window, obj.param.loco_rest_min_s);
            obj.time_table.loco_rest_puff = loco_rest_in_window(bout_beginsat(loco_rest_in_window, all_puff(:, 1)), :);

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

    end
end
