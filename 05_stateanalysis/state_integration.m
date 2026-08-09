classdef state_integration < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
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

        function obj = state_integration(base_path)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            [~, folder_name] = fileparts(base_path);
            obj.dir_struct.stateanalysis = fullfile(base_path,"state_analysis");
            if ~exist(obj.dir_struct.stateanalysis, 'dir')
                mkdir(obj.dir_struct.stateanalysis);
                disp(['Created directory: ', obj.dir_struct.stateanalysis]);
            end
            sleepscore_dir = fullfile(base_path,"peripheral","sleep_score.mat");
            analog_dir =     fullfile(base_path,"peripheral","analysis_analog.mat");
            if isfile(sleepscore_dir)
                fprintf('sleep_score.mat detected\n')
                obj.sleep_integration(sleepscore_dir);
            elseif contains(folder_name,'sleep','IgnoreCase',true)
                fprintf('sleepscoring is not done')
                return
            elseif contains(folder_name,'whis','IgnoreCase',true)
                if isfile(analog_dir)
                    fprintf('analog.mat detected with whisker stim folder')
                    obj.awake_integration(analog_dir)
                else
                    fprintf('peripheral analysis is not done')
                end
            else
                fprintf('Folder does not contain sleep or whis, check folder name')
                return
            end
            %%
        end


        function awake_integration(obj,analog_dir)
            %% Whisker stimulation time table
            obj.dir_struct.analog = analog_dir;
            loadstruct = load(analog_dir);
            primary_analog = loadstruct.primary_analog;
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
            absball = abs(primary_analog.ball.resampled_ball);
            %% 1 second smoothing
            kernelWidth = 1; % window (sec)
            smoothingKernel = gausswin(kernelWidth*fs)/sum(gausswin(kernelWidth*fs));
            ballSmooth = conv(absball,smoothingKernel,'same');
            %% Awake Integration: Continuous Movement Detection
            % 1. Detect binary movement from the SMOOTHED ball data first
            smooth_threshold = 0.02;
            binary_smooth    = abs(ballSmooth) > smooth_threshold;
            gap_tolerance    = 3 * fs;
            binary_move      = get_binaryball(absball, binary_smooth, gap_tolerance);

            % Convert binary to Nx2 [start, end] time table (min bout = 5 s)
            move_table = logic2timetable(taxis, binary_move, 5);
            %%
            obj.time_table.move = move_table;

            %% Per-duration merge: merged_durX = stim start + max(stim end, move end)
            all_air = [];
            all_merged = [];
            for idx = 1:numel(unique_dur)
                stim_dur  = unique_dur(idx);
                dur_field = strcat("dur", string(stim_dur));
                dur_air   = obj.time_table.(dur_field);  % Nx2 for this duration

                % Merge this duration's stim intervals with movement table
                merged_dur = merge_airtable_movtable(dur_air, move_table, 10);
                obj.time_table.(strcat("merged_", dur_field)) = merged_dur;  % ← stored per duration

                % Accumulate all stim intervals for pure_movement filtering
                all_air = [all_air; dur_air]; %#ok<AGROW>

                % Accumulate all merged intervals for all_merged_airpuff
                if ~isempty(merged_dur)
                    all_merged = [all_merged; merged_dur]; %#ok<AGROW>
                end
            end

            % Save all merged air puffs sorted by start time
            if ~isempty(all_merged)
                [~, sort_idx] = sort(all_merged(:, 1));
                all_merged = all_merged(sort_idx, :);
            end
            obj.time_table.all_merged_airpuff = all_merged;

            %% pure_movement: all movement bouts that DO NOT overlap with any merged airpuff
            is_pure = true(size(move_table, 1), 1);
            for m = 1:size(move_table, 1)
                m_start = move_table(m, 1);
                m_end   = move_table(m, 2);
                for a = 1:size(all_merged, 1)
                    a_start = all_merged(a, 1);
                    a_end   = all_merged(a, 2);
                    % Check overlap (strict overlap)
                    if m_start <= a_end && m_end >= a_start
                        is_pure(m) = false;
                        break;
                    end
                end
            end
            obj.time_table.pure_movement = move_table(is_pure, :);

            %% static: inverse of binary move (minimum 5 seconds of rest)
            obj.time_table.static = logic2timetable(taxis, ~binary_move, 5);

            %% Transitions
            if ~isfield(obj.param, 'transition_window')
                obj.param.transition_window = 10;
            end
            trans_win = obj.param.transition_window;

            % pure <-> static transitions
            obj.time_table.trans_static2pure = get_transition(trans_win, obj.time_table.static, obj.time_table.pure_movement, 0.5);
            obj.time_table.trans_pure2static = get_transition(trans_win, obj.time_table.pure_movement, obj.time_table.static, 0.5);

            % merged_durX <-> static transitions
            for idx = 1:numel(unique_dur)
                stim_dur = unique_dur(idx);
                merged_field = strcat("merged_dur", string(stim_dur));

                if isfield(obj.time_table, merged_field)
                    obj.time_table.(strcat("static_",merged_field,"_trans")) = ...
                        get_transition(trans_win, obj.time_table.static, obj.time_table.(merged_field), 0.5);

                    obj.time_table.(strcat(merged_field, "_static","_trans")) = ...
                        get_transition(trans_win, obj.time_table.(merged_field), obj.time_table.static, 0.5);
                end
            end

            %% Diagnostic plots
            figure()
            plot(taxis,binary_move)

            %%
            cla
            plot(taxis,ballSmooth)
            hold on
            plot(taxis,absball,'Color','k')
            %%
            plot(taxis, binary_move, 'y', 'LineWidth', 1.5)
            hold on
            plot(taxis, binary_move, 'r', 'LineWidth', 1.5)

            plot(taxis, absball > 0.02, 'k:', 'LineWidth', 1)
            legend('Interpolated (final_binary_movement)', 'Raw (absball > 0.02)')
            %%
            % Generate a color map for the different durations
            y_pos = 0;

            % Loop through each duration to draw its bars
            for idx = 1:numel(unique_dur)
                stim_dur = unique_dur(idx);

                % Retrieve the Nx2 table/matrix of [StartTime, EndTime]
                time_intervals = obj.time_table.(strcat("dur", string(stim_dur)));

                % Fast plotting strategy: Create vectors separated by NaNs
                % so we only call plot() once per duration type.
                num_intervals = size(time_intervals, 1);
                x_data = NaN(3, num_intervals);
                x_data(1, :) = time_intervals(:, 1)'; % StartTime
                x_data(2, :) = time_intervals(:, 2)'; % EndTime

                y_data = NaN(3, num_intervals);
                y_data(1:2, :) = y_pos; % Set the line height to 2.5

                % Plot all horizontal bars for this duration type
                plot(x_data(:), y_data(:), 'Color', 'g', ...
                    'LineWidth', 4, 'DisplayName', sprintf('Duration %d s', stim_dur));
            end
        end

        function sleep_integration(obj,sleepscore_dir)
            obj.dir_struct.sleep_score = sleepscore_dir;
            sleep_score = load(obj.dir_struct.sleep_score);
            obj.sleep_score = sleep_score;
            obj.param.transition_window = 25;
            obj.param.bigchunk_windowsize = 300;            obj.param.bigchunk_awake_weight = struct('w_awake', 10, 'w_nrem', -1, 'w_rem', -1000, 'w_drowsy', 0);
            obj.param.bigchunk_nrem_weight = struct('w_awake', -1, 'w_nrem', 10, 'w_rem', -1000, 'w_drowsy', 0);
            obj.param.bigchunk_rem_weight = struct('w_awake', -1, 'w_nrem', 0, 'w_rem', 10, 'w_drowsy', 0);
            obj.param.bigchunk_drowsy_weight = struct('w_awake', 0, 'w_nrem', -1, 'w_rem', -1000, 'w_drowsy', 10);

            % Binary array generation
            obj.binary_bin.awake = sleep_score.behavState == sleep_score.statecodes.NotSleep;
            obj.binary_bin.nrem = sleep_score.behavState == sleep_score.statecodes.NREM;
            obj.binary_bin.rem = sleep_score.behavState == sleep_score.statecodes.REM;
            obj.binary_bin.drowsy = sleep_score.behavState == sleep_score.statecodes.Drowsy;
            [obj.time_table.long_nrem, obj.info.bigchunk_nrem_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_nrem_weight,...
                obj.param.bigchunk_windowsize, sleep_score.binwidth_sec);
            [obj.time_table.long_awake, obj.info.bigchunk_awake_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_awake_weight,...
                obj.param.bigchunk_windowsize, sleep_score.binwidth_sec);
            [obj.time_table.long_drowsy, obj.info.bigchunk_drowsy_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_drowsy_weight,...
                obj.param.bigchunk_windowsize, sleep_score.binwidth_sec);
            [obj.time_table.long_rem, obj.info.bigchunk_rem_composition] = get_bigchunk(obj.binary_bin, obj.param.bigchunk_rem_weight,...
                obj.param.bigchunk_windowsize, sleep_score.binwidth_sec);
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

