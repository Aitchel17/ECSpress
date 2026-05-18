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
            plot(taxis,primary_analog.ball.resampled_ball)
            %%
            absball = abs(primary_analog.ball.resampled_ball);
            %% 1 second smoothing 
            kernelWidth = 1; % window (sec)
            smoothingKernel = gausswin(kernelWidth*fs)/sum(gausswin(kernelWidth*fs));
            ballSmooth = conv(absball,smoothingKernel,'same');    
            %% Awake Integration: Continuous Movement Detection
            % 1. Detect binary movement from the SMOOTHED ball data first
            smooth_threshold = 0.02; 
            binary_smooth = abs(ballSmooth) > smooth_threshold;
            gap_tolerance = 3 * fs;
            binary_move = get_binaryball(absball,binary_smooth,gap_tolerance);
            %%
            move_table = logic2timetable(taxis,binary_move,5);
            obj.time_table.awake_movement = [taxis(rise_idx(:))', taxis(fall_idx(:) - 1)'];
            obj.time_table.nr_trans = get_transition(obj.param.transition_window,obj.time_table.roughnrem, obj.time_table.rem);
            obj.time_table.na_trans = get_transition(obj.param.transition_window,obj.time_table.roughnrem, obj.time_table.roughawake);
            obj.time_table.ra_trans = get_transition(obj.param.transition_window,obj.time_table.rem, obj.time_table.roughawake);
            obj.time_table.an_trans = get_transition(obj.param.transition_window,obj.time_table.roughawake, obj.time_table.roughnrem);       

            %%
            figure()
            plot(taxis,binary_move)

            %%
            cla
            plot(taxis,ballSmooth)
            hold on
            plot(taxis,absball,'Color','k')
            %%
            plot(taxis, binary_move, 'r', 'LineWidth', 1.5)
            hold on
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

    end
end

