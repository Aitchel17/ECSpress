classdef state_linefwhm < handle
    %STATEANALYSIS_ABSTRACT Base class for sleep state analysis.
    %   Coordinated by sleep_integration object.

    properties
        sleep_obj   % Handle to sleep_integration object (Temporal Source of Truth)
        t_axis      % Time axis for the specific data being analyzed
        state_idx
        data_struct        % The raw data to analyze (Dimension varies by subclass)
        param        % Name of this analysis instance
    end


    properties
        state_summary      % Struct to store summary statistics
        powerdensity       % Struct to store power density
        band_decomposition % Struct to store band decomposition analysis
        transition          % Struct to store transition analysis
        peak_trough        % struct to store peak and trough analysis from decomposed data
    end

    methods
        function obj = state_linefwhm(sleep_obj)
            % Constructor
            obj.sleep_obj = sleep_obj;
        end

        function obj = get_state_indices(obj, t_axis,fs)
            obj.state_idx =  obj.sleep_obj.add_taxis(t_axis);
            obj.t_axis = t_axis;
            obj.param.fs = fs;
        end

        function obj = get_powerdensity(obj, name, data1d)
            sidx_fnames = fieldnames(obj.state_idx);
            nsidx = length(sidx_fnames);
            %%
            state_name = [];
            spectral_density = [];
            spectral_frequency = [];
            spectral_lognorm = [];
            bout_duration = [];
            bout_idx = [];
            total_bout = [];

            for sidx = 1:nsidx % loop here
                sidx_fname = sidx_fnames{sidx};
                nbouts =  size(obj.state_idx.(sidx_fname),1);
                for bout = 1:nbouts
                    start_idx = obj.state_idx.(sidx_fname)(bout,1);
                    end_idx = obj.state_idx.(sidx_fname)(bout,2);
                    duration = (end_idx - start_idx)/obj.param.fs;
                    crop_data = data1d(start_idx:end_idx);
                    spec = get_spectrogram(obj.param.fs,crop_data);
                    % specify
                    state_name = [state_name; string(sidx_fname)];
                    bout_idx = [bout_idx ; bout];
                    % spec contents
                    spectral_density = [spectral_density ; {spec.S}];
                    spectral_frequency = [spectral_frequency ; {spec.F}];
                    spectral_lognorm = [spectral_lognorm ; {spec.log_norm_spectrum}];
                    bout_duration = [bout_duration ; duration];
                    total_bout = [total_bout ; nbouts];
                end
            end
            obj.powerdensity.(name) = table(state_name,bout_idx,bout_duration,total_bout,...
                spectral_density,spectral_frequency,spectral_lognorm);

        end

        function obj = get_summary(obj, name, data1d)
            sidx_fnames = fieldnames(obj.state_idx);

            % 1. Pre-calculate total rows needed to pre-allocate
            total_bouts = 0;
            for i = 1:length(sidx_fnames)
                total_bouts = total_bouts + size(obj.state_idx.(sidx_fnames{i}), 1);
            end

            % 2. Pre-allocate arrays
            state_name = strings(total_bouts, 1);
            bout_idx = zeros(total_bouts, 1);
            bout_duration = zeros(total_bouts, 1);
            total_bout = zeros(total_bouts, 1);
            raw_data = cell(total_bouts, 1);
            raw_mean = zeros(total_bouts, 1);
            raw_median = zeros(total_bouts, 1);
            raw_q1 = zeros(total_bouts, 1);
            raw_q3 = zeros(total_bouts, 1);
            raw_var = zeros(total_bouts, 1);

            % 3. Fill Data
            row_counter = 1;
            for sidx = 1:length(sidx_fnames) % per state
                fname = sidx_fnames{sidx};
                indices = obj.state_idx.(fname); % Nx2 matrix
                nbouts = size(indices, 1);

                for b = 1:nbouts % bout
                    start_i = indices(b, 1);
                    end_i   = indices(b, 2);

                    crop_data = data1d(start_i:end_i);
                    duration = (end_i - start_i) / obj.param.fs;

                    % Store Metadata
                    state_name(row_counter) = fname;
                    bout_idx(row_counter) = b;
                    bout_duration(row_counter) = duration;
                    total_bout(row_counter) = nbouts;

                    % Store Data & Stats
                    raw_data{row_counter} = crop_data;

                    raw_mean(row_counter) = mean(crop_data, 'omitnan');
                    raw_median(row_counter) = median(crop_data, 'omitnan');
                    raw_var(row_counter) = var(crop_data, 'omitnan');

                    % Quartiles
                    qs = quantile(crop_data, [0.25, 0.75]);
                    raw_q1(row_counter) = qs(1);
                    raw_q3(row_counter) = qs(2);

                    row_counter = row_counter + 1;
                end
            end

            % 4. Create Table
            obj.state_summary.(name) = table(state_name, bout_idx, ...
                bout_duration, total_bout, raw_data, ...
                raw_mean, raw_median, raw_q1, raw_q3, raw_var);
        end

        function obj = decompose_signal(obj, name, data1d)
            % Get_band decomposition do three things which involves
            % 1. band decomposition (continous = cnt , Infraslow = isf, Verylow = vlf, Low = lf)
            % 2. peak2peak interval calculation and peak2trough amplitude calculation
            % 3. update summary_stat table
            % 4.
            decomposed = struct();
            %% Paramter setup
            bparam.dc_cutoff = 0.005;
            bparam.isf_range = [0.005, 0.1];
            bparam.vlf_range = [0.1, 0.3];
            bparam.lf_range = [0.3, 1.0];
            bparam.hf_cutoff = 1.0;
            bparam.order = 3;

            % Nyquist Frequency
            fn = obj.param.fs / 2;
            %% 0. Continuous (0 - 0.005 Hz) -> Lowpass
            f_cutoff = bparam.dc_cutoff;
            Wn = f_cutoff / fn;
            [b, a] = butter(bparam.order, Wn, 'low');
            decomposed.continuous = filtfilt(b, a, data1d);
            %% 1. Infraslow (0.005 - 0.1 Hz) -> Bandpass
            f_range = bparam.isf_range;
            Wn = f_range / fn;
            [b, a] = butter(bparam.order, Wn, 'bandpass');
            decomposed.isf = filtfilt(b, a, data1d);
            %% 2. VLF (0.1 - 0.3 Hz) -> Bandpass
            f_range = bparam.vlf_range;
            Wn = f_range / fn;
            [b, a] = butter(bparam.order, Wn, 'bandpass');
            decomposed.vlf = filtfilt(b, a, data1d);
            %% 3. LF (0.3 - 1 Hz) -> Bandpass
            f_range = bparam.lf_range;
            Wn = f_range / fn;
            [b, a] = butter(bparam.order, Wn, 'bandpass');
            decomposed.lf = filtfilt(b, a, data1d);
            %% 4. HF (1 - 1.5 Hz) -> Bandpass
            f_cutoff = bparam.hf_cutoff; % Just close to Nyquist
            Wn = f_cutoff / fn;
            [b, a] = butter(bparam.order, Wn, 'high');
            decomposed.hf_residual = filtfilt(b, a, data1d);
            obj.param.butter = bparam;
            obj.band_decomposition.(name) = decomposed;
        end

        function obj = get_pppt_decomposition(obj, name)
            % GET_PPPT_DECOMPOSITION Calculate peak-to-peak and peak-to-trough stats for decomposed signals

            if ~isfield(obj.band_decomposition, name)
                warning('Decomposition data for %s not found. Run get_band_decomposition first.', name);
                return;
            end

            decomposed = obj.band_decomposition.(name);
            bandnames = fieldnames(decomposed); % e.g. continuous, isf, vlf, lf, hf
            sidx_fnames = fieldnames(obj.state_idx);

            % 1. Pre-calculate total rows needed to pre-allocate
            total_bouts = 0;
            for i = 1:length(sidx_fnames)
                total_bouts = total_bouts + size(obj.state_idx.(sidx_fnames{i}), 1);
            end

            % 2. Pre-allocate arrays
            state_name = strings(total_bouts, 1);
            bout_idx = zeros(total_bouts, 1);
            bout_duration = zeros(total_bouts, 1);
            total_bout = zeros(total_bouts, 1);

            % Initialize struct for band data (will convert to table later)
            band_stats = struct();
            for k = 1:length(bandnames)
                bn = bandnames{k};
                band_stats.(bn).pp_time_mean = zeros(total_bouts, 1);
                band_stats.(bn).pp_time_std = zeros(total_bouts, 1);
                band_stats.(bn).pt_amp_mean = zeros(total_bouts, 1);
                band_stats.(bn).pt_amp_std = zeros(total_bouts, 1);
                band_stats.(bn).pp_time_raw = cell(total_bouts, 1);
                band_stats.(bn).pt_amp_raw = cell(total_bouts, 1);
            end

            % 3. Fill Data
            row_counter = 1;
            for sidx = 1:length(sidx_fnames) % 1. per state
                fname = sidx_fnames{sidx};
                indices = obj.state_idx.(fname); % Nx2 matrix
                nbouts = size(indices, 1);

                for b = 1:nbouts  % 2. Per Bout
                    start_i = indices(b, 1);
                    end_i   = indices(b, 2);
                    duration = (end_i - start_i) / obj.param.fs;

                    % Store Metadata
                    state_name(row_counter) = fname;
                    bout_idx(row_counter) = b;
                    bout_duration(row_counter) = duration;
                    total_bout(row_counter) = nbouts;

                    % Process each band
                    for k = 1:length(bandnames) % 3. per band
                        bn = bandnames{k};
                        signal_full = decomposed.(bn);

                        % Slicing check
                        if start_i > numel(signal_full) || end_i > numel(signal_full)
                            crop_sig = [];
                        else
                            crop_sig = signal_full(start_i:end_i);
                        end

                        if ~isempty(crop_sig) && length(crop_sig) > 3
                            % Calculate PP/PT stats
                            % Use local function or internal logic to keep it contained?
                            % Doing inline logic for P-P and P-T

                            % 1. Find Peaks
                            [pks, locs] = findpeaks(crop_sig);
                            % 2. Find Troughs
                            [troughs, t_locs] = findpeaks(-crop_sig);
                            troughs = -troughs;

                            % PP Time
                            if numel(locs) > 1
                                pp_times = diff(locs) / obj.param.fs;
                                band_stats.(bn).pp_time_mean(row_counter) = mean(pp_times, 'omitnan');
                                band_stats.(bn).pp_time_std(row_counter) = std(pp_times, 'omitnan');
                                band_stats.(bn).pp_time_raw{row_counter} = pp_times;
                            else
                                band_stats.(bn).pp_time_mean(row_counter) = NaN;
                                band_stats.(bn).pp_time_std(row_counter) = NaN;
                                band_stats.(bn).pp_time_raw{row_counter} = [];
                            end

                            % PT Amp
                            % Naive pairing: for each peak, find nearest subsequent trough?
                            % Or just general stats?
                            % "Peak to Trough Amplitude" usually implies Peak - Next Trough.
                            amps = [];
                            for p_i = 1:numel(locs)
                                c_loc = locs(p_i);
                                c_val = pks(p_i);
                                % Find first trough after this peak
                                next_t_idx = find(t_locs > c_loc, 1, 'first');
                                if ~isempty(next_t_idx)
                                    t_val = troughs(next_t_idx);
                                    amps(end+1) = c_val - t_val;
                                end
                            end

                            if ~isempty(amps)
                                band_stats.(bn).pt_amp_mean(row_counter) = mean(amps, 'omitnan');
                                band_stats.(bn).pt_amp_std(row_counter) = std(amps, 'omitnan');
                                band_stats.(bn).pt_amp_raw{row_counter} = amps;
                            else
                                band_stats.(bn).pt_amp_mean(row_counter) = NaN;
                                band_stats.(bn).pt_amp_std(row_counter) = NaN;
                                band_stats.(bn).pt_amp_raw{row_counter} = [];
                            end

                        else
                            % Empty or too short
                            band_stats.(bn).pp_time_mean(row_counter) = NaN;
                            band_stats.(bn).pp_time_std(row_counter) = NaN;
                            band_stats.(bn).pt_amp_mean(row_counter) = NaN;
                            band_stats.(bn).pt_amp_std(row_counter) = NaN;
                        end
                    end

                    row_counter = row_counter + 1;
                end
            end

            % 4. Create one giant table for all bands
            % Start with metadata columns
            giant_table = table(state_name, bout_idx, bout_duration, total_bout, ...
                'VariableNames', {'State', 'BoutIdx', 'Duration', 'TotalBouts'});

            % Iterate through bands and append standardized columns
            for k = 1:length(bandnames)
                bn = bandnames{k};
                % Create valid variable names (e.g. 'vlf_pp_time_mean')
                % Assuming 'bn' is a valid field name like 'vlf', 'continuous'

                band_tbl = table(band_stats.(bn).pp_time_mean, ...
                    band_stats.(bn).pp_time_std, ...
                    band_stats.(bn).pp_time_raw, ...
                    band_stats.(bn).pt_amp_mean, ...
                    band_stats.(bn).pt_amp_std, ...
                    band_stats.(bn).pt_amp_raw, ...
                    'VariableNames', ...
                    {strcat(bn, '_Mean_PP_Time'), ...
                    strcat(bn, '_Std_PP_Time'), ...
                    strcat(bn, '_Raw_PP_Time'), ...
                    strcat(bn, '_Mean_PT_Amp'), ...
                    strcat(bn, '_Std_PT_Amp'), ...
                    strcat(bn, '_Raw_PT_Amp')});

                giant_table = [giant_table, band_tbl];
            end

            % Save as requested: obj.band_decomposition.(name+'pppt')
            % Note: 'name' is a char/string.
            obj.peak_trough.(name) = giant_table;
        end

        function target_table = get_filtered_table(obj, prop_name, field_name, target_state)
            %GET_FILTERED_TABLE Access a specific table and filter by State
            %   tbl = obj.get_filtered_table('state_summary', 'bv_thickness', 'NREM')
            target_table = obj.(prop_name).(field_name);
            rows = target_table.state_name == string(target_state);
            target_table = target_table(rows, :);
        end

       function save2disk(obj,name,savepath)
            state_linefwhm = obj;
            save(fullfile(savepath,[name,'.mat']),'state_linefwhm')
        end


    end


end
