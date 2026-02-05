classdef StateAnalysis_1D < handle
    %STATEANALYSIS_ABSTRACT Base class for sleep state analysis.
    %   Coordinated by sleep_integration object.

    properties
        sleep_obj   % Handle to sleep_integration object (Temporal Source of Truth)
        t_axis      % Time axis for the specific data being analyzed
        state_idx
        data_struct        % The raw data to analyze (Dimension varies by subclass)
        power_analysis     % Struct to store analysis results
        summary_analysis   % Struct to store analysis results
        param        % Name of this analysis instance
    end

    methods
        function obj = StateAnalysis_1D(sleep_obj)
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
                nbouts = size(obj.state_idx.(sidx_fname),1);
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
            obj.power_analysis.(name) = table(state_name,bout_idx,bout_duration,total_bout,...
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
            val_mean = zeros(total_bouts, 1);
            val_median = zeros(total_bouts, 1);
            val_q1 = zeros(total_bouts, 1);
            val_q3 = zeros(total_bouts, 1);
            val_var = zeros(total_bouts, 1);

            % 3. Fill Data
            row_counter = 1;
            for sidx = 1:length(sidx_fnames)
                fname = sidx_fnames{sidx};
                indices = obj.state_idx.(fname); % Nx2 matrix
                nbouts = size(indices, 1);

                for b = 1:nbouts
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

                    val_mean(row_counter) = mean(crop_data, 'omitnan');
                    val_median(row_counter) = median(crop_data, 'omitnan');
                    val_var(row_counter) = var(crop_data, 'omitnan');

                    % Quartiles
                    qs = quantile(crop_data, [0.25, 0.75]);
                    val_q1(row_counter) = qs(1);
                    val_q3(row_counter) = qs(2);

                    row_counter = row_counter + 1;
                end
            end

            % 4. Create Table
            obj.summary_analysis.(name) = table(state_name, bout_idx, ...
                bout_duration, total_bout, raw_data, ...
                val_mean, val_median, val_q1, val_q3, val_var);
        end


    end


end
