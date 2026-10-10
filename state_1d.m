classdef state_1d < handle
    %STATE_1D  Centralize bout level analysis

    %   Caller: centralize_statesession, centralize_driftsession, stateintegration_test

    properties
        time_table         % 1x1 struct   table name -> N x 2 [start end] bouts in seconds, bout_sleep's or bout_puffmove's
        transition_window  % 1x1 float s  full width of a *_trans window, the bout function's param.transition_window
        t_axis      % Time axis for the specific data being analyzed
        state_idx
        param        % Name of this analysis instance
    end


    properties
        state_summary      % Struct to store summary statistics
        powerdensity       % Struct to store power density
        band_decomposition % Struct to store band decomposition analysis
        transition          % Struct to store transition analysis
        peak_trough        % struct to store peak and trough analysis from decomposed data
        eventlist          % 1 x N struct  the detected events, ALL contents share them
        event              % struct        one table per content name, rows = events
    end


    methods    % set-up and access
        function obj = state_1d(bout_table, window_s)
            % bout_table  1x1 struct   bout_sleep's or bout_puffmove's time_table, trimmed to this series' record
            % window_s    1x1 float s  the bout function's param.transition_window
            obj.time_table = bout_table;
            obj.transition_window = window_s;
        end

        function obj = get_state_indices(obj, t_axis,fs)
            obj.state_idx = bout_frames(t_axis, obj.time_table);   % every bout table on this series' frames
            obj.t_axis = t_axis; % copy time axis
            obj.param.fs = fs; % copy sampling rate
        end

        function obj = get_events(obj, eventlist, event_info)
        %GET_EVENTS  Hold the event list get_eventsummary builds its tables on; call it first.
        % IN   eventlist   1 x N struct   analysis_event's eventlist
        %      event_info  struct         its info; [] = do not record one
            arguments
                obj
                eventlist  (1,:) struct
                event_info = []
            end
            obj.eventlist = eventlist;
            if ~isempty(event_info)
                obj.param.event = event_info;
            end
        end

        function target_table = get_filtered_table(obj, prop_name, field_name, target_state)
            %GET_FILTERED_TABLE Access a specific table and filter by State
            %   tbl = obj.get_filtered_table('state_summary', 'bv_thickness', 'NREM')
            target_table = obj.(prop_name).(field_name);
            rows = target_table.state_name == string(target_state);
            target_table = target_table(rows, :);
        end

        function save2disk(obj,name,savepath)
            state_1d = obj;
            save(fullfile(savepath,[name,'.mat']),'state_1d')
        end
    end


    methods    % calculations: each fills one table per content name
        function obj = get_summary(obj, name, data1d)
            sidx_fnames = fieldnames(obj.state_idx);
            sidxfilt_logic = ~endsWith(sidx_fnames,'_trans');
            sidx_fnames = sidx_fnames(sidxfilt_logic);

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


        function obj = get_transitionsummary(obj, name, data1d)
            % GET_TRANSITIONSUMMARY Analyze transitions with center of bout as transition point
            % Splits each bout in half: before and after the transition center

            sidx_fnames = fieldnames(obj.state_idx);
            sidxfilt_logic = endsWith(sidx_fnames,'_trans');
            sidx_fnames = sidx_fnames(sidxfilt_logic);

            % 1. Pre-calculate total rows needed (one row per transition bout)
            total_bouts = 0;
            for i = 1:length(sidx_fnames)
                total_bouts = total_bouts + size(obj.state_idx.(sidx_fnames{i}), 1);
            end

            % 2. Pre-allocate arrays for metadata
            state_name = strings(total_bouts, 1);
            bout_idx = zeros(total_bouts, 1);
            bout_duration = zeros(total_bouts, 1);
            total_bout = zeros(total_bouts, 1);
            transition_point = zeros(total_bouts, 1);  % Center index

            % Pre-allocate for BEFORE transition period
            data = cell(total_bouts, 1);
            pre_mean = zeros(total_bouts, 1);
            pre_median = zeros(total_bouts, 1);
            pre_q1 = zeros(total_bouts, 1);
            pre_q3 = zeros(total_bouts, 1);
            pre_var = zeros(total_bouts, 1);

            % Pre-allocate for AFTER transition period
            post_mean = zeros(total_bouts, 1);
            post_median = zeros(total_bouts, 1);
            post_q1 = zeros(total_bouts, 1);
            post_q3 = zeros(total_bouts, 1);
            post_var = zeros(total_bouts, 1);

            % 3. Fill Data
            % half the transition window in samples, the same for every table so the
            % data column stacks
            onset = round(obj.transition_window * obj.param.fs / 2);

            row_counter = 1;
            for sidx = 1:length(sidx_fnames) % per state transition
                fname = sidx_fnames{sidx};
                indices = obj.state_idx.(fname); % Nx2 matrix
                nbouts = size(indices, 1);
                if ~isempty(indices)

                    for b = 1:nbouts % per bout
                        start_i = indices(b, 1);
                        end_i   = indices(b, 2);

                        % Calculate center point (transition point)
                        center_i = round((start_i + end_i) / 2);
                        duration = (end_i - start_i) / obj.param.fs;

                        % Store Metadata
                        state_name(row_counter) = fname;
                        bout_idx(row_counter) = b;
                        bout_duration(row_counter) = duration;
                        total_bout(row_counter) = nbouts;
                        transition_point(row_counter) = center_i;

                        % Split data at center
                        data(row_counter) = {data1d(center_i-onset:center_i+onset)};
                        pre_seg = data1d(center_i-onset:center_i);
                        post_seg = data1d(center_i:center_i+onset);

                        % BEFORE transition statistics
                        pre_mean(row_counter) = mean(pre_seg, 'omitnan');
                        pre_median(row_counter) = median(pre_seg, 'omitnan');
                        pre_var(row_counter) = var(pre_seg, 'omitnan');

                        qs_before = quantile(pre_seg, [0.25, 0.75]);
                        pre_q1(row_counter) = qs_before(1);
                        pre_q3(row_counter) = qs_before(2);

                        % AFTER transition statistics
                        post_mean(row_counter) = mean(post_seg, 'omitnan');
                        post_median(row_counter) = median(post_seg, 'omitnan');
                        post_var(row_counter) = var(post_seg, 'omitnan');

                        qs_after = quantile(post_seg, [0.25, 0.75]);
                        post_q1(row_counter) = qs_after(1);
                        post_q3(row_counter) = qs_after(2);

                        row_counter = row_counter + 1;
                    end
                end
            end

            % 4. Create Table with all columns
            obj.transition.(name) = table(state_name, bout_idx, bout_duration, total_bout, transition_point, ...
                data, pre_mean, pre_median, pre_q1, pre_q3, pre_var, ...
                post_mean, post_median, post_q1, post_q3, post_var);

        end

        function obj = get_powerdensity(obj, name, data1d)
            % a bout too short for the tapers is skipped, not shortened; bout_idx keeps its
            % original number so it still matches state_summary's
            min_samples = 5;        % dpss(N,2,3) needs N > 2*NW = 4
            n_skipped = 0;

            sidx_fnames = fieldnames(obj.state_idx);
            nsidx = length(sidx_fnames);
            state_name = [];
            spectral_density = [];
            spectral_frequency = [];
            spectral_lognorm = [];
            bout_duration = [];
            bout_idx = [];
            total_bout = [];

            for sidx = 1:nsidx
                sidx_fname = sidx_fnames{sidx};
                nbouts =  size(obj.state_idx.(sidx_fname),1);
                for bout = 1:nbouts
                    start_idx = obj.state_idx.(sidx_fname)(bout,1);
                    end_idx = obj.state_idx.(sidx_fname)(bout,2);
                    duration = (end_idx - start_idx)/obj.param.fs;
                    if (end_idx - start_idx + 1) < min_samples
                        n_skipped = n_skipped + 1;
                        continue
                    end
                    crop_data = data1d(start_idx:end_idx);
                    spec = get_spectrogram(obj.param.fs,crop_data);
                    state_name = [state_name; string(sidx_fname)];
                    bout_idx = [bout_idx ; bout];
                    spectral_density = [spectral_density ; {spec.S}];
                    spectral_frequency = [spectral_frequency ; {spec.F}];
                    spectral_lognorm = [spectral_lognorm ; {spec.log_norm_spectrum}];
                    bout_duration = [bout_duration ; duration];
                    total_bout = [total_bout ; nbouts];
                end
            end
            obj.powerdensity.(name) = table(state_name,bout_idx,bout_duration,total_bout,...
                spectral_density,spectral_frequency,spectral_lognorm);
            if n_skipped > 0
                fprintf('      %s : %d bout(s) under %d samples skipped\n', ...
                    name, n_skipped, min_samples);
            end

        end

        function obj = decompose_signal(obj, name, data1d)
            %DECOMPOSE_SIGNAL  The whole record in five bands, Butterworth order 3, zero phase:
            %   continuous < 0.005 Hz, isf 0.005-0.1, vlf 0.1-0.3, lf 0.3-1, hf > 1 Hz
            decomposed = struct();
            bparam.dc_cutoff = 0.005;
            bparam.isf_range = [0.005, 0.1];
            bparam.vlf_range = [0.1, 0.3];
            bparam.lf_range = [0.3, 1.0];
            bparam.hf_cutoff = 1.0;
            bparam.order = 3;

            % Nyquist Frequency
            fn = obj.param.fs / 2;
            % 0. Continuous (0 - 0.005 Hz) -> Lowpass
            f_cutoff = bparam.dc_cutoff;
            Wn = f_cutoff / fn;
            [b, a] = butter(bparam.order, Wn, 'low');
            decomposed.continuous = {filtfilt(b, a, data1d)};
            % 1. Infraslow (0.005 - 0.1 Hz) -> Bandpass
            f_range = bparam.isf_range;
            Wn = f_range / fn;
            [b, a] = butter(bparam.order, Wn, 'bandpass');
            decomposed.isf = {filtfilt(b, a, data1d)};
            % 2. VLF (0.1 - 0.3 Hz) -> Bandpass
            f_range = bparam.vlf_range;
            Wn = f_range / fn;
            [b, a] = butter(bparam.order, Wn, 'bandpass');
            decomposed.vlf = {filtfilt(b, a, data1d)};
            % 3. LF (0.3 - 1 Hz) -> Bandpass
            f_range = bparam.lf_range;
            Wn = f_range / fn;
            [b, a] = butter(bparam.order, Wn, 'bandpass');
            decomposed.lf = {filtfilt(b, a, data1d)};
            % 4. HF (> 1 Hz) -> Highpass
            f_cutoff = bparam.hf_cutoff;
            Wn = f_cutoff / fn;
            [b, a] = butter(bparam.order, Wn, 'high');
            decomposed.hf_residual = {filtfilt(b, a, data1d)};
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

                        if iscell(decomposed.(bn))
                            signal_full = decomposed.(bn){1};
                        else
                            signal_full = decomposed.(bn);
                        end

                        if start_i > numel(signal_full) || end_i > numel(signal_full)
                            crop_sig = [];
                        else
                            crop_sig = signal_full(start_i:end_i);
                        end

                        if ~isempty(crop_sig) && length(crop_sig) > 3
                            % 1. Find Peaks
                            [pks, locs] = findpeaks(crop_sig);
                            % 2. Find Troughs
                            [troughs, t_locs] = findpeaks(-crop_sig);
                            troughs = -troughs;

                            % peak-to-peak intervals (s)
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

                            % peak-to-trough amplitude: each peak against the first trough after it
                            amps = [];
                            for p_i = 1:numel(locs)
                                c_loc = locs(p_i);
                                c_val = pks(p_i);
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
            % the same four names the other nests use, so one flattener serves them all
            giant_table = table(state_name, bout_idx, bout_duration, total_bout, ...
                'VariableNames', {'state_name', 'bout_idx', 'bout_duration', 'total_bout'});

            % then six columns per band, named <band>_<statistic>
            for k = 1:length(bandnames)
                bn = bandnames{k};

                band_tbl = table(band_stats.(bn).pp_time_mean, ...
                    band_stats.(bn).pp_time_std, ...
                    band_stats.(bn).pp_time_raw, ...
                    band_stats.(bn).pt_amp_mean, ...
                    band_stats.(bn).pt_amp_std, ...
                    band_stats.(bn).pt_amp_raw, ...
                    'VariableNames', ...
                    {strcat(bn, '_pp_time_mean'), ...
                    strcat(bn, '_pp_time_std'), ...
                    strcat(bn, '_pp_time_raw'), ...
                    strcat(bn, '_pt_amp_mean'), ...
                    strcat(bn, '_pt_amp_std'), ...
                    strcat(bn, '_pt_amp_raw')});

                giant_table = [giant_table, band_tbl];
            end

            obj.peak_trough.(name) = giant_table;
        end

        function obj = get_eventsummary(obj, name, data1d, opt)
        %GET_EVENTSUMMARY  What one content did over each event, as a table; row = event, in eventlist order.
        % IN   name     char, the content's name, as get_summary uses it
        %      data1d   1 x T float px, that content's trace
        %      pre_s    float s, segment reach BEFORE the event start. Default 10
        %      post_s   float s, segment reach AFTER it. Default 30
        % OUT  obj.event.(name), a height x 10 table. The segment's time axis and
        %      the anchor sit in its Properties.UserData
            arguments
                obj
                name   (1,:) char
                data1d (1,:) double
                opt.pre_s  (1,1) double {mustBeNonnegative} = 10
                opt.post_s (1,1) double {mustBePositive}    = 30
            end
            if isempty(obj.eventlist)
                error('state_1d:noEvents', ...
                    'call get_events first; the tables are built on its list.');
            end
            n_row = numel(obj.eventlist);
            n_sample = numel(data1d);
            fs = obj.param.fs;
            % The segment window, in samples either side of the event's start
            n_pre  = round(opt.pre_s  * fs);
            n_post = round(opt.post_s * fs);
            n_seg  = n_pre + n_post + 1;
            seg_offset = -n_pre:n_post;

            state         = strings(n_row, 1);
            bout_idx      = zeros(n_row, 1);
            bout_duration = zeros(n_row, 1);
            total_bout    = zeros(n_row, 1);
            polarity      = strings(n_row, 1);
            from          = zeros(n_row, 1);
            to            = zeros(n_row, 1);
            rise_s        = zeros(n_row, 1);
            d_px     = nan(n_row, 1);
            range_px = nan(n_row, 1);
            sd_px    = nan(n_row, 1);
            mean_px  = nan(n_row, 1);
            seg_px   = nan(n_row, n_seg);

            % bouts per state, so total_bout means what it means in state_summary
            n_bout_of = struct();
            for f = string(fieldnames(obj.state_idx))'
                n_bout_of.(f) = size(obj.state_idx.(f), 1);
            end

            for k = 1:n_row
                row = obj.eventlist(k);
                state(k)         = string(row.state);
                polarity(k)      = string(row.pol);
                bout_idx(k)      = row.bout;
                bout_duration(k) = row.end_sec - row.start_sec;
                if isfield(n_bout_of, row.state)
                    total_bout(k) = n_bout_of.(row.state);
                else
                    total_bout(k) = NaN;
                end
                from(k)     = row.from;
                to(k)       = row.to;
                rise_s(k)   = row.rise_s;
                span        = row.from:row.to;
                span        = span(span >= 1 & span <= n_sample);
                if isempty(span)
                    continue
                end
                d_px(k)     = data1d(span(end)) - data1d(span(1));
                range_px(k) = max(data1d(span)) - min(data1d(span));
                sd_px(k)    = std(data1d(span));
                mean_px(k)  = mean(data1d(span));
                % the aligned segment; past either end of the record it stays NaN, not shifted
                idx  = row.from + seg_offset;
                have = idx >= 1 & idx <= n_sample;
                seg_px(k, have) = data1d(idx(have));
            end

            % "state", not "state_name": this column also holds transition labels
            tbl = table(state, bout_idx, bout_duration, total_bout, ...
                polarity, from, to, rise_s, ...
                d_px, range_px, sd_px, mean_px);
            tbl.seg_px = seg_px;
            tbl.Properties.VariableUnits = {'', '', 's', '', ...
                '', 'frame', 'frame', 's', ...
                'px', 'px', 'px', 'px', 'px'};
            tbl.Properties.UserData = struct( ...
                'seg_time_s', seg_offset / fs, ...   % 1 x n_seg, 0 = the anchor
                'seg_anchor', 'from', ...            % which endpoint index 0 is
                'pre_s',  opt.pre_s, ...
                'post_s', opt.post_s, ...
                'fs',     fs);
            obj.event.(name) = tbl;
        end
    end


end
