classdef analysis_radon < handle
    % ANALYSIS_RADON Wrapper class for executing and managing Radon analysis on an ROI.
    %   Extracts the 'radon' ROI from the provided dataset, executes analyze_radon,
    %   and provides methods to save the results.

    properties
        radon_result
        roi_image
        channel
    end

    methods
        function obj = analysis_radon(twophoton_processed, roilist, target_channel)
            % ANALYSIS_RADON Construct an instance of this class
            %   twophoton_processed: Struct containing ch1, ch2 data stacks
            %   roilist: roi_handle object containing 'radon' ROI
            %   target_channel: String (e.g., 'ch1' or 'ch2') indicating channel for Radon

            if nargin < 3
                target_channel = 'ch1'; % Default behavior
            end
            obj.channel = target_channel;

            % Extract the stack using the 'radon' ROI from Target Channel
            disp(['Extracting Radon ROI stack from ' obj.channel '...']);
            try
                obj.roi_image = roilist.applyvertices(twophoton_processed.(obj.channel), 'radon');
            catch ME
                error('Failed to extract Radon ROI. Ensure "radon" ROI exists. Error: %s', ME.message);
            end

            % Run the analysis function (Stats only)
            disp('Running Radon Analysis (Stats)...');
            obj.radon_result = analyze_radon(obj.roi_image);

            % Perform Event Selection
            disp('Selecting Events...');
            idx_struct = obj.select_events();

            % Perform Reconstruction (TIRS/IRTD)
            disp('Reconstructing Events...');
            obj.reconstruct_events(idx_struct);

            % Generate Median Projections
            obj.make_median_projection(twophoton_processed, roilist);
        end

        function idx_struct = select_events(obj)
            % SELECT_EVENTS Identify indices for different event types based on diameter variance.

            radon_result = obj.radon_result;

            % 1. Identify Highest Varying Angle & Smooth Trace
            [~, max_var_idx] = max(radon_result.var_diameter);
            raw_trace = radon_result.diameter(max_var_idx, :);
            if isa(raw_trace, 'gpuArray'), raw_trace = gather(raw_trace); end
            smooth_trace = movmedian(raw_trace, 5); % Smoothing
            med_val = median(smooth_trace);

            % 2. Define Indices for Events
            idx_struct = struct();
            % Thresholds based on 95th and 5th percentiles
            thresh_dilate = prctile(smooth_trace, 95);
            thresh_constrict = prctile(smooth_trace, 5);
            idx_struct.dilate = find(smooth_trace > thresh_dilate);
            idx_struct.constrict = find(smooth_trace < thresh_constrict);
            [~, sorted_diff] = sort(abs(smooth_trace - med_val));
            idx_struct.median = sorted_diff(1:min(200, end));

            % Max Event: Rising -> Peak -> Falling to Median
            [~, idx_peak] = max(smooth_trace);
            start_i = idx_peak;
            while start_i > 1 && smooth_trace(start_i) > med_val
                start_i = start_i - 1;
            end
            end_i = idx_peak;
            while end_i < length(smooth_trace) && smooth_trace(end_i) > med_val
                end_i = end_i + 1;
            end
            idx_struct.maxevent = start_i:end_i;
        end

        function reconstruct_events(obj, idx_struct)
            % RECONSTRUCT_EVENTS Generate TIRS and IRTD for identified events and store in obj.radon_result.events

            % Parameters from radon_result
            sz = obj.radon_result.radon_size;
            maxlocarray = obj.radon_result.idx_maxloc;
            upperboundary_idx = obj.radon_result.idx_uploc;
            bottomboundary_idx = obj.radon_result.idx_downloc;
            the_angles = 1:1:180; % Fixed angles

            % Initialize Output Struct Array
            event_types = {'dilate', 'constrict', 'median', 'maxevent'};
            events = struct('name', event_types, 'indices', [], 'tirs', [], 'irtd', []);

            for i = 1:length(event_types)
                ename = event_types{i};
                indices = idx_struct.(ename);

                events(i).indices = indices;

                if isempty(indices)
                    continue;
                end

                % Dimensions for subset
                n_frames = length(indices);
                sz_sub = [sz(1), sz(2), n_frames];

                % Reconstruct tirs and mask for this subset
                [sub_row_idx, ~, ~] = ndgrid(1:sz_sub(1), 1:sz_sub(2), 1:sz_sub(3));

                % Prepare boundaries for subset
                sub_max = reshape(maxlocarray(:, indices), [1, sz(2), n_frames]);
                sub_up  = reshape(upperboundary_idx(:, indices), [1, sz(2), n_frames]);
                sub_down = reshape(bottomboundary_idx(:, indices), [1, sz(2), n_frames]);

                % Expand to 3D
                max_3d = repmat(sub_max, [sz(1), 1, 1]);
                up_3d  = repmat(sub_up,  [sz(1), 1, 1]);
                down_3d = repmat(sub_down, [sz(1), 1, 1]);

                % Build TIRS (Visual Map)
                this_tirs = zeros(sz_sub);
                this_tirs(sub_row_idx == max_3d) = 2;
                this_tirs(sub_row_idx == up_3d)  = 1;
                this_tirs(sub_row_idx == down_3d) = 3;

                events(i).tirs = this_tirs;

                % Build Mask for Inverse Radon
                this_mask = false(sz_sub);
                this_mask(sub_row_idx >= up_3d) = 1;
                this_mask(sub_row_idx > down_3d) = 0;

                % Inverse Radon Transform
                recon_irtd = zeros([size(iradon(double(this_mask(:,:,1)), the_angles, 'linear', 'Hamming', 1, sz(1))), n_frames]);

                this_mask_gpu = gpuArray(this_mask);

                for f = 1:n_frames
                    recon_irtd(:,:,f) = iradon(this_mask_gpu(:,:,f), the_angles, 'linear', 'Hamming', 1, sz(1));
                end

                events(i).irtd = gather(recon_irtd); % Gather here to store in object
                events(i).tirs = gather(this_tirs);
            end

            obj.radon_result.events = events;
        end

        function make_median_projection(obj, twophoton_processed, roilist)
            % MAKE_MEDIAN_PROJECTION Computes median projections for events
            % using both channels from the original data.

            disp('Calculating Median Projections...');

            % Extract ROI stacks for both channels
            try
                stack_ch1 = roilist.applyvertices(twophoton_processed.ch1, 'radon');
            catch
                stack_ch1 = [];
            end
            try
                stack_ch2 = roilist.applyvertices(twophoton_processed.ch2, 'radon');
            catch
                stack_ch2 = [];
            end

            if isfield(obj.radon_result, 'events')
                for i = 1:length(obj.radon_result.events)
                    indices = obj.radon_result.events(i).indices;

                    if isempty(indices)
                        obj.radon_result.events(i).median_projected = [];
                        continue;
                    end

                    % Initialize 3D median projection (H, W, Ch)
                    % ch1 is index 1, ch2 is index 2
                    med_proj = [];

                    if ~isempty(stack_ch1)
                        img_ch1 = median(stack_ch1(:,:,indices), 3);
                        med_proj(:,:,1) = img_ch1;
                    end

                    if ~isempty(stack_ch2)
                        img_ch2 = median(stack_ch2(:,:,indices), 3);
                        med_proj(:,:,2) = img_ch2;
                    end

                    obj.radon_result.events(i).median_projected = med_proj;
                end
            end
            disp('Median Projection Complete.');
        end

        function save2disk(obj, save_dir)
            % SAVE2DISK Save the analysis results to a MAT file
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end

            radon_result = obj.radon_result; % Extract for saving variable name
            save_path = fullfile(save_dir, 'radon_result.mat');
            disp(['Saving Radon results to: ' save_path]);
            save(save_path, 'radon_result', '-v7.3');
        end
    end
end
