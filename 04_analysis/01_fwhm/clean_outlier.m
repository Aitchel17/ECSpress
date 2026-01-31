function idx_dataset = clean_outlier(idx_dataset, overwrite)
% CLEAN_OUTLIER Cleans outliers from index dataset by comparing with median filtered data.
arguments
    idx_dataset
    overwrite = false
end

idx_names = fieldnames(idx_dataset);
idx_names = string(idx_names');

% Identify already processed fields (those starting with 'clean_')
processed_nameidx = startsWith(idx_names,'clean_');
processed_names = idx_names(processed_nameidx);
processed_names = erase(processed_names,'clean_');

% Identify original fields
original_names = idx_names(~processed_nameidx);

if overwrite
    unprocessed_names = original_names;
else
    % Skip already processed fields if overwrite is false
    match_idx = ismember(processed_names,original_names);
    unprocess_logic = true(size(original_names));
    unprocess_logic(match_idx) = false;
    unprocessed_names = original_names(unprocess_logic);
end

for idx_name = unprocessed_names
    if isempty(idx_name)
        disp('all processed')
    else
        disp(idx_name{:})
        idx_data = idx_dataset.(idx_name{:});

        % Median filtering for reference
        ref_median = medfilt1(idx_data,9);
        diff_data = abs(idx_data - ref_median);
        above_noise = diff_data > 8;

        % Edge case handling
        edge_thresholded = above_noise;
        edge_thresholded(1) = 0;
        edge_thresholded(end) = 0;

        transit_thr = diff(edge_thresholded);

        if sum(transit_thr) == 0
            rise_logic = transit_thr == 1;
            rise_logic = [rise_logic, false];
            fall_logic = transit_thr == -1;
            fall_logic = [false, fall_logic];
        else
            disp(['number of rising and falling is different for ' idx_name{:}])
            % Handle mismatch if necessary or just warn
        end

        % Interpolate data between rising and falling edges
        interp_data = (idx_data(rise_logic) + idx_data(fall_logic)) / 2;

        % Map interpolation back to indices
        run = [1, abs(transit_thr)];
        run = cumsum(run);
        unique_run = run(edge_thresholded==1);
        unique_run = unique(unique_run);
        [~,runidx] = ismember(run,unique_run);

        idxdata_interp_idx = runidx > 0;
        interp_idx = runidx(idxdata_interp_idx);

        corrected_data = idx_data;

        % Apply correction
        % Check if there are points to interpolate to avoid indexing errors
        if ~isempty(interp_idx) && ~isempty(interp_data)
            corrected_data(idxdata_interp_idx) = interp_data(interp_idx);
        end

        idx_dataset.(['clean_' idx_name{:}]) = corrected_data;
    end
end
end
