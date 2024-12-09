function [interp_drifttable, zstack_corr] = pre_applymotion(zstack, drift_table)
    % Apply motion correction using drift table and preserve common area
    fprintf('Apply drift correction\n');
    
    % Validate inputs
    [height, width, num_frames] = size(zstack);
    assert(size(drift_table, 1) >= 4, 'drift_table must have at least 4 rows');
    
    % Extract and interpolate shifts
    original_frames = 1:size(drift_table, 2);
    interp_frames = linspace(1, size(drift_table, 2), num_frames);
    row_shifts_interp = interp1(original_frames, drift_table(3, :), interp_frames, 'linear');
    col_shifts_interp = interp1(original_frames, drift_table(4, :), interp_frames, 'linear');
    
    % Determine cropping boundaries for common area
    min_row_shift = floor(min(row_shifts_interp));
    max_row_shift = ceil(max(row_shifts_interp));
    min_col_shift = floor(min(col_shifts_interp));
    max_col_shift = ceil(max(col_shifts_interp));
    
    crop_top = max(1, 1 - min_row_shift);
    crop_bottom = min(height, height - max_row_shift);
    crop_left = max(1, 1 - min_col_shift);
    crop_right = min(width, width - max_col_shift);
    
    % Initialize corrected z-stack with the common area size
    cropped_height = crop_bottom - crop_top + 1;
    cropped_width = crop_right - crop_left + 1;
    zstack_corr = zeros(cropped_height, cropped_width, num_frames, 'like', zstack);
    
    % Process frames
    tic;
    for i = 1:num_frames
        disp(['Motion correction: ', num2str(i), '/', num2str(num_frames)]);
        frame = zstack(:, :, i);
        row_shift = round(row_shifts_interp(i));
        col_shift = round(col_shifts_interp(i));
        
        % Shift the frame
        corrected_frame = imtranslate(frame, [col_shift, row_shift], ...
                                      'OutputView', 'same', 'FillValues', 0);
        
        % Crop to the common area
        zstack_corr(:, :, i) = corrected_frame(crop_top:crop_bottom, crop_left:crop_right);
    end
    toc;
    fprintf('Pixel shift correction completed.\n');
    
    % Return results
    interp_drifttable = [row_shifts_interp; col_shifts_interp];
end
