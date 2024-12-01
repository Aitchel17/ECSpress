function [interp_drifttable, zstack_corr] = pre_applymotion(zstack, drift_table)
    % Apply motion correction using drift table and interpolate frame shifts
    fprintf('Apply drift correction');
    % Validate inputs
    [height, width, num_frames] = size(zstack);
    assert(size(drift_table, 1) >= 4, 'drift_table must have at least 4 rows');
    
    % Extract and interpolate shifts
    original_frames = 1:size(drift_table, 2);
    interp_frames = linspace(1, size(drift_table, 2), num_frames);
    row_shifts_interp = interp1(original_frames, drift_table(3, :), interp_frames, 'linear');
    col_shifts_interp = interp1(original_frames, drift_table(4, :), interp_frames, 'linear');
    
    % Prepare padded dimensions
    max_row_shift = ceil(max(abs(row_shifts_interp)));
    max_col_shift = ceil(max(abs(col_shifts_interp)));
    padded_height = height + 2 * max_row_shift;
    padded_width = width + 2 * max_col_shift;
    zstack_corr = zeros(padded_height, padded_width, num_frames, 'like', zstack);
    
    % Process frames in parallel

    tic;
    for i = 1:num_frames
        disp(['motion correction' num2str(i),'/',num2str(num_frames)])
        frame = zstack(:, :, i);
        row_shift = round(row_shifts_interp(i));
        col_shift = round(col_shifts_interp(i));
        
        % Pad and shift frame
        padded_frame = padarray(frame, [max_row_shift, max_col_shift], 0, 'both');
        corrected_frame = imtranslate(padded_frame, [col_shift, row_shift], ...
                                       'OutputView', 'same', 'FillValues', 0);
        zstack_corr(:, :, i) = corrected_frame;
    end
    toc;
    fprintf('Pixel shift correction completed.\n');
    
    % Return results
    interp_drifttable = [row_shifts_interp; col_shifts_interp];
end
