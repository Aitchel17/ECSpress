function final_binary_movement = get_binaryball(absball, binary_smooth,gap_tolerance)
%GET_BINARYBALL Summary of this function goes here
%   Detailed explanation goes here
% 3. Fill gaps where time between a falling edge and next rising edge is < 3 seconds

           % 2. Detect rising and falling edges of the smoothed binary movement
            diff_move = diff([0; binary_smooth(:); 0]);
            rise_idx = find(diff_move == 1);
            fall_idx = find(diff_move == -1);

            for i = 1:(length(rise_idx) - 1)
                gap_duration = rise_idx(i+1) - fall_idx(i);
                if gap_duration < gap_tolerance
                    binary_smooth(fall_idx(i):rise_idx(i+1)) = true;
                end
            end
            
            % 4. Recalculate edges of the continuous bouts after gap filling
            final_diff = diff([0; binary_smooth(:); 0]);
            final_rise_idx = find(final_diff == 1);
            final_fall_idx = find(final_diff == -1);
            
            % 5. Compensate for smoothing effect using absball > 0.02
            % For each filled bout, shrink or expand the start/end indices 
            % to perfectly match the raw absball > 0.02 boundaries.
            raw_movement = absball > 0.02;
            
            for i = 1:length(final_rise_idx)
                curr_rise = final_rise_idx(i);
                curr_fall = final_fall_idx(i);
                
                % --- Adjust Rising Edge (Onset) ---
                if ~raw_movement(curr_rise)
                    % SHRINK window forward until we hit raw movement
                    while curr_rise <= curr_fall && ~raw_movement(curr_rise)
                        curr_rise = curr_rise + 1;
                    end
                else
                    % EXPAND window backward as long as there is raw movement
                    while curr_rise > 1 && raw_movement(curr_rise - 1)
                        curr_rise = curr_rise - 1;
                    end
                end
                
                % --- Adjust Falling Edge (Offset) ---
                if ~raw_movement(curr_fall)
                    % SHRINK window backward until we hit raw movement
                    while curr_fall >= curr_rise && ~raw_movement(curr_fall)
                        curr_fall = curr_fall - 1;
                    end
                else
                    % EXPAND window forward as long as there is raw movement
                    while curr_fall < length(absball) && raw_movement(curr_fall + 1)
                        curr_fall = curr_fall + 1;
                    end
                end
                
                final_rise_idx(i) = curr_rise;
                final_fall_idx(i) = curr_fall;
            end
            
            % 6. Reconstruct the final compensated binary movement array
            final_binary_movement = false(size(absball));
            for i = 1:length(final_rise_idx)
                % Only add the bout if it contains actual raw movement 
                % (curr_rise <= curr_fall)
                if final_rise_idx(i) <= final_fall_idx(i)
                    final_binary_movement(final_rise_idx(i):final_fall_idx(i)) = true;
                end
            end
end

