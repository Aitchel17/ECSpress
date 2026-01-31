function time_idx = statebin2frame(sleepscore_timebin,target_taxis)
    %STATEBIN2FRAME Summary of this function goes here
    %   Detailed explanation goes here
    num_awake = size(sleepscore_timebin, 1);
    time_idx = [];
    
    for i = 1:num_awake
        [~, start_loc] = min(abs(target_taxis - sleepscore_timebin(i,1)));
        [~, end_loc] = min(abs(target_taxis - sleepscore_timebin(i,2)));
        time_idx = [time_idx, start_loc:end_loc];
    end
end

