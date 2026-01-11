
function trig_data = analyze_trigaverage(signal_array, start_idx,startoffset,duration)
    % arguments
    %     signal_array
    %     start_idx
    %     start_offset
    %     end_offset
    % end

    trig_data = [];
    trig_data.list = [];
    for idx = 1:length(start_idx)
        segment = medfilt1(signal_array(start_idx(idx)-startoffset:start_idx(idx)+duration), 5);
        baseline = mean(segment(1:30), 'all');
        norm_segment = (segment / baseline)-1;
        trig_data.list = [trig_data.list, norm_segment];
    end
    trig_data.mean = mean(trig_data.list, 2);
    sem_trace = std(trig_data.list, 0, 2) / sqrt(size(trig_data.list, 2));
    tval = tinv(0.975, size(trig_data.list, 2) - 1);
    trig_data.ci = tval * sem_trace;
end


