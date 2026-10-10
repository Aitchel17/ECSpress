function airtable = make_airtable(air_on, air_fs, gap_s)
%MAKE_AIRTABLE  The puffs on the air channel: pulses whose onsets lie within gap_s are one puff.
%   Caller: bout_puffmove
%   02_othersignal's analysis_analog.get_airtable, on the sample clock (k - 1) / air_fs.
% IN   air_on    1 x T bool         the air channel, on per sample
%      air_fs    1 x 1 float Hz
%      gap_s     1 x 1 float s      bout_puffmove's param.puff_gap_s
% OUT  airtable  N x 5 table        StartTime, EndTime, Duration (s), FrequencyHz, DutyCycle;
%                                   the last two the mean over the pulses that have a next
%                                   one in the same puff, NaN for a single pulse
    edge = diff(air_on(:));
    rise_idx = find(edge == 1);    % the last off sample before each pulse
    fall_idx = find(edge == -1);   % the last on sample of each pulse
    rise_s = (rise_idx - 1) / air_fs;
    fall_s = fall_idx / air_fs;    % the first off sample after it
    onset_gap = diff(rise_s);
    is_first = [true; onset_gap > gap_s];
    is_last = [onset_gap > gap_s; true];
    start_s = rise_s(is_first);
    end_s = fall_s(is_last);
    % pulse k against pulse k + 1, kept where both are in one puff
    pulse_duty = (fall_idx(1:end-1) - rise_idx(1:end-1)) ./ diff(rise_idx);
    has_next = ~is_last(1:end-1);
    puff_id = cumsum(is_first);
    pulse_puff = puff_id(has_next);
    n_puff = numel(start_s);
    mean_gap = accumarray(pulse_puff, onset_gap(has_next), [n_puff 1], @mean, NaN);
    mean_duty = accumarray(pulse_puff, pulse_duty(has_next), [n_puff 1], @mean, NaN);
    airtable = table(start_s, end_s, end_s - start_s, 1 ./ mean_gap, mean_duty, ...
        'VariableNames', {'StartTime', 'EndTime', 'Duration', 'FrequencyHz', 'DutyCycle'});
end
