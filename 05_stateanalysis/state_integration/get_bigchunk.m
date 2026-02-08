function [timetable, composition] = get_bigchunk(binary_bin, weights, window_thr, binwidth)
arguments
    binary_bin struct
    weights struct
    window_thr (1,1) double
    binwidth (1,1) double
end
%GET_BIGCHUNK Find optimally stable time windows for sleep states.
%   Pure function version.
%   Inputs:
%       binary_bin: struct with fields 'nrem', 'awake', 'rem', 'drowsy' (logical/double vectors)
%       weights:    struct with fields 'NREM', 'Awake', ... containing weights (.w_nrem, .w_awake...)
%       window_thr: window duration in seconds
%       binwidth:   duration of each bin in seconds

% 1. Setup Grid
% Assuming all binary vectors are same length
fields = fieldnames(binary_bin);
num_bins = length(binary_bin.(fields{1}));

% 2. Extract State Vectors
% Ensure they are double for calculation
vec_awake = double(binary_bin.awake)* weights.w_awake;
vec_nrem  = double(binary_bin.nrem)* weights.w_nrem;
vec_rem   = double(binary_bin.rem)* weights.w_rem;
vec_drowsy= double(binary_bin.drowsy)* weights.w_drowsy;

% 3. Define Analysis Window
window_bins = round(window_thr / binwidth);
kernel = ones(window_bins, 1);

% 4. Process Each Target State defined in weights
% Compute Score using passed weights
% Expecting w_cfg to have fields: w_nrem, w_drowsy, w_awake, w_rem
score_vec = vec_nrem + vec_drowsy + vec_awake + vec_rem;

% Convolve
score_conv = conv(score_vec, kernel, 'valid');

% Find Max
[~, max_idx] = max(score_conv);

% Calculate Times
start_bin_idx = max_idx;
start_time = (start_bin_idx - 1) * binwidth;
end_time = start_time + window_thr;

% Calculate Actual Composition
slice_s = start_bin_idx;
slice_e = min(start_bin_idx + window_bins - 1, num_bins);

timetable = [start_time, end_time];
composition.Awake = sum(vec_awake(slice_s:slice_e)) / window_bins;
composition.NREM = sum(vec_nrem(slice_s:slice_e)) / window_bins;
composition.REM = sum(vec_rem(slice_s:slice_e)) / window_bins;
composition.Drowsy = sum(vec_drowsy(slice_s:slice_e)) / window_bins;

end
