function stats = get_peak_analysis(signal, fs)
%GET_PEAK_ANALYSIS Calculate peak properties of a signal
%   Input:
%       signal: 1D vector
%       fs: sampling frequency
%   Output:
%       stats.peak2peak_time: array of time intervals between peaks (s)
%       stats.peak2trough_amp: array of amplitudes from Peak to Next Trough
%       stats.p2p_avg, stats.p2p_std
%       stats.p2t_avg, stats.p2t_std
%       stats.peak_idx, stats.peak_val
%       stats.trough_idx, stats.trough_val

stats = struct();
% Initialize all fields to empty/NaN
stats.peak2peak_time = [];
stats.peak2trough_amp = [];
stats.p2p_avg = NaN;
stats.p2p_std = NaN;
stats.p2t_avg = NaN;
stats.p2t_std = NaN;
stats.peak_idx = [];
stats.peak_val = [];
stats.trough_idx = [];
stats.trough_val = [];

if isempty(signal)
    return;
end

% 1. Find Peaks (Maxima)
[pks, locs] = findpeaks(signal);

% 2. Find Troughs (Minima)
[troughs, t_locs] = findpeaks(-signal);
troughs = -troughs; % Restore original values

% Store Raw Data
stats.peak_idx = locs;
stats.peak_val = pks;
stats.trough_idx = t_locs;
stats.trough_val = troughs;

% 3. Peak to Peak Time (seconds)
if numel(locs) > 1
    stats.peak2peak_time = diff(locs) / fs;
end

% 4. Peak to Trough Amplitude
amps = [];
for i = 1:numel(locs)
    current_peak_loc = locs(i);
    current_peak_val = pks(i);

    % Find first trough occurring after this peak
    next_trough_idx_in_list = find(t_locs > current_peak_loc, 1, 'first');

    if ~isempty(next_trough_idx_in_list)
        trough_val = troughs(next_trough_idx_in_list);
        amps(end+1) = current_peak_val - trough_val;
    end
end
stats.peak2trough_amp = amps;

% 5. Statistics
if ~isempty(stats.peak2peak_time)
    stats.p2p_avg = mean(stats.peak2peak_time);
    stats.p2p_std = std(stats.peak2peak_time);
end

if ~isempty(stats.peak2trough_amp)
    stats.p2t_avg = mean(stats.peak2trough_amp);
    stats.p2t_std = std(stats.peak2trough_amp);
end

end
