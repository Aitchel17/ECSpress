function stats = get_peak_analysis(signal, fs)
%GET_PEAK_ANALYSIS Calculate peak properties of a signal
%   Input:
%       signal: 1D vector
%       fs: sampling frequency
%   Output:
%       stats.peak2peak_time: array of time intervals between peaks (s)
%       stats.peak2trough_amp: array of amplitudes from Peak to Next Trough

stats = struct('peak2peak_time', [], 'peak2trough_amp', []);

if isempty(signal)
    return;
end

% 1. Find Peaks (Maxima)
% MinPeakProminence could be added if noisy, but for decomposed signals (filtered), usually clean.
[pks, locs] = findpeaks(signal);

% 2. Find Troughs (Minima)
[troughs, t_locs] = findpeaks(-signal);
troughs = -troughs;

% 3. Peak to Peak Time (seconds)
if numel(locs) > 1
    stats.peak2peak_time = diff(locs) / fs;
else
    stats.peak2peak_time = [];
end

% 4. Peak to Trough Amplitude
% Algorithm: For each peak, find the nearest subsequent trough.
% If a trough appears before the next peak, calculate amp.

amps = [];
for i = 1:numel(locs)
    current_peak_loc = locs(i);
    current_peak_val = pks(i);

    % Find first trough occurring after this peak
    next_trough_idx_in_list = find(t_locs > current_peak_loc, 1, 'first');

    if ~isempty(next_trough_idx_in_list)
        trough_loc = t_locs(next_trough_idx_in_list);

        % Check if this trough is 'valid' (i.e., before the next peak?)
        % Usually strictly alternating, but finding the *immediate* next is correct
        % for "Peak to Trough" drop.

        % Optional: Ensure we don't skip a peak (e.g. if Peak1 -> Peak2 -> Trough1,
        % then Peak1->Trough1 covers Peak2).
        % But findpeaks usually guarantees alternating check if prominence is handled,
        % but distinct local maxima can happen.
        % We will strictly take the nearest subsequent trough.

        trough_val = troughs(next_trough_idx_in_list);
        amps(end+1) = current_peak_val - trough_val;
    end
end

stats.peak2trough_amp = amps;

end
