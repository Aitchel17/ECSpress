function result = get_bigchunk_analysis(analysis_obj, sleep_score, state_name, fs, window_sec)
arguments
    analysis_obj
    sleep_score struct
    state_name string
    fs (1,1) double
    window_sec (1,1) double = 300
end
%GET_BIGCHUNK_ANALYSIS Identify 'big chunk' of state and run spectral analysis on it.
%   1. Finds best [window_sec] chunk for [state_name] using get_bigchunk.
%   2. Runs get_all_spectra on that specific time window.

result = struct();

% 1. Find the big chunk
bc_all = get_bigchunk(sleep_score, window_sec);
if isfield(bc_all, state_name)
    bc_info = bc_all.(state_name);
else
    error('State "%s" not found in get_bigchunk results.', state_name);
end

result.info = bc_info;

if isnan(bc_info.start_time)
    warning('No sufficient chunks found for state: %s', state_name);
    return;
end

% 2. Convert to indices
if ~isprop(analysis_obj, 't_axis') && ~isfield(analysis_obj, 't_axis')
    warning('Analysis object missing t_axis.');
    return;
end
t_axis = analysis_obj.t_axis;

idx_s = find(t_axis >= bc_info.start_time, 1, 'first');
idx_e = find(t_axis <= bc_info.end_time, 1, 'last');

if isempty(idx_s) || isempty(idx_e) || idx_s > idx_e
    warning('Chunk indices could not be mapped to t_axis.');
    return;
end

bouts = {idx_s:idx_e};

% 3. Run Spectrogram/Peak Analysis
% get_all_spectra expects a cell array of index vectors
result.spectral_analysis = get_all_spectra(bouts, analysis_obj, fs);
end
