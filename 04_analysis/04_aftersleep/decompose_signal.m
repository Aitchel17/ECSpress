function decomposed = decompose_signal(signal, fs)
%DECOMPOSE_SIGNAL Decompose signal into frequency bands using Butterworth filters
%   Bands:
%       Continuous (0-0.005 Hz)
%       Infraslow (0.005 - 0.1 Hz)
%       VLF (0.1-0.3 Hz)
%       LF (0.3-1 Hz)

%   Filter: Butterworth order 3, zero-phase (filtfilt)

decomposed = struct();

% Handle dimensions: ensure time is first dimension for filtfilt or handle usually
% standard: filtfilt operates along first non-singleton dim.
% If strict 1D vector, ensure column for consistency, or just let filtfilt handle.

% Nyquist Frequency
fn = fs / 2;

% Order
n = 3;

%% 0. Continuous (0 - 0.005 Hz) -> Lowpass
f_cutoff = 0.005;
Wn = f_cutoff / fn;
[b, a] = butter(n, Wn, 'low');
decomposed.continuous = filtfilt(b, a, signal);
%% 1. Infraslow (0.005 - 0.1 Hz) -> Bandpass
f_range = [0.005, 0.1];
Wn = f_range / fn;
[b, a] = butter(n, Wn, 'bandpass');
decomposed.isf = filtfilt(b, a, signal);
%% 2. VLF (0.1 - 0.3 Hz) -> Bandpass
f_range = [0.1, 0.3];
Wn = f_range / fn;
[b, a] = butter(n, Wn, 'bandpass');
decomposed.vlf = filtfilt(b, a, signal);

%% 3. LF (0.3 - 1 Hz) -> Bandpass
f_range = [0.3, 1.0];
Wn = f_range / fn;
[b, a] = butter(n, Wn, 'bandpass');
decomposed.lf = filtfilt(b, a, signal);

%% 4. HF (1 - 1.5 Hz) -> Bandpass
f_cutoff = 1.0; % Just close to Nyquist
Wn = f_cutoff / fn;
[b, a] = butter(n, Wn, 'high');
decomposed.hf_residual = filtfilt(b, a, signal);

end
