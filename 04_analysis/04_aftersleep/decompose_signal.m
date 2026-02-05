function decomposed = decompose_signal(signal, fs)
%DECOMPOSE_SIGNAL Decompose signal into frequency bands using Butterworth filters
%   Bands:
%       Continuous (0-0.1 Hz)
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

%% 0. Continuous (0 - 0.01 Hz) -> Lowpass
f_cutoff = 0.01;
Wn = f_cutoff / fn;
[b, a] = butter(n, Wn, 'low');
decomposed.continuous = filtfilt(b, a, signal);
%% 1. Infraslow (0.01 - 0.05) -> Bandpass
f_range = [0.01, 0.05];
Wn = f_range / fn;
[b, a] = butter(n, Wn, 'bandpass');
decomposed.isf = filtfilt(b, a, signal);
%% 2. VLF (0.05 - 0.3 Hz) -> Bandpass
f_range = [0.1, 0.3];
Wn = f_range / fn;
[b, a] = butter(n, Wn, 'bandpass');
decomposed.vlf = filtfilt(b, a, signal);

%% 3. LF (0.3 - 1 Hz) -> Bandpass
f_range = [0.3, 1.0];
Wn = f_range / fn;
[b, a] = butter(n, Wn, 'bandpass');
decomposed.lf = filtfilt(b, a, signal);

% Note: User mentioned 5 bands, but Fs=3Hz limits upper bands.
% HF (e.g. 1-4Hz) implies up to Nyquist 1.5Hz max.
% We can add a "High" band (1 - 1.5 Hz) if desired.
% f_range = [1.0, 1.49]; % Just close to Nyquist
% Wn = f_range / fn;
% [b, a] = butter(n, Wn, 'bandpass');
% decomposed.hf_residual = filtfilt(b, a, signal);

end
