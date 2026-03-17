function results = analyzeUltrasoundAudio(filepath)

% ANALYZE_ULTRASOUND
% Advanced ultrasonic diagnostics for AudioMoth v1.2.0 (Medium Gain)
%
% Outputs:
%   - Spectrogram
%   - Time-averaged spectrum
%   - Harmonic groups (>1 freq only)
%   - SPL of 20–100 kHz band (sound pressure level actual noise level, micropascal)
%   - SPL of strongest ultrasonic tone (peak amplitude)

%% ---- Input ----
if nargin < 1 || isempty(filepath)
    [file, path] = uigetfile('*.wav', 'Select AudioMoth WAV file');
    if isequal(file,0)
        error('No file selected.');
    end
    filepath = fullfile(path, file);
end

%% ---- Load ----
[x, Fs] = audioread(filepath);
if size(x,2) > 1
    x = mean(x,2);
end

nyquist = Fs/2;
maxFreq = min(100e3, nyquist);

fprintf('\nFile: %s\n', filepath);
fprintf('Sampling Rate: %.0f Hz\n', Fs);
fprintf('Duration: %.2f sec\n\n', length(x)/Fs);

%% ---- Spectrogram ----
windowLength = 2048;
overlap = 0.75;
nfft = 4096;

window = hann(windowLength);
noverlap = round(windowLength * overlap);

[S,F,T] = spectrogram(x, window, noverlap, nfft, Fs);
S_dB = 20*log10(abs(S) + eps);

figure;
imagesc(T, F/1000, S_dB);
axis xy;
ylim([0 maxFreq/1000]);
colormap turbo;
cb = colorbar;
ylabel(cb,'Magnitude (dB)');
xlabel('Time (s)');
ylabel('Frequency (kHz)');
title('Spectrogram (0–100 kHz)');
% clim([max(S_dB(:))-60 max(S_dB(:))]);
clim([max(S_dB(:))-70 max(S_dB(:))]);
% clim([max(S_dB(:))-80 max(S_dB(:))]);

%% ---- Time-Averaged Spectrum ----
meanSpectrum = mean(abs(S).^2, 2);
meanSpectrum_dB = 10*log10(meanSpectrum + eps);

figure;
plot(F/1000, meanSpectrum_dB, 'LineWidth',1);
xlim([0 maxFreq/1000]);
xlabel('Frequency (kHz)');
ylabel('Power (dB)');
title('Time-Averaged Power Spectrum');
grid on;

%% ---- Peak Detection ----
[pks, locs] = findpeaks(meanSpectrum_dB, ...
    'MinPeakProminence',5);

peakFreqs = F(locs);
peakPowers = meanSpectrum(locs);

hold on;
plot(peakFreqs/1000, pks, 'r.', 'MarkerSize',15);
hold off;

%% ---- 1️⃣ Harmonic Grouping (>1 frequency only) ----
fprintf('Harmonic Groups:\n');

tolerance = 0.02; % 2%
used = false(size(peakFreqs));
groups = {};

for i = 1:length(peakFreqs)
    if used(i), continue; end
    
    f0 = peakFreqs(i);
    groupIdx = i;
    
    for j = i+1:length(peakFreqs)
        ratio = peakFreqs(j)/f0;
        if abs(ratio - round(ratio)) < tolerance
            groupIdx(end+1) = j;
            used(j) = true;
        end
    end
    
    if length(groupIdx) > 1
        groups{end+1} = peakFreqs(groupIdx);
    end
end

if isempty(groups)
    fprintf('  No harmonic series detected.\n');
else
    for g = 1:length(groups)
        fprintf('  Series: ');
        fprintf('%.2f ', groups{g}/1000);
        fprintf('kHz\n');
    end
end

%% ---- SPL Calibration Constants (v1.2.0 Medium Gain) ----
S_effective = 0.100;     % V/Pa (includes +18 dB gain)
Vref = 0.707;            % Digital full-scale RMS volts
Pref = 20e-6;            % 20 µPa reference

%% ---- 2️⃣ SPL of 20–100 kHz Band ----
bandIdx = F >= 20e3 & F <= maxFreq;
bandPower = mean(meanSpectrum(bandIdx));
bandRMS_digital = sqrt(bandPower);

V_rms_band = bandRMS_digital * Vref;
P_rms_band = V_rms_band / S_effective;

SPL_band = 20*log10(P_rms_band / Pref);

fprintf('\nEstimated SPL (20–100 kHz): %.1f dB SPL\n', SPL_band);

%% ---- 3️⃣ SPL of Strongest Ultrasonic Tone ----
ultraIdx = peakFreqs >= 20e3 & peakFreqs <= maxFreq;

if any(ultraIdx)
    [~, strongestIdx] = max(peakPowers(ultraIdx));
    ultraPeaks = peakFreqs(ultraIdx);
    strongestFreq = ultraPeaks(strongestIdx);
    
    bw = 500; % 500 Hz integration bandwidth
    idx = F > (strongestFreq-bw) & F < (strongestFreq+bw);
    tonePower = mean(meanSpectrum(idx));
    toneRMS_digital = sqrt(tonePower);
    
    V_rms_tone = toneRMS_digital * Vref;
    P_rms_tone = V_rms_tone / S_effective;
    SPL_tone = 20*log10(P_rms_tone / Pref);
    
    fprintf('Strongest Ultrasonic Tone: %.2f kHz\n', strongestFreq/1000);
    fprintf('Estimated SPL at %.2f kHz: %.1f dB SPL\n', ...
        strongestFreq/1000, SPL_tone);
else
    fprintf('No ultrasonic peaks detected.\n');
    SPL_tone = NaN;
    strongestFreq = NaN;
end

%% ---- Return ----
results.harmonicGroups = groups;
results.SPL_20_100kHz = SPL_band;
results.strongestFrequency_kHz = strongestFreq/1000;
results.SPL_strongestTone = SPL_tone;

end