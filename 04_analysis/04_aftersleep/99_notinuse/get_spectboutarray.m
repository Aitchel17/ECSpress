function spectrogramlist = get_spectboutarray(bouts,signal,fs,eventfreq_thr)
arguments
    bouts cell % n*1 cell array
    signal double % 1D vector
    fs double {mustBePositive} % sampling frequency
    eventfreq_thr double {mustBePositive} % Hz
end

cnt = 1;

eventlen_thr = 2*1/eventfreq_thr; % seconds, twice long
eventlen_thr = eventlen_thr*fs;

for i = 1:numel(bouts)
    loc = bouts{i};
    if numel(loc)>eventlen_thr
        disp(i)
        % Preprocessing: Sgolay filter + Detrend
        cropped_signal = signal(loc);
        % bv_signal = detrend(sgolayfilt(bv_signal, 3, 5));
        spectrogram = get_spectrogram(fs,cropped_signal);
        spectrogram.boutidx = i;
        spectrogram.count = cnt;
        spectrogramlist(cnt) = spectrogram;
        cnt = cnt + 1;
    end
end
end