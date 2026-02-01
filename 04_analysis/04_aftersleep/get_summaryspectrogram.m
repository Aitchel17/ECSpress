function summary_spect = get_summaryspectrogram(spectrogramlist,fs)
summary_spect = struct();
common_F = 0:0.01:fs/2; % Common frequency axis
all_S = [];

for i = 1:numel(spectrogramlist)
    S_interp = interp1(spectrogramlist(i).F, spectrogramlist(i).S, common_F, 'linear', 'extrap');
    all_S = [all_S; S_interp(:)']; % Stack rows
end

summary_spect.all_S = all_S;
summary_spect.S = mean(all_S, 1, 'omitnan');
summary_spect.F = common_F;

% Calculate Confidence Interval (t-distribution)
n = size(all_S, 1);
if n > 1
    sem = std(all_S, 0, 1, 'omitnan') ./ sqrt(n);
    t_score = tinv(0.975, n-1); % 95% CI
    ci = t_score * sem;

    summary_spect.sem = sem;
    summary_spect.ci95 = ci;
else % REM...
    summary_spect.sem = NaN(size(summary_spect.S));
    summary_spect.ci95 = NaN(size(summary_spect.S));
end

end