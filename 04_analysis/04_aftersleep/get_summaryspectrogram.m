function spectrogram = get_summaryspectrogram(spectrogramlist)
common_F = 0:0.01:session.pax_fwhm.param.fs/2; % Common frequency axis
all_S = [];

for i = 1:numel(rem.thickness.bv_spectrogramlist)
    if isempty(rem.thickness.bv_spectrogramlist(i).S), continue; end
    % Interpolate to common axis
    S_interp = interp1(rem.thickness.bv_spectrogramlist(i).F, rem.thickness.bv_spectrogramlist(i).S, common_F, 'linear', 'extrap');
    all_S = [all_S; S_interp(:)']; % Stack rows
end

rem.thickness.bv_spectrogram.S = mean(all_S, 1, 'omitnan');
rem.thickness.bv_spectrogram.F = common_F;
end