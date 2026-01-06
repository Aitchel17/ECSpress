%%
baseecogspectrum = primary_datastruct.analog.ecogspectrum;
%%
baseecogspectrum.log_norm_spectrum = baseecogspectrum.log_norm_spectrum(:,100:200);
%%
baseecogspectrum.t_axis = baseecogspectrum.t_axis(:,100:200);
%%
normecogspectrum = primary_datastruct.analog.ecogspectrum;
%%
basespectrum = baseecogspectrum.log_norm_spectrum;
basespectrum = mean(basespectrum,2);
%%
normecogspectrum.log_norm_spectrum = normecogspectrum.log_norm_spectrum./basespectrum;
%%
figure()
%%
plot_ecogspectrum(normecogspectrum)