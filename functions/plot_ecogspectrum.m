function plot_ecogspectrum(ecog_spectrum)
    fig = figure(Name='ECoG spectrum');
    ax = axes('Parent', fig);
    ecog_spectrum.zdata = zeros(size(ecog_spectrum.log_norm_spectrum));
    %%
    surface(ax, ecog_spectrum.t_axis / 60, ...
        ecog_spectrum.f_axis, ...
        zeros(size(ecog_spectrum.log_norm_spectrum)), ...
        ecog_spectrum.log_norm_spectrum, ...
        'LineStyle', 'none');
    
    set(ax, 'YScale', 'log');
    xlim(ax, [0 30]);
    ylabel(ax, 'Freq (Hz)');
    xlabel(ax, 'min');
    set(ax, 'FontSize', 14);
end

