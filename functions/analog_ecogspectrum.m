function [ecog_spectrum] = pre_ecogspectrum(samplingfreq,ECoG)

% ECoG processing
    detrended_ECoG = detrend(ECoG,'constant');
    params5.tapers = [5,9];
    params5.Fs = samplingfreq;
    params5.fpass = [1,100];
    movingwin5 = [5,1/5];
    [S5,T5,F5] = mtspecgramc(detrended_ECoG,movingwin5,params5);
    % Output:
    %       S       (spectrum in form time x frequency x channels/trials if trialave=0; 
    %               in the form time x frequency if trialave=1)
    %       t       (times)
    %       f       (frequencies)
    %       Serr    (error bars) only for err(1)>=1
    
    normS5 = (S5-min(S5,[],'all'))./(max(S5,[],'all')-min(S5,[],'all'));
    % output
    ecog_spectrum.log_norm_spectrum = log(normS5)';
    ecog_spectrum.t_axis = T5;
    ecog_spectrum.f_axis = F5;

end

