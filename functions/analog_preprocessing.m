function [result] = pre_analogprocessing(analog_signal,parameter,downsample_rate,medianfilt)

analog_freq = parameter(1); % unit: Hz
analog_pixeldepth = parameter(2); % unit: bit
analog_inputrange = parameter(3); % unit: V

% Ball signal
downsample_signal = downsample(analog_signal,downsample_rate);
downsample_freq = analog_freq/downsample_rate;
medfilt_analog = medfilt1(downsample_signal,medianfilt);

% scale
rescale_signal = medfilt_analog*analog_inputrange/2^analog_pixeldepth;

% xaxis
downsample_duration = length(downsample_signal)/downsample_freq;

downsample_xaxis = linspace(0,downsample_duration,length(downsample_signal));
disp(analog_pixeldepth(1:end))
disp(class(downsample_xaxis))
result = cat(1,downsample_xaxis,rescale_signal);

end    


