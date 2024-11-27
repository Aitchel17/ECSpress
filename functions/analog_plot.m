function [fig] = analog_plot(analog)
%ANALOG_PLOT Summary of this function goes here
%   Detailed explanation goes here
% plot graph


% Make figure
fig = figure('Name','Analog signal','NumberTitle','off');
% Define the number of subplots
num_subplots = length(fieldnames(analog))/2;
current_subplot = 0; % Counter for subplots


% Ball plot
if isfield(analog,'raw_Ball')
current_subplot = current_subplot + 1;
ax= gca;
ax.Layer='top';
subplot(num_subplots, 1, current_subplot);
plot(analog.ds_ball(1,:)/60,analog.ds_ball(2,:))
ylabel('Ball (V)')
xlim([0 30])
xlabel('min')
set(gca,'fontsize',14)
end
% EMG plot
%
if isfield(analog,'raw_EMG')
xlabel('')
set(gca,'xtick',[])
current_subplot = current_subplot + 1;
subplot(num_subplots, 1, current_subplot);
plot(analog.ds_EMG(1,:)/60,analog.ds_EMG(2,:))
xlim([0 30])
ylabel('EMG (V)')
xlabel('min')
set(gca,'fontsize',14)
end
%
if isfield(analog,'raw_ECoG')
xlabel('')
set(gca,'xtick',[])
current_subplot = current_subplot + 1;
subplot(num_subplots, 1, current_subplot);
surface(analog.ecog_spectrum.t_axis/60,analog.ecog_spectrum.f_axis, ...
    zeros(size(analog.ecog_spectrum.log_norm_spectrum)),(analog.ecog_spectrum.log_norm_spectrum),'LineStyle','None')
set(gca,'YScale','log');
xlim([0 30])
ylabel('Freq (Hz)')
set(gca,'fontsize',14)
xlabel('min')
end
end

