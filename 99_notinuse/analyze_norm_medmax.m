function analyze_norm_medmax(image_ch1,x_scatterlabel,image_ch2,y_scatterlabel,image_x,offset,analog_x,w_stim)
%ANALYZE_NORM_MEDMAX Summary of this function goes here
%   Detailed explanation goes here

% x axis

image_x = image_x(offset:end);

% median filter
image_ch2 = medfilt1(image_ch2,10);
image_ch1 = medfilt1(image_ch1,10);

image_ch2= image_ch2(offset:end);
image_ch1= image_ch1(offset:end);


%% normalize by median value
%image_ch2 = image_ch2-median(image_ch2);
%image_ch1 = image_ch1-median(image_ch1);
image_ch2 = image_ch2/median(image_ch2);
image_ch1 = image_ch1/median(image_ch1);
image_ch2 = image_ch2 -1;
image_ch1 = image_ch1 -1;

[b, stats] = robustfit(image_ch1, image_ch2);
y_fit = b(1) + b(2)*image_ch1;
rsquare = corr(image_ch1,b(1)+b(2)*image_ch2)^2;

figure('Position', [200, 200, 1200, 400])
plot(image_x,image_ch1,Color='r');
hold on;
plot(image_x,image_ch2,Color='g');

%% Shade stimulus epochs clearly
stim_binary = w_stim > (0.5 * max(w_stim));
dStim = diff([0; stim_binary(:); 0]);
stim_onsets = analog_x(dStim(1:end-1) == 1);
stim_offsets = analog_x(dStim(1:end-1) == -1);

y_limits = ylim;
for k = 1:length(stim_onsets)
    fill([stim_onsets(k), stim_offsets(k), stim_offsets(k), stim_onsets(k)],...
         [y_limits(1), y_limits(1), y_limits(2), y_limits(2)],...
         [0.9 0.9 0.9], 'FaceAlpha', 0.9, 'EdgeColor', 'none');
end

xlim([image_x(offset),image_x(end)])
xlabel('Time (sec)');
ylabel('Normalized signal intensity F-Fmed/Fmed');
figure('Position', [200, 200, 400, 400])
scatter(image_ch1,image_ch2)
hold on;
% Overlay robust fit line
plot(image_ch1, y_fit, 'r-', 'LineWidth', 2);

% Optional: add labels and legend
xlabel(x_scatterlabel);
ylabel(y_scatterlabel);
title([x_scatterlabel '  vs  ' y_scatterlabel]);
fit_legend = ['slope: ',num2str(round(b(2),2)),' rsquare: ',num2str(round(rsquare,2))]

legend('At each time point', fit_legend, 'Location', 'Best');

grid on;
hold off;
end

