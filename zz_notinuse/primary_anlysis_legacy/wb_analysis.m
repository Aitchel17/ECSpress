format compact
marker_mw= [250,150,100,75,50,37,25,20,15];
% copy and paste ROI measure all, ROI order : 250kDA to 10 kDA, and last
% one to be target. 
imj_measure = [
1	11546362.722	17728.451	15958	18947	-90	4572.557
2	16716375.881	17607.022	15710	19014	-90	6619.970
3	25677732.023	17586.864	15710	19014	-90	10168.820
4	32915750.446	18608.847	15710	48830	-90	13035.199
5	52217132.906	19276.103	15710	51910	-90	20678.876
6	66693169.751	18943.997	15710	51910	-90	26411.634
7	87200888.615	19827.242	15710	51910	-90	34533.040
8	99781253.968	20007.932	15569	59181	-90	39515.080
9	107363939.935	19805.210	15569	59181	-90	42517.953
10	110293614.058	19732.177	15569	59181	-90	43678.154
11	86856221.071	18229.973	15866	31835	-90	34396.546
];
%%
marker_lw = imj_measure(1:end-2,7);
dye_front = imj_measure(end-1,7);
gfp_lw = imj_measure(end,7);
marker_rf = marker_lw./ dye_front;
marker_logmw = log10(marker_mw);
linearmodel = fitlm(marker_rf,marker_logmw);
%%
gfp_rf = gfp_lw/dye_front;
gfp_logmw = linearmodel.Coefficients.Estimate(2)*gfp_rf+linearmodel.Coefficients.Estimate(1);
gfp_mw = 10^gfp_logmw



%% 그래프
figure
hold on
% 데이터 포인트
plot(marker_rf, marker_logmw, 'o', 'MarkerFaceColor', 'k', 'DisplayName', 'Marker proteins')

% 트렌드라인
x_fit = linspace(min(marker_rf), max(marker_rf), 100);
y_fit = linearmodel.Coefficients.Estimate(2)*x_fit + linearmodel.Coefficients.Estimate(1);
plot(x_fit, y_fit, '-', 'LineWidth', 1.5, 'DisplayName', 'Linear fit')

% GFP 위치 빨간색
plot(gfp_rf, gfp_logmw, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'IgKL-mNG')

xlabel('Relative front (Rf)')
ylabel('MW (kDa)')
title('Standard Curve and IgKL-mNG estimation')
legend
%%
x_limits = xlim;
plot(x_limits, [gfp_logmw, gfp_logmw], 'r--', 'LineWidth', 1, 'HandleVisibility', 'off')
%%
hold off
%%
yticks_vals = get(gca, 'YTick');
ytick_labels = arrayfun(@(x) sprintf('%.0f', 10.^x), yticks_vals, 'UniformOutput', false);
set(gca, 'YTickLabel', ytick_labels)
%%
text_position_x = mean(x_limits); % 중앙에 표시
 a= text(text_position_x, gfp_logmw + 0.05, sprintf('Estimated MW: %.2f kDa', gfp_mw), ...
    'Color', 'r', 'FontWeight', 'normal', 'HorizontalAlignment', 'center')
%%
delete(a)