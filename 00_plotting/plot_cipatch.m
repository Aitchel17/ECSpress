function plot_cipatch(ax, f, data_matrix, color_char, label)
% Calculate stats
mu = mean(data_matrix, 2, 'omitnan')';
sigma = std(data_matrix, 0, 2, 'omitnan')';
sem = sigma ./ sqrt(size(data_matrix, 2));
ci = 1.96 * sem;

upper = mu + ci;
lower = mu - ci;

% Plot shaded CI

f = log(f);

fill(ax, [f, fliplr(f)], [upper, fliplr(lower)], color_char, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', [label ' 95% CI']);
% Plot mean
plot(ax, f, mu, color_char, 'LineWidth', 2, 'DisplayName', [label ' Mean']);
end