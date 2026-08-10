function gpt_testfunc(x_data, x_data2, y_data, scale, figtitle, cmap, cmap2)
    arguments
        x_data
        x_data2
        y_data
        scale
        figtitle
        cmap
        cmap2
    end

    %% Common bin edges (cover both datasets)
    all_x = [x_data; x_data2];
    scaled_x = all_x / scale;
    scaled_y = y_data / scale;

    max_x = ceil(max(scaled_x));
    max_y = ceil(max(scaled_y));
    bin_width = 1/scale;

    xedges = 0:bin_width:max_x;
    yedges = 0:bin_width:max_y;

    %% Histogram 1 (main)
    scaled_x1 = x_data / scale;
    [x_counts1, ~, ~] = histcounts2(scaled_x1, scaled_y, xedges, yedges);
    x_centers = (xedges(1:end-1) + xedges(2:end)) / 2;
    y_centers = (yedges(1:end-1) + yedges(2:end)) / 2;
    log_values1 = log10(x_counts1 + 1);

    %% Histogram 2 (overlay)
    scaled_x2 = x_data2 / scale;
    [x_counts2, ~, ~] = histcounts2(scaled_x2, scaled_y, xedges, yedges);
    log_values2 = log10(x_counts2 + 1);

    %% Plot on same axes
    figure()
    ax1 = gca;
    hold on

    % Plot first histogram
    img1 = imagesc(ax1, y_centers, x_centers, log_values1);
    set(ax1, 'YDir', 'normal')
    colormap(ax1, cmap)
    set(img1, 'AlphaData', log_values1 > 0)

    % Overlay second histogram
    img2 = imagesc(ax1, y_centers, x_centers, log_values2);
    colormap(ax1, cmap2)
    set(img2, 'AlphaData', 0.5 * (log_values2 > 0)) % Semi-transparent overlay

    % Final settings
    colorbar
    xlabel('y (scaled)')
    ylabel('x (scaled)')
    title(figtitle)
    axis equal
    xlim([0 max(y_centers)])
    ylim([0 max(x_centers)])
end