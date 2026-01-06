function plot_2dhistdistribution(x_data,y_data,scale,figtitle,cmap)
    arguments
        x_data
        y_data
        scale
        figtitle
        cmap
    end
%PLOT_2DHISTDISTRIBUTION Summary of this function goes here
%   Detailed explanation goes here
% 2d histogram plot
    % Create 2D histogram with fixed bin width
    scaled_x = x_data/scale;
    scaled_y = y_data/scale;
    [x_counts, xedges, yedges] = histcounts2(x_data/scale, y_data/scale, 'BinWidth', [1/scale, 1/scale]);
    %
    % Get bin centers
    x_centers = (xedges(1:end-1) + xedges(2:end)) / 2;
    y_centers = (yedges(1:end-1) + yedges(2:end)) / 2;
    
    % Apply log scale
    log_values = log10(x_counts + 1); % Avoid log(0)
    
    figure()
    imagesc(y_centers, x_centers, log_values)
    set(gca, 'YDir', 'normal')
    axis 
    colorbar
    xlabel('Vessel diameter changes (um)')
    ylabel('PVS diameter changes (um)')
    title(figtitle)
    colormap(cmap)
end

