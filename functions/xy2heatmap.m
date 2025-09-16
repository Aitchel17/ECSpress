function heatdata = xy2heatmap(x_data,y_data,resolution)
%HEATMAP_THICKNESS Summary of this function goes here
%   Detailed explanation goes here
[xy_counts, xedges, yedges] = histcounts2(x_data*resolution, y_data*resolution, 'BinWidth', [resolution, resolution],'Normalization','percentage');
heatdata.xy_counts = xy_counts';
% Get bin centers
heatdata.x_centers = (xedges(1:end-1) + xedges(2:end)) / 2;
heatdata.y_centers = (yedges(1:end-1) + yedges(2:end)) / 2;
% Apply log scale
heatdata.log_xycounts = log10(xy_counts' + 1); % Avoid log(0)
end

