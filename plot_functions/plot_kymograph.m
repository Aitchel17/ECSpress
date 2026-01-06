function plot_kymograph(fig_name,kymograph_data, opts)
arguments
    fig_name string
    kymograph_data
    opts.idx_data = [];
    opts.blackbackground logical = false
    opts.resolution = 1; % Default y axis unit: pixels
    opts.fps = 1; % Default x axis unit:frame
end

% Input: 
% 1. figure name : Name of figure and savename for eps file
% 2. kymograph data: m x time CData goes to imagesc
% 3. idx_data: 

%% Default setting
save_path = 'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\08_figures';
monitor_xyinch = [27 5];
xy_sizeinch = [5 2];
fontsize = 12;

%% calculate x axis and y axis

x_axis = linspace(0,(size(kymograph_data,2)-1)/opts.fps,size(kymograph_data,2)); % start from 0
y_axis = linspace(0,(size(kymograph_data,1)-1)*opts.resolution,size(kymograph_data,1)); % start from 0
%%



%%
fig = figure("Name",fig_name);
set(fig,'Units','inches',...
    "Position",[monitor_xyinch(1) monitor_xyinch(2) xy_sizeinch(1) xy_sizeinch(2)]) % 
ax = axes(fig);
%% make figure 
imagesc(ax,x_axis,y_axis,kymograph_data);
colormap('gray')
%%
hold on
for i = 1:2
plot(x_axis,opts.idx_data{i,1}*opts.resolution, opts.idx_data{i,2});
end


%%
xlabel(ax, 'Time (s)',  'FontName','Arial','FontSize', fontsize, 'FontWeight','normal', 'Color', 'k');
ylabel(ax, 'Diameter', 'FontName','Arial','FontSize', fontsize, 'FontWeight','normal', 'Color',  'k');
%%

set(ax,'Box','off', ...
        'TickDir','out', ...
        'LineWidth',1, ...
        'FontName','Arial', ...                  
        'FontSize',fontsize, ...
        'Layer','top')
%%
if opts.blackbackground
    fig_color = 'k';
    axis_color = 'w';
else
    fig_color = 'w';
    axis_color = 'k';
end

set(fig,'Color',fig_color)
set(ax, 'Color',fig_color, ...
        'XColor', axis_color, ...
        'YColor', axis_color);

end

