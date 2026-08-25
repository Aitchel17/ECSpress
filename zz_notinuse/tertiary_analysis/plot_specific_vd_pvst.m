function [outputArg1,outputArg2] = plot_specific_vd_pvst(secondary_struct,idx)
%PLOT_SPECIFIC_VD_PVST Summary of this function goes here
%   Detailed explanation goes here

colors = lines(size(secondary_struct,2));  % N개의 색 생성
for ssidx = 1:size(secondary_struct,2)
    hold on
    scatter(secondary_struct(ssidx).heatdata(idx).x_centers_aligned, secondary_struct(ssidx).heatdata(idx).modespvs, 100 ,"x",'MarkerEdgeColor', colors(ssidx,:), LineWidth=1)


    y_fit = secondary_struct(ssidx).heatdata(idx).intercept + secondary_struct(ssidx).heatdata(idx).slope*secondary_struct(ssidx).heatdata(idx).x_centers_aligned;
    plot(secondary_struct(ssidx).heatdata(3).x_centers_aligned,y_fit,'LineWidth',2,Color=colors(ssidx,:))
end
axis tight
ax = gca;              % Get current axes
ax.FontSize = 20;
xlim([-3 8])
ylim([-8 3])
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color',[0 0 0])
set(ax,'Color',[0 0 0])
end

