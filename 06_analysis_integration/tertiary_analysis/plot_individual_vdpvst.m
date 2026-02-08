function plot_individual_vdpvst(secondary_struct)
t = tiledlayout('flow');
color = ["r","g","w"];
for ssidx=1:8
nexttile
hold on
for idx = 1:3

        y_fit = secondary_struct(ssidx).heatdata(idx).intercept + secondary_struct(ssidx).heatdata(idx).slope*secondary_struct(ssidx).heatdata(idx).x_centers_aligned;
        scatter(secondary_struct(ssidx).heatdata(idx).x_centers_aligned, secondary_struct(ssidx).heatdata(idx).modespvs,"x",'MarkerEdgeColor', color(idx))
        plot(secondary_struct(ssidx).heatdata(idx).x_centers_aligned,y_fit,'LineWidth',2,color=color(idx))
end

axis tight equal
ax = gca;              % Get current axes
ax.FontSize = 14;

set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
    if secondary_struct(ssidx).heatdata(1).slope > secondary_struct(ssidx).heatdata(2).slope
        disp(ssidx)
        secondary_struct(ssidx).heatdata = secondary_struct(ssidx).heatdata([2,1,3]);
    end
hold off
end