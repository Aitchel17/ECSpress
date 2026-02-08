function plot_individual_thickness(secondary_struct,ssidx,opt)
%PLOT_INDIVIDUAL_THICKNESS Summary of this function goes here
%   Detailed explanation goes here
figure()
t =tiledlayout("flow",TileSpacing="tight")
for ssidx = 1:8
scale = strsplit(secondary_struct(ssidx).infodict("objpix"));
scale = str2double(scale(1));
nexttile
tmp.taxis = secondary_struct(ssidx).fwhmline.pax.t_axis(2:end);
plot(tmp.taxis,secondary_struct(ssidx).thickness(1).thickness*scale,'Color','c')
hold on
plot(tmp.taxis,secondary_struct(ssidx).thickness(2).thickness*scale,'Color','g')
plot(tmp.taxis,secondary_struct(ssidx).thickness(3).thickness*scale,'Color','r')
plot(tmp.taxis,secondary_struct(ssidx).thickness(4).thickness*scale,'Color','w')
hold off
ax = gca;
set(ax,'XColor',[0 0 0])
set(ax,'YColor',[0 0 0])
set(ax,'Color',[0 0 0])
axis tickaligned
axis tight
end
end

