figure()
t =tiledlayout("flow",TileSpacing="tight")
for ssidx = 8:8
pax = secondary_struct(8).fwhmline.pax;

scale = strsplit(secondary_struct(ssidx).infodict("objpix"));
scale = str2double(scale(1));

fps = str2double(secondary_struct(ssidx).infodict("savefps"))/3;
%%
fstart = 1;
fend = 3749;
nexttile
tmp.taxis = secondary_struct(ssidx).fwhmline.pax.t_axis(2:end);
plotdata_c = medfilt1(secondary_struct(ssidx).thickness(1).thickness(fstart:fend),30);
plotdata_w = medfilt1(secondary_struct(ssidx).thickness(4).thickness(fstart:fend),30);


tmp.taxis = tmp.taxis(fstart:fend);

plot(tmp.taxis,plotdata_c*scale,'Color','c')
hold on
%plot(tmp.taxis,secondary_struct(ssidx).thickness(2).thickness(fstart:fend)*scale,'Color','g')
%plot(tmp.taxis,secondary_struct(ssidx).thickness(3).thickness(fstart:fend)*scale,'Color','r')
plot(tmp.taxis, plotdata_w*scale,'Color','w')
hold off
ax = gca
set(ax,'XColor',[0 0 0])
set(ax,'YColor',[0 0 0])
set(ax,'Color',[0 0 0])
axis tickaligned
axis tight
end

%% kymograph construction

%%
csfkymograph = pax.kymograph.kgph_normcsf(:,851:end);
bvkymograph = pax.kymograph.kgph_normbv(:,851:end);
hmap_y = linspace(-size(bvkymograph,1)*scale/2,size(bvkymograph,1)*scale/2,size(bvkymograph,1));
hmap_x = linspace(1,size(bvkymograph,2)/fps,size(bvkymograph,2));

%%

figure()
imagesc(bvkymograph)