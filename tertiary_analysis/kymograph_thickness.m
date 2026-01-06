ssidx = 8;
clee = color_lee;
pax = secondary_struct(ssidx).fwhmline.pax;
scale = strsplit(secondary_struct(ssidx).infodict("objpix"));
scale = str2double(scale(1));
fps = str2double(secondary_struct(ssidx).infodict("savefps"))/3;
airtable = secondary_struct(ssidx).analog.airtable;
analog_t = secondary_struct(ssidx).analog.data.taxis;
%%

starttime = {airtable.StartTime};
starttime = starttime{1};

%%
endtime = {airtable.EndTime};
endtime = endtime{1};
%%

axsize = [1 1 8 3];
figsize = [1 1 12 5];
fig_bvkgph  = figure(Name='fig_bvkgph',NumberTitle='off');
set(gcf,'Units','inches','Position',figsize); 
fig_pvskgph = figure(Name='fig_pvskgph',NumberTitle='off');
set(gcf,'Units','inches','Position',figsize); 
fig_bvb     = figure(Name='fig_bvb',NumberTitle='off');
set(gcf,'Units','inches','Position',figsize); 
fig_pvsb    = figure(Name='fig_pvsb',NumberTitle='off');
set(gcf,'Units','inches','Position',figsize); 

%% whiskerstimulation point plot
fig = figure()

hold on
% Dashed magenta line, 50% transparent
xline(starttime, ...
    'LineStyle', '--', ...
    'Color', [1 1 1 0.5], ...   % [R G B Alpha]
    'LineWidth', 1);

% Dashed red line, 30% transparent
xline(endtime, ...
    'LineStyle', '--', ...
    'Color', [1 1 1 0.3], ...
    'LineWidth', 1);
set(ax,'XColor',[0 0 0])
set(ax,'YColor',[0 0 0])
set(ax,'Color',[0 0 0])

axis off
fig.Position = [100 100 1200 400];
htrigframe = getframe(gca);
htrigimg = frame2im(htrigframe);
%% kymograph plot
fig = figure()
xline(starttime, ...
    'LineStyle', '--', ...
    'Color', [1 0 1 0.5], ...   % [R G B Alpha]
    'LineWidth', 1);

% Dashed red line, 30% transparent
xline(endtime, ...
    'LineStyle', '--', ...
    'Color', [1 0 0 0.3], ...
    'LineWidth', 1);
hold on
imagesc(hmap_x,hmap_y, upbv_mask + uppvs_mask);
set(gca, 'YDir', 'normal');
colormap(clee.green_g);
xlim([0 1800])

fig.Position = [100 100 1200 400];
hbvframe = getframe(gca);
hbvimg = frame2im(hbvframe);

fig = figure() 
imagesc(hmap_x,hmap_y, downbv_mask + downpvs_mask);
set(gca, 'YDir', 'normal');
colormap(clee.red_g);
xlim([0 1800])

fig.Position = [100 100 1200 400];
hcsfframe = getframe(gca);
hcsfimg = frame2im(hcsfframe);

%% integration
figure()
imshow(hbvimg+hcsfimg)
%%

%% Thickness plot
fig = figure();
%%
hold on
ssidx = 8
tmp.taxis = secondary_struct(ssidx).fwhmline.pax.t_axis(2:end);
plot(tmp.taxis,secondary_struct(ssidx).thickness(1).thickness*scale, Color='c')
hold on
plot(tmp.taxis,secondary_struct(ssidx).thickness(2).thickness*scale, Color='g')
plot(tmp.taxis,secondary_struct(ssidx).thickness(3).thickness*scale, Color= 'r')
%plot(tmp.taxis,secondary_struct(ssidx).thickness(4).thickness*scale, Color='w')
axis tickaligned
axis tight
fig.Position = [100 100 1200 400];
ax = gca;
set(ax,'XColor',[0 0 0])
set(ax,'YColor',[0 0 0])
set(ax,'Color',[0 0 0])
%%
axis on



%%
axis off
xlim([400 1851])

%%

fig = figure() 
imagesc(hmap_x,hmap_y, uppvs_mask + downpvs_mask);
colormap(clee.magenta_g);
xlim([0 1800])

fig.Position = [100 100 1200 400];
%%
figure()
colormap(clee.red_g)



%%
medfilt_p = 19;
filtered_upbvidx = medfilt1(pax.idx.bv_upperboundary,medfilt_p);
upbv_mask = enlargedmask(filtered_upbvidx,pax.mask.rowidx);
filtered_downbvidx = medfilt1(pax.idx.bv_lowerboundary,medfilt_p);
downbv_mask = enlargedmask(filtered_downbvidx,pax.mask.rowidx);
bvrecon = pax.reconstruction(upbv_mask + downbv_mask);
%%
filtered_uppvs_idx = medfilt1(pax.idx.pvs_upperboundary,medfilt_p);
uppvs_mask = enlargedmask(filtered_uppvs_idx,pax.mask.rowidx);
filtered_downpvs_idx = medfilt1(pax.idx.pvs_lowerboundary,medfilt_p);
downpvs_mask = enlargedmask(filtered_downpvs_idx,pax.mask.rowidx);
pvsrecon = pax.reconstruction(uppvs_mask + downpvs_mask);

%%
csfkymograph = pax.kymograph.kgph_normcsf(:,851:end);
bvkymograph = pax.kymograph.kgph_normbv(:,851:end);
hmap_y = linspace(-size(bvkymograph,1)*scale/2,size(bvkymograph,1)*scale/2,size(bvkymograph,1));
hmap_x = linspace(1,size(bvkymograph,2)/fps,size(bvkymograph,2));
%%
figure()
imagesc(hmap_x,hmap_y, bvkymograph)
colormap(clee.red_g)
set(gca, 'YDir', 'normal','Units','inches','Position',axsize)
clim(gca,[0 0.7])
hrframe = getframe(gca);
hrimg = frame2im(hrframe);
%%
figure()
hold on
plot(filtered_upbvidx)
plot(filtered_downbvidx)
%%



%%
figure()
plot(filtered_upbvidx+filtered_downbvidx-median(filtered_upbvidx+filtered_downbvidx))
hold on
plot(-filtered_upbvidx+filtered_downbvidx-median(-filtered_upbvidx+filtered_downbvidx))


%%
% Preprocessing for fwhm based detection, 3 frame averaging, 5 frame median filter
afterprocessing_avgproj = 3;
[tmp.gpstack, tmp.ntailfr] =  analyze_grouproject(vivo_50um1.stackch2,afterprocessing_avgproj,"mean");
tmp.medcsf = medfilt3(tmp.gpstack,[1 1 5]);
gausscsf= imgaussfilt(tmp.medcsf,1);
tmp.gpstack = analyze_grouproject(vivo_50um1.stackch1,afterprocessing_avgproj,"mean");
tmp.medbv = medfilt3(tmp.gpstack,[1 1 5]);
gaussbv= imgaussfilt(tmp.medbv,1);
%%
util_checkstack(gaussbv.*~bvrecon)


%%
figure()
plot(medfilt1(pax.idx.bv_upperboundary,medfilt_p))
hold on
plot(medfilt1(pax.idx.bv_lowerboundary,medfilt_p))
%%
figure()
plot(medfilt1(pax.idx.bv_upperboundary,medfilt_p)+medfilt1(pax.idx.bv_lowerboundary,medfilt_p))
hold on

plot(medfilt1(pax.idx.bv_lowerboundary,medfilt_p))

%%

util_checkstack(bvrecon)



function mask = enlargedmask(idx,rowidx)
    mask = zeros(size(rowidx));
    mask(rowidx == idx) = 2^16-1;
    mask(rowidx == idx+1) = 2^15;
    mask(rowidx == idx-1) = 2^15;
end
