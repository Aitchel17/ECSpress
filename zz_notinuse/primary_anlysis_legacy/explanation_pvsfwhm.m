%% visualization testment
pax = tmp.loadstruct.line_fwhms.pax;
fps = str2double(vivo_50um1.info.savefps)/afterprocessing_avgproj;
scale = str2double(vivo_50um1.info.objpix(1:end-3));
%%
pax = secondary_struct(8).fwhmline.pax;

scale = strsplit(secondary_struct(ssidx).infodict("objpix"));
scale = str2double(scale(1));
fps = str2double(secondary_struct(ssidx).infodict("savefps"))/3;
clee = color_lee;
%%
arrayfun(@(s) s.infodict("mdfName"), filtered_secondarystruct) 
%%
filtered_secondarystruct.infodict("mdfName")

%%
csfkymograph = pax.kymograph.kgph_normcsf(:,851:end);
bvkymograph = pax.kymograph.kgph_normbv(:,851:end);
hmap_y = linspace(-size(bvkymograph,1)*scale/2,size(bvkymograph,1)*scale/2,size(bvkymograph,1));
hmap_x = linspace(1,size(bvkymograph,2)/fps,size(bvkymograph,2));
%%
bv_line = pax.kymograph.kgph_bv(:,100);
%%
fig_lineprofile = figure('Name','lineprofile',NumberTitle='off');

%%
figure(fig_lineprofile)
plot(hmap_y,(bv_line-min(bv_line))/max(bv_line-min(bv_line)),'r','LineWidth',3)

%%
figure(fig_lineprofile)

csf_line = pax.kymograph.kgph_csf(:,500);
figure()
plot(hmap_y,(csf_line-min(csf_line))/max(csf_line-min(csf_line)),'g','LineWidth',3)
hold on
plot(hmap_y,(bv_line-min(bv_line))/max(bv_line-min(bv_line)),'r','LineWidth',3)

ax = gca;              % Get current axes
ax.FontSize = 14;
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color',[0 0 0])
set(ax,'Color',[0 0 0])


%% Save
bv_mask = zeros(size(pax.mask.rowidx));
bv_mask(pax.mask.rowidx == medfilt1(pax.idx.bv_upperboundary,9)) = 2^16-1;
bv_mask(pax.mask.rowidx == medfilt1(pax.idx.bv_lowerboundary,9)) = 2^16-1;
bv_mask(pax.mask.rowidx == medfilt1(pax.idx.bv_upperboundary,9)+1) = 2^15;
bv_mask(pax.mask.rowidx == medfilt1(pax.idx.bv_lowerboundary,9)+1) = 2^15;
bv_mask(pax.mask.rowidx == medfilt1(pax.idx.bv_upperboundary,9)-1) = 2^15;
bv_mask(pax.mask.rowidx == medfilt1(pax.idx.bv_lowerboundary,9)-1) = 2^15;
imagesc(bv_mask)
bvrecon = pax.reconstruction(bv_mask);
analyze_savesimpletif(bvrecon,'bv_boundary.tif')
csf_mask = zeros(size(pax.mask.rowidx));
csf_mask(pax.mask.rowidx == medfilt1(pax.idx.pvs_upperboundary,9)) = 2^16-1;
csf_mask(pax.mask.rowidx == medfilt1(pax.idx.pvs_lowerboundary,9)) = 2^16-1;
csf_mask(pax.mask.rowidx == medfilt1(pax.idx.pvs_upperboundary,9)+1) = 2^15;
csf_mask(pax.mask.rowidx == medfilt1(pax.idx.pvs_lowerboundary,9)+1) = 2^15;
csf_mask(pax.mask.rowidx == medfilt1(pax.idx.pvs_upperboundary,9)-1) = 2^15;
csf_mask(pax.mask.rowidx == medfilt1(pax.idx.pvs_lowerboundary,9)-1) = 2^15;
imagesc(csf_mask)
csfrecon = pax.reconstruction(csf_mask);
analyze_savesimpletif(csfrecon,'csf_boundary.tif')

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

figure(fig_bvkgph)
imagesc(hmap_x,hmap_y, bvkymograph)
colormap(clee.red_g)
set(gca, 'YDir', 'normal','Units','inches','Position',axsize)
clim(gca,[0 0.7])
hrframe = getframe(gca);
hrimg = frame2im(hrframe);
figure(fig_pvskgph)
imagesc(hmap_x,hmap_y, csfkymograph)
colormap(clee.green_g)
set(gca, 'YDir', 'normal','Units','inches','Position',axsize)
clim(gca,[0.1 1.3])
hgframe = getframe(gca);
hgimg = frame2im(hgframe);

figure(fig_bvb)
imagesc(hmap_x,hmap_y, bv_mask(:,851:end))
colormap(clee.cyan_g)
set(gca, 'YDir', 'normal','Units','inches','Position',axsize)
hcyanframe = getframe(gca);
clim(gca,[0.1 0.7])
hcyaimg = frame2im(hcyanframe);

figure(fig_pvsb)
imagesc(hmap_x,hmap_y, csf_mask(:,851:end))
colormap(clee.magenta_g)
set(gca, 'YDir', 'normal','Units','inches','Position',axsize)
hmagframe = getframe(gca);
clim(gca,[0.1 0.7])
hmagimg = frame2im(hmagframe);
%

figure()
imshow(hgimg+hmagimg+hcyaimg)

clim(gca,[0 0.7])

%%

bvmidkymograph = maskrefine(round(pax.idx.bv_upperboundary/2+pax.idx.bv_lowerboundary/2),pax);


axsize = [1 1 8 3];
figsize = [1 1 12 5];
fig_bvkgph  = figure(Name='fig_bvkgph',NumberTitle='off');
set(gcf,'Units','inches','Position',figsize); 

figure(fig_bvkgph)
imagesc(hmap_x,hmap_y, bvkymograph)
colormap(clee.red_g)
set(gca, 'YDir', 'normal','Units','inches','Position',axsize)
clim(gca,[0 0.7])
hrframe = getframe(gca);
hrimg = frame2im(hrframe);

fig_bvmax  = figure(Name='fig_bvmax',NumberTitle='off');
set(gcf,'Units','inches','Position',figsize); 

figure(fig_bvmax)
imagesc(hmap_x,hmap_y, bvmidkymograph(:,851:end))
colormap(clee.yellow_g)
set(gca, 'YDir', 'normal','Units','inches','Position',axsize)
hgframe = getframe(gca);
hbvmaximg = frame2im(hgframe);
%%
function output = maskrefine(idx,pax)
output = zeros(size(pax.mask.rowidx));
output(pax.mask.rowidx == medfilt1(idx,9)) = 2^16-1;
output(pax.mask.rowidx == medfilt1(idx,9)+1) = 2^15;
output(pax.mask.rowidx == medfilt1(idx,9)-1) = 2^15;
end