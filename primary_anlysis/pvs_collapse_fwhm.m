clc,clear
%%
% loading file and info
vivo_50um1 = mdfExtractLoader();
afterprocessing_avgproj = 3;
%%
if isfile(fullfile(vivo_50um1.info.analyzefolder,'primary_analysis/line_fwhms.mat'))
    tmp.loadstruct = load(fullfile(vivo_50um1.info.analyzefolder,'primary_analysis/line_fwhms.mat'));
end

%%
savepath = fullfile(vivo_50um1.info.analyzefolder,'primary_analysis');
if ~exist(savepath, 'dir')
    mkdir(savepath);
end 
% load and process analog channel
vivo_50um1.analog = vivo_50um1.loadanalog;
%
primary_analog = analysis_analog(vivo_50um1.analog.info,vivo_50um1.analog.data);
primary_analog.airtable = primary_analog.get_airtable('raw_Air_puff1');
primary_analog.ecogspectrum = primary_analog.get_ecogspectrum('raw_ECoG');
% load imagestack
vivo_50um1.stackch1 = vivo_50um1.loadstack("ch1");
vivo_50um1.stackch2 = vivo_50um1.loadstack("ch2");
% original time axis
tmp.fps = str2double(vivo_50um1.info.savefps); % mdfExtractor saving frame
tmp.loadstart = str2double(vivo_50um1.info.loadstart)/tmp.fps; % load start frame --> load start seconds
tmp.taxis = linspace(tmp.loadstart,tmp.loadstart+size(vivo_50um1.stackch1,3)/tmp.fps,size(vivo_50um1.stackch1,3));
% Preprocessing for fwhm based detection, 3 frame averaging, 5 frame median filter
[tmp.gpstack, tmp.ntailfr] =  analyze_grouproject(vivo_50um1.stackch2,afterprocessing_avgproj,"mean");
tmp.medcsf = medfilt3(tmp.gpstack,[1 1 5]);
gausscsf= imgaussfilt(tmp.medcsf,1);
tmp.gpstack = analyze_grouproject(vivo_50um1.stackch1,afterprocessing_avgproj,"mean");
tmp.medbv = medfilt3(tmp.gpstack,[1 1 5]);
gaussbv= imgaussfilt(tmp.medbv,1);
% modified time axis
tmp.modified_taxis = linspace(tmp.taxis(1),tmp.taxis(end)-tmp.ntailfr/tmp.fps,size(gaussbv,3));

%% Create roi class
roilist = roi(gausscsf,'pax','line');
roilist = roilist.modifyroi(gausscsf,'pax');
roilist = roilist.modifyroi(gaussbv,'pax');

%% fwhm analysis class test
pax_fwhm = line_fwhm(gaussbv,gausscsf,roilist.vertices.pax);
pax_fwhm = pax_fwhm.bvanalysis(0.5,10);
pax_fwhm = pax_fwhm.csfanalysis(0.5,10);
bv_mask = pax_fwhm.reconstruction(pax_fwhm.mask.bv_upline+pax_fwhm.mask.bv_downline);
csf_mask = pax_fwhm.reconstruction(pax_fwhm.mask.pvs_upline+pax_fwhm.mask.pvs_downline);
pax_fwhm.t_axis = tmp.modified_taxis;

% kymograph with overlay
fig = figure;
ax = axes('Parent', fig);

% Simulated kymograph data
kymograph = pax_fwhm.kymograph.kgph_bv;
kymograph = (kymograph-min(kymograph,[],1))./max(kymograph,[],1);

% Display kymograph using imagesc with colormap
imagesc(ax, kymograph);
colormap(ax, parula);  % Or your inferno colormap
hold(ax, 'on');


%%
line_y = pax_fwhm.idx.bv_lowerboundary;
plot(ax, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
line_y = pax_fwhm.idx.bv_upperboundary;
plot(ax, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
%%
fig = figure;
ax = axes('Parent', fig);

% Simulated kymograph data
kymograph = pax_fwhm.kymograph.kgph_csf;
kymograph = (kymograph-min(kymograph,[],1))./max(kymograph,[],1);

% Display kymograph using imagesc with colormap
imagesc(ax, kymograph);
colormap(ax, parula);  % Or your inferno colormap
hold(ax, 'on');

% Draw horizontal line on top
line_y = pax_fwhm.idx.pvs_lowerboundary;
line_x = linspace(1,size(line_y,2),size(line_y,2));
plot(ax, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
line_y = pax_fwhm.idx.pvs_upperboundary;
plot(ax, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
line_y = pax_fwhm.idx.bv_lowerboundary;
plot(ax, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
line_y = pax_fwhm.idx.bv_upperboundary;
plot(ax, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);

%% saving
roilist.showroi('pax')
h2frame = getframe(gca);
h2img = frame2im(h2frame);

%% visualization testment
pax = tmp.loadstruct.line_fwhms.pax;
%%
tmp.reconkymomask = pax.mask.bv_upline+pax.mask.bv_downline; %+ pax.mask.bv_midline;
tmp.bvrecon = pax.reconstruction(tmp.reconkymomask)*65565;
analyze_savesimpletif(tmp.bvrecon,'bv_boundary.tif')
%%
tmp.reconkymomask = pax.mask.pvs_downline+pax.mask.pvs_upline;
tmp.csfrecon = pax.reconstruction(tmp.reconkymomask)*65565;
analyze_savesimpletif(tmp.csfrecon,'csf_boundary.tif')
%%
analyze_savesimpletif(gausscsf,'csf.tif')
analyze_savesimpletif(gaussbv,'bv.tif')


%%
figure()
%%
pax = tmp.loadstruct.line_fwhms.pax;
%%

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

%%

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

%%

util_checkstack(csfrecon)
%%
analyze_savesimpletif(gaussbv,'bv.tif')
%%
analyze_savesimpletif(gausscsf,'csf.tif')


%%
fps = str2double(vivo_50um1.info.savefps)/afterprocessing_avgproj;
clee = color_lee;
%%

%% make kymograph (hql080 002)
csfkymograph = pax.kymograph.kgph_normcsf(:,851:end);
bvkymograph = pax.kymograph.kgph_normbv(:,851:end);
hmap_x = linspace(0,size(bvkymograph,1)*0.57,size(bvkymograph,1));
hmap_y = linspace(1,size(bvkymograph,2)/fps,size(bvkymograph,2));
%%
figure()
imagesc(hmap_y,hmap_x, bvkymograph)
colormap(clee.red_g)
set(gca, 'YDir', 'normal')
hrframe = getframe(gca);
set(gcf,'Position',[100 100 1000 500])
hrimg = frame2im(hrframe);
%%
max(bvkymograph,[],"all")
%%
clim()

%%
figure()
imagesc(hmap_y,hmap_x, csfkymograph)
colormap(clee.green_g)
set(gca, 'YDir', 'normal')
hgframe = getframe(gca);
set(gcf,'Position',[100 100 1000 500])
hgimg = frame2im(hgframe);
%%
figure()
imagesc(hmap_y,hmap_x,bv_mask(:,851:end))
colormap(clee.cyan_g)
set(gca, 'YDir', 'normal')
hcyanframe = getframe(gca);
set(gcf,'Position',[100 100 1000 500])
hcyaimg = frame2im(hcyanframe);
%%
figure()
imagesc(hmap_y,hmap_x,csf_mask(:,851:end))
colormap(clee.magenta_g)
set(gca, 'YDir', 'normal')
hmagframe = getframe(gca);
set(gcf,'Position',[100 100 1000 500])
hmagimg = frame2im(hmagframe);

%%
figure()
imshow(hrimg+hgimg+hmagimg+hcyaimg)
%%


figure()
imshow()

%%
pax_fwhm.mask.roiimg = h2img;
%%
line_fwhms = struct();
line_fwhms.pax = pax_fwhm;
save(fullfile(savepath,'line_fwhms.mat'),'line_fwhms')
save(fullfile(savepath,'analog.mat'),'primary_analog')
save(fullfile(savepath,'roilist.mat'),'roilist')

%%
sliceViewer(csf_mask+csf_mask)
%%
util_checkstack(gaussbv)
%%
sliceViewer(~(csf_mask+csf_mask).*gausscsf)


%%
figure()
imagesc(double(pax_fwhm.mask.pvs_down)+double(pax_fwhm.mask.pvs_up))
merged = cat(4,box);


%%
figure()
plot(medfilt1(pax_fwhm.idx.pvs_upperboundary,5))
hold on
plot(medfilt1(pax_fwhm.idx.bv_upperboundary,5))


plot(medfilt1(pax_fwhm.idx.pvs_lowerboundary,5))
plot(medfilt1(pax_fwhm.idx.bv_lowerboundary,5))

%%
roilist = roilist.addroi(gaussbv,'radon','polygon');
%%
roilist.showroi('radon')
%%
stackfieldname ='radon';
%%
stack = roi_applymask(gaussbv,roilist.mask.radon);
stack(isnan(stack))= 0;
for i = 1:size(stack, 3)
    sli = medfilt2(stack(:, :, i), [3 3]);  % Apply 3x3 median filter (adjust size as needed)
    stack(:, :, i) = imgaussfilt(sli,1);
end
%%
radon_result = analyze_radon(stack); % do tirs
%%

util_checkstack(radon_result.irtd)
%%
radon_fwhm = radon_result.idx_downloc-radon_result.idx_uploc;
%%
rgb_img = zeros(size(radon_result.irtd(:,:,1),1),size(radon_result.irtd(:,:,1),2),3);
rgb_img(:,:,1) = mat2gray(radon_result.irtd(:,:,1))     >0.5;
rgb_img(:,:,2) = mat2gray(radon_result.irtd(:,:,500))   >0.5;
rgb_img(:,:,3) = mat2gray(radon_result.irtd(:,:,1))     >0.5;

figure()
imshow(rgb_img)

%%
figure(window1)
imagesc(radon_fwhm)
%%
figure()
plot(radon_fwhm(90,:))
hold on
plot(radon_fwhm(31,:))

%%
figure()
var(radon_fwhm,0,2)
plot(ans)
%%
[~, bottom_loc] = max(stacks.tirs_radon,[],1);
%%
stacks.tirs_radon(stacks.tirs_radon==0)=9999;

%%
[~, top_loc] = min(stacks.tirs_radon,[],1);

%%

bottom_loc = squeeze(bottom_loc);
top_loc = squeeze(top_loc);


%%
util_checkstack(stacks.irtd_radon)
%%
imadjust(stacks.irtd_radon(:,:,1),[0 1])
%%
irtd = stacks.irtd_radon(:,:,1);
%%

%%
window1=figure(Name='window1',NumberTitle='off');

figure(window1)
%%


%%
for i = 1:3000
irtd = stacks.irtd_radon(:,:,i);
img = mat2gray(stack(:,:,i));
resize_img = zeros(size(irtd));
resize_img(ceil(size(irtd,1)/2)+1-floor(size(img,1)/2):ceil(size(irtd,1)/2)+floor(size(img,1)/2),...
    ceil(size(irtd,1)/2)+1-floor(size(img,1)/2):ceil(size(irtd,1)/2)+floor(size(img,1)/2)) = img;
rgb_img = zeros(size(irtd,1),size(irtd,2),3);

rgb_img(:,:,1) = mat2gray(resize_img);
rgb_img(:,:,2) = mat2gray(irtd);
rgb_img(:,:,3) = mat2gray(resize_img);
imshow(rgb_img)
end

figure(window1)
imagesc(stacks.tirs_radon(:,:,1))

%%
% cat(4,irtd;resize_img;resize_img)
%%


