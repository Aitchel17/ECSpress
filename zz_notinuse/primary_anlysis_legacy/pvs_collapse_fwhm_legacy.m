clc,clear
% loading file and info
%%
vivo_50um1 = mdfExtractLoader();
loaded_data = util_load_primarydataset(vivo_50um1.info.analyzefolder);
primary_analog = add_analog(loaded_data,vivo_50um1);
roilist = add_roilist(loaded_data);
clee = color_lee;
%%
primaryprocess_fps = 3; % fps
% load imagestack
vivo_50um1.stackch1 = vivo_50um1.loadstack("ch1");
vivo_50um1.stackch2 = vivo_50um1.loadstack("ch2");
% original time axis
tmp.fps = str2double(vivo_50um1.info.savefps); % mdfExtractor saving fps
tmp.frameavg = floor(tmp.fps)/primaryprocess_fps;
%%
tmp.loadstart = str2double(vivo_50um1.info.loadstart)/tmp.fps; % load start frame --> load start seconds
tmp.taxis = linspace(tmp.loadstart,tmp.loadstart+size(vivo_50um1.stackch1,3)/tmp.fps,size(vivo_50um1.stackch1,3));
% Preprocessing for fwhm based detection, 3 frame averaging, 5 frame median filter
[tmp.gpstack, tmp.ntailfr] =  analyze_grouproject(vivo_50um1.stackch2,tmp.frameavg,"mean");
tmp.medcsf = medfilt3(tmp.gpstack,[1 1 5]);
gausscsf= imgaussfilt(tmp.medcsf,1);
tmp.gpstack = analyze_grouproject(vivo_50um1.stackch1,tmp.frameavg,"mean");
tmp.medbv = medfilt3(tmp.gpstack,[1 1 5]);
gaussbv= imgaussfilt(tmp.medbv,1);
% modified time axis
tmp.modified_taxis = linspace(tmp.taxis(1),tmp.taxis(end)-tmp.ntailfr/tmp.fps,size(gaussbv,3));

%% Create roi class
roilist = roi(gausscsf,'thick','line');
%%

roilist  = roilist.copyroi('thick','thin');
%%
roilist = roilist.modifyroi(gausscsf, 'thin');
%%
roilist = roilist.modifyroi(gausscsf, 'thick');

%%
get_bvoutter(gausscsf)
%%
roilist = loaded_data.roilist;
roilist = roilist.modifyroi(gausscsf,'pax');
roilist = roilist.modifyroi(gaussbv,'pax');
roilist = roilist.addroi(gaussbv, 'radon','polygon');
%% to kymograph
rc_bi30_thick = analyze_affine_rotate(gausscsf,roilist.vertices.thick(1:2,:), roilist.vertices.thick(3,1));
kgph_bi30_thick = squeeze(sum(rc_bi30_thick,1));
%
rc_bi30_thin = analyze_affine_rotate(gausscsf,roilist.vertices.thin(1:2,:), roilist.vertices.thin(3,1));
kgph_bi30_thin = squeeze(sum(rc_bi30_thin,1));
%
fwhm_thin = get_bvoutter(kgph_bi30_thin);
fwhm_thick = get_bvoutter(kgph_bi30_thick);
%
figure()
imagesc(kgph_bi30_thin)
colormap()
hold on
plot(fwhm_thin.bv_lowerboundary,'r')
plot(fwhm_thin.bv_upperboundary,'r')

figure()
imagesc(kgph_bi30_thick)
hold on
plot(fwhm_thick.bv_lowerboundary,'r')
plot(fwhm_thick.bv_upperboundary,'r')
%%

figure()
plot(fwhm_thick.bv_lowerboundary-fwhm_thick.bv_upperboundary,'r')
%
hold on
plot(fwhm_thin.bv_lowerboundary-fwhm_thin.bv_upperboundary,'g')


%%
figure()
imagesc(kgph_bi30_thick)
colormap(clee.gray_g)
ylim([0 50])
set(gcf,'Position',[100 100 1100 300])
%%
util_checkstack(gausscsf)
%%
figure()
imagesc(kgph_bi30_thin)
colormap(clee.gray_g)
ylim([0 50])
set(gcf,'Position',[100 100 1100 300])
%%
maskbv = roi_applymask(gaussbv,roilist.mask.radon);
%%
roilist.showroi('thick')


%%
radon = analyze_radon();


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

% % Simulated kymograph data
% kymograph = pax_fwhm.kymograph.kgph_csf;
% kymograph = (kymograph-min(kymograph,[],1))./max(kymograph,[],1);
%
% % Display kymograph using imagesc with colormap
% imagesc(ax, kymograph);
% colormap(ax, parula);  % Or your inferno colormap
% hold(ax, 'on');
hold on
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
line_fwhms = struct();
line_fwhms.pax = pax_fwhm;
save(fullfile(savepath,'line_fwhms.mat'),'line_fwhms')
save(fullfile(savepath,'analog.mat'),'primary_analog')
save(fullfile(savepath,'roilist.mat'),'roilist')
roilist.showroi('pax')
h2frame = getframe(gca);
h2img = frame2im(h2frame);
%%
tmp.reconkymomask = pax.mask.pvs_downline+pax.mask.pvs_upline;
tmp.csfrecon = pax.reconstruction(tmp.reconkymomask)*65565;
analyze_savesimpletif(tmp.csfrecon,'csf_boundary.tif')
%%
analyze_savesimpletif(gausscsf,'csf.tif')
analyze_savesimpletif(gaussbv,'bv.tif')

%%


figure()
imshow(double(pax.rotatecrop.rc_bv(:,:,100)))



%%
util_checkstack(double(pax.rotatecrop.rc_bv(:,:,:)))

%%
analyze_savesimpletif(gaussbv,'bv.tif')
%%
analyze_savesimpletif(gausscsf,'csf.tif')


%%

%% make kymograph (hql080 002)
csfkymograph = pax.kymograph.kgph_normcsf(:,851:end);
bvkymograph = pax.kymograph.kgph_normbv(:,851:end);
hmap_y = linspace(-size(bvkymograph,1)*scale/2,size(bvkymograph,1)*scale/2,size(bvkymograph,1));
hmap_x = linspace(1,size(bvkymograph,2)/fps,size(bvkymograph,2));
%%



%%
set(gcf,'Units','inches','Position',[1 1 7 5]);

set(gca,'Units','inches','Position',[1 1 5.6 3]);
%%


%%
figure()
imshow(hrimg+hgimg)
%%
figure()
imshow(hrimg+hgimg+hmagimg+hcyaimg)
%%
fig_lineprofile    = figure(Name='fig_lineprofile',NumberTitle='off')
set(gcf,'Position',[100 100 1000 400])

%%
pax_fwhm.mask.roiimg = h2img;
%%

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
imadjust(stacks.irtd_radon(:,:,1),[0 1])
irtd = stacks.irtd_radon(:,:,1);
%%
window1=figure(Name='window1',NumberTitle='off');
figure(window1)
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