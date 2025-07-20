clc,clear

% loading file and info
vivo_50um1 = mdfExtractLoader();
savepath = fullfile(vivo_50um1.info.analyzefolder,'primary_analysis');
if ~exist(savepath, 'dir')
    mkdir(savepath);
end 

% load and process analog channel
vivo_50um1.analog = vivo_50um1.loadanalog;
primary_analog = analysis_analog(vivo_50um1.analog.info,vivo_50um1.analog.data);
primary_analog.airtable = primary_analog.get_airtable('raw_Air_puff1');
primary_analog.ecogspectrum = primary_analog.get_ecogspectrum('raw_ECoG');
%% load imagestack
vivo_50um1.stackch1 = vivo_50um1.loadstack("ch1");
vivo_50um1.stackch2 = vivo_50um1.loadstack("ch2");
% original time axis
tmp.fps = str2double(vivo_50um1.info.savefps); % mdfExtractor saving frame
tmp.loadstart = str2double(vivo_50um1.info.loadstart)/tmp.fps; % load start frame --> load start seconds
tmp.taxis = linspace(tmp.loadstart,tmp.loadstart+size(vivo_50um1.stackch1,3)/tmp.fps,size(vivo_50um1.stackch1,3));
% Preprocessing for fwhm based detection, 3 frame averaging, 5 frame median filter
[tmp.gpstack, tmp.ntailfr] =  analyze_grouproject(vivo_50um1.stackch2,3,"mean");
tmp.medcsf = medfilt3(tmp.gpstack,[1 1 5]);
gausscsf= imgaussfilt(tmp.medcsf,1);
tmp.gpstack = analyze_grouproject(vivo_50um1.stackch1,3,"mean");
tmp.medbv = medfilt3(tmp.gpstack,[1 1 5]);
gaussbv= imgaussfilt(tmp.medbv,1);
% modified time axis
tmp.modified_taxis = linspace(tmp.taxis(1),tmp.taxis(end)-tmp.ntailfr/tmp.fps,size(gaussbv,3));

%% Create roi class
roilist = roi(gausscsf,'pax','line');
%%
roilist = roilist.modifyroi(gausscsf,'pax');
%%
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
kymograph = pax_fwhm.kymograph.kgph_csf;
kymograph = (kymograph-min(kymograph,[],1))./max(kymograph,[],1);

% Display kymograph using imagesc with colormap
imagesc(ax, kymograph);
colormap(ax, parula);  % Or your inferno colormap
hold(ax, 'on');

%% Draw horizontal line on top
fig = figure;
ax = axes('Parent', fig);
line_y = pax_fwhm.idx.pvs_lowerboundary;
line_x = linspace(1,size(line_y,2),size(line_y,2));
plot(ax, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
line_y = pax_fwhm.idx.pvs_upperboundary;
plot(ax, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
%%
line_y = pax_fwhm.idx.bv_lowerboundary;
plot(ax, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
line_y = pax_fwhm.idx.bv_upperboundary;
plot(ax, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
%%
fig = figure;
ax = axes('Parent', fig);

% Simulated kymograph data
kymograph = pax_fwhm.kymograph.kgph_bv;
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
line_fwhms = struct();
line_fwhms.pax = pax_fwhm;
save(fullfile(savepath,'line_fwhms.mat'),'line_fwhms')
save(fullfile(savepath,'analog.mat'),'primary_analog')



%%
sliceViewer(csf_mask+bv_mask)
%%

util_checkstack(gaussbv)
%%
sliceViewer(~(csf_mask+bv_mask).*gausscsf)


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



