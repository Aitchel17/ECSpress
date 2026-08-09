% This is main ecspress analysis file
% This file is intended to show top most level of pipeline
% Must be simple as possible
% 20251229: Currently pipeline consist of following step
% 1. Load mdfExtracted files using mdfExtractLoader object
%   the Loader object was intended to retrieve info and load data using
%   method within it. If the modification related mdf output to be made, the mdfExtractLoader must modified
%   not at analysis level
% 2. Directory tree management
% Perhaps its better to seperate analysis output from the
% mdfExtracted files.. as the data set itself become to big...
% 3. ROI setup
% This is also currently upstream process
% May be its better to seperate this to other processing as ROI
% selection loop untill the feature detection algorithm to work.

% 4. Analysis
% Currently planned primary analysis outputs are
% 1. Line full width halfmaximum (position-time) 3,4,5,7

% Path setup
addpath(genpath(pwd));
close all
% Directory setup
% sessiondir = 'G:\tmp\01_igkltdt\hql104\260607_hql104_sleep\HQL104_sleep260607_006';
sessiondir ='G:\tmp\01_igkltdt\hql080\250619_hql080_sleep\HQL080_sleep250619_002';

% 1. Load data & 3. Load processed data (Integrated via ECSSession)
session = ECSSession(sessiondir);
session = session.load_primary_results;
%%


session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
% 2. Twophoton data FPS matching & preprocessing
% Note: twophoton_preprocess expects a struct with stackch1/2 and img_param.
% ECSSession object works here as it has these properties.
twophoton_processed = twophoton_preprocess(session);
% Load Primary Results (replaces initialize_analysis_workspace/load_primaryresult)

%% 4.1 FWHM Analysis
% 4.1.1 FWHM analysis - ROI Setupw
session.roilist.addormodifyroi(twophoton_processed.ch2,'pax','line',twophoton_processed.ch1);
%% 4.1.2 Initialize Analysis Object & Lumen Analysis
session.pax_fwhm = line_fwhm(session.roilist.getvertices('pax'));
session.pax_fwhm.param.fs = twophoton_processed.outfps;
session.pax_fwhm.t_axis = twophoton_processed.t_axis;
% Lumen (vessel) processing
bvch = 'ch2';                                   % BV/lumen channel (recordings may swap ch1/ch2)
session.pax_fwhm.param.channel_lumen = bvch;    % single source: data + record both from bvch
session.pax_fwhm.addkymograph("lumen", twophoton_processed.(bvch),"max")
session.pax_fwhm.kymograph_afterprocess('lumen',[1 5])
session.pax_fwhm.fwhm("lumen");

%% 4.1.3 PVS Analysis
% PVS processing (using channel 2)
csfch = 'ch1';                                  % CSF/PVS channel
session.pax_fwhm.param.channel_pvs = csfch;     % single source: data + record both from csfch
session.pax_fwhm.addkymograph("pvs", twophoton_processed.(csfch),"mean")
session.pax_fwhm.kymograph_afterprocess('pvs',[1 5])
%%
session.pax_fwhm.pvsanalysis_inverted();
%%
session.pax_fwhm.clean_outlier(true)
session.pax_fwhm.getdiameter;
session.pax_fwhm.getdisplacement;
session.pax_fwhm.save2disk('paxfwhm',session.dir_struct.primary_analysis);
session.roilist.save2disk(session.dir_struct.primary_analysis)

%% 4.1.5 Manual ROI setup  (needs pax_fwhm; everything downstream reads these)
setup_rois(session.roilist, twophoton_processed, session.pax_fwhm, session.dir_struct.primary_analysis);
setup_rois_makefig(session.roilist, twophoton_processed, session.pax_fwhm, session.dir_struct.figures_roi);

%% 4.1.4 FWHM analysis figure generation
close all
pax_fig = analysis_pax_makefig(session.pax_fwhm, twophoton_processed.t_axis,...
    twophoton_processed.pixel2um, session.dir_struct.figures_fwhm);
%% Opt Output for FWHM boundary overlay video
session.pax_fwhm.reconstruction_overlay(twophoton_processed.ch1,twophoton_processed.ch2,...
    "SavePath",fullfile(session.dir_struct.primary_analysis,"overlay.tif"),"BlankBoundary",false);

%%
if isfile(fullfile(sessiondir,"peripheral/sleep_score.mat"))
    disp('Loading sleepscore file')
    sleepscore = load(fullfile(sessiondir,"peripheral/sleep_score.mat"));
else
    disp('no sleep_socre.mat')
end
fwhm = session.pax_fwhm.thickness;
%%
clee = color_lee;
fig = figure();
sgolayax = axes(fig);

cla(sgolayax)
plot(sgolayax,session.pax_fwhm.t_axis,sgolayfilt(fwhm.pvschanges_total*session.img_param.pixel2um,3,11),"Color",clee.clist.darkgreen,'LineWidth',1)
hold on
plot(sgolayax,session.pax_fwhm.t_axis,sgolayfilt(fwhm.epschanges*session.img_param.pixel2um,3,11),"Color",clee.clist.gold,'LineWidth',1)

plot(sgolayax,session.pax_fwhm.t_axis,sgolayfilt(fwhm.bvchanges*session.img_param.pixel2um,3,11),"Color",'red','LineWidth',1)
plot_sleep_patches(sgolayax,sleepscore)
%%
fig = figure();
sgolayax = axes(fig);

cla(sgolayax)
plot(sgolayax,session.pax_fwhm.t_axis,fwhm.epschanges*session.img_param.pixel2um,"Color",'b','LineWidth',2)
hold on
plot(sgolayax,session.pax_fwhm.t_axis,fwhm.bvchanges*session.img_param.pixel2um,"Color",'red','LineWidth',2)
plot_sleep_patches(sgolayax,sleepscore)
%%

%%
tmp.data= ans;
%%
figure()
imagesc(squeeze(tmp.data(:,:,1,2:4)))


%% 4.2 Cluster polar analysis
%% 5.1 Make cluster
session.polarcluster = analysis_clusterpolar(session.pax_fwhm, twophoton_processed, session.dir_struct.primary_analysis);
%% 5.2 Make cluster figure
analysis_cluster_makefig(session.polarcluster, session.roilist, session.pax_fwhm, twophoton_processed.t_axis, twophoton_processed.pixel2um, session.dir_struct.figures_polarcluster);
%% 5.3 Manual Contour Correction
session.polarcluster = analysis_clusterpolar_contour(session.polarcluster, session.roilist);
%% 5.4 Polar transform of contours (cartesian -> polar profiles)
session.polarcluster = analysis_clusterpolar_polarplot(session.polarcluster, session.roilist);
%% 5.4b Polar figure
analysis_polar_makefig(session.polarcluster, session.dir_struct.figures_polarcluster);
session.roilist.save2disk(session.dir_struct.primary_analysis);
% 5.5 save polar cluster
polarcluster = session.polarcluster;
save(fullfile(session.dir_struct.primary_analysis, 'polarcluster.mat'), "polarcluster");
%%
figure()

imagesc(polarcluster.dilated_medianimg)
%%
axis image
%%
colormap(clee.gradient.inferno)
%%
imcontrast
%% 6. Dynamic time warping based analysis
%% 7. PIV analysis

%% 8. Radon Analysis (Only capable for clear images without debris around artery)
bvch = session.pax_fwhm.param.channel_lumen;    % radon follows FWHM's BV/vessel channel
session.roilist.addormodifyroi(twophoton_processed.(bvch),'radon','rectangle');
%%
session.radon_analysis = analysis_radon(twophoton_processed, session.roilist, session.pax_fwhm.param.channel_lumen);



%% 8.2 Radon figures
analysis_radon_makefig(session.radon_analysis.radon_result, twophoton_processed.t_axis, session.dir_struct.figures_radon, twophoton_processed.pixel2um);
%%
session.radon_analysis.get_elipsoid
%%
figure()
%%




%%
radonresult = session.radon_analysis.radon_result;
tmp.normt_diameter =radonresult.diameter./mean(radonresult.diameter,2) ;
plot(tmp.normt_diameter(:,100))

%%
tmp.normd_diameter =radonresult.diameter./mean(radonresult.diameter,1) ;
%%
figure()
imagesc(tmp.normd_diameter)
%%
plot(tmp.normd_diameter(:,100))
hold on
plot(tmp.normd_diameter(:,1200))

%%
figure()
plot(tmp.normd_diameter(:,101))

hold on
plot(tmp.normd_diameter(:,880))

plot(tmp.normd_diameter(:,1204))
plot(tmp.normd_diameter(:,1242))


%%
figure()
plot(radonresult.elipsoidfit.AI_1minus)

%%

figure()
plot(radonresult.elipsoidfit.resnorm)






%%
session.radon_analysis.save2disk(session.dir_struct.primary_analysis);
session.roilist.save2disk
%%
util_checkstack(session.radon_analysis.radon_result.events(4).irtd)

