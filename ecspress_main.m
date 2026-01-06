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
% 1. Line full width halfmaximum (position-time)

% Path setup
addpath(genpath(pwd));

% Directory setup
base_path = 'G:\tmp\00_igkl\hql090\251020_hql090_whiskerb\HQL090_whiskerb251020_002';
directories = manage_directories(base_path);
%%

%% 1. Load data
mdfExtract = load_mdfextract(directories.load_dir);
% 2. Twophoton data FPS matching & preprocessing
twophoton_processed = twophoton_preprocess(mdfExtract);
% 3. ROIlist generation
roilist = roi_handle(fullfile(directories.primary_analysis,"roilist.mat"));

%% 4.1 FWHM Analysis
% 4.1.1 FWHM analysis - ROI Setup
roilist.addormodifyroi(twophoton_processed.ch2,'pax','line');
%% 4.1.2 Initialize Analysis Object & Lumen Analysis
pax_fwhm = line_fwhm(roilist.getvertices('pax'));

% Lumen (vessel) processing
pax_fwhm.addkymograph("lumen", twophoton_processed.ch1,"max")
pax_fwhm.kymograph_afterprocess('lumen',[1 3])
pax_fwhm.fwhm("lumen");
%% 4.1.3 PVS Analysis
% PVS processing (using channel 2)
pax_fwhm.addkymograph("pvs", twophoton_processed.ch2,"median")
pax_fwhm.kymograph_afterprocess('pvs',[3 5])
pax_fwhm.pvsanalysis();
pax_fwhm.clean_outlier(true)
pax_fwhm.getdiameter;
pax_fwhm.getdisplacement;
%% 4.1.4 FWHM analysis figure generation
analysis_pax_makefig(pax_fwhm, twophoton_processed.t_axis, twophoton_processed.pixel2um, directories.figures_fwhm);
roilist.save2disk

%% 4.2 Cluster polar analysis
% 5.1 Make cluster
pax_cluster = analysis_clusterpolar(pax_fwhm, twophoton_processed, directories.primary_analysis);
%% 5.2 Make cluster figure
analysis_clusterpolar_makefig(pax_cluster, pax_fwhm, twophoton_processed.t_axis, twophoton_processed.pixel2um, directories.figures_cluster);
%% 5.3 Manual Contour Correction
pax_cluster = analysis_clusterpolar_contour(pax_cluster, roilist);
%% 5.4 Polar Plot of Contours
analysis_clusterpolar_polarplot(pax_cluster, roilist, directories.figures_cluster);
%%
figure()
imagesc(pax_cluster.manual_roi(1).EdgeMask+pax_cluster.manual_roi(2).EdgeMask)
%%
figure()
imagesc(pax_cluster.manual_roi(3).EdgeMask+pax_cluster.manual_roi(4).EdgeMask)
%%
figure()
imagesc(pax_cluster.manual_roi(2).EdgeMask+pax_cluster.manual_roi(4).EdgeMask)
%%
roilist.save2disk();
save(fullfile(directories.primary_analysis, 'pax_cluster.mat'), "pax_cluster");
%%
figure()
imagesc(pax_cluster.cluster_bvcsf_constdil(:,:,1,1))
axis image
%%
imcontrast

%%
figure()
%%
imagesc(pax_cluster.manual_roi(1).SelectionMask)
axis image
%%

% 5.2 Polar analysis (with findimage edge, watershed algorithm)
%% 6. Dynamic time warping based analysis
%% 7. PIV analysis

%% 8. Radon Analysis (Only capable for clear images without debris around artery)
roilist.addormodifyroi(twophoton_processed.ch1,'radon','rectangle');
radon_analysis = analysis_radon(twophoton_processed, roilist, 'ch1');

%% 8.2 Radon figures
analysis_radon_makefig(radon_analysis, twophoton_processed.t_axis, directories.figures_radon);
%%
radon_analysis.save2disk(directories.primary_analysis);
roilist.save2disk

%% 7. ROI Setup
setup_rois(roilist, twophoton_processed);
