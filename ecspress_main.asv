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
base_path = 'G:\tmp\00_igkl\hql073\250626_hql073_whisker\HQL073_whisker250623_007';
directories = manage_directories(base_path);

%% 1. Load data
mdfExtract = load_mdfextract(directories.load_dir);
% 2. Twophoton data FPS matching & preprocessing
twophoton_processed = twophoton_preprocess(mdfExtract);
%% 3. Load processed data & ROIlist generation
[roilist, pax_fwhm, polarcluster, radon_analysis] = initialize_analysis_workspace(directories);

%% 4.1 FWHM Analysis
% 4.1.1 FWHM analysis - ROI Setupw
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
%%
pax_fwhm.clean_outlier(true)
pax_fwhm.getdiameter;
pax_fwhm.getdisplacement;
pax_fwhm.save2disk(directories.primary_analysis);
roilist.save2disk
%% 4.1.4 FWHM analysis figure generation
analysis_pax_makefig(pax_fwhm, twophoton_processed.t_axis, twophoton_processed.pixel2um, directories.figures_fwhm);


%% 4.2 Cluster polar analysis
%% 5.1 Make cluster
polarcluster = analysis_clusterpolar(pax_fwhm, twophoton_processed, directories.primary_analysis);
%% 5.2 Make cluster figure
analysis_clusterpolar_makefig(polarcluster, roilist, pax_fwhm, twophoton_processed.t_axis, twophoton_processed.pixel2um, directories.figures_polarcluster);
%% 5.3 Manual Contour Correction
polarcluster = analysis_clusterpolar_contour(polarcluster, roilist);
%% 5.4 Polar Plot of Contours
analysis_clusterpolar_polarplot(polarcluster, roilist, directories.figures_polarcluster);
roilist.save2disk();
save(fullfile(directories.primary_analysis, 'polarcluster.mat'), "polarcluster");


%% 6. Dynamic time warping based analysis
%% 7. PIV analysis

%% 8. Radon Analysis (Only capable for clear images without debris around artery)
roilist.addormodifyroi(twophoton_processed.ch1,'radon','rectangle');
radon_analysis = analysis_radon(twophoton_processed, roilist, 'ch1');

%% 8.2 Radon figures
analysis_radon_makefig(radon_analysis, twophoton_processed.t_axis, directories.figures_radon, twophoton_processed.pixel2um);
%%
radon_analysis.save2disk(directories.primary_analysis);
roilist.save2disk
%%
util_checkstack(radon_analysis.radon_result.events(4).irtd)

%% 7. ROI Setup
setup_rois(roilist, twophoton_processed);
%% 7.1 ROI Setup verification figures (Manual ROIs)
setup_rois_makefig(roilist, twophoton_processed.ch1, twophoton_processed.ch2, directories.figures_roi);
