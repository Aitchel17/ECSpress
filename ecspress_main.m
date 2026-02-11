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
sessiondir = 'G:\tmp\00_igkl\hql090\251012_hql090_sleep\HQL090_sleep251012_003';
% directories = manage_directories(base_path); % Removed, handled by ECSSession


%% 1. Load data & 3. Load processed data (Integrated via ECSSession)
session = ECSSession(sessiondir);
session = session.load_primary_results();
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
session.roilist.addormodifyroi(twophoton_processed.ch2,'pax','line');
%% 4.1.2 Initialize Analysis Object & Lumen Analysis
session.pax_fwhm = line_fwhm(session.roilist.getvertices('pax'));
session.pax_fwhm.param.fs = twophoton_processed.outfps;
session.pax_fwhm.t_axis = twophoton_processed.t_axis;
% Lumen (vessel) processing
session.pax_fwhm.addkymograph("lumen", twophoton_processed.ch1,"mean")
session.pax_fwhm.kymograph_afterprocess('lumen',[1 3])
session.pax_fwhm.fwhm("lumen");

%% 4.1.3 PVS Analysis
% PVS processing (using channel 2)
session.pax_fwhm.addkymograph("pvs", twophoton_processed.ch2,"mean")
session.pax_fwhm.kymograph_afterprocess('pvs',[1 5])
session.pax_fwhm.pvsanalysis();
%%
session.pax_fwhm.clean_outlier(true)
session.pax_fwhm.getdiameter;
session.pax_fwhm.getdisplacement;
session.pax_fwhm.save2disk('paxfwhm',session.dir_struct.primary_analysis);
session.roilist.save2disk(session.dir_struct.primary_analysis)
%% 4.1.4 FWHM analysis figure generation
analysis_pax_makefig(session.pax_fwhm, twophoton_processed.t_axis, twophoton_processed.pixel2um, session.dir_struct.figures_fwhm);


%% 4.2 Cluster polar analysis
%% 5.1 Make cluster
session.polarcluster = analysis_clusterpolar(session.pax_fwhm, twophoton_processed, session.dir_struct.primary_analysis);
%% 5.2 Make cluster figure
analysis_clusterpolar_makefig(session.polarcluster, session.roilist, session.pax_fwhm, twophoton_processed.t_axis, twophoton_processed.pixel2um, session.dir_struct.figures_polarcluster);
%% 5.3 Manual Contour Correction
session.polarcluster = analysis_clusterpolar_contour(session.polarcluster, session.roilist);
%% 5.4 Polar Plot of Contours
analysis_clusterpolar_polarplot(session.polarcluster, session.roilist, session.dir_struct.figures_polarcluster);
session.roilist.save2disk(session.dir_struct.primary_analysis);
% 5.5 save polar cluster
polarcluster = session.polarcluster;
save(fullfile(session.dir_struct.primary_analysis, 'polarcluster.mat'), "polarcluster");
%%








%% 6. Dynamic time warping based analysis
%% 7. PIV analysis

%% 8. Radon Analysis (Only capable for clear images without debris around artery)
session.roilist.addormodifyroi(twophoton_processed.ch1,'radon','rectangle');
session.radon_analysis = analysis_radon(twophoton_processed, session.roilist, 'ch1');

%% 8.2 Radon figures
analysis_radon_makefig(session.radon_analysis.radon_result, twophoton_processed.t_axis, session.dir_struct.figures_radon, twophoton_processed.pixel2um);
%%
session.radon_analysis.save2disk(session.dir_struct.primary_analysis);
session.roilist.save2disk
%%
util_checkstack(session.radon_analysis.radon_result.events(4).irtd)

%% 7. ROI Setup
setup_rois(session.roilist, twophoton_processed,session.dir_struct.primary_analysis);
%% 7.1 ROI Setup verification figures (Manual ROIs)
setup_rois_makefig(session.roilist, twophoton_processed.ch1, twophoton_processed.ch2, session.dir_struct.figures_roi);
