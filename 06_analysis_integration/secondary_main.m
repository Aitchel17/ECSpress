clc, clear, clean_editor
masterDirTable_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl\00_igkl_dirtable.xlsx';

%% 0. Initialize TableManager
% Initialize Manager (Loads Excel automatically)
mtable_FWHMsleep = tableManager(masterDirTable_path);
% 1. Table filtering
mtable_FWHMsleep.filter_refTable("Primary_paxFWHM","paxfwhm");
mtable_FWHMsleep.filter_refTable("State_PaxFWHM","paxfwhm");
mtable_FWHMsleep.parseParams();
% 2. Aggregate Data to Vessels
% Loads FWHM data from .mat files and populates Vessel objects
mtable_FWHMsleep.aggregateData("State_PaxFWHM");
nestname_arr = ["state_summary","transition","band_decomposition","powerdensity","peak_trough"];
mtable_FWHMsleep.addnest2subtable(nestname_arr)
% save (except subtables, auto reconstruct)
mtable_FWHMsleep.save2disk("mtable_FWHMsleep.mat")



