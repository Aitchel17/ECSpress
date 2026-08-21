clc, clear
masterDirTable_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl\00_igkl_dirtable.xlsx';
mtable_FWHMsleep = tableManager.load_recon(masterDirTable_path, "mtable_FWHMsleep.mat");
%%
transition_table = mtable_FWHMsleep.subTables.transition;
%% logic generation
% translogic = transition_table.state_name =="nr_trans";
vessellogic = contains(transition_table.VesselID, "PA", 'IgnoreCase', true);
depthlogic = transition_table.Depth <70;
combinedlogic = all([vessellogic, depthlogic],2);
%%
L1pa_transtable =   transition_table(combinedlogic,:);
data_col = "data";
numeric_cols = {'pre_mean','pre_median','pre_q1','pre_q3', 'post_mean', 'post_median', 'post_q1', 'post_q3'};
key_cols = ['MouseID', 'VesselID', "state_name"];
sessAve_transitionTable = average_table(L1pa_transtable,key_cols, numeric_cols,data_col);
%%
