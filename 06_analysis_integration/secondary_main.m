clc, clear
analysis_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis';
experiment_folder = 'G:\tmp\00_igkl';
%%

[~, tmp.exp_name] = fileparts(experiment_folder);
tmp.save_exppath = fullfile(analysis_path,tmp.exp_name);
tmp.dirtable_dir = fullfile(tmp.save_exppath,strcat(tmp.exp_name,'_dirtable.xlsx'));

tmp.opts = detectImportOptions(tmp.dirtable_dir, 'Sheet', 'reference');
ref_table = readtable(tmp.dirtable_dir, tmp.opts);

% Filter for Sleep Sessions
tmp.rows = ref_table.LineFWHM == "paxfwhm.mat";
fwhm_table = ref_table(tmp.rows, :);
tmp.rows = fwhm_table.SessionType == "sleep";
sleep_fwhm_table = fwhm_table(tmp.rows,:);

%% 1. Aggregate FWHM Data
% Helper function located in 00_mapdirectorystruct
% Ensure the folder is on the path (if not already)
addpath(fullfile(fileparts(which('map_directory.m')), '00_mapdirectorystruct'));

disp('Aggregating FWHM Data...');
master_fwhm_table = aggregate_fwhm(sleep_fwhm_table);

disp('Aggregation Complete. Master Table Preview:');
head(master_fwhm_table)


%%
