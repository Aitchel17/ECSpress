%% map all the subdirectories in mouse folder
clc,clear,clean_editor
save_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis';
experiment_folder = 'G:\tmp\00_igkl';

%% create primary analysis folder in primary analysis save folder

[~, exp_name] = fileparts(experiment_folder);
save_exppath = fullfile(save_path,exp_name);
if ~isfolder(save_exppath)
    mkdir(save_exppath)
end
dirtable_dir = fullfile(save_exppath,strcat(exp_name,'_dirtable.xlsx'));
%
primary_map = {
    'RadonResult', 'radon_result.mat';
    'RoiList', 'roilist.mat';
    'PolarCluster', 'polarcluster.mat';
    'PaxFWHM', 'paxfwhm.mat'
    };

peripheral_map = {
    'AnalogAnalysis', 'analysis_analog.mat';
    'BehaviorAnalysis', 'analysis_camera.mat';
    'SleepScore', 'sleep_score.mat'
    };

stateanalysis_map = {'paxfwhm_state', 'paxfwhm_state.mat'};

% Go down to exp_dir which contains mouse folder, ex.hql071 ...
dirstruct_table = mapdirstruct(experiment_folder, primary_map,peripheral_map,stateanalysis_map);
write_dirtable(dirstruct_table, dirtable_dir);
%% Read reference sheet (Example)
opts = detectImportOptions(dirtable_dir, 'Sheet', 'reference');
ref_table = readtable(dirtable_dir, opts);
% Read state analysis
rows = ref_table.LineFWHM == "paxfwhm.mat";

fwhm_table = ref_table(rows, :);
%
rows = fwhm_table.SessionType == "sleep";

sleep_fwhm_table = fwhm_table(rows,:);