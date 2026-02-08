function setup_workspace()
% SETUP_WORKSPACE Links dependencies from 01_ecspress and 02_othersignal
% Run this script before running analysis scripts in this repository.

% Get the root program directory (assuming specific depth)
% This file is in: g:\03_program\03_secondary_analysis
current_dir = fileparts(mfilename('fullpath'));
[root_program_dir, ~, ~] = fileparts(current_dir);

% Define Dependency Paths
path_ecspress = fullfile(root_program_dir, '01_ecspress');
path_othersignal = fullfile(root_program_dir, '02_othersignal');

% 1. Link ECSPRESS (Functions and Classes)
if isfolder(path_ecspress)
    addpath(genpath(fullfile(path_ecspress, 'functions')));
    addpath(genpath(fullfile(path_ecspress, 'class')));
    % Add primary analysis functions if needed?
    % Generally good to be explicit, but 'functions' is the main one.
    disp(['Linked: ' path_ecspress]);
else
    warning(['Not Found: ' path_ecspress]);
end

% 2. Link OTHERSIGNAL
if isfolder(path_othersignal)
    % Assuming 02_othersignal has a functions folder or similar structure
    % If structure is unknown, just add the root or known subfolders
    addpath(genpath(path_othersignal));
    disp(['Linked: ' path_othersignal]);
else
    warning(['Not Found: ' path_othersignal]);
end

% 3. Link This Repository (Recursive)
addpath(genpath(current_dir));
disp(['Linked Internal: ' current_dir]);

disp('Workspace Setup Complete.');
end
