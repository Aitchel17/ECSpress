function directories = manage_directories(base_path)
% MANAGE_DIRECTORIES Sets up the directory structure for analysis.
%   directories = manage_directories(base_path) returns a struct with fields:
%       - load_dir: The base path provided.
%       - primary_analysis: Path for primary analysis results.
%       - figures: Path for saving figures (unique per run to preserve results).
%   It also ensures the output directories exist.

if nargin < 1 || isempty(base_path)
    error('Base path must be provided.');
end

directories.load_dir = base_path;
directories.primary_analysis = fullfile(base_path, 'primary_analysis');

% Ensure primary_analysis directory exists
if ~exist(directories.primary_analysis, 'dir')
    mkdir(directories.primary_analysis);
    disp(['Created directory: ', directories.primary_analysis]);
end

% Create unique figure directory to preserve results
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
directories.figures = fullfile(directories.primary_analysis, ['figures_', timestamp]);
if ~exist(directories.figures, 'dir')
    mkdir(directories.figures);
    disp(['Created figure directory: ', directories.figures]);
end

% Create subdirectories for FWHM and Radon figures
directories.figures_fwhm = fullfile(directories.figures, 'fwhm');
if ~exist(directories.figures_fwhm, 'dir')
    mkdir(directories.figures_fwhm);
end

directories.figures_radon = fullfile(directories.figures, 'radon_figures');
if ~exist(directories.figures_radon, 'dir')
    mkdir(directories.figures_radon);
end

directories.figures_cluster = fullfile(directories.figures, 'cluster');
if ~exist(directories.figures_cluster, 'dir')
    mkdir(directories.figures_cluster);
end

end
