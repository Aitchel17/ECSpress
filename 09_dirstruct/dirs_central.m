function dirs = dirs_central()
%DIRS_CENTRAL  Where the centralized products live, on whatever machine this is.
%   dirs_ecspath finds the pipeline's CODE by deriving it from its own location. This
%   is the other half, and it cannot use the same trick: the products do not sit
%   beside the code, and their absolute path differs between machines.
%
%   It still needs no per-machine configuration. The path splits into a part Windows
%   already knows -- OneDrive publishes its own root as an environment variable on
%   every machine it is signed in to -- and a part that is a fact about this project
%   and therefore belongs in the repository. Neither half is a setting anyone has to
%   write down.
%
%   OUT  dirs  struct
%          .project_root    1x1 str  the project folder both of the next two sit under
%          .secondary_root  1x1 str  where the pipeline writes its products
%          .central         1x1 str  the centralized_*.mat files, under secondary_root
%          .figure_root     1x1 str  where finished figures go, beside secondary_root
%          .working_root    1x1 str  the raw acquisition tree. "" on a machine that
%                                    does not carry it
%          .source          1x1 str  which of the three rules resolved the root
%
%   WORKING_ROOT IS ALLOWED TO BE MISSING, and that is the point. Every stage from
%   centralize_primary onward reads only secondary_root, so the whole figure half runs
%   on a machine with no data tree at all -- and on that machine the read-only rule
%   cannot be broken, because there is nothing to write to. The two stages that do
%   need the tree check the field themselves.

% The project's own layout below whichever OneDrive this machine has, and the live
% acquisition tree. The first is a fact about the repository; the second is a fact
% about the rig, and its absence is not an error.
relative_project = "2023ecspress";
products_folder = "02_secondary_analysis";
figure_folder = "08_figures";
relative_products = fullfile(relative_project, products_folder);
working_candidate = "G:\tmp";

% 1. an explicit override wins. integrate_analysisresult hands the root to the
%    scripts it runs this way, and that handoff must keep working
secondary_root = string(getenv('ECSPRESS_ROOT'));
source = "ECSPRESS_ROOT";

% 2. otherwise ask OneDrive where it put itself. Commercial first: a machine signed
%    in to both reports a personal root under the plain name as well
if strlength(secondary_root) == 0
    for variable_name = ["OneDriveCommercial", "OneDrive"]
        onedrive_root = string(getenv(variable_name));
        if strlength(onedrive_root) == 0
            continue
        end
        candidate = fullfile(onedrive_root, relative_products);
        if isfolder(candidate)
            secondary_root = candidate;
            source = variable_name;
            break
        end
    end
end

% 3. a machine whose layout is neither. dirs_local is gitignored, so it is the one
%    place a path may be written down by hand
if strlength(secondary_root) == 0 && exist('dirs_local', 'file') == 2
    secondary_root = string(dirs_local());
    source = "dirs_local";
end

if strlength(secondary_root) == 0
    onedrive_report = getenv('OneDriveCommercial');
    if isempty(onedrive_report)
        onedrive_report = '(no OneDrive variable set)';
    end
    error('dirs_central:noRoot', ...
        ['cannot locate the products root.\n' ...
        '   OneDrive reports : %s\n' ...
        '   looked for       : %s under it\n' ...
        'Set ECSPRESS_ROOT, or put a dirs_local.m on the path returning that ' ...
        'folder as a string.'], onedrive_report, relative_products);
end

% the project folder is the parent of the products folder however the root was
% resolved, so figure_root follows an ECSPRESS_ROOT override as well
dirs.project_root = string(fileparts(secondary_root));
dirs.secondary_root = secondary_root;
dirs.central = fullfile(secondary_root, "centralized");
dirs.figure_root = fullfile(dirs.project_root, figure_folder);
dirs.working_root = "";
if isfolder(working_candidate)
    dirs.working_root = working_candidate;
end
dirs.source = source;
end
