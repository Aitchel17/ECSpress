function dirs = dirs_central()
%DIRS_CENTRAL  Where the centralized products live, on whatever machine this is.

%   OUT  dirs  struct
%          .project_root    1x1 str  the project folder both of the next two sit under
%          .secondary_root  1x1 str  where the pipeline writes its products
%          .central         1x1 str  the centralized_*.mat files, under secondary_root
%          .figure_root     1x1 str  where finished figures go, beside secondary_root
%          .work_root       1xK str  <disk>\tmp of every disk declaring work, in letter order
%          .archive_root    1xK str  <disk>\tmp of every disk declaring archive, dormant recordings
%          .lto_root        1xK str  <disk>\lto of those, raw .mdf and their 7z
%          .mirror_root     1xM str  the roots of the disks declaring mirror, never a pipeline root.
%                                    A role's disks are a pool: an item sits on one of them,
%                                    never two. Each is "" when no disk declares the role
%          .working_root    1xK str  work_root then archive_root, those that exist today.
%                                    Empty on a machine with none
%          .source          1x1 str  which of the three rules resolved the products root

relative_project = "2023ecspress";
products_folder = "02_secondary_analysis";
figure_folder = "08_figures";
relative_products = fullfile(relative_project, products_folder);

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
dirs.source = source;

% the data volumes by their beacon file, on whatever letter Windows gave them today
volume = dirs_volumeletter();
dirs.work_root = volume_folder(volume, "work", "tmp");
dirs.archive_root = volume_folder(volume, "archive", "tmp");
dirs.lto_root = volume_folder(volume, "archive", "lto");
dirs.mirror_root = volume_folder(volume, "mirror", "");
tmp_root = [dirs.work_root, dirs.archive_root];
dirs.working_root = tmp_root(isfolder(tmp_root));
end

function folder = volume_folder(volume, role, subfolder)
% <letter today>\<subfolder> of every disk in that role, in letter order; "" when none is mounted
letter = sort(volume.letter(volume.role == role));
if isempty(letter)
    folder = "";
    return
end
folder = fullfile(letter' + "\", subfolder);
end
