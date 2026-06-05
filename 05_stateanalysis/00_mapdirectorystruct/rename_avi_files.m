function rename_avi_files(root_dir)
% RENAME_AVI_FILES Recursively renames eye.avi and whisker.avi to match *info.txt prefix
%
% Usage: rename_avi_files('G:\path\to\start\folder')

if nargin < 1
    root_dir = uigetdir(pwd, 'Select Root Directory');
    if root_dir == 0
        return;
    end
end

% Find all subdirectories
subdirs = genpath(root_dir);

% Split into cell array of paths (semi-colon separated in Windows)
path_list = strsplit(subdirs, ';');

% Counter for reporting
count_renamed = 0;

for i = 1:length(path_list)
    current_dir = path_list{i};
    if isempty(current_dir)
        continue;
    end

    % Look for *info.txt file to get the prefix
    info_files = dir(fullfile(current_dir, '*info.txt'));

    if isempty(info_files)
        continue;
    end

    % Assume the first info file gives the correct prefix
    % Standard format: PREFIX_info.txt
    % We want to extract "PREFIX_"

    % Example: HQL088_sleep250927_009_info.txt -> HQL088_sleep250927_009_
    info_name = info_files(1).name;

    % Using regex or simple string replacement to get prefix
    % Pattern: Anything before "info.txt"
    idx = strfind(info_name, 'info.txt');
    if isempty(idx)
        continue;
    end

    prefix = info_name(1:idx-1);

    % Check for eye.avi
    process_file(current_dir, 'eye.avi', prefix);

    % Check for whisker.avi
    process_file(current_dir, 'whisker.avi', prefix);

    count_renamed = count_renamed + 1;
end

fprintf('Finished checking directories.\n');
end

function process_file(folder, target_filename, prefix)
target_path = fullfile(folder, target_filename);

if exist(target_path, 'file')

    % Construct new name
    new_filename = [prefix target_filename];
    new_path = fullfile(folder, new_filename);

    % Avoid overwriting if it already exists (or is the same)
    if strcmp(target_path, new_path)
        return;
    end

    if exist(new_path, 'file')
        fprintf('Skipping: %s already exists.\n', new_filename);
    else
        movefile(target_path, new_path);
        fprintf('Renamed: %s -> %s\n', target_filename, new_filename);
    end
end
end
