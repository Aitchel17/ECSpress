function found_files = scan_analysisfolder(analysis_dir, expected_files_map)
% SCAN_ANALYSISFOLDER Scan a specific analysis folder for expected files
%
% Inputs:
%   analysis_dir       - Full path to the analysis directory (or 'NA')
%   expected_files_map - Nx2 cell array: {field_name, file_name_to_look_for}
%
% Output:
%   found_files        - Struct with field_names as keys and found filenames (or 'NA')

found_files = struct();

% Initialize all expected fields to 'NA'
for i = 1:size(expected_files_map, 1)
    field_name = expected_files_map{i, 1};
    found_files.(field_name) = 'NA';
end

% If directory is invalid, return immediately with all 'NA'
if strcmp(analysis_dir, 'NA') || ~isfolder(analysis_dir)
    return;
end

% Get directory contents
contents = dir(analysis_dir);
file_names = {contents.name};

% Check for each expected file
for i = 1:size(expected_files_map, 1)
    field_name = expected_files_map{i, 1};
    target_file = expected_files_map{i, 2};

    if any(strcmp(file_names, target_file))
        found_files.(field_name) = target_file;
    end
end
end
