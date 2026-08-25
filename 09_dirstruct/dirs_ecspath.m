function dirs_ecspath()
%DIRS_ECSPATH  Put the pipeline and its two siblings on the MATLAB path.
%   The one file the path cannot find for you, so every script opens by adding
%   this folder by hand and calling it. It derives the roots from its own
%   location, so moving it one level changes what it finds.

sibling = ["00_mdfExtractor",  "02_othersignal"];
exclude_name = "zz_notinuse"; % Exclude legacy

this_file = mfilename('fullpath'); % where I am?
dirstruct_dir = fileparts(this_file); % where is mother
ecspress_root = fileparts(dirstruct_dir); % where is grandmother
program_root = fileparts(ecspress_root); % who is ancester

candidate_path = string(ecspress_root);
candidate_name = "01_ecspress";
%%
for k = 1:numel(sibling)
    sibling_path = fullfile(program_root, sibling(k));
    candidate_path(end+1) = sibling_path; 
    candidate_name(end+1) = sibling(k); 
end
%% Scan folder
keep_folder = strings(1, 0);
for k = 1:numel(candidate_path)
    root_path = candidate_path(k);
    folder_string = genpath(root_path);
    folder_list = string(strsplit(folder_string, pathsep));
    folder_list = folder_list(strlength(folder_list) > 0);
    is_excluded = contains(folder_list, exclude_name);
    keep_folder = [keep_folder, folder_list(~is_excluded)]; 
end
%%

addpath(strjoin(keep_folder, pathsep));



end
