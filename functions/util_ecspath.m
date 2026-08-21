function info = util_ecspath()
%UTIL_ECSPATH  The three ECSpress roots on the MATLAB path, without zz_notinuse.
%   Works out where the roots are from where THIS FILE sits, so a clone at another
%   location needs no edit. The three are siblings of each other:
%
%     00_mdfExtractor   mdfExtractLoader, without which ECSSession loads a
%                       recording's analysis_analog.mat as a bare uint32
%     01_ecspress       this repository, the only one guaranteed to be here
%     02_othersignal    analysis_analog, without which the analog and sleep
%                       products cannot be opened as objects
%
%   Only 01_ecspress is certain -- it is the folder this file lives two levels
%   inside. The other two are SEPARATE repositories that happen to sit beside it,
%   so cloning 01_ecspress on its own still works, minus what those two provide.
%   A missing one is named along with what stops working, rather than failing
%   later inside a load.
%
%   zz_notinuse is LEFT OFF, and that is the point of the function. genpath takes
%   every subfolder it can see, so a retired script moved into zz_notinuse stays
%   callable by name -- which is the one thing moving it was meant to stop.
%   map_directory was the case that mattered: it still writes paxfwhm_state.mat
%   into the session folders. see CLAUDE_LOG.md
%
%   OUT  info  struct
%          .root        1xN str    roots that were found and added
%          .missing     1xM str    roots that are not beside this one
%          .n_added     double     folders put on the path
%          .n_excluded  double     folders genpath offered that were left off
%          .n_removed   double     zz_notinuse folders something else had already
%                                  added, taken back off

exclude_name = "zz_notinuse";

this_file = mfilename('fullpath');
functions_dir = fileparts(this_file);
ecspress_root = fileparts(functions_dir);
program_root = fileparts(ecspress_root);

% <folder name>, <what is lost when it is not there>
sibling = [ ...
    "00_mdfExtractor", "ECSSession cannot open a recording -- mdfExtractLoader"; ...
    "02_othersignal",  "analog and sleep products load as raw data -- analysis_analog"];

% string, not char: fileparts returns char, and appending to a char row vector
% grows it by one CHARACTER instead of by one path
candidate_path = string(ecspress_root);
candidate_name = "01_ecspress";
for k = 1:size(sibling, 1)
    sibling_path = fullfile(program_root, sibling(k, 1));
    candidate_path(end+1) = sibling_path; %#ok<AGROW>
    candidate_name(end+1) = sibling(k, 1); %#ok<AGROW>
end

found_root = strings(1, 0);
missing_root = strings(1, 0);
keep_folder = strings(1, 0);
n_excluded = 0;
for k = 1:numel(candidate_path)
    root_path = candidate_path(k);
    if ~isfolder(root_path)
        missing_root(end+1) = candidate_name(k); %#ok<AGROW>
        continue
    end
    found_root(end+1) = candidate_name(k); %#ok<AGROW>

    folder_string = genpath(root_path);
    folder_list = string(strsplit(folder_string, pathsep));
    folder_list = folder_list(strlength(folder_list) > 0);
    is_excluded = contains(folder_list, exclude_name);
    n_excluded = n_excluded + sum(is_excluded);
    keep_folder = [keep_folder, folder_list(~is_excluded)]; %#ok<AGROW>
end

if isempty(keep_folder)
    error('util_ecspath:noRoot', 'no root folder found beside %s', functions_dir);
end
addpath(strjoin(keep_folder, pathsep));

% Leaving them out is only half of it: addpath adds and never removes, so an
% earlier genpath in the same session leaves zz_notinuse on the path and the
% retired file is still what a name resolves to. Take them off, and say so -- a
% session that was deliberately working on a legacy file needs to see this line.
current_list = string(strsplit(path, pathsep));
stale_list = current_list(contains(current_list, exclude_name));
if ~isempty(stale_list)
    rmpath(strjoin(stale_list, pathsep));
end

info = struct();
info.root = found_root;
info.missing = missing_root;
info.n_added = numel(keep_folder);
info.n_excluded = n_excluded;
info.n_removed = numel(stale_list);

fprintf('ecspath : %d folders, %d roots', info.n_added, numel(info.root));
if info.n_excluded > 0
    fprintf(', %d left off under %s', info.n_excluded, exclude_name);
end
if info.n_removed > 0
    fprintf(', %d already on the path REMOVED', info.n_removed);
end
fprintf('\n');
for k = 1:numel(missing_root)
    hit = strcmp(sibling(:, 1), missing_root(k));
    fprintf('   MISSING %s -- %s\n', missing_root(k), sibling(hit, 2));
end
end
