function delete_emptydir(target_path)
    % folder_path: folder_path fullfile( folderpath , *)
    % 1. get directory info
    % 2.1 if subdirectory exist from each sub directory (except . and ..) go back to 1.
    % 2.2 if subdirectory does not exist remove directory
    
    dir_info = dir(target_path); %1.
    for i = 1:length(dir_info)
        % Check if it is a directory and not '.' or '..'
        if dir_info(i).isdir && ~strcmp(dir_info(i).name, '.') && ~strcmp(dir_info(i).name, '..')
            subfolder = fullfile(dir_info(i).folder, dir_info(i).name); % 2.1
            delete_emptydir(subfolder);
        end
    end
    try
        rmdir(target_path); % 2.2
        disp(['Removed empty directory: ', target_path]);
    catch
    end
end