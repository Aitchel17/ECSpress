function [info] = analyze_init()
% Read _info.txt file
    info = struct();
    info.analyzefolder = uigetdir;
    info.infoname = dir(fullfile(info.analyzefolder, '*_info.txt'));
    if length(info.infoname) ~= 1
        disp('_info.txt file not exist or plural num info.txt: ' + string(length(info.infoname)));
        info = rmfield(info, 'info_dir');
    else
        info.infoname = info.infoname.name;
        tmp.infopath = fullfile(info.analyzefolder, info.infoname);
        tmp.info_table = readtable(tmp.infopath);
        
        for i = 1:height(tmp.info_table)
            info.(tmp.info_table.Field{i}) = tmp.info_table.Value{i}; % Add the field and value to the struct
        end
    end
end

