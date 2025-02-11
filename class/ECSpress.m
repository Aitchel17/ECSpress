classdef ECSpress
    % ANALYZE Summary of this class goes here
    %   Detailed explanation goes here
    properties
        info
        analog
        stackch1
        stackch2
    end
    
    methods
        function obj = ECSpress(savepath)

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
        obj.info.analysis_savepath = savepath;

        end


        function stack = loadstack(obj, channel)
            arguments
            obj
            channel
            end

            disp('analyze class constructor activated')
            obj.info = analyze_init();
            obj.analog = analyze_readanalog(obj.info.analyzefolder);
            obj.stackch1 = analyze_readtiff(obj.info.analyzefolder, '*ch1.tif');
            obj.stackch2 = analyze_readtiff(obj.info.analyzefolder, '*ch2.tif');


        end
    end
end



