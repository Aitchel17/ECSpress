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
        info.infoname = info.infoname.name;
        tmp.infopath = fullfile(info.analyzefolder, info.infoname);
        tmp.info_table = readtable(tmp.infopath);
        
        for i = 1:height(tmp.info_table)
            info.(tmp.info_table.Field{i}) = tmp.info_table.Value{i}; % Add the field and value to the struct
        end
        
        obj.info = info;
        obj.info.analysis_savepath = savepath;

        end
        % function analog = loadanalog(obj)
        % 
        % end


        function stack = loadstack(obj, channel)
            arguments
            obj
            channel (1,:) char {mustBeMember(channel, ["ch1","ch2"])}
            end

            disp('analyze class constructor activated')

            stack = analyze_readtiff(obj.info.analyzefolder, ['*',channel,'.tif']);
        end
    
        function analog = loadanalog(obj)
            analogname = dir(fullfile(obj.info.analyzefolder,'*_analog.txt')).name;
            filename = fullfile(obj.info.analyzefolder,analogname);
            % Open the file
            fid = fopen(filename, 'r');
            if fid == -1
                error('Could not open the file.');
            end
        
            % Initialize output structs
            analog = struct();
            analog.info = struct();
            analog.data = struct();
            
            % Read header information
            section = 'header'; % Track which section we are in
            while ~feof(fid)
                line = strtrim(fgetl(fid)); % Read line and trim whitespace
                
                % Check for section headers
                if contains(line, '--- Analog Info')
                    section = 'header';
                    continue;
                elseif contains(line, '--- Analog Data')
                    section = 'data';
                    continue;
                end
        
                % Process header info
                if strcmp(section, 'header') && contains(line, ':')
                    tokens = split(line, ':'); % Split by colon
                    key = strtrim(tokens{1}); 
                    value = strtrim(tokens{2});
                    % Store in info struct
                    analog.info.(key) = value;
                end
        
                % Process data section
                if strcmp(section, 'data') && contains(line, ':')
                    tokens = split(line, ':'); % Split by colon
                    key = strtrim(tokens{1}); 
                    value = strtrim(tokens{2});
                    % Convert to numeric array
                    value = str2num(value); %#ok<ST2NM> 
                    % Store in analog data struct
                    analog.data.(key) = value;
                end
            end
        
            % Close the file
            fclose(fid);
        end
    end
end



