function [data] = analyze_readanalog(folder_path)
    % Initialize an empty structure for the output data
    data = struct();
    
    % Directory file searches
    io.ball_dir = dir(fullfile(folder_path, '*_ball.txt'));
    io.emg_dir = dir(fullfile(folder_path, '*_emg.txt'));
    io.ecog_dir = dir(fullfile(folder_path, '*_ecog.xlsx'));
        
    % Read _ball.txt file
    if length(io.ball_dir) ~= 1
        disp('_ball.txt file not exist or plural num ball.txt: ' + string(length(io.ball_dir)));
        io = rmfield(io, 'ball_dir');
    else
        io.ball_dir = fullfile(folder_path, io.ball_dir.name);
        tmp.ball_opts = detectImportOptions(io.ball_dir, 'FileType', 'text','VariableNamingRule','preserve');
        data.ball_table = readtable(io.ball_dir, tmp.ball_opts);   
    end
    
    % Read _emg.txt file
    if length(io.emg_dir) ~= 1
        disp('_emg.txt file not exist or plural num emg.txt: ' + string(length(io.emg_dir)));
        io = rmfield(io, 'emg_dir');
    else
        io.emg_dir = fullfile(folder_path, io.emg_dir.name);
        tmp.emg_opts = detectImportOptions(io.emg_dir, 'FileType', 'text','VariableNamingRule','preserve');
        data.emg_table = readtable(io.emg_dir, tmp.emg_opts);
    end
    
    % Read _ecog.xlsx file
    if length(io.ecog_dir) ~= 1
        disp('_ecog.xlsx file not exist or plural num ecog.xlsx: ' + string(length(io.ecog_dir)));
        io = rmfield(io, 'ecog_dir');
    else
        io.ecog_dir = fullfile(folder_path, io.ecog_dir.name);
        tmp.table_ecog = readtable(io.ecog_dir, 'Sheet', 'ecog_mtspec','VariableNamingRule','preserve');
        
        % Extract f_ axis (frequency)
        tmp.ftable = table2array(tmp.table_ecog(:, 1));  
        data.ecog.faxis = flip(tmp.ftable');
        
        % Extract t_ axis (time)
        tmp.ttable = tmp.table_ecog.Properties.VariableNames(2:end);
        for i = 1:numel(tmp.ttable)
            tmp.ttable{i} = str2double(strrep(tmp.ttable{i}(2:end),'_','.'));
        end
        data.ecog.taxis = cell2mat(tmp.ttable);
        
        % Extract log_norm_spectrum
        data.ecog.spectrum = table2array(tmp.table_ecog(:, 2:end));
    end
end