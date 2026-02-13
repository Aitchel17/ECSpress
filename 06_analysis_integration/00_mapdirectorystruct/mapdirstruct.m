function dirstruct_table = mapdirstruct(experiment_folder, primary_map, peri_map, state_map)
%MAPDIRECTORY Summary of this function goes here
%   Detailed explanation goes here
    exp_dir = dir(experiment_folder);
    exp_dir = exp_dir([exp_dir.isdir] & ~ismember({exp_dir.name}, {'.', '..'}));
    % Go down to mouse_dir which contains imaging date folder, ex. 251012_hql090_sleep_....
    % Determine current mouse id
    % Double check that mouse folder contains that mouse
    % If labeled detect third delimeter as type of imaging session, ex. sleep
    % Initialize columns structure
    cols = struct();
    cols.MouseID = {};
    cols.Date = {};
    cols.SessionType = {};
    cols.SessionID = {};
    cols.Directory = {};
    
    % Metadata
    cols.VesselID = {};
    cols.Power = {};
    cols.Depth = {};
    cols.Notes = {};
    cols.xlocation = {};
    cols.ylocation = {};
    cols.zlocation = {};
    cols.PercentPower = {};
    cols.Wavelength = {};
    
    % Main Files
    cols.Analog = {};
    cols.Ch1 = {};
    cols.Ch2 = {};
    cols.Eye = {};
    cols.Info = {};
    cols.Motion = {};
    cols.Whisker = {};
    cols.Peripheral = {};
    cols.PrimaryAnalysis = {};
    cols.StateAnalysis = {};
    
    % Primary Analysis Files
    cols.Primary_RadonResult = {};
    cols.Primary_RoiList = {};
    cols.Primary_PolarCluster = {};
    cols.Primary_paxFWHM = {};
    
    % Peripheral Files
    cols.Peripheral_AnalogAnalysis = {};
    cols.Peripheral_BehaviorAnalysis = {};
    cols.Peripheral_SleepScore = {};
    
    % State Analysis Files
    cols.State_PaxFWHM = {};
    
    for mouse_idx = 1:numel(exp_dir)
        current.mouseid = exp_dir(mouse_idx).name;
        mouse.directory = fullfile(exp_dir(mouse_idx).folder, current.mouseid);
        mouse.dirs = dir(mouse.directory);
        mouse.dirs = mouse.dirs([mouse.dirs.isdir] & ~ismember({mouse.dirs.name}, {'.', '..', 'zstack'}));
    
        for datefolder_idx = 1:numel(mouse.dirs)
            current.dateid = mouse.dirs(datefolder_idx).name;
            datefolder.directory = fullfile(mouse.dirs(datefolder_idx).folder, current.dateid);
            datefolder.dirs = dir(datefolder.directory);
            datefolder.dirs = datefolder.dirs([datefolder.dirs.isdir] & ~ismember({datefolder.dirs.name}, {'.', '..'}));
    
            current.datefoldername = strsplit(current.dateid, '_');
            current.date = current.datefoldername{1};
            current.mouseid_datefolder = current.datefoldername{2};
    
            if numel(current.datefoldername) >= 3
                current.sessiontype = current.datefoldername{3};
            else
                current.sessiontype = 'unknown';
            end
    
            for session_idx = 1:numel(datefolder.dirs)
                session.name = datefolder.dirs(session_idx).name;
                s_split = strsplit(session.name, '_');
                current_session_id = s_split{3};
    
                % --- Scan Session Folder ---
                session.directory = fullfile(datefolder.dirs(session_idx).folder, session.name);
                session_data = scan_sessionfolder(session.directory);
    
                % --- Scan Primary Analysis Folder ---
                if ~strcmp(session_data.PrimaryAnalysis, 'NA')
                    p_dir = fullfile(session.directory, session_data.PrimaryAnalysis);
                    found_primary = scan_analysisfolder(p_dir, primary_map);
                else
                    found_primary = scan_analysisfolder('NA', primary_map);
                end
    
                % --- Scan Peripheral Folder ---
                if ~strcmp(session_data.Peripheral, 'NA')
                    p_dir = fullfile(session.directory, session_data.Peripheral);
                    found_peri = scan_analysisfolder(p_dir, peri_map);
                else
                    found_peri = scan_analysisfolder('NA', peri_map);
                end
    
                % --- Scan State Analysis Folder ---
                if isfield(session_data, 'State_analysis') && ~strcmp(session_data.State_analysis, 'NA')
                    p_dir = fullfile(session.directory, session_data.State_analysis);
                    found_state = scan_analysisfolder(p_dir, state_map);
                else
                    found_state = scan_analysisfolder('NA', state_map);
                end
    
                % --- Append Data to Columns ---
                % Basic Info
                cols.MouseID{end+1,1} = current.mouseid;
                cols.Date{end+1,1} = current.date;
                cols.SessionType{end+1,1} = current.sessiontype;
                cols.SessionID{end+1,1} = current_session_id;
                cols.Directory{end+1,1} = session.directory;
    
                % Metadata (maps to specific column names)
                cols.VesselID{end+1,1} = session_data.vesselID;
                cols.Power{end+1,1} = session_data.power;
                cols.Depth{end+1,1} = session_data.depth;
                cols.Notes{end+1,1} = session_data.comments; % 'comments' maps to 'Notes'
                cols.xlocation{end+1,1} = session_data.objx;
                cols.ylocation{end+1,1} = session_data.objy;
                cols.zlocation{end+1,1} = session_data.objz;
                cols.PercentPower{end+1,1} = session_data.power_percent;
                cols.Wavelength{end+1,1} = session_data.wavelength;
    
                % Main Files
                cols.Analog{end+1,1} = session_data.Analog;
                cols.Ch1{end+1,1} = session_data.Ch1;
                cols.Ch2{end+1,1} = session_data.Ch2;
                cols.Eye{end+1,1} = session_data.Eye;
                cols.Info{end+1,1} = session_data.Info;
                cols.Motion{end+1,1} = session_data.Motion;
                cols.Whisker{end+1,1} = session_data.Whisker;
                cols.Peripheral{end+1,1} = session_data.Peripheral;
                cols.PrimaryAnalysis{end+1,1} = session_data.PrimaryAnalysis;
                if isfield(session_data, 'State_analysis')
                    cols.StateAnalysis{end+1,1} = session_data.State_analysis;
                else
                    cols.StateAnalysis{end+1,1} = 'NA';
                end
    
                % Primary Analysis Files
                cols.Primary_RadonResult{end+1,1} = found_primary.RadonResult;
                cols.Primary_RoiList{end+1,1} = found_primary.RoiList;
                cols.Primary_PolarCluster{end+1,1} = found_primary.PolarCluster;
                cols.Primary_paxFWHM{end+1,1} = found_primary.PaxFWHM; % Maps 'PaxFWHM' to 'LineFWHM' column
    
                % Peripheral Files
                cols.Peripheral_AnalogAnalysis{end+1,1} = found_peri.AnalogAnalysis;
                cols.Peripheral_BehaviorAnalysis{end+1,1} = found_peri.BehaviorAnalysis;
                cols.Peripheral_SleepScore{end+1,1} = found_peri.SleepScore;
    
                % State Analysis Files
                cols.State_PaxFWHM{end+1,1} = found_state.paxfwhm_state;
            end
        end
    end
    
    % Create Table
    dirstruct_table = struct2table(cols);
    
    % Sort by date
    if ~isempty(dirstruct_table)
        dirstruct_table = sortrows(dirstruct_table, 'Date');
    end
end

