function fwhm_table = aggregate_fwhm(ref_table)
%AGGREGATE_FWHM Aggregates FWHM analysis data from primary and state analysis files.
%   Iterates through the provided reference table (filtered for sessions of interest),
%   loads 'paxfwhm.mat' and 'paxfwhm_state.mat', and reorganizes the data into
%   a master table with metadata columns and cell/struct data columns.

    % Initialize storage
    % We will build a struct array, but we need to pre-allocate or grow it carefully.
    % To be safe with different types, we'll initialize the first element if we have data.
    data_struct = struct([]); 
    
    row_cnt = 0;

    % Loop through each session
    for i = 1:height(ref_table)
        
        % 1. Get Session Info
        session_dir = ref_table.Directory{i};
        mouse_id = ref_table.MouseID(i); 
        date = ref_table.Date(i);
        session_id = ref_table.SessionID(i);
        
        % Metadata (Handle potential missing values/types)
        vessel_id = string(ref_table.VesselID(i));
        depth = ref_table.Depth(i);
        if iscell(depth), depth = str2double(depth); end 
        
        % File Names from Table
        p_file_name = ref_table.LineFWHM(i);     % e.g., 'paxfwhm.mat'
        s_file_name = ref_table.PaxFWHM_state(i); % e.g., 'paxfwhm_state.mat'
        
        % Construct Full Paths (handle string/cell)
        if iscell(p_file_name), p_file_name = p_file_name{1}; end
        if iscell(s_file_name), s_file_name = s_file_name{1}; end
        if iscell(session_dir), session_dir = session_dir{1}; end
        
        p_path = fullfile(session_dir, "primary_analysis",p_file_name);
        s_path = fullfile(session_dir,"state_analysis" ,s_file_name);
        
        if ~isfile(p_path)
            continue; 
        end

        row_cnt = row_cnt + 1;
        
        % --- Metadata Columns ---
        data_struct(row_cnt).MouseID = string(mouse_id);
        data_struct(row_cnt).Date = string(date);
        data_struct(row_cnt).SessionID = string(session_id);
        data_struct(row_cnt).VesselID = vessel_id;
        data_struct(row_cnt).Depth = depth;
        data_struct(row_cnt).Directory = string(session_dir);
        
        % --- 3. Load Primary Analysis (Time Series) ---
        try
            loaded_p = load(p_path, 'line_fwhm');
            if isfield(loaded_p, 'line_fwhm')
                obj_p = loaded_p.line_fwhm; 
                
                % Extract Time Series (Always present)
                data_struct(row_cnt).Time = {obj_p.t_axis};
                
                % Primary Thickness
                if isfield(obj_p, 'thickness') && ~isempty(obj_p.thickness)
                    data_struct(row_cnt).Primary_Thickness = obj_p.thickness;
                else
                    data_struct(row_cnt).Primary_Thickness = struct();
                end
                
                % Primary Displacement
                if isfield(obj_p, 'displacement') && ~isempty(obj_p.displacement)
                    data_struct(row_cnt).Primary_Displacement = obj_p.displacement;
                else
                    data_struct(row_cnt).Primary_Displacement = struct();
                end
                
                % Primary Index (Boundaries)
                if isfield(obj_p, 'idx') && ~isempty(obj_p.idx)
                    data_struct(row_cnt).Primary_Idx = obj_p.idx;
                else
                    data_struct(row_cnt).Primary_Idx = struct();
                end

            else
                 warning('Variable line_fwhm not found in %s', p_path);
                 data_struct(row_cnt).Time = {[]};
            end
             
        catch ME
            warning('Failed to load primary analysis for %s: %s', session_dir, ME.message);
            data_struct(row_cnt).Time = {[]};
        end
        
        % --- 4. Load State Analysis (Subdivided Info) ---
        % Extract EVERYTHING from State Analysis as requested
        
        if isfile(s_path)
            try
                loaded_s = load(s_path, 'state_linefwhm');
                if isfield(loaded_s, 'state_linefwhm')
                    obj_s = loaded_s.state_linefwhm;
                    
                    % 1. State Summary
                    if ~isempty(obj_s.state_summary)
                       data_struct(row_cnt).State_Summary = obj_s.state_summary;
                    else
                       data_struct(row_cnt).State_Summary = struct();
                    end

                    % 2. Power Density
                    if ~isempty(obj_s.powerdensity)
                        data_struct(row_cnt).State_PowerDensity = obj_s.powerdensity;
                    else
                        data_struct(row_cnt).State_PowerDensity = struct();
                    end
                    
                    % 3. Band Decomposition
                    if ~isempty(obj_s.band_decomposition)
                        data_struct(row_cnt).State_BandDecomp = obj_s.band_decomposition;
                    else
                        data_struct(row_cnt).State_BandDecomp = struct();
                    end
                    
                    % 4. Transitions
                    if ~isempty(obj_s.transition)
                        data_struct(row_cnt).State_Transition = obj_s.transition;
                    else
                        data_struct(row_cnt).State_Transition = struct();
                    end
                    
                    % 5. Peak/Trough
                    if ~isempty(obj_s.peak_trough)
                        data_struct(row_cnt).State_PeakTrough = obj_s.peak_trough;
                    else
                         data_struct(row_cnt).State_PeakTrough = struct();
                    end
                    
                    % 6. State Indices
                     if ~isempty(obj_s.state_idx)
                        data_struct(row_cnt).State_Idx = obj_s.state_idx;
                    else
                         data_struct(row_cnt).State_Idx = struct();
                    end
                end
            catch ME
                warning('Failed to load state analysis for %s: %s', session_dir, ME.message);
            end
        end
    end

    % Convert to Table - Handling missing fields is crucial
    % If one row has 'Thickness_newmetric' and others don't, struct2table might error
    % Use our 'struct2table_robust' approach or fill missing fields
    
    if row_cnt > 0
        % Get all unique field names
        all_fields = fieldnames(data_struct);
        % Ensure every struct in the array has all fields
        for i = 1:numel(data_struct)
            missing = setdiff(all_fields, fieldnames(data_struct(i)));
            for m = 1:numel(missing)
                data_struct(i).(missing{m}) = {[]}; % Fill with empty cell
            end
        end
        fwhm_table = struct2table(data_struct);
    else
        fwhm_table = table();
    end
end
