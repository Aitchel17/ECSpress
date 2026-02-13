classdef tableManager < handle
    % TABLEMANAGER Centralized class for managing FWHM analysis data tables.
    % Handles loading from Excel, aggregating analysis files, parsing metadata,
    % and organizing data by vessel/depth groups.
    
    properties
        secondary_analysisPath        % Base path for analysis results
        experimentFolder    % Path to the experiment folder
        refTable            % The raw reference table loaded from Excel
        primaryTable         % The aggregated data table (with loaded .mat data)
        stateTable
        subTables
        depth_thr = 70 % 70 um as boundary between L1 and L2/3
        key_names = {'MouseID', 'Date', 'VesselID', 'Depth'};
    end
    
    methods
        function obj = tableManager(experiment_folder, analysis_base_path)
            % Constructor: Initializes paths and loads the Excel reference table         
            obj.experimentFolder = experiment_folder;
            obj.secondary_analysisPath = analysis_base_path;
            obj.refTable = obj.loadExcel();
        end
        
        function refTable = loadExcel(obj)
            % Loads the reference table from the standard directory structure
            [~, exp_name] = fileparts(obj.experimentFolder);
            save_exppath = fullfile(obj.secondary_analysisPath, exp_name);
            dirtable_path = fullfile(save_exppath, strcat(exp_name, '_dirtable.xlsx'));
            
            if ~isfile(dirtable_path)
                error('Directory table not found: %s', dirtable_path);
            end
            
            opts = detectImportOptions(dirtable_path, 'Sheet', 'reference');
            refTable = readtable(dirtable_path, opts);
            fprintf('Loaded reference table with %d rows.\n', height(obj.refTable));
        end

        function withkey_table = addkey(obj,data_table)
            % Merge with key metadata columns
            key_cols = obj.stateTable(:, obj.key_names);
            withkey_table = [key_cols, data_table];
        end
        
        function filter_refTable(obj, column, pattern)
            % Filters the RefTable for valid FWHM files and specific criteria
            % column: Column name to filter (string)
            % pattern: String pattern to search for (using contains)            
            target_data = string(obj.refTable.(column)); 
            rows = contains(target_data, pattern, 'IgnoreCase', true);
            obj.refTable = obj.refTable(rows, :);
            fprintf('Filtered reference table with %d rows.\n', height(obj.refTable));
        end
        
        function data_table = aggregateData(obj, column_name)
            % Aggregates data from individual session files based on a target column
            % column_name: The column in refTable containing filenames (e.g., 'Primary_LineFWHM')
            % returns: A table with key metadata and the loaded data
            
            refTable = obj.refTable;
            if isempty(refTable)
                warning('RefTable is empty. Cannot aggregate.');
                data_table = table();
                return;
            end
            
            % 1. Parse Column Name to Determine Folder Type
            % Expecting format like "Primary_LineFWHM" or "State_PaxFWHM"
            name_split = strsplit(column_name, '_');
            parent_type = name_split{1}; 
            
            switch lower(parent_type)
                case {'primary', 'primaryanalysis'}
                    parent_folder_col = 'PrimaryAnalysis';
                case {'peripheral'}
                    parent_folder_col = 'Peripheral';
                case {'state', 'stateanalysis'}
                    parent_folder_col = 'StateAnalysis';
                otherwise
                    warning('Unknown parent type parsed from "%s". Defaulting to empty.', column_name);
                    parent_folder_col = '';
            end
            
            % 2. Initialize Storage
            n_sessions = height(refTable);
            data_struct = struct([]);
            
            fprintf('Aggregating %d sessions from column "%s"...\n', n_sessions, column_name);
            
            % 3. Loop through sessions
            for i = 1:n_sessions
                % -- Key Metadata --
                data_struct(i).MouseID = string(refTable.MouseID{i});
                data_struct(i).Date = string(refTable.Date{i});
                data_struct(i).VesselID = string(refTable.VesselID{i});
                data_struct(i).Depth = refTable.NumericDepth(i);
                
                % -- File Loading --
                target_file = refTable.(column_name){i};
                % Check if parent folder column exists and get folder name
                parent_folder = refTable.(parent_folder_col){i};
                session_dir = refTable.Directory{i};
                % Construct full path
                full_path = fullfile(session_dir, parent_folder, target_file);
                % Load Data
                loaded_content = [];
                % Load Data
                loaded_content = [];
                raw_load = load(full_path);
                fields = fieldnames(raw_load);
                loaded_content = raw_load.(fields{1});

                % Unfold loaded content into data_struct
                fnames = fieldnames(loaded_content);
                for f = 1:numel(fnames)
                     data_struct(i).(fnames{f}) = loaded_content.(fnames{f});
                end
                
                data_struct(i).SourceFile = string(full_path);
            end
            
            % 4. Convert to Table
            data_table = struct2table(data_struct, 'AsArray', true);
            fprintf('Aggregation complete. Returned table with %d rows.\n', height(data_table));
        end
        
        function parseDepths(obj)
            % Parses 'Depth' column to numeric 'NumericDepth'
            refTable = obj.refTable;
            raw_depth = refTable.Depth;
            numeric_depth = zeros(height(refTable), 1);
            
            for i = 1:height(refTable)
                d_val = raw_depth{i};
                tmp_str = strsplit(d_val, "um");
                numeric_depth(i) = str2double(tmp_str(1));                
            end
            
            refTable.NumericDepth = numeric_depth;
            
            % Assign DepthState label
            depth_logic = numeric_depth > obj.depth_thr;
            depth_name = ["L1", "L2"];
            % +1 because logic 0->1, 1->2
            refTable.DepthLayer = depth_name(depth_logic + 1)';
            obj.refTable = refTable;
            fprintf('Depths parsed. Threshold used: %d um.\n', obj.depth_thr);
        end

        
        function vessels = aggregateVessels(obj)
            % Groups all data by Vessel and returns an array of Vessel objects   
            % Unique Vessel Identifier: MouseID + VesselID
            % Combine strings to create unique keys
            vessel_keys = string(obj.refTable.MouseID) + "_" + string(obj.refTable.VesselID);
            [unique_keys, unique_idx, key_idx] = unique(vessel_keys);
            
            % Initialize Vessel Array
            vessels = Vessel.empty(0, numel(unique_keys));
            
            fprintf('aggregating %d unique vessels from %d sessions...\n', numel(unique_keys), height(obj.refTable));

            for key_order = 1:numel(unique_keys)
                % Metadata from the first occurrence
                row_idx = unique_idx(key_order);
                
                vid = obj.refTable.VesselID{row_idx};
                mid = obj.refTable.MouseID{row_idx};
                depth_layer = obj.refTable.DepthLayer(row_idx);
                % Create Vessel Object
                newVessel = Vessel(vid, mid,depth_layer);
                
                % Find all sessions for this vessel
                session_indices = find(key_idx == key_order);
                
                for session_order = 1:numel(session_indices)
                    s_idx = session_indices(session_order);
                    
                    % Get Session Info
                    sess_id = obj.refTable.SessionID{s_idx};
                    date_val = obj.refTable.Date{s_idx};
                    depth = obj.refTable.NumericDepth(s_idx); % Ensure parseDepths ran!

                    % Extract Data from Primary/State Tables (if they exist)
                    % Requires that aggregateData() was run first to populate primaryTable/stateTable!
                    primary_struct = table2struct(obj.primaryTable(s_idx, :));    
                    state_struct = table2struct(obj.stateTable(s_idx, :));
                    
                    % Add to Vessel
                    newVessel.addSession(sess_id, date_val, depth,primary_struct, state_struct);
                end
                
                vessels(key_order) = newVessel;
            end
        end
    end
end
