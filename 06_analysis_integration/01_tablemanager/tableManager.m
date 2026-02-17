classdef tableManager < handle
    % TABLEMANAGER Centralized class for managing FWHM analysis data tables.
    % Handles loading from Excel, aggregating analysis files, parsing metadata,
    % and organizing data by vessel/depth groups.
    
    properties
        masterdirtable_path        % Base path for analysis results
        refTable            % The raw reference table loaded from Excel
        aggregatedTable
        subTables
        depth_thr = 70 % 70 um as boundary between L1 and L2/3
        action_log
    end
    
    methods
        function obj = tableManager(masterdirtable_path)
            % Constructor: Initializes paths and loads the Excel reference table         
            obj.masterdirtable_path = masterdirtable_path;
            obj.refTable = obj.loadExcel();
            obj.action_log.key_names = ["MouseID","Date","VesselID","NumericDepth","NumericResolution"];
        end
        
        function refTable = loadExcel(obj)
            % Loads the reference table from the standard directory structure
            
            if ~isfile(obj.masterdirtable_path)
                error('Directory table not found: %s', dirtable_path);
            end
            
            opts = detectImportOptions(obj.masterdirtable_path, 'Sheet', 'reference');
            refTable = readtable(obj.masterdirtable_path, opts);
            fprintf('Loaded reference table with %d rows.\n', height(obj.refTable));
        end

        function withkey_table = addkey(obj,data_table)
            % Merge with key metadata columns
            key_cols = obj.stateTable(:, obj.action_log.key_names);
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
            obj.action_log.(column) = pattern; 
        end


        function addnest2subtable(obj, nestname_arr)
            key_columns = obj.aggregatedTable(:, obj.action_log.key_names);
            nesttable = obj.aggregatedTable;
            for name_idx = 1:numel(nestname_arr)
                obj.subTables.(nestname_arr(name_idx)) = explode_nest(nesttable.(nestname_arr(name_idx)),key_columns);
            end
            obj.action_log.nestname_arr = nestname_arr;
        end
        
        function aggregateData(obj, column_name)
            % Aggregates data from individual session files based on a target column
            % column_name: The column in refTable containing filenames (e.g., 'Primary_LineFWHM')
            % returns: A table with key metadata and the loaded data
            
            refTable = obj.refTable;
            if isempty(obj.refTable)
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
                key_names = obj.action_log.key_names;
                for key_idx = 1:numel(key_names)
                    target_col = refTable.(key_names{key_idx});
                    data_struct(i).(key_names{key_idx}) = target_col(i);
                end
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
            obj.aggregatedTable = data_table;
            obj.action_log.aggregatedTable = column_name;
            fprintf('Aggregation complete. Returned table with %d rows.\n', height(data_table));
        end
        
        function parseParams(obj)
            % Parses 'Depth' column to numeric 'NumericDepth'
            refTable = obj.refTable;
            n_session = height(refTable);
            % depth init
            raw_depth = refTable.Depth;
            numeric_depth = zeros(n_session, 1);
            % resolution init
            raw_resolution = refTable.Resolution;
            numeric_resolution = zeros(n_session,1);
            
            for sidx = 1:n_session
                d_val = raw_depth{sidx};
                tmp_str = strsplit(d_val, "um");
                numeric_depth(sidx) = str2double(tmp_str(1));   
                r_val = raw_resolution{sidx};
                tmp_str = strsplit(r_val, "µm");
                numeric_resolution(sidx) = str2double(tmp_str(1));
            end
            
            refTable.NumericDepth = numeric_depth;
            refTable.NumericResolution = numeric_resolution;
            % Assign DepthState label
            depth_logic = numeric_depth > obj.depth_thr;
            depth_name = ["L1", "L2"];
            % +1 because logic 0->1, 1->2
            refTable.DepthLayer = depth_name(depth_logic + 1)';
            obj.refTable = refTable;
            fprintf('Depths parsed. Threshold used: %d um.\n', obj.depth_thr);
            fprintf('Resolution parsed. \n');

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

        function save2disk(obj,savename)
            tableManager = obj;
            tableManager.subTables = [];
            savedir = fileparts(obj.masterdirtable_path);
            savepath = fullfile(savedir, savename);
            save(savepath,"tableManager")
        end
        
    end
    methods (Static)
        function loaded_table = load_recon(masterDirTable_path,mtable_name)
            load_folderdir = fileparts(masterDirTable_path);
            load_filedir = fullfile(load_folderdir, mtable_name);
            load_str = load(load_filedir);
            nestname_arr = load_str.tableManager.action_log.nestname_arr;
            load_str.tableManager.addnest2subtable(nestname_arr)
            loaded_table = load_str.tableManager; 
        end
    end
end
