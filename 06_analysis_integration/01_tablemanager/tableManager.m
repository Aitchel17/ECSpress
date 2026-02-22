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
        analysis_table
        filtered_table
        filtLogics
        numeric_tables
        data_tables
    end
    
    methods % Methods for generating First Normal form of analysis
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
    
    methods % Methods for query and additional data column generation
        function apply_filter(obj)
           combinedlogic = all(table2array(struct2table(obj.filtLogics)), 2);
           fprintf('Filtered %d from % d rows.\n', sum(combinedlogic), height(obj.analysis_table));
           obj.filtered_table = obj.analysis_table(combinedlogic,:);
           obj.action_log.resolution_applied = false;
        end

        function meanFrom2(obj,datacolName,NewcolName,start_fraction,end_fraction)
            dlength = height(obj.filtered_table);
            numeric_arrs = obj.filtered_table.(datacolName);
            result_col = zeros(dlength,1);
            for arr_idx = 1:dlength
                numeric_arr = numeric_arrs{arr_idx};
                arrlength = length(numeric_arr);
                start_idx = floor(arrlength * start_fraction);
                end_idx = floor(arrlength * end_fraction);
                result_col(arr_idx) = mean(numeric_arr(start_idx:end_idx));
            end

            obj.filtered_table.(NewcolName) = result_col;

            current_cols = string(obj.action_log.numeric_colnames);
            chk_duplication = ismember(current_cols, string(NewcolName));
            if ~any(chk_duplication)
                obj.action_log.numeric_colnames = [current_cols, string(NewcolName)];
            end
            fprintf("mean from %s to %s calculated and added as %s column\n", start_fraction,end_fraction, NewcolName)
        end
 
        function addPrctilecol(obj,datacolName,NewcolName,percentile)
            dlength = height(obj.filtered_table);
            numeric_arrs = obj.filtered_table.(datacolName);
            result_col = zeros(dlength,1);
            for arr_idx = 1:dlength
                numeric_arr = numeric_arrs{arr_idx};
                result_col(arr_idx) = prctile(numeric_arr,percentile);  
            end
            obj.filtered_table.(NewcolName) = result_col;
            
            current_cols = string(obj.action_log.numeric_colnames);
            chk_duplication = ismember(current_cols, string(NewcolName));
            if ~any(chk_duplication)
                obj.action_log.numeric_colnames = [current_cols, string(NewcolName)];
            end
            fprintf("%d Percentile calculated and added as %s column\n", percentile, NewcolName)
        end


        function apply_resolution(obj,scale_colname,data_colnames,numeric_colnames)
            % Applies resolution scaling (multiplication) to specified columns
            if obj.action_log.resolution_applied
                fprintf("Resolution already applied")
            else
                resolution_vec = obj.filtered_table.(scale_colname);
                numeric_data = obj.filtered_table{:, numeric_colnames};
                obj.filtered_table{:, numeric_colnames} = numeric_data .* resolution_vec;
                % 2. Scale Data Columns (Cell arrays of vectors)
                % Explicit loop over rows as requested
                
                for i = 1:numel(data_colnames)
                    col_name = data_colnames{i};
                    % Loop over rows
                    current_col = obj.filtered_table.(col_name);
                    for j = 1:height(obj.filtered_table)
                        % Multiply content
                        if ~isempty(current_col{j})
                            current_col{j} = current_col{j} * resolution_vec(j);
                        end
                    end
                    obj.filtered_table.(col_name) = current_col;                
                end
                obj.action_log.resolution_applied = true;
                obj.action_log.numeric_colnames = numeric_colnames;
                obj.action_log.data_colnames = data_colnames;

                obj.action_log.key_names = ["DataType","state_name","MouseID","VesselID","Date","bout_idx"];
                fprintf("Units are scaled to resolution column\n")
            end
        end

        function get_numericsummary(obj, new_keyname, target_tablename)
            var_names = string(obj.action_log.numeric_colnames);
            
            [target_table, summary_key, aveTable_name, ci95Table_name] = obj.prepare_summary("numeric_tables", target_tablename, new_keyname);
           
            stats = groupsummary(target_table, summary_key, {'mean', 'std'}, var_names);  % Calculate Mean, Std

                
            temp_table = stats(:, [summary_key, "mean_" + var_names]);  %             % Select columns starting with mean_ and keys and store as Mean Table
            
            % Rename 'mean_' columns back to original variable names
            for i = 1:numel(var_names)
                current_name = "mean_" + var_names(i);
                if ismember(current_name, temp_table.Properties.VariableNames)
                    temp_table.Properties.VariableNames(current_name) = var_names(i);
                end
            end
            obj.numeric_tables.(aveTable_name) = temp_table;
            
            % Calculate 95% CI Table (Margin of Error)
            ci_table = stats(:, summary_key);
            
            for i = 1:length(var_names)
                vname = var_names(i);
                m_col = "mean_" + vname;
                s_col = "std_" + vname;
                
                % Standard Error and Margin of Error
                n = stats.GroupCount;
                ci95 = obj.calculate_ci95(stats.(s_col), n);
                
                % Store in table using original variable name (e.g., raw_mean)
                ci_table.(vname) = ci95;
            end
            
            obj.save_summary("numeric_tables", aveTable_name, temp_table, ci95Table_name, ci_table, new_keyname);
        end

        function get_datasummary(obj, new_keyname, target_tablename)
            % Summarizes cell-array columns containing 1D data arrays,
            % computing element-wise mean (and CI) across groups.
            % Analogous to get_numericsummary but for data_colnames.
            var_name = obj.action_log.data_colnames{1};

            [target_table, summary_key, aveTable_name, ci95Table_name] = obj.prepare_summary("data_tables", target_tablename, new_keyname);

            % Find unique group combinations
            group_table = target_table(:, summary_key);
            [unique_groups, ~, group_idx] = unique(group_table, 'rows', 'stable');
            n_groups = height(unique_groups);

            % Preallocate output tables
            ave_table = unique_groups;
            ci_table  = unique_groups;
            ave_table.(var_name) = cell(n_groups, 1);
            ci_table.(var_name)  = cell(n_groups, 1);
            col_data = target_table.(var_name);   % cell array of 1D arrays
            % Loop over groups and compute element-wise mean / CI
            for g = 1:n_groups
                row_idx = find(group_idx == g);
                n = numel(row_idx);
                
                if n == 1                     % Only one observation – bypass mean, copy directly
                    ave_table.(var_name){g} = col_data{row_idx};
                    ci_table.(var_name){g}  = zeros(size(col_data{row_idx}));
                    continue;
                end

                % Stack: ensure result is [n x len]
                mat = cell2mat(col_data(row_idx));
                mu  = mean(mat, 1, 'omitnan');
                sig = std(mat,  0, 1, 'omitnan');
                ci95 = obj.calculate_ci95(sig, n);

                ave_table.(var_name){g} = mu;
                ci_table.(var_name){g}  = ci95;
                
            end

            obj.save_summary("data_tables", aveTable_name, ave_table, ci95Table_name, ci_table, new_keyname);
        end

    end 

    methods (Access = private)
        function ci95 = calculate_ci95(obj, std_val, n)
            sem = std_val ./ sqrt(n);
            df = max(n - 1, 1); % Ensure df is at least 1 for tinv to avoid NaN if n=1
            t_val = tinv(0.975, df);
            ci95 = t_val .* sem;
        end

        function [target_table, summary_key, aveTable_name, ci95Table_name] = prepare_summary(obj, table_typename, target_tablename, new_keyname)
            aveTable_name = strcat(new_keyname, "_ave");
            ci95Table_name = strcat(new_keyname, "_ci95");
            [target_table, summary_key] = obj.retrieve_keytable(table_typename, target_tablename, new_keyname);
        end

        function save_summary(obj, table_typename, aveTable_name, ave_table, ci95Table_name, ci_table, new_keyname)
            obj.(table_typename).(aveTable_name) = ave_table;
            obj.(table_typename).(ci95Table_name) = ci_table;
            name_parts = strsplit(table_typename, "_");
            fprintf("%s %s summarytable added\n", name_parts{1}, new_keyname);
        end

        function [target_table,summary_key] = retrieve_keytable(obj,table_typename,target_tablename,new_keyname)
                        % Key name array retrieval
            newkey_name = strcat(new_keyname, "_key");
            if strcmp(target_tablename, 'filtered_table')
                target_table = obj.filtered_table;
                key_names = obj.action_log.key_names;
            else
                target_table = obj.(table_typename).(target_tablename);
                actionLog_keyfield = strsplit(target_tablename,"_");
                actionLog_keyfield = strjoin(actionLog_keyfield(1:end-1),"_");
                actionLog_keyfield = strcat(actionLog_keyfield,"_key");
                key_names = obj.action_log.(actionLog_keyfield);
            end
             % exclude
            logic_key = ismember(key_names,new_keyname);
            if  any(logic_key)
                summary_key = key_names(~logic_key);
                obj.action_log.(newkey_name) = summary_key;
            else
                fprintf("target_keyname is not in key_names at actionlog")
            end
        end
    end

    methods (Static) % Methods for load and reconstruct tables
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
