function write_dirtable(dirstruct_table, dirtable_dir)
%WRITE_DIRTABLE Update directory mapping Excel with smart merging and logging.
%   Orchestrates the loading, updating, and saving of directory tables.

    if ~isfile(dirtable_dir)
        create_new_tables(dirstruct_table, dirtable_dir);
    else
        load_and_update_tables(dirstruct_table, dirtable_dir);
    end
end

function create_new_tables(dirstruct_table, dirtable_dir)
%CREATE_NEW_TABLES Initialize tables if file doesn't exist.
    
    % Initialize log table with current time
    log_table = dirstruct_table;
    var_names = log_table.Properties.VariableNames;
    for i = 1:numel(var_names)
        log_table.(var_names{i}) = repmat(datetime("now"), height(log_table), 1);
    end

    writetable(dirstruct_table, dirtable_dir, "Sheet", "reference");
    writetable(dirstruct_table, dirtable_dir, "Sheet", "auto_mapped");
    writetable(log_table, dirtable_dir, "Sheet", "log");
end

function load_and_update_tables(dirstruct_table, dirtable_dir)
%LOAD_AND_UPDATE_TABLES Load existing tables and merge new data.

    try
        % Load 'reference' sheet
        opts = detectImportOptions(dirtable_dir, 'Sheet', 'reference');
        % opts.VariableTypes = repmat({'string'}, 1, numel(opts.VariableNames)); 
        ref_table = readtable(dirtable_dir, opts);

        % Load 'log' sheet (create if missing)
        try
            opts_log = detectImportOptions(dirtable_dir, 'Sheet', 'log');
            log_table = readtable(dirtable_dir, opts_log);
        catch
            % If log sheet likely doesn't exist, initialize it based on ref_table
            log_table = ref_table;
            var_names = log_table.Properties.VariableNames;
            for i = 1:numel(var_names)
                log_table.(var_names{i}) = repmat(datetime("now"), height(log_table), 1);
            end
        end
        
        % Define Columns
        meta_cols = {'VesselID', 'Power', 'Depth', 'Notes', 'Wavelength', 'xlocation', 'ylocation', 'zlocation', 'PercentPower'};
        all_cols = dirstruct_table.Properties.VariableNames;
        file_cols = setdiff(all_cols, meta_cols);

        % Iterate and Update
        for i = 1:height(dirstruct_table)
            new_row = dirstruct_table(i, :);
            [ref_table, log_table] = process_session_update(ref_table, log_table, new_row, file_cols);
        end

        % Sort by Date
        if ismember('Date', ref_table.Properties.VariableNames)
            [ref_table, sort_idx] = sortrows(ref_table, 'Date');
            log_table = log_table(sort_idx, :); % Keep log in sync
        end

        % Save
        writetable(ref_table, dirtable_dir, "Sheet", "reference");
        writetable(dirstruct_table, dirtable_dir, "Sheet", "auto_mapped");
        writetable(log_table, dirtable_dir, "Sheet", "log");
        
    catch ME
        warning('Failed to update directory table: %s', ME.message);
    end
end

function [ref_table, log_table] = process_session_update(ref_table, log_table, new_row, file_cols)
%PROCESS_SESSION_UPDATE Update single session row in reference and log tables.

    key = new_row.Directory; % Unique Key
    
    % Find match in reference table
    if iscell(key)
        key = key{1};
    end
    
    match_idx = find(strcmpi(ref_table.Directory, key), 1);

    if ~isempty(match_idx)
        % --- EXISTING SESSION ---
        % Update only file columns if changed
        for c = 1:numel(file_cols)
            col_name = file_cols{c};
            if ismember(col_name, ref_table.Properties.VariableNames)
                
                old_val = ref_table.(col_name)(match_idx);
                new_val = new_row.(col_name);
                
                % Handle cell/string/numeric comparison
                has_changed = false;
                if iscell(old_val) || isstring(old_val)
                    if ~strcmp(string(old_val), string(new_val))
                        has_changed = true;
                    end
                elseif isnumeric(old_val)
                    if old_val ~= new_val
                        has_changed = true;
                    end
                else
                    % Fallback for other types
                    if ~isequal(old_val, new_val)
                         has_changed = true;
                    end
                end

                if has_changed
                    ref_table.(col_name)(match_idx) = new_val;
                    log_table.(col_name)(match_idx) = datetime("now");
                    fprintf('Updated: %s [%s]\n', key, col_name);
                end
            end
        end
    else
        % --- NEW SESSION ---
        % Append to ref_table
        ref_table = append_row(ref_table, new_row);
        
        % Create new log row with current timestamps
        new_log_row = new_row;
        var_names = new_log_row.Properties.VariableNames;
        for i = 1:numel(var_names)
             new_log_row.(var_names{i}) = datetime("now");
        end
        log_table = append_row(log_table, new_log_row);
        
        fprintf('New Session Added: %s\n', key);
    end
end

function tbl = append_row(tbl, new_row)
%APPEND_ROW Helper to robustly append a row, handling missing/extra columns.

    % Add missing columns to table
    missing_cols = setdiff(new_row.Properties.VariableNames, tbl.Properties.VariableNames);
    if ~isempty(missing_cols)
         for mc = 1:numel(missing_cols)
             tbl.(missing_cols{mc}) = repmat("", height(tbl), 1); % Default empty string
         end
    end

    % Add missing columns to new_row
    missing_in_new = setdiff(tbl.Properties.VariableNames, new_row.Properties.VariableNames);
     if ~isempty(missing_in_new)
         for mc = 1:numel(missing_in_new)
             new_row.(missing_in_new{mc}) = ""; % Default empty string
         end
    end
    
    % Reorder and append
    new_row = new_row(:, tbl.Properties.VariableNames);
    tbl = [tbl; new_row];
end
