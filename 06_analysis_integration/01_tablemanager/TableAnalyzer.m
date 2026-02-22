classdef TableAnalyzer < handle
    % TABLEANALYZER Specialized class for calculating summaries and applying mathematical logic to data tables.
    % Receives the managed tables from a TableManager.
    properties
        filtered_table
        action_log
        numeric_tables
        data_tables
    end
    
    methods
        function obj = TableAnalyzer(input_table, input_log)
            % Initializes TableAnalyzer with data loaded/organized by TableManager
            obj.filtered_table = input_table;
            if nargin > 1
                obj.action_log = input_log;
            else
                obj.action_log = struct(); % Initialize if empty
            end
        end
        
        function meanFrom2(obj, datacolName, NewcolName, start_fraction, end_fraction)
            dlength = height(obj.filtered_table);
            numeric_arrs = obj.filtered_table.(datacolName);
            result_col = zeros(dlength, 1);
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
            fprintf("mean from %s to %s calculated and added as %s column\n", start_fraction, end_fraction, NewcolName)
        end
 
        function addPrctilecol(obj, datacolName, NewcolName, percentile)
            dlength = height(obj.filtered_table);
            numeric_arrs = obj.filtered_table.(datacolName);
            result_col = zeros(dlength, 1);
            for arr_idx = 1:dlength
                numeric_arr = numeric_arrs{arr_idx};
                result_col(arr_idx) = prctile(numeric_arr, percentile);  
            end
            obj.filtered_table.(NewcolName) = result_col;
            
            current_cols = string(obj.action_log.numeric_colnames);
            chk_duplication = ismember(current_cols, string(NewcolName));
            if ~any(chk_duplication)
                obj.action_log.numeric_colnames = [current_cols, string(NewcolName)];
            end
            fprintf("%d Percentile calculated and added as %s column\n", percentile, NewcolName)
        end


        function scale_table(obj, scale_colname, data_colnames, numeric_colnames)
            % Applies resolution scaling (multiplication) to specified columns
            scale_logname = strcat(scale_colname, "_applied");
            if isfield(obj.action_log, scale_logname)
                fprintf("%s already applied\n", scale_colname)
            else
                scale_vec = obj.filtered_table.(scale_colname);
                numeric_data = obj.filtered_table{:, numeric_colnames};
                obj.filtered_table{:, numeric_colnames} = numeric_data .* scale_vec;
                % 2. Scale Data Columns (Cell arrays of vectors)
                
                for i = 1:numel(data_colnames)
                    col_name = data_colnames{i};
                    % Loop over rows
                    current_col = obj.filtered_table.(col_name);
                    for j = 1:height(obj.filtered_table)
                        % Multiply content
                        if ~isempty(current_col{j})
                            current_col{j} = current_col{j} * scale_vec(j);
                        end
                    end
                    obj.filtered_table.(col_name) = current_col;                
                end
                obj.action_log.(scale_logname) = true;
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

            temp_table = stats(:, [summary_key, "mean_" + var_names]);           
            
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
            var_name = obj.action_log.data_colnames{1};

            [target_table, summary_key, aveTable_name, ci95Table_name] = obj.prepare_summary("data_tables", target_tablename, new_keyname);

            % Find unique group combinations
            group_table = target_table(:, summary_key);
            [unique_groups, ~, group_idx] = unique(group_table, 'rows', 'stable');
            n_groups = height(unique_groups);

            ave_table = unique_groups;
            ci_table  = unique_groups;
            ave_table.(var_name) = cell(n_groups, 1);
            ci_table.(var_name)  = cell(n_groups, 1);
            col_data = target_table.(var_name);
            for g = 1:n_groups
                row_idx = find(group_idx == g);
                n = numel(row_idx);
                
                if n == 1                     
                    ave_table.(var_name){g} = col_data{row_idx};
                    ci_table.(var_name){g}  = zeros(size(col_data{row_idx}));
                    continue;
                end

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
            df = max(n - 1, 1); 
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

        function [target_table, summary_key] = retrieve_keytable(obj, table_typename, target_tablename, new_keyname)
            newkey_name = strcat(new_keyname, "_key");
            if strcmp(target_tablename, 'filtered_table')
                target_table = obj.filtered_table;
                key_names = obj.action_log.key_names;
            else
                target_table = obj.(table_typename).(target_tablename);
                actionLog_keyfield = strsplit(target_tablename, "_");
                actionLog_keyfield = strjoin(actionLog_keyfield(1:end-1), "_");
                actionLog_keyfield = strcat(actionLog_keyfield, "_key");
                key_names = obj.action_log.(actionLog_keyfield);
            end
            logic_key = ismember(key_names, new_keyname);
            if  any(logic_key)
                summary_key = key_names(~logic_key);
                obj.action_log.(newkey_name) = summary_key;
            else
                fprintf("target_keyname is not in key_names at actionlog\n")
            end
        end
    end
end
