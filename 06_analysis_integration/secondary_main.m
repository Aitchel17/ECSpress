clc, clear
analysis_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis';
experiment_folder = 'G:\tmp\00_igkl';

% metadata table filtering
% load table
[~, tmp.exp_name] = fileparts(experiment_folder);
save_exppath = fullfile(analysis_path,tmp.exp_name);
tmp.dirtable_dir = fullfile(save_exppath,strcat(tmp.exp_name,'_dirtable.xlsx'));
tmp.opts = detectImportOptions(tmp.dirtable_dir, 'Sheet', 'reference');
ref_table = readtable(tmp.dirtable_dir, tmp.opts);

% Filter for Sleep Sessions
tmp.rows = ref_table.LineFWHM == "paxfwhm.mat";
fwhm_table = ref_table(tmp.rows, :);
tmp.rows = fwhm_table.SessionType == "sleep";
sleep_fwhm_table = fwhm_table(tmp.rows,:);

% 1. Aggregate FWHM Data
% Helper function located in 00_mapdirectorystruct
disp('Aggregating FWHM Data...');
master_fwhm_table = aggregate_fwhm(sleep_fwhm_table);
% Amend table for depth logic
% 1. Parse Depth Logic securely
raw_depth = master_fwhm_table.Depth;
numeric_depth = zeros(height(master_fwhm_table), 1);
for i = 1:height(master_fwhm_table)
    d_val = raw_depth(i);
    tmp_str = strsplit(d_val,'um');
    tmp_str = tmp_str(1);
    numeric_depth(i) = str2double(tmp_str);  
end
depth_logic = numeric_depth > 70;
depth_name = ["L1", "L2"];
depthn_array = depth_name(depth_logic+1);
master_fwhm_table.Depthstate = depthn_array';

save(fullfile(save_exppath,"sleep_paxfwhm_table.mat"),"master_fwhm_table")





























%% 2. Split into L1 and L2 Tables
% L1: <= 70um (Shallow)
% L2: > 70um (Deep)
L1_logic = master_fwhm_table.NumericDepth <= 70;
L2_logic = master_fwhm_table.NumericDepth > 70;

L1_data = master_fwhm_table(L1_logic, :);
L2_data = master_fwhm_table(L2_logic, :);

disp(['L1 Sets: ', num2str(height(L1_data))]);
disp(['L2 Sets: ', num2str(height(L2_data))]);

%% 3. Group and Aggregate (Average Multi-Session Vessels)
disp('Grouping L1 Vessels...');
L1_Vessels_Struct = group_vessels(L1_data);

disp('Grouping L2 Vessels...');
L2_Vessels_Struct = group_vessels(L2_data);


%% Combine the same vessel
target_table = L1_Vessels_Struct(2).RawTable;
trans_logic = target_table.State_Transition(2).thickness_bv.state_name == "ra_trans";
data = target_table.State_Transition(2).thickness_bv(trans_logic,:);
%%




%%
function vessel_struct = group_vessels(input_table)
    unique_mice = unique(lower(input_table.MouseID));
    vessel_struct = struct([]);
    count = 0;
    
    for m = 1:numel(unique_mice)
        mid = unique_mice(m);
        m_rows = input_table(lower(input_table.MouseID) == mid, :);
        
        unique_vids = unique(m_rows.VesselID);
        for v = 1:numel(unique_vids)
            vid = unique_vids(v);
            v_rows = m_rows(m_rows.VesselID == vid, :);
            
            count = count + 1;
            vessel_struct(count).MouseID = mid;
            vessel_struct(count).VesselID = vid;
            vessel_struct(count).AvgDepth = mean(v_rows.NumericDepth);
            vessel_struct(count).SessionCount = height(v_rows);
            vessel_struct(count).RawTable = v_rows; % Store all sessions for this vessel
            
            % Setup for Averaging (Placeholder for now)
            % Here you would average time series or transition metrics
            % e.g. vessel_struct(count).AvgTransition = ...
        end
    end
end
