clc, clear
masterDirTable_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl\00_igkl_dirtable.xlsx';
mtable_FWHMsleep = tableManager.load_recon(masterDirTable_path, "mtable_FWHMsleep.mat");
mtable_FWHMsleep.analysis_table = mtable_FWHMsleep.subTables.transition; % State summary to be target table
%%
mtable_FWHMsleep.filtLogics = [];
mtable_FWHMsleep.filtLogics.vType = contains(mtable_FWHMsleep.analysis_table.VesselID, "PA", 'IgnoreCase', true); % Artery filter
mtable_FWHMsleep.filtLogics.depth = mtable_FWHMsleep.analysis_table.NumericDepth <70; % L1 filter
mtable_FWHMsleep.apply_filter

%%
data_colnames = {"data"};
numeric_colnames = {'pre_mean','pre_median','pre_q1','pre_q3', 'pre_var',...
                'post_mean','post_median','post_q1','post_q3', 'post_var'};

myTransAnalyzer = TableAnalyzer(mtable_FWHMsleep.filtered_table, mtable_FWHMsleep.action_log);

myTransAnalyzer.scale_table("NumericResolution",data_colnames,numeric_colnames);

%%


myTransAnalyzer.get_numericsummary("Date","filtered_table")


%% Mouse table
ra_transitions = [];

unique_mids = unique(l1_table.MouseID);
for midx = 1:size(unique_mids,1)
    mouse_name = unique_mids(midx);
    mid_logic = l1_table.MouseID == mouse_name;
    mouse_table = l1_table(mid_logic,:);
    % vessel table
    unique_vids = unique(mouse_table.VesselID);
    n_vids = length(unique_vids);
    for vidx = 1:n_vids
        vessel_id = unique_vids(vidx);
        if contains(vessel_id, "PA")
        
        vid_logic = mouse_table.VesselID == vessel_id;
        vessel_table = mouse_table(vid_logic,:);
        
        % Get all RA transition
        % merge session
        merged_table = [];
        for idx = 1:size(vessel_table.State_Transition,1)
            trans_table = vessel_table.State_Transition(idx).thickness_totalpvs;
            merged_table = [merged_table ; trans_table];
        end
        
        % RA transition
        trans_logic = merged_table.state_name == 'ra_trans';
        ra_table = merged_table(trans_logic,:);
        averaged_ratransition = [];
        n_trans = size(ra_table.data,1);
        for trans_idx = 1:n_trans
            averaged_ratransition = [averaged_ratransition;ra_table.data{trans_idx}];
        end
        ra_transition = mean(averaged_ratransition,1);
        ra_transitions = [ra_transitions; ra_transition];
        end
    end
end



%%

d1 = ra_table.data{1};
d2 = ra_table.data{2};
%%

transition_length = (size(ra_transitions,2)-1)/2;

pre_transit = ra_transitions(:,1:transition_length);

max_transit = max(pre_transit,[],2);

norm_transit = ra_transitions./max_transit;

%%
cla
plot(mean(norm_transit,1))
hold on
for i = 1:size(ra_transitions,1)
    plot(norm_transit(i,:),Color=[0.7 0.7 0.7])
end