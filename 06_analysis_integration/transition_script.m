clc, clear
analysis_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis';
experiment_folder = 'G:\tmp\00_igkl';
target_transition = 'sleep_paxfwhm_table.mat';
[~, tmp.exp_name] = fileparts(experiment_folder);
save_exppath = fullfile(analysis_path,tmp.exp_name);
load(fullfile(save_exppath,target_transition))
%% L1 Table
depth_logic = master_fwhm_table.Depthstate == "L1";
l1_table = master_fwhm_table(depth_logic,:);

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