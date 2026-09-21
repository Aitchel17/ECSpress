function polarcluster = analysis_clusterpolar(pax_fwhm, twophoton_processed, output_dir)
% ANALYSIS_PAX_CLUSTER Performs cluster analysis on PAX FWHM data.
%   1. Clusters the kymograph based on pairwise correlation.
%   2. Generates median filtered images for each cluster.
%   3. Prompts user to identify constricted/dilated clusters.
%   4. Saves output to disk.

%% 1. Initialize parameter
polarcluster.cluster_num = 20;

%% 2. Make cluster
[polarcluster.kgph_lumen_columnidx, polarcluster.clusterboundary] = ...
    analysis_cluster_kymograph(pax_fwhm.kymograph.kgph_lumen_processed,polarcluster.cluster_num);

%% 3. Make median filter of images at each cluster
% Follow the channel FWHM used (recorded on pax_fwhm.param by addkymograph);
% default BV=ch1, CSF=ch2 when not recorded.
bvch = 'ch1';  csfch = 'ch2';
if isfield(pax_fwhm.param, 'channel_lumen'); bvch  = char(pax_fwhm.param.channel_lumen); end
if isfield(pax_fwhm.param, 'channel_pvs');   csfch = char(pax_fwhm.param.channel_pvs);   end
polarcluster.clust_med_bv  = zeros([size(twophoton_processed.(bvch),1,2),polarcluster.cluster_num]);
polarcluster.clust_med_csf = zeros([size(twophoton_processed.(csfch),1,2),polarcluster.cluster_num]);

tmp.rearanged_bvstack  = twophoton_processed.(bvch)(:,:,polarcluster.kgph_lumen_columnidx);
tmp.rearanged_csfstack = twophoton_processed.(csfch)(:,:,polarcluster.kgph_lumen_columnidx);

for cluster_idx = 1:polarcluster.cluster_num
    tmp.clusterstart = polarcluster.clusterboundary(cluster_idx,1);
    tmp.clusterend = polarcluster.clusterboundary(cluster_idx,2);
    polarcluster.clust_med_bv(:,:,cluster_idx) = median(tmp.rearanged_bvstack(:,:,tmp.clusterstart:tmp.clusterend),3);
    polarcluster.clust_med_csf(:,:,cluster_idx) = median(tmp.rearanged_csfstack(:,:,tmp.clusterstart:tmp.clusterend),3);
end

% 4. Make summary
tmp.bv_diameter = pax_fwhm.idx.lowerBVboundary-pax_fwhm.idx.upperBVboundary;
tmp.lowidxchange = pax_fwhm.idx.lowerBVboundary - median(pax_fwhm.idx.lowerBVboundary);
tmp.upidxchange = median(pax_fwhm.idx.upperBVboundary) -  pax_fwhm.idx.upperBVboundary;
tmp.clustersort_bv_diameter = tmp.bv_diameter(polarcluster.kgph_lumen_columnidx);
tmp.clustersort_lowidxchange = tmp.lowidxchange(polarcluster.kgph_lumen_columnidx);
tmp.clustersort_upidxchange  = tmp.upidxchange(polarcluster.kgph_lumen_columnidx);

% 1. mean, 2. std, 3. length, 4. bv_lowidx changes, 5. bv_upidx change
% 1. mean, 2. std, 3. length, 4. bv_lowidx changes, 5. bv_upidx change
polarcluster.cluster_summary = struct('id',{},'mean_diameter',{},'std_diameter',{},'duration',{},'mean_lowchange',{},'mean_upchange',{});
for cluster_idx = 1: size(polarcluster.clusterboundary,1)
    polarcluster.cluster_summary(cluster_idx).id = cluster_idx;
    polarcluster.cluster_summary(cluster_idx).mean_diameter = mean(tmp.clustersort_bv_diameter(polarcluster.clusterboundary(cluster_idx,1):polarcluster.clusterboundary(cluster_idx,2)));
    polarcluster.cluster_summary(cluster_idx).std_diameter = std(tmp.clustersort_bv_diameter(polarcluster.clusterboundary(cluster_idx,1):polarcluster.clusterboundary(cluster_idx,2)));
    polarcluster.cluster_summary(cluster_idx).duration = polarcluster.clusterboundary(cluster_idx,2)-polarcluster.clusterboundary(cluster_idx,1)+1;
    polarcluster.cluster_summary(cluster_idx).mean_lowchange = mean(tmp.clustersort_lowidxchange(polarcluster.clusterboundary(cluster_idx,1):polarcluster.clusterboundary(cluster_idx,2)));
    polarcluster.cluster_summary(cluster_idx).mean_upchange = mean(tmp.clustersort_upidxchange (polarcluster.clusterboundary(cluster_idx,1):polarcluster.clusterboundary(cluster_idx,2)));
end

%% 5. Filter cluster (length > 10 frames)
% Filter length below 10 frames 5fps, 2 sec average
polarcluster.filtered_clusteridx = [polarcluster.cluster_summary.id];
polarcluster.filtered_clusteridx = polarcluster.filtered_clusteridx([polarcluster.cluster_summary.duration]>10);

%% 6. Manually find constricted dilated cluster (Only from filtered clusters)
% Pass filtered stack to util_checkstack
filtered_csf_stack = polarcluster.clust_med_csf(:,:,polarcluster.filtered_clusteridx);
temp_const_idx = util_checkstack(filtered_csf_stack, 'Select the constricted image.');
temp_dil_idx = util_checkstack(filtered_csf_stack, 'Select the dilated image.');

% Map back to original cluster index
polarcluster.constricted_cluster_idx = polarcluster.filtered_clusteridx(temp_const_idx);
polarcluster.dilated_cluster_idx = polarcluster.filtered_clusteridx(temp_dil_idx);

polarcluster.constricted_medianimg = polarcluster.clust_med_csf(:,:,polarcluster.constricted_cluster_idx);
polarcluster.dilated_medianimg = polarcluster.clust_med_csf(:,:,polarcluster.dilated_cluster_idx);

%% 6b. Intermediate cluster, chosen rather than picked
% Two states fix one area ratio; a third says whether the contour area follows the FWHM
% diameter squared. It is the cluster nearest the median diameter, among those whose median
% image averages at least as many frames as the thinner of the two picks.
cluster_diameter = [polarcluster.cluster_summary.mean_diameter];
cluster_duration = [polarcluster.cluster_summary.duration];
picked_idx = [polarcluster.constricted_cluster_idx, polarcluster.dilated_cluster_idx];
frame_floor = min(cluster_duration(picked_idx));
target_diameter = median(cluster_diameter(polarcluster.filtered_clusteridx));
is_long_enough = cluster_duration(polarcluster.filtered_clusteridx) >= frame_floor;
is_unpicked = ~ismember(polarcluster.filtered_clusteridx, picked_idx);
keep_candidate = is_long_enough & is_unpicked;
candidate_idx = polarcluster.filtered_clusteridx(keep_candidate);
[~, nearest] = min(abs(cluster_diameter(candidate_idx) - target_diameter));
polarcluster.intermediate_cluster_idx = candidate_idx(nearest);
fprintf('intermediate cluster : %d of %d candidates, %d frames, %.2f px against a median of %.2f\n', ...
    polarcluster.intermediate_cluster_idx, numel(candidate_idx), ...
    cluster_duration(polarcluster.intermediate_cluster_idx), ...
    cluster_diameter(polarcluster.intermediate_cluster_idx), target_diameter);

%% 7. (x,y,bv/csf,state), state 1 constricted, 2 intermediate, 3 dilated
polarcluster.state_cluster_idx = [polarcluster.constricted_cluster_idx, ...
    polarcluster.intermediate_cluster_idx, polarcluster.dilated_cluster_idx];
polarcluster.cluster_bvcsf_state = zeros([size(tmp.rearanged_bvstack,1,2),2,numel(polarcluster.state_cluster_idx)]);
for state_idx = 1:numel(polarcluster.state_cluster_idx)
    cluster_idx = polarcluster.state_cluster_idx(state_idx);
    polarcluster.cluster_bvcsf_state(:,:,1,state_idx) = polarcluster.clust_med_bv(:,:,cluster_idx);
    polarcluster.cluster_bvcsf_state(:,:,2,state_idx) = polarcluster.clust_med_csf(:,:,cluster_idx);
end

%% Output
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save(fullfile(output_dir, 'polarcluster.mat'),"polarcluster")

end
