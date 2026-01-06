function pax_cluster = analysis_clusterpolar(pax_fwhm, twophoton_processed, output_dir)
% ANALYSIS_PAX_CLUSTER Performs cluster analysis on PAX FWHM data.
%   1. Clusters the kymograph based on pairwise correlation.
%   2. Generates median filtered images for each cluster.
%   3. Prompts user to identify constricted/dilated clusters.
%   4. Saves output to disk.

%% 1. Initialize parameter
pax_cluster.cluster_num = 20;

%% 2. Make cluster
[pax_cluster.kgph_lumen_columnidx, pax_cluster.clusterboundary] = ...
    analysis_cluster_kymograph(pax_fwhm.kymograph.kgph_lumen_processed,pax_cluster.cluster_num);

%% 3. Make median filter of images at each cluster
pax_cluster.clust_med_bv = zeros([size(twophoton_processed.ch1,1,2),pax_cluster.cluster_num]);
pax_cluster.clust_med_csf = zeros([size(twophoton_processed.ch2,1,2),pax_cluster.cluster_num]);

tmp.rearanged_bvstack = twophoton_processed.ch1(:,:,pax_cluster.kgph_lumen_columnidx);
tmp.rearanged_csfstack = twophoton_processed.ch2(:,:,pax_cluster.kgph_lumen_columnidx);

for cluster_idx = 1:pax_cluster.cluster_num
    tmp.clusterstart = pax_cluster.clusterboundary(cluster_idx,1);
    tmp.clusterend = pax_cluster.clusterboundary(cluster_idx,2);
    pax_cluster.clust_med_bv(:,:,cluster_idx) = median(tmp.rearanged_bvstack(:,:,tmp.clusterstart:tmp.clusterend),3);
    pax_cluster.clust_med_csf(:,:,cluster_idx) = median(tmp.rearanged_csfstack(:,:,tmp.clusterstart:tmp.clusterend),3);
end

% 4. Make summary
tmp.bv_diameter = pax_fwhm.idx.lowerBVboundary-pax_fwhm.idx.upperBVboundary;
tmp.lowidxchange = pax_fwhm.idx.lowerBVboundary - median(pax_fwhm.idx.lowerBVboundary);
tmp.upidxchange = median(pax_fwhm.idx.upperBVboundary) -  pax_fwhm.idx.upperBVboundary;
tmp.clustersort_bv_diameter = tmp.bv_diameter(pax_cluster.kgph_lumen_columnidx);
tmp.clustersort_lowidxchange = tmp.lowidxchange(pax_cluster.kgph_lumen_columnidx);
tmp.clustersort_upidxchange  = tmp.upidxchange(pax_cluster.kgph_lumen_columnidx);

% 1. mean, 2. std, 3. length, 4. bv_lowidx changes, 5. bv_upidx change
% 1. mean, 2. std, 3. length, 4. bv_lowidx changes, 5. bv_upidx change
pax_cluster.cluster_summary = struct('id',{},'mean_diameter',{},'std_diameter',{},'duration',{},'mean_lowchange',{},'mean_upchange',{});
for cluster_idx = 1: size(pax_cluster.clusterboundary,1)
    pax_cluster.cluster_summary(cluster_idx).id = cluster_idx;
    pax_cluster.cluster_summary(cluster_idx).mean_diameter = mean(tmp.clustersort_bv_diameter(pax_cluster.clusterboundary(cluster_idx,1):pax_cluster.clusterboundary(cluster_idx,2)));
    pax_cluster.cluster_summary(cluster_idx).std_diameter = std(tmp.clustersort_bv_diameter(pax_cluster.clusterboundary(cluster_idx,1):pax_cluster.clusterboundary(cluster_idx,2)));
    pax_cluster.cluster_summary(cluster_idx).duration = pax_cluster.clusterboundary(cluster_idx,2)-pax_cluster.clusterboundary(cluster_idx,1)+1;
    pax_cluster.cluster_summary(cluster_idx).mean_lowchange = mean(tmp.clustersort_lowidxchange(pax_cluster.clusterboundary(cluster_idx,1):pax_cluster.clusterboundary(cluster_idx,2)));
    pax_cluster.cluster_summary(cluster_idx).mean_upchange = mean(tmp.clustersort_upidxchange (pax_cluster.clusterboundary(cluster_idx,1):pax_cluster.clusterboundary(cluster_idx,2)));
end

%% 5. Filter cluster (length > 10 frames)
% Filter length below 10 frames 5fps, 2 sec average
pax_cluster.filtered_clusteridx = [pax_cluster.cluster_summary.id];
pax_cluster.filtered_clusteridx = pax_cluster.filtered_clusteridx([pax_cluster.cluster_summary.duration]>10);

%% 6. Manually find constricted dilated cluster (Only from filtered clusters)
% Pass filtered stack to util_checkstack
filtered_csf_stack = pax_cluster.clust_med_csf(:,:,pax_cluster.filtered_clusteridx);
temp_const_idx = util_checkstack(filtered_csf_stack, 'Select the constricted image.');
temp_dil_idx = util_checkstack(filtered_csf_stack, 'Select the dilated image.');

% Map back to original cluster index
pax_cluster.constricted_cluster_idx = pax_cluster.filtered_clusteridx(temp_const_idx);
pax_cluster.dilated_cluster_idx = pax_cluster.filtered_clusteridx(temp_dil_idx);

pax_cluster.constricted_medianimg = pax_cluster.clust_med_csf(:,:,pax_cluster.constricted_cluster_idx);
pax_cluster.dilated_medianimg = pax_cluster.clust_med_csf(:,:,pax_cluster.dilated_cluster_idx);

%% 7. (x,y,bv/csf,constricted/dilated)
pax_cluster.cluster_bvcsf_constdil = zeros([size(tmp.rearanged_bvstack,1,2),2,2]); % (x,y,bv/csf,constricted/dilated)
pax_cluster.cluster_bvcsf_constdil(:,:,1,1) = pax_cluster.clust_med_bv(:,:,pax_cluster.constricted_cluster_idx);
pax_cluster.cluster_bvcsf_constdil(:,:,2,1) = pax_cluster.clust_med_csf(:,:,pax_cluster.constricted_cluster_idx);
pax_cluster.cluster_bvcsf_constdil(:,:,1,2) = pax_cluster.clust_med_bv(:,:,pax_cluster.dilated_cluster_idx);
pax_cluster.cluster_bvcsf_constdil(:,:,2,2) = pax_cluster.clust_med_csf(:,:,pax_cluster.dilated_cluster_idx);

%% Output
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
save(fullfile(output_dir, 'pax_cluster.mat'),"pax_cluster")

end
