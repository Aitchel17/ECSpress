% MAKE CLUSTER
% 1. Make cluster from pax kymograph based on pairwise correlation distance
% 2. Using cluster idx, make median filtered representation of cluster
% 3. Manually find clean representative constricted and dilated cluster

%% 1. Make cluster
%%
if isfield(primary_datastruct,'pax_cluster')
    pax_cluster = primary_datastruct.pax_cluster;
end
%%
pax_cluster.cluster_num = 20;
[pax_cluster.kgph_lumen_columnidx, pax_cluster.clusterboundary] = analysis_cluster_kymograph(pax_fwhm.kymograph.kgph_lumen_processed,pax_cluster.cluster_num);
% 2. Make median filter of images at each cluster
pax_cluster.clust_med_bv = zeros([size(roianalysis.preprocessed_ch1,1,2),pax_cluster.cluster_num]);
pax_cluster.clust_med_csf = zeros([size(roianalysis.preprocessed_ch2,1,2),pax_cluster.cluster_num]);
tmp.rearanged_bvstack = roianalysis.preprocessed_ch1(:,:,pax_cluster.kgph_lumen_columnidx);
tmp.rearanged_csfstack = roianalysis.preprocessed_ch2(:,:,pax_cluster.kgph_lumen_columnidx);
for cluster_idx = 1:pax_cluster.cluster_num
    tmp.clusterstart = pax_cluster.clusterboundary(cluster_idx,1);
    tmp.clusterend = pax_cluster.clusterboundary(cluster_idx,2);
    pax_cluster.clust_med_bv(:,:,cluster_idx) = median(tmp.rearanged_bvstack(:,:,tmp.clusterstart:tmp.clusterend),3);
    pax_cluster.clust_med_csf(:,:,cluster_idx) = median(tmp.rearanged_csfstack(:,:,tmp.clusterstart:tmp.clusterend),3);
end
% 3. Make summary
tmp.bv_diameter = pax_fwhm.idx.lowerboundary-pax_fwhm.idx.upperboundary;
tmp.lowidxchange = pax_fwhm.idx.lowerboundary - median(pax_fwhm.idx.lowerboundary);
tmp.upidxchange = median(pax_fwhm.idx.upperboundary) -  pax_fwhm.idx.upperboundary;
tmp.clustersort_bv_diameter = tmp.bv_diameter(pax_cluster.kgph_lumen_columnidx);
tmp.clustersort_lowidxchange = tmp.lowidxchange(pax_cluster.kgph_lumen_columnidx);
tmp.clustersort_upidxchange  = tmp.upidxchange(pax_cluster.kgph_lumen_columnidx);
% 1. mean, 2. std, 3. length, 4. bv_lowidx changes, 5. bv_upidx change
pax_cluster.cluster_summary = zeros([size(pax_cluster.clusterboundary,1),6]);
for cluster_idx = 1: size(pax_cluster.clusterboundary,1)
    pax_cluster.cluster_summary(cluster_idx,1) = cluster_idx;
    pax_cluster.cluster_summary(cluster_idx,2) = mean(tmp.clustersort_bv_diameter(pax_cluster.clusterboundary(cluster_idx,1):pax_cluster.clusterboundary(cluster_idx,2)));
    pax_cluster.cluster_summary(cluster_idx,3) = std(tmp.clustersort_bv_diameter(pax_cluster.clusterboundary(cluster_idx,1):pax_cluster.clusterboundary(cluster_idx,2)));
    pax_cluster.cluster_summary(cluster_idx,4) = pax_cluster.clusterboundary(cluster_idx,2)-pax_cluster.clusterboundary(cluster_idx,1)+1;
    pax_cluster.cluster_summary(cluster_idx,5) = mean(tmp.clustersort_lowidxchange(pax_cluster.clusterboundary(cluster_idx,1):pax_cluster.clusterboundary(cluster_idx,2)));
    pax_cluster.cluster_summary(cluster_idx,6) = mean(tmp.clustersort_upidxchange (pax_cluster.clusterboundary(cluster_idx,1):pax_cluster.clusterboundary(cluster_idx,2)));
end

%% 3. Manually find constricted dilated cluster
pax_cluster.constricted_cluster_idx = util_checkstack(pax_cluster.clust_med_csf);
%%
pax_cluster.dilated_cluster_idx = util_checkstack(pax_cluster.clust_med_csf);
%% (x,y,bv/csf,constricted/dilated)
pax_cluster.cluster_bvcsf_constdil = zeros([size(tmp.rearanged_bvstack,1,2),2,2]); % (x,y,bv/csf,constricted/dilated)
pax_cluster.cluster_bvcsf_constdil(:,:,1,1) = pax_cluster.clust_med_bv(:,:,pax_cluster.constricted_cluster_idx);
pax_cluster.cluster_bvcsf_constdil(:,:,2,1) = pax_cluster.clust_med_csf(:,:,pax_cluster.constricted_cluster_idx);
pax_cluster.cluster_bvcsf_constdil(:,:,1,2) = pax_cluster.clust_med_bv(:,:,pax_cluster.dilated_cluster_idx);
pax_cluster.cluster_bvcsf_constdil(:,:,2,2) = pax_cluster.clust_med_csf(:,:,pax_cluster.dilated_cluster_idx);

% Filter length below 10 frames 5fps, 2 sec average (opt... may be I need range of vascular dilation for 
% PIV --> in case porosity change only occur during initial dilation)
pax_cluster.filtered_clusteridx = pax_cluster.cluster_summary(:,1);
pax_cluster.filtered_clusteridx = pax_cluster.filtered_clusteridx(pax_cluster.cluster_summary(:,4)>10);
%%
save(fullfile(directories.save_dir,'pax_cluster.mat'),"pax_cluster")

%% Check sorting result
fig.paxsortedbv_fig = make_fig('paxsortedBV_figure');
%%
fig.paxsortedbv_fig.bring_fig
fig.paxsortedbv_fig.update_figsize([8 3])
fig.paxsortedbv_fig.reset_axis()
fig.paxsortedbv_fig.resolution = primary_datastruct.img_param.pixel2um;

fig.paxsortedbv_fig.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed(:,pax_cluster.kgph_lumen_columnidx),roianalysis.taxis)
fig.paxsortedbv_fig.plot_line(pax_fwhm.idx.clean_upperboundary(pax_cluster.kgph_lumen_columnidx),'r');
fig.paxsortedbv_fig.plot_line(pax_fwhm.idx.clean_lowerboundary(pax_cluster.kgph_lumen_columnidx),'r');
% fig.paxsortedbv_fig.plot_line(pax_fwhm.idx.clean_pvsupedge_idx(pax_cluster.kgph_lumen_columnidx),'g');
% fig.paxsortedbv_fig.plot_line(pax_fwhm.idx.clean_pvsdownedge_idx(pax_cluster.kgph_lumen_columnidx),'g');
fig.paxsortedbv_fig.plot_xline(pax_cluster.clusterboundary(:,2),'m')
fig.paxsortedbv_fig.put_yaxistitle('Length (\mum)')
fig.paxsortedbv_fig.put_xaxistitle('Total duration (sec)')
text_x_dilate =fig.paxsortedbv_fig.loc.x(pax_cluster.clusterboundary(pax_cluster.dilated_cluster_idx));
text_x_const =fig.paxsortedbv_fig.loc.x(pax_cluster.clusterboundary(pax_cluster.constricted_cluster_idx));

text_y =fig.paxsortedbv_fig.loc.y(end)+1;

pax_cluster.clusterboundary(pax_cluster.constricted_cluster_idx);
dilatedcluster_startframe = pax_cluster.clusterboundary(pax_cluster.dilated_cluster_idx);
text(fig.paxsortedbv_fig.ax,text_x_const,text_y,'Constricted')
text(fig.paxsortedbv_fig.ax,text_x_dilate,text_y,'Dilated')
fig.paxsortedbv_fig.save2svg(directories.save_dir)
%%

brain_distortionstruct = sturct();



%%
fig.cluster_cbv = make_fig('cluster_constricted_bv');
fig.cluster_cbv.showimg(pax_cluster.cluster_bvcsf_constdil(:,:,1,1))
%%
fig.cluster_cbv.save2svg(directories.save_dir)
%%
fig.cluster_dbv = make_fig('cluster_dilated_bv');
fig.cluster_dbv.showimg(pax_cluster.cluster_bvcsf_constdil(:,:,1,2))
%%
fig.cluster_dbv.save2svg(directories.save_dir)
%%
fig.cluster_cpvs = make_fig('cluster_constricted_pvs');
fig.cluster_cpvs.showimg(pax_cluster.cluster_bvcsf_constdil(:,:,2,1))
%%
fig.cluster_cpvs.save2svg(directories.save_dir)
%%
fig.cluster_dpvs = make_fig('cluster_dilated_pvs');
fig.cluster_dpvs.showimg(pax_cluster.cluster_bvcsf_constdil(:,:,2,2))
%%
fig.cluster_dpvs.save2svg(directories.save_dir)
%%
fig.roi_pax.showrois(roilist,3,["pax","dipin","dipout","constricted_bv","dilated_bv"],["-c","-y","-y","y","y"])
