% MAKE CLUSTER
% 1. Make cluster from pax kymograph based on pairwise correlation distance
% 2. Using cluster idx, make median filtered representation of cluster
% 3. Manually find clean representative constricted and dilated cluster

%% 1. Make cluster
%%
sleep.awakeidx = (MergedData.sleep.logicals.Manual.awakeLogical);
sleep.nremidx = (MergedData.sleep.logicals.Manual.nremLogical);
%%

nrem_kymograph = pax_fwhm.kymograph.kgph_lumen(:,logical(sleep.nremidx));
nrem_intracellular = roianalysis.preprocessed_ch1(:,:,logical(sleep.nremidx));
nrem_bv = roianalysis.preprocessed_ch2(:,:,logical(sleep.nremidx));
nrem_upperidx = pax_fwhm.idx.upperboundary(logical(sleep.nremidx));
nrem_loweridx =  pax_fwhm.idx.lowerboundary(logical(sleep.nremidx));
%%
awake_kymograph = pax_fwhm.kymograph.kgph_lumen(:,logical(sleep.awakeidx));
awake_intracellular = roianalysis.preprocessed_ch1(:,:,logical(sleep.awakeidx));
awake_bv = roianalysis.preprocessed_ch2(:,:,logical(sleep.awakeidx));
awake_upperidx = pax_fwhm.idx.upperboundary(logical(sleep.awakeidx));
awake_loweridx =  pax_fwhm.idx.lowerboundary(logical(sleep.awakeidx));

%%
nrem_cluster = pax_clustering(nrem_kymograph,nrem_bv,nrem_intracellular,nrem_upperidx,nrem_loweridx,25);

%%
awake_cluster = pax_clustering(awake_kymograph,awake_bv,awake_intracellular,awake_upperidx,awake_loweridx,25);

%%
save(fullfile(directories.save_dir,'nrem_cluster.mat'),"nrem_cluster")
%% Check sorting result
fig.paxsortedbv_fig = make_fig('paxsortedBV_figure_awake');
fig.paxsortedbv_fig.update_figsize([8 3])
fig.paxsortedbv_fig.reset_axis()
fig.paxsortedbv_fig.resolution = primary_datastruct.img_param.pixel2um;
fig.paxsortedbv_fig.plot_kymograph(awake_kymograph(:,awake_cluster.kgph_lumen_columnidx))
%%
fig.paxsortedbv_fig.plot_line(awake_upperidx(awake_cluster.kgph_lumen_columnidx),'r');
fig.paxsortedbv_fig.plot_line(awake_loweridx(awake_cluster.kgph_lumen_columnidx),'r');
%fig.paxsortedbv_fig.plot_line(pax_fwhm.idx.pvs_upperboundary(pax_cluster.kgph_lumen_columnidx),'g');
%fig.paxsortedbv_fig.plot_line(pax_fwhm.idx.pvs_lowerboundary(pax_cluster.kgph_lumen_columnidx),'g');
fig.paxsortedbv_fig.plot_xline(awake_cluster.clusterboundary(:,2),'m')
fig.paxsortedbv_fig.put_yaxistitle('Length (\mum)')
fig.paxsortedbv_fig.put_xaxistitle('Time (sec)')


%% Check sorting result
fig.paxsortedbv_fig = make_fig('paxsortedBV_figure_nrem');
fig.paxsortedbv_fig.update_figsize([8 3])
fig.paxsortedbv_fig.reset_axis()
fig.paxsortedbv_fig.resolution = primary_datastruct.img_param.pixel2um;
fig.paxsortedbv_fig.plot_kymograph(nrem_kymograph(:,nrem_cluster.kgph_lumen_columnidx))
%%
fig.paxsortedbv_fig.plot_line(nrem_upperidx(nrem_cluster.kgph_lumen_columnidx),'r');
fig.paxsortedbv_fig.plot_line(nrem_loweridx(nrem_cluster.kgph_lumen_columnidx),'r');
%fig.paxsortedbv_fig.plot_line(pax_fwhm.idx.pvs_upperboundary(pax_cluster.kgph_lumen_columnidx),'g');
%fig.paxsortedbv_fig.plot_line(pax_fwhm.idx.pvs_lowerboundary(pax_cluster.kgph_lumen_columnidx),'g');
fig.paxsortedbv_fig.plot_xline(nrem_cluster.clusterboundary(:,2),'m')
fig.paxsortedbv_fig.put_yaxistitle('Length (\mum)')
fig.paxsortedbv_fig.put_xaxistitle('Time (sec)')
%%
fig.paxsortedbv_fig.save2svg(directories.save_dir);
%%
fig.cluster_cbv = make_fig('cluster_constricted_bv');
fig.cluster_cbv.showimg(nrem_cluster.cluster_bvcsf_constdil(:,:,1,1))
%%
fig.cluster_cbv.save2svg(directories.save_dir)
%%
fig.cluster_dbv = make_fig('cluster_dilated_bv');
fig.cluster_dbv.showimg(nrem_cluster.cluster_bvcsf_constdil(:,:,1,2))
%%
fig.cluster_dbv.save2svg(directories.save_dir)
%%
fig.cluster_cpvs = make_fig('cluster_constricted_pvs_nrem');
fig.cluster_cpvs.showimg(nrem_cluster.cluster_bvcsf_constdil(:,:,2,1))
%%
fig.cluster_cpvs = make_fig('cluster_constricted_pvs_awake');
fig.cluster_cpvs.showimg(awake_cluster.cluster_bvcsf_constdil(:,:,2,1))
%%
fig.cluster_cpvs.save2svg(directories.save_dir)
%%
fig.cluster_dpvs = make_fig('cluster_dilated_pvs_nrem');
fig.cluster_dpvs.showimg(nrem_cluster.cluster_bvcsf_constdil(:,:,2,2))
%%
fig.cluster_dpvs.save2svg(directories.save_dir)
%%
fig.cluster_dpvs = make_fig('cluster_dilated_pvs_awake');
fig.cluster_dpvs.showimg(awake_cluster.cluster_bvcsf_constdil(:,:,2,2))
%%

function pax_cluster = pax_clustering(lumen_kymograph,bv_ch,pvs_ch,upperboundary,lowerboundary,cluster_num)
    pax_cluster = struct();
    pax_cluster.cluster_num = cluster_num;
    [pax_cluster.kgph_lumen_columnidx, pax_cluster.clusterboundary] = analysis_cluster_kymograph(lumen_kymograph,pax_cluster.cluster_num);
    % 2. Make median filter of images at each cluster
    pax_cluster.clust_med_bv = zeros([size(bv_ch,1,2),pax_cluster.cluster_num]);
    pax_cluster.clust_med_csf = zeros([size(pvs_ch,1,2),pax_cluster.cluster_num]);
    tmp.rearanged_bvstack = bv_ch(:,:,pax_cluster.kgph_lumen_columnidx);
    tmp.rearanged_csfstack = pvs_ch(:,:,pax_cluster.kgph_lumen_columnidx);
    for cluster_idx = 1:pax_cluster.cluster_num
        tmp.clusterstart = pax_cluster.clusterboundary(cluster_idx,1);
        tmp.clusterend = pax_cluster.clusterboundary(cluster_idx,2);
        pax_cluster.clust_med_bv(:,:,cluster_idx) = median(tmp.rearanged_bvstack(:,:,tmp.clusterstart:tmp.clusterend),3);
        pax_cluster.clust_med_csf(:,:,cluster_idx) = median(tmp.rearanged_csfstack(:,:,tmp.clusterstart:tmp.clusterend),3);
    end
    % 3. Make summary
    tmp.bv_diameter = lowerboundary-upperboundary;
    tmp.lowidxchange = lowerboundary - median(lowerboundary);
    tmp.upidxchange = median(upperboundary) -  upperboundary;
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
    pax_cluster.constricted_cluster_idx = util_checkstack(pax_cluster.clust_med_bv);
    pax_cluster.dilated_cluster_idx = util_checkstack(pax_cluster.clust_med_bv);
    %% (x,y,bv/csf,constricted/dilated)
    pax_cluster.cluster_bvcsf_constdil = zeros([size(tmp.rearanged_bvstack,1,2),2,2]); % (x,y,bv/csf,constricted/dilated)
    pax_cluster.cluster_bvcsf_constdil(:,:,1,1) = pax_cluster.clust_med_bv(:,:,pax_cluster.constricted_cluster_idx);
    pax_cluster.cluster_bvcsf_constdil(:,:,2,1) = pax_cluster.clust_med_csf(:,:,pax_cluster.constricted_cluster_idx);
    pax_cluster.cluster_bvcsf_constdil(:,:,1,2) = pax_cluster.clust_med_bv(:,:,pax_cluster.dilated_cluster_idx);
    pax_cluster.cluster_bvcsf_constdil(:,:,2,2) = pax_cluster.clust_med_csf(:,:,pax_cluster.dilated_cluster_idx);
    
    
    
    %% Filter length below 10 frames 5fps, 2 sec average (opt... may be I need range of vascular dilation for 
    % PIV --> in case porosity change only occur during initial dilation)
    pax_cluster.filtered_clusteridx = pax_cluster.cluster_summary(:,1);
    pax_cluster.filtered_clusteridx = pax_cluster.filtered_clusteridx(pax_cluster.cluster_summary(:,4)>10);
end
