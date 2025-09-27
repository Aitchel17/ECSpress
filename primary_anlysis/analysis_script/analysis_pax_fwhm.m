% Pre requisite: Run main_primary to determine 'pax' region of interest

%% PAX BV analysis
pax_fwhm = line_fwhm(roilist.getvertices('pax'));
pax_fwhm.addkymograph("lumen", roianalysis.preprocessed_ch1)
pax_fwhm.kymograph_afterprocess('lumen',[3 5])
pax_fwhm.fwhm("lumen",0.3,20);
%% PAX BV summary fig
paxbv_fig = make_fig('paxBV_figure');
paxbv_fig.update_figsize([8 3])
%%
paxbv_fig.reset_axis()
paxbv_fig.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed)
%%
paxbv_fig.plot_line(pax_fwhm.idx.upperboundary,'r');
paxbv_fig.plot_line(pax_fwhm.idx.lowerboundary,'r');
paxbv_fig.plot_line(pax_fwhm.idx.max_idx,'c')
paxbv_fig.plot_line(pax_fwhm.idx.downoffsetloc,'y')
paxbv_fig.plot_line(pax_fwhm.idx.upoffsetloc,'y')

paxbv_fig.plot_line(medfilt1(pax_fwhm.idx.upperboundary,5),'g');
paxbv_fig.plot_line(medfilt1(pax_fwhm.idx.lowerboundary,5),'g');
%%
bv_diameter = pax_fwhm.idx.lowerboundary-pax_fwhm.idx.upperboundary;

%%
lowidxchange = pax_fwhm.idx.lowerboundary - median(pax_fwhm.idx.lowerboundary);
upidxchange = median(pax_fwhm.idx.upperboundary)-  pax_fwhm.idx.upperboundary;

%% MAKE CLUSTER
[kgph_lumen_columnidx, clusterboundary] = analysis_cluster_kymograph(pax_fwhm.kymograph.kgph_lumen_processed,25);
%% Check sorting result
paxsortedbv_fig = make_fig('paxsortedBV_figure');
paxsortedbv_fig.update_figsize([8 3])
paxsortedbv_fig.reset_axis()
paxsortedbv_fig.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed(:,kgph_lumen_columnidx))
paxsortedbv_fig.plot_line(pax_fwhm.idx.upperboundary(kgph_lumen_columnidx),'r');
paxsortedbv_fig.plot_line(pax_fwhm.idx.lowerboundary(kgph_lumen_columnidx),'r');
%%
clustersort_bv_diameter = bv_diameter(kgph_lumen_columnidx);
clustersort_lowidxchange = lowidxchange(kgph_lumen_columnidx);
clustersort_upidxchange = upidxchange(kgph_lumen_columnidx);
%%
cluster_diameter_meanstd = zeros([size(clusterboundary,1),5]);
cluster_mean_ch1 = zeros([size,,size(clusterboundary,1)]);
cluster_mean_ch2 = zeros([size,,size(clusterboundary,1)]);
for cluster_idx = 1: size(clusterboundary,1)
    cluster_diameter_meanstd(cluster_idx,1) = mean(clustersort_bv_diameter(clusterboundary(cluster_idx,1):clusterboundary(cluster_idx,2)));
    cluster_diameter_meanstd(cluster_idx,2) = std(clustersort_bv_diameter(clusterboundary(cluster_idx,1):clusterboundary(cluster_idx,2)));
    cluster_diameter_meanstd(cluster_idx,3) = clusterboundary(cluster_idx,2)-clusterboundary(cluster_idx,1)+1;
    cluster_diameter_meanstd(cluster_idx,4) = mean(clustersort_lowidxchange(clusterboundary(cluster_idx,1):clusterboundary(cluster_idx,2)));
    cluster_diameter_meanstd(cluster_idx,5) = mean(clustersort_upidxchange(clusterboundary(cluster_idx,1):clusterboundary(cluster_idx,2)));
end
%%
cluster_idx = 3;
cluster_columnidx = kgph_lumen_columnidx(clusterboundary(cluster_idx,1):clusterboundary(cluster_idx,2));
cluster3_bv = roianalysis.preprocessed_ch1(:,:,cluster_columnidx);
%%
cluster_idx = 67;
cluster_columnidx = kgph_lumen_columnidx(clusterboundary(cluster_idx,1):clusterboundary(cluster_idx,2));
cluster67_bv = roianalysis.preprocessed_ch1(:,:,cluster_columnidx);
%%
util_checkstack(cluster3_bv)

%% PAX PVS analysis
pax_fwhm.addkymograph("pvs", roianalysis.preprocessed_ch2)
pax_fwhm.kymograph_afterprocess('lumen',[5 2])
pax_fwhm.addkymograph("pvs", roianalysis.preprocessed_ch2)
pax_fwhm.kymograph_afterprocess('pvs',[3 5])

%%
paxpvs_fig = make_fig('paxPVS_figure');
paxpvs_fig.update_figsize([8 3])
paxpvs_fig.reset_axis()
paxpvs_fig.plot_kymograph(pax_fwhm.kymograph.kgph_pvs_processed)
paxpvs_fig.plot_line(pax_fwhm.idx.upperboundary,'r');
paxpvs_fig.plot_line(pax_fwhm.idx.lowerboundary,'r');
%%
paxsortedbv_fig = make_fig('paxsortedBV_figure');
paxsortedbv_fig.update_figsize([8 3])
paxsortedbv_fig.reset_axis()
paxsortedbv_fig.plot_kymograph(pax_fwhm.kymograph.kgph_pvs_processed(:,kgph_lumen_columnidx))
paxsortedbv_fig.plot_line(pax_fwhm.idx.upperboundary(kgph_lumen_columnidx),'r');
paxsortedbv_fig.plot_line(pax_fwhm.idx.lowerboundary(kgph_lumen_columnidx),'r');
