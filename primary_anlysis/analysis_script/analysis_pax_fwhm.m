% Pre requisite: Run main_primary to determine 'pax' region of interest

%% PAX BV analysis
pax_fwhm = line_fwhm(roilist.getvertices('pax'));
pax_fwhm.t_axis = roianalysis.taxis;
pax_fwhm.addkymograph("lumen", roianalysis.preprocessed_ch1)
pax_fwhm.kymograph_afterprocess('lumen',[3 5])
pax_fwhm.fwhm("lumen");

%% PAX PVS analysis
pax_fwhm.addkymograph("pvs", roianalysis.preprocessed_ch2)
pax_fwhm.kymograph_afterprocess('pvs',[3 5])
%%
pax_fwhm.pvsanalysis(0.5,5)
%%
pax_fwhm.save2disk(directories.save_dir);
%% PAX BV summary fig
fig.paxbv_fig = make_fig('paxBV_figure');
fig.paxbv_fig.update_figsize([8 3])
fig.paxbv_fig.resolution = primary_datastruct.img_param.pixel2um;
%%
fig.paxbv_fig.bring_fig
fig.paxbv_fig.reset_axis()
fig.paxbv_fig.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed,roianalysis.taxis);
%%
fig.paxbv_fig.plot_line(pax_fwhm.idx.upperboundary,'r');
fig.paxbv_fig.plot_line(pax_fwhm.idx.lowerboundary,'r');
%%
fig.paxbv_fig.bring_fig
fig.paxbv_fig.plot_line(pax_fwhm.idx.pvs_upperboundary,'g')
fig.paxbv_fig.plot_line(pax_fwhm.idx.pvs_lowerboundary,'g')
%%
fig.paxbv_fig.put_yaxistitle('Length (\mum)')
fig.paxbv_fig.put_xaxistitle('Time (sec)')
fig.paxbv_fig.save2svg(directories.save_dir);

%%
fig.paxbv_fig.bring_fig
%% PAX PVS summary fig
fig.paxpvs_fig = make_fig('paxPVS_figure');
fig.paxpvs_fig.resolution = primary_datastruct.img_param.pixel2um;
fig.paxpvs_fig.update_figsize([8 3])
%%
fig.paxpvs_fig.reset_axis()
fig.paxpvs_fig.plot_kymograph(pax_fwhm.kymograph.kgph_pvs_processed,roianalysis.taxis)
%%
fig.paxpvs_fig.plot_line(pax_fwhm.idx.upperboundary,'r');
fig.paxpvs_fig.plot_line(pax_fwhm.idx.lowerboundary,'r');
%%
fig.paxpvs_fig.plot_line(pax_fwhm.idx.pvs_upperboundary,'g')
fig.paxpvs_fig.plot_line(pax_fwhm.idx.pvs_lowerboundary,'g')
fig.paxpvs_fig.put_yaxistitle('Length (\mum)')
fig.paxpvs_fig.put_xaxistitle('Time (sec)')
fig.paxpvs_fig.save2svg(directories.save_dir);

%%
fig.paxpvs_fig.bring_fig
