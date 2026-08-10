%% load stack, vertices from main_primary script

%%
polar_analysis = struct();

%% load pax_cluster and roilist
if exist('pax_cluster','var')
    paxc = pax_cluster;
elseif isfield(primary_datastruct,'pax_cluster')
    paxc = primary_datastruct.pax_cluster;
else
    disp('Prerequisite run cluster_analysis_script.mat')
end
if isfield(primary_datastruct,'roilist')
    roilist = primary_datastruct.roilist;
end
%% make filtered (>10frames) cluster idx and clustered stack for polar analysis
polar_analysis.filtered_dilated_cluster_idx = find(paxc.filtered_clusteridx==paxc.dilated_cluster_idx);
polar_analysis.filtered_constricted_cluster_idx = find(paxc.filtered_clusteridx==paxc.constricted_cluster_idx);
tmp.bvclusterstack = paxc.clust_med_bv(:,:,paxc.filtered_clusteridx);
tmp.pvsclusterstack = paxc.clust_med_csf(:,:,paxc.filtered_clusteridx);
%% %% make filtered (>10frames) cluster idx for polar analysis
tmp.pax_vertices = roilist.getvertices('pax');
tmp.pax_vertices = tmp.pax_vertices(1:2,:);
% Primary axis used for 1d calculation
tmp.pax_angle = tmp.pax_vertices(2,:)-tmp.pax_vertices(1,:);
polar_analysis.parm.pax_angle = atan2d(tmp.pax_angle(2),tmp.pax_angle(1));
% calculate centerposition by projecting center point to PAX line
tmp.pax_vertices = roilist.getvertices('pax'); %(1:2,:);
tmp.pax_vertices = tmp.pax_vertices(1:2,:);
polar_analysis.parm.constricted_center = (min(roilist.getvertices('constricted_bv'),[],1) + max(roilist.getvertices('constricted_bv'),[],1))/2;
polar_analysis.parm.pax_center = analyze_dpoint2line(polar_analysis.parm.constricted_center,tmp.pax_vertices);

%% Polar calculation for sinogram-kymograph
polar_analysis.parm.n_angles = 24;
polar_analysis.bvcluster = analyze_polar(tmp.bvclusterstack,...
    polar_analysis.parm.constricted_center,...
    polar_analysis.parm.pax_angle, polar_analysis.parm.n_angles ,true);
polar_analysis.pvscluster = analyze_polar(tmp.pvsclusterstack,...
    polar_analysis.parm.constricted_center,...
    polar_analysis.parm.pax_angle, polar_analysis.parm.n_angles ,true);

% Trim the kymograph at each angle and merge to make sinogram-kymograph array
polar_analysis.parm.minlength = min(arrayfun(@(s) size(s.kymograph, 1), polar_analysis.bvcluster));
% BV trimming
tmp.trimmed_polarbv = arrayfun(@(s) s.kymograph(1:polar_analysis.parm.minlength,:), ...
    polar_analysis.bvcluster,'UniformOutput', false);
polar_analysis.sgram_kgph_bv = cat(3,tmp.trimmed_polarbv{:}); % [rho,cluster,angle]
% PVS trimming
tmp.trimmed_polarpvs = arrayfun(@(s) s.kymograph(1:polar_analysis.parm.minlength,:), ...
    polar_analysis.pvscluster,'UniformOutput', false);
polar_analysis.sgram_kgph_pvs = cat(3,tmp.trimmed_polarpvs{:});
% Normalization at each column
% bv sinogram-kymograph normalization
polar_analysis.norm_sgram_kgph_bv = polar_analysis.sgram_kgph_bv-min(polar_analysis.sgram_kgph_bv,[],1);
polar_analysis.norm_sgram_kgph_bv = polar_analysis.norm_sgram_kgph_bv./max(polar_analysis.norm_sgram_kgph_bv,[],1);
% pvs sinogram-kymograph normalization
polar_analysis.norm_sgram_kgph_pvs = polar_analysis.sgram_kgph_pvs-min(polar_analysis.sgram_kgph_pvs,[],1);
polar_analysis.norm_sgram_kgph_pvs = polar_analysis.norm_sgram_kgph_pvs./max(polar_analysis.norm_sgram_kgph_pvs,[],1);
% Make input array for bv/pvs constricted/dilated sinogram
polar_analysis.constdial_sgram = zeros([size(polar_analysis.norm_sgram_kgph_bv,1,3),2,2]);
%
polar_analysis.constdial_sgram(:,:,1,1) = squeeze(polar_analysis.norm_sgram_kgph_bv(:,polar_analysis.filtered_constricted_cluster_idx,:)); % bv - const
polar_analysis.constdial_sgram(:,:,1,2) = squeeze(polar_analysis.norm_sgram_kgph_bv(:,polar_analysis.filtered_dilated_cluster_idx,:));% bv - dial
polar_analysis.constdial_sgram(:,:,2,1) = squeeze(polar_analysis.norm_sgram_kgph_pvs(:,polar_analysis.filtered_constricted_cluster_idx,:)); % pvs - const
polar_analysis.constdial_sgram(:,:,2,2) = squeeze(polar_analysis.norm_sgram_kgph_pvs(:,polar_analysis.filtered_dilated_cluster_idx,:)); % pvs - dial





%%
polar_analysis.parm.bvthr = 0.25;
polar_analysis.parm.bvoffset = 5;
polar_analysis.parm.pvsthr = 0.5;
polar_analysis.parm.pvsoffset = 5;

polar_analysis.polar_boundary = polar_halfmax(polar_analysis.constdial_sgram,...
    polar_analysis.parm.bvthr,polar_analysis.parm.bvoffset,...
    polar_analysis.parm.pvsthr,polar_analysis.parm.pvsoffset);
polar_analysis.bvdilation =polar_analysis.polar_boundary.dilate_bvboundary-polar_analysis.polar_boundary.const_bvboundary;
polar_analysis.compression = polar_analysis.polar_boundary.dilate_pvsboundary-polar_analysis.polar_boundary.const_pvsboundary;
polar_analysis.dilate_pvsthickness = polar_analysis.polar_boundary.dilate_pvsboundary-polar_analysis.polar_boundary.dilate_bvboundary;
polar_analysis.const_pvsthickness = polar_analysis.polar_boundary.const_pvsboundary-polar_analysis.polar_boundary.const_bvboundary;
polar_analysis.pvschange = polar_analysis.const_pvsthickness- polar_analysis.dilate_pvsthickness;
%%
save(fullfile(directories.save_dir,'polar_analysis.mat'),"polar_analysis")

%%
%% confirm by figure
fig.polar_sgram_cbv = make_fig('polar_sinogram_bv_const');
%%
fig.polar_sgram_cbv.bring_fig
fig.polar_sgram_cbv.update_figsize([8,5])
fig.polar_sgram_cbv.reset_axis
fig.polar_sgram_cbv.plot_kymograph(polar_analysis.constdial_sgram(:,:,1,1)) % constricted pvs
fig.polar_sgram_cbv.plot_line(polar_analysis.polar_boundary.const_pvsboundary,'g');
fig.polar_sgram_cbv.plot_line(polar_analysis.polar_boundary.const_bvboundary,'r');
fig.polar_sgram_cbv.plot_line(polar_analysis.polar_boundary.dilate_pvsboundary,clee.clist.darkgreen);
fig.polar_sgram_cbv.plot_line(polar_analysis.polar_boundary.dilate_bvboundary,clee.clist.magenta);
fig.polar_sgram_cbv.save2svg(directories.save_dir)
%%
fig.polar_sgram_dbv = make_fig('polar_sinogram_bv_dilate');
%%
fig.polar_sgram_dbv.bring_fig
fig.polar_sgram_dbv.update_figsize([8,5])
fig.polar_sgram_dbv.reset_axis
fig.polar_sgram_dbv.plot_kymograph(polar_analysis.constdial_sgram(:,:,1,2)) % constricted pvs
fig.polar_sgram_dbv.plot_line(polar_analysis.polar_boundary.const_pvsboundary,'g');
fig.polar_sgram_dbv.plot_line(polar_analysis.polar_boundary.const_bvboundary,'r');
fig.polar_sgram_dbv.plot_line(polar_analysis.polar_boundary.dilate_pvsboundary,clee.clist.darkgreen);
fig.polar_sgram_dbv.plot_line(polar_analysis.polar_boundary.dilate_bvboundary,clee.clist.magenta);
fig.polar_sgram_dbv.save2svg(directories.save_dir)

%%
fig.polar_sgram_cpvs = make_fig('polar_sinogram_pvs_const');
fig.polar_sgram_cpvs.update_figsize([8,5])
%%
fig.polar_sgram_cpvs.bring_fig

fig.polar_sgram_cpvs.reset_axis
fig.polar_sgram_cpvs.plot_kymograph(polar_analysis.constdial_sgram(:,:,2,1)) % constricted pvs
fig.polar_sgram_cpvs.plot_line(polar_analysis.polar_boundary.const_pvsboundary,'g');
fig.polar_sgram_cpvs.plot_line(polar_analysis.polar_boundary.const_bvboundary,'r');
fig.polar_sgram_cpvs.plot_line(polar_analysis.polar_boundary.dilate_pvsboundary,clee.clist.darkgreen);
fig.polar_sgram_cpvs.plot_line(polar_analysis.polar_boundary.dilate_bvboundary,clee.clist.magenta);
fig.polar_sgram_cpvs.save2svg(directories.save_dir)

%%
fig.polar_sgram_dpvs = make_fig('polar_sinogram_pvs_dilate');
fig.polar_sgram_dpvs.update_figsize([8,5])
%%
fig.polar_sgram_dpvs.bring_fig
fig.polar_sgram_dpvs.reset_axis
fig.polar_sgram_dpvs.plot_kymograph(polar_analysis.constdial_sgram(:,:,2,2)) % constricted pvs
fig.polar_sgram_dpvs.plot_line(polar_analysis.polar_boundary.const_pvsboundary,'g');
fig.polar_sgram_dpvs.plot_line(polar_analysis.polar_boundary.const_bvboundary,'r');
fig.polar_sgram_dpvs.plot_line(polar_analysis.polar_boundary.dilate_pvsboundary,clee.clist.darkgreen);
fig.polar_sgram_dpvs.plot_line(polar_analysis.polar_boundary.dilate_bvboundary,clee.clist.magenta);
fig.polar_sgram_dpvs.save2svg(directories.save_dir)
fig.polar_sgram_dpvs.save2svg(directories.save_dir)



%% Draw all boundary (constrict, dilate bv and pvs outter boundary)
fig.polar_boundary = make_fig('polar_boundary_all','polar');
fig.polar_boundary.update_figsize([6,6])
fig.polar_boundary.plot_polar(polar_analysis.polar_boundary.const_bvboundary,clee.clist.red,'*')
hold(fig.polar_boundary.ax,"on")
fig.polar_boundary.plot_polar(polar_analysis.polar_boundary.dilate_bvboundary,clee.clist.magenta,'*')
fig.polar_boundary.plot_polar(polar_analysis.polar_boundary.const_pvsboundary,clee.clist.darkgreen,'*')
fig.polar_boundary.plot_polar(polar_analysis.polar_boundary.dilate_pvsboundary,clee.clist.green,'*')
fig.polar_boundary.save2svg(directories.save_dir)

%% Draw all thickness
fig.polar_thickness = make_fig('polar_thickness','polar');
%%
fig.polar_thickness.bring_fig
fig.polar_thickness.reset_axis
fig.polar_thickness.update_figsize([6,6])
%
fig.polar_thickness.plot_polar(polar_analysis.polar_boundary.const_bvboundary,clee.clist.red,'*')
hold(fig.polar_thickness.ax,"on")
fig.polar_thickness.plot_polar(polar_analysis.polar_boundary.dilate_bvboundary,clee.clist.magenta,'*')
fig.polar_thickness.plot_polar(polar_analysis.const_pvsthickness,clee.clist.darkgreen,'*')
fig.polar_thickness.plot_polar(polar_analysis.dilate_pvsthickness,clee.clist.green,'*')

%%
fig.polar_diffboundary = make_fig('polar_boundary_differential','polar');
%%
fig.polar_diffboundary.reset_axis
fig.polar_diffboundary.update_figsize([6,6])
%
fig.polar_diffboundary.plot_polar(polar_analysis.bvdilation,clee.clist.red,'*')
hold(fig.polar_diffboundary.ax,"on")
fig.polar_diffboundary.plot_polar(polar_analysis.pvschange,clee.clist.darkgreen,'*')
%%
fig.polar_diffboundary.plot_polar(polar_analysis.compression,clee.clist.black,'*')
fig.polar_diffboundary.save2svg(directories.save_dir)

%% show center position
fig.roi_polarpaxcenter = make_fig('polar_paxcenter');
fig.roi_polarpaxcenter.update_figsize([8 6])
%%
fig.roi_polarpaxcenter.reset_axis
fig.roi_polarpaxcenter.showrois(roilist,2,["constricted_bv","dilated_bv"],["-y","-y"])
hold on
ax = gca;
plot(ax,[tmp.pax_vertices(:, 1); tmp.pax_vertices(1, 1)], [tmp.pax_vertices(:, 2); tmp.pax_vertices(1, 2)], 'r-', 'LineWidth', 2);
plot(ax, polar_analysis.parm.constricted_center(1), polar_analysis.parm.constricted_center(2), 'g+', 'MarkerSize', 10, 'LineWidth', 1);
plot(ax, polar_analysis.parm.pax_center(1), polar_analysis.parm.pax_center(2), 'g+', 'MarkerSize', 10, 'LineWidth', 1);

%%
fig.roi_polarpaxcenter.save2svg(directories.save_dir)