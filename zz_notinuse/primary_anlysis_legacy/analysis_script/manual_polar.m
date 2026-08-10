

manual_polar_roilist = roi_handle(fullfile(directories.save_dir,"manual_polar.mat"));
% load pax_cluster and roilist
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
polar_analysis.bvclusterstack = paxc.clust_med_bv(:,:,paxc.filtered_clusteridx);
polar_analysis.pvsclusterstack = paxc.clust_med_csf(:,:,paxc.filtered_clusteridx);
% For rgb comparison make rgb format stack
polar_analysis.bv_pvs_cluster = cat(4,mat2gray(polar_analysis.bvclusterstack),mat2gray(polar_analysis.pvsclusterstack),zeros(size(polar_analysis.bvclusterstack)));
%% Just extract constricted and dilated portion
polar_analysis.bv_pvs_clusterrgb_constricted = squeeze(polar_analysis.bv_pvs_cluster(:,:,polar_analysis.filtered_constricted_cluster_idx,:));
polar_analysis.bv_pvs_clusterrgb_dilated = squeeze(polar_analysis.bv_pvs_cluster(:,:,polar_analysis.filtered_dilated_cluster_idx,:));

%%
savestruct = struct();
%%
savestruct.bv_constricted = polar_analysis.bv_pvs_clusterrgb_constricted(:,:,1);
savestruct.pvs_constricted = polar_analysis.bv_pvs_clusterrgb_constricted(:,:,2);
savestruct.bv_dilated = polar_analysis.bv_pvs_clusterrgb_dilated(:,:,1);
savestruct.pvs_dilated = polar_analysis.bv_pvs_clusterrgb_dilated(:,:,2);
%%
save('250722_hql080savestruct.mat',"savestruct")

%%
figure()
imshow(savestruct.bv_constricted)


%% Manual constricted vessel boundary selection
try
manual_polar_roilist.addrgb(polar_analysis.bv_pvs_clusterrgb_constricted,'constricted_bv','polygon');
catch ME
    disp(ME.message)
    manual_polar_roilist.modifyrgb(polar_analysis.bv_pvs_clusterrgb_constricted,'constricted_bv');
end
%% Manual constricted PVS boundary selection
try
manual_polar_roilist.addrgb(polar_analysis.bv_pvs_clusterrgb_constricted,'constricted_pvs','polygon');
catch ME
    disp(ME.message)
    manual_polar_roilist.modifyrgb(polar_analysis.bv_pvs_clusterrgb_constricted,'constricted_pvs');
end
%% Manual dilated vessel boundary selection
try
manual_polar_roilist.addrgb(polar_analysis.bv_pvs_clusterrgb_dilated,'dilated_bv','polygon');
catch ME
    disp(ME.message)
    manual_polar_roilist.modifyrgb(polar_analysis.bv_pvs_clusterrgb_dilated,'dilated_bv');
end
%% Manual dilated pvs boundary selection
try
manual_polar_roilist.addrgb(polar_analysis.bv_pvs_clusterrgb_dilated,'dilated_pvs','polygon');
catch ME
    disp(ME.message)
    manual_polar_roilist.modifyrgb(polar_analysis.bv_pvs_clusterrgb_dilated,'dilated_pvs');
end
%% Save roilist
manual_polar_roilist.save2disk

%% Put
pax_vertices = roilist.getvertices('pax');
pax_vertices = pax_vertices(1:2,:);
% Primary axis used for 1d calculation
pax_angle = pax_vertices(2,:)-pax_vertices(1,:);
polar_analysis.parm.pax_angle = atan2d(pax_angle(2),pax_angle(1));
% calculate centerposition by projecting center point to PAX line
pax_vertices = roilist.getvertices('pax'); %(1:2,:);
pax_vertices = pax_vertices(1:2,:);
polar_analysis.parm.constricted_center = (min(roilist.getvertices('constricted_bv'),[],1) + max(roilist.getvertices('constricted_bv'),[],1))/2;
polar_analysis.parm.pax_center = analyze_dpoint2line(polar_analysis.parm.constricted_center,pax_vertices);

% Polar calculation for sinogram-kymograph
polar_analysis.parm.n_angles = 24;
[~,radius_map,angle_map] = analyze_polar(manual_polar_roilist.ROIs(1).Mask,...
    polar_analysis.parm.constricted_center,...
    polar_analysis.parm.pax_angle, polar_analysis.parm.n_angles ,true);

% figure()
% imagesc(radius_map)
% axis image
% hold on
% plot([cbv_vertex(:,1);cbv_vertex(1,1)],[cbv_vertex(:,2);cbv_vertex(1,2)],'r')
% Convert manualy drawn vertices into polar coordination
manual_polarstruct.dbv_thetaradi = roivtx2polar(manual_polar_roilist.getvertices('dilated_bv'),angle_map,radius_map);
manual_polarstruct.cbv_thetaradi = roivtx2polar(manual_polar_roilist.getvertices('constricted_bv'),angle_map,radius_map);
manual_polarstruct.dpvs_thetaradi = roivtx2polar(manual_polar_roilist.getvertices('dilated_pvs'),angle_map,radius_map);
manual_polarstruct.cpvs_thetaradi = roivtx2polar(manual_polar_roilist.getvertices('constricted_pvs'),angle_map,radius_map);
% Interpolate the manual points to calculate  thickness and etc..
manual_polarstruct.interp_dbv_thetaradi = sinogram_interpolation(manual_polarstruct.dbv_thetaradi);
manual_polarstruct.interp_cbv_thetaradi = sinogram_interpolation(manual_polarstruct.cbv_thetaradi);
manual_polarstruct.interp_dpvs_thetaradi = sinogram_interpolation(manual_polarstruct.dpvs_thetaradi);
manual_polarstruct.interp_cpvs_thetaradi = sinogram_interpolation(manual_polarstruct.cpvs_thetaradi);
save(fullfile(directories.save_dir,'manual_polarstruct.mat'),"manual_polarstruct")


%% Draw all interpolated boundary (constrict, dilate bv and pvs outter boundary)
fig.polar_boundary = make_fig('polar_boundary_all','polar');
%%
fig.polar_boundary.bring_fig
fig.polar_boundary.reset_axis
fig.polar_boundary.update_figsize([6,6])
fig.polar_boundary.plot_polar(deg2rad(manual_polarstruct.interp_dbv_thetaradi(1,:)),manual_polarstruct.interp_dbv_thetaradi(end,:),clee.clist.magenta,'none')

hold(fig.polar_boundary.ax,"on")
fig.polar_boundary.plot_polar(deg2rad(manual_polarstruct.interp_cbv_thetaradi(1,:)),manual_polarstruct.interp_cbv_thetaradi(end,:),clee.clist.red,'none')
fig.polar_boundary.plot_polar(deg2rad(manual_polarstruct.interp_dpvs_thetaradi(1,:)),manual_polarstruct.interp_dpvs_thetaradi(end,:),clee.clist.green,'none')
fig.polar_boundary.plot_polar(deg2rad(manual_polarstruct.interp_cpvs_thetaradi(1,:)),manual_polarstruct.interp_cpvs_thetaradi(end,:),clee.clist.darkgreen,'none')
fig.polar_boundary.save2svg(directories.save_dir);


%% %% make filtered (>10frames) cluster idx for polar analysis
pax_vertices = roilist.getvertices('pax');
pax_vertices = pax_vertices(1:2,:);
% Primary axis used for 1d calculation
pax_angle = pax_vertices(2,:)-pax_vertices(1,:);
polar_analysis.parm.pax_angle = atan2d(pax_angle(2),pax_angle(1));
% calculate centerposition by projecting center point to PAX line

%%
polar_analysis.constdial_sgram(:,:,1,1) = squeeze(polar_analysis.norm_sgram_kgph_bv(:,polar_analysis.filtered_constricted_cluster_idx,:)); % bv - const
polar_analysis.constdial_sgram(:,:,1,2) = squeeze(polar_analysis.norm_sgram_kgph_bv(:,polar_analysis.filtered_dilated_cluster_idx,:));% bv - dial
polar_analysis.constdial_sgram(:,:,2,1) = squeeze(polar_analysis.norm_sgram_kgph_pvs(:,polar_analysis.filtered_constricted_cluster_idx,:)); % pvs - const
polar_analysis.constdial_sgram(:,:,2,2) = squeeze(polar_analysis.norm_sgram_kgph_pvs(:,polar_analysis.filtered_dilated_cluster_idx,:)); % pvs - dial


function theta_radius = roivtx2polar(vertex,angle_map,radius_map)
    imgsize = size(angle_map);
    angle = angle_map(sub2ind(imgsize,vertex(:,2),vertex(:,1)));
    radius = radius_map(sub2ind(imgsize,vertex(:,2),vertex(:,1)));
    [deg_angle, order] = sort(angle,'ascend');
    radius = radius(order);
    angle = deg2rad(deg_angle)';
    radius = radius';
    theta_radius = [angle ; deg_angle' ;radius];
end

