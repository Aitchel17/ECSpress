function analysis_pax_makefig(pax_fwhm, t_axis, pixel2um, save_dir)
% ANALYSIS_PAX_MAKEFIG Generates and saves figures for PAX analysis.
%   pax_fwhm: line_fwhm object containing analysis results.
%   t_axis: Time axis vector.
%   pixel2um: Conversion factor (microns per pixel).
%   save_dir: Directory where figures will be saved (SVGs).

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% Ensure pax_fwhm has t_axis for consistency
if isempty(pax_fwhm.t_axis) && ~isempty(t_axis)
    pax_fwhm.t_axis = t_axis;
elseif ~isempty(pax_fwhm.t_axis)
    t_axis = pax_fwhm.t_axis;
end
%%
clee = color_lee;

%% 1. PAX BV summary fig
fig_bv = make_fig('paxBV_figure');
fig_bv.update_figsize([8 3]);
fig_bv.resolution = pixel2um;
fig_bv.reset_axis();

fig_bv.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed, t_axis);

% Plot boundaries - use new names 'upperBVboundary' / 'lowerBVboundary'
% Also check for legacy names just in case
if isfield(pax_fwhm.idx, 'clean_lowerBVboundary')
    lower = pax_fwhm.idx.clean_lowerBVboundary;
    upper = pax_fwhm.idx.clean_upperBVboundary;
elseif isfield(pax_fwhm.idx, 'clean_lowerboundary')
    lower = pax_fwhm.idx.clean_lowerboundary;
    upper = pax_fwhm.idx.clean_upperboundary;
else
    warning('BV boundaries not found for plotting.');
    lower = []; upper = [];
end

if ~isempty(lower)
    fig_bv.plot_line(lower, 'r');
end
if ~isempty(upper)
    fig_bv.plot_line(upper, 'r');
end

fig_bv.put_yaxistitle('Length (\mum)');
fig_bv.put_xaxistitle('Time (sec)');
fig_bv.save2svg(save_dir);
disp(['Saved paxBV_figure to ', save_dir]);


%% 2. PVS innerboundary and outter boundary (paxPVSonly_figure)
fig_pvsonly = make_fig('paxPVSonly_figure');
fig_pvsonly.update_figsize([8 3]);
fig_pvsonly.resolution = pixel2um;
fig_pvsonly.reset_axis();

fig_pvsonly.plot_kymograph(pax_fwhm.kymograph.kgph_pvs_processed, t_axis);

% Plot lines if they exist
if isfield(pax_fwhm.idx, 'clean_pvsupcent_idx')
    fig_pvsonly.plot_line(pax_fwhm.idx.clean_pvsupcent_idx, 'r')
end
if isfield(pax_fwhm.idx, 'clean_pvsdowncent_idx')
    fig_pvsonly.plot_line(pax_fwhm.idx.clean_pvsdowncent_idx, 'r')
end
if isfield(pax_fwhm.idx, 'clean_pvsupedge_idx')
    fig_pvsonly.plot_line(pax_fwhm.idx.clean_pvsupedge_idx, 'r')
end
if isfield(pax_fwhm.idx, 'clean_pvsdownedge_idx')
    fig_pvsonly.plot_line(pax_fwhm.idx.clean_pvsdownedge_idx, 'r')
end

fig_pvsonly.put_yaxistitle('Length (\mum)');
fig_pvsonly.put_xaxistitle('Time (sec)');
fig_pvsonly.save2svg(save_dir);
disp(['Saved paxPVSonly_figure to ', save_dir]);

%% 3. PAX PVS summary fig (paxPVS_figure)
fig_pvs = make_fig('paxPVS_figure');
fig_pvs.resolution = pixel2um;
fig_pvs.update_figsize([8 3]);
fig_pvs.reset_axis();

% Prepare RGB Kymograph
rgb_kymo = plot_make_rgb(pax_fwhm.kymograph.kgph_lumen_processed, pax_fwhm.kymograph.kgph_pvs_processed);

% Plot RGB Kymograph
fig_pvs.plot_kymograph(rgb_kymo, t_axis);

% Plot boundaries (filtered)
% BV Boundaries: Red Background + Yellow Dotted Overlay (Wide)
if ~isempty(upper)
    fig_pvs.plot_line(upper, clee.clist.yellow, 'none', '-', 3);
    fig_pvs.plot_line(upper, clee.clist.red, 'none', '-');

end
if ~isempty(lower)
    fig_pvs.plot_line(lower, clee.clist.yellow, 'none', '-', 3);
    fig_pvs.plot_line(lower, clee.clist.red, 'none', '-');

end

% PVS Boundaries: DarkGreen Background + Yellow Dotted Overlay (Wide)
if isfield(pax_fwhm.idx, 'clean_pvsupedge_idx')
    fig_pvs.plot_line(pax_fwhm.idx.clean_pvsupedge_idx, clee.clist.yellow, 'none', '-', 3);
    fig_pvs.plot_line(pax_fwhm.idx.clean_pvsupedge_idx, clee.clist.darkgreen, 'none', '-');

end
if isfield(pax_fwhm.idx, 'clean_pvsdownedge_idx')
    fig_pvs.plot_line(pax_fwhm.idx.clean_pvsdownedge_idx, clee.clist.yellow, 'none', '-', 3);
    fig_pvs.plot_line(pax_fwhm.idx.clean_pvsdownedge_idx, clee.clist.darkgreen, 'none', '-');
end

fig_pvs.put_yaxistitle('Length (\mum)');
fig_pvs.put_xaxistitle('Time (sec)');
fig_pvs.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);

%% 4. Vessel and perivascular space (total) thickness
fig_pvsthickness = make_fig('pvsTotalthickness_figure');
%%
fig_pvsthickness.reset_axis
%%
fig_pvsthickness.update_figsize([8 3]);
fig_pvsthickness.resolution = pixel2um;
%%
fig_pvsthickness.loc.x = t_axis;
%%
fig_pvsthickness.plot_line(pax_fwhm.thickness.bv,clee.clist.magenta)
%%
hold(fig_pvsthickness.ax)
fig_pvsthickness.plot_line(pax_fwhm.thickness.totalpvs,clee.clist.darkgreen)
%%
maxyrange = max([pax_fwhm.thickness.totalpvs,pax_fwhm.thickness.bv])*pixel2um;
minyrange = min([pax_fwhm.thickness.totalpvs,pax_fwhm.thickness.bv])*pixel2um;
%%
ylim(fig_pvsthickness.ax,[minyrange,maxyrange])
%%
fig_pvsthickness.put_yaxistitle('Length (\mum)');
fig_pvsthickness.put_xaxistitle('Time (sec)');
fig_pvsthickness.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);



%% 5. Vessel diameter changes and PVS thicness changes
fig_totalpvs_changes = make_fig('pvsTotalthicknesschages_figure');

fig_totalpvs_changes.reset_axis

fig_totalpvs_changes.update_figsize([8 3]);
fig_totalpvs_changes.resolution = pixel2um;

fig_totalpvs_changes.loc.x = t_axis;
hold(fig_totalpvs_changes.ax)

fig_totalpvs_changes.plot_line(pax_fwhm.thickness.bvchanges,clee.clist.magenta);
fig_totalpvs_changes.plot_line(pax_fwhm.thickness.ecschanges_residual,clee.clist.coloredorange)
fig_totalpvs_changes.plot_line(pax_fwhm.thickness.pvschanges_total,clee.clist.darkgreen);

maxyrange = max([pax_fwhm.thickness.bvchanges,pax_fwhm.thickness.ecschanges_residual,pax_fwhm.thickness.pvschanges_total])*pixel2um;
minyrange = min([pax_fwhm.thickness.bvchanges,pax_fwhm.thickness.ecschanges_residual,pax_fwhm.thickness.pvschanges_total])*pixel2um;
ylim(fig_totalpvs_changes.ax,[minyrange,maxyrange])
fig_totalpvs_changes.put_yaxistitle('Length (\mum)');
fig_totalpvs_changes.put_xaxistitle('Time (sec)');
fig_totalpvs_changes.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);

%% 6. Dynamic and Static PVS thickness changes
fig_dspvs_changes = make_fig('pvsDynamicandStaticthicknesschages_figure');
hold(fig_dspvs_changes.ax)
fig_dspvs_changes.reset_axis
fig_dspvs_changes.update_figsize([8 3]);
fig_dspvs_changes.resolution = pixel2um;
fig_dspvs_changes.loc.x = t_axis;
fig_dspvs_changes.plot_line(pax_fwhm.thickness.bvchanges,clee.clist.magenta);
fig_dspvs_changes.plot_line(pax_fwhm.thickness.pvschanges_dynamic,clee.clist.darkgreen,'none','-');
fig_dspvs_changes.plot_line(pax_fwhm.thickness.pvschanges_static,clee.clist.lightgreen,'none','-');
maxyrange = max([pax_fwhm.thickness.bvchanges,pax_fwhm.thickness.pvschanges_dynamic,pax_fwhm.thickness.pvschanges_static])*pixel2um;
minyrange = min([pax_fwhm.thickness.bvchanges,pax_fwhm.thickness.pvschanges_dynamic,pax_fwhm.thickness.pvschanges_static])*pixel2um;
ylim(fig_dspvs_changes.ax,[minyrange,maxyrange])
fig_dspvs_changes.put_yaxistitle('Length (\mum)');
fig_dspvs_changes.put_xaxistitle('Time (sec)');
fig_dspvs_changes.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);

%% 7. Displacement of dynamic boundary
fig_dynamicdisplacement = make_fig('pvsDynamicDisplacement_figure');
fig_dynamicdisplacement.reset_axis
fig_dynamicdisplacement.update_figsize([8 3]);
fig_dynamicdisplacement.resolution = pixel2um;
fig_dynamicdisplacement.loc.x = t_axis;
hold(fig_dynamicdisplacement.ax)
fig_dynamicdisplacement.plot_line(pax_fwhm.displacement.dynamicpvs,clee.clist.orange);
fig_dynamicdisplacement.plot_line(pax_fwhm.displacement.dynamicbv,clee.clist.magenta);
fig_dynamicdisplacement.put_yaxistitle('Length (\mum)');
fig_dynamicdisplacement.put_xaxistitle('Time (sec)');
fig_dynamicdisplacement.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);

%% 8. Displacement of static boundary
fig_staticdisplacement = make_fig('pvsStaticDisplacement_figure');

fig_staticdisplacement.reset_axis
fig_staticdisplacement.update_figsize([8 3]);
fig_staticdisplacement.resolution = pixel2um;
fig_staticdisplacement.loc.x = t_axis;
hold(fig_staticdisplacement.ax)
fig_staticdisplacement.plot_line(pax_fwhm.displacement.staticpvs,clee.clist.darkgreen);
fig_staticdisplacement.plot_line(pax_fwhm.displacement.staticbv,clee.clist.magenta);
fig_staticdisplacement.put_yaxistitle('Length (\mum)');
fig_staticdisplacement.put_xaxistitle('Time (sec)');
fig_staticdisplacement.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);
%% 9. Comparison of static and dynamic pvs displacement
fig_pvsdisplacement = make_fig('pvsDisplacement_figure');

fig_pvsdisplacement.reset_axis
fig_pvsdisplacement.update_figsize([8 3]);
fig_pvsdisplacement.resolution = pixel2um;
fig_pvsdisplacement.loc.x = t_axis;
hold(fig_pvsdisplacement.ax)
fig_pvsdisplacement.plot_line(pax_fwhm.displacement.staticpvs,clee.clist.darkgreen);
fig_pvsdisplacement.plot_line(pax_fwhm.displacement.dynamicpvs,clee.clist.orange);
fig_pvsdisplacement.put_yaxistitle('Length (\mum)');
fig_pvsdisplacement.put_xaxistitle('Time (sec)');
fig_pvsdisplacement.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);
end
