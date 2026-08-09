function fig = analysis_pax_makefig(pax_fwhm, t_axis, pixel2um, save_dir)
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
fig.bv = make_fig('paxBV_figure');
fig.bv.update_figsize([8 3]);
fig.bv.resolution = pixel2um;
fig.bv.reset_axis();

fig.bv.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed, t_axis);

%%

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
    fig.bv.plot_line(lower, 'r');
end
if ~isempty(upper)
    fig.bv.plot_line(upper, 'r');
end
if isfield(pax_fwhm.idx, 'clean_pvsupedge_idx')
    fig.bv.plot_line(pax_fwhm.idx.clean_pvsupedge_idx, 'g')
end

if isfield(pax_fwhm.idx, 'clean_pvsdownedge_idx')
    fig.bv.plot_line(pax_fwhm.idx.clean_pvsdownedge_idx, 'g')
end


%%
fig.bv.put_yaxistitle('Length (\mum)');
fig.bv.put_xaxistitle('Time (sec)');
fig.bv.save2svg(save_dir);
disp(['Saved paxBV_figure to ', save_dir]);


%% 2. PVS innerboundary and outter boundary (paxPVSonly_figure)
fig.pvsonly = make_fig('paxPVSonly_figure');
fig.pvsonly.update_figsize([8 3]);
fig.pvsonly.resolution = pixel2um;
fig.pvsonly.reset_axis();

fig.pvsonly.plot_kymograph(pax_fwhm.kymograph.kgph_pvs_processed, t_axis);

% Plot lines if they exist
if isfield(pax_fwhm.idx, 'clean_pvsupcent_idx')
    fig.pvsonly.plot_line(pax_fwhm.idx.clean_pvsupcent_idx, 'r')
end

if isfield(pax_fwhm.idx, 'clean_pvsdowncent_idx')
    fig.pvsonly.plot_line(pax_fwhm.idx.clean_pvsdowncent_idx, 'r')
end

if isfield(pax_fwhm.idx, 'clean_pvsupedge_idx')
    fig.pvsonly.plot_line(pax_fwhm.idx.clean_pvsupedge_idx, 'r')
end

if isfield(pax_fwhm.idx, 'clean_pvsdownedge_idx')
    fig.pvsonly.plot_line(pax_fwhm.idx.clean_pvsdownedge_idx, 'r')
end

fig.pvsonly.put_yaxistitle('Length (\mum)');
fig.pvsonly.put_xaxistitle('Time (sec)');
fig.pvsonly.save2svg(save_dir);
disp(['Saved paxPVSonly_figure to ', save_dir]);

%% 3. PAX PVS summary fig (paxPVS_figure)
fig.pvs = make_fig('paxPVS_figure');
fig.pvs.resolution = pixel2um;
fig.pvs.update_figsize([8 3]);
fig.pvs.reset_axis();

% Prepare RGB Kymograph
rgb_kymo = plot_make_rgb(pax_fwhm.kymograph.kgph_lumen, pax_fwhm.kymograph.kgph_pvs);

% Plot RGB Kymograph
fig.pvs.plot_kymograph(rgb_kymo, t_axis);

% Plot boundaries (filtered)
% BV Boundaries: Red Background + Yellow Dotted Overlay (Wide)
if ~isempty(upper)
    fig.pvs.plot_line(upper, clee.clist.yellow, 'none', '-', 1.5);
    fig.pvs.plot_line(upper, clee.clist.red, 'none', '-',0.5);

end
if ~isempty(lower)
    fig.pvs.plot_line(lower, clee.clist.yellow, 'none', '-', 1.5);
    fig.pvs.plot_line(lower, clee.clist.red, 'none', '-',0.5);

end

% PVS Boundaries: DarkGreen Background + Yellow Dotted Overlay (Wide)
if isfield(pax_fwhm.idx, 'clean_pvsupedge_idx')
    fig.pvs.plot_line(pax_fwhm.idx.clean_pvsupedge_idx, clee.clist.yellow, 'none', '-', 1.5);
    fig.pvs.plot_line(pax_fwhm.idx.clean_pvsupedge_idx, clee.clist.darkgreen, 'none', '-',0.5);

end
if isfield(pax_fwhm.idx, 'clean_pvsdownedge_idx')
    fig.pvs.plot_line(pax_fwhm.idx.clean_pvsdownedge_idx, clee.clist.yellow, 'none', '-', 1.5);
    fig.pvs.plot_line(pax_fwhm.idx.clean_pvsdownedge_idx, clee.clist.darkgreen, 'none', '-',0.5);
end

fig.pvs.put_yaxistitle('Length (\mum)');
fig.pvs.put_xaxistitle('Time (sec)');
fig.pvs.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);

%% 4. Vessel and perivascular space (total) thickness
fig.pvsthickness = make_fig('pvsTotalthickness_figure');
%%
fig.pvsthickness.reset_axis
%%
fig.pvsthickness.update_figsize([8 3]);
fig.pvsthickness.resolution = pixel2um;
%%
fig.pvsthickness.loc.x = t_axis;
%%
fig.pvsthickness.plot_line(pax_fwhm.thickness.bv,clee.clist.magenta)
%%
hold(fig.pvsthickness.ax)
fig.pvsthickness.plot_line(pax_fwhm.thickness.totalpvs,clee.clist.darkgreen)
%%
maxyrange = max([pax_fwhm.thickness.totalpvs,pax_fwhm.thickness.bv])*pixel2um;
minyrange = min([pax_fwhm.thickness.totalpvs,pax_fwhm.thickness.bv])*pixel2um;
%%
ylim(fig.pvsthickness.ax,[minyrange,maxyrange])
%%
fig.pvsthickness.put_yaxistitle('Length (\mum)');
fig.pvsthickness.put_xaxistitle('Time (sec)');
fig.pvsthickness.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);



%% 5. Vessel diameter changes and PVS thicness changes
fig.totalpvs_changes = make_fig('pvsTotalthicknesschages_figure');

fig.totalpvs_changes.reset_axis

fig.totalpvs_changes.update_figsize([8 3]);
fig.totalpvs_changes.resolution = pixel2um;

fig.totalpvs_changes.loc.x = t_axis;
hold(fig.totalpvs_changes.ax)

fig.totalpvs_changes.plot_line(pax_fwhm.thickness.bvchanges,clee.clist.magenta);
fig.totalpvs_changes.plot_line(pax_fwhm.thickness.epschanges,clee.clist.coloredorange)
fig.totalpvs_changes.plot_line(pax_fwhm.thickness.pvschanges_total,clee.clist.darkgreen);

maxyrange = max([pax_fwhm.thickness.bvchanges,pax_fwhm.thickness.epschanges,pax_fwhm.thickness.pvschanges_total])*pixel2um;
minyrange = min([pax_fwhm.thickness.bvchanges,pax_fwhm.thickness.epschanges,pax_fwhm.thickness.pvschanges_total])*pixel2um;
ylim(fig.totalpvs_changes.ax,[minyrange,maxyrange])
fig.totalpvs_changes.put_yaxistitle('Length (\mum)');
fig.totalpvs_changes.put_xaxistitle('Time (sec)');
fig.totalpvs_changes.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);

%% 6. Dynamic and Static PVS thickness changes
fig.dspvs_changes = make_fig('pvsDynamicandStaticthicknesschages_figure');
hold(fig.dspvs_changes.ax)
fig.dspvs_changes.reset_axis
fig.dspvs_changes.update_figsize([8 3]);
fig.dspvs_changes.resolution = pixel2um;
fig.dspvs_changes.loc.x = t_axis;
fig.dspvs_changes.plot_line(pax_fwhm.thickness.bvchanges,clee.clist.magenta);
fig.dspvs_changes.plot_line(pax_fwhm.thickness.pvschanges_dynamic,clee.clist.darkgreen,'none','-');
fig.dspvs_changes.plot_line(pax_fwhm.thickness.pvschanges_static,clee.clist.lightgreen,'none','-');
maxyrange = max([pax_fwhm.thickness.bvchanges,pax_fwhm.thickness.pvschanges_dynamic,pax_fwhm.thickness.pvschanges_static])*pixel2um;
minyrange = min([pax_fwhm.thickness.bvchanges,pax_fwhm.thickness.pvschanges_dynamic,pax_fwhm.thickness.pvschanges_static])*pixel2um;
ylim(fig.dspvs_changes.ax,[minyrange,maxyrange])
fig.dspvs_changes.put_yaxistitle('Thickness changes (um)');
fig.dspvs_changes.put_xaxistitle('Time (sec)');
fig.dspvs_changes.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);

%% 7. Displacement of dynamic boundary
fig.dynamicdisplacement = make_fig('pvsDynamicDisplacement_figure');
fig.dynamicdisplacement.reset_axis
fig.dynamicdisplacement.update_figsize([8 3]);
fig.dynamicdisplacement.resolution = pixel2um;
fig.dynamicdisplacement.loc.x = t_axis;
hold(fig.dynamicdisplacement.ax)
fig.dynamicdisplacement.plot_line(pax_fwhm.displacement.dynamicpvs,clee.clist.orange);
fig.dynamicdisplacement.plot_line(pax_fwhm.displacement.dynamicbv,clee.clist.magenta);
fig.dynamicdisplacement.put_yaxistitle('Length (\mum)');
fig.dynamicdisplacement.put_xaxistitle('Time (sec)');
fig.dynamicdisplacement.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);

%% 8. Displacement of static boundary
fig.staticdisplacement = make_fig('pvsStaticDisplacement_figure');

fig.staticdisplacement.reset_axis
fig.staticdisplacement.update_figsize([8 3]);
fig.staticdisplacement.resolution = pixel2um;
fig.staticdisplacement.loc.x = t_axis;
hold(fig.staticdisplacement.ax)
fig.staticdisplacement.plot_line(pax_fwhm.displacement.staticpvs,clee.clist.darkgreen);
fig.staticdisplacement.plot_line(pax_fwhm.displacement.staticbv,clee.clist.magenta);
fig.staticdisplacement.put_yaxistitle('Length (\mum)');
fig.staticdisplacement.put_xaxistitle('Time (sec)');
fig.staticdisplacement.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);
%% 9. Comparison of static and dynamic pvs displacement
fig.pvsdisplacement = make_fig('pvsDisplacement_figure');

fig.pvsdisplacement.reset_axis
fig.pvsdisplacement.update_figsize([8 3]);
fig.pvsdisplacement.resolution = pixel2um;
fig.pvsdisplacement.loc.x = t_axis;
hold(fig.pvsdisplacement.ax)
fig.pvsdisplacement.plot_line(pax_fwhm.displacement.staticpvs,clee.clist.darkgreen);
fig.pvsdisplacement.plot_line(pax_fwhm.displacement.dynamicpvs,clee.clist.orange);
fig.pvsdisplacement.put_yaxistitle('Length (\mum)');
fig.pvsdisplacement.put_xaxistitle('Time (sec)');
fig.pvsdisplacement.save2svg(save_dir);
disp(['Saved paxPVS_figure to ', save_dir]);
end
