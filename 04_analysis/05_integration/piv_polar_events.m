function [profile, fields] = piv_polar_events(piv_run, coremask, pixel2um, opt)
%PIV_POLAR_EVENTS  Read every event's displacement field in polar coordinates.
%   One settings set for all events, so the rows stay comparable, and one
%   vfield_profile holding them all against each other. The stages are called in
%   order here, not through a helper. It replaced a ring/sector dv_r/dr estimator
%   and a box cross-check in 2026-08; why both went : see CLAUDE_LOG.md
%
% IN   piv_run    output of piv_run_events. uv is already a TOTAL displacement.
%                 Each row's gate_name travels into the profile, which refuses to
%                 pool rows gated differently
%      coremask   H x W bool, traced at maximum dilation. Radius is measured from
%                 its wall by bwdist, wedges from its centroid
%      pixel2um   float, microns per pixel; the bins are microns
%      n_wedge        int, angular bins
%      bin_edges_um   1 x nB+1 float, radial bin edges in MICRONS from the wall
%      exclmask       H x W bool, true = never interpolate (the PIV mask)
%      min_tri_wedge  int, triangles a wedge needs before it is used at all
% OUT  profile    vfield_profile, one row per event, ready for difference()
%      fields     1 x N cell of vfield_polar, kept for the display layer
%
%   see also VFIELD_POLAR, VFIELD_PROFILE, PIV_RUN_EVENTS, PIV_FIGURE

arguments
    piv_run  struct
    coremask logical
    pixel2um (1,1) double {mustBePositive}
    opt.n_wedge       (1,1) double {mustBePositive, mustBeInteger} = 12
    opt.bin_edges_um  (1,:) double = 0:1.5:40
    opt.exclmask      logical = []
    opt.min_tri_wedge (1,1) double = 10
    opt.verbose       (1,1) logical = false
end

if ~any(coremask, 'all')
    error('piv_polar_events:noCore', 'coremask is empty (manual_dilated_pvs is blank).');
end

n_event = numel(piv_run);
profile = vfield_profile(opt.bin_edges_um, opt.n_wedge);
profile.param.verbose = opt.verbose;
fields = cell(1, n_event);

for k = 1:n_event
    vp = vfield_polar(piv_run(k).xyuv, coremask, pixel2um, ...
        n_wedge = opt.n_wedge, ...
        bin_edges_um = opt.bin_edges_um, ...
        exclmask = opt.exclmask, ...
        gate_name = string(piv_run(k).gate_name));
    vp.param.min_tri_wedge = opt.min_tri_wedge;
    vp.param.verbose = false;

    % measure -> gate_wedge -> accumulate; accumulate takes the verdict as an argument
    vp.applydelaunay();
    vp.trifilter();
    vp.placetri();
    vp.measure();
    wedge_mask = vp.gate_wedge();
    cells = vp.accumulate(wedge_mask);

    profile.append(cells, ...
        event_idx = piv_run(k).id, ...
        pol = string(piv_run(k).pol), ...
        state = string(piv_run(k).state), ...
        dD = piv_run(k).diameter_change, ...
        gate_name = string(piv_run(k).gate_name));
    fields{k} = vp;

    if opt.verbose
        fprintf('.');
    end
end

if opt.verbose
    fprintf(' %d events\n', n_event);
end
end
