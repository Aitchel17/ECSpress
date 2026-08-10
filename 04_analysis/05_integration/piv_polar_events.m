function [profile, fields] = piv_polar_events(piv_run, coremask, pixel2um, opt)
%PIV_POLAR_EVENTS  Read every event's displacement field in polar coordinates.
%   One settings set for all events, so the rows stay comparable, and one
%   vfield_profile holding them all against each other.
%
%   WHAT CHANGED, AND WHY THE OUTPUT IS NOT THE OLD ONE.
%   This used to return a cell of dense maps built from two estimators at once:
%   a ring/sector dv_r/dr, and a cartesian du/dx+dv/dy over equal-area boxes as a
%   cross-check. Both are gone.
%
%     the ring/sector estimator averaged over a wedge set that changed with the
%     ring, because a wedge that runs off the frame stops sooner than its
%     neighbours. vfield_polar bins the triangles instead, and vfield_profile
%     refuses to pool a wedge that does not span the whole comparison interval
%
%     dv_r/dr is ONE COMPONENT of the strain, and its null is not zero: under
%     volume conservation strain_radial = -disp_radial/r, so it is already
%     negative before anything is absorbed. The trace is the quantity with a
%     zero null, and it is what divergence means
%
%     the box cross-check used a per-cell least-squares plane, so the set being
%     fitted changed with cell coverage -- the same fault in a second place
%
%   The order below is spelled out rather than bundled into one call. Which stage
%   runs when belongs to vfield_polar, and a helper that hid it would take that
%   away from the class that has to own it.
%
% IN   piv_run    output of piv_run_events. uv is already a TOTAL displacement
%      coremask   H x W bool, traced at maximum dilation. Radius is measured from
%                 its wall by bwdist, wedges from its centroid
%      pixel2um   float, microns per pixel. Required now: the bins are microns,
%                 so that rows from different recordings can meet
%      n_wedge        int, angular bins
%      bin_edges_um   1 x nB+1 float, radial bin edges in MICRONS from the wall
%      exclmask       H x W bool, true = never interpolate (the PIV mask)
%      min_tri_wedge  int, triangles a wedge needs before it is used at all
%      gated          bool, whether piv_run's uv came through the PIV gates.
%                     Required. The gates take a radius-dependent share of the
%                     LARGER vectors, so a gated row and an ungated one are not
%                     the same measurement. See CLAUDE_LOG.md
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
    opt.gated         (1,1) logical
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
    vp = vfield_polar(piv_run(k).uv, coremask, pixel2um, ...
        n_wedge = opt.n_wedge, ...
        bin_edges_um = opt.bin_edges_um, ...
        exclmask = opt.exclmask, ...
        gated = opt.gated);
    vp.param.min_tri_wedge = opt.min_tri_wedge;
    vp.param.verbose = false;

    % measure -> gate_wedge -> accumulate. accumulate takes the verdict as a
    % required argument, so a gate cannot be formed here and then forgotten
    vp.measure();
    wedge_mask = vp.gate_wedge();
    cells = vp.accumulate(wedge_mask);

    profile.append(cells, ...
        event_idx = piv_run(k).id, ...
        pol = string(piv_run(k).pol), ...
        state = string(piv_run(k).state), ...
        dD = piv_run(k).diameter_change, ...
        gated = opt.gated);
    fields{k} = vp;

    if opt.verbose
        fprintf('.');
    end
end

if opt.verbose
    fprintf(' %d events\n', n_event);
end
end
