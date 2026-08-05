function polar = piv_polar_events(piv_run, coremask, opt)
%PIV_POLAR_EVENTS  Perivascular divergence maps for every event in piv_run.
%   One settings set for all events, so the maps stay comparable. Two measures:
%     radial (cells/pdiv)  mean v_r per (ring,sector), then dv_r/dr along the radius.
%                          < 0 = PVS compressed, > 0 = PVS refilling.
%     plain  (div)         cartesian du/dx+dv/dy from a block plane fit. Mixes in
%                          tangential shear -- the cross-check, not the headline.
%   uv is already a total displacement, so both come out as areal strain directly.
%
% IN   piv_run    output of piv_run_events
%      coremask   logical HxW blob rings start at (manual_dilated_pvs). Measured
%                 from its wall (bwdist), sectors from its centroid
%      exclmask, n_sectors (8), min_valid (4, both partitions),
%      box_step (3, box side in vector-grid steps), verbose (false)
% OUT  polar      1xN cell of struct: cells, rcen, vr, vt (+ = counter-clockwise),
%                 pdiv, vrmap, vtmap, div, info

arguments
    piv_run struct
    coremask logical
    opt.exclmask  logical = []
    opt.n_sectors (1,1) double = 8
    opt.min_valid (1,1) double = 4
    opt.box_step  (1,1) double = 3
    opt.verbose   (1,1) logical = false
end

% 0. Setup
% 0.1 The blob the rings start at
if ~any(coremask, 'all')
    error('piv_polar_events:noCore', 'coremask is empty (manual_dilated_pvs is blank).');
end
% 0.2 Vector-grid step: both cell shapes are sized from it, and it is a property
%     of the PIV run, not of the geometry, so it is measured once here
step = piv_grid_step(piv_run(1).uv);

% 0.3 Box partition -- geometry only, so it is built once and reused for every
%     event. The wedge partition is built inside piv_divergence_polar
[tiles, tinfo] = sector_block(coremask, opt.box_step * step);
if ~isempty(opt.exclmask); tiles(opt.exclmask) = 0; end

% 1. Both divergence measures, event by event
polar = cell(1, numel(piv_run));
for k = 1:numel(piv_run)
    uv = piv_run(k).uv;                 % already a total displacement -- no * nfr
    % 1.1 Radial: ring/sector v_r, then dv_r/dr
    [pdiv, cells, pinfo] = piv_divergence_polar(uv, coremask, ...
        'ring_width', 2*step, 'n_sectors', opt.n_sectors, ...
        'min_valid', opt.min_valid, 'exclmask', opt.exclmask, 'method', 'radial');
    % 1.2 Plain: cartesian du/dx + dv/dy over equal-area boxes
    tk  = sector_validate(tiles, ~isnan(uv(:,:,1)) & ~isnan(uv(:,:,2)), ...
        tinfo.n_labels, opt.min_valid);
    div = sector_paint(tk, piv_divergence(uv, tk, tinfo.n_labels));
    polar{k} = struct('cells', cells, 'rcen', pinfo.rcen, ...
        'vr', pinfo.vr_cells, 'vt', pinfo.vt_cells, 'pdiv', pdiv, ...
        'vrmap', pinfo.vr_map, 'vtmap', pinfo.vt_map, 'div', div, 'info', pinfo);
    if opt.verbose; fprintf('.'); end
end
if opt.verbose; fprintf(' %d events\n', numel(piv_run)); end
end
