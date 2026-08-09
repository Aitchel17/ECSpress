function frame = polar_frame(img_size, center, opt)
%POLAR_FRAME  Shared polar frame for cross-analysis co-registration.
%   The one polar coordinate system that clusterpolar (02), radon (03) and PIV (04)
%   all reference, so their per-sector outputs land in the same frame. Frames with
%   matching img_size and center stay reconcilable whatever pax_angle is: theta_raw
%   is pax-independent and bin_start_deg records the offset, so a differing PAX
%   alignment is only a sector-label rotation. Rings are deliberately not binned
%   here - each analysis owns its own radial scheme.
%
% IN   img_size   [H W] pixel size of the image the maps span
%      center     [x y] polar origin (px), e.g. constricted-BV centroid
%      pax_angle  reference axis angle (deg), default 0 = raw image axes
%      n_sectors  angular sector count, default 24 (matches clusterpolar)
% OUT  frame      struct; img_size, center, pax_angle echoed for provenance, plus
%                 .angle_range    sector width (deg)
%                 .radius_map     HxW center distance (px), NOT a wall distance;
%                                 PIV owns its own bwdist shells
%                 .angle_map      HxW PAX-aligned azimuth (deg) in [0,360)
%                 .sector_idx     HxW sector id 1..n_sectors
%                 .theta_raw      HxW azimuth before recentering (deg), the
%                                 pax-independent reconciliation key
%                 .bin_start_deg  screen angle = angle_map + bin_start_deg

arguments
    img_size      (1,2) double {mustBePositive}
    center        (1,2) double
    opt.pax_angle (1,1) double = 0
    opt.n_sectors (1,1) double {mustBePositive, mustBeInteger} = 24
end

% 0. Image size and sector width
H = img_size(1);  W = img_size(2);
angle_range = 360 / opt.n_sectors;

% 1. Polar coordinates about the center (analyze_polar step 1)
[meshx, meshy] = meshgrid(1:W, 1:H);
meshx = meshx - center(1);
meshy = meshy - center(2);
[theta, radius_map] = cart2pol(meshx, meshy);      % 1.1 signed rad, center distance (px)

% 2. Azimuth recentered on the PAX axis (analyze_polar step 2)
theta_raw = mod(-rad2deg(theta), 360);             % 2.1 azimuth (deg) in [0,360)
bin_start = mod(-opt.pax_angle, 180) - angle_range / 2;   % 2.2 recentering offset
angle_map = mod(theta_raw - bin_start, 360);       % 2.3 raw axes when pax_angle = 0

% 3. Sector id (analyze_polar step 3)
sector_idx = floor(angle_map / angle_range) + 1;   % 3.1 1..n_sectors
sector_idx = min(max(sector_idx, 1), opt.n_sectors);      % 3.2 clamp to range

% 4. Pack the frame
frame = struct( ...
    'img_size',      img_size, ...
    'center',        center, ...
    'pax_angle',     opt.pax_angle, ...
    'n_sectors',     opt.n_sectors, ...
    'angle_range',   angle_range, ...
    'radius_map',    radius_map, ...
    'angle_map',     angle_map, ...
    'sector_idx',    sector_idx, ...
    'theta_raw',     theta_raw, ...
    'bin_start_deg', bin_start);
end
