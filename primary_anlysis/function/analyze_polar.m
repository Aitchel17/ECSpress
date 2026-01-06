function [polarstruct,radius_map,centered_anglemap,angleid_map] = analyze_polar(Stack, Centervertice,Start_angle,n_angles,equidistance_mode)
%ANALYZE_POLAR Summary of this function goes here
%   Detailed explanation goes here
    % Input: Stack, Centervertice,Start_angle,equidistance_mode
    % Stack  (x,y,t)
    % Centervertices 1x2 double (x,y)
    % Start_angle 1x1 double , made from atan2d
    % Equidistance mode logical, True: equidistance radius binning
        % False: equalarea radius binning
% Core steps
% 1. Polar coodrination mask generation
    % 1.1 Using first frame to generate grid array of x and y
    % 1.2 Make center point (0,0) by subtracting each center coordinate
        % ^^^^^^^^^^Center vertices used in this step^^^^^^^^^^
    % 1.3 Get Polar mask of angle and radius (pixels)
            % cart2pol from step 1.1 provide flipped signed radian map 
% 2. Anglemask modification
    % 2.1 Make Anglemask from radian to degree
    % 2.2 Signed angle 0, -179 , 180 ,0 to Azimuth 0, 180, 180, 360 
    % 2.3 Calculate center angle from signed angle to Azimuth
        % ^^^^^^^^ Start angle used in this step ^^^^^^^^
    % 2.4 Recenter Anglemask

% 3. Get Angleid map and id list
    % 3.1 Divide by angle_range (angle bin size) and floor
    % 3.2 Get idlist
% 4.
%% Hard coded parameter
angle_range = 360/n_angles; % 30 12 section
bin_pixel = 15; % 15 pixels for binedge calculation
bin_distance = 2; % 5 pixel distance for equidistance binedge calculation
nframes = size(Stack,3);
%% Initialize data container
polarstruct = struct('angle_id', {}, 'kymograph', {}, 'radius_bin', {},'angle_range', {});

%% 1. Polar coodrination generation
[meshx,meshy] = meshgrid(1:size(Stack(:,:,1),2),1:size(Stack(:,:,1),1)); % 1.1
meshx = meshx - Centervertice(1); % 1.2
meshy = meshy - Centervertice(2); % 1.2
[theta_map, radius_map] = cart2pol(meshx,meshy); % 1.3
%% 

%% 2. Anglemask Modification
theta_map = rad2deg(theta_map); % 2.1
theta_map = mod(-theta_map,360); % 2.2
bin_start = mod(-Start_angle,180) - angle_range/2; % 2.3
centered_anglemap = mod(theta_map - bin_start, 360); % 2.4
%% 3. Get Angleid
angleid_map = floor( centered_anglemap / angle_range ) + 1; % 3.1
angle_idlist = unique(angleid_map); % 3.2
%% 4. Access to each segment
collapsed_stack = reshape(Stack, [], nframes);  % (H*W) x T % 4.1
for angle_idx = angle_idlist' % 4.2
    polarstruct(angle_idx).angle_id = angle_idx; % 4.2.1
    angle_mask = angleid_map == angle_idx; % 4.2.2
    polarstruct(angle_idx).angle_range = [min(theta_map(angle_mask)), max(theta_map(angle_mask))]; % 4.2.3
    %% 5. Radius
    radius_idxvalue = radius_map(angle_mask); 
    if equidistance_mode
        start_radius = rem(max(radius_idxvalue),bin_distance);
        radius_nbins = floor(max(radius_idxvalue)/bin_distance);
        radius_edges = linspace(bin_distance+start_radius,max(radius_idxvalue),radius_nbins);
    else % otherwise equal area
        radius_idxvalue = sort(radius_idxvalue,'ascend');
        radius_nbins = floor(length(radius_idxvalue)/bin_pixel);
        radius_edges = radius_idxvalue(bin_pixel*(1:radius_nbins));
    end
    radius_edges = [0, radius_edges];
    polarstruct(angle_idx).radius_bin = radius_edges(2:end)-radius_edges(1:end-1);
    %% 6. Preallocate
    polarstruct(angle_idx).kymograph = zeros([radius_nbins, nframes]);
    %% 7. Access to each Radius
    for radius_idx = 1:radius_nbins
        % the final submask
        submask = (angleid_map == angle_idx) & (radius_map > radius_edges(radius_idx)) & (radius_map <= radius_edges(radius_idx+1));
        % apply final submask to get average values over time
        M = reshape(submask,[],1);
        %
        ts = mean(collapsed_stack(M, :), 1, 'omitnan');  % 1 x T, 프레임별 평균
        polarstruct(angle_idx).kymograph(radius_idx,:) = squeeze(ts);
    end

end


end

