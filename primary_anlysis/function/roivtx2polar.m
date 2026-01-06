function theta_radius = roivtx2polar(vertex,angle_map,radius_map)
% ROIVTX2POLAR Converts ROI vertices to polar coordinates using pre-calculated maps.
%   vertex: [x, y] coordinates of vertices/pixels.
%   angle_map: 2D array of angles (result from analyze_polar).
%   radius_map: 2D array of radii (result from analyze_polar).
%
%   Returns:
%       theta_radius: [3 x N] matrix
%           Row 1: Angle in Radians
%           Row 2: Angle in Degrees
%           Row 3: Radius

imgsize = size(angle_map);
% Ensure vertex indices are within bounds
vertex = round(vertex);
idx = sub2ind(imgsize, vertex(:,2), vertex(:,1));

angle = angle_map(idx);
radius = radius_map(idx);

[deg_angle, order] = sort(angle,'ascend');
radius = radius(order);
angle = deg2rad(deg_angle)';
radius = radius';
theta_radius = [angle ; deg_angle' ;radius];
end
