function interpolated_sinogram = sinogram_interpolation(original_sinogram)
%SINOGRAM_INTERPOLATION Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = original_sinogram;
%%
radius = original_sinogram(end,:);
angle = original_sinogram(2,:);
%% Extended radius and angle to avoid padding issue at the edge
ex_radius = repmat(radius,1,3);
ex_angle = [angle-360,angle,angle+360];
%%  Interpolate
interp_angle = [-90:1:450];
interp_radius = interp1(ex_angle,ex_radius,interp_angle,"makima");

interpolated_sinogram = zeros([2,361]);
interpolated_sinogram(1,:) =interp_angle(91:451);
interpolated_sinogram(2,:) = interp_radius(91:451);


%% Debugging 
% cla(gca)
% plot(angle,radius,'x')
% hold on 
% plot(ex_angle,ex_radius,'-r')
% plot(interp_angle,interp_radius,'-ko');

end

