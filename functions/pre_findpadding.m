function [start_x,end_x] = pre_findpadding(frames)
%% find padding caused by sinusoidal correction
mean_x = mean(frames,[1,3]); % calculate mean value of y, z axis (y,x,z)
tmp.nzloc = find(mean_x~=-2048); % find location of value not -2048
start_x = tmp.nzloc(1); % start point of non zero
end_x = tmp.nzloc(end); % end point of non zero
end

