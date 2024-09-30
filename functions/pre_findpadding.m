function [start_x,end_x,mean_x] = pre_findpadding(mobj,img_ch,tpsm_info,start_z)
%% find padding caused by sinusoidal correction
frames = io_readframes(mobj,img_ch,[1,tpsm_info.fw],[1,tpsm_info.fh],[start_z,start_z+99]); % load first 100 frames
mean_x = mean(frames,[1,3]); % calculate mean value of y, z axis (y,x,z)
tmp.nzloc = find(mean_x~=-2048); % find location of value not -2048

start_x = tmp.nzloc(1); % start point of non zero
end_x = tmp.nzloc(end); % end point of non zero

end

