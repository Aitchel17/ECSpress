% roilist = roilist.addroi(preprocessed_ch1,'test','line');
% roilist = roilist.modifyroi(preprocessed_ch1,'test');
% roilist = roilist.removeroi('test');

%% load stack, vertices from main_primary script
bv_stack = preprocessed_ch1;
pax_vertices = roilist.getvertices('pax');
pax_vertices = pax_vertices(1:2,:);
%% Primary axis used for 1d calculation
pax_angle = pax_vertices(2,:)-pax_vertices(1,:);
pax_angle = atan2d(pax_angle(2),pax_angle(1));

%% bin angle
angle_range = 30;
% section area
bin_pixel = 15;
%% calculate centerposition by projecting center point to PAX line
pax_vertices = roilist.getvertices('pax'); %(1:2,:);
pax_vertices = pax_vertices(1:2,:);
exp_center = (min(roilist.getvertices('extraparenchyma'),[],1) + max(roilist.getvertices('extraparenchyma'),[],1))/2;
pax_center = analyze_dpoint2line(exp_center,pax_vertices);
%% show center position
roilist.showroi('extraparenchyma')
hold on
ax = gca;
plot(ax,[pax_vertices(:, 1); pax_vertices(1, 1)], [pax_vertices(:, 2); pax_vertices(1, 2)], 'r-', 'LineWidth', 2);
plot(ax, exp_center(1), exp_center(2), 'g+', 'MarkerSize', 10, 'LineWidth', 1);
plot(ax, pax_center(1), pax_center(2), 'g+', 'MarkerSize', 10, 'LineWidth', 1);
%%
x = analyze_polar(bv_stack,pax_center,pax_angle,true);
%%
figure()
imagesc(x(6).kymograph)



%% bring center position to cropped mask
exp_center = exp_center-min(roilist.getvertices('extraparenchyma'),[],1);
pax_center = pax_center-min(roilist.getvertices('extraparenchyma'),[],1);
pax_vertices = pax_vertices-min(roilist.getvertices('extraparenchyma'),[],1);
%% show center position on cropped mask
fig = figure();
ax  = axes('Parent',fig);
imshow(tmp.nonnan2d, 'Parent', ax);
hold on
ax = gca;
plot(ax,[pax_vertices(:, 1); pax_vertices(1, 1)], [pax_vertices(:, 2); pax_vertices(1, 2)], 'r-', 'LineWidth', 2);
plot(ax, exp_center(1), exp_center(2), 'g+', 'MarkerSize', 10, 'LineWidth', 2);
plot(ax, pax_center(1), pax_center(2), 'g+', 'MarkerSize', 10, 'LineWidth', 2);
%% polar coordination
% input mask, output polar coordinate with origin from center of mask
% center coordinate of PVS
bv_holdstack = roi_applyvertices(bv_stack,roilist.getvertices('extraparenchyma'));

nframes = size(bv_stack,3);
tmp.nonnan2d = ~isnan(bv_stack(:,:,1));
[~,center.loc_x] = max(sum(tmp.nonnan2d,2));
[~,center.loc_y] = max(sum(tmp.nonnan2d,1));
[tmp.meshx,tmp.meshy] = meshgrid(1:size(tmp.nonnan2d,2),1:size(tmp.nonnan2d,1));
tmp.meshx = tmp.meshx - center.loc_x;
tmp.meshy = tmp.meshy - center.loc_y;
% get polar coordination
[theta_map, radius_map] = cart2pol(tmp.meshx,tmp.meshy);
theta_map_deg = rad2deg(theta_map);
theta_map_deg = mod(-theta_map_deg,360); % clockwise

%%
pax_anglemap = mod(theta_map_deg - bin_start, 360);
%% divide and floor to get angular section
angleid_map = floor( pax_anglemap / angle_range ) + 1;
angle_idlist = unique(angleid_map);
%%
idx = find(angleid_map ==1); % pixel position corresopond to angular id 1
radius_idxvalue = radius_map(idx);
radius_idxvalue = sort(radius_idxvalue,'ascend');
%%
radius_nbins = floor(length(radius_idxvalue)/bin_pixel);
%% radius edge
radius_edges = zeros(1, radius_nbins+1);
radius_edges(1) = 0;
for bincount = 1:radius_nbins
    radius_edges(1+bincount) = radius_idxvalue(bincount*bin_pixel); 
end 
%%
wedge_t = zeros([radius_nbins, nframes]);
%%
clc
for radius_idx = 1:radius_nbins
    % the final submask
    submask = (angleid_map == 1) & (radius_map > radius_edges(radius_idx)) & (radius_map <= radius_edges(radius_idx+1));
    % apply final submask to get average values over time
    X = reshape(bv_holdstack, [], nframes);  % (H*W) x T
    M = reshape(submask,[],1);
    %
    ts = mean(X(M, :), 1, 'omitnan');  % 1 x T, 프레임별 평균
    wedge_t(radius_idx,:) = squeeze(ts);
end
%%
util_checkstack(bv_holdstack)

%%
t = [1:nframes];
%%
figure()
surface(t,radius_edges(2:end),zeros(size(wedge_t)),(wedge_t),'LineStyle','none');
q = gca;
q.Layer = 'top'; % put the axes/ticks on the top layer
%%
figure()
imagesc(line_fwhms.pax.kymograph.kgph_bv)
%%
set(gca,'YDir','normal')
%%
figure()
plot(radius_edges(2:end),wedge_t(:,2500))
