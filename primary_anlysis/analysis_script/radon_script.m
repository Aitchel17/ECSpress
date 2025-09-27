
%%



%%

%%
hold on
plot(ax, pax_center(1), pax_center(2), 'g+', 'MarkerSize', 10, 'LineWidth', 2);
plot(ax,[pax_vertices(:, 1); pax_vertices(1, 1)], [pax_vertices(:, 2); pax_vertices(1, 2)], 'r-', 'LineWidth', 2);

%% util_checkstack(bv_stack)

%%
bv_holdstack = roi_applyvertices(bv_stack,roilist.vertices.extraparenchyma);
%%
util_checkstack(bv_holdstack);
%% polar coordination
% input mask, output polar coordinate with origin from center of mask
% center coordinate of PVS
tmp.nonnan2d = ~isnan(bv_holdstack(:,:,1));
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



%%

%%
analyze_radon();