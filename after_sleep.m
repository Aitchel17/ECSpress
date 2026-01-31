sessiondir = 'G:\tmp\00_igkl\hql088\250927_hql088_sleep\HQL088_sleep250927_010';
session = ECSSession(sessiondir);
session = session.load_primary_results();
%
session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
% Two photon analysis time axis calculation
session.pax_fwhm.t_axis = linspace(session.img_param.imgstarttime,session.img_param.imgendtime,numel(session.pax_fwhm.idx.upperBVboundary));
sleep_score = load(fullfile(sessiondir,"peripheral","sleep_score.mat"));
t_axis = linspace(session.img_param.imgstarttime,session.img_param.imgendtime,numel(session.pax_fwhm.idx.upperBVboundary));
% Calculate
% Find indices for all REM epochs
fwhmloc.rem = statebin2frame(sleep_score.REMTimes,t_axis);
fwhmloc.nrem = statebin2frame(sleep_score.NREMTimes,t_axis);
fwhmloc.awake = statebin2frame(sleep_score.NREMTimes,t_axis);
fwhmloc.drowsy = statebin2frame(sleep_score.DrowsyTimes,t_axis);
sleep_score.behavState
% Image frame location
imgloc.rem = statebin2frame(sleep_score.REMTimes,session.img_param.taxis);
imgloc.nrem = statebin2frame(sleep_score.NREMTimes,session.img_param.taxis);
imgloc.awake = statebin2frame(sleep_score.NREMTimes,session.img_param.taxis);
imgloc.drowsy = statebin2frame(sleep_score.DrowsyTimes,session.img_param.taxis);
%%

fig = figure('Name','Thickness plot');
t = tiledlayout(1,1);
%%
ax.bvthickness = axes(t);
%%
plot(ax.bvthickness,t_axis(fwhmloc.nrem),session.pax_fwhm.thickness.pvschanges_dynamic(fwhmloc.nrem))






%%

channel = ["bv","csf"];
state = ["rem","nrem","awake"];

%%
fig = figure();
fig.Units = 'inches';
fig.Position = [21 1 18 7]; % [left bottom width height]
t = tiledlayout(numel(channel),numel(state), 'TileSpacing', 'none', 'Padding', 'none');

img_handles = gobjects(numel(channel),numel(state));
%
ax = struct();
for ch =1:numel(channel)
    for st = 1:numel(state)
        fieldname =  strcat(channel(ch),"_",state(st));
        % Use nexttile to get the next tile in the layout
        ax.(fieldname) = nexttile(t);
        current_locs = imgloc.(state(st));
        if isempty(current_locs)
            title(ax.(fieldname), [fieldname " (Empty)"], 'Interpreter', 'none');
            continue;
        end

        if ch == 1
            img_data = mean(session.stackch1(:,:,current_locs), 3);
        elseif ch == 2
            img_data = mean(session.stackch2(:,:,current_locs), 3);
        end

        % Contrast adjustment: 10th to 90th percentile
        clim_val = prctile(img_data(:), [10 99]);
        img_handles(ch, st) = imagesc(ax.(fieldname), img_data, clim_val);
        axis(ax.(fieldname), 'image'); % Fix aspect ratio
        title(ax.(fieldname), fieldname, 'Interpreter', 'none');
    end
end
%% Sync all axes for simultaneous panning/zooming
all_axes = struct2cell(ax);
all_axes = [all_axes{:}]; % Convert cell array to object array
linkaxes(all_axes, 'xy');
%%
%% Interactive Marker Control (Click on Figure)
fprintf('--- Interactive Marker Mode ---\n');
fprintf('Click on any image to mark position on all axes.\n');
fprintf('Press Enter (key) to stop.\n');

marker_handles = [];
while true
    try
        [x_val, y_val, button] = ginput(1);
    catch
        break; % Handle figure closure
    end

    % Break if Enter (empty button) or any key press (button > 3 typically)
    % Mouse clicks are 1 (Left), 2 (Middle), 3 (Right)
    if isempty(button) || ~ismember(button, [1, 2, 3])
        fprintf('Exiting marker mode.\n');
        break;
    end

    % Remove old markers
    delete(marker_handles);
    marker_handles = gobjects(0);

    % Add new markers to all axes
    for i = 1:numel(all_axes)
        try
            curr_ax = all_axes(i);
            hold(curr_ax, 'on');
            h = plot(curr_ax, x_val, y_val, 'rx', 'MarkerSize', 8, 'LineWidth', 1);
            marker_handles(end+1) = h;
        catch
            % Handle closed axes
        end
    end
    drawnow;
end

%%

imagesc(ans)
axis image
%%




%%
session.pax_fwhm.kymograph.kgph_lumen_processed(:,drowsy_locs)


%%
clee = color_lee;
%%
thickness_fig = make_fig('Thickness figure');
%%
thickness_fig.hold_axis(true)
thickness_fig.resolution = session.img_param.pixel2um;
thickness_fig.update_figsize([15,3])
thickness_fig.monitor_xyinch = [21, 5];
thickness_fig.reset_axis
thickness_fig.loc.x = t_axis;
xlim(thickness_fig.ax,[0 max(t_axis)])
thickness_fig.plot_line(session.pax_fwhm.thickness.bvchanges,'r')
thickness_fig.plot_line(session.pax_fwhm.thickness.pvschanges_total, clee.clist.darkgreen)
thickness_fig.put_xaxistitle('Time (sec)')
thickness_fig.put_yaxistitle('Thickness change (um)')
%%
kymograph_fig = make_fig('Kymograph figure');
kymograph_fig.resolution = session.img_param.pixel2um;
kymograph_fig.update_figsize([15,3])
kymograph_fig.monitor_xyinch = [21, 5];
kymograph_fig.reset_axis
kymograph_fig.loc.x = t_axis;

kymograph_fig.plot_kymograph(session.pax_fwhm.kymograph.kgph_pvs_processed(:,fwhmloc.rem),t_axis(fwhmloc.rem))
xlim(kymograph_fig.ax,[min(t_axis(fwhmloc.rem)) max(t_axis(fwhmloc.rem))])


%%
scope = fwhmloc.nrem;
kymograph_fig = make_fig('Kymograph figure');
kymograph_fig.resolution = session.img_param.pixel2um;
kymograph_fig.update_figsize([15,3])
kymograph_fig.monitor_xyinch = [21, 5];
kymograph_fig.reset_axis
kymograph_fig.loc.x = t_axis;

kymograph_fig.plot_kymograph(session.pax_fwhm.kymograph.kgph_pvs_processed(:,scope),t_axis(scope))
xlim(kymograph_fig.ax,[min(t_axis(scope)) max(t_axis(scope))])







