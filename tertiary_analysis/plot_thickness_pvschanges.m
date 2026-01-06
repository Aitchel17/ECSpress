uppvs_changes = zeros([size(secondary_struct,2),1]);
uppvs_thickness = zeros([size(secondary_struct,2),1]);

for ssidx = 1:size(secondary_struct,2)
    scale = strsplit(secondary_struct(ssidx).infodict("objpix"));
    scale = str2double(scale(1));
    idxarray = strcmp({secondary_struct(ssidx).heatdata.type}, 'uppvs');
    tmp.pvs_thicknessarray = secondary_struct(ssidx).thickness(idxarray).thickness; % ('uppvs')
    uppvs_thickness(ssidx) = max(tmp.pvs_thicknessarray*scale);
    uppvs_changes(ssidx) = max(secondary_struct(ssidx).heatdata(idxarray).modespvs)-min(secondary_struct(ssidx).heatdata(idxarray).modespvs);
end
%%

downpvs_changes = zeros([size(secondary_struct,2),1]);
downpvs_thickness = zeros([size(secondary_struct,2),1]);

for ssidx = 1:size(secondary_struct,2)
    scale = strsplit(secondary_struct(ssidx).infodict("objpix"));
    scale = str2double(scale(1));
    idxarray = strcmp({secondary_struct(ssidx).heatdata.type}, 'downpvs');
    tmp.pvs_thicknessarray = secondary_struct(ssidx).thickness(idxarray).thickness; % ('downpvs')
    downpvs_thickness(ssidx) = max(tmp.pvs_thicknessarray*scale);
    downpvs_changes(ssidx) = max(secondary_struct(ssidx).heatdata(idxarray).modespvs)-min(secondary_struct(ssidx).heatdata(idxarray).modespvs);
end
%%


%%
fig = figure()
hold on
scatter(uppvs_thickness,uppvs_changes,100,"x",'MarkerEdgeColor', 'w')
scatter(downpvs_thickness,downpvs_changes,100,"x",'MarkerEdgeColor', 'w')
plot([0 25], [0 25])
axis tight equal
ax = gca;              % Get current axes
ax.FontSize = 14;
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
fig.Color = [0 0 0];

%%
all_pvsthickness = cat(1,uppvs_thickness,downpvs_thickness);
all_pvschanges = cat(1,uppvs_changes,downpvs_changes);
[b2, stat] = robustfit(all_pvsthickness,all_pvschanges);

yfit = b2(1) + b2(2).*all_pvsthickness;

fig = figure()
hold on
scatter(all_pvsthickness,all_pvschanges,100,"x",'MarkerEdgeColor', 'w')
plot([0 25], [0 25],'Color','r')
plot(all_pvsthickness,yfit,'Color','w')


axis tight equal
ax = gca;              % Get current axes
ax.FontSize = 14;
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
fig.Color = [0 0 0];

%%
totalthickness = uppvs_thickness+downpvs_thickness;
eccentricity = abs(0.5-uppvs_thickness./totalthickness);
totalpvs_changes = uppvs_changes+downpvs_changes;
totalpvs_changes = totalpvs_changes./totalthickness;
[b2, stat] = robustfit(eccentricity,totalpvs_changes);
yfit = b2(1) + b2(2).*eccentricity;
fig = figure()
hold on
scatter(eccentricity,totalpvs_changes,100,"x",'MarkerEdgeColor', 'w')
%plot([0 25], [0 25],'Color','r')
plot(eccentricity,yfit,'Color','w')

axis tight equal
ax = gca;              % Get current axes
ax.FontSize = 14;
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
fig.Color = [0 0 0];

%% not normalized

totalthickness = uppvs_thickness+downpvs_thickness;
eccentricity = abs(totalthickness/2-uppvs_thickness);
totalpvs_changes = uppvs_changes+downpvs_changes;
totalpvs_changes = totalpvs_changes;
[b2, stat] = robustfit(eccentricity,totalpvs_changes);
yfit = b2(1) + b2(2).*eccentricity;
fig = figure()
hold on
scatter(eccentricity,totalpvs_changes,100,"x",'MarkerEdgeColor', 'w')
%plot([0 25], [0 25],'Color','r')
plot(eccentricity,yfit,'Color','w')


axis tight equal
ax = gca;              % Get current axes
ax.FontSize = 14;
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
fig.Color = [0 0 0];
