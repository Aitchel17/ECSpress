clc, clear
folderpath = 'G:\tmp\**';
secondary_struct = secondary_integration(folderpath);
%%
secondary_struct = secondary_afterproceesing(secondary_struct);
%%
slopes.dynamic = arrayfun(@(s) s.heatdata(1).slope, secondary_struct);
slopes.static = arrayfun(@(s) s.heatdata(2).slope, secondary_struct);
slopes.total = arrayfun(@(s) s.heatdata(3).slope, secondary_struct);
%%
plot_static_dynamic_slope = make_fig('static-dynamic_slope');
%%
plot_static_dynamic_slope.update_figsize([7 7])
%%
plot_static_dynamic_slope.plot_scatter(slopes.dynamic,slopes.static,'go')
%%
plot_static_dynamic_slope.change_xylim([-1 0],[-1 0])
%%
plot_static_dynamic_slope.initialize_axis
%%
plot_static_dynamic_slope.put_xaxistitle('Dynamic slope')
%%
plot_static_dynamic_slope.put_yaxistitle('Static slope')
%%
plot_static_dynamic_slope.fontsize = 20;

%%
plot_static_dynamic_slope.convert_background(true)
%%



axis(plot_static_dynamic_slope.ax,'equal')
%%
figure()
scatter(slopes.dynamic,slopes.static,'rx')
%%
hold on
plot(slopes.static,slopes.static)
%%
plot([-0.7, 0.1],[-0.7, 0.1])


%% check total modespvs value
figure()
clf
plot_specific_vd_pvst(filtered_secondarystruct,1)
%%
figure()
clf
plot_individual_vdpvst(filtered_secondarystruct)
%% Does the toal pvs thickness changes nonlinearly ( steep firsthalf gentle second half )
ssidx = 8;
idx = 3;

%%

firsthalf = arrayfun(@(s)  s.heatdata(2).slope_firsthalf, filtered_secondarystruct);
secondhalf = arrayfun(@(s)  s.heatdata(2).slope_secondhalf, filtered_secondarystruct);

%%
figure;
hold on
histogram(firsthalf, 'BinWidth', 0.1, 'FaceAlpha', 0.5);
histogram(secondhalf,'BinWidth', 0.1, 'FaceAlpha', 0.5);
legend('firsthalf','secondhalf');

%%
figure()
t = tiledlayout('flow')
idx= 3

for ssidx = 1:8
    nexttile
    tmp.secondarystruct = filtered_secondarystruct(ssidx).heatdata(idx);
scatter(tmp.secondarystruct.x_centers_aligned, tmp.secondarystruct.modespvs,"x",'MarkerEdgeColor', 'w')
hold on 
bv_length = size(tmp.secondarystruct.x_centers_aligned,2);
midpoint = floor(bv_length/2);
tmp.yfit = tmp.secondarystruct.slope.*tmp.secondarystruct.x_centers_aligned(1:end) +tmp.secondarystruct.intercept;
% plot(tmp.secondarystruct.x_centers_aligned(1:end),tmp.yfit(1:end),'LineWidth',2,color='k')
tmp.yfit = tmp.secondarystruct.slope_firsthalf.*tmp.secondarystruct.x_centers_aligned(1:midpoint) +tmp.secondarystruct.intercept_firsthalf;
plot(tmp.secondarystruct.x_centers_aligned(1:midpoint),tmp.yfit,'LineWidth',2,color='m')

tmp.yfit = tmp.secondarystruct.slope_secondhalf.*tmp.secondarystruct.x_centers_aligned(midpoint:end) +tmp.secondarystruct.intercept_secondhalf;
plot(tmp.secondarystruct.x_centers_aligned(midpoint:end),tmp.yfit,'LineWidth',2,color='c')
end

%% make figure background as black

fig = gcf;  % 현재 figure (또는 figure handle)
ax = findall(fig, 'type', 'axes');  % figure 내 모든 axes 핸들 찾기
set(ax, 'XLim', [-3 8], 'YLim', [-6 5]);  % 모든 축에 동일하게 적용
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
set(ax,'Color', [0 0 0])
% figure 안의 모든 axes 찾기
% 1) 우선 모두 axis equal 적용
arrayfun(@(a) axis(a, 'equal'), ax);
% 2) 모든 axes를 link해서 확대/이동이 같이 적용되도록
linkaxes(ax, 'xy');   % 'x', 'y', 'xy' 중 선택 가능
%%
secondary_struct.heatdata
%%
arrayfun(@(s) length(s.heatdata(3).x_centers_aligned), filtered_secondarystruct);

%%
mean(dynamic_slopes)
std(dynamic_slopes)
%%
d1 = dynamic_slopes./total_slopes;
d2 = static_slopes./total_slopes;
%%
d3 = d1+d2

%%
figure()
scatter(dynamic_slopes,static_slopes)
hold on
%% make Total pvs - Dynamic pvs plot
figure()
hold on
scatter(slopes.dynamic,slopes.total, 100, "x",MarkerEdgeColor='w',LineWidth=2)
plot([-1,0],[-1,0], color='r',LineWidth=2)
plot([-0.5,0],[-1,0],color='g',LineWidth=2)
axis tight equal
ax = gca;              % Get current axes
ax.FontSize = 18;
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
set(ax,'Color', [0 0 0])
xlim([-1.1 0])
ylim([-1.1 0])

%%
data = dynamic_slopes;
mu = mean(data)     % 평균
se = std(data) / sqrt(length(data));  % 표준오차 (standard error)
n = length(data);
alpha = 0.05;  % 95% CI
tval = tinv(1 - alpha/2, n - 1);  % t-score

ci95 = tval * se

%%
hist(slopes)
first_vals-second_vals
%%
figure()
plot_individual_vdpvst(secondary_struct)


figure('Position', [200, 200, 1200, 400])
plot(image_x,image_ch1,Color='r');
hold on;
plot(image_x,image_ch2,Color='g');

secondary_struct(ssidx).heatdata(3).modespvs


%%
figure()
plot(secondary_struct(1).heatdata(2).modespvs)
plot(demo.upbvdata,demo.upbvdata*(demo.upbvdata'\demo.uppvsdata))
%%

heatdata = secondary_struct(9).heatdata(3);
%%
tmp.baseline_pvs = sum(heatdata.xy_counts_aligned,2,'omitmissing');
[~,tmp.maxpvsloc]= max(tmp.baseline_pvs);
%%
heatdata.y_baseceneters = heatdata.y_centers_aligned - heatdata.y_centers_aligned(tmp.maxpvsloc);

%%
% Mode of PVS for each x point
[~,result.modepvslocs] = max(heatdata.xy_counts_aligned,[],1);
%%
heatdata.y_centers_aligned
%%
result.modepvs = zeros([length(result.modepvslocs),1]);
    for idx = 1:length(result.modepvslocs)
        result.modepvs(idx) = result.y_baseceneters(result.modepvslocs(idx));
    end
%%
demo.downbvdata = secondary_struct(1).heatdata(2).x_centers_aligned;
demo.downpvsdata = secondary_struct(1).heatdata(2).modespvs;
%%
demo.bvdata'\demo.pvsdata
%%
figure()
scatter(demo.bvdata,demo.pvsdata)
hold on
plot(demo.bvdata,demo.bvdata*(demo.bvdata'\demo.pvsdata))

%% Plot thickness changes of all example
figure()
t =tiledlayout("flow",TileSpacing="tight")
for idx = 1:8
nexttile
tmp.taxis = secondary_struct(idx).fwhmline.pax.t_axis(2:end);
plot(tmp.taxis,secondary_struct(idx).thickness(1).thickness)
hold on
plot(tmp.taxis,secondary_struct(idx).thickness(2).thickness)
plot(tmp.taxis,secondary_struct(idx).thickness(3).thickness)
plot(tmp.taxis,secondary_struct(idx).thickness(4).thickness)
hold off
axis tickaligned
axis tight
end
%%
plot_individual_thickness(secondary_struct,8)
%% plot one of the example

%%
airtable = secondary_struct(idx).analog.airtable;
%%
starttime = airtable.StartTime;
endtime = airtable.EndTime;
%%

secondary_struct(idx).analog.data.taxis(end)
%%

secondary_struct(idx).analog.data.raw_Ball
%%
figure()
plot(secondary_struct(idx).analog.data.taxis(end),secondary_struct(idx).analog.data.raw_Ball)

%%
figure()
plot(tmp.taxis,secondary_struct(idx).thickness(1).thickness*scale,'Color','k')

hold on
xline(starttime, 'Color', 'g')
xline(endtime,'Color','r')
plot(secondary_struct(idx).analog.data.raw_Ball,secondary_struct(idx).analog.data.taxis(end))


%%
plot(heatdata.x_baseceneters,heatdata.modepvs)



