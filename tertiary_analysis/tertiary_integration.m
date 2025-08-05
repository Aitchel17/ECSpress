clc, clear
folderpath = 'G:\tmp\**';
secondary_struct = secondary_integration(folderpath);

%%
window1=figure(Name='window1',NumberTitle='off');
%%
ssidx_list = [];
figure()
hold on
for ssidx = 1:size(secondary_struct,2)
    tmp.minpvschanges =  abs(min(secondary_struct(ssidx).heatdata(1).modespvs,[],'all'));
    tmp.maxpvschanges =  abs(max(secondary_struct(ssidx).heatdata(1).modespvs,[],'all'));

    if tmp.minpvschanges > tmp.maxpvschanges
    plot(secondary_struct(ssidx).heatdata(3).x_centers_aligned,secondary_struct(ssidx).heatdata(3).modespvs)
    ssidx_list = [ssidx_list, ssidx];
    else
        disp(ssidx)
    end
end

axis tight
%%
secondary_struct = secondary_struct(ssidx_list);
%%



%%
ssidx = 8;
demo.bvdata = secondary_struct(ssidx).heatdata(3).x_centers_aligned;
demo.pvsdata = secondary_struct(ssidx).heatdata(3).modespvs;
figure()
scatter(demo.bvdata,demo.pvsdata)
hold on
plot(demo.bvdata,demo.bvdata*(demo.bvdata'\demo.pvsdata))
%%
demo.upbvdata = secondary_struct(ssidx).heatdata(1).x_centers_aligned;
demo.uppvsdata = secondary_struct(ssidx).heatdata(1).modespvs;
figure()
scatter(demo.upbvdata,demo.uppvsdata)
hold on
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
%%
figure()
scatter(demo.upbvdata,demo.uppvsdata)
hold on
plot(demo.bvdata,demo.bvdata*(demo.bvdata'\demo.pvsdata))
%%
t = tiledlayout()
%%
figure()
plot(seco)




%%
figure()
%%
for k = 1:length(secondary_struct.goodsessionlist)
    disp(k)

hold on
plot(secondary_struct.(secondary_struct.goodsessionlist(k)).heatdata(1).x_centers_aligned,...
    secondary_struct.(secondary_struct.goodsessionlist(k)).heatdata(1).modespvs,'r')
plot(secondary_struct.(secondary_struct.goodsessionlist(k)).heatdata(2).x_centers_aligned,...
    secondary_struct.(secondary_struct.goodsessionlist(k)).heatdata(2).modespvs,'g')
pause(1)
end





%%
figure()
plot(tmp.result_up.x_baseceneters,tmp.result_up.modepvs,'r')
hold on
plot(tmp.result_down.x_baseceneters,tmp.result_down.modepvs,'r')
axis tight

%%

plot(heatdata.x_baseceneters,heatdata.modepvs)

%%
heatdata = secondary_struct.(secondary_struct.goodsessionlist(3)).heatdata;
spattern = 'flat';
hidx = 1
figure()
s = pcolor(heatdata(hidx).x_centers_aligned,heatdata(hidx).y_centers_aligned, heatdata(hidx).xy_counts_aligned);
s.FaceColor = spattern;
set(s,'EdgeColor','none');
set(gca, 'YDir', 'normal')
axis equal
axis tight
axis on
colormap(clee.gray_g)

%% pvs빈도 계산
 

figure()
bar(heatdata.y_centers,baseline_pvs,'hist')
xlim([0 max(heatdata.y_centers)])




%%
primary_session.infodict("mdfName")



%% plot function start from here




%


%%
window1=figure(Name='window1',NumberTitle='off');
window2=figure(Name='window2',NumberTitle='off');
%%

figure(window2)
hold off
plot(tmp.downbv,'r')
hold on
plot(tmp.downcell,'k')
plot(tmp.downpvs,'g')

plot(tmp.upbv,'r')
hold on
plot(tmp.upcell,'k')
plot(tmp.uppvs,'g')
figure(window1)
plot(tmp.downpvs-tmp.downbv)
plot(tmp.uppvs-tmp.upbv)
%%


% 3.1 원점보정 히스토그램
figure()
hold on

b1 = bar(heatdata.x_baseceneters, tmp.baseline_bv, 'hist');
b2 = bar(-heatdata.y_baseceneters, tmp.baseline_pvs, 'hist');

% 색상 설정
b1.FaceColor = [0.2 0.6 1];   % 파랑
b2.FaceColor = [1 0.4 0.4];   % 빨강

% 투명도 설정
b1.FaceAlpha = 0.5;
b2.FaceAlpha = 0.5;

% 공통 범위로 x축 설정
all_x = [heatdata.x_baseceneters(:); heatdata.y_baseceneters(:)];
xlim([min(all_x), max(all_x)]);

xlabel('Offset from peak')
ylabel('Counts')
legend({'BV', 'PVS'})

%%
test.snippet_fname = 'HQL080_whiskerb_250722_002';


