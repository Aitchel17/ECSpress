%%
ssidx = 8;
clee = color_lee;
heatdata = secondary_struct(ssidx).heatdata;
spattern = 'flat';
hidx = 2;
fig = figure();
s = pcolor(heatdata(hidx).x_centers_aligned,heatdata(hidx).y_centers_aligned, log(heatdata(hidx).xy_counts_aligned));
s.FaceColor = spattern;
set(s,'EdgeColor','none');
set(gca, 'YDir', 'normal')
axis equal
ax = gca;              % Get current axes
ax.FontSize = 14;
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
fig.Color = [0 0 0];
colormap(clee.green_g)

%% pvs빈도 계산
%%
sum(heatdata(hidx).xy_counts,"all")
%%
tmp.bv_freq = sum(heatdata(hidx).xy_counts_aligned,1,"omitmissing");
%%
fig = figure()
b1 = bar(heatdata(hidx).x_centers_aligned,tmp.bv_freq);
b1.FaceColor = [1 0 0];   % 파랑
ax = gca;              % Get current axes
ax.FontSize = 14;
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
fig.Color = [0 0 0];
%%
tmp.bv_freq = sum(heatdata(hidx).xy_counts_aligned,2,"omitmissing");
%%
fig = figure()
b1 = bar(heatdata(hidx).y_centers_aligned,tmp.bv_freq);
b1.FaceColor = [0 0 1];   % 파랑
ax = gca;              % Get current axes
ax.FontSize = 14;
set(ax,'XColor',[1 1 1])
set(ax,'YColor',[1 1 1])
set(ax,'Color', [0 0 0])
fig.Color = [0 0 0];


%% UNFINISHED -- commented out so the file parses. see CLAUDE_LOG.md
% Overlaid BV and PVS marginals on one axis. Three names it needs do not exist
% anywhere in this file or in what feeds it:
%   tmp.baseline_bv    the b1 argument is literally empty below
%   tmp.baseline_pvs
%   primary_session
% The cells above it work and draw the heatmap and each marginal on its own.
%
% primary_session.infodict("mdfName")
%
% % 3.1 원점보정 히스토그램
% figure()
% hold on
%
% b1 = bar(heatdata.x_centers_aligned, tmp.baseline_bv, 'hist');
% b2 = bar(-heatdata.y_centers_aligned, tmp.baseline_pvs, 'hist');
%
% % 색상 설정
% b1.FaceColor = [0.2 0.6 1];   % 파랑
% b2.FaceColor = [1 0.4 0.4];   % 빨강
%
% % 투명도 설정
% b1.FaceAlpha = 0.5;
% b2.FaceAlpha = 0.5;
%
% % 공통 범위로 x축 설정
% all_x = [heatdata.x_baseceneters(:); heatdata.y_baseceneters(:)];
% xlim([min(all_x), max(all_x)]);
%
% xlabel('Offset from peak')
% ylabel('Counts')
% legend({'BV', 'PVS'})

