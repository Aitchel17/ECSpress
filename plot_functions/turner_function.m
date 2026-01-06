
%%











Fig4A = figure('Name','Figure Panel 4 - Turner et al. 2023','Units','Normalized','OuterPosition',[0,0,1,1]);
%% histogram of interblink interval
ax1 = subplot(2,3,1);
[~,edges] = histcounts(log10(data.allInterBlink));
histogram(data.allInterBlink,10.^edges,'Normalization','probability','EdgeColor',colors('vegas gold'),'FaceColor','k')
set(gca,'xscale','log')
xlabel('Interblink interval (IBI) (s)')
ylabel('Probability')
title('Interblink interval histogram')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% mean interblink interval
ax5 = subplot(2,3,2);
scatter(ones(1,length(data.interblink))*1,data.interblink,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('vegas gold'),'jitter','on','jitterAmount',0);
hold on
e1 = errorbar(1,data.meanInterblink,data.stdInterblink,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
ylabel('Mean IBI (s)')
title('Interblink interval')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,2])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% blinking post-stimulus
ax2 = subplot(2,3,3);
stimTimeVec = 0.5:0.5:5;
plot(stimTimeVec,data.meanBinProb,'color',colors('magenta'),'LineWidth',2)
hold on;
plot(stimTimeVec,data.meanBinProb + data.stdBinProb,'color',colors('magenta')','LineWidth',0.5)
plot(stimTimeVec,data.meanBinProb - data.stdBinProb,'color',colors('magenta')','LineWidth',0.5)
title('blink location post-stimulus')
xlabel('Post-stim time (s)')
ylabel('Probability')
set(gca,'box','off')
xlim([0.5,5]);
ylim([0,0.45])
axis square
ax2.TickLength = [0.03,0.03];
%% percentage blinks post-puff
ax6 = subplot(2,3,4);
scatter(ones(1,length(data.stimPerc))*1,data.stimPerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colors('magenta'),'jitter','on','jitterAmount',0);
hold on
e1 = errorbar(1,data.meanStimPerc,data.stdStimPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
ylabel('Percentage (%)')
title('Probabilty of blinking after a whisker puff')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,2])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% blink transitions
ax3 = subplot(2,3,5);
p1 = plot(awakeProbability,'color',colors('black'),'LineWidth',2);
hold on
p2 = plot(nremProbability,'color',colors('cyan'),'LineWidth',2);
p3 = plot(remProbability,'color',colors('candy apple red'),'LineWidth',2);
x1 = xline(7,'color',colors('magenta'),'LineWidth',2);
title('Peri-blink state probability')
xlabel('Peri-blink time (s)')
ylabel('Probability')
legend([p1,p2,p3,x1],'Awake','NREM','REM','Blink','Location','NorthWest')
xticks([1,3,5,7,9,11,13])
xticklabels({'-30','-20','-10','0','10','20','30'})
xlim([1,13])
ylim([0,1])
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% whisking before/after a blink
ax4 = subplot(2,3,6);
p1 = plot(timeVector,data.Awake.meanWhisk,'color',colors('black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanWhisk + data.Awake.stdWhisk,'color',colors('black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanWhisk - data.Awake.stdWhisk,'color',colors('black'),'LineWidth',0.5)
p2 = plot(timeVector,data.Asleep.meanWhisk,'color',colors('royal purple'),'LineWidth',2);
plot(timeVector,data.Asleep.meanWhisk + data.Asleep.stdWhisk,'color',colors('royal purple'),'LineWidth',0.5)
plot(timeVector,data.Asleep.meanWhisk - data.Asleep.stdWhisk,'color',colors('royal purple'),'LineWidth',0.5)
x1 = xline(0,'color',colors('magenta'),'LineWidth',2);
xlabel('Peri-blink Time (s)')
ylabel('Probability')
title('Peri-blink whisk probability')
legend([p1,p2,x1],'Awake','Asleep','Blink','Location','NorthEast')
set(gca,'box','off')
axis square
ax4.TickLength = [0.03,0.03];
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig4A,[dirpath 'Fig4A_JNeurosci2023']);
    set(Fig4A,'PaperPositionMode','auto');
    print('-vector','-dpdf','-bestfit',[dirpath 'Fig4A_JNeurosci2023'])
end
%% awake blink triggered averages
Fig4B = figure('Name','Figure Panel 4 - Turner et al. 2023','Units','Normalized','OuterPosition',[0,0,1,1]);
sgtitle('Awake high vs. low whisk blink triggered averages')
%% awake diameter
ax1 = subplot(2,3,1);
p1 = plot(timeVector,data.Awake.meanDiameter_highWhisk,'color',colors('black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanDiameter_highWhisk + data.Awake.stdDiameter_highWhisk,'color',colors('black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanDiameter_highWhisk - data.Awake.stdDiameter_highWhisk,'color',colors('black'),'LineWidth',0.5)
p2 = plot(timeVector,data.Awake.meanDiameter_lowWhisk,'color',colors('ash grey'),'LineWidth',2);
plot(timeVector,data.Awake.meanDiameter_lowWhisk + data.Awake.stdDiameter_lowWhisk,'color',colors('ash grey'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanDiameter_lowWhisk - data.Awake.stdDiameter_lowWhisk,'color',colors('ash grey'),'LineWidth',0.5)
xline(0,'color',colors('magenta'),'LineWidth',2);
ylabel('Diameter (z-units)')
xlabel('Peri-blink time (s)')
legend([p1,p2],'whisking blink','low whisk blink')
set(gca,'box','off')
xlim([-10,10])
axis square
ax1.TickLength = [0.03,0.03];
%% awake Hbt
ax2 = subplot(2,3,2);
plot(timeVector,data.Awake.meanHbT_highWhisk,'color',colors('black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanHbT_highWhisk + data.Awake.stdHbT_highWhisk,'color',colors('black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanHbT_highWhisk - data.Awake.stdHbT_highWhisk,'color',colors('black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanHbT_lowWhisk,'color',colors('battleship grey'),'LineWidth',2);
plot(timeVector,data.Awake.meanHbT_lowWhisk + data.Awake.stdHbT_lowWhisk,'color',colors('battleship grey'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanHbT_lowWhisk - data.Awake.stdHbT_lowWhisk,'color',colors('battleship grey'),'LineWidth',0.5)
xline(0,'color',colors('magenta'),'LineWidth',2);
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax2.TickLength = [0.03,0.03];
%% awake EMG
ax3 = subplot(2,3,3);
plot(timeVector,data.Awake.meanEMG_highWhisk,'color',colors('black'),'LineWidth',2);
hold on
plot(timeVector,data.Awake.meanEMG_highWhisk + data.Awake.stdEMG_highWhisk,'color',colors('black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanEMG_highWhisk - data.Awake.stdEMG_highWhisk,'color',colors('black'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanEMG_lowWhisk,'color',colors('ash grey'),'LineWidth',2);
plot(timeVector,data.Awake.meanEMG_lowWhisk + data.Awake.stdEMG_lowWhisk,'color',colors('ash grey'),'LineWidth',0.5)
plot(timeVector,data.Awake.meanEMG_lowWhisk - data.Awake.stdEMG_lowWhisk,'color',colors('ash grey'),'LineWidth',0.5)
xline(0,'color',colors('magenta'),'LineWidth',2);
ylabel('EMG Power (a.u.)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
axis square
ax3.TickLength = [0.03,0.03];
%% awake cortical low whisk
ax5 = subplot(2,4,5);
imagesc(T,F,data.Awake.meanCort_lowWhisk)
hold on
xline(0,'color',colors('magenta'),'LineWidth',2);
title('LW Cort LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-20,20])
axis square
axis xy
ax5.TickLength = [0.03,0.03];
%% awake cortical high whisk
ax6 = subplot(2,4,6);
imagesc(T,F,data.Awake.meanCort_highWhisk)
hold on
xline(0,'color',colors('magenta'),'LineWidth',2);
title('HW Cort LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-20,20])
axis square
axis xy
ax6.TickLength = [0.03,0.03];
%% awake hip low whisk
ax7 = subplot(2,4,7);
imagesc(T,F,data.Awake.meanHip_lowWhisk)
hold on
xline(0,'color',colors('magenta'),'LineWidth',2);
title('LW Hip LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-20,20])
axis square
axis xy
ax7.TickLength = [0.03,0.03];
%% awake hip high whisk
ax8 = subplot(2,4,8);
imagesc(T,F,data.Awake.meanHip_highWhisk)
hold on
xline(0,'color',colors('magenta'),'LineWidth',2);
title('HW Hip LFP')
ylabel('Freq (Hz)')
xlabel('Peri-blink time (s)')
set(gca,'box','off')
xlim([-10,10])
caxis([-20,20])
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
axis square
axis xy
ax8.TickLength = [0.03,0.03];
% axes properties
ax5Pos = get(ax5,'position');
ax8Pos = get(ax8,'position');
ax8Pos(3:4) = ax5Pos(3:4);
set(ax8,'position',ax8Pos);
%% save figure(s)
if saveFigs == true
    dirpath = [rootFolder delim 'MATLAB Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(Fig4B,[dirpath 'Fig4B_JNeurosci2023']);
    set(Fig4B,'PaperPositionMode','auto');
    print('-vector','-dpdf','-bestfit',[dirpath 'Fig4B_JNeurosci2023'])
end