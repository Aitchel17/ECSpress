directories.load_dir = 'G:\tmp\01_igkltdt\hql080\250722_hql080_whiskerb\HQL080_whiskerb_250722_007';
directories.save_dir = fullfile(directories.load_dir,"primary_analysis");
directories.savefig_dir =  'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\08_figures';

figure.clee = color_lee;
primary_datastruct = primary_loaddata(directories.load_dir);
directories.savepath = fullfile(primary_datastruct.mdfextract.info.analyzefolder,'primary_analysis');
% fix frame rate to 5 fps
[roianalysis.preprocessed_ch2, roianalysis.outfps_ch2] = analyze_resample(primary_datastruct.stackch2,primary_datastruct.img_param.fps,5);
[roianalysis.preprocessed_ch1, roianalysis.outfps_ch1] = analyze_resample(primary_datastruct.stackch1,primary_datastruct.img_param.fps,5);
% load roi.mat if exist
roilist = roi_handle(directories.save_dir);

%% ROI Level 1.1, PAX : 
% Primary axis (pax) - Where successful Kymograph generation to be guarented
try
roilist.addroi(roianalysis.preprocessed_ch1,'pax','line');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch2,'pax');
end
%% ROI Level 1.2, Extraparenchyma : 
% Extraparenchyma - For rough calculation of constricted, dilated state
try
roilist.addroi(roianalysis.preprocessed_ch2,'extraparenchyma','polygon');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'extraparenchyma');
end
%% ROI Level 1.3, Dip in angle : 
% Extraparenchyma - For rough calculation of constricted, dilated state
try
roilist.addroi(roianalysis.preprocessed_ch1,'dipin','line');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'dipin');
end
%% ROI Level 1.4, Dip out angle : 
% Extraparenchyma - For rough calculation of constricted, dilated state
try
roilist.addroi(roianalysis.preprocessed_ch2,'dipout','line');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'dipout');
end

%% Put image slice for Level 1 ROI list
roilist.addimgchannel(cat(4,roianalysis.preprocessed_ch1,roianalysis.preprocessed_ch2),...
    ["pax","extraparenchyma","dipin","dipout"]); % put image slices

%%


%% Rough calculation of BV area using multi sigma blending after thresholding 250920-250924
roianalysis.exp_bv = roilist.applyvertices(roianalysis.preprocessed_ch1,'extraparenchyma');
[roianalysis.bv_psudoarea,roianalysis.bv_stack] = analyze_quant_sigmablending(roianalysis.exp_bv);
%%

util_checkstack(roianalysis.bv_stack)
%% Plot area
plot_bvarea = make_fig('bvarea');
plot_bvarea.update_figsize([8 4])
plot_bvarea.update_position([21 1.5])
plot_bvarea.convert_background(true);
plot_bvarea.bring_fig
plot_bvarea.reset_axis
plot_bvarea.plot_line(roianalysis.bv_psudoarea,'w');
plot_bvarea.update_figsize([15 4])
plot_bvarea.convert_background(true);
%% Find constricted and dilated frame chunk
[roianalysis.dilated_framechunks, roianalysis.dilatedlength] = analyze_prctileframecluster(roianalysis.bv_psudoarea ,[95 99]);
[roianalysis.constricted_framechunks, roianalysis.constrictedlength] = analyze_prctileframecluster(roianalysis.bv_psudoarea ,[5 10]);
%%
roianalysis.dilated_ch1 = roianalysis.preprocessed_ch1(:,:,cat(1,roianalysis.dilated_framechunks{:}));
roianalysis.dilated_ch2 = roianalysis.preprocessed_ch2(:,:,cat(1,roianalysis.dilated_framechunks{:}));
roianalysis.constricted_ch1 = roianalysis.preprocessed_ch1(:,:,cat(1,roianalysis.constricted_framechunks{:}));
roianalysis.constricted_ch2 = roianalysis.preprocessed_ch2(:,:,cat(1,roianalysis.constricted_framechunks{:}));
% Using clustered frames to make constricted and dilated ROI 250925
%% ROI level 2 
try
    roilist.copyroi('extraparenchyma','constricted_bv');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'constricted_bv');
end
%%
try
    roilist.copyroi('extraparenchyma','constricted_bv');
catch ME
    disp(ME.message)
    roilist.modifyroi(cluster3_bv,'constricted_bv');
end

%%
try
    roilist.copyroi('extraparenchyma','dilated_bv');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch2(:,:,kgph_lumen_columnidx),'dilated_bv');
end
%%
try
    roilist.copyroi('extraparenchyma','dilated_bv');
catch ME
    disp(ME.message)
    roilist.modifyroi(cluster67_bv,'dilated_bv');
end

%%
try
    roilist.copyroi('extraparenchyma','constricted_bv');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1(:,:,cat(1,roianalysis.constricted_framechunks{:})),'constricted_bv');
end
%%
try
    roilist.copyroi('extraparenchyma','dilated_bv');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1(:,:,cat(1,roianalysis.dilated_framechunks{:})),'dilated_bv');
end
%% 
roilist.addimgchannel(cat(4,roianalysis.preprocessed_ch1,roianalysis.preprocessed_ch2),...
    ["dilated_bv","constricted_bv"]); % put image slices

%%
roilist.addimgchannel(cat(4,roianalysis.dilated_ch1, roianalysis.dilated_ch2),"dilated_bv"); % put image slices
roilist.addimgchannel(cat(4,roianalysis.constricted_ch1, roianalysis.constricted_ch2),"constricted_bv"); % put image slices

%%
roilist.save2disk
%%
roi_fig = make_fig('roi_fig');
roi_fig.update_figsize([8 6])

%%
roi_fig.showroi(roilist,'dilated_bv',1)
%%
roi_fig.reset_axis

roi_fig.showrois(roilist,2,["dilated_bv","constricted_bv"])
%%
roi_fig.reset_axis

roi_fig.showrois(roilist,2,["constricted_bv","dilated_bv"])

%%
roi_fig.showrois(roilist,1,["dipin","dipout"])
%%

%%
%% 2. Do polar analysis and get sinogram like things
img = mean(roianalysis.dilated_ch1,3);
figure()
imshow(mat2gray(img))
%%
img = mean(roianalysis.preprocessed_ch2(:,:,constricted_frames{5}),3);
figure()
imshow(mat2gray(img))

%%
img = median(roianalysis.preprocessed_ch1(:,:,constricted_frames{2}),3);
figure()
imshow(mat2gray(img))

%% sinogram generation






%%
plot_bvpaxdiam = make_fig('bv_paxdiam');
plot_bvpaxdiam.update_figsize([8 4])
plot_bvpaxdiam.update_position([21 1.5])
plot_bvpaxdiam.convert_background(true);
plot_bvpaxdiam.bring_fig
plot_bvpaxdiam.reset_axis
plot_bvpaxdiam.plot_line(ans,'w');
plot_bvpaxdiam.update_figsize([15 4])
plot_bvpaxdiam.convert_background(true);
constricted = find(bv_psudoarea < prctile(bv_psudoarea,10) & bv_psudoarea > prctile(bv_psudoarea,5));
dilated = find(bv_psudoarea > prctile(bv_psudoarea,90) & bv_psudoarea < prctile(bv_psudoarea,95));
%%
plot_frame_bv = make_fig('frames_bv');
%%
plot_frame_bv.bring_fig
%%
plot_frame_bv.plot_line(constricted,'g');
%%
plot_frame_bv.hold_axis(true);
%%
plot_frame_bv.plot_line(dilated,'r')

%%
plot_frame_bv.update_position([21 1.5])

plot_frame_bv.update_figsize([15 4])





%%
try
roilist.copyroi('extraparenchyma','constricted');
catch ME
    disp(ME.message)
    roilist.modifyroi(preprocessed_ch1(:,:,constricted),'constricted');
end
%%
analyze_savesimpletif(preprocessed_ch1(:,:,constricted),fullfile(directories.savepath,"constricted_tif"))
%%
analyze_savesimpletif(preprocessed_ch2(:,:,constricted),fullfile(directories.savepath,"constricted_csftif"))


%%
%%
function [] = Semilog_ImageSC(x,y,C,logaxis)
% 9/2018 Patrick Drew
% make a surface at points x,y, of height 0 and with colors given by the matrix C
% logaxis - which axis to plot logarithmically: 'x', 'y' or 'xy'
surface(x,y,zeros(size(C)),(C),'LineStyle','none');
q = gca;
q.Layer = 'top'; % put the axes/ticks on the top layer
if strcmp(logaxis,'y') == 1
    set(gca,'YScale','log');
elseif strcmp(logaxis,'x') == 1
    set(gca,'XScale','log');
elseif strcmp(logaxis,'xy') == 1
    set(gca,'XScale','log');
    set(gca,'YScale','log');
end
axis xy
axis tight
end
