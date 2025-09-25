directories.load_dir = 'G:\tmp\01_igkltdt\hql080\250722_hql080_whiskerb\HQL080_whiskerb_250722_002';
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

%%

% ROI Level1 : 
% Primary axis (pax) - Where successful Kymograph generation to be guarented
% Extraparenchyma - For rough calculation of constricted, dilated state
try
roilist.addroi(roianalysis.preprocessed_ch1,'pax','line');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'pax');
end
%%
try
roilist.addroi(roianalysis.preprocessed_ch2,'extraparenchyma','polygon');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'extraparenchyma');
end
%% PAX BV analysis
pax_fwhm = line_fwhm(roilist.getvertices('pax'));
pax_fwhm.addkymograph("lumen", roianalysis.preprocessed_ch1)
pax_fwhm.kymograph_afterprocess('lumen',[5 2])
pax_fwhm.fwhm("lumen",0.5,10);


%% PAX BV summary fig
paxbv_fig = make_fig('paxBV_figure');
paxbv_fig.update_figsize([8 3])
paxbv_fig.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed)
%%
paxbv_fig.plot_line(medfilt1(pax_fwhm.idx.upperboundary,3),'r');
paxbv_fig.plot_line(medfilt1(pax_fwhm.idx.lowerboundary,3),'r');
%% PAX PVS analysis
pax_fwhm.addkymograph("pvs", roianalysis.preprocessed_ch2)
pax_fwhm.kymograph_afterprocess('lumen',[5 2])
pax_fwhm.addkymograph("pvs", roianalysis.preprocessed_ch2)
%%
paxpvs_fig = make_fig('paxPVS_figure');
paxpvs_fig.update_figsize([8 3])
paxpvs_fig.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed)
%%
plot_kymograph('kymograph',pax_fwhm.kymograph.kgph_lumen_processed,...
    'resolution', primary_datastruct.img_param.pixel2um,'fps',outfps_ch1,'blackbackground',true ...
    ,'idx_data',{pax_fwhm.idx.upperboundary,'r';...
    pax_fwhm.idx.lowerboundary,'g'});
%%
kgphfig = make_fig('kymograph');
kgphfig.plot_kymograph(pax_fwhm.kymograph.kgph_lumen_processed)
%%
kgphfig.reset_axis
kgphfig.plot_kymograph(x(7).kymograph)
%% load stack, vertices from main_primary script
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
x = analyze_polar(bv_stack,pax_center,pax_angle,true);

%%
try
roilist.addroi(preprocessed_ch2,'extraparenchyma','polygon');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'extraparenchyma');
end
%%
roianalysis.exp_bv = roilist.applyvertices(roianalysis.preprocessed_ch1,'extraparenchyma');
roianalysis.bv_psudoarea = analyze_quant_sigmablending(roianalysis.exp_bv);
%%
plot_bvarea = make_fig('bvarea');
plot_bvarea.update_figsize([8 4])
plot_bvarea.update_position([21 1.5])
plot_bvarea.convert_background(true);
plot_bvarea.bring_fig
plot_bvarea.reset_axis
plot_bvarea.plot_line(roianalysis.bv_psudoarea,'w');
plot_bvarea.update_figsize([15 4])
plot_bvarea.convert_background(true);
%%
constricted = find(roianalysis.bv_psudoarea < prctile(roianalysis.bv_psudoarea,10) & roianalysis.bv_psudoarea > prctile(roianalysis.bv_psudoarea,5));
dilated = find(roianalysis.bv_psudoarea > prctile(roianalysis.bv_psudoarea,90) & roianalysis.bv_psudoarea < prctile(roianalysis.bv_psudoarea,95));
%%

pax_fwhm.idx.lowerboundary - pax_fwhm.idx.upperboundary;

%%
figure()
plot(constricted)
%%
figure()
plot(dilated)
%%
figure()
plot(diff(constricted))
%% To

%% constricted : x axis frame count, y axis frame
%% Diff_ dim - 1, choose right side, difference should be 2 frame jump 
dilated = find(roianalysis.bv_psudoarea > prctile(roianalysis.bv_psudoarea,90) & roianalysis.bv_psudoarea < prctile(roianalysis.bv_psudoarea,95));
diff_dilated = diff(dilated);
diff_dilated = diff_dilated<3;
idxs_dilated = bwconncomp(diff_dilated);
length_dilated = cellfun(@numel,idxs_dilated.PixelIdxList);
[fifth_dilatedlegnth, fifth_dilatedidx] = maxk(length_dilated(2:end),5);
fith_idxs = idxs_dilated.PixelIdxList{fifth_dilatedidx(5)};
dilated_frames = dilated(fith_idxs);
%%

%%
[dilated_frames, dilatedlength] = analyze_prctileframecluster(roianalysis.bv_psudoarea ,[95 99]);
%%
[constricted_frames, constrictedlength] = analyze_prctileframecluster(roianalysis.bv_psudoarea ,[5 10]);
%%
img = mean(roianalysis.preprocessed_ch2(:,:,dilated_frames{5}),3);
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
roilist.save2disk
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
