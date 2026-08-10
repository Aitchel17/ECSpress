
% ROI class:
% roi related: roimod, vertices, refslices, roislice, stacks, mean, mask
% saving purpose:
% plot purpose: clist (color code)
% create roi: obj_name = roi(vivo_50um1.stackch2,'epecs','polygon',vivo_50um1);
% add roi:
% remove roi: roi = roi.removeroi('name');
% modify roi:

% apply mask: roi_applymask(vivo_50um1.stackch2,roi.mask.pax);

% load file
savepath = 'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\01_primary_analysis\01_pvscollapse';
%% load extracted mdf file
vivo_50um1 = mdfExtractLoader();
vivo_50um1.analog = vivo_50um1.loadanalog;
vivo_50um1.stackch1 = vivo_50um1.loadstack("ch1");
vivo_50um1.stackch2 = vivo_50um1.loadstack("ch2");


%% preprocessing
gpstack =  analyze_grouproject(vivo_50um1.stackch2,3,"mean");
medcsf = medfilt3(gpstack,[1 1 5]);
gausscsf= imgaussfilt(medcsf,1);
gpstack = analyze_grouproject(vivo_50um1.stackch1,3,"mean");
medbv = medfilt3(gpstack,[1 1 5]);
gaussbv= imgaussfilt(medbv,1);


%% making Roi object
roi1 = roi(gausscsf,'pax','line');
%%
roi1 = roi1.modifyroi(gausscsf,'pax');
%%

roi1 = roi1.addroi(gaussbv, 'aaa', 'rectangle');
%%

roi1.showroi()
%% dlz lsdlg


%%
roi1 = roi1.copyroi('epecs','cbv');
roi1 = roi1.modifyroi(vivo_50um1.stackch1,'cbv');
roi1 = roi1.addroi(vivo_50um1.stackch2,'pax','line');
roi1.mask.pvs = roi1.mask.epecs-roi1.mask.cbv;
%%
roi1 = roi1.removeroi('pax');
roi1 = roi1.addroi(vivo_50um1.stackch1,'pax','line');
%%
roi1 = roi1.modifyroi(vivo_50um1.stackch2,'pax');
% radon based method
% pvsstack = roi_applymask(gausscsf,roi1.mask.pax);
% [lprofile,sinogram] = analyze_fwhm_radon(pvsstack,roi1.vertices.pax(1:2,1:2),5);

%% affinebased method
roi1 = roi1.modifyroi(gausscsf,'pax');
%%
roi1 = roi1.removeroi('pax');
%%
roi1 = roi1.addroi(gaussbv,'pax','line');
%%
util_checkstack(gaussbv)
rc_bv=analyze_affine_rotate(gaussbv,roi1.vertices.pax(1:2,:), roi1.vertices.pax(3,1));
rc_csf=analyze_affine_rotate(gausscsf,roi1.vertices.pax(1:2,:), roi1.vertices.pax(3,1));
csf_projected = squeeze(sum(rc_csf,1));
bv_projected = squeeze(sum(rc_bv,1));
% normalize projection
a = bv_projected;
a = a-min(a,[],1);
a = a./max(a,[], 1);
%%
[bvpax_idx, bvpax_kymomask] = get_bvoutter(bv_projected);
%%
thresholded = false(size(bv_projected));
%%
bvpax_kymomask.boundary
bvpax_fwhm = bvpax_idx.lowerboundary_idx-bvpax_idx.upperboundary_idx;
%% reconstruction
v_thr = repmat(bvpax_kymomask.boundary,[1,1,size(rc_bv,1)]);
v_thr = permute(v_thr,[3,1,2]);
reconstructed = analyze_affine_reverse(v_thr,size(gaussbv),roi1.vertices.pax(1:2,:));
%%
util_checkstack(reconstructed)

%%
util_checkstack(~reconstructed.*gaussbv)
%%
[csfpax_idx, csf_paxmask] = get_pvsoutter(csf_projected,bvpax_idx.upperboundary_idx,bvpax_idx.lowerboundary_idx);
%%

figure()
imagesc(csf_paxmask.up)

%%
bvpax_idx.upperboundary_idx-csfidx.upboundary;

figure()
plot(ans)
hold on
plot(bvpax_fwhm)
%%
csfidx.downboundary-bvpax_idx.lowerboundary_idx ;

figure()
plot(ans)
hold on
plot(bvpax_fwhm)



%%
util_checkstack(rc_csf)
util_checkstack(rc_bv)


ch1 = mat2gray(csf_projected);
ch2 = mat2gray(bv_projected);
R = ch2;
G = ch1;
B = ch2;
mergedRGB = cat(3, R, G, B);

figure()
imshow(mergedRGB(:,500:1500,:))
%%

figure()
plot(csf_projected(:,2369)/max(csf_projected(:,2369),[],'all'))
hold on
plot(csf_projected(:,2510)/max(csf_projected(:,2510),[],'all'))

figure()

plot(bv_projected(:,2369)/max(bv_projected(:,2369),[],'all'))
hold on
plot(bv_projected(:,2510)/max(bv_projected(:,2510),[],'all'))

%%
roi1 = roi1.copyroi('pax','sax');
%%
roi1 = roi1.modifyroi(gausscsf,'sax');
%% stable csf
rc_bv2=analyze_affine_rotate(gaussbv,roi1.vertices.sax(1:2,:),roi1.vertices.sax(3,1));
rc_csf2=analyze_affine_rotate(gausscsf,roi1.vertices.sax(1:2,:),roi1.vertices.sax(3,1));
csf_projected2 = squeeze(sum(rc_csf2,1));
bv_projected2 = squeeze(sum(rc_bv2,1));
%%

%%
figure()
imagesc(sinogram)

%%
figure()
plot(bv_projected2(:,2380)/max(bv_projected2(:,2380),[],'all'))
hold on
plot(bv_projected2(:,2520)/max(bv_projected2(:,2520),[],'all'))
%%

%%
figure()
plot(a(:,2369))
hold on
plot(a(:,2520))
%
util_checkstack(rc_bv)

%%

% reconstruction
v_thr = repmat(thresholded,[1,1,size(rc_bv,1)]);
v_thr = permute(v_thr,[3,1,2]);


rotatedCrop = v_thr;
originalSize = size(gaussbv);
vertices = roi1.vertices.pax(1:2,:);

%%

reconstructed = analyze_affine_reverse(v_thr,size(gaussbv),roi1.vertices.pax(1:2,:));
%%
figure()
sliceViewer(gausscsf.*~reconstructed)

%%
figure()
sliceViewer(invRotated)
%%
roi1.modifyroi(gaussbv,'pax')
%%
c = csf_projected;
c = c-min(c,[],1);
c = c./max(c,[], 1);
figure()
plot(c(:,2369))
hold on
plot(c(:,2520))


%%
mask = mask.*a;



%%
% 0. seperate up and down while each does not contain center of vessel for brain-csf boundary detection by min()
% 1. find outter boundary of CSF (adjacent minimum value)
% 2. masking extrapaenchymal part
% 3. find maximum part

% 0.
%% top process
upcsf = false(sz);
upcsf(row_idx < upperboundary_idx) =1; % upper half
% 0. result
c_up = c;
c_up(~upcsf) = NaN;
%% 1. outter boundary detection
up_csfoffset = prctile(c_up,10,1); % bottom 10 percentile signal intensity
mask = false(sz); % initialize mask
mask(c_up<up_csfoffset) = 1; % brain to CSF max region, below offset become 1
[~,upoffsetloc] = max(mask.*row_idx,[],1); % bottom most offset location
%% FWHM calculation
upcsf = false(sz);
upcsf(row_idx < bv_centeridx) =1; % upper half
c_up = c;
c_up(~upcsf) = NaN;
c_up(row_idx<upoffsetloc) = NaN;
%
c_up = c_up-up_csfoffset;
c_up=c_up./max(c_up,[],1);
thr_cup = c_up >threshold;
csf_thridx = row_idx.*thr_cup;
csf_thridx(csf_thridx==0) =9999;
up_csf_boundary =min(csf_thridx,[],1);
%% bottom process
botcsf = false(sz);
botcsf(row_idx > bottomboundary_idx) = 1; % bottom half
c_bot = c;
c_bot(~botcsf) = NaN;
%%
bot_csfoffset = prctile(c_bot,10,1); % bottom 10 percentile signal intensity
mask = false(sz); % initialize mask
mask(c_bot<bot_csfoffset) = 1; % brain to CSF max region, below offset become 1
mask = mask.*row_idx;
mask(mask==0) = Inf;
[~,botoffsetloc] = min(mask.*row_idx,[],1); % bottom most offset location
%%
botcsf = false(sz);
downcsf(row_idx > bv_centeridx) =1; % bottom half
%%
c_down =c;
c_down(~downcsf) = NaN;
c_down(row_idx>botoffsetloc) = NaN;
%
c_down = c_down - bot_csfoffset;
c_down = c_down./max(c_down,[],1);
thr_cdown = c_down > threshold;
csf_thrdownidx = thr_cdown.*row_idx;
down_csf_Boundary = max(csf_thrdownidx,[],1);
%%
cb_thresholded = thresholded;
cb_thresholded(row_idx == down_csf_Boundary) = 1;
cb_thresholded(row_idx == up_csf_boundary) = 1;



%% PVS calculation require 1. vessel top center bottom boundary








%%
v_thr = repmat(cb_thresholded,[1,1,size(rc_bv,1)]);
v_thr = permute(v_thr,[3,1,2]);

rotatedCrop = v_thr;
originalSize = size(gaussbv);
vertices = roi1.vertices.pax(1:2,:);



reconstructed = analyze_affine_reverse(v_thr,size(gaussbv),roi1.vertices.pax(1:2,:));

%%
figure()
sliceViewer(gausscsf.*~reconstructed)
%%
util_checkstack(gausscsf.*~reconstructed)


%% bvside

bvsidecsf = false(sz);


%%
mask = false(sz); % initialize mask
mask(row_idx < upcsfmaxloc) = 1;  % brain -> CSF max position
up_csf = c; % initialize matrix
up_csf(~mask) = NaN; % masking brain to CSF max region
up_csfoffset = prctile(up_csf,10,1); % bottom 10 percentile signal intensity
%%
mask = false(sz); % initialize mask
mask(up_csf<up_csfoffset) = 1; % brain to CSF max region, below offset become 1
[~,upoffsetloc] = max(mask.*row_idx,[],1); % bottom most offset location
%%

up_thrcsf = c; % initialize value
up_thrcsf(row_idx<upoffsetloc) = NaN;
up_thrcsf(row_idx > upperboundary_idx) = NaN;
%%

figure()
imagesc(up_thrcsf)

%%
upcsfmax-up_csfoffset;
%%

mask.*row_idx;


%%


bv_csf(row_idx > botcsfmaxloc) = 1;
%%


%%

bv_csf = (c).*(~bv_csf);


%%
figure()
plot(bv_csf(:,3627))
hold on
plot(a(:,3627))

%%
% 1. maximum point =1
% 2. half max, offset = 10% offset value + 0.5

halfmax_value = prctile(c,10,1)+0.5;
[~,halfmax_pos]=min(abs(c-halfmax_value));

%%

figure()
plot(halfmax_pos);


%%

roi1 = roi1.modifyroi(vivo_50um1.stackch2,'epecs');



%%

roi1 = roi1.modifyroi(vivo_50um1.stackch2, 'pax');
%%
cmask = roi_applymask(vivo_50um1.stackch2,roi1.mask.epecs);
bmask = roi_applymask(vivo_50um1.stackch2,roi1.mask.cbv);
lmask = roi_applymask(vivo_50um1.stackch2,roi1.mask.pax);
%%
figure("Name",'a','NumberTitle','off')
imshow(roi1.mask.pax)
%%
roi1 = roi1.modifyroi(vivo_50um1.stackch1,'epecs');
%%
roi1 = roi1.copyroi('epecs','ipecs');
roi1 = roi1.modifyroi(vivo_50um1.stackch2,'ipecs');
%%
roi1 = roi1.copyroi('epecs','cbv');
roi1 = roi1.modifyroi(vivo_50um1.stackch1,'cbv');
%%

roi1.mean(parameter)
%%

roi1.removeroi('cellbody')


%%
roi1 = roi1.modifyroi(vivo_50um1.stackch1,'cbv');

%%
% put abstract image to roi object
roi1.stacks.medfilt = 3;
[roi1.stacks.ch1, roi1.stacks.ch1deprciate] = analyze_grouproject(vivo_50um1.stackch1,roi1.stacks.medfilt,'median');
[roi1.stacks.ch2, roi1.stacks.ch2deprciate] = analyze_grouproject(vivo_50um1.stackch2,roi1.stacks.medfilt,'median');
%% put stimulation, and image time axis
roi1.info.imgfpserror = 0.022; % fps correction
roi1.info.correctedfps = str2double(roi1.info.savefps)+roi1.info.imgfpserror; % corrected fps
roi1.mean.t = linspace(str2double(vivo_50um1.info.loadstart)/roi1.info.correctedfps,size(vivo_50um1.stackch1,3)/roi1.info.correctedfps,size(vivo_50um1.stackch1,3))'; % make time axis for image
%%
%%
% mask generation
[tmp.cols, tmp.rows] = meshgrid(1:size(vivo_50um1.stackch2, 2), 1:size(vivo_50um1.stackch2, 1));
roi1.mask.epecs = inpolygon(tmp.cols, tmp.rows,roi1.vertices.epecs(:, 1), roi1.vertices.epecs(:, 2));
roi1.mask.ipecs = inpolygon(tmp.cols, tmp.rows,roi1.vertices.ipecs(:, 1), roi1.vertices.ipecs(:, 2));
roi1.mask.cbv = double(inpolygon(tmp.cols, tmp.rows,roi1.vertices.cbv(:, 1), roi1.vertices.cbv(:, 2)));
roi1.mask.ip = double((roi1.mask.ipecs - roi1.mask.epecs)>0);
roi1.mask.pvs = double((roi1.mask.epecs - roi1.mask.cbv)>0);

% mask average calculation
% ip_mask: concentric ring, intraparenchymal space, space just inside of parenchyma
% pvs_mask: concentric ring, perivascular space, between the space where vessel constricted and border of parencyma
% epecs_mask: circle, extraparenchymal space, BV+Endothelial+Mural cell
% Intraparenchymal space
roi1.mask.ip((roi1.mask.ip == 0)) = NaN;
roi1.mask.pvs((roi1.mask.pvs == 0)) = NaN;
roi1.mask.cbv((roi1.mask.cbv == 0)) = NaN;
roi1.mean.ip_ch2 = squeeze(mean(double(vivo_50um1.stackch2).*roi1.mask.ip,[1,2],"omitmissing"));
roi1.mean.ip_ch1 = squeeze(mean(double(vivo_50um1.stackch1).*roi1.mask.ip,[1,2],"omitmissing"));
roi1.mean.pvs_ch2 = squeeze(mean(double(vivo_50um1.stackch2).*roi1.mask.pvs,[1,2],"omitmissing"));
roi1.mean.pvs_ch1 = squeeze(mean(double(vivo_50um1.stackch1).*roi1.mask.pvs,[1,2],"omitmissing"));
roi1.mean.cbv_ch2 = squeeze(mean(double(vivo_50um1.stackch2).*roi1.mask.cbv,[1,2],"omitmissing"));
roi1.mean.cbv_ch1 = squeeze(mean(double(vivo_50um1.stackch1).*roi1.mask.cbv,[1,2],"omitmissing"));
%%

figure('Name','a','NumberTitle','off')
plot(roi1.mean.t/60,medfilt1(roi1.mean.pvs_ch1/median(roi1.mean.pvs_ch1),10),Color=color.red)

hold on

plot(roi1.mean.t/60,medfilt1(roi1.mean.pvs_ch2/median(roi1.mean.pvs_ch2),10),Color=color.darkgreen)
plot(roi1.mean.t/60,medfilt1(roi1.mean.ip_ch2/median(roi1.mean.ip_ch2),10),Color=color.green)


%%

figure('Name','a','NumberTitle','off')
plot(medfilt1(roi1.mean.pvs_ch1/median(roi1.mean.pvs_ch1),10),Color=color.red)

hold on

plot(medfilt1(roi1.mean.ip_ch2/median(roi1.mean.ip_ch2),10),Color=color.darkgreen)




%% test underdura
roi1 = roi1(vivo_50um1.stackch1,'underduralvessel','polygon',vivo_50um1);
tmp.global_ch1 = squeeze(mean(vivo_50um1.stackch1,[1,2]))/median(squeeze(mean(vivo_50um1.stackch1,[1,2])));
% mask generation

[tmp.cols, tmp.rows] = meshgrid(1:size(vivo_50um1.stackch2, 2), 1:size(vivo_50um1.stackch2, 1));
roi1.mask.underdura = inpolygon(tmp.cols, tmp.rows,roi1.vertices.underduralvessel(:, 1), roi1.vertices.underduralvessel(:, 2));
%roi.mask.underdura((roi.mask.underdura == 0)) = NaN;
roi1.mean.underdura_ch1 = squeeze(mean(double(vivo_50um1.stackch1).*roi1.mask.underdura,[1,2],"omitmissing"));
roi1.mean.underdura_ch2 = squeeze(mean(double(vivo_50um1.stackch2).*roi1.mask.underdura,[1,2],"omitmissing"));
tmp.normch1_underdura = medfilt1(roi1.mean.underdura_ch1/mean(roi1.mean.underdura_ch1),10);
tmp.normch2_underdura = medfilt1(roi1.mean.underdura_ch2/mean(roi1.mean.underdura_ch2),10);
tmp.ratioIgkl_tdT_underdura = tmp.normch2_underdura./tmp.normch1_underdura;
%%
tmp.global_ch1 = squeeze(mean(vivo_50um1.stackch1,[1,2]))/median(squeeze(mean(vivo_50um1.stackch1,[1,2])));
%% test bv
roi1 = roi1.copyroi('underduralvessel','cbv');
roi1 = roi1.modifyroi(vivo_50um1.stackch1,'cbv');
%% mask generation
[tmp.cols, tmp.rows] = meshgrid(1:size(vivo_50um1.stackch2, 2), 1:size(vivo_50um1.stackch2, 1));
roi1.mask.cbv = inpolygon(tmp.cols, tmp.rows,roi1.vertices.cbv(:, 1), roi1.vertices.cbv(:, 2));
roi1.mean.cbv_ch1 = squeeze(mean(double(vivo_50um1.stackch1).*roi1.mask.cbv,[1,2],"omitmissing"));
roi1.mean.cbv_ch2 = squeeze(mean(double(vivo_50um1.stackch2).*roi1.mask.cbv,[1,2],"omitmissing"));
tmp.normch1_bv = medfilt1(roi1.mean.cbv_ch1/mean(roi1.mean.cbv_ch1),10);
tmp.normch2_bv = medfilt1(roi1.mean.cbv_ch2/mean(roi1.mean.cbv_ch2),10);
%%
tmp.ratioIgkl_tdT = tmp.normch2./tmp.normch1;
%%
figure()
plot(medfilt1(roi1.mean.cbv_ch1,3));

%%
figure()
plot(tmp.normch1_underdura,'Color',color.magenta);
hold on
plot(tmp.normch2_underdura,'Color',color.green);
plot(tmp.normch1_bv,'Color',color.black)
%%
tmp.ratioIgkl_tdT = tmp.normch2_underdura./tmp.normch1_underdura;
%%
figure()
plot(tmp.ratioIgkl_tdT,'color',color.black)
hold on
plot(tmp.normch1_bv,'Color',color.red)


%%

%%
svd_result = svd_filter([(roi1.mean.ip_ch1).',roi1.mean.ip_ch2.'],100);
%%
figure()
plot(((roi1.mean.ip_ch1./median(roi1.mean.ip_ch1)))-tmp.global_ch1)
hold on
plot(svd_result(:,2)*0.01/median(svd_result(:,2)))
%%
plot(((roi1.mean.ip_ch1./median(roi1.mean.ip_ch1)))-tmp.global_ch1)
hold on

plot(roi1.mean.ip_ch2./median(roi1.mean.ip_ch2)-tmp.global_ch1)


%%
% put stimulation, and image time axis
roi1.info.imgfpserror = 0.022; % fps correction
roi1.info.correctedfps = str2double(vivo_50um1.info.savefps)+roi1.info.imgfpserror; % corrected fps
roi1.mean.t = linspace(str2double(vivo_50um1.info.loadstart)/roi1.info.correctedfps,size(vivo_50um1.stackch1,3)/roi1.info.correctedfps,size(vivo_50um1.stackch1,3))'; % make time axis for image
%%
roi1.analog.data.air_puff_table = vivo_50um1.air_puff_extract('raw_Air_puff1'); % change the name to match analog air puff field name

tmp.air_durationlist = unique(round(roi1.analog.data.air_puff_table.Duration))';
tmp.airpuff_error_tol = 0.3;
% calculate triggered average
roi1.analog.data.airpuff_startframe = [];
for duration = tmp.air_durationlist
    tmp.fieldname = ['dur',num2str(duration)];
    disp(tmp.fieldname)
    tmp.airpuff_starttime = roi1.analog.data.air_puff_table(abs(roi1.analog.data.air_puff_table.Duration-duration) < tmp.airpuff_error_tol,:).StartTime;
    [~,roi1.analog.data.airpuff_startframe.(tmp.fieldname)] = min(abs(roi1.mean.t'-tmp.airpuff_starttime),[],2);
end
%% filtering successful trial
roi1.mean.cbv_ch1_airtrig15=analyze_trigaverage(roi1.mean.cbv_ch1,roi1.analog.data.airpuff_startframe.dur15,30,90);
roi1.mean.pvs_ch1_airtrig15=analyze_trigaverage(roi1.mean.pvs_ch1,roi1.analog.data.airpuff_startframe.dur15,30,90);
roi1.mean.filtered_idx = find(max(roi1.mean.pvs_ch1_airtrig15.list(30:end,:),[],1)>0.1);
%% filtering and recalculate
roi1.mean.pvs_ch1_airtrig15=analyze_trigaverage(roi1.mean.pvs_ch1,roi1.analog.data.airpuff_startframe.dur15(roi1.mean.filtered_idx),30,90);
roi1.mean.pvs_ch2_airtrig15=analyze_trigaverage(roi1.mean.pvs_ch2,roi1.analog.data.airpuff_startframe.dur15(roi1.mean.filtered_idx),30,90);
roi1.mean.ip_ch2_airtrig15=analyze_trigaverage(roi1.mean.ip_ch2,roi1.analog.data.airpuff_startframe.dur15(roi1.mean.filtered_idx),30,90);
roi1.mean.ip_ch1_airtrig15=analyze_trigaverage(roi1.mean.ip_ch1,roi1.analog.data.airpuff_startframe.dur15(roi1.mean.filtered_idx),30,90);
roi1.analog.data.airpuff_startframe.dur15(roi1.mean.filtered_idx)

%% Plot mean with shaded CI ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
t_array = linspace(-30/roi1.info.correctedfps,120/roi1.info.correctedfps,length(roi1.mean.pvs_ch2_airtrig15.mean));
figure()
hold on
analyze_shaded_graph(t_array,roi1.mean.ip_ch1_airtrig15.mean,roi1.mean.ip_ch1_airtrig15.ci,color.magenta)
analyze_shaded_graph(t_array,roi1.mean.ip_ch2_airtrig15.mean,roi1.mean.ip_ch1_airtrig15.ci,color.green)
analyze_shaded_graph(t_array,roi1.mean.pvs_ch1_airtrig15.mean,roi1.mean.ip_ch1_airtrig15.ci,color.red)
analyze_shaded_graph(t_array,roi1.mean.pvs_ch2_airtrig15.mean,roi1.mean.ip_ch2_airtrig15.ci,color.darkgreen)
xlabel('time (sec)')
ylabel('dF/F')
title('')
xlim([-1,8])
ylim([-0.10,0.3])



%% plot each stim
figure()
hold on
for i = 1:10
    plot(roi1.mean.pvs_ch1_airtrig15.list(:,i))
    pause(0.5)
end


%%
figure()
for i = 1:10
    hold off
    plot(roi1.mean.ip_ch2_airtrig15.list(:,i))
    hold on
    plot(roi1.mean.pvs_ch2_airtrig15.list(:,i))
    pause(2)
end
%%
roi1.analog.data.t = vivo_50um1.analog.data.t;
%% plot stim position

tmp.normalized_pvssignal = medfilt1((roi1.mean.pvs_ch1/median(roi1.mean.pvs_ch1,"all"))-1,20);
tmp.normalized_pvssignal2 = medfilt1((roi1.mean.pvs_ch2/median(roi1.mean.pvs_ch2,"all"))-1,20);
tmp.normalized_ipsignal2 = medfilt1((roi1.mean.ip_ch2/median(roi1.mean.ip_ch2,"all"))-1,20);
tmp.normalized_ipsignal1 =  medfilt1((roi1.mean.ip_ch1/median(roi1.mean.ip_ch1,"all"))-1,20);


tmp.normalized_raw_ball = roi1.analog.data.raw_Ball/max(roi1.analog.data.raw_Ball,[],'all');
tmp.thresholded_ball=zeros(size(tmp.normalized_raw_ball));
tmp.thresholded_ball(find(tmp.normalized_raw_ball>0.1))=0.1;

figure()
hold on
plot(roi1.mean.t, tmp.normalized_pvssignal,'color',color.red)
plot(roi1.mean.t, tmp.normalized_pvssignal2,'color',color.darkgreen)
plot(roi1.mean.t, tmp.normalized_ipsignal2,'color',color.green)
plot(roi1.mean.t, tmp.normalized_ipsignal1,'color',color.magenta)

plot(roi1.analog.data.t,tmp.thresholded_ball,'color',[1 1 1 0.2])
for idx = roi1.analog.data.airpuff_startframe.dur15
    t_start = roi1.mean.t(idx);
    t_end = t_start + 10;

    % Define x and y for the patch
    x_fill = [t_start, t_end, t_end, t_start];
    y_fill = [min(tmp.normalized_pvssignal), min(tmp.normalized_pvssignal), ...
        max(tmp.normalized_pvssignal), max(tmp.normalized_pvssignal)];

    % Draw the shaded area
    fill(x_fill, y_fill, [0.8 0.8 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none')

    % Optional: plot xline for reference
    % xline(t_start, 'k');
end
xlim([roi1.mean.t(1) roi1.mean.t(end)])
%ylim([0 1])
xlabel('Time (s)')
ylabel('dF/F')
title('PVS Signal with Air Puff Intervals')
%% for ppt
set(gcf, 'Color', 'k')
set(gca, 'XColor', 'w', 'YColor', 'w')        % Axis lines and ticks to white
set(gca, 'Color', 'k')
set(findall(gca, 'Type', 'line'), 'LineWidth',1);
set(gca, 'FontSize', 14);   % Set tick label font size (X, Y, Z axes)

%% Generate triggered average video
tmp.normstack = vivo_50um1.stackch1(:,:,idx-50:idx+499);
tmp.triggeredstacklist = zeros([size(tmp.normstack),length(roi1.mean.filtered_idx)]);
tmp.triggeredstacklist2 = zeros([size(tmp.normstack),length(roi1.mean.filtered_idx)]);

for i = 1:length(roi1.mean.filtered_idx)
    idx_list = roi1.analog.data.airpuff_startframe.dur10(roi1.mean.filtered_idx);
    idx = idx_list(i);
    tmp.baseimg = squeeze(median(double(vivo_50um1.stackch1(:,:,idx-50:idx)),3));
    tmp.normstack = double(vivo_50um1.stackch1(:,:,idx-50:idx+499))./tmp.baseimg;
    tmp.triggeredstacklist(:,:,:,i) = tmp.normstack;
    % ch2
    tmp.baseimg = squeeze(mean(double(vivo_50um1.stackch2(:,:,idx-50:idx)),3));
    tmp.normstack = double(vivo_50um1.stackch2(:,:,idx-50:idx+499))./tmp.baseimg;
    tmp.triggeredstacklist2(:,:,:,i) = tmp.normstack;
end
%%

mean(tmp.triggeredstacklist,4);
%%

analyze_savesimpletif((tmp.triggeredstacklist-1)*2048+1024,'ch1.tif')

%%

analyze_savesimpletif((tmp.triggeredstacklist2-1)*2048+1024,'ch2.tif')









%%

figure()
imshow(roi1.mask.ip)
%%


figure()
imshow(roi1.mask.cbv)

%% save
save(fullfile(roi1.info.analyzefolder,['pvs_collapse_',roi1.info.mdfName(1:end-4)]),"roi1")



% first make function that gets roi.vertices and
%%

roi1.showroi('cellbody1')




%%

% total imaging time
tmp.totalimg = str2num(vivo_50um1.info.fcount)*str2num(vivo_50um1.info.fduration(1:end-1));
% total analog time
tmp.totalanalog = str2num(vivo_50um1.analog.info.analogcount)/str2num(vivo_50um1.analog.info.analogfreq(1:end-2));

tmp.totalimg-tmp.totalanalog

%% IF CELL BODY DETECTED
roi1 = roi1.addroi(vivo_50um1.stackch1,'cellbody4','polygon');

%%
roi1 = roi1.modifyroi(vivo_50um1.stackch1,'cellbody1');
%%
roi1.mask.cellbody1 = double(inpolygon(tmp.cols, tmp.rows,roi1.vertices.cellbody1(:, 1), roi1.vertices.cellbody1(:, 2)));
roi1.mask.cellbody1((roi1.mask.cellbody1 == 0)) = NaN;
roi1.mean.cellbody1_ch2 = squeeze(mean(double(vivo_50um1.stackch2).*roi1.mask.cellbody1,[1,2],"omitmissing"));
roi1.mean.cellbody1_ch1 = squeeze(mean(double(vivo_50um1.stackch1).*roi1.mask.cellbody1,[1,2],"omitmissing"));
roi1.mean.cellbody1_ch1_airtrig15=analyze_trigaverage(roi1.mean.cellbody1_ch1,roi1.analog.data.airpuff_startframe.dur10(roi1.mean.filtered_idx));
roi1.mean.cellbody1_ch2_airtrig15=analyze_trigaverage(roi1.mean.cellbody1_ch2,roi1.analog.data.airpuff_startframe.dur10(roi1.mean.filtered_idx));

figure()
hold on
fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(roi1.mean.cellbody1_ch1_airtrig15.mean+roi1.mean.cellbody1_ch1_airtrig15.ci)' fliplr((roi1.mean.cellbody1_ch1_airtrig15.mean-roi1.mean.cellbody1_ch1_airtrig15.ci)')], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, roi1.mean.cellbody1_ch1_airtrig15.mean, 'r', 'LineWidth', 2)
fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(roi1.mean.cellbody1_ch2_airtrig15.mean+roi1.mean.cellbody1_ch2_airtrig15.ci)' fliplr((roi1.mean.cellbody1_ch2_airtrig15.mean-roi1.mean.cellbody1_ch2_airtrig15.ci)')], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, roi1.mean.cellbody1_ch2_airtrig15.mean, 'g', 'LineWidth', 2)
%% IF CELL BODY DETECTED, surrounding cellbody
roi1 = roi1.copyroi('cellbody1','aroundcellbody1');
%%
roi1 = roi1.modifyroi(vivo_50um1.stackch2,'aroundcellbody1');
%%
roi1.mask.aroundcellbody1 = double(inpolygon(tmp.cols, tmp.rows,roi1.vertices.aroundcellbody1(:, 1), roi1.vertices.aroundcellbody1(:, 2)));
roi1.mask.aroundcellbody1((roi1.mask.aroundcellbody1 == 0)) = NaN;
roi1.mean.aroundcellbody1_ch2 = squeeze(mean(double(vivo_50um1.stackch2).*roi1.mask.aroundcellbody1,[1,2],"omitmissing"));
roi1.mean.aroundcellbody1_ch1 = squeeze(mean(double(vivo_50um1.stackch1).*roi1.mask.aroundcellbody1,[1,2],"omitmissing"));
roi1.mean.aroundcellbody1_ch1_airtrig15=analyze_trigaverage(roi1.mean.aroundcellbody1_ch1,roi1.analog.data.airpuff_startframe.dur10(roi1.mean.filtered_idx));
roi1.mean.aroundcellbody1_ch2_airtrig15=analyze_trigaverage(roi1.mean.aroundcellbody1_ch2,roi1.analog.data.airpuff_startframe.dur10(roi1.mean.filtered_idx));

figure()
hold on
fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(roi1.mean.aroundcellbody1_ch1_airtrig15.mean+roi1.mean.aroundcellbody1_ch1_airtrig15.ci)' fliplr((roi1.mean.aroundcellbody1_ch1_airtrig15.mean-roi1.mean.aroundcellbody1_ch1_airtrig15.ci)')], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, roi1.mean.aroundcellbody1_ch1_airtrig15.mean, 'r', 'LineWidth', 2)
fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(roi1.mean.aroundcellbody1_ch2_airtrig15.mean+roi1.mean.aroundcellbody1_ch2_airtrig15.ci)' fliplr((roi1.mean.aroundcellbody1_ch2_airtrig15.mean-roi1.mean.aroundcellbody1_ch2_airtrig15.ci)')], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, roi1.mean.aroundcellbody1_ch2_airtrig15.mean, 'g', 'LineWidth', 2)
%%

%% IF CELL BODY2 DETECTED
roi1 = roi1.copyroi('cellbody1','cellbody2');
roi1 = roi1.modifyroi(vivo_50um1.stackch1,'cellbody2');
%%
roi1.mask.cellbody2 = double(inpolygon(tmp.cols, tmp.rows,roi1.vertices.cellbody1(:, 1), roi1.vertices.cellbody1(:, 2)));
roi1.mask.cellbody2((roi1.mask.cellbody2 == 0)) = NaN;
roi1.mean.cellbody2_ch2 = squeeze(mean(double(vivo_50um1.stackch2).*roi1.mask.cellbody2,[1,2],"omitmissing"));
roi1.mean.cellbody2_ch1 = squeeze(mean(double(vivo_50um1.stackch1).*roi1.mask.cellbody2,[1,2],"omitmissing"));
roi1.mean.cellbody2_ch1_airtrig15=analyze_trigaverage(roi1.mean.cellbody2_ch1,roi1.analog.data.airpuff_startframe.dur10(roi1.mean.filtered_idx));
roi1.mean.cellbody2_ch2_airtrig15=analyze_trigaverage(roi1.mean.cellbody2_ch2,roi1.analog.data.airpuff_startframe.dur10(roi1.mean.filtered_idx));

figure()
hold on
fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(roi1.mean.cellbody2_ch1_airtrig15.mean+roi1.mean.cellbody2_ch1_airtrig15.ci)' fliplr((roi1.mean.cellbody2_ch1_airtrig15.mean-roi1.mean.cellbody2_ch1_airtrig15.ci)')], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, roi1.mean.cellbody2_ch1_airtrig15.mean, 'r', 'LineWidth', 2)
fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(roi1.mean.cellbody2_ch2_airtrig15.mean+roi1.mean.cellbody2_ch2_airtrig15.ci)' fliplr((roi1.mean.cellbody2_ch2_airtrig15.mean-roi1.mean.cellbody2_ch2_airtrig15.ci)')], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, roi1.mean.cellbody2_ch2_airtrig15.mean, 'g', 'LineWidth', 2)
%% IF CELL BODY DETECTED, surrounding cellbody
roi1 = roi1.copyroi('cellbody2','aroundcellbody2');
%%
roi1 = roi1.modifyroi(vivo_50um1.stackch1,'aroundcellbody2');
%%
roi1.mask.aroundcellbody2 = double(inpolygon(tmp.cols, tmp.rows,roi1.vertices.aroundcellbody2(:, 1), roi1.vertices.aroundcellbody2(:, 2)));
roi1.mask.aroundcellbody2((roi1.mask.aroundcellbody2 == 0)) = NaN;
roi1.mean.aroundcellbody2_ch2 = squeeze(mean(double(vivo_50um1.stackch2).*roi1.mask.aroundcellbody2,[1,2],"omitmissing"));
roi1.mean.aroundcellbody2_ch1 = squeeze(mean(double(vivo_50um1.stackch1).*roi1.mask.aroundcellbody2,[1,2],"omitmissing"));
roi1.mean.aroundcellbody2_ch1_airtrig15=analyze_trigaverage(roi1.mean.aroundcellbody2_ch1,roi1.analog.data.airpuff_startframe.dur10(roi1.mean.filtered_idx));
roi1.mean.aroundcellbody2_ch2_airtrig15=analyze_trigaverage(roi1.mean.aroundcellbody2_ch2,roi1.analog.data.airpuff_startframe.dur10(roi1.mean.filtered_idx));

figure()
hold on
fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(roi1.mean.aroundcellbody2_ch1_airtrig15.mean+roi1.mean.aroundcellbody2_ch1_airtrig15.ci)' fliplr((roi1.mean.aroundcellbody2_ch1_airtrig15.mean-roi1.mean.aroundcellbody2_ch1_airtrig15.ci)')], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, roi1.mean.aroundcellbody2_ch1_airtrig15.mean, 'r', 'LineWidth', 2)
fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(roi1.mean.aroundcellbody2_ch2_airtrig15.mean+roi1.mean.aroundcellbody2_ch2_airtrig15.ci)' fliplr((roi1.mean.aroundcellbody2_ch2_airtrig15.mean-roi1.mean.aroundcellbody2_ch2_airtrig15.ci)')], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, roi1.mean.aroundcellbody2_ch2_airtrig15.mean, 'g', 'LineWidth', 2)



%% Correction
figure()
hold on
for i = 1:size(trig_data.list,2)
    plot(trig_data.list(:,i))
end
%%
figure()
hold on
for i = 1:size(trig_data.list,2)
    plot(trig.data.list(:,i))
end


%%
figure()
hold on
trig_data.list = roi1.mean.aroundcellbody1_ch1_airtrig15.list-roi1.mean.cellbody1_ch1_airtrig15.list;
trig_data.mean = mean(trig_data.list, 2);
tmp.sem_trace = std(trig_data.list, 0, 2) / sqrt(size(trig_data.list, 2));
tmp.tval = tinv(0.975, size(trig_data.list, 2) - 1);
trig_data.ci = tmp.tval * tmp.sem_trace;
fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(trig_data.mean+trig_data.ci)' fliplr((trig_data.mean-trig_data.ci)')], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, trig_data.mean, 'r', 'LineWidth', 2)
%%
trig.data.list = roi1.mean.aroundcellbody1_ch2_airtrig15.list-roi1.mean.cellbody1_ch2_airtrig15.list;
trig.data.mean = mean(trig.data.list, 2);
tmp.sem_trace = std(trig.data.list, 0, 2) / sqrt(size(trig.data.list, 2));
tmp.tval = tinv(0.975, size(trig.data.list, 2) - 1);
trig.data.ci = tmp.tval * tmp.sem_trace;

fill([roi1.mean.airtrig_t fliplr(roi1.mean.airtrig_t)], [(trig.data.mean+trig.data.ci)' fliplr((trig.data.mean-trig.data.ci)')], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.4) % closed polygon generation l->r, r->l
plot(roi1.mean.airtrig_t, trig.data.mean, 'g', 'LineWidth', 2)
%%
figure()
hold on
plot(roi1.mean.aroundcellbody1_ch2_airtrig15.mean)
plot(roi1.mean.cellbody1_ch2_airtrig15.mean)


%%
plot(roi1.mean.aroundcellbody1_ch1_airtrig15.mean-roi1.mean.cellbody1_ch1_airtrig15.mean);
%%
plot(roi1.mean.aroundcellbody1_ch2_airtrig15.mean-roi1.mean.cellbody1_ch2_airtrig15.mean);


%% show mask
roi1.mask.colormask = zeros([size(roi1.mask.cbv),3],"uint8");
%
roi1.mask.colormask(:,:,1) = uint8(roi1.mask.cbv)*255+uint8(roi1.mask.cellbody2)*255;
roi1.mask.colormask(:,:,2) = uint8(roi1.mask.pvs)*255+uint8(roi1.mask.ip)*255;
roi1.mask.colormask(:,:,3) = uint8(roi1.mask.ip)*255+uint8(roi1.mask.cellbody2)*255;

figure()
imshow(roi1.mask.colormask)

%%
roi1.showroi('epecs')


%%
figure()
plot()
hold on
%%
plot(mean(airtrig5sec.csf,2))
plot(mean(airtrig5sec.isf,2))

%%
pvs.tdT_filt = medfilt1(pvs.tdT_mean,100);
pvs.tdT_diff = diff(pvs.tdT_filt);
pvs.tdT_diff = movmean(pvs.tdT_diff(2:end),10);
%%
pvs.tdT_diff = pvs.tdT_diff-median(pvs.tdT_diff);
pvs.tdT_diff = pvs.tdT_diff/median(pvs.tdT_diff);


%%
figure()
plot(roi1.mean.t(3:end),pvs.tdT_diff)
hold on
plot(roi1.mean.t(2:end),(pvs.tdT_filt(2:end)/median(pvs.tdT_filt(2:end))-1))

hold on
plot(analog_x,air_puff/400000)


%%
figure()
plot(roi1.mean.t(2:end),(pvs.tdT_filt(2:end)/median(pvs.tdT_filt(2:end))-1))
hold on
plot(roi1.mean.t(3:end),pvs.tdT_diff(2:end)/max(pvs.tdT_diff(1:end)))
%%
hold on
plot(roi1.mean.t(2:end),(pvs.igkl_mean(2:end)/median(pvs.igkl_mean(1:end))-1))

%%
figure()
plot(roi1.mean.t(2:end),movmean(pvs.tdT_diff,10))
%%
hold on
plot(analog_x,air_puff)
%%
x = pvs.tdT_diff(2:end);
y = pvs.igkl_mean(4:end);
%%
x-mean(x)
%%
figure()
[r,lags] = xcorr(x-mean(x),y-mean(y));
plot(lags,r)
%%
figure()
scatter(x(1:end-50),y(51:end))
%%
ipecs.modifyroi(vivo_50um1.stackch2);
%%
analyze_norm_medmax(pvs.tdT_mean,'pvs TRITC',ip.igkl_mean,'ip IgKL',roi1.mean.t, 5000,analog_x,air_puff)
analyze_norm_medmax(pvs.tdT_mean,'pvs TRITC',pvs.igkl_mean,'pvs IgKL',roi1.mean.t, 5000,analog_x,air_puff)
analyze_norm_medmax(ip.tdT_mean,'ip tdT',ip.igkl_mean,'ip EGFP',roi1.mean.t, 5000,analog_x,air_puff)
analyze_norm_medmax(pvs.tdT_mean,'pvs tdT',pvs.igkl_mean,'pvs IgKL',roi1.mean.t, 5000,analog_x,air_puff)
analyze_norm_medmax(pvs.tdT_mean,'pvs tdT',ip.tdT_mean,'ip tdT',roi1.mean.t, 5000,analog_x,air_puff)
analyze_norm_medmax(pvs.tdT_mean,'pvs tdT',pvs.igkl_mean,'pvs IgKL',roi1.mean.t, 5000,analog_x,air_puff)


%%
%% conversion from old roi struct
disp('convert')
roi1.roimod.cbv = cbv.roimod;
roi1.vertices.cbv = cbv.vertices;
roi1.refslice.cbv = cbv.refslice;
roi1.stacks.cbv = cbv.stacks;
roi1.mask.cbv = cbv.mask;
roi1.info    = cbv.info;
roi1.analog = cbv.analog;

roi1.roimod.epecs = epecs.roimod;
roi1.vertices.epecs = epecs.vertices;
roi1.refslice.epecs = epecs.refslice;
roi1.stacks.epecs = epecs.stacks;
roi1.mask.epecs = epecs.mask;
roi1.info    = epecs.info;
roi1.analog = epecs.analog;

roi1.roimod.ipecs = ipecs.roimod;
roi1.vertices.ipecs = ipecs.vertices;
roi1.refslice.ipecs = ipecs.refslice;
roi1.stacks.ipecs = ipecs.stacks;
roi1.mask.ipecs = ipecs.mask;
roi1.info    = ipecs.info;
roi1.analog = ipecs.analog;
