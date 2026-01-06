
% Photobleach experiment analysis script



    % extra-mNG + intra-tdT
    % intra-mNG + intra-tdT (control)

% 1. Fixed sample
    % 1.1 Brain vivo (proof of concept)
        % 1.1.1 Neuropile
        % 1.1.2 Cellbody
        % 1.1.3 vessel lumen (void)
            
    % 1.2 Ex-vivo, perfusion fix + decapitate (fixed control)
        % 1.1.1 Neuropile
        % 1.1.2 Cellbody
        % 1.1.3 pvs (void)


% 2. In vivo alive
        % 2.1.1
        % 2.1.2
        
% load file
savepath = 'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\01_primary_analysis\01_pvscollapse';
%%
vivo_50um1 = mdfExtractLoader();
vivo_50um1.analog = vivo_50um1.loadanalog;
vivo_50um1.stackch1 = vivo_50um1.loadstack("ch1");
vivo_50um1.stackch2 = vivo_50um1.loadstack("ch2");


%%

%epecs : extraparenchymal extracellularspace
%ipecs : intraparenchymal extracellularspace
%cbv : when vesel constricted

%%
epecs = dilation(vivo_50um1.stackch2,"polygon",vivo_50um1);
%%
epecs = epecs.modifyroi(vivo_50um1.stackch1);
%%
ipecs = epecs;
%%
ipecs = ipecs.modifyroi(vivo_50um1.stackch2);
%%
cbv = epecs;
cbv = cbv.modifyroi(vivo_50um1.stackch2);

%%
tmp.filename = vivo_50um1.info.mdfName(1:end-4);

%%
save(fullfile(savepath,tmp.filename),'epecs','cbv','ipecs')



%%
figure()
imshow(ip_mask+cbv_mask)

%%
epecs = epecs.modifyroi(vivo_50um1.stackch2);


%% Mask generation
stack = vivo_50um1.stackch2;
[cols, rows] = meshgrid(1:size(stack, 2), 1:size(stack, 1));
cbv_mask = inpolygon(cols, rows,cbv.vertices(:, 1), cbv.vertices(:, 2));
epecs_mask = inpolygon(cols, rows,epecs.vertices(:, 1), epecs.vertices(:, 2));
ipecs_mask = inpolygon(cols, rows,ipecs.vertices(:, 1), ipecs.vertices(:, 2));

ip_mask = ipecs_mask - epecs_mask;
pvs_mask = epecs_mask - cbv_mask;
%
% ip_mask: concentric ring, intraparenchymal space, space just inside of parenchyma
% pvs_mask: concentric ring, perivascular space, between the space where vessel constricted and border of parencyma
% epecs_mask: circle, extraparenchymal space, BV+Endothelial+Mural cell

% Intraparenchymal space
ip_mask((ip_mask == 0)) = nan;

ip_igkl = double(vivo_50um1.stackch2).*ip_mask;
ip_igkl_mean = squeeze(mean(ip_igkl,[1,2],"omitmissing"));


ip_tdT = double(vivo_50um1.stackch1).*ip_mask;
ip_tdT_mean = squeeze(mean(ip_tdT,[1,2],"omitmissing"));


pvs_mask((pvs_mask == 0)) = nan;

pvs_igkl = double(vivo_50um1.stackch2).*pvs_mask;
pvs_igkl_mean = squeeze(mean(pvs_igkl,[1,2],"omitmissing"));

pvs_tdT = double(vivo_50um1.stackch1).*pvs_mask;
pvs_tdT_mean = squeeze(mean(pvs_tdT,[1,2],"omitmissing"));

cbv_igkl = double(vivo_50um1.stackch2).*cbv_mask;
cbv_igkl_mean = squeeze(mean(cbv_igkl,[1,2],"omitmissing"));

cbv_tdT = double(vivo_50um1.stackch1).*cbv_mask;
cbv_tdT_mean = squeeze(mean(cbv_tdT,[1,2],"omitmissing"));
%%






%%
y = ip_igkl_mean(5000:10000);
x = ip_tdT_mean(5000:10000);
x = medfilt1(x,10);
y = medfilt1(y,10);

[b, stats] = robustfit(x, y);
y_fit = b(1) + b(2)*x;
rsquare = corr(y,b(1)+b(2)*x)^2;


figure()
plot(x);
hold on;

plot(y);

figure('Position', [200, 200, 600, 800])
scatter(x,y)
hold on;
% Overlay robust fit line
plot(x, y_fit, 'r-', 'LineWidth', 2);
% Optional: add labels and legend
xlabel('IP tdT');
ylabel('IP IgKL-mNG');
title('IP IgKL-mNG vs IP tdT');
legend('Data', 'Robust Fit', 'Location', 'Best');

grid on;
hold off;


%% PVS

y = pvs_igkl_mean(5000:10000);
x = pvs_tdT_mean(5000:10000);
x = medfilt1(x,10);
y = medfilt1(y,10);
[b, stats] = robustfit(x, y);
y_fit = b(1) + b(2)*x;
rsquare = corr(y,b(1)+b(2)*x)^2;

figure()
plot(pvs_igkl_mean);
hold on;


plot(pvs_tdT_mean);

figure('Position', [200, 200, 600, 800])
scatter(x,y)
hold on;
% Overlay robust fit line
plot(x, y_fit, 'r-', 'LineWidth', 2);
% Optional: add labels and legend
xlabel('PVS TRITC');
ylabel('PVS IgKL-mNG');
title('PVS TRITC vs PVS IgKL-mNG');
legend('Data', 'Robust Fit', 'Location', 'Best');

grid on;
hold off;


%% PVS vs intraparenchymal space


y = ip_igkl_mean(5000:10000);
x = pvs_tdT_mean(5000:10000);
x = medfilt1(x,10);
y = medfilt1(y,10);
[b, stats] = robustfit(x, y);
y_fit = b(1) + b(2)*x;
rsquare = corr(y,b(1)+b(2)*x)^2;

% Scatter plot with robust fit overlay
figure('Position', [200, 200, 600, 800]);
scatter(x, y, 'filled');
hold on;

% Overlay robust fit line
plot(x, y_fit, 'r-', 'LineWidth', 2);

% Optional: add labels and legend
xlabel('pvs TRITC');
ylabel('ip IgKL-mNG');
title('PVS TRITC vs IP IgKL-mNG');
legend('Data', 'Robust Fit', 'Location', 'Best');

grid on;
hold off;

%% PVS vs intraparenchymal space


y = ip_igkl_mean(5000:10000);
x = pvs_igkl_mean(5000:10000);
x = medfilt1(x,10);
y = medfilt1(y,10);
[b, stats] = robustfit(x, y);
y_fit = b(1) + b(2)*x;
rsquare = corr(y,b(1)+b(2)*x)^2;

% Scatter plot with robust fit overlay
figure('Position', [200, 200, 600, 800]);
scatter(x, y, 'filled');
hold on;

% Overlay robust fit line
plot(x, y_fit, 'r-', 'LineWidth', 2);

% Optional: add labels and legend
xlabel('pvs FITC');
ylabel('ip EGFP');
title('PVS FITC vs IP EGFP');
legend('Data', 'Robust Fit', 'Location', 'Best');

grid on;
hold off;

%%  Inside of vessel
%cbv_mask((cbv_mask == 0)) = nan;
y = ip_cbv_mean(5000:10000);
x = pvs_cbv_mean(5000:10000);
figure()
plot(cbv_igkl_mean);
hold on;

plot(cbv_tdT_mean);
figure()
scatter(cbv_igkl_mean(5000:10000),cbv_tdT_mean(5000:10000))
set(gcf, 'defaultFigurePosition', [200, 200, 600, 800]);



%%

figHandles = findall(0, 'Type', 'figure');
for k = 1:length(figHandles)
    set(figHandles(k), 'Position', [100, 100, 600, 800]);
end






%% second example
%
%
%%
vivo_50um2 = mdfExtractLoader();
vivo_50um2.analog = vivo_50um2.loadanalog;
vivo_50um2.stackch1 = vivo_50um2.loadstack("ch1");
vivo_50um2.stackch2 = vivo_50um2.loadstack("ch2");


%%

%epecs : extraparenchymal extracellularspace
%ipecs : intraparenchymal extracellularspace
%cbv : when vesel constricted

%%
epecs = dilation(vivo_50um2.stackch1,"polygon",vivo_50um2);
epecs = epecs.modifyroi(vivo_50um2.stackch1);

ipecs = epecs;
ipecs = ipecs.modifyroi(vivo_50um2.stackch1);

cbv = epecs;
cbv = cbv.modifyroi(vivo_50um2.stackch1);

%%
tmp.filename = vivo_50um2.info.mdfName(1:end-4);

%%
save(fullfile(savepath,tmp.filename),'epecs','cbv','ipecs')



%%
figure()
imshow(ip_mask+cbv_mask)

%%
epecs = epecs.modifyroi(vivo_50um2.stackch2);


%% Mask generation
stack = vivo_50um2.stackch2;
[cols, rows] = meshgrid(1:size(stack, 2), 1:size(stack, 1));
cbv_mask = inpolygon(cols, rows,cbv.vertices(:, 1), cbv.vertices(:, 2));
epecs_mask = inpolygon(cols, rows,epecs.vertices(:, 1), epecs.vertices(:, 2));
ipecs_mask = inpolygon(cols, rows,ipecs.vertices(:, 1), ipecs.vertices(:, 2));
%% 
% ip_mask: concentric ring, intraparenchymal space, space just inside of parenchyma
% pvs_mask: concentric ring, perivascular space, between the space where vessel constricted and border of parencyma
% epecs_mask: circle, extraparenchymal space, BV+Endothelial+Mural cell

%% Intraparenchymal space
ip_mask = ipecs_mask - epecs_mask;
ip_mask((ip_mask == 0)) = nan;

ip_igkl = double(vivo_50um2.stackch2).*ip_mask;
ip_igkl_mean = squeeze(mean(ip_igkl,[1,2],"omitmissing"));

figure()
plot(ip_igkl_mean);
hold on;

ip_tdT = double(vivo_50um2.stackch1).*ip_mask;
ip_tdT_mean = squeeze(mean(ip_tdT,[1,2],"omitmissing"));

plot(ip_tdT_mean);

figure()
scatter(ip_tdT_mean(5000:10000),ip_igkl_mean(5000:10000))

%% PVS
pvs_mask = epecs_mask - cbv_mask;
pvs_mask((pvs_mask == 0)) = nan;

pvs_igkl = double(vivo_50um2.stackch2).*pvs_mask;
pvs_igkl_mean = squeeze(mean(pvs_igkl,[1,2],"omitmissing"));

figure()
plot(pvs_igkl_mean);
hold on;

pvs_tdT = double(vivo_50um2.stackch1).*pvs_mask;
pvs_tdT_mean = squeeze(mean(pvs_tdT,[1,2],"omitmissing"));

plot(pvs_tdT_mean);

figure()
scatter(pvs_tdT_mean(5000:10000),pvs_igkl_mean(5000:10000))

%%  Inside of vessel
%cbv_mask((cbv_mask == 0)) = nan;
cbv_igkl = double(vivo_50um2.stackch2).*cbv_mask;
cbv_igkl_mean = squeeze(mean(cbv_igkl,[1,2],"omitmissing"));
figure()
plot(cbv_igkl_mean);
hold on;
cbv_tdT = double(vivo_50um2.stackch1).*cbv_mask;
cbv_tdT_mean = squeeze(mean(cbv_tdT,[1,2],"omitmissing"));
plot(cbv_tdT_mean);
figure()
scatter(cbv_tdT_mean(5000:10000),cbv_igkl_mean(5000:10000))

%%

%%
figure()
scatter(pvs_tdT_mean(5000:10000),ip_tdT_mean(5000:10000))
%%
cbv.showroi
%%
epecs.showroi
%%
ipecs.showroi
%%



t = 0.05:0.5:2*pi;
x1 = cos(t);
y1 = sin(t);
x2 = 0.5*cos(t);
y2 = 0.5*sin(t);
pgon = polyshape({x1,x2},{y1,y2})










%% Draw three ring 

% User defined mask
% cbv_mask: constricted vessel reflect always BV portion
% epecs_mask: maximally dialated vessel + shade where intraparenchymal tdT to be not included and modify it to expand by seeing
% igkl signal contrast
% ipecs_mask: using epecs mask, expand furthermore to include some of brain
% parenchyma, shaded region..?

% Auto generation mask
% pvs_mask = epecs_mask - cbv_mask
% ip_mask = ipecs_mask - epecs_mask





%% Draw three ring







%%

figure()
scatter(pvs_tdT_mean(5000:10000),ip_igkl_mean(5000:10000))




%%
figure('Name','pvs_mask')
imshow(cbv_mask)
figure('Name','pvs_mask')
imshow(cbv_mask)
%%













%%
cbv_igkl = vivo_50um2.stackch2.*uint16(cbv_mask);
%%




%%
ip_tdt = vivo_50um2.stackch1.*uint16(cbv_mask);
%%
bv_tdt = vivo_50um2.stackch1.*uint16(epecs_mask);


%%

figure()
plot(squeeze(mean(cbv_igkl,[1,2],"omitmissing")))
hold on 
plot(squeeze(mean(cbv_igkl,[1,2],'omitmissing')))
plot(squeeze(mean(ip_tdt,[1,2],'omitmissing')))
plot(squeeze(mean(bv_tdt,[1,2],'omitmissing')))

%%














%%

cbv.modifyroi(vivo_50um2.stackch2)



%%
vertices = pvsbv4.vertices;
%%

stack = vivo_50um2.stackch2;

%%

[cols, rows] = meshgrid(1:size(stack, 2), 1:size(stack, 1));
mask4 = inpolygon(cols, rows, vertices(:, 1), vertices(:, 2));

%%

figure()
imshow(mask2-mask4)

%%







%%
pvsbv1 = pvsbv1.modifyroi(vivo_50um2.stackch2);
pvsbv1.stacks.tritc_cagtdT = pvsbv1.addstack(vivo_50um2.stackch1);
pvsbv1.stacks.igklmNG = pvsbv1.addstack(vivo_50um2.stackch2);

%%
pvsbv1.showroi
%%

pvsbv1 = pvsbv1.modifyroi(vivo_50um2.stackch1);

%%

pvsbv1.showroi

%%

pvsbv1.showstack('tritc_cagtdT')
pvsbv1.showstack('igklmNG')

%%

pvsbv1.data.mean3digklmNG = mean(pvsbv1.stacks.igklmNG,"all","omitmissing");
pvsbv1.data.mean2digklmNG = squeeze(mean(pvsbv1.stacks.igklmNG,[1,2],"omitmissing"));
pvsbv1.data.normigklmNG = pvsbv1.data.mean2digklmNG/pvsbv1.data.mean3digklmNG;
%%
figure()
plot(pvsbv1.data.normigklmNG);

%%

pvsbv1.data.mean3dtritc = mean(pvsbv1.stacks.tritc_cagtdT,"all","omitmissing");
pvsbv1.data.mean2dtritc = squeeze(mean(pvsbv1.stacks.tritc_cagtdT,[1,2],"omitmissing"));
pvsbv1.data.normtritc = pvsbv1.data.mean2dtritc/pvsbv1.data.mean3dtritc;
%%
hold on
plot(pvsbv1.data.normtritc)

%% aaaaaaaaaa


pvs1 = dilation(vivo_50um2.stackch1,"polygon",vivo_50um2);
%%
pvs1 = pvs1.modifyroi(vivo_50um2.stackch2);
%%
pvs1.stacks.tritc_cagtdT = pvs1.addstack(vivo_50um2.stackch1);
pvs1.stacks.igklmNG = pvs1.addstack(vivo_50um2.stackch2);

%%
pvs1.showroi
%%

pvs1 = pvs1.modifyroi(vivo_50um2.stackch2);

%%

pvs1.showroi

%%

pvs1.showstack('tritc_cagtdT')
pvs1.showstack('igklmNG')

%%

pvs1.data.mean3digklmNG = mean(pvs1.stacks.igklmNG,"all","omitmissing");
pvs1.data.mean2digklmNG = squeeze(mean(pvs1.stacks.igklmNG,[1,2],"omitmissing"));
pvs1.data.normigklmNG = pvs1.data.mean2digklmNG/pvs1.data.mean3digklmNG;
%%
figure()
plot(pvs1.data.normigklmNG);

%%
pvs1.showstack('igklmNG')
%%

pvs1.showstack('tritc_cagtdT')

%%
pvs1.data.mean3dtritc = mean(pvs1.stacks.tritc_cagtdT,"all","omitmissing");
pvs1.data.mean2dtritc = squeeze(mean(pvs1.stacks.tritc_cagtdT,[1,2],"omitmissing"));
pvs1.data.normtritc = pvs1.data.mean2dtritc/pvs1.data.mean3dtritc;
%%
hold on
plot(pvs1.data.normtritc)






















%%
pvs.stacks.syn1gfap_tdt_inverted = max(pvs.stacks.cag_egfp,[],'all')-pvs.stacks.syn1gfap_tdt;
pvs.stacks.syn1gfap_tdt_inverted = pre_groupaverage(pvs.stacks.syn1gfap_tdt_inverted,5);
%%

pvs = pvs.radonthresholding('syn1gfap_tdt_inverted');

%%

figure("Name",'iradon pvs')
sliceViewer(pvs.stacks.irtd_syn1gfap_tdt_inverted)

%%

figure()
sliceViewer(pvs.stacks.syn1gfap_tdt_inverted)
%%

tmp.groupaverage = pre_groupaverage(pvs.stacks.syn1gfap_tdt_inverted,5);

figure()
sliceViewer(tmp.groupaverage)


%%

pvs = pvs.radonthresholding('cag_egfp_inverted');

%%
figure()
sliceViewer(ans)


%%
pvs = dilation(vivo_50um2.stackch1,'rectangle',vivo_50um2);
%%
pvs.stacks.syn1gfap_tdt = pvs.addstack(vivo_50um2.stackch1);
%%
pvs.stacks.cag_egfp = pvs.addstack(vivo_50um2.stackch2);

%%
figure()
plot(squeeze(max(pvs.stacks.syn1gfap_tdt,[],[1,2])))
%%
figure()
plot(squeeze(max(pvs.stacks.cag_egfp,[],[1,2])))


%%

% saveprefix
tmp.wavelength = vivo_50um2.info.excitation(1:end-3);
tmp.power = vivo_50um2.info.laserpower(1:end-1);
tmp.filename = vivo_50um2.info.mdfName(1:end-4);
saveprefix = [tmp.filename,'_lp',tmp.power,'_ex',tmp.wavelength];



figure("Name",'Ch1')
vivoViewer(vivo_50um2.stackch1)
figure("Name",'Ch2')
vivoViewer(vivo_50um2.stackch2)


ch1name = 'tdT';
ch2name = 'mNG';

%% Get User Input for ROI Name
roiName = input('Enter ROI name (e.g., cellbody1, neuropile1, void1, pvs1): ', 's');

% Create ROI dynamically
bleachData.(roiName) = bleach(vivo_50um2.stackch2, 'polygon', vivo_50um2);
%%
bleachData.(roiName) = bleachData.(roiName).modifyroi(vivo_50um2.stackch1);


% Assign stacks
bleachData.(roiName).stacks.(ch1name) = bleachData.(roiName).addstack(vivo_50um2.stackch1);
bleachData.(roiName).stacks.(ch2name) = bleachData.(roiName).addstack(vivo_50um2.stackch2);
%
% Compute Decay
bleachData.(roiName).data = bleachData.(roiName).decay(ch2name);
bleachData.(roiName).data = bleachData.(roiName).decay(ch1name);
%
bleachData.(roiName).showroi

% Save Data Using ROI Name
saveroi = bleachData.(roiName);
%
save(fullfile(savepath, [saveprefix, '_', roiName, '.mat']), "saveroi",'-v7.3');

disp(['Data saved: ', saveprefix, '_', roiName, '.mat']);

%%
figure()
plot(bleachData.(roiName).data.norm_mNG)
hold on
plot(bleachData.(roiName).data.norm_tdT)
%%
figure()
plot(bleachData.cellbody.data.exponential2_mNG)
hold on
plot(bleachData.cellbody.data.exponential2_tdT)

