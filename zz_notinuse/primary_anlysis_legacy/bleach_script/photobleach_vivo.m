
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
        
%% load file
savepath = 'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\01_primary_analysis\00_bleach\vivo';


vivo = mdfExtractLoader();
vivo.analog = vivo.loadanalog;
vivo.stackch1 = vivo.loadstack("ch1");
vivo.stackch2 = vivo.loadstack("ch2");


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
pvs = dilation(vivo.stackch1,'rectangle',vivo);
%%
pvs.stacks.syn1gfap_tdt = pvs.addstack(vivo.stackch1);
%%
pvs.stacks.cag_egfp = pvs.addstack(vivo.stackch2);

%%
figure()
plot(squeeze(max(pvs.stacks.syn1gfap_tdt,[],[1,2])))
%%
figure()
plot(squeeze(max(pvs.stacks.cag_egfp,[],[1,2])))


%%

% saveprefix
tmp.wavelength = vivo.info.excitation(1:end-3);
tmp.power = vivo.info.laserpower(1:end-1);
tmp.filename = vivo.info.mdfName(1:end-4);
saveprefix = [tmp.filename,'_lp',tmp.power,'_ex',tmp.wavelength];



figure("Name",'Ch1')
vivoViewer(vivo.stackch1)
figure("Name",'Ch2')
vivoViewer(vivo.stackch2)


ch1name = 'tdT';
ch2name = 'mNG';

%% Get User Input for ROI Name
roiName = input('Enter ROI name (e.g., cellbody1, neuropile1, void1, pvs1): ', 's');

% Create ROI dynamically
bleachData.(roiName) = bleach(vivo.stackch2, 'polygon', vivo);
%%
bleachData.(roiName) = bleachData.(roiName).modifyroi(vivo.stackch1);


% Assign stacks
bleachData.(roiName).stacks.(ch1name) = bleachData.(roiName).addstack(vivo.stackch1);
bleachData.(roiName).stacks.(ch2name) = bleachData.(roiName).addstack(vivo.stackch2);
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

