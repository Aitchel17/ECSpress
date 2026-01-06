
% Photobleach experiment analysis script



    % extra-mNG + intra-tdT
    % intra-mNG + intra-tdT (control)

% 1. Fixed sample
    % 1.1 Brain slice (proof of concept)
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
savepath = 'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\01_primary_analysis\00_bleach\slice';


slice = mdfExtractLoader();
slice.analog = slice.loadanalog;
slice.stackch1 = slice.loadstack("ch1");
slice.stackch2 = slice.loadstack("ch2");

% saveprefix
tmp.wavelength = slice.info.excitation(1:end-3);
tmp.power = slice.info.laserpower(1:end-1);
tmp.filename = slice.info.mdfName(1:end-4);
saveprefix = [tmp.filename,'_lp',tmp.power,'_ex',tmp.wavelength];



figure("Name",'Ch1')
sliceViewer(slice.stackch1)
figure("Name",'Ch2')
sliceViewer(slice.stackch2)


ch1name = 'tdT';
ch2name = 'mNG';

%% Get User Input for ROI Name
roiName = input('Enter ROI name (e.g., cellbody1, neuropile1, void1, pvs1): ', 's');

% Create ROI dynamically
bleachData.(roiName) = bleach(slice.stackch2, 'polygon', slice);
%%
bleachData.(roiName) = bleachData.(roiName).modifyroi(slice.stackch1);


% Assign stacks
bleachData.(roiName).stacks.(ch1name) = bleachData.(roiName).addstack(slice.stackch1);
bleachData.(roiName).stacks.(ch2name) = bleachData.(roiName).addstack(slice.stackch2);
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

