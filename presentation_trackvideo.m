% This script requires twophoton_processed (the input for FWHM) and FWHM info.
% ch1: BV, ch2: CSF, ch3: BV boundary, ch4: CSF boundary.

fwhm = session.pax_fwhm;

%% Reconstruct + overlay boundaries, save 4-channel TIFF
overlay = fwhm.reconstruction_overlay( ...
    twophoton_processed.ch1, ...   % -> output channel 1 (BV)
    twophoton_processed.ch2, ...   % -> output channel 2 (CSF)
    BoundaryIntensity = 2^16, ...
    Pixel2um          = twophoton_processed.pixel2um, ...
    Fps               = twophoton_processed.fps, ...
    SavePath          = fullfile(session.dir_struct.primary_analysis, "fwhm_boundary.tif"));

%%
util_checkstack(overlay(:, :, :, 1))

%% Reference views
util_checkstack(fwhm.rotatecrop.rc_lumen)
%%
figure(); imagesc(fwhm.rotatecrop.rc_lumen(:,:,1618)); axis image; colormap gray; axis off
%%
figure(); plot(fwhm.kymograph.kgph_lumen_processed(:,1618),Color='k',LineWidth=2)
hold on; yline(0.5,LineWidth=1, Color='r'),axis off
