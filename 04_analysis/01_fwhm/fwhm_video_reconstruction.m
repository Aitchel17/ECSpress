% FWHM_VIDEO_RECONSTRUCTION  Overlay reconstructed FWHM boundaries on the 2-photon
% stack and save an ImageJ hyperstack TIFF.
%
% Requires in the workspace:
%   session            - with .pax_fwhm (line_fwhm) and .img_param.pixel2um
%   twophoton_processed - with .ch1, .ch2 and .outfps
%   sessiondir         - output folder
%
% The overlay/reconstruction/save logic now lives in the line_fwhm method
% reconstruction_overlay.

overlay = session.pax_fwhm.reconstruction_overlay( ...
    twophoton_processed.ch2, ...   % -> output channel 1
    twophoton_processed.ch1, ...   % -> output channel 2
    BoundaryIntensity = 2048, ...
    Pixel2um          = session.img_param.pixel2um, ...
    Fps               = twophoton_processed.outfps, ...
    SavePath          = fullfile(sessiondir, "fwhm_tracking.tif"));

util_checkoverlay(overlay)
