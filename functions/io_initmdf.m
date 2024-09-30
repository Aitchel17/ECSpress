%%%%%%%%%%%%%%%%%%%%%%%%
% Prerequisite: MCSX 
%                       >> https://www.sutter.com/MICROSCOPES/mcs.html

% FUNCTION NAME:    io_tpsm
%
% DESCRIPTION:      Sutter .mdf Data input - output function
% INPUT:            .mdf
%
% NOTES:            If wavelength and objective length 
%
% WRITTEN BY:       C. Hyunseok Lee 2024-09-14
%
%%%%%%%%%%%%%%%%%%%%%%%%
function [info, tpsm_info, mode_info, mobj] = io_initmdf()
    %%
    [info.mdfName, info.mdfPath] = uigetfile({'*.tif ; *.mdf'}); % select file by UI
    mdfPath = [info.mdfPath, info.mdfName];
    mobj = actxserver('MCSX.Data'); % Create Component Object Model (COM)
    mobj.invoke('OpenMCSFile', mdfPath); % Using COM open .mdf file
    
    % General Info
    info.User = mobj.ReadParameter('Created by');
    info.Date = mobj.ReadParameter('Created on');
    
    % two photon scanning microscope info
    tpsm_info.scanmode      = mobj.ReadParameter('Scan Mode');
    tpsm_info.pclock        = mobj.ReadParameter('Pixel Clock');
    tpsm_info.yoffset       = mobj.ReadParameter('Y Frame Offset');
    
    % Objective info
    tpsm_info.objname       = mobj.ReadParameter('Objective');
    tpsm_info.objpix        = mobj.ReadParameter('Microns per Pixel');
    tpsm_info.objx          = mobj.ReadParameter('X Position');
    tpsm_info.objy          = mobj.ReadParameter('Y Position');
    tpsm_info.objz          = mobj.ReadParameter('Z Position');
    
    % Scan info
    tpsm_info.zoom          = mobj.ReadParameter('Magnification');
    tpsm_info.tpex          = strcat(mobj.ReadParameter('Laser Wavelength (nm)'),' nm');
    tpsm_info.fbit          = mobj.ReadParameter('Frame Bit Depth');
    tpsm_info.fdur          = mobj.ReadParameter('Frame Duration (s)');
    tpsm_info.fcount        = str2double(mobj.ReadParameter('Frame Count'));
    tpsm_info.fint          = mobj.ReadParameter('Frame Interval (ms)');
    tpsm_info.fps           = 1/str2double(tpsm_info.fdur(1:end-1)); % Hz
    tpsm_info.fh            = str2double(mobj.ReadParameter('Frame Height'));
    tpsm_info.fw            = str2double(mobj.ReadParameter('Frame Width'));
    tpsm_info.lpower        = mobj.ReadParameter('Laser intensity');
    % Imaging Channel info
    tpsm_info.imgch0range   = mobj.ReadParameter('Scanning Ch 0 Input Range');
    tpsm_info.imgch1range   = mobj.ReadParameter('Scanning Ch 1 Input Range');



    % Imaging channel info
    %tpsm_info.image0        =


%% Scan mode specific info
    if strcmp(tpsm_info.scanmode, 'Image Stack')
        disp('Image stack loaded')
        mode_info.fave      = mobj.ReadParameter('Average Count');
        mode_info.pinit     = mobj.ReadParameter('Initial Intensity');
        % mode_info.pfinl     = mobj.ReadParameter('Final Intensity');
        % final intensity activex control has bug the .ocx file should be
        % editted
        mode_info.zinter    = mobj.ReadParameter('Z- interval');
        
    elseif strcmp(tpsm_info.scanmode, 'XY Movie')
        mode_info.analogfreq    = mobj.ReadParameter('Analog Acquisition Frequency (Hz)');
        mode_info.analogfreq    = str2double(mode_info.analogfreq(1:end -3));
        mode_info.analogcount   = str2double(mobj.ReadParameter('Analog Sample Count'));
        mode_info.analogresolution   = mobj.ReadParameter('Analog Resolution');

        mode_info.analog0   = mobj.ReadParameter('Analog Ch 0 Name');
        mode_info.analog1   = mobj.ReadParameter('Analog Ch 1 Name');
        mode_info.analog2   = mobj.ReadParameter('Analog Ch 2 Name');
        mode_info.analog3   = mobj.ReadParameter('Analog Ch 3 Name');
        mode_info.analog4   = mobj.ReadParameter('Analog Ch 4 Name');
        mode_info.analogch4range = mobj.ReadParameter('Analog Ch 4 Input Range');
        mode_info.analog5   = mobj.ReadParameter('Analog Ch 5 Name');
        mode_info.analog6   = mobj.ReadParameter('Analog Ch 6 Name');
        mode_info.analog7   = mobj.ReadParameter('Analog Ch 7 Name');
    end

end