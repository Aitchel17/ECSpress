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
function [info, mobj] = io_initmdf()
    %%
    [info.mdfName, info.mdfPath] = uigetfile({'*.mdf'}); % select file by UI
    mdfPath = [info.mdfPath, info.mdfName];
    mobj = actxserver('MCSX.Data'); % Create Component Object Model (COM)
    mobj.invoke('OpenMCSFile', mdfPath); % Using COM open .mdf file
    
    % General Info
    info.User = mobj.ReadParameter('Created by');
    info.Date = mobj.ReadParameter('Created on');
    
    % two photon scanning microscope info
    info.scanmode      = mobj.ReadParameter('Scan Mode');
    info.pclock        = mobj.ReadParameter('Pixel Clock');
    info.yoffset       = mobj.ReadParameter('Y Frame Offset');
    
    % Objective info
    info.objname       = mobj.ReadParameter('Objective');
    info.objpix        = mobj.ReadParameter('Microns per Pixel');
    info.objx          = mobj.ReadParameter('X Position');
    info.objy          = mobj.ReadParameter('Y Position');
    info.objz          = mobj.ReadParameter('Z Position');
    
    % Scan info
    info.zoom          = mobj.ReadParameter('Magnification');
    info.tpex          = strcat(mobj.ReadParameter('Laser Wavelength (nm)'),' nm');
    info.fbit          = mobj.ReadParameter('Frame Bit Depth');
    info.fdur          = mobj.ReadParameter('Frame Duration (s)');
    info.fcount        = str2double(mobj.ReadParameter('Frame Count'));
    info.fint          = mobj.ReadParameter('Frame Interval (ms)');
    info.fps           = 1/str2double(info.fdur(1:end-1)); % Hz
    info.fh            = str2double(mobj.ReadParameter('Frame Height'));
    info.fw            = str2double(mobj.ReadParameter('Frame Width'));
    info.lpower        = mobj.ReadParameter('Laser intensity');
    % Imaging Channel info
    info.imgch0name    = mobj.ReadParameter('Scanning Ch 0 Name');
    info.imgch0range   = mobj.ReadParameter('Scanning Ch 0 Input Range');
    info.imgch1name    = mobj.ReadParameter('Scanning Ch 1 Name');
    info.imgch1range   = mobj.ReadParameter('Scanning Ch 1 Input Range');



    % Imaging channel info
    %tpsm_info.image0        =


%% Scan mode specific info
    if strcmp(info.scanmode, 'Image Stack')
        disp('Image stack loaded')
        info.fave      = mobj.ReadParameter('Average Count');
        info.pinit     = mobj.ReadParameter('Initial Intensity');
        % mode_info.pfinl     = mobj.ReadParameter('Final Intensity');
        % final intensity activex control has bug the .ocx file should be
        % editted
        info.zinter    = mobj.ReadParameter('Z- interval');
        
    elseif strcmp(info.scanmode, 'XY Movie')
        info.analogfreq    = mobj.ReadParameter('Analog Acquisition Frequency (Hz)');
        info.analogfreq    = str2double(info.analogfreq(1:end -3));
        info.analogcount   = str2double(mobj.ReadParameter('Analog Sample Count'));
        info.analogresolution   = mobj.ReadParameter('Analog Resolution');
        info.analog0   = mobj.ReadParameter('Analog Ch 0 Name');
        info.analog1   = mobj.ReadParameter('Analog Ch 1 Name');
        info.analog2   = mobj.ReadParameter('Analog Ch 2 Name');
        info.analog3   = mobj.ReadParameter('Analog Ch 3 Name');
        info.analog4   = mobj.ReadParameter('Analog Ch 4 Name');
        info.analog5   = mobj.ReadParameter('Analog Ch 5 Name');
        info.analog6   = mobj.ReadParameter('Analog Ch 6 Name');
        info.analog7   = mobj.ReadParameter('Analog Ch 7 Name');
        for analog_ch = string(0:1:7)
            field_name = strcat('analog',analog_ch);
            if strcmp(info.(field_name),"")
                info = rmfield(info,field_name);
            else
                range_fname = strcat('analog',analog_ch,'range');
                range_call = strcat('Analog ', 'ch ',analog_ch,' Input Range');
                info.(range_fname) = mobj.ReadParameter(range_call);
            end
        end
    end



end