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
function [info, analog, mobj] = io_initmdf()
    %%
    [info.mdfName, info.mdfPath] = uigetfile({'*.mdf'}); % select file by UI
    mdfPath = [info.mdfPath, info.mdfName];
    disp([info.mdfName,'is loaded'])
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
        info.analogcount   = str2double(mobj.ReadParameter('Analog Sample Count'));
        info.analogresolution   = mobj.ReadParameter('Analog Resolution');

        for analog_ch = string(0:1:7)
            field_name = mobj.ReadParameter(sprintf('Analog Ch %s Name',analog_ch));

            if strcmp(field_name,"")
                fprintf('\nAnalog Ch %s not detected',analog_ch);
            else
                channel_fname = [field_name,'_channel'];
                info.(channel_fname) = str2double(analog_ch);
                range_fname = [field_name,'inputrange'];
                range_call = sprintf('Analog Ch %s Input Range',analog_ch);
                info.(range_fname) = mobj.ReadParameter(range_call);

               
                fprintf('\nAnalog Ch %s is %s\n',analog_ch,field_name);

                analog.(['raw_',field_name]) = double(mobj.ReadAnalog(str2double(analog_ch)+1,info.analogcount,0));

            end
        end
    end



end