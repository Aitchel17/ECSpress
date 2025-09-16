%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION NAME:    Calibration_PSF
%
% DESCRIPTION:      Calculate FWHM of point spreadfunction
%
% INPUT:            Stack images of point source
%
% VARIABLES:        fileNames: (string) full path to image file
%                   beadSize: (double) size of bead (um)
%                   : (double) width of grid holes
%                   minPeakDistanceAuto: (double) minimum distance between peaks when finding intensity peaks during automation
%                   objMag: (double) magnification of objective used
%                   digMag: (double) digital magnification used
%                   turnabout: (double) turnabout value of microscope
%                   commentString: (string) general comments for record keeping
%
% OUTPUT:
%
% FUNCTIONS USED:   [width,ttrail,tlead,Pol] = fwhm_calibGrid(x,y);

%                   Adapted from: Patrick Egan (2024). fwhm (https://www.mathworks.com/matlabcentral/fileexchange/10590-fwhm), MATLAB Central File Exchange. Retrieved April 3, 2024. 
%
% LIBARIES USED:
%
% NOTES:            Any grid can be used that has a sufficiently small grid pattern for the objective. For a 16x objective on a
%                   two-photon microscope, the following grids were used:
%                   1000 (lines per inch) Mesh Copper Grid (2145C from SPI supplies): bar size = 6 microns, hole size = 19 microns
%                   2000 (lines per inch) Mesh Copper Grid (2155C from SPI supplies): bar size = 5 microns, hole size = 7.5 microns
%                   
%                   Draw box for each full width at half maximum in only one row or column but as wide as possible for best
%                   smoothing of intensity curve with averaging.
%
%                   As of now, manual selection of points on intensity curve can calculate full width at half maximum in
%                   either polarity (bar or hole width) depending on points chosen. Auto selection will be limited to single polarity
%                   due to relying on finding peaks for calculation (hole width only).
%
% WRITTEN BY:       C. Hyunseok Lee 2024-09-14
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%