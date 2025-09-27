

% This function is to load data set for primary analysis
% The mission is not very solid at this point, it can be 
% 1. fwhm calculation for BV, and PVS thickness n for
% 2. Radon based calculation of all angle vessel expension
% 3. Dip direction angle
% all angles

function raw_datastruct = primary_loaddata(test_dir)
% Primary_loaddata initialize data structure required for roi
% based analysis with matlab
% primary init generate folder which will contain
    % 1. analog.mat - analysis_analog class
    % 2. roilist.mat - roi class
    % 3. behavior.mat - behav class
    % 4. 
    %% 
    % 4. linefwhm.mat  this is for kymograph based analysis
    % 5. radon.mat
%   Detailed explanation goes here

    if nargin == 0
        mdfextract = mdfExtractLoader();
    elseif nargin == 1
        mdfextract = mdfExtractLoader(test_dir);
    end
    raw_datastruct = struct();
    stackch1 = mdfextract.loadstack("ch1");
    stackch2 = mdfextract.loadstack("ch2");
    loaded_data = util_load_primarydataset(mdfextract.info.analyzefolder);
    primary_analog = add_analog(loaded_data,mdfextract);
    roilist = add_roilist(loaded_data);
        
    % original time axis
    img_param.fps = str2double(mdfextract.info.savefps); % mdfExtractor saving fps
    img_param.pixel2um = str2double(mdfextract.info.objpix(1:end-2));
    img_param.loadstart = str2double(mdfextract.info.loadstart)/img_param.fps; % load start frame --> load start seconds
    img_param.taxis = linspace(img_param.loadstart,img_param.loadstart+size(stackch1,3)/img_param.fps,size(stackch1,3));

    raw_datastruct.mdfextract = mdfextract;
    raw_datastruct.primaryanalog = primary_analog;
    raw_datastruct.roilist = roilist;
    raw_datastruct.img_param = img_param;
    raw_datastruct.stackch1 = stackch1;
    raw_datastruct.stackch2 = stackch2;

end

function loadstruct = util_load_primarydataset(extractfolder_path)
%UTIL_LOAD_PRIMARYDATASET Summary of this function goes here
%   Detailed explanation goes here
%   check existence of primary anlysis and make folder
    savepath = fullfile(extractfolder_path,'primary_analysis');
    loadstruct = struct();
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end

    % check existence of line_fwhms 
    if isfile(fullfile(extractfolder_path,'primary_analysis/line_fwhms.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/line_fwhms.mat'));
        loadstruct.line_fwhms = tmp.load.line_fwhms;
    end

    % check existence of line_fwhms 
    if isfile(fullfile(extractfolder_path,'primary_analysis/analog.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/analog.mat'));
        loadstruct.analog = tmp.load.primary_analog;
    end

    % check existence of line_fwhms 
    if isfile(fullfile(extractfolder_path,'primary_analysis/roilist.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/roilist.mat'));
        loadstruct.roilist = tmp.load.roilist;
    end
end

