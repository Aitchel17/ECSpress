

% This function is to load data set for primary analysis
% The mission is not very solid at this point, it can be 
% 1. fwhm calculation for BV, and PVS thickness n for
% 2. Radon based calculation of all angle vessel expension
% 3. Dip direction angle
% all angles

function output_datastruct = primary_loaddata(test_dir)
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
    %%


    %%
    output_datastruct = util_load_primarydataset(mdfextract.info.analyzefolder, mdfextract);

    stackch1 = mdfextract.loadstack("ch1");
    stackch2 = mdfextract.loadstack("ch2");
        
    % original time axis
    img_param.save_fps = str2double(mdfextract.info.savefps); % mdfExtractor saving fps
    img_param.record_fps = str2double(mdfextract.info.fps); %
    img_param.imgstartframe = str2double(mdfextract.info.loadstart);
    img_param.imgstarttime = img_param.imgstartframe/img_param.record_fps;
    img_param.imgendframe = str2double(mdfextract.info.loadend);
    img_param.imgendtime = img_param.imgendframe/img_param.record_fps;
    img_param.pixel2um = str2double(mdfextract.info.objpix(1:end-2));
    img_param.taxis = linspace(img_param.imgstarttime,img_param.imgendtime,size(stackch1,3));


    output_datastruct.mdfextract = mdfextract;
    output_datastruct.img_param = img_param;
    output_datastruct.stackch1 = stackch1;
    output_datastruct.stackch2 = stackch2;

end

function loadstruct = util_load_primarydataset(extractfolder_path, mdfextract)
%UTIL_LOAD_PRIMARYDATASET Summary of this function goes here
%   Detailed explanation goes here
%   check existence of primary anlysis and make folder
    savepath = fullfile(extractfolder_path,'primary_analysis');
    loadstruct = struct();
    if ~exist(savepath, 'dir')
        mkdir(savepath);
    end
    % check existence of line_fwhms 
    if isfile(fullfile(extractfolder_path,'primary_analysis/analog.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/analog.mat'));
        loadstruct.analog = tmp.load.primary_analog;
    else % if not exist make it
            disp('analog file does not exist, generate new one')
            analog = mdfextract.loadanalog;
            primary_analog =  analysis_analog(analog.info,analog.data);
            if isfield(primary_analog.data,'raw_Air_puff1')
                primary_analog.airtable = primary_analog.get_airtable('raw_Air_puff1');
            end
            if isfield(primary_analog.data,'raw_ECoG')
                primary_analog.ecogspectrum = primary_analog.get_ecogspectrum('raw_ECoG');
            end

            save(fullfile(mdfextract.info.analyzefolder,"primary_analysis","analog.mat"),"primary_analog")
    end

    % check existence of line_fwhms 
    if isfile(fullfile(extractfolder_path,'primary_analysis/line_fwhm.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/line_fwhm.mat'));
        loadstruct.line_fwhms = tmp.load.line_fwhm;
    end

    

    % check existence of line_fwhms 
    if isfile(fullfile(extractfolder_path,'primary_analysis/roilist.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/roilist.mat'));
        loadstruct.roilist = tmp.load.roilist;
    end

    % check existence of line_fwhms 
    if isfile(fullfile(extractfolder_path,'primary_analysis/pax_cluster.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/pax_cluster.mat'));
        loadstruct.pax_cluster = tmp.load.pax_cluster;
    end

    if isfile(fullfile(extractfolder_path,'primary_analysis/pax_cluster.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/pax_cluster.mat'));
        loadstruct.pax_cluster = tmp.load.pax_cluster;
    end

    if isfile(fullfile(extractfolder_path,'primary_analysis/manual_polar.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/manual_polar.mat'));
        loadstruct.manual_polar_roilist = tmp.load.roilist;
    end

    if isfile(fullfile(extractfolder_path,'primary_analysis/manual_polarstruct.mat'))
        tmp.load = load(fullfile(extractfolder_path,'primary_analysis/manual_polarstruct.mat'));
        loadstruct.manual_polarstruct = tmp.load.manual_polarstruct;
    end
end

