

% This function is to load data set for primary analysis
% The mission is not very solid at this point, it can be 
% 1. fwhm calculation for BV, and PVS thickness n for
% 2. Radon based calculation of all angle vessel expension
% 3. Dip direction angle
% all angles

function output_datastruct = load_mdfextract(test_dir)
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



