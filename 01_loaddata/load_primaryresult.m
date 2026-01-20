function loadstruct = load_primaryresult(extractfolder_path)

savepath = fullfile(extractfolder_path,'primary_analysis');
loadstruct = struct();
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

% check existence of roilist
if isfile(fullfile(extractfolder_path,'primary_analysis/roilist.mat'))
    tmp.load = load(fullfile(extractfolder_path,'primary_analysis/roilist.mat'));
    loadstruct.roilist = tmp.load.roilist;
end

% check existence of line_fwhms
if isfile(fullfile(extractfolder_path,'primary_analysis/line_fwhm.mat'))
    tmp.load = load(fullfile(extractfolder_path,'primary_analysis/line_fwhm.mat'));
    loadstruct.line_fwhms = tmp.load.line_fwhm;
end




% check existence of polarcluster
if isfile(fullfile(extractfolder_path,'primary_analysis/polarcluster.mat'))
    tmp.load = load(fullfile(extractfolder_path,'primary_analysis/polarcluster.mat'));
    loadstruct.polarcluster = tmp.load.polarcluster;
end


if isfile(fullfile(extractfolder_path,'primary_analysis/radon_result.mat'))
    tmp.load = load(fullfile(extractfolder_path,'primary_analysis/radon_result.mat'));
    loadstruct.radon_result = tmp.load.radon_result;
end

end


