function loadstruct = load_primaryresult(extractfolder_path)

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


