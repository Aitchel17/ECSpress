function primary_struct = primary_integration(folderpath,targetfile)
    %FWHM_INTEGRATION Summary of this function goes here
    %   Detailed explanation goes here
    folderpath = fullfile(folderpath,'**');
    tmp.dir_list = struct2table(dir(folderpath));
    
    tmp.line_dirloc = matches(tmp.dir_list.name,targetfile);

    tmp.line_list = tmp.dir_list(tmp.line_dirloc,:);
    
    tmp.info_dirloc = endsWith(tmp.dir_list.name, '_info.txt');
    tmp.info_list = tmp.dir_list(tmp.info_dirloc,:);
    
    tmp.analog_dirloc = matches(tmp.dir_list.name, 'analog.mat');
    tmp.analog_list = tmp.dir_list(tmp.analog_dirloc,:);
    
    %% path
    primary_struct = repmat(struct(), 1, height(tmp.line_list));
    for idx = 1:height(tmp.line_list)
        tmp.linedata_path = fullfile(tmp.line_list.folder{idx},tmp.line_list.name{idx});
        % using linedata folder information extract parent folder name to
        % access info.txt or current folder directories analog.nat
        tmp.info_listloc = matches(tmp.info_list.folder, fileparts(tmp.line_list.folder{idx}));
        tmp.info_dir = tmp.info_list(tmp.info_listloc,:);
        tmp.info_path = fullfile(tmp.info_dir.folder{1},tmp.info_dir.name{1});
        disp(tmp.info_path)
        tmp.analog_listloc = matches(tmp.analog_list.folder,tmp.line_list.folder{idx});
        tmp.analog_dir = tmp.analog_list(tmp.analog_listloc,:);
        tmp.analog_path = fullfile(tmp.analog_dir.folder{1},tmp.analog_dir.name{1});
            
        tmp.info_table = readtable(tmp.info_path);
        tmp.mdforigin = tmp.info_table(matches(tmp.info_table.Field,'mdfName'),:).Value{1};
    
        tmp.fieldname = tmp.mdforigin(1:end-4);
        disp(tmp.fieldname)
        % information dictionary construction
        tmp.infodict = dictionary(string(tmp.info_table.Field),string(tmp.info_table.Value));
        tmp.comment = tmp.infodict('Comments');
        tmp.split_comment = strsplit(tmp.comment, ' ');
        tmp.vesselid = tmp.split_comment{1};
        tmp.depth = regexp(tmp.comment, '\d+\s*um', 'match', 'once', 'ignorecase');
        tmp.depth = regexp(tmp.depth, '\d+', 'match', 'once');
        tmp.mwpower = regexp(tmp.comment, '\d+\s*mw', 'match', 'once', 'ignorecase');
        tmp.mwpower = regexp(tmp.mwpower, '\d+', 'match', 'once');
        tmp.rcv = regexp(tmp.comment, 'rcv\s*\.?(\d+)', 'tokens', 'once', 'ignorecase');
        tmp.rcv = tmp.rcv{1};
        tmp.gcv = regexp(tmp.comment, 'gcv\s*\.?(\d+)', 'tokens', 'once', 'ignorecase');
        tmp.gcv = tmp.gcv{1};
        
        tmp.infodict = insert(tmp.infodict, ["vesselid", "depth","power","rcv","gcv"],...
            [tmp.vesselid, tmp.depth, tmp.mwpower, tmp.rcv, tmp.gcv]);
        primary_struct(idx).sessionid = tmp.fieldname;
        primary_struct(idx).infodict = tmp.infodict;
        % load primaryanlysis data
        primary_struct(idx).analog = load(tmp.analog_path).primary_analog;
        % primary_struct(idx).fwhmline = load(tmp.linedata_path).line_fwhm;
        primary_struct(idx).manualpolar = load(tmp.linedata_path).manual_polarstruct;    

    end

end

