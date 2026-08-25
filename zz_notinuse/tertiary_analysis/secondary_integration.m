function secondary_struct = secondary_integration(folderpath)
%SECONDARY_INTEGRATION Summary of this function goes here
%   Detailed explanation goes here
tmp.dir_list = struct2table(dir(folderpath));
tmp.secondary_dirloc = matches(tmp.dir_list.name,'paxfwhm_struct.mat');
tmp.secondary_list = tmp.dir_list(tmp.secondary_dirloc,:);
tmp.primary_dirloc = matches(tmp.dir_list.name,'primarystruct.mat');
tmp.primary_list = tmp.dir_list(tmp.primary_dirloc,:);

%
secondary_struct=repmat(struct(),1,2);
tmp.paxsecondarypath = fullfile(tmp.secondary_list.folder{1},tmp.secondary_list.name{1});

tidx = 1;
for midx = 1: height(tmp.secondary_list)
    tmp.paxsecondarypath = fullfile(tmp.secondary_list.folder{midx},tmp.secondary_list.name{midx});
    z = load(tmp.paxsecondarypath);
    for sidx = 1:size(z.secondary_paxfwhm,2)
        secondary_struct(tidx).session_id = z.secondary_paxfwhm(sidx).session_id;
        secondary_struct(tidx).thickness = z.secondary_paxfwhm(sidx).thickness;
        secondary_struct(tidx).heatdata = z.secondary_paxfwhm(sidx).heatdata;
        tidx = tidx+1;
    end

end
%
for midx = 1: height(tmp.secondary_list)
    tmp.primary_matchloc = matches(tmp.primary_list.folder,tmp.secondary_list.folder{midx});
    
    tmp.primary_dirinfo = tmp.primary_list(tmp.primary_matchloc,:);
    % load primarystruct
    tmp.paxprimarypath = fullfile(tmp.primary_dirinfo.folder{1},tmp.primary_dirinfo.name{1});
    tmp.loadstruct = load(tmp.paxprimarypath);
    primarystruct = tmp.loadstruct.primarystruct;
    % find matching session
    tmp.primary_sessionidlist = string({primarystruct.sessionid});
    tmp.p2s_sessionidloc = find(matches(tmp.primary_sessionidlist,[secondary_struct.session_id]));
    for psidx = tmp.p2s_sessionidloc
        sidx = find(matches([secondary_struct.session_id],primarystruct(psidx).sessionid));
        secondary_struct(sidx).infodict  = primarystruct(psidx).infodict;
        secondary_struct(sidx).analog  = primarystruct(psidx).analog;
        secondary_struct(sidx).fwhmline  = primarystruct(psidx).fwhmline;
    end
end

end

