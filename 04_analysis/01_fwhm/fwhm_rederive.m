function report = fwhm_rederive(root, opt)
%FWHM_REDERIVE  Recompute thickness and displacement from the boundaries already stored.
%   Walks a tree of paxfwhm.mat, finds the ones whose thickness struct is missing a
%   field the current getdiameter produces, and rebuilds thickness and displacement
%   from obj.idx. Nothing is re-measured: the boundaries the kymograph analysis
%   found are read straight off the saved object, so there is no manual step and no
%   parameter to choose.
%
%   why  fields have been added and renamed over time. Some stored objects still
%        carry 'ecschanges_residual' and no 'eps', the naming from before
%        getdiameter gained eps/epschanges. see FINDINGS.md for the count
%   err  one such session is enough to break a merge. struct2table wraps EVERY row
%        in a cell when one row's field set differs, and explode_nest then fails on
%        a cell where it wants a struct array (line 14). That is how one 2025 file
%        stopped a 41-session table from building
%   caution  thickness and displacement are REPLACED. kymograph, idx, mask and param
%            are untouched, and a .bak copy is written first unless backup = false
%
% IN   root      char/str     tree to walk, e.g. 'G:\tmp\00_igkl'
%      require   1 x K str    thickness fields that must exist; a file missing any
%                             of them is rebuilt. Default "eps"
%      dryrun    bool         list what would change and write nothing. Default true
%      backup    bool         copy to paxfwhm_backup_<yyyymmdd>.mat first. Default true
%      stamp     char         the date that backup name carries. Default today
% OUT  report    table        one row per candidate: Directory, had, rebuilt, error
%
% USAGE
%   fwhm_rederive('G:\tmp\00_igkl')                     % see what would change
%   fwhm_rederive('G:\tmp\00_igkl', dryrun = false)     % do it

arguments
    root        (1,:) char
    opt.require (1,:) string = "eps"
    opt.dryrun  (1,1) logical = true
    opt.backup  (1,1) logical = true
    opt.stamp   (1,:) char = datestr(now, 'yyyymmdd') %#ok<TNOW1,DATST>
end

files = dir(fullfile(root, '**', 'primary_analysis', 'paxfwhm.mat'));
fprintf('%s : %d paxfwhm.mat\n', root, numel(files));

Directory = strings(0,1);   had = strings(0,1);
rebuilt   = false(0,1);     errmsg = strings(0,1);

for k = 1:numel(files)
    p = fullfile(files(k).folder, files(k).name);
    loaded = load(p, 'line_fwhm');
    if ~isfield(loaded, 'line_fwhm')
        continue
    end
    obj = loaded.line_fwhm;

    % 1. Does it already carry everything asked for
    present = string(fieldnames(obj.thickness))';
    if all(ismember(opt.require, present))
        continue
    end

    Directory(end+1,1) = string(files(k).folder); %#ok<AGROW>
    had(end+1,1)       = strjoin(present, ',');   %#ok<AGROW>
    rebuilt(end+1,1)   = false;                   %#ok<AGROW>
    errmsg(end+1,1)    = "";                      %#ok<AGROW>
    fprintf('  %-56s missing %s\n', ...
        extractBefore(extractAfter(string(files(k).folder), string(root) + filesep), ...
                      filesep + "primary_analysis"), ...
        strjoin(setdiff(opt.require, present), ','));
    if opt.dryrun
        continue
    end

    % 2. Rebuild. getdiameter reads obj.idx.clean_* only, getdisplacement the same,
    %    so both are pure re-derivations of what is already on disk
    try
        if opt.backup
            copyfile(p, fullfile(files(k).folder, ...
                sprintf('paxfwhm_backup_%s.mat', opt.stamp)));
        end
        obj.getdiameter();
        obj.getdisplacement();
        obj.save2disk('paxfwhm', files(k).folder);
        rebuilt(end,1) = true;
    catch ME
        errmsg(end,1) = string(ME.message);
        fprintf('    FAILED : %s\n', ME.message);
    end
end

report = table(Directory, had, rebuilt, errmsg);
tail = '';
if opt.dryrun
    tail = '  (dryrun, nothing written)';
end
fprintf('%d candidate(s), %d rebuilt%s\n', height(report), nnz(rebuilt), tail);
end
