function stats = sweep_empty_figuredirs(root, opts)
%SWEEP_EMPTY_FIGUREDIRS  Remove empty figure folders under every primary_analysis below ROOT.
%
%   A candidate figure folder is removed ONLY if its entire tree contains no
%   files (only empty subfolders). No file is ever deleted, so this is safe to
%   run over a whole data drive.
%
%   stats = sweep_empty_figuredirs(root)                 % DRY RUN (default): reports only
%   stats = sweep_empty_figuredirs(root, dryrun=false)   % actually delete
%
%   Targets both the legacy nested run folders (figures_*) and the new flat
%   category folders (*_fwhm, *_radon_figures, *_polarcluster, *_roi_fig).
%
%   Output stats has fields: dryrun, primary_analysis_scanned, empty_found,
%   nonempty_kept, list (cellstr of the empty folders found/removed).
arguments
    root (1,:) char
    opts.dryrun (1,1) logical = true
end

patterns = {'figures_*', '*_fwhm', '*_radon_figures', '*_polarcluster', '*_roi_fig'};

% Find every primary_analysis folder under root
pa = dir(fullfile(root, '**', 'primary_analysis'));
pa = pa([pa.isdir]);
pa_paths = unique(fullfile({pa.folder}, {pa.name}));

removed = {};
nonempty_kept = 0;
for p = 1:numel(pa_paths)
    for k = 1:numel(patterns)
        cand = dir(fullfile(pa_paths{p}, patterns{k}));
        cand = cand([cand.isdir]);
        for c = 1:numel(cand)
            if strcmp(cand(c).name, '.') || strcmp(cand(c).name, '..'), continue; end
            folder = fullfile(cand(c).folder, cand(c).name);
            if local_is_tree_empty(folder)
                removed{end+1} = folder; %#ok<AGROW>
                if ~opts.dryrun
                    rmdir(folder, 's');   % safe: verified no files in the tree
                end
            else
                nonempty_kept = nonempty_kept + 1;
            end
        end
    end
end

stats.dryrun = opts.dryrun;
stats.primary_analysis_scanned = numel(pa_paths);
stats.empty_found = numel(removed);
stats.nonempty_kept = nonempty_kept;
stats.list = removed(:);

tag = '[DRY RUN] would remove';
if ~opts.dryrun, tag = '[DELETED] removed'; end
fprintf('%s %d EMPTY figure folders; kept %d non-empty; scanned %d primary_analysis folders under %s\n', ...
    tag, stats.empty_found, stats.nonempty_kept, stats.primary_analysis_scanned, root);
end

function tf = local_is_tree_empty(folder)
    % Empty == no files anywhere in the tree (subfolders may exist but must be empty too).
    entries = dir(fullfile(folder, '**', '*'));
    tf = ~any(~[entries.isdir]);
end
