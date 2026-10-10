function [sheet_path, n_scanned] = dirs_scantable(working_root, secondary_root, cohort)
%DIRS_SCANTABLE  Scan one cohort under every root and write its own dirtable.
%   Caller: dirtable_main
%   READ ONLY under working_root. The only writes are the cohort workbook and its MDF
%   metadata cache, both under <secondary_root>\<cohort>. The cache is keyed by the
%   session's full path, so one file serves every root. dirs_mergetable joins the sheets
%   this writes; a cohort that stands alone is scanned and never joined.
%
% IN   working_root    1 x K str     the mounted data trees, each one folder per cohort
%      secondary_root  1 x 1 str     where the sheet goes. NOT the data
%      cohort          1 x 1 str     the cohort folder name
% OUT  sheet_path      1 x 1 str     <secondary_root>\<cohort>\<cohort>_dirtable.xlsx
%      n_scanned       1 x 1 double  sessions the scan found, every root together
arguments
    working_root (1, :) string
    secondary_root (1, 1) string
    cohort (1, 1) string
end

cohort_dir = fullfile(working_root, cohort);        % 1xK, one per root
cohort_out = fullfile(secondary_root, cohort);
has_cohort = isfolder(cohort_dir);
if ~any(has_cohort)
    error('dirs_scantable:noCohort', 'cohort folder %s under none of: %s', ...
        cohort, strjoin(working_root, ', '));
end
if ~isfolder(cohort_out)
    mkdir(cohort_out)
end
sheet_path = fullfile(cohort_out, cohort + "_dirtable.xlsx");
cache_path = fullfile(cohort_out, cohort + "_dircache.mat");

scan_start = tic;
scanned = cell(1, numel(cohort_dir));
for root_idx = find(has_cohort)
    scanned{root_idx} = dirs_mapstruct(char(cohort_dir(root_idx)), char(cache_path));
end
scanned_table = vertcat(scanned{has_cohort});
repeated = rootoverlap(string(scanned_table.Directory));
if ~isempty(repeated)
    error('dirs_scantable:overlap', 'session folders under more than one root:\n%s', ...
        strjoin(repeated, newline));
end
dirs_writetable(scanned_table, char(sheet_path));
n_scanned = height(scanned_table);
fprintf('%-12s scanned %d sessions under %d root(s) in %.1f s\n', cohort, ...
    n_scanned, nnz(has_cohort), toc(scan_start));
end

%% ------------------------------------------------------------------------
function repeated = rootoverlap(directory)
%ROOTOVERLAP  Session folders the scan found under more than one root, the letter removed.
%   Caller: dirs_scantable
%
% IN   directory  N x 1 str   Directory of every scanned row, all roots together
% OUT  repeated   M x 1 str   the rows whose path below the letter appeared already
below_letter = lower(extractAfter(directory, ':'));
[~, first_idx] = unique(below_letter, 'stable');
is_repeat = true(numel(directory), 1);
is_repeat(first_idx) = false;
repeated = directory(is_repeat);
end
