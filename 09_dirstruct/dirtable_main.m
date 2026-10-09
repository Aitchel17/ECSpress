%DIRTABLE_MAIN  The directory tables: each cohort's own sheet, then the joined dataset.
%   Reads the data tree and writes each cohort's workbook and MDF metadata cache under
%   <secondary_root>\<cohort>, and the dataset workbook under <secondary_root>\<dataset>. A cohort in
%   param.alone_cohort gets its own sheet and is never joined into a dataset.
%   The only place the dataset's membership is decided: centralize_primary and centralize_drift
%   read the sheet this writes and hold no cohort list of their own.
clc
clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;

%% 0. Settings
resolved = dirs_central();
dirs.working_root   = resolved.working_root;     % 1 x K str  every mounted work and archive tmp, by beacon
dirs.secondary_root = resolved.secondary_root;   % 1 x 1 str  where the sheets go

param.dataset        = "merged";                                 % 1 x 1 str  the joined dataset
param.dataset_cohort = ["00_igkl", "01_igkltdt", "02_egfptdt"];  % 1 x N str  its cohorts, in the order they join
param.alone_cohort   = strings(1, 0);                            % 1 x M str  cohorts that keep their own sheet only; none today

fprintf('roots   %s\nsheets  %s\n', strjoin(dirs.working_root, ' | '), dirs.secondary_root);

%% 1. Each cohort's own sheet, scanned under every root
for cohort = [param.dataset_cohort, param.alone_cohort]
    dirs_scantable(dirs.working_root, dirs.secondary_root, cohort);
end

%% 2. The dataset, joined from its cohorts' sheets; alone_cohort is never joined
dirs_mergetable(dirs.working_root, dirs.secondary_root, param.dataset_cohort, param.dataset);
