%ANALYSIS_MAIN_PIVJUDGE  One mouse's sessions judged for PIV by eye, a marker file per session
%   Needs no FWHM. A short stretch of both stacks, preprocessed as analysis_main does, under the pax line dialog,
%   then y/N : pivquality_pass.mat or pivquality_fail.mat in primary_analysis. Sessions holding either are skipped.
clc, clear
addpath('g:/03_program/01_ecspress/09_dirstruct');    % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

%% settings
dirs.mouse = 'K:\tmp\02_egfptdt\hql066';    % char  one mouse : every session folder under it
param.channel_lumen = 'ch2';                 % the vessel channel of this mouse (char)
param.channel_parenchyma   = 'ch1';                 % the tissue channel of this mouse (char)
param.n_page        = 1000;                  % pages read from the middle of each recording (count)

%% 1. Sessions : every XY movie under the mouse with no verdict yet; a z-stack holds stacks too but no frame rate
info_files = dir(fullfile(dirs.mouse, '**', '*_info.txt'));
info_paths = fullfile({info_files.folder}, {info_files.name});
is_xymovie = cellfun(@(p) contains(fileread(p), "scanmode,XY Movie"), info_paths);
session_dirs = unique(string({info_files(is_xymovie).folder}));
pass_files = fullfile(session_dirs, 'primary_analysis', 'pivquality_pass.mat');
fail_files = fullfile(session_dirs, 'primary_analysis', 'pivquality_fail.mat');
is_judged = isfile(pass_files) | isfile(fail_files);
fprintf('%s : %d sessions, %d judged, %d to go\n', dirs.mouse, numel(session_dirs), nnz(is_judged), nnz(~is_judged));
session_dirs = session_dirs(~is_judged);

%% 2. Judge, one session at a time : the pax line on a short stretch, then the verdict
for k = 1:numel(session_dirs)
    session = ECSSession(char(session_dirs(k)));
    session.roilist = roi_handle(fullfile(session.dir_struct.primary_analysis, 'roilist.mat'));   % loads the file if there is one; FWHM and cluster are not read
    % the middle of the recording, the stretch the person looks at
    n_page = numel(session.img_param.taxis);
    first_page = max(1, floor((n_page - param.n_page) / 2) + 1);
    page_range = [first_page, min(n_page, first_page + param.n_page - 1)];
    session.stackch1 = session.loadstack('ch1', page_range);
    session.stackch2 = session.loadstack('ch2', page_range);
    twophoton_processed = twophoton_preprocess(session);
    % the pax dialog as a viewer only : the line is not saved, analysis_main checks it against its FWHM output
    session.roilist.addormodifyroi(twophoton_processed.(param.channel_parenchyma), 'pax', 'line', ...
        twophoton_processed.(param.channel_lumen));
    verdict_reply = input(sprintf('[%d/%d] %s : pass to PIV? [y/N]: ', k, numel(session_dirs), session_dirs(k)), 's');
    pivquality = struct();
    pivquality.pass       = strcmpi(verdict_reply, 'y');   % logical  the person's verdict, also the file name
    pivquality.judged_at  = string(datetime('now'));       % str      when
    pivquality.page_range = page_range;                    % 1 x 2 int  the stretch that was shown, raw pages
    marker_name = 'pivquality_fail.mat';
    if pivquality.pass
        marker_name = 'pivquality_pass.mat';
    end
    marker_path = fullfile(session.dir_struct.primary_analysis, marker_name);
    save(marker_path, 'pivquality');
    fprintf('-> %s\n', marker_path);
end
