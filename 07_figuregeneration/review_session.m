% Path setup
clc, clear
addpath("G:\03_program\01_ecspress\09_dirstruct"); dirs_ecspath;

%close all
% Directory setup
% sessiondir = 'G:\tmp\01_igkltdt\hql104\260607_hql104_sleep\HQL104_sleep260607_006';
sessiondir ='G:\tmp\00_igkl\hql090\251016_hql090_sleep\HQL090_sleep251016_005';

% 1. Load data & 3. Load processed data (Integrated via ECSSession)
session = ECSSession(sessiondir);
session = session.load_primary_results;

if isfile(fullfile(sessiondir,"peripheral/sleep_score.mat"))
    disp('Loading sleepscore file')
    sleepscore = load(fullfile(sessiondir,"peripheral/sleep_score.mat"));
else
    disp('no sleep_socre.mat')
end

% what this session has been graded as, if it has. Read before the figures so the
% previous verdict is on screen while you look, not after you have formed a new one
ledger_quality(sessiondir);

%%
pax_fig = analysis_pax_makefig(session.pax_fwhm, session.pax_fwhm.t_axis,...
    session.img_param.pixel2um, session.dir_struct.figures_fwhm);

%% the cluster panels. No stack is read: load_primary_results already brought
% polarcluster.mat in, and t_axis and pixel2um stand in for twophoton_processed the
% same way the pax call above does it
analysis_cluster_makefig(session.polarcluster, session.roilist, session.pax_fwhm,...
    session.pax_fwhm.t_axis, session.img_param.pixel2um,...
    session.dir_struct.figures_polarcluster);

%% the polar profiles. Only for a session that went through
% analysis_clusterpolar.transform -- this reads pax_angle, polar_theta,
% polar_profiles and polar_binstart_deg, and a polarcluster.mat saved before that
% step carries none of the four
analysis_polar_makefig(session.polarcluster, session.dir_struct.figures_polarcluster);

%% the awake bouts, for a session with no sleep score. Puff x locomotion from the ball and
% the airtable in peripheral/analysis_analog.mat, which peripheral_batch writes
if ~exist('sleepscore', 'var')
    primary_analog = analysis_analog.load_object(fullfile(sessiondir, "peripheral"));
    awake_integrate = state_integration([], primary_analog);
    awake_fig = analysis_awake_makefig(awake_integrate.time_table, session.pax_fwhm.t_axis, ...
        session.pax_fwhm.thickness.bv, session.pax_fwhm.thickness.totalpvs, ...
        session.img_param.pixel2um, primary_analog.ball.rs_taxis, primary_analog.ball.resampled_ball, ...
        session.dir_struct.figures_fwhm);
end
%%

%% Grade it. The note is in ENGLISH and says what the columns cannot see -- where the
% kymograph drifts, which boundary the FWHM loses, whether the vessel sits at an angle.
% bv_range and cond_sd are already in the table and do not belong in the note.
%   S  goes in a representative figure. Clean from start to end, boundaries steady
%   A  used in the analysis as it is. Local defects that do not move the measurement
%   B  used, never representative. The defect is visible and the note says where
%   C  dropped once more data arrives. Kept now because the sample is worth having
% ledger_quality(sessiondir, "A", "")
