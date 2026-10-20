%E_TABLEGENERATION  The centralized state files -> the tables the figure scripts read.
%   Run d_centralize first : every function here reads a centralized_*_state.mat. Nothing is
%   read from or written to the data tree; the tables land in the dataset folder.
%   paxfwhm       centralized_paxfwhm_state -> state_summary, transition
%   drift         centralized_drift_state, sessions in centralized_paxfwhm_state -> drift_summary, drift_transition
%   fwhmrelation  centralized_paxfwhm_state -> fwhmrelation
%   Every row filter is set here and nowhere else. The figure scripts only choose what to
%   draw from these tables, and each function's quality gates stay in that function.
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

%% settings
param.dataset = "merged";   % str   the dataset folder the tables are written to
resolved = dirs_central();
dirs.central = resolved.central;                                % 1 x 1 str  the centralized_*.mat files
dirs.out = fullfile(resolved.secondary_root, param.dataset);    % 1 x 1 str  where the tables go
fprintf('central %s\nout     %s\n', dirs.central, dirs.out);

%% filters : which sessions the tables are about
% the state tables, paxfwhm and drift share them
state_filt.vessel_prefix = "PA";                                 % str   arteries only
state_filt.depth_thr     = 70;                                   % um    L1 / L2-3 boundary
state_filt.min_bout_s    = 30;                                   % s     shorter is not a state
state_filt.state_names   = ["awake", "drowsy", "nrem", "rem"];   % str   the scored states kept
drift_filt.require_fwhm  = true;                                 % bool  only recordings the FWHM pipeline has judged
% the heatmap tables
fwhmrelation_filt.vessel_prefix = "PA";      % str    arteries only
fwhmrelation_filt.depth_thr     = 70;        % um     cortical depth kept, shallower than this
fwhmrelation_filt.session_type  = "sleep";   % str    session type kept, "" = every type
fwhmrelation_filt.drop_flagged  = true;      % bool   drop a folder whose name says notanalyzable
fwhmrelation_filt.max_pxsize    = 0.3;       % um/px  coarsest pixel size kept, Inf = every size

%% tables : each reads its own centralized state file, none reads another's output
info.paxfwhm = tablegeneration_paxfwhm(dirs.central, dirs.out, ...
    state_filt.vessel_prefix, state_filt.depth_thr, state_filt.min_bout_s, state_filt.state_names);
info.drift = tablegeneration_drift(dirs.central, dirs.out, ...
    state_filt.vessel_prefix, state_filt.depth_thr, state_filt.min_bout_s, state_filt.state_names, ...
    drift_filt.require_fwhm);
info.fwhmrelation = tablegeneration_fwhmrelation(dirs.central, dirs.out, ...
    fwhmrelation_filt.vessel_prefix, fwhmrelation_filt.depth_thr, fwhmrelation_filt.session_type, ...
    fwhmrelation_filt.drop_flagged, fwhmrelation_filt.max_pxsize);
