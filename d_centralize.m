%D_CENTRALIZE  Every session's products into the centralized files, then each series split by its bouts.
%   Run dirtable_main first : cell 1 reads the sheet it writes. Nothing is read from or written to
%   the data tree; everything lands under dirs.central, one file per product.
%   1  gather    the sheet's sessions -> centralized_{paxfwhm, polarcluster, roilist, sleep_score, analysis_analog, drift}
%   2  by bout   centralized_{sleep_score, analysis_analog} -> centralized_bout, then
%                centralized_{paxfwhm, drift} x centralized_bout -> centralized_{paxfwhm_state, drift_state}
clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

%% settings
param.dataset = "merged";   % str   the dirtable whose reference sheet lists the sessions
param.rebuild = false;      % bool  true recomputes every row of every file
resolved = dirs_central();
dirs.central = resolved.central;                                                                    % 1 x 1 str  the centralized_*.mat files
dirs.sheet = fullfile(resolved.secondary_root, param.dataset, param.dataset + "_dirtable.xlsx");   % 1 x 1 str  dirtable_main's workbook
if ~isfolder(dirs.central)
    mkdir(dirs.central);
end
fprintf('sheet   %s\ncentral %s\n', dirs.sheet, dirs.central);

%% 1. gather : one centralized file per session product, read from the session folders the sheet lists
info.mat = centralize_mat(dirs.central, dirs.sheet, param.rebuild);
info.drift = centralize_drift(dirs.central, dirs.sheet, param.rebuild);

%% 2. by bout : the scorings as bout tables, then every series over them
info.bout = centralize_bout(dirs.central, param.rebuild);
info.paxfwhm_state = bout_paxfwhm(dirs.central, param.rebuild);
info.drift_state = bout_drift(dirs.central, param.rebuild);
