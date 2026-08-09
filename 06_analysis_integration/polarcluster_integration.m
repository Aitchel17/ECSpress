%POLARCLUSTER_INTEGRATION  Collect every session's polar contour profiles into one file.
%   Walks a data root, rebuilds the cartesian->polar transform for every session
%   that carries the four manual contours, measures the max-dilation axis, and
%   writes polar_ave.mat. Nothing here is aligned, normalised or averaged --
%   polar_ave_fig does that, so a colour or normalisation change never re-runs
%   the edge detection.
%
%   why      the transform is the expensive half (edge detect x4 per session over
%            the whole set); alignment and normalisation are arithmetic on 360
%            numbers. Splitting there is what makes the figure script cheap to
%            iterate on
%   caution  almost no polarcluster.mat on disk carries polar_profiles, while
%            nearly every one carries all four contours in roilist.mat. The
%            contours are the manual work; the profiles are recomputed here and
%            NOT saved back -- see FINDINGS.md
%   note     nothing under the data root is written. The only output is
%            polar_ave.mat under dirs.save_dir

clc, clear
addpath(genpath('g:\03_program\01_ecspress'));

dirs.data_root = 'G:\tmp';
dirs.save_dir  = ['E:\OneDrive - The Pennsylvania State University\2023ecspress' ...
    '\02_secondary_analysis\polar_ave'];

% One place where every field of a session row is declared, so its shape is
% readable without tracing which line wrote what.
%   caution  index k of every r_*_px is the bin whose CENTER is (k-0.5) + bin_start
%            degrees. bin_start is the same for every session (-7.5, pax-independent
%            since analysis_clusterpolar_polarplot stopped baking PAX in), so the
%            offset is common and cancels once each session is rotated to its own
%            max-dilation axis
session_proto = struct( ...
    'session_id',       "",  ...   % string           session folder name
    'mouse_id',         "",  ...   % string           hql###
    'vessel_id',        "",  ...   % string           PA01 ... , first token of Comments
    'session_type',     "",  ...   % string           sleep | whiskerb | whisker
    'depth_um',         NaN, ...   % float um         from Comments, NaN if absent
    'pixel2um',         NaN, ...   % float um/px      from objpix, NaN if absent
    'r_cbv_px',         [],  ...   % 1 x 360 float px constricted BV, NaN = empty bin
    'r_cpvs_px',        [],  ...   % 1 x 360 float px constricted PVS
    'r_dbv_px',         [],  ...   % 1 x 360 float px dilated BV
    'r_dpvs_px',        [],  ...   % 1 x 360 float px dilated PVS
    'theta_maxdil_deg', NaN, ...   % float deg        dilation axis, 1st circ harmonic
    'theta_argmax_deg', NaN, ...   % float deg        smoothed peak, for comparison only
    'maxdil_px',        NaN, ...   % float px         peak of the smoothed dilation
    'dil_anisotropy',   NaN, ...   % float 0-1        0 = swells evenly, axis meaningless
    'n_empty_bin',      [],  ...   % 1 x 4 int        empty bins per profile, of 360
    'maxgap_deg',       [],  ...   % 1 x 4 int deg    longest contiguous empty run
    'folder',           "");       % string           primary_analysis folder

param.smooth_deg = 15;    % circular moving-average width for the alignment measurement

%% 1. Find every session that has both files
file_list = dir(fullfile(dirs.data_root, '**', 'polarcluster.mat'));
fprintf('polarcluster.mat found : %d under %s\n', numel(file_list), dirs.data_root);

sessions = session_proto([]);
skipped  = strings(0, 2);

for file_idx = 1:numel(file_list)
    primary_dir = file_list(file_idx).folder;
    roi_path    = fullfile(primary_dir, 'roilist.mat');
    session_dir = fileparts(primary_dir);
    [~, session_name] = fileparts(session_dir);

    if ~isfile(roi_path)
        skipped(end+1, :) = [string(session_name), "no roilist.mat"];  %#ok<SAGROW>
        continue
    end

    loaded_pc = load(fullfile(primary_dir, 'polarcluster.mat'));
    loaded_rl = load(roi_path);
    polarcluster = loaded_pc.polarcluster;
    roilist      = loaded_rl.roilist;

    % 1.1 the transform itself -- the SAME function the per-session figure uses,
    %     so a row here and that session's polar figure cannot drift apart
    polarcluster = analysis_clusterpolar_polarplot(polarcluster, roilist);
    if ~isfield(polarcluster, 'polar_profiles') || isempty(polarcluster.polar_profiles)
        skipped(end+1, :) = [string(session_name), "transform declined"];  %#ok<SAGROW>
        continue
    end

    % 1.2 unpack by label, not by position
    profile_r = nan(4, 360);
    want = {'constrictedBV_contour','constrictedPVS_contour', ...
            'dilatedBV_contour','dilatedPVS_contour'};
    labels_found = string({polarcluster.polar_profiles.roi_label});
    for want_idx = 1:4
        hit = find(labels_found == want{want_idx}, 1);
        if ~isempty(hit)
            profile_r(want_idx, :) = polarcluster.polar_profiles(hit).binned_r;
        end
    end

    % 1.3 the dilation axis. Gaps are filled FIRST : a contour is a closed curve,
    %     so an empty bin is a sampling gap, not a missing radius
    r_cbv_filled = fill_circular(profile_r(1, :));
    r_dbv_filled = fill_circular(profile_r(3, :));
    dilation_px  = r_dbv_filled - r_cbv_filled;
    [theta_maxdil, theta_argmax, maxdil, anisotropy] = dilation_axis(dilation_px, param.smooth_deg);

    % 1.4 provenance from the session's own _info.txt, the same source
    %     primary_integration and the dirtable read
    info = read_sessioninfo(session_dir);

    row = session_proto;
    row.session_id       = string(session_name);
    row.mouse_id         = info.mouse_id;
    row.vessel_id        = info.vessel_id;
    row.session_type     = info.session_type;
    row.depth_um         = info.depth_um;
    row.pixel2um         = info.pixel2um;
    row.r_cbv_px         = profile_r(1, :);
    row.r_cpvs_px        = profile_r(2, :);
    row.r_dbv_px         = profile_r(3, :);
    row.r_dpvs_px        = profile_r(4, :);
    row.theta_maxdil_deg = theta_maxdil;
    row.theta_argmax_deg = theta_argmax;
    row.maxdil_px        = maxdil;
    row.dil_anisotropy   = anisotropy;
    row.n_empty_bin      = sum(isnan(profile_r), 2)';
    row.maxgap_deg       = arrayfun(@(k) longest_gap(profile_r(k, :)), 1:4);
    row.folder           = string(primary_dir);
    sessions(end+1) = row;  %#ok<SAGROW>

    fprintf('%-44s %-8s %-6s %-9s theta %5.1f (argmax %5.1f) aniso %.2f  maxdil %5.1f px  gap %s\n', ...
        session_name, row.mouse_id, row.vessel_id, row.session_type, theta_maxdil, ...
        theta_argmax, anisotropy, maxdil, num2str(row.maxgap_deg, '%4d'));
end

fprintf('\ncollected %d sessions, skipped %d\n', numel(sessions), size(skipped, 1));
for skip_idx = 1:size(skipped, 1)
    fprintf('  skip %-44s %s\n', skipped(skip_idx, 1), skipped(skip_idx, 2));
end

%% 2. Save
polar_ave = struct();
polar_ave.sessions = sessions;
polar_ave.param    = struct( ...
    'data_root',    string(dirs.data_root), ...
    'smooth_deg',   param.smooth_deg, ...
    'n_bin',        360, ...
    'bin_start_deg', -7.5, ...     % see the caution on session_proto
    'built',        datetime('now'));

if ~exist(dirs.save_dir, 'dir')
    mkdir(dirs.save_dir);
end
save(fullfile(dirs.save_dir, 'polar_ave.mat'), 'polar_ave');
fprintf('\nwrote %s\n', fullfile(dirs.save_dir, 'polar_ave.mat'));

% a table view for eyeballing, rows side by side
summary_table = struct2table(rmfield(sessions, ...
    {'r_cbv_px','r_cpvs_px','r_dbv_px','r_dpvs_px','folder'}));
disp(summary_table)


%% ── Local functions ────────────────────────────────────────────────────

function filled = fill_circular(r_deg)
%FILL_CIRCULAR  Interpolate empty angular bins around the circle.
%   why  the contour is closed, so a radius exists at every angle -- an empty bin
%        means no edge pixel landed in it, not that the wall is absent
%   err  a plain interp1 breaks at the wrap, leaving the 359->0 gap unfilled
    filled = r_deg;
    valid  = ~isnan(r_deg);
    if nnz(valid) < 2
        return
    end
    theta = 0:359;
    % tile three turns so the wrap is interior to the interpolation
    theta_tiled = [theta - 360, theta, theta + 360];
    r_tiled     = [r_deg, r_deg, r_deg];
    valid_tiled = ~isnan(r_tiled);
    filled = interp1(theta_tiled(valid_tiled), r_tiled(valid_tiled), theta, 'linear');
end

function [theta_center, theta_argmax, peak, anisotropy] = dilation_axis(dilation_px, smooth_deg)
%DILATION_AXIS  Which way the vessel swells, as the first circular harmonic.
%   why      a single argmax bin is one number off a profile that was two thirds
%            interpolated. The harmonic weights the whole lobe instead
%   err      the first try weighted max(dilation,0), which keeps the DC term : a
%            vessel that swells everywhere then hands back the centroid of the
%            whole ring, which on a real session put the axis most of a quadrant
%            off the visible lobe (see FINDINGS.md). sum(dil .* exp(i*theta))
%            has no DC
%            term (sum of exp(i*theta) over a full circle is 0), so only the
%            ANISOTROPIC part steers it
%   note     anisotropy near 0 means the dilation was radially even and the axis
%            is arbitrary -- the caller decides whether such a session belongs in
%            an aligned average
    smoothed = smooth_circular(dilation_px, smooth_deg);
    theta    = deg2rad(0:359);
    if all(isnan(smoothed))
        theta_center = NaN;  theta_argmax = NaN;  peak = NaN;  anisotropy = NaN;  return
    end
    smoothed(isnan(smoothed)) = 0;
    harmonic     = sum(smoothed .* exp(1i * theta));
    theta_center = mod(rad2deg(angle(harmonic)), 360);
    anisotropy   = abs(harmonic) / sum(abs(smoothed));
    [peak, argmax_idx] = max(smoothed);
    theta_argmax = argmax_idx - 1;
end

function smoothed = smooth_circular(r_deg, width_deg)
%SMOOTH_CIRCULAR  Moving average that wraps at 0/360.
    tiled    = [r_deg, r_deg, r_deg];
    tiled    = movmean(tiled, width_deg, 'omitnan');
    smoothed = tiled(361:720);
end

function gap = longest_gap(r_deg)
%LONGEST_GAP  Longest run of empty bins, wrapping at 0/360 (degrees).
    empty_bin = isnan(r_deg);
    if all(empty_bin)
        gap = 360;
        return
    end
    if ~any(empty_bin)
        gap = 0;
        return
    end
    tiled = [empty_bin, empty_bin];
    gap   = 0;  run = 0;
    for k = 1:numel(tiled)
        if tiled(k); run = run + 1; gap = max(gap, run); else; run = 0; end
    end
    gap = min(gap, 360);
end

function info = read_sessioninfo(session_dir)
%READ_SESSIONINFO  Vessel / depth / resolution from the session's _info.txt.
%   why  this is the same file primary_integration and the dirtable read, so a
%        vessel id here matches the one the transition figures counted
    info = struct('mouse_id', "", 'vessel_id', "", 'session_type', "", ...
        'depth_um', NaN, 'pixel2um', NaN);

    % .../<mouse>/<date>_<mouse>_<type>/<session>
    day_dir = fileparts(session_dir);
    [mouse_dir, day_name] = fileparts(day_dir);
    [~, mouse_name] = fileparts(mouse_dir);
    info.mouse_id = string(mouse_name);
    day_parts = split(string(day_name), '_');
    if numel(day_parts) >= 3
        info.session_type = day_parts(end);
    end

    info_file = dir(fullfile(session_dir, '*_info.txt'));
    if isempty(info_file)
        return
    end
    info_table = readtable(fullfile(info_file(1).folder, info_file(1).name));
    fields = string(info_table.Field);
    values = string(info_table.Value);

    if any(fields == "Comments")
        comment = values(fields == "Comments");
        comment = comment(1);
        tokens  = split(strtrim(comment), ' ');
        % err  the first token carries whatever punctuation the comment had --
        %      "PA01," and "pv01" both showed up, and either breaks a string
        %      compare against the dirtable's VesselID
        info.vessel_id = upper(regexprep(tokens(1), '[^A-Za-z0-9]', ''));
        depth_txt = regexp(comment, '(\d+)\s*um', 'tokens', 'once', 'ignorecase');
        if ~isempty(depth_txt)
            info.depth_um = str2double(depth_txt{1});
        end
    end
    if any(fields == "objpix")
        objpix = values(fields == "objpix");
        info.pixel2um = str2double(regexprep(objpix(1), '[^0-9.]', ''));
    end
end
