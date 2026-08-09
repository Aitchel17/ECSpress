%POLAR_AVE_FIG  Vessel-averaged polar contours, aligned on the dilation axis.
%   Reads polar_ave.mat (polarcluster_integration), rotates every session so its
%   own max-dilation direction sits at 0 deg, divides by that session's mean
%   constricted-BV radius, averages sessions within a vessel and then draws every
%   vessel faint with the across-vessel mean bold.
%
%   why      alignment and normalisation are arithmetic on 360 numbers, so they
%            live here : changing either re-runs in a second, where the transform
%            in polarcluster_integration takes minutes
%   caution  the radial unit is "x mean constricted BV radius". Constricted BV
%            therefore sits near 1 BY CONSTRUCTION -- it is a reference circle,
%            not a result
%
%   NOTHING IS INTERPOLATED HERE. Empty angular bins stay NaN, every mean is
%   omitnan, and both the per-vessel lines and the mean break where the contour
%   was not measured.
%     why      a contour at radius r px is drawn from ~2*pi*r pixels, so 360 bins
%              oversamples it -- interpolation only resamples, it does not add
%              information, and it hid how thin the support was
%     see FINDINGS.md
%     caution  rotation is a whole-degree circshift for the same reason. theta is
%              a circular centroid and lands between bins, so interpolating to it
%              would smear NaNs into measured bins. The cost is <= 0.5 deg

clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
clee = color_lee();

dirs.save_dir = ['E:\OneDrive - The Pennsylvania State University\2023ecspress' ...
    '\02_secondary_analysis\polar_ave'];
load_struct = load(fullfile(dirs.save_dir, 'polar_ave.mat'));
sessions = load_struct.polar_ave.sessions;

param.session_type   = "sleep";   % "" = every type
param.vessel_prefix  = "PA";      % arteries only, the same filter tablegeneration_main uses
param.maxgap_deg     = 120;       % drop a session whose contour has a longer empty arc
param.min_anisotropy = 0.10;      % below this the dilation was radially even -- see below
param.min_vessel_n   = 3;         % blank the mean at angles fewer vessels than this reach
param.sector_half    = 30;        % half width of the reporting sectors, deg
param.fontsize       = 14;        % matched to transition_fig / transition_prepost

% Angular bin width. Pooling MEASURED 1-deg bins, still no interpolation.
%   why      1 deg oversamples the pixel grid : a contour at r px has only
%            ~2*pi*r pixels to spread over 360 bins, so most stay empty
%   see FINDINGS.md
%   note     the answer does not move across that range -- only the fragmentation
%            of the drawn line does. 360 must divide by this
param.bin_deg = 5;

% Same hues as the transition figures so the two sets read as one story
figconfig.bv.faint  = clee.lch(80, 130,  10);
figconfig.bv.bold   = clee.lch(45, 150,  15);
figconfig.pvs.faint = clee.lch(80, 130, 110);
figconfig.pvs.bold  = clee.lch(45, 150, 115);

%% 1. Select
% caution  aligning a session whose dilation was radially even puts an arbitrary
%          direction at 0 deg, which does not average out -- it flattens the mean
%          lobe of every OTHER vessel. That is what param.min_anisotropy is for
keep = true(1, numel(sessions));
if param.session_type ~= ""
    keep = keep & [sessions.session_type] == param.session_type;
end
if param.vessel_prefix ~= ""
    keep = keep & startsWith([sessions.vessel_id], param.vessel_prefix);
end
keep = keep & arrayfun(@(s) max(s.maxgap_deg) <= param.maxgap_deg, sessions);
keep = keep & ~isnan([sessions.theta_maxdil_deg]);
keep = keep & [sessions.dil_anisotropy] >= param.min_anisotropy;
selected = sessions(keep);
fprintf('%d of %d sessions kept (type=%s, %s*, maxgap<=%d deg, aniso>=%.2f)\n', ...
    numel(selected), numel(sessions), param.session_type, param.vessel_prefix, param.maxgap_deg, param.min_anisotropy);

%% 2. Align + normalise, one row per session, gaps left as NaN
n_bin = 360 / param.bin_deg;
theta_deg = (0:n_bin-1) * param.bin_deg + param.bin_deg/2;      % bin centres
n_sel = numel(selected);
profile_fields = ["r_cbv", "r_cpvs", "r_dbv", "r_dpvs"];
source_fields  = ["r_cbv_px", "r_cpvs_px", "r_dbv_px", "r_dpvs_px"];
aligned = struct();
for field_idx = 1:4
    aligned.(profile_fields(field_idx)) = nan(n_sel, n_bin);
end
for sel_idx = 1:n_sel
    s = selected(sel_idx);
    % the divisor is the mean over MEASURED bins only
    scale_px = mean(s.r_cbv_px, 'omitnan');
    for field_idx = 1:4
        rotated = rotate_bins(s.(source_fields(field_idx)), s.theta_maxdil_deg);
        aligned.(profile_fields(field_idx))(sel_idx, :) = ...
            mean(reshape(rotated, param.bin_deg, n_bin), 1, 'omitnan') / scale_px;
    end
end

%% 3. Sessions -> vessels
% why  the same vessel imaged on two days is one sample, not two. This is the
%      unit the transition figures count (Date -> VesselID), so n matches
vessel_key = [selected.mouse_id] + "_" + [selected.vessel_id];
[vessel_names, ~, vessel_group] = unique(vessel_key, 'stable');
n_vessel = numel(vessel_names);
byvessel = struct();
for field_idx = 1:4
    f = profile_fields(field_idx);
    byvessel.(f) = nan(n_vessel, n_bin);
    for vessel_idx = 1:n_vessel
        byvessel.(f)(vessel_idx, :) = ...
            mean(aligned.(f)(vessel_group == vessel_idx, :), 1, 'omitnan');
    end
end
fprintf('%d vessels from %d mice\n', n_vessel, numel(unique([selected.mouse_id])));
for field_idx = 1:4
    f = profile_fields(field_idx);
    per_angle_n = sum(~isnan(byvessel.(f)), 1);
    fprintf('  %-6s vessels reaching an angle : median %d, %.0f%% of angles have >= %d\n', ...
        f, median(per_angle_n), 100*mean(per_angle_n >= param.min_vessel_n), param.min_vessel_n);
end

%% 4. Figures
draw_polar_pair('polar_ave_BV', 'Blood vessel', theta_deg, ...
    byvessel.r_cbv, byvessel.r_dbv, figconfig.bv, param.min_vessel_n, param.fontsize, dirs.save_dir);
draw_polar_pair('polar_ave_PVS', 'Perivascular space', theta_deg, ...
    byvessel.r_cpvs, byvessel.r_dpvs, figconfig.pvs, param.min_vessel_n, param.fontsize, dirs.save_dir);

% caution  a change is NaN unless BOTH walls were measured in that same 1-deg bin,
%          so this figure is the thinnest-supported one of the three
change_bv  = byvessel.r_dbv  - byvessel.r_cbv;
change_pvs = byvessel.r_dpvs - byvessel.r_cpvs;
draw_polar_change('polar_ave_dilation', theta_deg, {change_bv, change_pvs}, ...
    {figconfig.bv, figconfig.pvs}, ["BV", "PVS"], param.min_vessel_n, param.fontsize, dirs.save_dir);

report_sector('BV ', change_bv,  theta_deg, param.sector_half);
report_sector('PVS', change_pvs, theta_deg, param.sector_half);


%% ── Local functions ────────────────────────────────────────────────────

function rotated = rotate_bins(r_deg, theta0_deg)
%ROTATE_BINS  Put theta0 at index 1 by whole-degree circshift.
%   why  circshift moves measured and empty bins together. interp1 to a fractional
%        theta would blend a NaN neighbour into every measured bin next to a gap,
%        turning a 30-bin gap into a 32-bin one
    rotated = circshift(r_deg, -round(theta0_deg));
end

function [mean_r, n_r] = mean_supported(r_matrix, min_n)
%MEAN_SUPPORTED  Column mean over measured entries, blank where too few reach.
    n_r = sum(~isnan(r_matrix), 1);
    mean_r = mean(r_matrix, 1, 'omitnan');
    mean_r(n_r < min_n) = NaN;
end

function draw_polar_pair(fig_name, title_text, theta_deg, r_constricted, r_dilated, cfg, min_n, fontsize, save_dir)
%DRAW_POLAR_PAIR  Every vessel faint, the across-vessel mean bold.
%   note  NaN bins break polarplot's line, which is the point : a gap in the
%         drawing is a gap in the data, not a straight segment across it
    fig = make_fig(fig_name, 'polar');
    fig.bring_fig();
    fig.update_figsize([4.5 4.5]);
    ax = fig.ax;
    hold(ax, 'on');
    theta = deg2rad([theta_deg, theta_deg(1)]);
    closed = @(r) [r, r(1)];

    for vessel_idx = 1:size(r_dilated, 1)
        polarplot(ax, theta, closed(r_constricted(vessel_idx, :)), ...
            'Color', cfg.faint, 'LineWidth', 0.5, 'LineStyle', ':');
        polarplot(ax, theta, closed(r_dilated(vessel_idx, :)), ...
            'Color', cfg.faint, 'LineWidth', 0.5);
    end
    [mean_constricted, n_constricted] = mean_supported(r_constricted, min_n);
    [mean_dilated,     n_dilated]     = mean_supported(r_dilated,     min_n);
    polarplot(ax, theta, closed(mean_constricted), ...
        'Color', cfg.bold, 'LineWidth', 2, 'LineStyle', ':');
    polarplot(ax, theta, closed(mean_dilated), ...
        'Color', cfg.bold, 'LineWidth', 2.5);

    % 0 deg is the max-dilation direction by construction -- mark it
    polarplot(ax, [0 0], [0 ax.RLim(2)], 'k--', 'LineWidth', 0.8);
    title(ax, sprintf('%s  (n = %d vessels)', title_text, size(r_dilated, 1)));
    subtitle(ax, sprintf('n \\geq %d at %.0f%% / %.0f%% of angles (constr. / dil.)', ...
        min_n, 100*mean(n_constricted >= min_n), 100*mean(n_dilated >= min_n)));
    ax.ThetaTick = 0:45:315;
    set(ax, 'FontName', 'Arial', 'FontSize', fontsize, 'LineWidth', 1);
    % caution  make_fig's polaraxes fills the figure, so a title plus subtitle is
    %          clipped by the top edge unless the axes is pulled down
    ax.Position = [0.13 0.06 0.74 0.78];
    fig.save2svg(save_dir);
end

function draw_polar_change(fig_name, theta_deg, change_list, cfg_list, name_list, min_n, fontsize, save_dir)
%DRAW_POLAR_CHANGE  Dilated minus constricted. Offset so a negative lobe still plots.
%   caution  polarplot cannot show negative radius. The curves are drawn on
%            1 + change so the unit circle IS zero change; the reported numbers
%            below the figure are the un-offset values
    fig = make_fig(fig_name, 'polar');
    fig.bring_fig();
    fig.update_figsize([4.5 4.5]);
    ax = fig.ax;
    hold(ax, 'on');
    theta = deg2rad([theta_deg, theta_deg(1)]);
    closed = @(r) [r, r(1)];
    support = zeros(1, numel(change_list));

    for list_idx = 1:numel(change_list)
        change = change_list{list_idx};
        cfg    = cfg_list{list_idx};
        for vessel_idx = 1:size(change, 1)
            polarplot(ax, theta, closed(1 + change(vessel_idx, :)), ...
                'Color', cfg.faint, 'LineWidth', 0.5);
        end
        [mean_change, n_change] = mean_supported(change, min_n);
        polarplot(ax, theta, closed(1 + mean_change), ...
            'Color', cfg.bold, 'LineWidth', 2.5);
        support(list_idx) = 100 * mean(n_change >= min_n);
    end
    polarplot(ax, deg2rad(0:360), ones(1, 361), 'k--', 'LineWidth', 0.8);   % zero change
    polarplot(ax, [0 0], [0 ax.RLim(2)], 'k--', 'LineWidth', 0.8);          % alignment axis
    title(ax, sprintf('%s  (1 = no change)', strjoin(name_list, ' / ')));
    subtitle(ax, sprintf('both walls at n \\geq %d : %s', min_n, ...
        strjoin(compose('%s %.0f%%', name_list(:), support(:)), ', ')));
    ax.ThetaTick = 0:45:315;
    set(ax, 'FontName', 'Arial', 'FontSize', fontsize, 'LineWidth', 1);
    ax.Position = [0.13 0.06 0.74 0.78];
    fig.save2svg(save_dir);
end

function report_sector(tag, change, theta_deg, half_deg)
%REPORT_SECTOR  The aligned sector against its opposite, measured bins only.
%   why  a single bin is not a measurement here -- with no fill, almost no vessel
%        has a measured point at exactly 0 or exactly 180 deg. A sector pools
%        whatever that vessel did measure nearby
    offset_0   = mod(theta_deg + 180, 360) - 180;        % signed distance from 0
    offset_180 = mod(theta_deg,       360) - 180;        % signed distance from 180
    sector_0   = abs(offset_0)   <= half_deg;
    sector_180 = abs(offset_180) <= half_deg;
    at_0   = mean(change(:, sector_0),   2, 'omitnan');
    at_180 = mean(change(:, sector_180), 2, 'omitnan');
    paired = ~isnan(at_0) & ~isnan(at_180);
    n = nnz(paired);
    t_val = tinv(0.975, n - 1);
    ci = @(v) t_val * std(v, 'omitnan') / sqrt(n);
    [~, p] = ttest(at_0(paired), at_180(paired));
    fprintf(['%s  0+/-%d deg %+.3f +/- %.3f | 180+/-%d deg %+.3f +/- %.3f | ' ...
             'ratio %.2f | %d of %d vessels higher at 0 | p = %.2g\n'], ...
        tag, half_deg, mean(at_0(paired)), ci(at_0(paired)), ...
        half_deg, mean(at_180(paired)), ci(at_180(paired)), ...
        mean(at_0(paired)) / mean(at_180(paired)), ...
        nnz(at_0(paired) > at_180(paired)), n, p);
end
