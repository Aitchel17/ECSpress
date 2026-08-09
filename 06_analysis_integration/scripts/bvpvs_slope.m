%BVPVS_SLOPE  Per-vessel BV vs total-PVS slope, averaged across vessels by state.
%   The population version of bv_pvs.m, which draws one vessel on one date. Each
%   vessel gets its own robust slope per state; the slopes are what is averaged,
%   not the points.
%
%   why      pooling raw points would weight whichever vessel has the most frames
%            and the widest caliber range. A slope per vessel makes every vessel
%            one sample, the same unit the transition figures count
%   why      both axes are divided by the SAME scalar, that vessel's awake BV
%            median, so the slope stays d(PVS)/d(BV) and is directly comparable
%            to the transition table. Dividing each by its own awake value would
%            turn it into an elasticity that moves with the PVS/BV size ratio
%
%   CAUTION -- the -1 floor. line_fwhm builds
%       bv       = lowerBV - upperBV
%       totalpvs = (upperBV - pvsupedge) + (pvsdownedge - lowerBV) = eps - bv
%   so the BV boundaries enter totalpvs with the OPPOSITE sign. With independent
%   edge noise cov(bv, totalpvs) = cov(bv, eps) - var(bv) = -var(bv), i.e. pure
%   measurement noise drives this slope to -1, which is also what perfect
%   compensation (a fixed outer wall) looks like. The two are not separable on
%   this plot.
%     note   the noise-clean pair is eps vs bv : pvsupedge/pvsdownedge and
%            upperBV/lowerBV are different edges, so their noise is independent
%            and the null slope there is 0. It is exactly this slope + 1, printed
%            alongside, so nothing extra has to be computed
%     note   param.slow_sec re-runs the fit on time-smoothed traces. Boundary noise is
%            frame-to-frame; if the slope moves away from -1 once the fast
%            component is gone, the -1 pull was noise and the remainder is real

clc, clear
addpath(genpath('g:\03_program\01_ecspress\06_analysis_integration\01_tablemanager'));
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
clee = color_lee();

param.dataset = 'merged_igkl_igkltdt';
masterDirTable_path = fullfile( ...
    'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis', ...
    param.dataset, [param.dataset '_dirtable.xlsx']);

param.states      = ["awake", "nrem", "rem"];  % drowsy dropped : too few bouts per vessel
param.min_frame   = 100;                       % frames a (vessel,state) needs to be fitted
param.slow_sec    = 3;                         % moving-median width for the timescale check
param.split_state = "nrem";                    % caliber that divides the constricted/dilated sectors
param.fontsize    = 14;

state_color = struct( ...
    'awake', clee.lch(55,  0,   0), ...   % grey
    'nrem',  clee.lch(55, 80, 240), ...   % skyblue
    'rem',   clee.lch(55, 80,  10));      % scarlet

%% 1. Load and filter, same gates as tablegeneration_main
mtable = tableManager.load_recon(masterDirTable_path, "mtable_FWHMsleep.mat");
mtable.analysis_table = mtable.subTables.state_summary;
mtable.filtLogics = [];
mtable.filtLogics.vType    = contains(mtable.analysis_table.VesselID, "PA", 'IgnoreCase', true);
mtable.filtLogics.depth    = mtable.analysis_table.NumericDepth < 70;
mtable.filtLogics.bout_dur = mtable.analysis_table.bout_duration > 30;
mtable.filtLogics.state    = ismember(mtable.analysis_table.state_name, ["awake","drowsy","nrem","rem"]);
mtable.apply_filter

analyzer = tableAnalyzer(mtable.filtered_table, mtable.action_log);
analyzer.scale_table("NumericResolution", {"raw_data"}, ...
    {'raw_mean','raw_median','raw_q1','raw_q3','raw_var'}, "division");
T = analyzer.filtered_table;

%% 2. Pair the two traces bout by bout
T_bv  = T(string(T.DataType) == "thickness_bv", :);
T_pvs = T(string(T.DataType) == "thickness_totalpvs", :);
key = @(t) string(t.MouseID) + "|" + string(t.Date) + "|" + string(t.VesselID) + "|" + ...
           string(t.state_name) + "|" + string(t.bout_idx);
key_bv  = key(T_bv);
key_pvs = key(T_pvs);
[has_pair, pair_row] = ismember(key_bv, key_pvs);
fprintf('bouts: %d bv rows, %d paired\n', height(T_bv), nnz(has_pair));

paired = struct('vessel_key', {}, 'mouse_id', {}, 'vessel_id', {}, ...
                'state', {}, 'bv', {}, 'pvs', {});
bv_rows = find(has_pair);
for k = 1:numel(bv_rows)
    i = bv_rows(k);
    bv_trace  = T_bv.raw_data{i}(:)';
    pvs_trace = T_pvs.raw_data{pair_row(i)}(:)';
    n = min(numel(bv_trace), numel(pvs_trace));
    paired(end+1) = struct( ...
        'vessel_key', string(T_bv.MouseID(i)) + "_" + string(T_bv.VesselID(i)), ...
        'mouse_id',   string(T_bv.MouseID(i)), ...
        'vessel_id',  string(T_bv.VesselID(i)), ...
        'state',      string(T_bv.state_name(i)), ...
        'bv',         bv_trace(1:n), ...
        'pvs',        pvs_trace(1:n));  %#ok<SAGROW>
end

%% 3. Per vessel : awake BV median is the one divisor for BOTH axes
vessel_names = unique([paired.vessel_key], 'stable');
% caution  a vessel with no awake bout has no baseline and is dropped, not
%          silently rescaled to something else
awake_bv_um = nan(1, numel(vessel_names));
for v = 1:numel(vessel_names)
    sel = [paired.vessel_key] == vessel_names(v) & [paired.state] == "awake";
    if ~any(sel); continue; end
    awake_bv_um(v) = median([paired(sel).bv], 'omitnan');
end
% the caliber that divides the two sectors, same units as the awake baseline
split_bv_um = nan(1, numel(vessel_names));
for v = 1:numel(vessel_names)
    sel = [paired.vessel_key] == vessel_names(v) & [paired.state] == param.split_state;
    if ~any(sel); continue; end
    split_bv_um(v) = median([paired(sel).bv], 'omitnan');
end
fprintf('%d vessels, %d with an awake baseline, %d with a %s split point\n', ...
    numel(vessel_names), nnz(~isnan(awake_bv_um)), nnz(~isnan(split_bv_um)), param.split_state);

%% 3b. Is the sector split just the state split wearing a different name?
% caution  awake / nrem / rem sit at different median calibers, so cutting at the
%          nrem median could put awake entirely left and rem entirely right. If it
%          does, a per-state fit inside a sector is fitting one state
fprintf('\nframes falling LEFT of the %s split, by state\n', param.split_state);
for s = 1:numel(param.states)
    n_left = 0; n_all = 0;
    for v = 1:numel(vessel_names)
        if isnan(split_bv_um(v)); continue; end
        sel = [paired.vessel_key] == vessel_names(v) & [paired.state] == param.states(s);
        bv_all = [paired(sel).bv];
        n_left = n_left + nnz(bv_all < split_bv_um(v));
        n_all  = n_all  + numel(bv_all);
    end
    fprintf('  %-6s %6.1f%% left  (%d frames)\n', param.states(s), 100*n_left/max(n_all,1), n_all);
end

%% 4. Robust slope per (vessel, state), fast and slow
row_proto = struct('mouse_id', "", 'vessel_id', "", 'state', "", ...
    'n_frame', 0, 'awake_bv_um', NaN, 'slope', NaN, 'intercept', NaN, ...
    'slope_slow', NaN, 'bv_pct_iqr', NaN);
slopes = row_proto([]);
fs_guess = 3;                                    % Hz, every sleep session in this set
slow_win = max(3, round(param.slow_sec * fs_guess));

for v = 1:numel(vessel_names)
    if isnan(awake_bv_um(v)); continue; end
    scale = awake_bv_um(v);
    for s = 1:numel(param.states)
        sel = find([paired.vessel_key] == vessel_names(v) & [paired.state] == param.states(s));
        if isempty(sel); continue; end
        bv_pct  = [];  pvs_pct  = [];
        bv_slow = [];  pvs_slow = [];
        for k = sel
            b = paired(k).bv  / scale * 100;
            p = paired(k).pvs / scale * 100;
            bv_pct  = [bv_pct,  b];   pvs_pct  = [pvs_pct,  p];   %#ok<AGROW>
            % smooth WITHIN a bout : a moving median across a bout join would
            % blend two unrelated stretches of recording
            bv_slow  = [bv_slow,  movmedian(b, slow_win)];        %#ok<AGROW>
            pvs_slow = [pvs_slow, movmedian(p, slow_win)];        %#ok<AGROW>
        end
        ok = isfinite(bv_pct) & isfinite(pvs_pct);
        if nnz(ok) < param.min_frame; continue; end

        coef      = robustfit(bv_pct(ok),  pvs_pct(ok));
        ok_slow   = isfinite(bv_slow) & isfinite(pvs_slow);
        coef_slow = robustfit(bv_slow(ok_slow), pvs_slow(ok_slow));

        row = row_proto;
        parts = split(vessel_names(v), "_");
        row.mouse_id    = parts(1);
        row.vessel_id   = parts(2);
        row.state       = param.states(s);
        row.n_frame     = nnz(ok);
        row.awake_bv_um = scale;
        row.intercept   = coef(1);
        row.slope       = coef(2);
        row.slope_slow  = coef_slow(2);
        row.bv_pct_iqr  = diff(prctile(bv_pct(ok), [25 75]));
        slopes(end+1) = row;  %#ok<SAGROW>
    end
end

%% 4b. Sector slopes : constricted vs dilated relative to the param.split_state caliber
% why      one line per vessel per state made too many short segments to read as
%          anything but noise. The question underneath is whether the compensation
%          ratio depends on how dilated the vessel already is, which is a caliber
%          split, not a state split
% see FINDINGS.md
% caution  rem sits almost entirely in the dilated sector, so a sector difference
%          could be a rem difference wearing a caliber label. SECTOR_STATES runs
%          the whole thing
%          a second time on awake+nrem only -- if the split survives there it is
%          caliber, not state
sector_proto = struct('mouse_id', "", 'vessel_id', "", 'sector', "", 'stateset', "", ...
    'n_frame', 0, 'slope', NaN, 'slope_slow', NaN, 'span_pct', NaN);
sectors = sector_proto([]);
statesets = {param.states, ["awake","nrem"]};
statesetnames = ["all", "awake+nrem"];
for ss = 1:numel(statesets)
for v = 1:numel(vessel_names)
    if isnan(awake_bv_um(v)) || isnan(split_bv_um(v)); continue; end
    scale = awake_bv_um(v);
    split_pct = split_bv_um(v) / scale * 100;
    sel = find([paired.vessel_key] == vessel_names(v) & ismember([paired.state], statesets{ss}));
    bv_pct = [];  pvs_pct = [];  bv_slow = [];  pvs_slow = [];
    for k = sel
        b = paired(k).bv / scale * 100;
        p = paired(k).pvs / scale * 100;
        bv_pct = [bv_pct, b];  pvs_pct = [pvs_pct, p];                 %#ok<AGROW>
        bv_slow = [bv_slow, movmedian(b, slow_win)];                   %#ok<AGROW>
        pvs_slow = [pvs_slow, movmedian(p, slow_win)];                 %#ok<AGROW>
    end
    for side = ["constricted", "dilated"]
        if side == "constricted"
            in_sector = bv_pct < split_pct;
        else
            in_sector = bv_pct >= split_pct;
        end
        ok = in_sector & isfinite(bv_pct) & isfinite(pvs_pct);
        if nnz(ok) < param.min_frame; continue; end
        coef = robustfit(bv_pct(ok) - split_pct, pvs_pct(ok));
        ok_slow = in_sector & isfinite(bv_slow) & isfinite(pvs_slow);
        coef_slow = robustfit(bv_slow(ok_slow) - split_pct, pvs_slow(ok_slow));

        parts = split(vessel_names(v), "_");
        row = sector_proto;
        row.mouse_id  = parts(1);   row.vessel_id = parts(2);
        row.sector    = side;       row.n_frame   = nnz(ok);
        row.stateset  = statesetnames(ss);
        row.slope     = coef(2);    row.slope_slow = coef_slow(2);
        row.span_pct  = prctile(abs(bv_pct(ok) - split_pct), 90);
        sectors(end+1) = row;  %#ok<SAGROW>
    end
end
end   % ss, the state-set control loop

for ss = 1:numel(statesetnames)
    fprintf('\nstates = %s\n', statesetnames(ss));
    fprintf('%-12s %4s %10s %12s %12s %14s %10s\n', ...
        'sector', 'n', 'slope', '95%% CI', 'p vs -1', 'slope+1 (eps)', 'slow');
    for side = ["constricted", "dilated"]
        sel = [sectors.sector] == side & [sectors.stateset] == statesetnames(ss);
        v = [sectors(sel).slope];  v_slow = [sectors(sel).slope_slow];
        n = numel(v);
        if n < 2; continue; end
        ci = tinv(0.975, n-1) * std(v) / sqrt(n);
        [~, p1] = ttest(v, -1);
        fprintf('%-12s %4d %10.3f %12s %12.2g %14.3f %10.3f\n', ...
            side, n, mean(v), sprintf('+/-%.3f', ci), p1, mean(v) + 1, mean(v_slow));
    end
    in_set = [sectors.stateset] == statesetnames(ss);
    c_sel = in_set & [sectors.sector] == "constricted";
    d_sel = in_set & [sectors.sector] == "dilated";
    kc = [sectors(c_sel).mouse_id] + "_" + [sectors(c_sel).vessel_id];
    kd = [sectors(d_sel).mouse_id] + "_" + [sectors(d_sel).vessel_id];
    vc = [sectors(c_sel).slope];   vd = [sectors(d_sel).slope];
    [~, ic, id] = intersect(kc, kd, 'stable');
    if numel(ic) < 2; continue; end
    [~, p_sector] = ttest(vc(ic), vd(id));
    fprintf('paired constricted vs dilated  n=%d  %+.3f vs %+.3f  d=%+.3f  p=%.3g\n', ...
        numel(ic), mean(vc(ic)), mean(vd(id)), mean(vd(id) - vc(ic)), p_sector);
end

%% 5. Average the SLOPES, one vessel one vote
fprintf('\n%-6s %4s %10s %12s %12s %14s %10s\n', ...
    'state', 'n', 'slope', '95%% CI', 'p vs -1', 'slope+1 (eps)', 'slow');
summary = struct('state', {}, 'n', {}, 'mean_slope', {}, 'ci', {}, 'p_vs_minus1', {});
for s = 1:numel(param.states)
    sel = [slopes.state] == param.states(s);
    v = [slopes(sel).slope];
    v_slow = [slopes(sel).slope_slow];
    n = numel(v);
    if n < 2; continue; end
    ci = tinv(0.975, n-1) * std(v) / sqrt(n);
    [~, p1] = ttest(v, -1);
    fprintf('%-6s %4d %10.3f %12s %12.2g %14.3f %10.3f\n', ...
        param.states(s), n, mean(v), sprintf('+/-%.3f', ci), p1, mean(v) + 1, mean(v_slow));
    summary(end+1) = struct('state', param.states(s), 'n', n, 'mean_slope', mean(v), ...
        'ci', ci, 'p_vs_minus1', p1);  %#ok<SAGROW>
end

% paired state comparisons on the vessels that have both
fprintf('\npaired state comparisons (vessels with both states)\n');
for a = 1:numel(param.states)
    for b = a+1:numel(param.states)
        ka = [slopes([slopes.state] == param.states(a)).mouse_id] + "_" + ...
             [slopes([slopes.state] == param.states(a)).vessel_id];
        kb = [slopes([slopes.state] == param.states(b)).mouse_id] + "_" + ...
             [slopes([slopes.state] == param.states(b)).vessel_id];
        va = [slopes([slopes.state] == param.states(a)).slope];
        vb = [slopes([slopes.state] == param.states(b)).slope];
        [common, ia, ib] = intersect(ka, kb, 'stable');
        if numel(common) < 2; continue; end
        [~, p] = ttest(va(ia), vb(ib));
        fprintf('  %-6s vs %-6s  n=%2d  %+.3f vs %+.3f  d=%+.3f  p=%.3g\n', ...
            param.states(a), param.states(b), numel(common), mean(va(ia)), mean(vb(ib)), ...
            mean(vb(ib) - va(ia)), p);
    end
end

%% 6. Figure : every vessel's fitted line faint, the mean slope bold
% why  each line is drawn through ITS OWN (BV, PVS) operating point, i.e. as a
%      deviation. Plotting the raw intercepts spread the panel over a range that
%      buried the slopes (see FINDINGS.md) -- the PVS/BV size ratio differs
%      between vessels and that vertical offset is not what was averaged.
%      Centring leaves only the slope visible
% caution  close old figures first : make_fig opens a NEW window per run and
%          several windows sharing a Name silently overwrite each other on export
close(findobj(groot, 'Type', 'figure', 'Name', 'bvpvs_slope'));
fig = make_fig('bvpvs_slope', 'normal');
fig.bring_fig(); fig.update_figsize([4.0 4.0]);
ax = fig.ax; cla(ax); hold(ax, 'on');

% One line per vessel per SECTOR, each spanning only its own side of the split, so
% the two fans separate instead of 68 segments overlapping through the origin
sides      = ["constricted", "dilated"];
side_faint = {clee.lch(80, 60, 250), clee.lch(80, 90,  40)};
side_bold  = {clee.lch(45, 80, 250), clee.lch(45, 110, 40)};
for si = 1:2
    for k = find([sectors.sector] == sides(si) & [sectors.stateset] == "all")
        if sides(si) == "constricted"
            dx = [-sectors(k).span_pct, 0];
        else
            dx = [0, sectors(k).span_pct];
        end
        plot(ax, dx, sectors(k).slope*dx, ...
            'Color', side_faint{si}, 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end
mean_handles = gobjects(1, 2);
for si = 1:2
    sel = [sectors.sector] == sides(si) & [sectors.stateset] == "all";
    if ~any(sel); continue; end
    span = mean([sectors(sel).span_pct]);
    if sides(si) == "constricted"; dx = [-span, 0]; else; dx = [0, span]; end
    mean_handles(si) = plot(ax, dx, mean([sectors(sel).slope])*dx, ...
        'Color', side_bold{si}, 'LineWidth', 3, ...
        'DisplayName', sprintf('%s  %.2f  (n=%d)', sides(si), ...
            mean([sectors(sel).slope]), nnz(sel)));
end
% the two references this plot is read against
plot(ax, [-40 40], [40 -40], 'k--', 'LineWidth', 1, ...
    'DisplayName', 'slope -1 (fixed outer wall)');
yline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
xline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');

xlabel(ax, sprintf('\\DeltaBV from %s caliber (%% of awake BV)', param.split_state));
ylabel(ax, '\DeltaTotal PVS (% of awake BV)');
xlim(ax, [-40 40]);  ylim(ax, [-40 40]);
set(ax, 'Box','off', 'TickDir','out', 'LineWidth',1, ...
    'FontName','Arial', 'FontSize',param.fontsize, 'Layer','top');
legend(ax, [mean_handles(isgraphics(mean_handles)), ax.Children(1)], ...
    'Location', 'southwest', 'Box', 'off', 'FontSize', param.fontsize-3);

slope_table = struct2table(slopes);
disp(slope_table)
